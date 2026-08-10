#!/usr/bin/env python3
"""
Controlled Romanian TTS -> SpeD ASR benchmark.

This script:
1. Synthesizes the same Romanian prompts with several online TTS services
   and two local baselines.
2. Converts all audio to 16 kHz, mono, PCM16 WAV.
3. Transcribes it with the exact CTC-greedy ONNX export of
   gabrielpirlo/SpeD_ParakeetRo_110M_TDT-CTC.
4. Computes normalized WER/CER, accent-insensitive metrics, and CPU latency.
5. Writes CSV, JSON, Markdown, and a self-contained HTML report with audio players.

Synthetic speech is a controlled pronunciation test, not a substitute for
real-microphone evaluation.
"""

from __future__ import annotations

import asyncio
import csv
import hashlib
import html
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
import traceback
import unicodedata
import wave
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import edge_tts
import onnx_asr
import requests
from gtts import gTTS
from huggingface_hub import snapshot_download


ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = ROOT / "results"
AUDIO_DIR = RESULTS_DIR / "audio"
SOURCE_AUDIO_DIR = AUDIO_DIR / "source"
WAV_AUDIO_DIR = AUDIO_DIR / "wav_16k_mono"
MODEL_DIR = ROOT / ".cache" / "sped-ro-onnx"
PIPER_DIR = ROOT / ".cache" / "piper"
LOG_PATH = RESULTS_DIR / "benchmark.log"

MODEL_REPO = "AlinClaudiu/SpeD-ParakeetRo-110M-onnx"
MODEL_TYPE = "nemo-conformer-ctc"
EXPECTED_SHA256 = {
    "model.onnx": "c5ffc434d9f16d9ba257f380b3634b7ca19ef9b5b7420489705669539e2a640d",
    "model.onnx.data": "81e884781377a1365802777295dec05c068e4005be9f657a5519994a5071544b",
}
PIPER_VOICE = "ro_RO-mihai-medium"

PROMPTS: list[dict[str, str]] = [
    {
        "id": "clean_prose",
        "category": "Clean Romanian prose",
        "text": (
            "În această dimineață, cerul este senin, iar temperatura va crește "
            "treptat până la douăzeci și șapte de grade."
        ),
    },
    {
        "id": "diacritics",
        "category": "Diacritics and difficult phonemes",
        "text": (
            "Ștefan și Țiței au cumpărat pâine proaspătă, brânză și câteva "
            "căpșuni pentru călătoria de mâine."
        ),
    },
    {
        "id": "numbers_dates",
        "category": "Numbers, date, and time",
        "text": (
            "Întâlnirea este programată pe paisprezece august două mii douăzeci "
            "și șase, la ora cincisprezece și treizeci de minute, în sala două sute patru."
        ),
    },
    {
        "id": "technical",
        "category": "Romanian technical dictation",
        "text": (
            "Serverul Proxmox rulează trei mașini virtuale, iar baza de date "
            "PostgreSQL ascultă pe portul cinci mii patru sute treizeci și doi."
        ),
    },
    {
        "id": "code_switching",
        "category": "Romanian-English code-switching",
        "text": (
            "Trebuie să facem rollback la deployment și să verificăm health check-ul "
            "din Kubernetes pentru PhotonSpark."
        ),
    },
    {
        "id": "self_correction",
        "category": "Spoken self-correction",
        "text": (
            "Trimite raportul joi după-amiază, de fapt vineri dimineață, înainte de ora nouă."
        ),
    },
]

SERVICE_METADATA: dict[str, dict[str, str]] = {
    "edge_alina": {
        "label": "Microsoft Edge — Alina",
        "kind": "online neural",
        "voice": "ro-RO-AlinaNeural",
    },
    "edge_emil": {
        "label": "Microsoft Edge — Emil",
        "kind": "online neural",
        "voice": "ro-RO-EmilNeural",
    },
    "google_gtts": {
        "label": "Google Translate TTS (gTTS)",
        "kind": "online",
        "voice": "Romanian / ro",
    },
    "kiprio": {
        "label": "Kiprio demo API",
        "kind": "online",
        "voice": "Romanian / ro",
    },
    "piper_mihai": {
        "label": "Piper — Mihai",
        "kind": "local open-source neural",
        "voice": PIPER_VOICE,
    },
    "espeak_ng": {
        "label": "eSpeak NG — Romanian",
        "kind": "local formant baseline",
        "voice": "ro",
    },
}


@dataclass
class GenerationRecord:
    service: str
    prompt_id: str
    status: str
    source_path: str | None = None
    wav_path: str | None = None
    error: str | None = None


@dataclass
class BenchmarkRecord:
    service: str
    service_label: str
    service_kind: str
    voice: str
    prompt_id: str
    category: str
    reference: str
    hypothesis_raw: str
    reference_normalized: str
    hypothesis_normalized: str
    reference_accentless: str
    hypothesis_accentless: str
    wer: float
    cer: float
    wer_accentless: float
    cer_accentless: float
    duration_s: float
    inference_s: float
    rtf: float
    rtfx: float
    wav_path: str
    source_path: str
    status: str = "ok"
    error: str | None = None


def log(message: str) -> None:
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    line = f"[{stamp}] {message}"
    print(line, flush=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def ensure_directories() -> None:
    for path in (RESULTS_DIR, AUDIO_DIR, SOURCE_AUDIO_DIR, WAV_AUDIO_DIR, MODEL_DIR, PIPER_DIR):
        path.mkdir(parents=True, exist_ok=True)


def run_command(
    command: list[str],
    *,
    input_text: str | None = None,
    check: bool = True,
    timeout: int = 300,
) -> subprocess.CompletedProcess[str]:
    log("$ " + " ".join(command))
    return subprocess.run(
        command,
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=check,
        timeout=timeout,
    )


def convert_to_wav(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    result = run_command(
        [
            "ffmpeg",
            "-hide_banner",
            "-loglevel",
            "error",
            "-y",
            "-i",
            str(source),
            "-ac",
            "1",
            "-ar",
            "16000",
            "-c:a",
            "pcm_s16le",
            str(destination),
        ],
        timeout=180,
    )
    if result.stdout.strip():
        log(result.stdout.strip())
    if not destination.exists() or destination.stat().st_size < 100:
        raise RuntimeError(f"ffmpeg did not produce a valid WAV: {destination}")


def wav_duration(path: Path) -> float:
    with wave.open(str(path), "rb") as handle:
        frames = handle.getnframes()
        rate = handle.getframerate()
        channels = handle.getnchannels()
        width = handle.getsampwidth()
    if rate != 16000 or channels != 1 or width != 2:
        raise ValueError(
            f"Unexpected normalized WAV format for {path}: "
            f"rate={rate}, channels={channels}, sample_width={width}"
        )
    return frames / rate


async def generate_edge_voice(voice: str, text: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: BaseException | None = None
    for attempt in range(1, 4):
        try:
            communicate = edge_tts.Communicate(text=text, voice=voice)
            await communicate.save(str(destination))
            if destination.exists() and destination.stat().st_size > 100:
                return
            raise RuntimeError("edge-tts returned an empty audio file")
        except BaseException as exc:
            last_error = exc
            log(f"edge-tts {voice} attempt {attempt}/3 failed: {exc}")
            await asyncio.sleep(2**attempt)
    raise RuntimeError(f"edge-tts failed after retries: {last_error}")


def generate_google_tts(text: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: BaseException | None = None
    for attempt in range(1, 4):
        try:
            gTTS(text=text, lang="ro", slow=False).save(str(destination))
            if destination.exists() and destination.stat().st_size > 100:
                return
            raise RuntimeError("gTTS returned an empty audio file")
        except BaseException as exc:
            last_error = exc
            log(f"gTTS attempt {attempt}/3 failed: {exc}")
            time.sleep(2**attempt)
    raise RuntimeError(f"gTTS failed after retries: {last_error}")


def generate_kiprio(text: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    url = "https://kiprio.com/v1/tts/demo"
    last_error: BaseException | None = None
    for attempt in range(1, 5):
        try:
            response = requests.post(
                url,
                headers={
                    "Content-Type": "application/json",
                    "Accept": "audio/mpeg,application/octet-stream;q=0.9,*/*;q=0.1",
                    "User-Agent": "SpeD-Romanian-TTS-Benchmark/1.0",
                },
                json={"text": text, "lang": "ro", "slow": False},
                timeout=90,
            )
            if response.status_code == 429:
                retry_after = int(response.headers.get("Retry-After", "20"))
                log(f"Kiprio rate-limited; sleeping {retry_after}s")
                time.sleep(max(20, retry_after))
                continue
            response.raise_for_status()
            content_type = response.headers.get("Content-Type", "")
            if len(response.content) < 100:
                raise RuntimeError(
                    f"Kiprio response too small ({len(response.content)} bytes, {content_type})"
                )
            if "json" in content_type.lower():
                payload = response.json()
                audio_b64 = payload.get("audio") or payload.get("audio_base64")
                if not audio_b64:
                    raise RuntimeError(f"Unexpected Kiprio JSON response keys: {list(payload)}")
                import base64

                destination.write_bytes(base64.b64decode(audio_b64))
            else:
                destination.write_bytes(response.content)
            return
        except BaseException as exc:
            last_error = exc
            log(f"Kiprio attempt {attempt}/4 failed: {exc}")
            time.sleep(min(30, 3 * attempt))
    raise RuntimeError(f"Kiprio failed after retries: {last_error}")


def ensure_piper_voice() -> None:
    model_path = PIPER_DIR / f"{PIPER_VOICE}.onnx"
    config_path = PIPER_DIR / f"{PIPER_VOICE}.onnx.json"
    if model_path.exists() and config_path.exists():
        return
    result = run_command(
        [
            sys.executable,
            "-m",
            "piper.download_voices",
            "--data-dir",
            str(PIPER_DIR),
            PIPER_VOICE,
        ],
        check=False,
        timeout=600,
    )
    if result.stdout.strip():
        log(result.stdout.strip())
    if not model_path.exists() or not config_path.exists():
        raise RuntimeError(
            f"Piper voice download did not produce {model_path.name} and {config_path.name}"
        )


def generate_piper(text: str, destination: Path) -> None:
    ensure_piper_voice()
    destination.parent.mkdir(parents=True, exist_ok=True)
    model_path = PIPER_DIR / f"{PIPER_VOICE}.onnx"
    result = run_command(
        [
            sys.executable,
            "-m",
            "piper",
            "--model",
            str(model_path),
            "--output_file",
            str(destination),
            "--",
            text,
        ],
        check=False,
        timeout=180,
    )
    if result.stdout.strip():
        log(result.stdout.strip())
    if result.returncode != 0 or not destination.exists() or destination.stat().st_size < 100:
        raise RuntimeError(f"Piper generation failed with exit code {result.returncode}")


def generate_espeak(text: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    result = run_command(
        [
            "espeak-ng",
            "-v",
            "ro",
            "-s",
            "165",
            "-w",
            str(destination),
            text,
        ],
        check=False,
        timeout=60,
    )
    if result.stdout.strip():
        log(result.stdout.strip())
    if result.returncode != 0 or not destination.exists() or destination.stat().st_size < 100:
        raise RuntimeError(f"eSpeak generation failed with exit code {result.returncode}")


def source_extension(service: str) -> str:
    if service in {"piper_mihai", "espeak_ng"}:
        return ".wav"
    return ".mp3"


async def generate_all_audio() -> list[GenerationRecord]:
    records: list[GenerationRecord] = []

    for service, voice in (
        ("edge_alina", "ro-RO-AlinaNeural"),
        ("edge_emil", "ro-RO-EmilNeural"),
    ):
        semaphore = asyncio.Semaphore(2)

        async def one_edge(prompt: dict[str, str]) -> GenerationRecord:
            source = SOURCE_AUDIO_DIR / service / f"{prompt['id']}.mp3"
            wav = WAV_AUDIO_DIR / service / f"{prompt['id']}.wav"
            try:
                async with semaphore:
                    await generate_edge_voice(voice, prompt["text"], source)
                convert_to_wav(source, wav)
                return GenerationRecord(
                    service=service,
                    prompt_id=prompt["id"],
                    status="ok",
                    source_path=str(source.relative_to(RESULTS_DIR)),
                    wav_path=str(wav.relative_to(RESULTS_DIR)),
                )
            except BaseException as exc:
                return GenerationRecord(
                    service=service,
                    prompt_id=prompt["id"],
                    status="failed",
                    error=f"{type(exc).__name__}: {exc}",
                )

        records.extend(await asyncio.gather(*(one_edge(prompt) for prompt in PROMPTS)))

    for prompt in PROMPTS:
        service = "google_gtts"
        source = SOURCE_AUDIO_DIR / service / f"{prompt['id']}.mp3"
        wav = WAV_AUDIO_DIR / service / f"{prompt['id']}.wav"
        try:
            generate_google_tts(prompt["text"], source)
            convert_to_wav(source, wav)
            records.append(
                GenerationRecord(
                    service=service,
                    prompt_id=prompt["id"],
                    status="ok",
                    source_path=str(source.relative_to(RESULTS_DIR)),
                    wav_path=str(wav.relative_to(RESULTS_DIR)),
                )
            )
        except BaseException as exc:
            records.append(
                GenerationRecord(
                    service=service,
                    prompt_id=prompt["id"],
                    status="failed",
                    error=f"{type(exc).__name__}: {exc}",
                )
            )

    for index, prompt in enumerate(PROMPTS):
        service = "kiprio"
        source = SOURCE_AUDIO_DIR / service / f"{prompt['id']}.mp3"
        wav = WAV_AUDIO_DIR / service / f"{prompt['id']}.wav"
        try:
            generate_kiprio(prompt["text"], source)
            convert_to_wav(source, wav)
            records.append(
                GenerationRecord(
                    service=service,
                    prompt_id=prompt["id"],
                    status="ok",
                    source_path=str(source.relative_to(RESULTS_DIR)),
                    wav_path=str(wav.relative_to(RESULTS_DIR)),
                )
            )
        except BaseException as exc:
            records.append(
                GenerationRecord(
                    service=service,
                    prompt_id=prompt["id"],
                    status="failed",
                    error=f"{type(exc).__name__}: {exc}",
                )
            )
        if index < len(PROMPTS) - 1:
            time.sleep(13)

    for service, generator in (
        ("piper_mihai", generate_piper),
        ("espeak_ng", generate_espeak),
    ):
        for prompt in PROMPTS:
            extension = source_extension(service)
            source = SOURCE_AUDIO_DIR / service / f"{prompt['id']}{extension}"
            wav = WAV_AUDIO_DIR / service / f"{prompt['id']}.wav"
            try:
                generator(prompt["text"], source)
                convert_to_wav(source, wav)
                records.append(
                    GenerationRecord(
                        service=service,
                        prompt_id=prompt["id"],
                        status="ok",
                        source_path=str(source.relative_to(RESULTS_DIR)),
                        wav_path=str(wav.relative_to(RESULTS_DIR)),
                    )
                )
            except BaseException as exc:
                records.append(
                    GenerationRecord(
                        service=service,
                        prompt_id=prompt["id"],
                        status="failed",
                        error=f"{type(exc).__name__}: {exc}",
                    )
                )
                if service == "piper_mihai":
                    break

    (RESULTS_DIR / "generation.json").write_text(
        json.dumps([asdict(item) for item in records], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return records


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_model() -> dict[str, Any]:
    log(f"Downloading exact ONNX export: {MODEL_REPO}")
    snapshot_download(
        repo_id=MODEL_REPO,
        local_dir=str(MODEL_DIR),
        allow_patterns=["model.onnx", "model.onnx.data", "vocab.txt", "config.json", "README.md"],
    )

    checksums: dict[str, Any] = {}
    for filename, expected in EXPECTED_SHA256.items():
        path = MODEL_DIR / filename
        if not path.exists():
            raise FileNotFoundError(f"Required model file missing: {path}")
        actual = sha256_file(path)
        checksums[filename] = {
            "expected_sha256": expected,
            "actual_sha256": actual,
            "matched": actual == expected,
            "bytes": path.stat().st_size,
        }
        if actual != expected:
            raise RuntimeError(
                f"Checksum mismatch for {filename}: expected {expected}, got {actual}"
            )

    for filename in ("vocab.txt", "config.json"):
        path = MODEL_DIR / filename
        checksums[filename] = {
            "actual_sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        }

    (RESULTS_DIR / "model_checksums.json").write_text(
        json.dumps(checksums, indent=2),
        encoding="utf-8",
    )
    return checksums


_LEGACY_DIACRITICS = str.maketrans({"ş": "ș", "Ş": "Ș", "ţ": "ț", "Ţ": "Ț"})


def normalize_text(text: str, *, strip_diacritics: bool = False) -> str:
    text = unicodedata.normalize("NFC", text.translate(_LEGACY_DIACRITICS)).lower()
    text = text.replace("‑", "-").replace("–", "-").replace("—", "-")
    text = text.replace("-", " ")
    text = re.sub(r"[^\w\săâîșț]", " ", text, flags=re.UNICODE)
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if strip_diacritics:
        text = "".join(
            char
            for char in unicodedata.normalize("NFD", text)
            if not unicodedata.combining(char)
        )
        text = unicodedata.normalize("NFC", text)
    return text


def levenshtein_distance(reference: list[str], hypothesis: list[str]) -> int:
    if len(reference) < len(hypothesis):
        reference, hypothesis = hypothesis, reference
    previous = list(range(len(hypothesis) + 1))
    for i, ref_item in enumerate(reference, start=1):
        current = [i]
        for j, hyp_item in enumerate(hypothesis, start=1):
            substitution = previous[j - 1] + (ref_item != hyp_item)
            insertion = current[j - 1] + 1
            deletion = previous[j] + 1
            current.append(min(substitution, insertion, deletion))
        previous = current
    return previous[-1]


def error_rate(reference: Iterable[str], hypothesis: Iterable[str]) -> float:
    ref = list(reference)
    hyp = list(hypothesis)
    if not ref:
        return 0.0 if not hyp else 1.0
    return levenshtein_distance(ref, hyp) / len(ref)


def compute_metrics(reference: str, hypothesis: str) -> dict[str, Any]:
    ref_norm = normalize_text(reference)
    hyp_norm = normalize_text(hypothesis)
    ref_accentless = normalize_text(reference, strip_diacritics=True)
    hyp_accentless = normalize_text(hypothesis, strip_diacritics=True)

    return {
        "reference_normalized": ref_norm,
        "hypothesis_normalized": hyp_norm,
        "reference_accentless": ref_accentless,
        "hypothesis_accentless": hyp_accentless,
        "wer": error_rate(ref_norm.split(), hyp_norm.split()),
        "cer": error_rate(list(ref_norm.replace(" ", "")), list(hyp_norm.replace(" ", ""))),
        "wer_accentless": error_rate(ref_accentless.split(), hyp_accentless.split()),
        "cer_accentless": error_rate(
            list(ref_accentless.replace(" ", "")),
            list(hyp_accentless.replace(" ", "")),
        ),
    }


def recognize_all(
    generation: list[GenerationRecord],
) -> tuple[list[BenchmarkRecord], float, dict[str, Any]]:
    available = [item for item in generation if item.status == "ok" and item.wav_path]
    if not available:
        raise RuntimeError("No audio was generated successfully; cannot run ASR benchmark")

    checksums = prepare_model()

    log(f"Loading {MODEL_TYPE} from {MODEL_DIR}")
    load_started = time.perf_counter()
    model = onnx_asr.load_model(
        MODEL_TYPE,
        str(MODEL_DIR),
        providers=["CPUExecutionProvider"],
    )
    model_load_s = time.perf_counter() - load_started
    log(f"Model loaded in {model_load_s:.3f}s")

    prompt_by_id = {item["id"]: item for item in PROMPTS}
    records: list[BenchmarkRecord] = []

    warmup_wav = RESULTS_DIR / available[0].wav_path
    log(f"Warming up ASR with {warmup_wav}")
    _ = model.recognize(str(warmup_wav))

    for item in available:
        prompt = prompt_by_id[item.prompt_id]
        wav = RESULTS_DIR / str(item.wav_path)
        duration = wav_duration(wav)
        started = time.perf_counter()
        try:
            hypothesis = model.recognize(str(wav))
            inference_s = time.perf_counter() - started
            hypothesis = str(hypothesis).strip()
            metrics = compute_metrics(prompt["text"], hypothesis)
            metadata = SERVICE_METADATA[item.service]
            records.append(
                BenchmarkRecord(
                    service=item.service,
                    service_label=metadata["label"],
                    service_kind=metadata["kind"],
                    voice=metadata["voice"],
                    prompt_id=item.prompt_id,
                    category=prompt["category"],
                    reference=prompt["text"],
                    hypothesis_raw=hypothesis,
                    duration_s=duration,
                    inference_s=inference_s,
                    rtf=inference_s / duration if duration else float("inf"),
                    rtfx=duration / inference_s if inference_s else float("inf"),
                    wav_path=str(item.wav_path),
                    source_path=str(item.source_path),
                    **metrics,
                )
            )
            log(
                f"{item.service}/{item.prompt_id}: "
                f"WER={metrics['wer']:.3f}, CER={metrics['cer']:.3f}, "
                f"duration={duration:.2f}s, inference={inference_s:.3f}s"
            )
        except BaseException as exc:
            inference_s = time.perf_counter() - started
            metadata = SERVICE_METADATA[item.service]
            records.append(
                BenchmarkRecord(
                    service=item.service,
                    service_label=metadata["label"],
                    service_kind=metadata["kind"],
                    voice=metadata["voice"],
                    prompt_id=item.prompt_id,
                    category=prompt["category"],
                    reference=prompt["text"],
                    hypothesis_raw="",
                    reference_normalized=normalize_text(prompt["text"]),
                    hypothesis_normalized="",
                    reference_accentless=normalize_text(prompt["text"], strip_diacritics=True),
                    hypothesis_accentless="",
                    wer=1.0,
                    cer=1.0,
                    wer_accentless=1.0,
                    cer_accentless=1.0,
                    duration_s=duration,
                    inference_s=inference_s,
                    rtf=inference_s / duration if duration else float("inf"),
                    rtfx=duration / inference_s if inference_s else float("inf"),
                    wav_path=str(item.wav_path),
                    source_path=str(item.source_path),
                    status="failed",
                    error=f"{type(exc).__name__}: {exc}",
                )
            )
            log(f"Recognition failed for {item.service}/{item.prompt_id}: {exc}")

    metadata = {
        "model_repo": MODEL_REPO,
        "model_type": MODEL_TYPE,
        "decoder": "CTC greedy (onnx-asr)",
        "model_load_s": model_load_s,
        "checksums": checksums,
    }
    return records, model_load_s, metadata


def mean_or_none(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def aggregate_records(records: list[BenchmarkRecord], generation: list[GenerationRecord]) -> list[dict[str, Any]]:
    aggregates: list[dict[str, Any]] = []
    generation_by_service: dict[str, list[GenerationRecord]] = {}
    for item in generation:
        generation_by_service.setdefault(item.service, []).append(item)

    for service, metadata in SERVICE_METADATA.items():
        valid = [item for item in records if item.service == service and item.status == "ok"]
        generated = generation_by_service.get(service, [])
        aggregates.append(
            {
                "service": service,
                "service_label": metadata["label"],
                "service_kind": metadata["kind"],
                "voice": metadata["voice"],
                "generated_ok": sum(item.status == "ok" for item in generated),
                "generated_failed": sum(item.status != "ok" for item in generated),
                "recognized_ok": len(valid),
                "mean_wer": mean_or_none([item.wer for item in valid]),
                "mean_cer": mean_or_none([item.cer for item in valid]),
                "mean_wer_accentless": mean_or_none([item.wer_accentless for item in valid]),
                "mean_cer_accentless": mean_or_none([item.cer_accentless for item in valid]),
                "mean_inference_s": mean_or_none([item.inference_s for item in valid]),
                "mean_rtf": mean_or_none([item.rtf for item in valid]),
                "mean_rtfx": mean_or_none([item.rtfx for item in valid]),
            }
        )
    return aggregates


def format_pct(value: float | None) -> str:
    return "—" if value is None else f"{100 * value:.2f}%"


def format_number(value: float | None, digits: int = 3) -> str:
    return "—" if value is None else f"{value:.{digits}f}"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(
    records: list[BenchmarkRecord],
    generation: list[GenerationRecord],
    model_metadata: dict[str, Any],
) -> None:
    record_rows = [asdict(item) for item in records]
    aggregates = aggregate_records(records, generation)

    write_csv(RESULTS_DIR / "results.csv", record_rows)
    write_csv(RESULTS_DIR / "summary_by_voice.csv", aggregates)

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "methodology": {
            "scope": "Synthetic Romanian TTS pronunciation stress test",
            "warning": (
                "Synthetic TTS speech is not representative of spontaneous speech, "
                "microphone noise, room acoustics, disfluencies, or real users."
            ),
            "audio_format": "16 kHz mono PCM16 WAV",
            "metric_normalization": (
                "Unicode NFC, lowercase, punctuation removed, hyphens treated as token "
                "boundaries; both diacritic-sensitive and accent-insensitive WER/CER are reported."
            ),
        },
        "model": model_metadata,
        "environment": collect_environment(),
        "prompts": PROMPTS,
        "generation": [asdict(item) for item in generation],
        "results": record_rows,
        "summary_by_voice": aggregates,
    }
    (RESULTS_DIR / "results.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    write_markdown_report(payload)
    write_html_report(payload)
    write_environment_file(payload["environment"])


def collect_environment() -> dict[str, Any]:
    env: dict[str, Any] = {
        "python": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "github_run_id": os.getenv("GITHUB_RUN_ID"),
        "github_run_attempt": os.getenv("GITHUB_RUN_ATTEMPT"),
        "github_sha": os.getenv("GITHUB_SHA"),
        "github_ref": os.getenv("GITHUB_REF"),
    }
    try:
        from importlib import metadata

        env["onnx_asr"] = metadata.version("onnx-asr")
    except BaseException:
        env["onnx_asr"] = None
    try:
        env["onnxruntime"] = __import__("onnxruntime").__version__
    except BaseException:
        env["onnxruntime"] = None
    for command_name in ("ffmpeg", "espeak-ng"):
        try:
            result = run_command([command_name, "--version"], check=False, timeout=20)
            env[command_name] = result.stdout.splitlines()[0] if result.stdout else None
        except BaseException as exc:
            env[command_name] = f"unavailable: {exc}"
    return env


def write_environment_file(environment: dict[str, Any]) -> None:
    lines = [f"{key}: {value}" for key, value in environment.items()]
    (RESULTS_DIR / "environment.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_markdown_report(payload: dict[str, Any]) -> None:
    aggregates = payload["summary_by_voice"]
    results = payload["results"]
    generation = payload["generation"]

    lines = [
        "# SpeD Romanian Parakeet 110M — Romanian TTS benchmark",
        "",
        f"Generated: `{payload['generated_at']}`",
        "",
        "## Important limitation",
        "",
        (
            "This is a **controlled synthetic-speech pronunciation test**, not a real-world "
            "dictation benchmark. Neural TTS generally has clean acoustics, regular pacing, "
            "and fewer hesitations than human microphone audio. Results can identify model "
            "vocabulary/pronunciation weaknesses, but cannot establish production accuracy."
        ),
        "",
        "## Model/runtime",
        "",
        f"- Model export: `{MODEL_REPO}`",
        "- Original checkpoint: `gabrielpirlo/SpeD_ParakeetRo_110M_TDT-CTC`",
        f"- Decoder used here: **{payload['model']['decoder']}**",
        "- Audio: 16 kHz mono PCM16",
        f"- Model load time: {payload['model']['model_load_s']:.3f} s",
        "",
        "The original checkpoint also has a TDT decoder and stronger CTC beam-search/KenLM "
        "paths. This run uses the compact ONNX CTC-greedy export and should not be presented "
        "as the checkpoint's best attainable accuracy.",
        "",
        "## Aggregate results",
        "",
        "| Voice/service | Type | Clips | WER | CER | Accentless WER | CPU RTFx |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for item in aggregates:
        lines.append(
            "| {label} | {kind} | {clips} | {wer} | {cer} | {awer} | {rtfx} |".format(
                label=item["service_label"],
                kind=item["service_kind"],
                clips=item["recognized_ok"],
                wer=format_pct(item["mean_wer"]),
                cer=format_pct(item["mean_cer"]),
                awer=format_pct(item["mean_wer_accentless"]),
                rtfx=format_number(item["mean_rtfx"], 1),
            )
        )

    lines.extend(["", "## Per-clip transcriptions", ""])
    for item in results:
        lines.extend(
            [
                f"### {item['service_label']} — {item['category']}",
                "",
                f"- **Reference:** {item['reference']}",
                f"- **SpeD raw:** {item['hypothesis_raw'] or '∅'}",
                f"- **WER/CER:** {format_pct(item['wer'])} / {format_pct(item['cer'])}",
                (
                    f"- **Accentless WER/CER:** {format_pct(item['wer_accentless'])} / "
                    f"{format_pct(item['cer_accentless'])}"
                ),
                (
                    f"- **Latency:** {item['inference_s']:.3f} s for "
                    f"{item['duration_s']:.2f} s audio ({item['rtfx']:.1f}× real time)"
                ),
                f"- **Audio:** `{item['wav_path']}`",
                "",
            ]
        )

    failures = [item for item in generation if item["status"] != "ok"]
    lines.extend(["## Generation failures", ""])
    if failures:
        for item in failures:
            lines.append(
                f"- `{item['service']}/{item['prompt_id']}`: {item.get('error', 'unknown error')}"
            )
    else:
        lines.append("None.")

    lines.extend(
        [
            "",
            "## Metric interpretation",
            "",
            "- WER and CER are computed after lowercasing and removing punctuation.",
            "- Romanian diacritics are preserved in the primary metrics.",
            "- Accentless metrics additionally strip diacritics to isolate lexical recognition.",
            "- The model's raw output is retained so formatting and diacritic behavior remain visible.",
            "",
        ]
    )
    (RESULTS_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def write_html_report(payload: dict[str, Any]) -> None:
    aggregates = payload["summary_by_voice"]
    results = payload["results"]
    generation = payload["generation"]

    aggregate_rows = "\n".join(
        (
            "<tr>"
            f"<td>{html.escape(item['service_label'])}</td>"
            f"<td>{html.escape(item['service_kind'])}</td>"
            f"<td>{item['recognized_ok']}</td>"
            f"<td>{format_pct(item['mean_wer'])}</td>"
            f"<td>{format_pct(item['mean_cer'])}</td>"
            f"<td>{format_pct(item['mean_wer_accentless'])}</td>"
            f"<td>{format_number(item['mean_rtfx'], 1)}</td>"
            "</tr>"
        )
        for item in aggregates
    )

    detail_cards: list[str] = []
    for item in results:
        wav_path = html.escape(item["wav_path"])
        detail_cards.append(
            f"""
            <article class="card">
              <div class="eyebrow">{html.escape(item['service_label'])}</div>
              <h3>{html.escape(item['category'])}</h3>
              <p><strong>Reference</strong><br>{html.escape(item['reference'])}</p>
              <p><strong>SpeD raw</strong><br>{html.escape(item['hypothesis_raw'] or '∅')}</p>
              <div class="metrics">
                <span>WER {format_pct(item['wer'])}</span>
                <span>CER {format_pct(item['cer'])}</span>
                <span>Accentless WER {format_pct(item['wer_accentless'])}</span>
                <span>{item['rtfx']:.1f}× realtime</span>
              </div>
              <audio controls preload="none" src="{wav_path}"></audio>
              <p class="small">{item['duration_s']:.2f}s audio · {item['inference_s']:.3f}s CPU inference</p>
            </article>
            """
        )

    failures = [item for item in generation if item["status"] != "ok"]
    if failures:
        failure_html = "<ul>" + "".join(
            f"<li><code>{html.escape(item['service'])}/{html.escape(item['prompt_id'])}</code>: "
            f"{html.escape(item.get('error') or 'unknown error')}</li>"
            for item in failures
        ) + "</ul>"
    else:
        failure_html = "<p>None.</p>"

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SpeD Romanian TTS benchmark</title>
<style>
  :root {{ color-scheme: light dark; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
  body {{ max-width: 1180px; margin: 0 auto; padding: 32px 20px 80px; line-height: 1.5; }}
  h1 {{ margin-bottom: 4px; }}
  .subtitle, .small {{ opacity: .72; }}
  .warning {{ border-left: 4px solid #c98000; padding: 12px 16px; background: color-mix(in srgb, #c98000 10%, transparent); }}
  table {{ width: 100%; border-collapse: collapse; margin: 18px 0 36px; }}
  th, td {{ border-bottom: 1px solid color-mix(in srgb, currentColor 20%, transparent); padding: 10px 8px; text-align: left; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(330px, 1fr)); gap: 16px; }}
  .card {{ border: 1px solid color-mix(in srgb, currentColor 18%, transparent); border-radius: 14px; padding: 18px; }}
  .card h3 {{ margin-top: 4px; }}
  .eyebrow {{ font-size: .78rem; font-weight: 700; text-transform: uppercase; letter-spacing: .05em; opacity: .68; }}
  .metrics {{ display: flex; flex-wrap: wrap; gap: 8px; margin: 12px 0; }}
  .metrics span {{ padding: 4px 8px; border-radius: 999px; background: color-mix(in srgb, currentColor 9%, transparent); font-size: .82rem; }}
  audio {{ width: 100%; margin-top: 8px; }}
  code {{ overflow-wrap: anywhere; }}
</style>
</head>
<body>
<h1>SpeD Romanian Parakeet 110M</h1>
<p class="subtitle">Romanian TTS → CTC-greedy ONNX ASR benchmark · {html.escape(payload['generated_at'])}</p>

<div class="warning">
<strong>Controlled test only.</strong> Synthetic speech is cleaner and more regular than real microphone dictation.
These numbers identify pronunciation and vocabulary weaknesses; they do not establish real-world production accuracy.
</div>

<h2>Method</h2>
<ul>
  <li>Exact export: <code>{MODEL_REPO}</code></li>
  <li>Decoder: <strong>{html.escape(payload['model']['decoder'])}</strong></li>
  <li>Input normalized to 16 kHz mono PCM16</li>
  <li>Primary WER/CER preserve Romanian diacritics; accentless metrics strip them</li>
  <li>Punctuation and case are excluded from metric scoring; raw output is shown</li>
</ul>

<h2>Aggregate results</h2>
<table>
<thead><tr><th>Voice/service</th><th>Type</th><th>Clips</th><th>WER</th><th>CER</th><th>Accentless WER</th><th>CPU RTFx</th></tr></thead>
<tbody>{aggregate_rows}</tbody>
</table>

<h2>Listen and inspect</h2>
<div class="grid">{''.join(detail_cards)}</div>

<h2>Generation failures</h2>
{failure_html}

<h2>Interpretation caveats</h2>
<p>
The original SpeD checkpoint includes both TDT and CTC heads, and its published best path uses stronger decoding than this
portable CTC-greedy ONNX run. TTS voices can also align unusually well—or poorly—with a model's training distribution.
A production decision requires recordings from multiple real Romanian speakers, microphones, rooms, speaking rates,
accents, whispered speech, disfluencies, and Romanian-English technical dictation.
</p>
</body>
</html>
"""
    (RESULTS_DIR / "report.html").write_text(document, encoding="utf-8")


def write_failure_bundle(exc: BaseException) -> None:
    ensure_directories()
    failure = {
        "status": "failed",
        "error_type": type(exc).__name__,
        "error": str(exc),
        "traceback": traceback.format_exc(),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "environment": collect_environment(),
    }
    (RESULTS_DIR / "fatal_error.json").write_text(
        json.dumps(failure, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


async def main() -> int:
    ensure_directories()
    LOG_PATH.write_text("", encoding="utf-8")
    log("Starting Romanian TTS -> SpeD benchmark")
    generation = await generate_all_audio()
    success_count = sum(item.status == "ok" for item in generation)
    failure_count = sum(item.status != "ok" for item in generation)
    log(f"Audio generation finished: {success_count} succeeded, {failure_count} failed")

    records, model_load_s, model_metadata = recognize_all(generation)
    log(f"Recognition completed for {len(records)} clips; model load {model_load_s:.3f}s")
    write_outputs(records, generation, model_metadata)
    log("Reports written to results/")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(asyncio.run(main()))
    except SystemExit:
        raise
    except BaseException as exc:
        print(traceback.format_exc(), file=sys.stderr, flush=True)
        write_failure_bundle(exc)
        raise
