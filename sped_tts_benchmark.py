#!/usr/bin/env python3
"""Benchmark the actual SpeD Romanian Parakeet 110M ONNX/CTC checkpoint.

Creates Romanian speech with several TTS engines, normalizes all clips to
16 kHz mono PCM WAV, transcribes them with SpeD, and exports audio plus
per-clip/aggregate WER, CER, and latency results.
"""

from __future__ import annotations

import asyncio
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import unicodedata
import wave
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Callable

import edge_tts
import onnx_asr
import soundfile as sf
from gtts import gTTS
from huggingface_hub import snapshot_download
from jiwer import cer, wer
from piper import PiperVoice

ROOT = Path("benchmark_artifact")
AUDIO_RAW = ROOT / "audio_raw"
AUDIO_WAV = ROOT / "audio_wav_16k"
MODEL_DIR = ROOT / "model"
PIPER_DIR = ROOT / "piper"

MODEL_REPO = "AlinClaudiu/SpeD-ParakeetRo-110M-onnx"
MODEL_FILES = ["model.onnx", "model.onnx.data", "vocab.txt", "config.json"]
EXPECTED_SHA256 = {
    "model.onnx": "c5ffc434d9f16d9ba257f380b3634b7ca19ef9b5b7420489705669539e2a640d",
    "model.onnx.data": "81e884781377a1365802777295dec05c068e4005be9f657a5519994a5071544b",
}
PIPER_REPO = "rhasspy/piper-voices"
PIPER_MODEL_PATH = "ro/ro_RO/mihai/medium/ro_RO-mihai-medium.onnx"
PIPER_CONFIG_PATH = PIPER_MODEL_PATH + ".json"

UTTERANCES: list[dict[str, str]] = [
    {
        "id": "clean",
        "category": "clean_ro",
        "text": "Bună ziua. Acesta este un test de recunoaștere automată a vorbirii în limba română.",
    },
    {
        "id": "diacritics",
        "category": "diacritics",
        "text": "Ștefan și Țăndărică au găsit o căsuță lângă râul înghețat, într-o dimineață liniștită.",
    },
    {
        "id": "date_time",
        "category": "numbers_dates",
        "text": "Întâlnirea este programată pe douăzeci și trei august, la ora cincisprezece și treizeci.",
    },
    {
        "id": "technical",
        "category": "technical_ro_en",
        "text": "Serverul rulează pe portul opt mii optzeci și folosește o bază de date PostgreSQL.",
    },
    {
        "id": "mixed_devops",
        "category": "technical_ro_en",
        "text": "Trebuie să repornim deployment-ul PhotonSpark și să verificăm logurile din Kubernetes.",
    },
    {
        "id": "networking",
        "category": "technical_ro_en",
        "text": "Configurează adresa IPv6, certificatul Cloudflare și tunelul WireGuard pe nodul Proxmox.",
    },
    {
        "id": "self_correction",
        "category": "self_correction",
        "text": "Trimite raportul joi, de fapt vineri dimineață, după verificarea copiilor de siguranță.",
    },
    {
        "id": "short_command",
        "category": "short_command",
        "text": "Pornește serverul și verifică starea serviciului.",
    },
    {
        "id": "product_names",
        "category": "proper_names",
        "text": "Cesar testează SearXNG, Qdrant, FastAPI și Pterodactyl pe stația cu placă NVIDIA.",
    },
    {
        "id": "long_dictation",
        "category": "long_dictation",
        "text": "Înainte să publicăm noua versiune, compară logurile de pe toate nodurile, verifică latența bazei de date și creează o copie de siguranță completă, astfel încât serviciul să poată fi restaurat rapid dacă apare o eroare.",
    },
]

PROVIDERS: list[dict[str, str]] = [
    {
        "id": "edge_alina",
        "kind": "edge",
        "label": "Microsoft Edge Alina Neural",
        "voice": "ro-RO-AlinaNeural",
        "rate": "+0%",
    },
    {
        "id": "edge_emil",
        "kind": "edge",
        "label": "Microsoft Edge Emil Neural",
        "voice": "ro-RO-EmilNeural",
        "rate": "+0%",
    },
    {
        "id": "edge_alina_fast",
        "kind": "edge",
        "label": "Microsoft Edge Alina Neural (+25% rate)",
        "voice": "ro-RO-AlinaNeural",
        "rate": "+25%",
    },
    {
        "id": "google_gtts",
        "kind": "gtts",
        "label": "Google Romanian gTTS",
    },
    {
        "id": "piper_mihai",
        "kind": "piper",
        "label": "Piper Mihai medium",
    },
    {
        "id": "espeak_ro",
        "kind": "espeak",
        "label": "eSpeak NG Romanian baseline",
    },
]


def run(cmd: list[str], *, input_text: str | None = None) -> None:
    subprocess.run(
        cmd,
        input=None if input_text is None else input_text.encode("utf-8"),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def retry(fn: Callable[[], None], attempts: int = 4, base_delay: float = 2.0) -> None:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            fn()
            return
        except Exception as exc:  # network TTS can fail transiently
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(base_delay * (2**attempt))
    assert last_error is not None
    raise last_error


def ffmpeg_to_wav(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    run(
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
            str(target),
        ]
    )


async def edge_save(text: str, path: Path, voice: str, rate: str) -> None:
    communicate = edge_tts.Communicate(text=text, voice=voice, rate=rate)
    await communicate.save(str(path))


def synthesize_edge(text: str, target: Path, provider: dict[str, str]) -> None:
    retry(lambda: asyncio.run(edge_save(text, target, provider["voice"], provider["rate"])))


def synthesize_gtts(text: str, target: Path) -> None:
    retry(lambda: gTTS(text=text, lang="ro", slow=False, tld="com").save(str(target)))


def synthesize_espeak(text: str, target: Path) -> None:
    run(["espeak-ng", "-v", "ro", "-s", "165", "-w", str(target), text])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_text(text: str) -> str:
    text = unicodedata.normalize("NFKC", text)
    text = text.translate(str.maketrans({"ş": "ș", "Ş": "Ș", "ţ": "ț", "Ţ": "Ț"}))
    text = text.lower().replace("-", " ")
    text = re.sub(r"[^0-9a-zăâîșț\s]", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def audio_duration(path: Path) -> float:
    info = sf.info(path)
    return float(info.frames / info.samplerate)


def safe_metric(metric: Callable[[str, str], float], reference: str, hypothesis: str) -> float:
    if not reference:
        return 0.0 if not hypothesis else 1.0
    return float(metric(reference, hypothesis))


def markdown_escape(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def aggregate_rows(rows: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["status"] == "ok":
            groups[str(row[key])].append(row)

    output: list[dict[str, Any]] = []
    for name, group in sorted(groups.items()):
        ref = " ".join(str(item["reference_normalized"]) for item in group)
        hyp = " ".join(str(item["hypothesis_normalized"]) for item in group)
        latencies = [float(item["inference_seconds"]) for item in group]
        rtfx_values = [float(item["rtfx"]) for item in group]
        output.append(
            {
                key: name,
                "clips": len(group),
                "duration_seconds": round(sum(float(item["duration_seconds"]) for item in group), 3),
                "wer": round(safe_metric(wer, ref, hyp), 6),
                "cer": round(safe_metric(cer, ref, hyp), 6),
                "median_inference_seconds": round(median(latencies), 6),
                "median_rtfx": round(median(rtfx_values), 3),
            }
        )
    return output


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    rows: list[dict[str, Any]],
    provider_summary: list[dict[str, Any]],
    category_summary: list[dict[str, Any]],
    model_manifest: dict[str, Any],
) -> None:
    lines = [
        "# SpeD Romanian Parakeet 110M — synthetic TTS benchmark",
        "",
        "This run transcribes generated Romanian speech with the actual `AlinClaudiu/SpeD-ParakeetRo-110M-onnx` CTC export. It uses greedy CTC decoding through `onnx-asr`, not the original checkpoint's TDT decoder or the published 6-gram language-model beam search.",
        "",
        "## Aggregate by TTS profile",
        "",
        "| Profile | Clips | WER | CER | Median inference | Median RTFx |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    labels = {provider["id"]: provider["label"] for provider in PROVIDERS}
    for item in provider_summary:
        lines.append(
            f"| {markdown_escape(labels.get(item['provider'], item['provider']))} | {item['clips']} | "
            f"{100 * item['wer']:.2f}% | {100 * item['cer']:.2f}% | "
            f"{item['median_inference_seconds']:.3f} s | {item['median_rtfx']:.1f}× |"
        )

    lines += [
        "",
        "## Aggregate by sentence category",
        "",
        "| Category | Clips | WER | CER | Median inference | Median RTFx |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for item in category_summary:
        lines.append(
            f"| {markdown_escape(item['category'])} | {item['clips']} | "
            f"{100 * item['wer']:.2f}% | {100 * item['cer']:.2f}% | "
            f"{item['median_inference_seconds']:.3f} s | {item['median_rtfx']:.1f}× |"
        )

    ok_rows = [row for row in rows if row["status"] == "ok"]
    worst = sorted(ok_rows, key=lambda row: (float(row["wer"]), float(row["cer"])), reverse=True)[:15]
    lines += [
        "",
        "## Fifteen highest-error clips",
        "",
        "| TTS profile | Case | WER | Reference | SpeD output |",
        "|---|---|---:|---|---|",
    ]
    for row in worst:
        lines.append(
            f"| {markdown_escape(labels.get(row['provider'], row['provider']))} | {row['utterance_id']} | "
            f"{100 * float(row['wer']):.1f}% | {markdown_escape(row['reference'])} | "
            f"{markdown_escape(row['hypothesis'])} |"
        )

    lines += [
        "",
        "## Complete clip-level results",
        "",
        "| TTS profile | Case | Duration | Latency | WER | CER | SpeD output |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        if row["status"] == "ok":
            lines.append(
                f"| {markdown_escape(labels.get(row['provider'], row['provider']))} | {row['utterance_id']} | "
                f"{float(row['duration_seconds']):.2f} s | {float(row['inference_seconds']):.3f} s | "
                f"{100 * float(row['wer']):.1f}% | {100 * float(row['cer']):.1f}% | "
                f"{markdown_escape(row['hypothesis'])} |"
            )
        else:
            lines.append(
                f"| {markdown_escape(labels.get(row['provider'], row['provider']))} | {row['utterance_id']} | "
                f"— | — | — | — | ERROR: {markdown_escape(str(row['error']))} |"
            )

    lines += [
        "",
        "## Model integrity",
        "",
        "```json",
        json.dumps(model_manifest, ensure_ascii=False, indent=2),
        "```",
        "",
        "## Interpretation limits",
        "",
        "Synthetic TTS is clean, deterministic, close-mic audio. It is useful for catching pronunciation, vocabulary, decoder, and speed failures, but it is not a substitute for Romanian microphone recordings with room acoustics, accents, hesitations, clipped word starts, and background noise.",
    ]
    (ROOT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    for directory in (ROOT, AUDIO_RAW, AUDIO_WAV, MODEL_DIR, PIPER_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
    }
    (ROOT / "environment.json").write_text(json.dumps(environment, indent=2), encoding="utf-8")

    print(f"Downloading exact SpeD ONNX model from {MODEL_REPO}", flush=True)
    snapshot_download(repo_id=MODEL_REPO, local_dir=MODEL_DIR, allow_patterns=MODEL_FILES)

    model_manifest: dict[str, Any] = {"repo": MODEL_REPO, "files": {}}
    for filename in MODEL_FILES:
        path = MODEL_DIR / filename
        if not path.exists():
            raise FileNotFoundError(f"Required model file missing: {path}")
        digest = sha256(path)
        expected = EXPECTED_SHA256.get(filename)
        model_manifest["files"][filename] = {
            "bytes": path.stat().st_size,
            "sha256": digest,
            "expected_sha256": expected,
            "matches_expected": None if expected is None else digest == expected,
        }

    print("Downloading Piper Romanian Mihai voice", flush=True)
    piper_snapshot = Path(
        snapshot_download(
            repo_id=PIPER_REPO,
            local_dir=PIPER_DIR,
            allow_patterns=[PIPER_MODEL_PATH, PIPER_CONFIG_PATH],
        )
    )
    piper_model = piper_snapshot / PIPER_MODEL_PATH
    piper_config = piper_snapshot / PIPER_CONFIG_PATH
    piper_voice: PiperVoice | None = None
    try:
        piper_voice = PiperVoice.load(str(piper_model), config_path=str(piper_config))
    except Exception as exc:
        print(f"Piper initialization failed; provider will be recorded as failed: {exc}", flush=True)

    synthesis_manifest: list[dict[str, Any]] = []
    for provider in PROVIDERS:
        print(f"Synthesizing provider: {provider['label']}", flush=True)
        for utterance in UTTERANCES:
            stem = f"{provider['id']}__{utterance['id']}"
            raw_suffix = ".mp3" if provider["kind"] in {"edge", "gtts"} else ".wav"
            raw_path = AUDIO_RAW / f"{stem}{raw_suffix}"
            wav_path = AUDIO_WAV / f"{stem}.wav"
            item: dict[str, Any] = {
                "provider": provider["id"],
                "provider_label": provider["label"],
                "utterance_id": utterance["id"],
                "category": utterance["category"],
                "reference": utterance["text"],
                "raw_path": str(raw_path),
                "wav_path": str(wav_path),
            }
            try:
                if provider["kind"] == "edge":
                    synthesize_edge(utterance["text"], raw_path, provider)
                elif provider["kind"] == "gtts":
                    synthesize_gtts(utterance["text"], raw_path)
                elif provider["kind"] == "espeak":
                    synthesize_espeak(utterance["text"], raw_path)
                elif provider["kind"] == "piper":
                    if piper_voice is None:
                        raise RuntimeError("Piper voice could not be loaded")
                    with wave.open(str(raw_path), "wb") as wav_file:
                        piper_voice.synthesize_wav(utterance["text"], wav_file)
                else:
                    raise ValueError(f"Unknown provider kind: {provider['kind']}")

                ffmpeg_to_wav(raw_path, wav_path)
                item["status"] = "ok"
                item["duration_seconds"] = round(audio_duration(wav_path), 6)
            except Exception as exc:
                item["status"] = "error"
                item["error"] = f"{type(exc).__name__}: {exc}"
                print(f"Synthesis error for {stem}: {item['error']}", flush=True)
            synthesis_manifest.append(item)

    successful = [item for item in synthesis_manifest if item["status"] == "ok"]
    if not successful:
        raise RuntimeError("No TTS clips were generated")

    print("Loading SpeD NeMo FastConformer CTC model", flush=True)
    model = onnx_asr.load_model(
        "nemo-conformer-ctc",
        MODEL_DIR,
        providers=["CPUExecutionProvider"],
        preprocessor_config={"max_concurrent_workers": 1, "use_numpy_preprocessors": True},
    )

    print("Warming up inference", flush=True)
    _ = model.recognize(successful[0]["wav_path"])

    rows: list[dict[str, Any]] = []
    synthesis_lookup = {(item["provider"], item["utterance_id"]): item for item in synthesis_manifest}
    for provider in PROVIDERS:
        for utterance in UTTERANCES:
            synth = synthesis_lookup[(provider["id"], utterance["id"])]
            base: dict[str, Any] = {
                "provider": provider["id"],
                "provider_label": provider["label"],
                "utterance_id": utterance["id"],
                "category": utterance["category"],
                "reference": utterance["text"],
                "reference_normalized": normalize_text(utterance["text"]),
                "audio_path": synth["wav_path"],
            }
            if synth["status"] != "ok":
                base.update({"status": "error", "error": synth["error"]})
                rows.append(base)
                continue

            try:
                print(f"Transcribing {provider['id']} / {utterance['id']}", flush=True)
                started = time.perf_counter()
                hypothesis = str(model.recognize(synth["wav_path"]))
                elapsed = time.perf_counter() - started
                duration = float(synth["duration_seconds"])
                ref_norm = base["reference_normalized"]
                hyp_norm = normalize_text(hypothesis)
                base.update(
                    {
                        "status": "ok",
                        "hypothesis": hypothesis,
                        "hypothesis_normalized": hyp_norm,
                        "duration_seconds": round(duration, 6),
                        "inference_seconds": round(elapsed, 6),
                        "rtf": round(elapsed / duration, 6),
                        "rtfx": round(duration / elapsed, 3),
                        "wer": round(safe_metric(wer, ref_norm, hyp_norm), 6),
                        "cer": round(safe_metric(cer, ref_norm, hyp_norm), 6),
                    }
                )
            except Exception as exc:
                base.update({"status": "error", "error": f"{type(exc).__name__}: {exc}"})
                print(f"Inference error: {base['error']}", flush=True)
            rows.append(base)

    provider_summary = aggregate_rows(rows, "provider")
    category_summary = aggregate_rows(rows, "category")

    write_csv(ROOT / "results.csv", rows)
    write_csv(ROOT / "summary_by_provider.csv", provider_summary)
    write_csv(ROOT / "summary_by_category.csv", category_summary)
    (ROOT / "results.json").write_text(
        json.dumps(
            {
                "environment": environment,
                "model": model_manifest,
                "utterances": UTTERANCES,
                "providers": PROVIDERS,
                "synthesis": synthesis_manifest,
                "results": rows,
                "summary_by_provider": provider_summary,
                "summary_by_category": category_summary,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    write_report(rows, provider_summary, category_summary, model_manifest)

    # Model weights are reproducibly downloadable and would inflate the user artifact.
    shutil.rmtree(MODEL_DIR, ignore_errors=True)
    shutil.rmtree(PIPER_DIR, ignore_errors=True)

    print("\nProvider summary:", flush=True)
    print(json.dumps(provider_summary, indent=2, ensure_ascii=False), flush=True)
    print(f"\nArtifact prepared at {ROOT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
