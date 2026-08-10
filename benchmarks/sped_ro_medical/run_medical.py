#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import csv
import html
import importlib.util
import json
import statistics
import sys
import time
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASE_PATH = ROOT / "benchmarks" / "sped_ro_tts" / "benchmark.py"
CORPUS_PATH = HERE / "medical_prompts.json"
RESULTS = ROOT / "medical_results"

spec = importlib.util.spec_from_file_location("sped_base", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot import {BASE_PATH}")
base = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = base
spec.loader.exec_module(base)

PRIMARY = {"edge_alina", "edge_emil", "google_gtts", "piper_mihai"}
base.SERVICE_METADATA = {
    key: value for key, value in base.SERVICE_METADATA.items()
    if key in PRIMARY or key == "espeak_ng"
}
base.RESULTS_DIR = RESULTS
base.AUDIO_DIR = RESULTS / "audio"
base.SOURCE_AUDIO_DIR = base.AUDIO_DIR / "source"
base.WAV_AUDIO_DIR = base.AUDIO_DIR / "wav_16k_mono"
base.LOG_PATH = RESULTS / "benchmark.log"


def load_prompts() -> list[dict[str, Any]]:
    prompts = json.loads(CORPUS_PATH.read_text(encoding="utf-8"))
    seen: set[str] = set()
    for p in prompts:
        if p["id"] in seen:
            raise ValueError(f"duplicate id: {p['id']}")
        seen.add(p["id"])
        text = base.normalize_text(p["text"])
        for term in p["medical_terms"]:
            if base.normalize_text(term) not in text:
                raise ValueError(f"term {term!r} not found in prompt {p['id']}")
    return prompts


async def generate_audio(prompts: list[dict[str, Any]]) -> list[Any]:
    records: list[Any] = []
    for service, voice in (("edge_alina", "ro-RO-AlinaNeural"), ("edge_emil", "ro-RO-EmilNeural")):
        sem = asyncio.Semaphore(3)
        async def one(p: dict[str, Any]) -> Any:
            src = base.SOURCE_AUDIO_DIR / service / f"{p['id']}.mp3"
            wav = base.WAV_AUDIO_DIR / service / f"{p['id']}.wav"
            try:
                async with sem:
                    await base.generate_edge_voice(voice, p["text"], src)
                base.convert_to_wav(src, wav)
                return base.GenerationRecord(service, p["id"], "ok", str(src.relative_to(RESULTS)), str(wav.relative_to(RESULTS)))
            except BaseException as exc:
                return base.GenerationRecord(service, p["id"], "failed", error=f"{type(exc).__name__}: {exc}")
        records.extend(await asyncio.gather(*(one(p) for p in prompts)))

    for p in prompts:
        service = "google_gtts"
        src = base.SOURCE_AUDIO_DIR / service / f"{p['id']}.mp3"
        wav = base.WAV_AUDIO_DIR / service / f"{p['id']}.wav"
        try:
            base.generate_google_tts(p["text"], src)
            base.convert_to_wav(src, wav)
            records.append(base.GenerationRecord(service, p["id"], "ok", str(src.relative_to(RESULTS)), str(wav.relative_to(RESULTS))))
        except BaseException as exc:
            records.append(base.GenerationRecord(service, p["id"], "failed", error=f"{type(exc).__name__}: {exc}"))
        time.sleep(0.2)

    for service, generator in (("piper_mihai", base.generate_piper), ("espeak_ng", base.generate_espeak)):
        for p in prompts:
            ext = base.source_extension(service)
            src = base.SOURCE_AUDIO_DIR / service / f"{p['id']}{ext}"
            wav = base.WAV_AUDIO_DIR / service / f"{p['id']}.wav"
            try:
                generator(p["text"], src)
                base.convert_to_wav(src, wav)
                records.append(base.GenerationRecord(service, p["id"], "ok", str(src.relative_to(RESULTS)), str(wav.relative_to(RESULTS))))
            except BaseException as exc:
                records.append(base.GenerationRecord(service, p["id"], "failed", error=f"{type(exc).__name__}: {exc}"))
                if service == "piper_mihai":
                    break
    (RESULTS / "generation.json").write_text(json.dumps([asdict(x) for x in records], ensure_ascii=False, indent=2), encoding="utf-8")
    return records


def contains(tokens: list[str], phrase: list[str]) -> bool:
    return any(tokens[i:i + len(phrase)] == phrase for i in range(max(0, len(tokens) - len(phrase) + 1))) if phrase else True


def best_span(ref: list[str], hyp: list[str]) -> tuple[int, str]:
    if not hyp:
        return len(ref), ""
    best = (10**9, "")
    for size in range(max(1, len(ref) - 2), min(len(hyp), len(ref) + 3) + 1):
        for i in range(len(hyp) - size + 1):
            span = hyp[i:i + size]
            d = base.levenshtein_distance(ref, span)
            if d < best[0]:
                best = (d, " ".join(span))
    return best


def clip_counts(record: Any) -> dict[str, int]:
    rw = record.reference_normalized.split(); hw = record.hypothesis_normalized.split()
    rc = list(record.reference_normalized.replace(" ", "")); hc = list(record.hypothesis_normalized.replace(" ", ""))
    return {"we": base.levenshtein_distance(rw, hw), "rw": len(rw), "ce": base.levenshtein_distance(rc, hc), "rc": len(rc)}


def aggregate(rows: list[dict[str, Any]], fields: tuple[str, ...]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[tuple(row[f] for f in fields)].append(row)
    out = []
    for key, items in sorted(groups.items(), key=lambda x: tuple(str(v) for v in x[0])):
        r = {f: v for f, v in zip(fields, key)}
        r.update({
            "clips": len(items), "audio_seconds": sum(x["duration_s"] for x in items),
            "reference_words": sum(x["rw"] for x in items), "word_errors": sum(x["we"] for x in items),
            "reference_chars": sum(x["rc"] for x in items), "char_errors": sum(x["ce"] for x in items),
            "term_words": sum(x["term_words"] for x in items), "term_errors": sum(x["term_errors"] for x in items),
            "terms": sum(x["terms"] for x in items), "exact_hits": sum(x["exact_hits"] for x in items),
            "accentless_hits": sum(x["accentless_hits"] for x in items), "fuzzy_hits": sum(x["fuzzy_hits"] for x in items),
            "mean_rtfx": statistics.fmean(x["rtfx"] for x in items),
        })
        r["micro_wer"] = r["word_errors"] / r["reference_words"]
        r["micro_cer"] = r["char_errors"] / r["reference_chars"]
        r["mter"] = r["term_errors"] / r["term_words"]
        r["exact_term_recall"] = r["exact_hits"] / r["terms"]
        r["accentless_exact_recall"] = r["accentless_hits"] / r["terms"]
        r["fuzzy_term_recall"] = r["fuzzy_hits"] / r["terms"]
        out.append(r)
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8"); return
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def medical_outputs(prompts: list[dict[str, Any]], records: list[Any]) -> dict[str, Any]:
    prompt_map = {p["id"]: p for p in prompts}
    clips: list[dict[str, Any]] = []
    terms: list[dict[str, Any]] = []
    for rec in records:
        if rec.status != "ok": continue
        p = prompt_map[rec.prompt_id]; counts = clip_counts(rec)
        hyp = rec.hypothesis_normalized.split(); hyp_a = rec.hypothesis_accentless.split()
        te = tw = exact = exact_a = fuzzy = 0
        for term in p["medical_terms"]:
            t = base.normalize_text(term).split(); ta = base.normalize_text(term, strip_diacritics=True).split()
            d, span = best_span(t, hyp); er = d / len(t)
            hit = contains(hyp, t); hit_a = contains(hyp_a, ta)
            te += d; tw += len(t); exact += hit; exact_a += hit_a; fuzzy += er <= .25
            terms.append({"service": rec.service, "prompt_id": rec.prompt_id, "suite": p["suite"], "category": p["category"], "term": term, "exact_hit": hit, "accentless_exact_hit": hit_a, "best_span": span, "term_errors": d, "term_words": len(t), "term_error_rate": er})
        clips.append({"service": rec.service, "service_label": rec.service_label, "service_kind": rec.service_kind, "primary_neural": rec.service in PRIMARY, "prompt_id": rec.prompt_id, "suite": p["suite"], "category": p["category"], "reference": rec.reference, "hypothesis": rec.hypothesis_raw, "duration_s": rec.duration_s, "inference_s": rec.inference_s, "rtfx": rec.rtfx, **counts, "term_errors": te, "term_words": tw, "terms": len(p["medical_terms"]), "exact_hits": exact, "accentless_hits": exact_a, "fuzzy_hits": fuzzy})

    primary = [x for x in clips if x["primary_neural"]]
    by_voice = aggregate(clips, ("service", "service_label", "service_kind", "primary_neural"))
    by_suite = aggregate(primary, ("suite",)); by_category = aggregate(primary, ("category",))
    overall = aggregate(primary, ("primary_neural",))[0]; overall.pop("primary_neural", None)

    tg: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for x in terms:
        if x["service"] in PRIMARY: tg[x["term"]].append(x)
    term_summary = []
    for term, items in tg.items():
        term_summary.append({"term": term, "instances": len(items), "exact_recall": sum(x["exact_hit"] for x in items) / len(items), "accentless_exact_recall": sum(x["accentless_exact_hit"] for x in items) / len(items), "mter": sum(x["term_errors"] for x in items) / sum(x["term_words"] for x in items), "recognized_spans": " | ".join(dict.fromkeys(x["best_span"] for x in items if x["best_span"]))[:500]})
    term_summary.sort(key=lambda x: (x["exact_recall"], -x["mter"], x["term"].lower()))

    write_csv(RESULTS / "medical_clip_metrics.csv", clips)
    write_csv(RESULTS / "medical_term_details.csv", terms)
    write_csv(RESULTS / "medical_summary_by_voice.csv", by_voice)
    write_csv(RESULTS / "medical_summary_by_suite.csv", by_suite)
    write_csv(RESULTS / "medical_summary_by_category.csv", by_category)
    write_csv(RESULTS / "medical_term_summary.csv", term_summary)
    payload = {"generated_at": base.datetime.now(base.timezone.utc).isoformat(), "corpus": {"prompts": len(prompts), "annotated_terms_per_voice": sum(len(p["medical_terms"]) for p in prompts), "categories": len({p["category"] for p in prompts})}, "overall_primary_neural": overall, "by_voice": by_voice, "by_suite": by_suite, "by_category": by_category, "term_summary": term_summary}
    (RESULTS / "medical_metrics.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload


def pct(x: float) -> str: return f"{100*x:.2f}%"


def write_medical_report(data: dict[str, Any]) -> None:
    o = data["overall_primary_neural"]
    lines = ["# SpeD Romanian 110M — medical TTS benchmark", "", "Controlled synthetic-speech benchmark; not clinical validation and not a substitute for real clinician microphone recordings.", "", "## Headline result — four neural voices", "", f"- Prompts per voice: **{data['corpus']['prompts']}**", f"- Annotated medical terms per voice: **{data['corpus']['annotated_terms_per_voice']}**", f"- Clips: **{o['clips']}**", f"- Audio: **{o['audio_seconds']/60:.1f} minutes**", f"- Micro-WER: **{pct(o['micro_wer'])}**", f"- Micro-CER: **{pct(o['micro_cer'])}**", f"- Medical-term error rate: **{pct(o['mter'])}**", f"- Exact medical-term recall: **{pct(o['exact_term_recall'])}**", f"- Accentless exact recall: **{pct(o['accentless_exact_recall'])}**", f"- Fuzzy term recall, TER at most 25%: **{pct(o['fuzzy_term_recall'])}**", "", "## By voice", "", "| Voice | Clips | WER | CER | MTER | Exact terms | Fuzzy terms | CPU RTFx |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for x in data["by_voice"]: lines.append(f"| {x['service_label']} | {x['clips']} | {pct(x['micro_wer'])} | {pct(x['micro_cer'])} | {pct(x['mter'])} | {pct(x['exact_term_recall'])} | {pct(x['fuzzy_term_recall'])} | {x['mean_rtfx']:.1f} |")
    lines += ["", "## By suite — neural voices", "", "| Suite | Clips | WER | CER | MTER | Exact terms |", "|---|---:|---:|---:|---:|---:|"]
    for x in data["by_suite"]: lines.append(f"| {x['suite']} | {x['clips']} | {pct(x['micro_wer'])} | {pct(x['micro_cer'])} | {pct(x['mter'])} | {pct(x['exact_term_recall'])} |")
    lines += ["", "## Hardest categories — neural voices", "", "| Category | Clips | WER | CER | MTER | Exact terms |", "|---|---:|---:|---:|---:|---:|"]
    for x in sorted(data["by_category"], key=lambda z: z["micro_wer"], reverse=True): lines.append(f"| {x['category']} | {x['clips']} | {pct(x['micro_wer'])} | {pct(x['micro_cer'])} | {pct(x['mter'])} | {pct(x['exact_term_recall'])} |")
    lines += ["", "## Hardest medical terms — neural voices", "", "| Term | Exact recall | MTER | Recognized spans |", "|---|---:|---:|---|"]
    for x in data["term_summary"][:50]: lines.append(f"| {x['term']} | {pct(x['exact_recall'])} | {pct(x['mter'])} | {x['recognized_spans'][:180]} |")
    lines += ["", "## Metrics", "", "WER and CER ignore case and punctuation. MTER is the micro-averaged best-window word edit rate over annotated medical terms. Exact recall requires a contiguous exact normalized term; fuzzy recall allows at most twenty-five percent term error.", "", "The tested model is the deployable ONNX CTC-greedy export. The original hybrid SpeD checkpoint supports stronger decoding paths."]
    (RESULTS / "summary.md").write_text("\n".join(lines), encoding="utf-8")

    def tr(rows: list[dict[str, Any]], key: str) -> str:
        return "".join(f"<tr><td>{html.escape(str(x[key]))}</td><td>{x['clips']}</td><td>{pct(x['micro_wer'])}</td><td>{pct(x['micro_cer'])}</td><td>{pct(x['mter'])}</td><td>{pct(x['exact_term_recall'])}</td></tr>" for x in rows)
    terms = "".join(f"<tr><td>{html.escape(x['term'])}</td><td>{pct(x['exact_recall'])}</td><td>{pct(x['mter'])}</td><td>{html.escape(x['recognized_spans'][:240])}</td></tr>" for x in data["term_summary"])
    doc = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>SpeD Romanian medical benchmark</title><style>body{{max-width:1200px;margin:auto;padding:28px;font-family:system-ui;line-height:1.45}}.k{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px}}.k div{{border:1px solid #8885;border-radius:12px;padding:12px}}.k b{{display:block;font-size:1.35rem}}table{{width:100%;border-collapse:collapse;margin:12px 0 28px}}th,td{{border-bottom:1px solid #8885;padding:8px;text-align:left}}.w{{border-left:4px solid #c80;padding:12px;background:#c801}}a{{font-weight:700}}</style></head><body><h1>SpeD Romanian 110M — medical TTS benchmark</h1><p class="w"><b>Controlled synthetic test.</b> Not clinical validation; real clinician and microphone recordings are still required.</p><div class="k"><div>Neural clips<b>{o['clips']}</b></div><div>Micro-WER<b>{pct(o['micro_wer'])}</b></div><div>Micro-CER<b>{pct(o['micro_cer'])}</b></div><div>MTER<b>{pct(o['mter'])}</b></div><div>Exact terms<b>{pct(o['exact_term_recall'])}</b></div><div>Fuzzy terms<b>{pct(o['fuzzy_term_recall'])}</b></div></div><p><a href="all_clips_report.html">Open the complete report with every audio clip and transcription</a></p><h2>By voice</h2><table><tr><th>Voice</th><th>Clips</th><th>WER</th><th>CER</th><th>MTER</th><th>Exact terms</th></tr>{tr(data['by_voice'],'service_label')}</table><h2>By suite — neural voices</h2><table><tr><th>Suite</th><th>Clips</th><th>WER</th><th>CER</th><th>MTER</th><th>Exact terms</th></tr>{tr(data['by_suite'],'suite')}</table><h2>By category — neural voices</h2><table><tr><th>Category</th><th>Clips</th><th>WER</th><th>CER</th><th>MTER</th><th>Exact terms</th></tr>{tr(sorted(data['by_category'],key=lambda z:z['micro_wer'],reverse=True),'category')}</table><h2>Medical terms — neural voices</h2><table><tr><th>Term</th><th>Exact recall</th><th>MTER</th><th>Recognized spans</th></tr>{terms}</table></body></html>'''
    (RESULTS / "report.html").write_text(doc, encoding="utf-8")


async def main() -> int:
    prompts = load_prompts(); base.PROMPTS = prompts; base.ensure_directories(); base.LOG_PATH.write_text("", encoding="utf-8")
    base.log(f"Starting medical benchmark: {len(prompts)} prompts and {sum(len(p['medical_terms']) for p in prompts)} terms per voice")
    generation = await generate_audio(prompts)
    base.log(f"Generation: {sum(x.status=='ok' for x in generation)} succeeded, {sum(x.status!='ok' for x in generation)} failed")
    records, _, metadata = base.recognize_all(generation)
    base.write_outputs(records, generation, metadata)
    (RESULTS / "report.html").replace(RESULTS / "all_clips_report.html")
    (RESULTS / "summary.md").replace(RESULTS / "base_summary.md")
    data = medical_outputs(prompts, records); write_medical_report(data)
    base.log("Medical benchmark complete")
    return 0

if __name__ == "__main__":
    raise SystemExit(asyncio.run(main()))
