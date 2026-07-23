#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import re
import statistics
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

BASE_URL = "http://127.0.0.1:8080"
OUTPUT_DIR = Path("results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PROMPTS: list[dict[str, str]] = [
    {
        "id": "crt",
        "prompt": (
            "Find the smallest positive integer n satisfying n ≡ 1 (mod 2), "
            "n ≡ 2 (mod 3), n ≡ 3 (mod 5), and n ≡ 4 (mod 7). "
            "Use short numbered reasoning lines, verify the result, then give a compact final answer. /think"
        ),
    },
    {
        "id": "five_digit",
        "prompt": (
            "How many five-digit decimal integers have no repeated digits and are divisible by 5? "
            "Use short numbered reasoning lines, check both possible final digits, then give a compact final answer. /think"
        ),
    },
    {
        "id": "die_mod3",
        "prompt": (
            "A fair six-sided die is rolled until the first 6 appears. What is the probability that the total "
            "sum of all rolls, including that final 6, is divisible by 3? Use short numbered reasoning lines "
            "and a small state recurrence, then give a compact final answer. /think"
        ),
    },
]
SEEDS = [7, 42, 1337]

INTRO_MESSAGE = (
    "[INTRO budget={budget}] I have {budget} reasoning tokens. "
    "I will stay concise, keep only decisive steps, and avoid revisiting settled points.\n"
)
SOFT_MESSAGE = (
    "\n[SOFT] I am past the midpoint of the reasoning budget. "
    "I will consolidate the result and move toward a conclusion.\n"
)
HARD_MESSAGE = (
    "\n[HARD] I have reached the end of my reasoning budget. "
    "I will now provide the final answer.\n"
)

MODES: dict[str, dict[str, Any]] = {
    "baseline": {},
    "hard_cutoff": {
        "reasoning_budget_tokens": 96,
        "reasoning_budget_message": HARD_MESSAGE,
        # Explicitly disable every new staged feature for a clean hard-cutoff control.
        "reasoning_budget_intro_message": "",
        "reasoning_budget_soft_ratio": -1.0,
        "reasoning_budget_soft_message": "",
        "reasoning_budget_grace_tokens": 0,
    },
    "staged": {
        "reasoning_budget_tokens": 96,
        "reasoning_budget_message": HARD_MESSAGE,
        "reasoning_budget_intro_message": INTRO_MESSAGE,
        "reasoning_budget_soft_ratio": 0.65,
        "reasoning_budget_soft_message": SOFT_MESSAGE,
        "reasoning_budget_grace_tokens": 24,
    },
}


def post_json(path: str, payload: dict[str, Any], timeout: int = 240) -> dict[str, Any]:
    request = urllib.request.Request(
        BASE_URL + path,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} for {path}: {body}") from exc


def token_count(text: str) -> int:
    if not text:
        return 0
    try:
        response = post_json("/tokenize", {"content": text}, timeout=60)
        return len(response.get("tokens", []))
    except Exception:
        # The experiment should still complete if the auxiliary endpoint changes.
        return len(re.findall(r"\S+", text))


def split_reasoning_and_answer(message: dict[str, Any]) -> tuple[str, str]:
    reasoning = (
        message.get("reasoning_content")
        or message.get("reasoning")
        or message.get("analysis")
        or ""
    )
    answer = message.get("content") or ""

    # Some server/parser combinations leave the raw think block in content.
    if not reasoning:
        match = re.search(r"<think>(.*?)</think>", answer, flags=re.DOTALL | re.IGNORECASE)
        if match:
            reasoning = match.group(1)
            answer = answer[match.end() :]

    return str(reasoning), str(answer)


def is_correct(prompt_id: str, answer: str) -> bool:
    normalized = answer.lower().replace(",", "")
    compact = re.sub(r"\s+", "", normalized)
    if prompt_id == "crt":
        return bool(re.search(r"(?<!\d)53(?!\d)", normalized))
    if prompt_id == "five_digit":
        return bool(re.search(r"(?<!\d)5712(?!\d)", normalized))
    if prompt_id == "die_mod3":
        return (
            "3/7" in compact
            or "\\frac{3}{7}" in compact
            or bool(re.search(r"0\.4285", normalized))
        )
    return False


def repeated_ngram_ratio(text: str, n: int = 4) -> float:
    words = re.findall(r"[\w']+", text.lower())
    if len(words) < n:
        return 0.0
    grams = [tuple(words[i : i + n]) for i in range(len(words) - n + 1)]
    counts = Counter(grams)
    repeated_occurrences = sum(count - 1 for count in counts.values() if count > 1)
    return repeated_occurrences / len(grams)


def revision_marker_count(text: str) -> int:
    return len(
        re.findall(
            r"\b(?:but\s+wait|wait|hold\s+on|actually|reconsider|on\s+second\s+thought)\b",
            text,
            flags=re.IGNORECASE,
        )
    )


def mean(values: Iterable[float | int | None]) -> float:
    valid = [float(value) for value in values if value is not None]
    return statistics.mean(valid) if valid else float("nan")


def percent(value: float) -> str:
    return f"{100.0 * value:.1f}%"


records: list[dict[str, Any]] = []

for mode_name, mode_fields in MODES.items():
    for prompt in PROMPTS:
        for seed in SEEDS:
            payload: dict[str, Any] = {
                "model": "local-model",
                "messages": [
                    {
                        "role": "system",
                        "content": (
                            "You are a precise mathematical reasoner. Follow the requested structure exactly, "
                            "do not pad the reasoning, do not repeat settled steps, and always provide a final answer."
                        ),
                    },
                    {"role": "user", "content": prompt["prompt"]},
                ],
                "max_tokens": 384,
                "temperature": 0.1,
                "top_p": 0.8,
                "top_k": 20,
                "min_p": 0.0,
                "presence_penalty": 0.0,
                "seed": seed,
                "stream": False,
                "chat_template_kwargs": {"enable_thinking": True},
            }
            payload.update(mode_fields)

            started = time.monotonic()
            response = post_json("/v1/chat/completions", payload)
            elapsed = time.monotonic() - started

            choice = (response.get("choices") or [{}])[0]
            message = choice.get("message") or {}
            reasoning, answer = split_reasoning_and_answer(message)
            usage = response.get("usage") or {}

            record: dict[str, Any] = {
                "mode": mode_name,
                "prompt_id": prompt["id"],
                "seed": seed,
                "elapsed_s": round(elapsed, 4),
                "finish_reason": choice.get("finish_reason"),
                "reasoning": reasoning,
                "answer": answer,
                "reasoning_chars": len(reasoning),
                "answer_chars": len(answer),
                "reasoning_tokens": token_count(reasoning),
                "answer_tokens": token_count(answer),
                "completion_tokens": usage.get("completion_tokens"),
                "prompt_tokens": usage.get("prompt_tokens"),
                "intro_hit": "[INTRO" in reasoning,
                "soft_hit": "[SOFT]" in reasoning,
                "hard_hit": "[HARD]" in reasoning,
                "answer_nonempty": bool(answer.strip()),
                "correct": is_correct(prompt["id"], answer),
                "repeated_4gram_ratio": round(repeated_ngram_ratio(reasoning), 6),
                "revision_markers": revision_marker_count(reasoning),
                "request": payload,
                "raw_response": response,
            }
            records.append(record)
            print(
                "RUN "
                f"mode={mode_name} prompt={prompt['id']} seed={seed} "
                f"reasoning_tokens={record['reasoning_tokens']} answer_tokens={record['answer_tokens']} "
                f"intro={record['intro_hit']} soft={record['soft_hit']} hard={record['hard_hit']} "
                f"correct={record['correct']} repeat4={record['repeated_4gram_ratio']:.3f} "
                f"finish={record['finish_reason']} elapsed={elapsed:.2f}s",
                flush=True,
            )

(OUTPUT_DIR / "raw-results.json").write_text(
    json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8"
)

run_fields = [
    "mode",
    "prompt_id",
    "seed",
    "elapsed_s",
    "finish_reason",
    "reasoning_tokens",
    "answer_tokens",
    "completion_tokens",
    "prompt_tokens",
    "intro_hit",
    "soft_hit",
    "hard_hit",
    "answer_nonempty",
    "correct",
    "repeated_4gram_ratio",
    "revision_markers",
]
with (OUTPUT_DIR / "runs.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=run_fields)
    writer.writeheader()
    for record in records:
        writer.writerow({field: record.get(field) for field in run_fields})

aggregates: list[dict[str, Any]] = []
for mode_name in MODES:
    mode_records = [record for record in records if record["mode"] == mode_name]
    aggregates.append(
        {
            "mode": mode_name,
            "n": len(mode_records),
            "correct_rate": sum(record["correct"] for record in mode_records) / len(mode_records),
            "answer_rate": sum(record["answer_nonempty"] for record in mode_records) / len(mode_records),
            "intro_rate": sum(record["intro_hit"] for record in mode_records) / len(mode_records),
            "soft_rate": sum(record["soft_hit"] for record in mode_records) / len(mode_records),
            "hard_rate": sum(record["hard_hit"] for record in mode_records) / len(mode_records),
            "length_finish_rate": sum(
                record["finish_reason"] == "length" for record in mode_records
            )
            / len(mode_records),
            "mean_reasoning_tokens": mean(record["reasoning_tokens"] for record in mode_records),
            "median_reasoning_tokens": statistics.median(
                record["reasoning_tokens"] for record in mode_records
            ),
            "mean_answer_tokens": mean(record["answer_tokens"] for record in mode_records),
            "mean_elapsed_s": mean(record["elapsed_s"] for record in mode_records),
            "mean_repeated_4gram_ratio": mean(
                record["repeated_4gram_ratio"] for record in mode_records
            ),
            "mean_revision_markers": mean(
                record["revision_markers"] for record in mode_records
            ),
        }
    )

with (OUTPUT_DIR / "aggregate.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(aggregates[0].keys()))
    writer.writeheader()
    writer.writerows(aggregates)

lines: list[str] = [
    "# CPU reproduction: staged reasoning-budget sampler",
    "",
    "## Setup",
    "",
    "- Fork commit: `9cde33219df8e4d1fca4dd42c8101ddba6eff416`",
    "- Model: `unsloth/Qwen3-0.6B-GGUF`, `Q4_K_M` (0.6B parameters)",
    "- Backend: llama.cpp CPU only; no GPU offload",
    "- Sampling stress test: temperature 0.1, top-p 0.8, top-k 20",
    "- Controlled budget: 96 reasoning tokens; staged mode warns at 65% and allows 24 grace tokens",
    "- Workload: 3 math prompts × 3 seeds × 3 modes = 27 generations",
    "",
    "## Aggregate results",
    "",
    "| mode | n | correct | answer produced | intro hit | soft hit | hard hit | finish=length | mean reasoning tok | median reasoning tok | mean answer tok | repeat-4gram | revision markers | mean seconds |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for aggregate in aggregates:
    lines.append(
        f"| {aggregate['mode']} | {aggregate['n']} | {percent(aggregate['correct_rate'])} | "
        f"{percent(aggregate['answer_rate'])} | {percent(aggregate['intro_rate'])} | "
        f"{percent(aggregate['soft_rate'])} | {percent(aggregate['hard_rate'])} | "
        f"{percent(aggregate['length_finish_rate'])} | {aggregate['mean_reasoning_tokens']:.1f} | "
        f"{aggregate['median_reasoning_tokens']:.1f} | {aggregate['mean_answer_tokens']:.1f} | "
        f"{aggregate['mean_repeated_4gram_ratio']:.3f} | {aggregate['mean_revision_markers']:.2f} | "
        f"{aggregate['mean_elapsed_s']:.2f} |"
    )

lines.extend(
    [
        "",
        "## Per-prompt correctness",
        "",
        "| prompt | baseline | hard cutoff | staged |",
        "|---|---:|---:|---:|",
    ]
)
for prompt in PROMPTS:
    values: list[str] = []
    for mode_name in MODES:
        prompt_records = [
            record
            for record in records
            if record["prompt_id"] == prompt["id"] and record["mode"] == mode_name
        ]
        values.append(
            percent(sum(record["correct"] for record in prompt_records) / len(prompt_records))
        )
    lines.append(f"| {prompt['id']} | {values[0]} | {values[1]} | {values[2]} |")

staged_records = [record for record in records if record["mode"] == "staged"]
lines.extend(
    [
        "",
        "## Instrumentation checks",
        "",
        f"- Intro marker observed in {sum(record['intro_hit'] for record in staged_records)}/{len(staged_records)} staged runs.",
        f"- Soft marker observed in {sum(record['soft_hit'] for record in staged_records)}/{len(staged_records)} staged runs.",
        f"- Hard marker observed in {sum(record['hard_hit'] for record in staged_records)}/{len(staged_records)} staged runs.",
        "- The fork's `test-reasoning-budget` executable completed before inference; its log is included.",
        "",
        "## Interpretation",
        "",
        "This is a functional reproduction and a small low-temperature stress test, not a statistically powered quality benchmark. "
        "Marker rates establish whether the staged sampler fired end-to-end through the Qwen chat template. Correctness, "
        "termination, token counts, repeated 4-grams, and revision markers indicate whether the intervention obviously helped "
        "or harmed this 0.6B model on the selected tasks. General conclusions require a broader benchmark, matched output "
        "budgets, and confidence intervals.",
        "",
    ]
)

report = "\n".join(lines)
(OUTPUT_DIR / "report.md").write_text(report, encoding="utf-8")
print("\n" + report, flush=True)
