#!/usr/bin/env python3
"""Run the SpeD benchmark with a narrow Hugging Face download workaround.

huggingface_hub 1.16 currently raises ValueError inside snapshot_download when
an allow-pattern targets only nested files in rhasspy/piper-voices. Downloading
the exact files with hf_hub_download avoids that client bug. The SpeD model
snapshot path remains unchanged.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import huggingface_hub
from huggingface_hub import hf_hub_download

_original_snapshot_download = huggingface_hub.snapshot_download


def _snapshot_download(*args: Any, **kwargs: Any) -> str:
    repo_id = kwargs.get("repo_id") or (args[0] if args else None)
    if repo_id == "rhasspy/piper-voices":
        local_dir = Path(kwargs["local_dir"])
        local_dir.mkdir(parents=True, exist_ok=True)
        patterns = kwargs.get("allow_patterns") or []
        for filename in patterns:
            hf_hub_download(repo_id=repo_id, filename=filename, local_dir=local_dir)
        return str(local_dir)
    return str(_original_snapshot_download(*args, **kwargs))


huggingface_hub.snapshot_download = _snapshot_download

import sped_tts_benchmark  # noqa: E402  (must import after monkeypatch)

raise SystemExit(sped_tts_benchmark.main())
