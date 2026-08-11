"""Repair one transport-damaged frozen benchmark chunk before Python jobs run."""
from __future__ import annotations

import hashlib
from pathlib import Path

_PATH = Path(__file__).resolve().parent / "benchmarks" / "sped_ro_medical_robust" / "bootstrap.b64.03a"
_TARGET = "f8747e71fc43164437e6e52549b4476877fb0929b5f3d115fa8544cd4af79937"
_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

if _PATH.exists():
    _text = _PATH.read_text(encoding="utf-8").strip()
    if hashlib.sha256(_text.encode()).hexdigest() != _TARGET:
        if len(_text) != 3499:
            raise RuntimeError(f"unexpected frozen chunk length: {len(_text)}")
        _matches: list[tuple[int, str, str]] = []
        for _position in range(len(_text) + 1):
            for _character in _ALPHABET:
                _candidate = _text[:_position] + _character + _text[_position:]
                if hashlib.sha256(_candidate.encode()).hexdigest() == _TARGET:
                    _matches.append((_position, _character, _candidate))
        if len(_matches) != 1:
            raise RuntimeError(f"unable to uniquely repair frozen chunk: {len(_matches)} matches")
        _position, _character, _repaired = _matches[0]
        _PATH.write_text(_repaired, encoding="utf-8")
        print(f"Repaired frozen benchmark chunk: inserted {_character!r} at {_position}")
