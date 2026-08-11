"""Verify and repair frozen benchmark transport chunks before CI.

Large source chunks 02 and 06 are rebuilt from independently transferred
1,000-character fragments. Every reconstructed chunk is then checked against a
frozen SHA-256. A final one-character repair is permitted only when it produces
the unique expected hash. The corpus, lexicon, thresholds, and metrics are never
modified here.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent / "benchmarks" / "sped_ro_medical_robust"
ALPHABET = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
EXPECTED: dict[str, tuple[int, str]] = {
    "bootstrap.b64.00": (7000, "7c8b0b267f9464daaa60b93bf7938f70657ed2a1e1aa88620e0cb89000558d5a"),
    "bootstrap.b64.01": (7000, "5e84de689ee780e82385c2fd22b7e7a31bd43a958e3faaf5aa8283575ea45fc8"),
    "bootstrap.b64.02": (7000, "a353cb2998f0b153d0e440a52455ae3baf34fa1f33c72be8ed95ea925316a593"),
    "bootstrap.b64.03a": (3500, "f8747e71fc43164437e6e52549b4476877fb0929b5f3d115fa8544cd4af79937"),
    "bootstrap.b64.03b": (3500, "5c31ce0ca94513041d679de1f910d047e5f3d66cf5a5298bd53c3fd73b96e6a0"),
    "bootstrap.b64.04": (7000, "d112d62894199a17ff080256bc8d6082a2e657b1ed0733f0a42d5ec3e2583135"),
    "bootstrap.b64.05": (7000, "ab1714b4efb8ee8c37af9c9274690afef3bfa382e33ad13842d9b685df0d6285"),
    "bootstrap.b64.06": (6796, "85f071379117a638ce195de0b74d1b6f05e61c42c29dd0d65d88a27ac31dd8d2"),
}
FRAGMENT_GROUPS: dict[str, list[str]] = {
    "bootstrap.b64.02": [f"bootstrap.b64.02.{index:02d}" for index in range(7)],
    "bootstrap.b64.06": [f"bootstrap.b64.06.{index:02d}" for index in range(7)],
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def repair_one_edit(data: bytes, expected_len: int, expected_hash: str) -> tuple[bytes, str]:
    matches: list[tuple[bytes, str]] = []

    if len(data) == expected_len - 1:
        for pos in range(len(data) + 1):
            prefix, suffix = data[:pos], data[pos:]
            for char in ALPHABET:
                candidate = prefix + bytes((char,)) + suffix
                if digest(candidate) == expected_hash:
                    matches.append((candidate, f"inserted {chr(char)!r} at {pos}"))

    elif len(data) == expected_len + 1:
        for pos in range(len(data)):
            candidate = data[:pos] + data[pos + 1 :]
            if digest(candidate) == expected_hash:
                matches.append((candidate, f"deleted {chr(data[pos])!r} at {pos}"))

    elif len(data) == expected_len:
        mutable = bytearray(data)
        for pos, old in enumerate(data):
            for char in ALPHABET:
                if char == old:
                    continue
                mutable[pos] = char
                candidate = bytes(mutable)
                if digest(candidate) == expected_hash:
                    matches.append((candidate, f"replaced {chr(old)!r} with {chr(char)!r} at {pos}"))
            mutable[pos] = old
    else:
        raise RuntimeError(
            f"length differs by more than one byte: got {len(data)}, expected {expected_len}"
        )

    if len(matches) != 1:
        raise RuntimeError(f"found {len(matches)} matching one-edit repairs")
    return matches[0]


if ROOT.exists():
    for target_name, fragment_names in FRAGMENT_GROUPS.items():
        fragment_paths = [ROOT / name for name in fragment_names]
        if all(path.exists() for path in fragment_paths):
            assembled = b"".join(path.read_bytes().strip() for path in fragment_paths)
            (ROOT / target_name).write_bytes(assembled)
            print(
                f"Rebuilt {target_name} from {len(fragment_paths)} fragments: "
                f"len={len(assembled)} sha256={digest(assembled)}"
            )

    errors: list[str] = []
    for name, (expected_len, expected_hash) in EXPECTED.items():
        path = ROOT / name
        if not path.exists():
            errors.append(f"{name}: missing")
            continue
        data = path.read_text(encoding="utf-8").strip().encode("ascii")
        actual_hash = digest(data)
        if len(data) == expected_len and actual_hash == expected_hash:
            print(f"OK {name}: len={len(data)} sha256={actual_hash}")
            continue

        print(
            f"MISMATCH {name}: len={len(data)} expected_len={expected_len} "
            f"sha256={actual_hash} expected_sha256={expected_hash}"
        )
        try:
            repaired, description = repair_one_edit(data, expected_len, expected_hash)
        except Exception as exc:
            errors.append(f"{name}: {exc}")
            continue
        path.write_bytes(repaired)
        print(f"REPAIRED {name}: {description}")

    for name, (expected_len, expected_hash) in EXPECTED.items():
        path = ROOT / name
        if not path.exists():
            continue
        data = path.read_bytes().strip()
        actual_hash = digest(data)
        if len(data) != expected_len or actual_hash != expected_hash:
            errors.append(
                f"{name}: verification failed len={len(data)} sha256={actual_hash}"
            )

    if errors:
        raise RuntimeError("frozen chunk errors:\n- " + "\n- ".join(dict.fromkeys(errors)))
    print(f"Verified {len(EXPECTED)} frozen benchmark chunks")
