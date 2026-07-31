from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path


SUPPORTED_CHIPS = ("HG-U219", "HG-U133A")


def _normalized(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def detect_chip_type(header: bytes) -> str:
    """Detect supported Affymetrix chip identifiers in text or Calvin CEL headers."""

    for chip in SUPPORTED_CHIPS:
        variants = {
            chip.encode("ascii"),
            chip.encode("utf-16-be"),
            chip.encode("utf-16-le"),
        }
        if any(variant in header for variant in variants):
            return chip
    raise ValueError("No supported Affymetrix chip identifier found in CEL header")


def validate_chip_type(detected: str, expected: str) -> None:
    if _normalized(detected) != _normalized(expected):
        raise ValueError(f"CEL chip mismatch: expected {expected}, detected {detected}")


def detect_archive_chip_type(tar_path: Path, read_bytes: int = 2_000_000) -> str:
    with tarfile.open(tar_path, "r") as archive:
        member = next(
            item
            for item in archive.getmembers()
            if item.isfile() and item.name.lower().endswith(".cel.gz")
        )
        raw = archive.extractfile(member)
        if raw is None:
            raise ValueError(f"Could not read {member.name} from {tar_path}")
        with gzip.GzipFile(fileobj=raw) as handle:
            header = handle.read(read_bytes)
    return detect_chip_type(header)
