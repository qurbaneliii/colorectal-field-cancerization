from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
import tarfile
from pathlib import Path

from src.data.platform import detect_chip_type, validate_chip_type
from src.data.provenance import file_sha256


GEO_PLATFORMS = {"GSE44076": ("GPL13667", "HG-U219"), "GSE41258": ("GPL96", "HG-U133A")}


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def archive_spec(root: Path, accession: str) -> dict:
    metadata = root / "data/metadata"
    manifest = csv_rows(metadata / "data_manifest.csv")
    row, = [
        entry for entry in manifest
        if entry["accession"] == accession and entry["filename"] == f"{accession}_RAW.tar"
    ]
    platform, chip = GEO_PLATFORMS[accession]
    samples = csv_rows(metadata / f"{accession.lower()}_samples.csv")
    members = csv_rows(metadata / f"{accession}_tar_members.csv")
    sample_by_file = {entry["cel_filename"]: entry["geo_accession"] for entry in samples}
    order = [entry["member"] for entry in members]
    expected_count = int(row["archive_member_count"])
    if (
        len(samples) != expected_count
        or len(sample_by_file) != expected_count
        or len(set(sample_by_file.values())) != expected_count
        or len(order) != expected_count
        or set(order) != set(sample_by_file)
        or {entry["platform"] for entry in samples} != {platform}
    ):
        raise ValueError(f"{accession}: historical CEL manifest and GEO sample metadata disagree")
    return {
        "accession": accession,
        "filename": row["filename"],
        "bytes": int(row["size_bytes"]),
        "historical_sha256": row["sha256"],
        "url": (
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/"
            f"{accession}/suppl/{row['filename']}"
        ),
        "chip": chip,
        "sample_by_file": sample_by_file,
        "historical_order": order,
    }


def cel_payload(compressed: bytes) -> tuple[str, int, str]:
    """Read the entire gzip stream, including its CRC and length trailer."""
    sha = hashlib.sha256()
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as source:
        header = source.read(256 * 1024)
        chip = detect_chip_type(header)
        sha.update(header)
        size = len(header)
        while chunk := source.read(1024 * 1024):
            sha.update(chunk)
            size += len(chunk)
    return sha.hexdigest(), size, chip


def inspect_cel_archive(path: Path, spec: dict) -> tuple[dict, list[dict]]:
    """Verify a GEO CEL tar and fingerprint the archive envelope and each CEL."""
    expected = spec["sample_by_file"]
    seen: set[str] = set()
    rows = []
    with tarfile.open(path, "r:") as archive, path.open("rb") as raw:
        for order, member in enumerate(archive, start=1):
            if not member.isfile() or member.name not in expected or member.name in seen:
                raise ValueError(f"Unexpected or duplicate CEL member {member.name!r} in {path}")
            gsm = re.match(r"^(GSM\d+)_", member.name)
            if gsm is None or gsm[1] != expected[member.name]:
                raise ValueError(f"CEL filename and metadata disagree for {member.name}")
            seen.add(member.name)
            source = archive.extractfile(member)
            if source is None:
                raise ValueError(f"Cannot read {member.name} from {path}")
            with source:
                compressed = source.read()
            payload_hash, payload_bytes, chip = cel_payload(compressed)
            validate_chip_type(chip, spec["chip"])
            raw.seek(member.offset)
            tar_header = raw.read(member.offset_data - member.offset)
            rows.append({
                "accession": spec["accession"],
                "order": order,
                "geo_accession": gsm[1],
                "member": member.name,
                "platform_chip": chip,
                "tar_offset": member.offset,
                "tar_header_sha256": hashlib.sha256(tar_header).hexdigest(),
                "tar_mtime_unix": member.mtime,
                "gzip_bytes": len(compressed),
                "gzip_sha256": hashlib.sha256(compressed).hexdigest(),
                "cel_bytes": payload_bytes,
                "cel_sha256": payload_hash,
            })
        if seen != set(expected):
            raise ValueError(f"Missing CEL members: {sorted(set(expected) - seen)[:10]}")
        end = member.offset_data + ((member.size + 511) // 512) * 512
        raw.seek(end)
        trailer = raw.read()
        if len(trailer) < 1024 or any(trailer):
            raise ValueError(f"{path} lacks a complete zero-padded TAR footer")
    observed_sha = file_sha256(path)
    names = [row["member"] for row in rows]
    return {
        "accession": spec["accession"],
        "source_url": spec["url"],
        "archive_bytes": path.stat().st_size,
        "historical_archive_bytes": spec["bytes"],
        "archive_sha256": observed_sha,
        "historical_archive_sha256": spec["historical_sha256"],
        "historical_sha256_matches": observed_sha == spec["historical_sha256"],
        "member_count": len(rows),
        "historical_member_order_matches": names == spec["historical_order"],
        "tar_zero_trailer_bytes": len(trailer),
        "all_members_match_metadata_and_platform": True,
    }, rows


def compare_archives(current: list[dict], historical: list[dict]) -> dict:
    old = {row["member"]: row for row in historical}
    new = {row["member"]: row for row in current}
    shared = old.keys() & new.keys()
    gzip_differences = sorted(name for name in shared if old[name]["gzip_sha256"] != new[name]["gzip_sha256"])
    cel_differences = sorted(name for name in shared if old[name]["cel_sha256"] != new[name]["cel_sha256"])
    return {
        "removed_members": sorted(old.keys() - new.keys()),
        "added_members": sorted(new.keys() - old.keys()),
        "member_order_matches": list(old) == list(new),
        "tar_header_differences": sorted(
            name for name in shared if old[name]["tar_header_sha256"] != new[name]["tar_header_sha256"]
        ),
        "gzip_differences": gzip_differences,
        "cel_differences": cel_differences,
        "interpretation": (
            "member_set_changed" if old.keys() != new.keys()
            else "cel_payload_changed" if cel_differences
            else "gzip_recompression_or_header_change" if gzip_differences
            else "cel_members_identical_archive_envelope_may_differ"
        ),
    }
