from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_geo_metadata import samples
from src.data.geo_archive import (
    GEO_PLATFORMS,
    archive_spec,
    cel_payload,
    compare_archives,
    inspect_cel_archive,
)
from src.data.resumable_download import resumable_download


def compare_live_sample_files(rows: list[dict], geo_text: Path, count: int) -> list[dict]:
    """Compare evenly spaced archive members with their independently served GEO CELs."""
    records = samples(geo_text)
    indices = (
        range(len(rows))
        if count >= len(rows)
        else sorted({round(i * (len(rows) - 1) / max(count - 1, 1)) for i in range(count)})
    )
    checks = []
    for index in indices:
        row = rows[index]
        url = records[row["geo_accession"]]["supplementary_file"].replace(
            "ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/", 1
        )
        if not url.startswith("https://ftp.ncbi.nlm.nih.gov/geo/samples/"):
            raise ValueError(f"Unexpected sample file endpoint for {row['geo_accession']}")
        with urlopen(url, timeout=45) as response:
            data = response.read()
        compressed_sha = hashlib.sha256(data).hexdigest()
        payload_sha, _, chip = cel_payload(data)
        checks.append({
            "geo_accession": row["geo_accession"],
            "source_url": url,
            "source_gzip_sha256": compressed_sha,
            "source_cel_sha256": payload_sha,
            "chip_matches": row["platform_chip"] == chip,
            "gzip_matches": row["gzip_sha256"] == compressed_sha,
            "cel_matches": row["cel_sha256"] == payload_sha,
        })
    return checks


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("accession", choices=GEO_PLATFORMS)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "results/ai_audit")
    parser.add_argument("--accept-sha256", help="Explicit alternate archive pin, separate from historical")
    parser.add_argument("--inspect-only", type=Path, help="Audit an existing archive without promoting it")
    parser.add_argument("--historical-archive", type=Path)
    parser.add_argument("--geo-text", type=Path)
    parser.add_argument("--live-samples", type=int, default=0)
    args = parser.parse_args()
    spec = archive_spec(ROOT, args.accession)
    pin = args.accept_sha256 or spec["historical_sha256"]
    if not re.fullmatch(r"[0-9a-f]{64}", pin):
        parser.error("Expected a lowercase 64-character SHA-256")
    if args.live_samples < 0 or (args.live_samples and not args.geo_text):
        parser.error("--live-samples requires a nonnegative count and --geo-text")
    archive = args.inspect_only or args.cache / spec["filename"]
    transfer = {"mode": "inspect_only"}
    if not args.inspect_only:
        transfer = resumable_download(spec["url"], archive, spec["bytes"], pin)
    report, members = inspect_cel_archive(archive, spec)
    report["transfer"] = transfer
    if args.historical_archive:
        historical, old_members = inspect_cel_archive(args.historical_archive, spec)
        if not historical["historical_sha256_matches"]:
            raise ValueError("The historical comparison archive does not match the historical digest")
        report["historical_byte_comparison"] = compare_archives(members, old_members)
    args.output.mkdir(parents=True, exist_ok=True)
    stem = args.accession.lower() + "_raw_integrity"
    write_csv(args.output / f"{stem}_members.csv", members)
    if args.live_samples:
        checks = compare_live_sample_files(members, args.geo_text, args.live_samples)
        write_csv(args.output / f"{stem}_live_samples.csv", checks)
        report["live_sample_checks"] = {
            "count": len(checks),
            "gzip_matches": sum(row["gzip_matches"] for row in checks),
            "cel_matches": sum(row["cel_matches"] for row in checks),
        }
    (args.output / f"{stem}.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
