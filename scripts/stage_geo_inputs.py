from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.geo_archive import archive_spec, csv_rows, inspect_cel_archive
from src.data.provenance import file_sha256
from src.data.resumable_download import resumable_download


def download_inputs(cache: Path, alternate_pins: dict[str, str]) -> None:
    for row in csv_rows(ROOT / "data/metadata/data_manifest.csv"):
        accession, filename = row["accession"], row["filename"]
        geo_root = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}"
        folder = "matrix" if filename.endswith("_series_matrix.txt.gz") else "suppl"
        path = cache / filename
        transfer = resumable_download(
            f"{geo_root}/{folder}/{filename}", path, int(row["size_bytes"]),
            alternate_pins.get(filename, row["sha256"]),
        )
        if filename.endswith("_RAW.tar"):
            summary, _ = inspect_cel_archive(path, archive_spec(ROOT, accession))
            print(f"{accession}: verified {summary['member_count']} CEL streams", flush=True)
        print(f"{filename}: {transfer['mode']} {transfer['sha256']}", flush=True)


def stage(cache: Path, worktree: Path, alternate_pins: dict[str, str]) -> dict:
    if worktree.resolve() == ROOT or not (worktree / ".git").exists():
        raise ValueError("Use a separate Git worktree for raw-data reproduction")
    historical = csv_rows(ROOT / "data/metadata/data_manifest.csv")
    unknown = alternate_pins.keys() - {row["filename"] for row in historical}
    if unknown:
        raise ValueError(f"Unknown alternate inputs: {sorted(unknown)}")
    files = []
    for row in historical:
        source = cache / row["filename"]
        target = worktree / row["relative_path"]
        expected = alternate_pins.get(row["filename"], row["sha256"])
        size = int(row["size_bytes"])
        if source.stat().st_size != size or file_sha256(source) != expected:
            raise ValueError(f"Input failed its explicit size or SHA-256 pin: {source.name}")
        if target.exists() and (
            target.stat().st_size != size or file_sha256(target) != expected
        ):
            raise ValueError(f"Existing worktree input differs: {target}")
        files.append((row, source, target, expected))

    for _, source, target, _ in files:
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.link(source, target)
            except OSError:
                with source.open("rb") as raw, target.open("xb") as staged:
                    shutil.copyfileobj(raw, staged, length=1024 * 1024)
    return {
        "mode": "alternate_pinned_inputs" if alternate_pins else "historical_inputs",
        "inputs": [
            {
                "accession": row["accession"],
                "relative_path": row["relative_path"],
                "bytes": int(row["size_bytes"]),
                "historical_sha256": row["sha256"],
                "verified_sha256": expected,
            }
            for row, _, _, expected in files
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--worktree", type=Path, required=True)
    parser.add_argument("--alternate-pins", type=Path)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--download", action="store_true", help="Acquire missing pinned files first")
    args = parser.parse_args()
    pins = (
        json.loads(args.alternate_pins.read_text(encoding="utf-8"))["archives"]
        if args.alternate_pins
        else {}
    )
    if args.download:
        download_inputs(args.cache, pins)
    evidence = stage(args.cache, args.worktree, pins)
    args.evidence.parent.mkdir(parents=True, exist_ok=True)
    args.evidence.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(f"Verified and staged {len(evidence['inputs'])} inputs under {evidence['mode']}")


if __name__ == "__main__":
    main()
