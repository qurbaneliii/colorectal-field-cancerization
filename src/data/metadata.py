from __future__ import annotations

import csv
import gzip
import re
import tarfile
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


def discover_accession_files(raw_root: Path, accession: str) -> dict[str, list[Path]]:
    """Discover inputs robustly by accession and extension."""
    files = [p for p in raw_root.rglob("*") if p.is_file() and accession.lower() in p.name.lower()]
    return {
        "raw_tar": sorted(p for p in files if p.name.lower().endswith("_raw.tar")),
        "series_matrix": sorted(
            p for p in files if "series_matrix" in p.name.lower() and p.name.lower().endswith(".gz")
        ),
        "clinical": sorted(
            p for p in files if "clinical" in p.name.lower() and p.name.lower().endswith(".gz")
        ),
    }


def parse_series_metadata(path: Path) -> dict[str, list[list[str]]]:
    """Parse all repeated GEO series-matrix metadata rows without loading expression."""
    result: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            row = next(csv.reader([line], delimiter="\t"))
            result.setdefault(row[0], []).append(row[1:])
    return result


def _sample_characteristics(meta: dict[str, list[list[str]]], n: int) -> list[dict[str, str]]:
    output = [dict() for _ in range(n)]
    for characteristic_row in meta.get("!Sample_characteristics_ch1", []):
        if len(characteristic_row) != n:
            raise ValueError("Sample characteristic row length does not match sample count")
        for idx, value in enumerate(characteristic_row):
            if ":" not in value:
                continue
            key, val = value.split(":", 1)
            output[idx][key.strip().lower()] = val.strip()
    return output


def _single_sample_row(meta: dict[str, list[list[str]]], key: str, n: int) -> list[str]:
    rows = meta.get(key, [])
    if len(rows) != 1 or len(rows[0]) != n:
        raise ValueError(f"Expected exactly one complete {key} row")
    return rows[0]


def _supplementary_basename(value: str) -> str:
    return Path(urlparse(value).path).name


def tar_members(tar_path: Path) -> list[str]:
    with tarfile.open(tar_path, "r") as archive:
        return archive.getnames()


def build_gse44076_metadata(series_path: Path, tar_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_series_metadata(series_path)
    accessions = meta["!Sample_geo_accession"][0]
    n = len(accessions)
    titles = _single_sample_row(meta, "!Sample_title", n)
    sources = _single_sample_row(meta, "!Sample_source_name_ch1", n)
    platforms = _single_sample_row(meta, "!Sample_platform_id", n)
    supplementary = _single_sample_row(meta, "!Sample_supplementary_file", n)
    chars = _sample_characteristics(meta, n)
    members = set(tar_members(tar_path))

    source_to_class = {
        "healthy colon mucosa cells": "healthy",
        "normal distant colon mucosa cells": "adjacent_normal",
        "primary colon adenocarcinoma cells": "tumor",
    }
    rows: list[dict[str, object]] = []
    for gsm, title, source, platform, supp, char in zip(
        accessions, titles, sources, platforms, supplementary, chars, strict=True
    ):
        tissue = source_to_class.get(source.strip().lower())
        patient_id = char.get("individual id")
        cel = _supplementary_basename(supp)
        resolved = bool(tissue and patient_id)
        rows.append(
            {
                "geo_accession": gsm,
                "sample_title": title,
                "original_tissue": source,
                "tissue_class": tissue or "unresolved",
                "patient_id": patient_id or "",
                "donor_or_patient_group": patient_id or f"unresolved_{gsm}",
                "pairing_status": "not_applicable",
                "platform": platform,
                "cel_filename": cel,
                "cel_present_in_tar": cel in members,
                "inclusion_status": "included" if resolved and cel in members else "excluded",
                "exclusion_reason": "" if resolved and cel in members else "unresolved metadata or CEL missing",
                "metadata_parsing_source": "GEO source_name_ch1 and individual id characteristic",
                "metadata_parsing_confidence": "high" if resolved else "low",
                "stage": char.get("stage", ""),
                "location": char.get("location", ""),
                "sex": char.get("gender", ""),
                "age": char.get("age", ""),
            }
        )
    frame = pd.DataFrame(rows)
    pair_counts = (
        frame[frame["tissue_class"].isin(["adjacent_normal", "tumor"])]
        .groupby(["patient_id", "tissue_class"])["geo_accession"]
        .agg(list)
        .unstack(fill_value=[])
    )
    for col in ["adjacent_normal", "tumor"]:
        if col not in pair_counts:
            pair_counts[col] = [[] for _ in range(len(pair_counts))]
    pair_audit = pd.DataFrame(
        {
            "patient_id": pair_counts.index,
            "adjacent_normal_sample": pair_counts["adjacent_normal"].map(
                lambda x: ";".join(x) if isinstance(x, list) else ""
            ),
            "tumor_sample": pair_counts["tumor"].map(
                lambda x: ";".join(x) if isinstance(x, list) else ""
            ),
            "both_present": [
                isinstance(a, list) and len(a) == 1 and isinstance(t, list) and len(t) == 1
                for a, t in zip(pair_counts["adjacent_normal"], pair_counts["tumor"], strict=True)
            ],
            "duplicate_samples": [
                (isinstance(a, list) and len(a) > 1) or (isinstance(t, list) and len(t) > 1)
                for a, t in zip(pair_counts["adjacent_normal"], pair_counts["tumor"], strict=True)
            ],
        }
    )
    paired_ids = set(pair_audit.loc[pair_audit["both_present"], "patient_id"])
    paired_mask = frame["patient_id"].isin(paired_ids) & frame["tissue_class"].isin(
        ["adjacent_normal", "tumor"]
    )
    frame.loc[paired_mask, "pairing_status"] = "paired"
    frame.loc[frame["tissue_class"].eq("healthy"), "pairing_status"] = "healthy_unique_donor"
    return frame, pair_audit


def _canonical_external_tissue(value: str) -> str:
    normalized = re.sub(r"\s+", " ", value.strip().lower())
    mapping = {
        "primary tumor": "primary_tumor",
        "normal colon": "normal_colon",
        "polyp": "polyp",
        "polyp, high grade": "polyp",
        "microadenoma": "polyp",
        "liver metastasis": "liver_metastasis",
        "lung metastasis": "lung_metastasis",
        "normal liver": "normal_liver",
        "normal lung": "normal_lung",
        "cell line": "cell_line",
    }
    return mapping.get(normalized, "other")


def build_gse41258_metadata(
    series_path: Path, tar_path: Path, technical_pattern: str
) -> pd.DataFrame:
    meta = parse_series_metadata(series_path)
    accessions = meta["!Sample_geo_accession"][0]
    n = len(accessions)
    titles = _single_sample_row(meta, "!Sample_title", n)
    platforms = _single_sample_row(meta, "!Sample_platform_id", n)
    supplementary = _single_sample_row(meta, "!Sample_supplementary_file", n)
    chars = _sample_characteristics(meta, n)
    members = set(tar_members(tar_path))
    tech_re = re.compile(technical_pattern)
    rows: list[dict[str, object]] = []
    for gsm, title, platform, supp, char in zip(
        accessions, titles, platforms, supplementary, chars, strict=True
    ):
        original_tissue = char.get("tissue", "")
        tissue = _canonical_external_tissue(original_tissue)
        patient = char.get("patient id", "")
        author_included = char.get("included in analysis", "").lower() == "yes"
        technical = bool(tech_re.search(title))
        cel = _supplementary_basename(supp)
        reasons: list[str] = []
        if tissue not in {"primary_tumor", "normal_colon"}:
            reasons.append("not primary external-validation tissue")
        if not author_included:
            reasons.append("GEO author flag included in analysis is not Yes")
        if technical:
            reasons.append("technical-replicate indicator")
        if cel not in members:
            reasons.append("CEL missing from TAR")
        if not patient:
            reasons.append("missing explicit patient id")
        rows.append(
            {
                "geo_accession": gsm,
                "sample_title": title,
                "original_tissue": original_tissue,
                "tissue_class": tissue,
                "patient_id": patient,
                "donor_or_patient_group": patient or f"unresolved_{gsm}",
                "pairing_status": "patient_grouped" if patient else "unresolved",
                "platform": platform,
                "cel_filename": cel,
                "cel_present_in_tar": cel in members,
                "author_included": author_included,
                "technical_replicate_candidate": technical,
                "technical_replicate_rule": technical_pattern,
                "inclusion_status": "included" if not reasons else "excluded",
                "exclusion_reason": "; ".join(reasons),
                "metadata_parsing_source": "GEO tissue, patient id, and included characteristics",
                "metadata_parsing_confidence": "high" if original_tissue and patient else "low",
            }
        )
    return pd.DataFrame(rows)
