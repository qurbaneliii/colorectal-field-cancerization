from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_LABELS = {
    "Healthy colon mucosa cells": "healthy",
    "Normal distant colon mucosa cells": "adjacent_normal",
    "Primary colon adenocarcinoma cells": "tumor",
}
EXTERNAL_LABELS = {"Primary Tumor": "primary_tumor", "Normal Colon": "normal_colon"}
CLINICAL_FIELDS = {
    "Sample Date": "sample_date",
    "Age": "age",
    "Gender": "sex",
    "T": "t_stage",
    "N": "n_stage",
    "M": "m_stage",
    "Group Stage": "stage",
    "Microsattelite instability": "msi_status",
    "P53 mutation status": "p53_status",
    "Anatomic Location": "location",
}


def samples(path: Path) -> dict[str, dict[str, str]]:
    """Read GEO text records independently of the series-matrix parser."""
    records: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if " = " not in line:
            continue
        field, value = line.split(" = ", 1)
        if field == "^SAMPLE":
            if value in records:
                raise ValueError(f"Duplicate GEO sample {value}")
            current = {}
            records[value] = current
        elif field.startswith("^"):
            current = None
        elif current is not None:
            if field == "!Sample_characteristics_ch1":
                name, _, characteristic = value.partition(":")
                current[name.lower().strip()] = characteristic.strip()
            elif field in {
                "!Sample_title",
                "!Sample_source_name_ch1",
                "!Sample_platform_id",
                "!Sample_supplementary_file",
            }:
                current[field.removeprefix("!Sample_")] = value
    if not records:
        raise ValueError(f"No GEO sample records in {path}")
    return records


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as source:
        return list(csv.DictReader(source))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def clinical_values(path: Path) -> dict:
    fields: defaultdict = defaultdict(lambda: defaultdict(set))
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as source:
        for record in csv.DictReader(source, delimiter="\t"):
            row = {" ".join(key.split()): (value or "").strip() for key, value in record.items()}
            for original, name in CLINICAL_FIELDS.items():
                if row.get(original, "").lower() not in {"", "--", "na", "n/a"}:
                    fields[row["Patient ID"]][name].add(row[original])
    return {
        patient: {name: next(iter(values)) for name, values in columns.items() if len(values) == 1}
        for patient, columns in fields.items()
    }


def reconcile(
    root: Path, accession: str, source: Path, update_scan_dates: bool, clinical: dict | None = None
) -> tuple[list[dict], list[dict], list[dict]]:
    records = samples(source)
    metadata_path = root / "data/metadata" / f"{accession.lower()}_samples.csv"
    metadata = read_csv(metadata_path)
    primary = accession == "GSE44076"
    recorded = {row["geo_accession"]: row for row in metadata}
    if len(metadata) != len(recorded) or set(records) != set(recorded):
        raise ValueError(f"{accession} GEO accessions disagree with the repository")
    disagreements: list[dict] = []
    labels: Counter = Counter()
    source_labels: Counter = Counter()
    exclusion_counts: Counter = Counter()
    donor_tissue: defaultdict = defaultdict(set)
    archive_cel = set()
    for row in read_csv(root / "data/metadata" / f"{accession}_tar_members.csv"):
        for value in row.values():
            if value.lower().endswith(".cel.gz"):
                archive_cel.add(value.replace("\\", "/").rsplit("/", 1)[-1])
    for gsm, record in records.items():
        stored = recorded[gsm]
        cel = record.get("supplementary_file", "").rsplit("/", 1)[-1]
        expected = {
            "sample_title": record["title"],
            "platform": record["platform_id"],
            "cel_filename": cel,
            "cel_present_in_tar": str(cel in archive_cel),
        }
        if primary:
            tissue = record["source_name_ch1"]
            label = PRIMARY_LABELS[tissue]
            donor = record["individual id"]
            eligible = bool(donor and cel in archive_cel)
            for field in ("age", "location", "stage"):
                expected[field] = record.get(field, "")
            expected["sex"] = record.get("gender", "")
        else:
            tissue = record.get("tissue", "")
            label = EXTERNAL_LABELS.get(tissue)
            donor = record.get("patient id", "")
            author_included = record.get("included in analysis", "").lower() == "yes"
            technical = bool(re.search(r"_ez|rehyb", record["title"], re.IGNORECASE))
            reasons = [
                name
                for name, excluded in (
                    ("outside_endpoint", not label),
                    ("geo_author_not_included", not author_included),
                    ("technical_replicate_indicator", technical),
                    ("cel_missing_from_tracked_archive_inventory", cel not in archive_cel),
                    ("missing_patient_id", not donor),
                )
                if excluded
            ]
            exclusion_counts.update(reasons)
            if label and technical:
                exclusion_counts["technical_indicator_in_endpoint_tissues"] += 1
            eligible = not reasons
            expected.update(
                author_included=str(author_included),
                technical_replicate_candidate=str(technical),
            )
            scan_date = record.get("scan date", "")
            if update_scan_dates:
                stored["scan_date"] = scan_date
            expected["scan_date"] = scan_date
            if clinical is not None:
                expected.update({
                    field: clinical.get(donor, {}).get(field, "")
                    for field in CLINICAL_FIELDS.values()
                })
        expected["original_tissue"] = tissue
        expected["patient_id"] = donor
        expected["inclusion_status"] = "included" if eligible else "excluded"
        if label is not None:
            expected["tissue_class"] = label
        for field, value in expected.items():
            if stored[field] != value:
                disagreements.append(
                    dict(accession=accession, sample_id=gsm, field=field, geo=value, stored=stored[field])
                )
        source_labels[tissue] += 1
        if eligible:
            labels[label] += 1
            donor_tissue[donor].add(label)

    pair = {"tumor", "adjacent_normal"} if primary else {"primary_tumor", "normal_colon"}
    counts = {
        "total_arrays": len(records),
        "eligible_arrays": sum(labels.values()),
        "excluded_arrays": len(records) - sum(labels.values()),
        "eligible_donor_or_patient_groups": len(donor_tissue),
        "groups_with_both_endpoint_tissues": sum(pair <= labels for labels in donor_tissue.values()),
        "compared_fields_with_disagreements": len(disagreements),
        **{f"eligible_{k}": v for k, v in sorted(labels.items())},
        **{f"source_tissue_{k}": v for k, v in sorted(source_labels.items())},
        **{f"exclusion_{k}": v for k, v in sorted(exclusion_counts.items())},
    }
    missingness = []
    for label in sorted(labels):
        stratum = [
            row for row in metadata
            if row["tissue_class"] == label and row["inclusion_status"] == "included"
        ]
        for field in (
            "age", "sex", "location", "stage", "msi_status",
            "scan_date", "processing_date", "center", "batch",
        ):
            missing = sum(row[field].strip().lower() in {"", "--", "na", "n/a"} for row in stratum)
            missingness.append(
                dict(accession=accession, tissue=label, field=field, arrays=len(stratum), missing=missing)
            )
    if update_scan_dates and not primary and not disagreements:
        write_csv(metadata_path, metadata)
    return (
        [dict(accession=accession, measure=k, value=v) for k, v in counts.items()],
        missingness,
        disagreements,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse44076", required=True, type=Path)
    parser.add_argument("--gse41258", required=True, type=Path)
    parser.add_argument("--gse41258-clinical", type=Path)
    parser.add_argument("--update-external-scan-dates", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT / "results/ai_audit")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summary, missingness, disagreements, sources = [], [], [], []
    for accession, source in (("GSE44076", args.gse44076), ("GSE41258", args.gse41258)):
        clinical = (
            clinical_values(args.gse41258_clinical)
            if accession == "GSE41258" and args.gse41258_clinical
            else None
        )
        counts, missing, differences = reconcile(
            ROOT, accession, source, args.update_external_scan_dates, clinical
        )
        summary.extend(counts)
        missingness.extend(missing)
        disagreements.extend(differences)
        sources.append(
            dict(
                accession=accession,
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=all&form=text&view=brief",
                bytes=source.stat().st_size,
                sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            )
        )
    if args.gse41258_clinical:
        source = args.gse41258_clinical
        sources.append(dict(
            accession="GSE41258_clinical",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41258/suppl/GSE41258_clinical_data.txt.gz",
            bytes=source.stat().st_size,
            sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        ))
    write_csv(args.output / "geo_cohort_reconciliation.csv", summary)
    write_csv(args.output / "geo_covariate_missingness.csv", missingness)
    write_csv(args.output / "geo_metadata_source_hashes.csv", sources)
    if disagreements:
        write_csv(args.output / "geo_metadata_disagreements.csv", disagreements)
        raise SystemExit(f"{len(disagreements)} sample fields disagree. Inspect the reconciliation CSV.")
    (args.output / "geo_metadata_disagreements.csv").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
