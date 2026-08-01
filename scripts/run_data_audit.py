from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import numpy as np
import sklearn
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.metadata import (
    build_gse41258_metadata,
    build_gse44076_metadata,
    discover_accession_files,
    tar_members,
)
from src.data.validation import validate_metadata


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_one(discovery: dict[str, list[Path]], key: str, accession: str) -> Path:
    observed = discovery[key]
    if len(observed) != 1:
        names = [str(p) for p in observed]
        raise FileNotFoundError(f"{accession}: expected one {key}, observed {len(observed)}: {names}")
    return observed[0]


def command_version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=20)
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return f"unavailable: {exc}"
    text = (result.stdout or result.stderr).strip().replace("\n", " | ")
    return text or f"exit {result.returncode}"


def discover_rscript() -> str | None:
    configured = os.environ.get("RSCRIPT")
    candidates = [
        configured,
        shutil.which("Rscript"),
        r"D:\R\R-4.5.1\bin\Rscript.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    return None


def main() -> None:
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    for path in paths.values():
        (ROOT / path).mkdir(parents=True, exist_ok=True)

    raw_root = ROOT / paths["raw"]
    primary_files = discover_accession_files(raw_root, "GSE44076")
    external_files = discover_accession_files(raw_root, "GSE41258")
    p_matrix = require_one(primary_files, "series_matrix", "GSE44076")
    p_tar = require_one(primary_files, "raw_tar", "GSE44076")
    e_matrix = require_one(external_files, "series_matrix", "GSE41258")
    e_tar = require_one(external_files, "raw_tar", "GSE41258")
    e_clinical = require_one(external_files, "clinical", "GSE41258")

    processed_at = datetime.now(timezone.utc).isoformat()
    manifest_rows = []
    for accession, files in (("GSE44076", primary_files), ("GSE41258", external_files)):
        for candidates in files.values():
            for path in candidates:
                manifest_rows.append(
                    {
                        "accession": accession,
                        "filename": path.name,
                        "relative_path": path.relative_to(ROOT).as_posix(),
                        "size_bytes": path.stat().st_size,
                        "sha256": sha256(path),
                        "archive_member_count": (
                            len(tar_members(path)) if path.name.lower().endswith("_raw.tar") else ""
                        ),
                        "processed_at_utc": processed_at,
                        "software_environment": (
                            f"Python {platform.python_version()}; {platform.platform()}"
                        ),
                        "data_provenance": "NCBI GEO deposited file",
                    }
                )
    manifest = pd.DataFrame(manifest_rows).sort_values(["accession", "filename"])
    manifest.to_csv(ROOT / paths["metadata"] / "data_manifest.csv", index=False)

    primary, pairs = build_gse44076_metadata(p_matrix, p_tar)
    external = build_gse41258_metadata(
        e_matrix,
        e_tar,
        config["metadata"]["technical_replicate_pattern"],
        e_clinical,
    )
    validate_metadata(primary, {"healthy", "adjacent_normal", "tumor"})
    validate_metadata(external, {"primary_tumor", "normal_colon"})
    primary.to_csv(ROOT / paths["metadata"] / "gse44076_samples.csv", index=False)
    external.to_csv(ROOT / paths["metadata"] / "gse41258_samples.csv", index=False)
    pairs.to_csv(ROOT / paths["metadata"] / "gse44076_pairing_audit.csv", index=False)

    dictionary = pd.DataFrame(
        [
            ["geo_accession", "GEO GSM identifier", "GEO series matrix", "string"],
            ["sample_title", "Original GEO sample title", "GEO series matrix", "string"],
            ["original_tissue", "Unmodified deposited tissue category", "GEO metadata", "string"],
            ["tissue_class", "Canonical analysis label", "explicit mapping rule", "category"],
            ["patient_id", "Explicit deposited individual/patient ID", "GEO characteristic", "string"],
            ["donor_or_patient_group", "Group used for every CV split", "patient_id", "string"],
            ["pairing_status", "Tumor/adjacent pair availability", "pair audit", "category"],
            ["technical_replicate_candidate", "Title matches _ez or rehyb", "deterministic regex", "bool"],
            ["inclusion_status", "Eligibility for the primary analysis of that accession", "rules", "category"],
            ["exclusion_reason", "Semicolon-separated deterministic reasons", "rules", "string"],
        ],
        columns=["field", "definition", "derivation", "type"],
    )
    dictionary.to_csv(ROOT / paths["metadata"] / "metadata_dictionary.csv", index=False)

    covariates = [
        "batch",
        "processing_date",
        "scan_date",
        "sex",
        "age",
        "location",
        "stage",
        "msi_status",
        "molecular_subtype",
        "center",
    ]
    missingness_rows = []
    association_rows = []
    for accession, frame in (("GSE44076", primary), ("GSE41258", external)):
        included = frame[frame["inclusion_status"].eq("included")].copy()
        for covariate in covariates:
            values = (
                included[covariate].fillna("").astype(str).str.strip()
                if covariate in included
                else pd.Series([""] * len(included), index=included.index)
            )
            available = values.ne("") & values.str.lower().ne("nan")
            unique_values = values[available].nunique()
            missingness_rows.append(
                {
                    "accession": accession,
                    "covariate": covariate,
                    "included_samples": len(included),
                    "available_values": int(available.sum()),
                    "missing_values": int((~available).sum()),
                    "unique_nonmissing_values": int(unique_values),
                    "usable_for_association": bool(available.sum() >= 10 and unique_values >= 2),
                    "source": "deposited metadata only",
                }
            )
            if available.sum() >= 10 and unique_values >= 2:
                numeric = pd.to_numeric(values[available], errors="coerce")
                if numeric.notna().mean() >= 0.9:
                    groups = [
                        numeric[included.loc[available, "tissue_class"].eq(label)]
                        for label in sorted(included.loc[available, "tissue_class"].unique())
                    ]
                    groups = [group.dropna().to_numpy() for group in groups if len(group.dropna())]
                    statistic, p_value = stats.kruskal(*groups) if len(groups) >= 2 else (np.nan, np.nan)
                    method = "Kruskal-Wallis"
                else:
                    contingency = pd.crosstab(
                        included.loc[available, "tissue_class"], values[available]
                    )
                    statistic, p_value, _, _ = stats.chi2_contingency(contingency)
                    method = "chi-square"
                association_rows.append(
                    {
                        "accession": accession,
                        "covariate": covariate,
                        "method": method,
                        "statistic": statistic,
                        "p_value": p_value,
                        "samples": int(available.sum()),
                        "interpretation": "screening association only; no automatic batch correction",
                    }
                )
    missingness = pd.DataFrame(missingness_rows)
    associations = pd.DataFrame(association_rows)
    missingness.to_csv(ROOT / paths["metadata"] / "covariate_missingness_usability.csv", index=False)
    associations.to_csv(ROOT / paths["metadata"] / "covariate_tissue_associations.csv", index=False)
    confounding_report = f"""# Batch and confounding audit

Only deposited covariates were used; absent fields were not invented. Screening
associations do not trigger automatic ComBat or sample exclusion.

## Missingness and usability

{missingness.to_markdown(index=False)}

## Tissue-label association screens

{associations.to_markdown(index=False) if not associations.empty else 'No covariate met the minimum usability rule.'}

Any estimable covariate must be incorporated through a prespecified inferential
design or learned inside training folds. GSE41258 labels must never be used to
harmonize the external cohort with GSE44076.
"""
    (ROOT / paths["reports"] / "batch_confounding_audit.md").write_text(
        confounding_report, encoding="utf-8"
    )

    rscript = discover_rscript()
    r_status = command_version([rscript, "--version"]) if rscript else "unavailable"
    env = {
        "processed_at_utc": processed_at,
        "python": sys.version,
        "platform": platform.platform(),
        "scikit_learn": sklearn.__version__,
        "Rscript": r_status,
        "Rscript_path": rscript,
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False
        ).stdout.strip(),
        "raw_data_immutable_policy": True,
    }
    (ROOT / paths["metadata"] / "software_environment.json").write_text(
        json.dumps(env, indent=2), encoding="utf-8"
    )

    expected = config["expected"]["GSE44076"]
    primary_counts = primary["tissue_class"].value_counts().to_dict()
    observed_pairs = int(pairs["both_present"].sum())
    discrepancies = []
    for key in ("healthy", "adjacent_normal", "tumor"):
        if primary_counts.get(key, 0) != expected[key]:
            discrepancies.append(f"{key}: expected {expected[key]}, observed {primary_counts.get(key, 0)}")
    if len(primary) != expected["total"]:
        discrepancies.append(f"total: expected {expected['total']}, observed {len(primary)}")
    if observed_pairs != expected["paired_patients"]:
        discrepancies.append(
            f"paired patients: expected {expected['paired_patients']}, observed {observed_pairs}"
        )

    r_blocker = (
        "- R was discovered and the raw-CEL, Bioconductor QC, limma, composition, "
        "and enrichment stages have executable repository routes."
        if rscript
        else "- Rscript was not discovered on PATH, via RSCRIPT, or at the documented Windows path."
    )
    repository_audit = f"""# Repository audit

## Audited state

The audit preserves the existing hybrid R/Python implementation and reads the
five deposited GEO inputs without modifying `data/raw/`.

## Data discovered

{manifest[['accession', 'filename', 'size_bytes']].to_markdown(index=False)}

- GSE44076 raw TAR members: {len(tar_members(p_tar))}
- GSE41258 raw TAR members: {len(tar_members(e_tar))}
- Raw files were read but not modified.

## Runtime

- Python: {platform.python_version()}
- Rscript: {env['Rscript']}
- Git: {env['git_commit']}

## Missing or blocked requirements

{r_blocker}
- Docker verification is a separate host-level gate and is not inferred from
  the presence of a Dockerfile.

## Reusable components and modifications

The existing hybrid R/Python pipeline is reused and extended with deterministic
metadata parsing, provenance capture, independent platform annotation,
group-aware modeling, external validation, tests, and manuscript/report
generation.
"""
    (ROOT / paths["reports"] / "repository_audit.md").write_text(
        repository_audit, encoding="utf-8"
    )

    technical = external[external["technical_replicate_candidate"]]
    metadata_audit = f"""# Metadata audit

## GSE44076

Class counts:

{primary['tissue_class'].value_counts().rename_axis('class').reset_index(name='n').to_markdown(index=False)}

- Unique patient/donor IDs: {primary['patient_id'].nunique()}
- Complete adjacent-normal/tumor pairs: {observed_pairs}
- Duplicate pair entries: {int(pairs['duplicate_samples'].sum())}
- Missing CEL members: {int((~primary['cel_present_in_tar']).sum())}
- Unmatched or excluded samples: {int(primary['inclusion_status'].eq('excluded').sum())}
- Expected-count discrepancies: {('; '.join(discrepancies)) if discrepancies else 'none'}

## GSE41258

Preserved canonical category counts:

{external['tissue_class'].value_counts().rename_axis('class').reset_index(name='n').to_markdown(index=False)}

- Unique explicit patient IDs: {external.loc[external['patient_id'].ne(''), 'patient_id'].nunique()}
- Technical-replicate candidates: {len(technical)}
- Main external-validation included samples: {int(external['inclusion_status'].eq('included').sum())}
- Missing CEL members: {int((~external['cel_present_in_tar']).sum())}
- Excluded samples: {int(external['inclusion_status'].eq('excluded').sum())}

Technical-replicate candidates:

{technical[['geo_accession', 'sample_title', 'tissue_class', 'patient_id']].to_markdown(index=False)}

## Parsing rules

1. GSE44076 tissue labels are mapped only from explicit `source_name_ch1`;
   patient IDs come only from the `individual id` characteristic.
2. GSE41258 preserves the explicit `tissue` characteristic before canonical
   mapping. Only Primary Tumor and Normal Colon are eligible for main external
   validation.
3. GSE41258 samples must also carry the author flag `included in analysis: Yes`.
4. Titles matching `{config['metadata']['technical_replicate_pattern']}` are
   excluded deterministically as technical-replicate candidates.
5. Missing or unresolved metadata causes exclusion, never imputation.
"""
    (ROOT / paths["reports"] / "metadata_audit.md").write_text(
        metadata_audit, encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "GSE44076_counts": primary_counts,
                "GSE44076_pairs": observed_pairs,
                "GSE41258_included_external": int(
                    external["inclusion_status"].eq("included").sum()
                ),
                "count_discrepancies": discrepancies,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
