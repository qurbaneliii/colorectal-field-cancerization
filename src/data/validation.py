from __future__ import annotations

import numpy as np
import pandas as pd


def validate_metadata(frame: pd.DataFrame, allowed_labels: set[str]) -> None:
    required = {
        "geo_accession",
        "tissue_class",
        "patient_id",
        "platform",
        "inclusion_status",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise AssertionError(f"Missing metadata columns: {sorted(missing)}")
    if frame["geo_accession"].duplicated().any():
        raise AssertionError("Duplicate sample identifiers")
    included = frame["inclusion_status"].eq("included")
    if frame.loc[included, "tissue_class"].isna().any():
        raise AssertionError("Included samples contain missing target labels")
    invalid = set(frame.loc[included, "tissue_class"]).difference(allowed_labels)
    if invalid:
        raise AssertionError(f"Unexpected included labels: {sorted(invalid)}")
    if frame.loc[included, "patient_id"].astype(str).str.strip().eq("").any():
        raise AssertionError("Included samples contain missing patient/donor IDs")


def validate_expression(
    expression: pd.DataFrame, metadata: pd.DataFrame, gene_column: str = "gene_symbol"
) -> None:
    if gene_column not in expression:
        raise AssertionError(f"Missing {gene_column} column")
    if expression[gene_column].duplicated().any():
        raise AssertionError("Gene identifiers are not unique after aggregation")
    sample_columns = expression.columns.drop(gene_column)
    if sample_columns.duplicated().any():
        raise AssertionError("Duplicate expression sample columns")
    if set(sample_columns) != set(metadata["geo_accession"]):
        missing_expr = set(metadata["geo_accession"]).difference(sample_columns)
        missing_meta = set(sample_columns).difference(metadata["geo_accession"])
        raise AssertionError(
            f"Expression/metadata mismatch: missing expression={sorted(missing_expr)[:5]}, "
            f"missing metadata={sorted(missing_meta)[:5]}"
        )
    values = expression[sample_columns].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise AssertionError("Expression matrix contains missing or infinite values")
    if values.shape[0] <= values.shape[1]:
        raise AssertionError("Expected genes as rows and samples as columns")


def assert_no_group_overlap(train_index: np.ndarray, test_index: np.ndarray, groups: np.ndarray) -> None:
    overlap = set(groups[train_index]).intersection(groups[test_index])
    if overlap:
        raise AssertionError(f"Patient/donor leakage detected: {sorted(overlap)[:5]}")


def validate_fold_assignments(assignments: pd.DataFrame) -> None:
    required = {"task", "repeat", "outer_fold", "split", "patient_id"}
    if not required.issubset(assignments):
        raise AssertionError(f"Fold assignments missing {sorted(required.difference(assignments))}")
    for key, fold in assignments.groupby(["task", "repeat", "outer_fold"]):
        train = set(fold.loc[fold["split"].eq("train"), "patient_id"])
        test = set(fold.loc[fold["split"].eq("test"), "patient_id"])
        overlap = train.intersection(test)
        if overlap:
            raise AssertionError(f"Group leakage in {key}: {sorted(overlap)[:5]}")
