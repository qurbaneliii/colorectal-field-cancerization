from __future__ import annotations

import numpy as np
import pandas as pd

from src.data.validation import assert_no_group_overlap
from src.data.validation import validate_fold_assignments
from src.modeling.splitters import stratified_group_splits


def _primary_arrays(root):
    frame = pd.read_csv(root / "data/metadata/gse44076_samples.csv", dtype={"patient_id": str})
    return (
        frame["tissue_class"].to_numpy(),
        frame["donor_or_patient_group"].astype(str).to_numpy(),
    )


def test_no_patient_overlap_in_every_fold(root):
    y, groups = _primary_arrays(root)
    for train, test in stratified_group_splits(y, groups, 5, 44076):
        assert_no_group_overlap(train, test, groups)
        assert set(groups[train]).isdisjoint(groups[test])


def test_fold_generation_is_reproducible(root):
    y, groups = _primary_arrays(root)
    first = stratified_group_splits(y, groups, 5, 44076)
    second = stratified_group_splits(y, groups, 5, 44076)
    for (train_a, test_a), (train_b, test_b) in zip(first, second, strict=True):
        np.testing.assert_array_equal(train_a, train_b)
        np.testing.assert_array_equal(test_a, test_b)


def test_saved_nested_assignments_have_no_overlap(root):
    assignments = pd.read_csv(
        root / "results/metrics/fold_assignments.csv", dtype={"patient_id": str}
    )
    validate_fold_assignments(assignments)
