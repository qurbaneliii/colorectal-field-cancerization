from __future__ import annotations

import pandas as pd

from src.data.metadata import discover_accession_files
from src.data.validation import validate_metadata


def test_expected_data_file_discovery(root):
    for accession in ("GSE44076", "GSE41258"):
        found = discover_accession_files(root / "data/raw", accession)
        assert len(found["raw_tar"]) == 1
        assert len(found["series_matrix"]) == 1
    assert len(discover_accession_files(root / "data/raw", "GSE41258")["clinical"]) == 1


def test_primary_metadata_labels_and_ids(root):
    frame = pd.read_csv(root / "data/metadata/gse44076_samples.csv", dtype={"patient_id": str})
    validate_metadata(frame, {"healthy", "adjacent_normal", "tumor"})
    assert frame["tissue_class"].value_counts().to_dict() == {
        "adjacent_normal": 98,
        "tumor": 98,
        "healthy": 50,
    }
    assert frame["patient_id"].nunique() == 148


def test_tumor_adjacent_pairing(root):
    pairs = pd.read_csv(root / "data/metadata/gse44076_pairing_audit.csv")
    assert len(pairs) == 98
    assert pairs["both_present"].all()
    assert not pairs["duplicate_samples"].any()


def test_external_exclusions(root):
    frame = pd.read_csv(root / "data/metadata/gse41258_samples.csv", dtype={"patient_id": str})
    included = frame[frame["inclusion_status"].eq("included")]
    assert set(included["tissue_class"]) == {"primary_tumor", "normal_colon"}
    forbidden = {
        "liver_metastasis",
        "lung_metastasis",
        "normal_liver",
        "normal_lung",
        "cell_line",
    }
    assert forbidden.isdisjoint(included["tissue_class"])
    assert not included["technical_replicate_candidate"].any()
    assert included["author_included"].all()
