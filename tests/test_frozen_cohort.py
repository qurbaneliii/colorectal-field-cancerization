import json

import pytest

from scripts.evaluate_frozen_cohort import evaluate
from src.data.provenance import file_sha256


def test_frozen_cohort_checks_lock_groups_and_outcome_alignment(tmp_path):
    lock, model, mapping, matrix, metadata = [tmp_path / name for name in
                                                ("lock.json", "model.json", "map.csv", "x.csv", "y.csv")]
    model.write_text(json.dumps({"features": ["A"], "coefficients": [1], "intercept": 0}))
    mapping.write_text("feature,column,annotation\nA,A.1,synthetic\n")
    matrix.write_text("sample_id,A.1\n001,-2\n002,-1\n003,1\n004,2\n")
    header = "sample_id,donor_id,tissue,site,acquisition,source_tissue,source_record\n"
    rows = ("004,d4,cancer_associated_normal,left,biopsy,cancer,fixture\n"
            "002,d2,healthy,left,biopsy,control,fixture\n"
            "003,d3,cancer_associated_normal,left,biopsy,cancer,fixture\n"
            "001,d1,healthy,left,biopsy,control,fixture\n")
    metadata.write_text(header + rows)
    lock.write_text(json.dumps({
        "mode": "synthetic", "training_accession": "training", "evaluation_accession": "evaluation",
        "model_sha256": file_sha256(model), "mapping_sha256": file_sha256(mapping),
        "normalization": "frozen_single_sample", "threshold": 0.5, "negative_tissue": "healthy",
        "positive_tissue": "cancer_associated_normal", "site": "left", "acquisition": "biopsy",
        "bootstrap_seed": 3, "bootstrap_iterations": 12,
    }))
    args = lock, model, mapping, matrix, metadata
    result = evaluate(*args, synthetic=True)
    assert result["mode"] == "synthetic"
    assert result["eligible_donors"] == 4
    assert result["metrics"]["tp"] == result["metrics"]["tn"] == 2
    assert result["conditional_donor_bootstrap_intervals"]["roc_auc"] == [1, 1]
    with pytest.raises(ValueError, match="synthetic mode"):
        evaluate(*args)
    metadata.write_text(header + rows.replace("d3", "d4"))
    with pytest.raises(ValueError, match="one prespecified specimen"):
        evaluate(*args, synthetic=True)
    mapping.write_text(mapping.read_text().replace("A.1", "changed"))
    with pytest.raises(ValueError, match="mapping digest"):
        evaluate(*args, synthetic=True)
