import pandas as pd
import pytest

from src.data.independent_eligibility import screen


def metadata():
    return pd.DataFrame([{ "sample_id": f"SYNTHETIC_{i}", "donor_id": f"donor{i}",
        "source_record": "synthetic fixture", "source_tissue": "normal mucosa",
        "cohort": "healthy" if i < 2 else "crc", "histology": "normal_confirmed",
        "site": "sigmoid", "acquisition": "endoscopic_biopsy", "pretreatment": "yes",
        "tumor_distance_cm": "" if i < 2 else "12", "qc_pass": "yes",
        "consent_compatible": "yes", "biopsy_order": "1", "age": "60", "sex": "female",
        "batch": "fixture", "platform": "synthetic", "overlap_review": "independent"}
        for i in range(4)])


def test_explicit_distance_and_donor_rules():
    data = metadata()
    assert screen(data).eligible.all()
    data.loc[2, "tumor_distance_cm"] = ""
    data.loc[3, "tumor_distance_cm"] = "2"
    result = screen(data)
    assert result.eligible.tolist() == [True, True, False, False]
    assert screen(metadata(), development_donors=["donor0"]).eligible.tolist() == [False, True, True, True]
    data = metadata()
    second = data.iloc[[0]].copy()
    second["sample_id"], second["biopsy_order"] = "SYNTHETIC_repeat", "2"
    result = screen(pd.concat([data, second], ignore_index=True))
    assert result.eligible.sum() == 4
    assert result.iloc[-1].exclusion_reason == "secondary_biopsy_same_donor"
    second["cohort"] = "crc"
    with pytest.raises(ValueError, match="Conflicting"):
        screen(pd.concat([data, second], ignore_index=True))


@pytest.mark.parametrize("field,value", [("site", "ileum"), ("qc_pass", "unknown"),
    ("histology", "unknown"), ("acquisition", "surgical"), ("pretreatment", "no"),
    ("age", "nan"), ("consent_compatible", "unknown"), ("overlap_review", "unknown")])
def test_uncertain_metadata_never_silently_eligible(field, value):
    data = metadata()
    data.loc[0, field] = value
    assert not screen(data).iloc[0].eligible
