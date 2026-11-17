"""Explicit sample eligibility for a prospective acquisition-matched Task B study.

Schema validation is not legal authorization. Real records stay in approved compute.
"""
import numpy as np
import pandas as pd

REQUIRED = {"sample_id", "donor_id", "source_record", "source_tissue", "cohort",
            "histology", "site", "acquisition", "pretreatment", "tumor_distance_cm",
            "qc_pass", "consent_compatible", "biopsy_order", "age", "sex", "batch",
            "platform", "overlap_review"}
ENUMS = {
    "cohort": {"healthy", "crc", "adenoma", "unknown"},
    "histology": {"normal_confirmed", "tumor", "adenoma", "unknown"},
    "acquisition": {"endoscopic_biopsy", "surgical", "unknown"},
    "pretreatment": {"yes", "no", "unknown"},
    "qc_pass": {"yes", "no", "unknown"},
    "consent_compatible": {"yes", "no", "unknown"},
    "overlap_review": {"independent", "overlap", "unknown"},
    "site": {"cecum", "ascending", "transverse", "descending", "sigmoid", "rectum", "ileum", "unknown"},
}


def screen(metadata, *, site="sigmoid", minimum_distance_cm=10.0, development_donors=()):
    if REQUIRED.difference(metadata):
        raise ValueError(f"Missing fields: {sorted(REQUIRED.difference(metadata))}")
    if site not in ENUMS["site"] - {"ileum", "unknown"} or not np.isfinite(minimum_distance_cm) or minimum_distance_cm < 0:
        raise ValueError("Invalid prespecified site or distance")
    data = metadata.copy().reset_index(drop=True).fillna("").astype(str)
    for name in ["sample_id", "donor_id", "source_record", "source_tissue", "batch", "platform"]:
        if data[name].str.strip().eq("").any():
            raise ValueError(f"Explicit {name} required")
    if data.sample_id.duplicated().any():
        raise ValueError("Duplicate sample identity")
    for name, choices in ENUMS.items():
        if not set(data[name]).issubset(choices):
            raise ValueError(f"Unknown {name} values must use the explicit unknown category")
    if data.groupby("donor_id").cohort.nunique().gt(1).any():
        raise ValueError("Conflicting donor cohort labels")
    distance = pd.to_numeric(data.tumor_distance_cm, errors="coerce")
    age = pd.to_numeric(data.age, errors="coerce")
    order = pd.to_numeric(data.biopsy_order, errors="coerce")
    if (order.isna().any() or (~np.isfinite(order)).any() or (order < 1).any()
            or (order % 1 != 0).any() or data.assign(order=order).duplicated(["donor_id", "order"]).any()):
        raise ValueError("Unique positive integer source biopsy order is required per donor")
    reasons = [[] for _ in range(len(data))]
    rules = {
        "endpoint_not_eligible": ~data.cohort.isin(["healthy", "crc"]),
        "histology_not_confirmed": ~data.histology.eq("normal_confirmed"),
        "site_not_primary": ~data.site.eq(site),
        "acquisition_not_matched": ~data.acquisition.eq("endoscopic_biopsy"),
        "pretreatment_not_confirmed": ~data.pretreatment.eq("yes"),
        "distance_not_remote_or_unknown": data.cohort.eq("crc") & (~np.isfinite(distance) | (distance < minimum_distance_cm)),
        "qc_not_passed": ~data.qc_pass.eq("yes"),
        "consent_not_confirmed": ~data.consent_compatible.eq("yes"),
        "independence_not_confirmed": ~data.overlap_review.eq("independent") | data.donor_id.isin(development_donors),
        "age_missing_or_invalid": ~np.isfinite(age) | (age < 18) | (age > 120),
        "sex_missing": ~data.sex.isin(["female", "male", "other"]),
    }
    for reason, excluded in rules.items():
        for i in np.flatnonzero(excluded.to_numpy()):
            reasons[i].append(reason)
    data["exclusion_reason"] = [";".join(row) for row in reasons]
    data["_order"] = order
    candidates = data.loc[data.exclusion_reason.eq("")].sort_values(["donor_id", "_order"])
    duplicates = candidates.index[candidates.donor_id.duplicated()]
    data.loc[duplicates, "exclusion_reason"] = "secondary_biopsy_same_donor"
    data["eligible"] = data.exclusion_reason.eq("")
    data["tissue"] = data.cohort.map({"healthy": "healthy", "crc": "cancer_associated_normal"}).fillna("excluded")
    return data.drop(columns="_order")
