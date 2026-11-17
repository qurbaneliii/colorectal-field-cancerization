"""Training-only quantile reference for single-array median PM percentiles.

This research comparator omits learned probe-affinity effects and is not fRMA.
No labels are accepted. Mapping is a fixed, uniquely annotated manufacturer map.
"""
from dataclasses import dataclass
import hashlib
import json

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class FrozenQuantile:
    target: tuple[float, ...]
    training_ids: tuple[str, ...]

    def __post_init__(self):
        values = np.asarray(self.target)
        if (values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all()
                or (values <= 0).any() or (np.diff(values) < 0).any()
                or not self.training_ids or len(set(self.training_ids)) != len(self.training_ids)):
            raise ValueError("Invalid frozen quantile reference")

    @classmethod
    def fit(cls, grids, training_ids, held_out_ids=()):
        ids = tuple(training_ids)
        if not ids or len(set(ids)) != len(ids) or set(ids) & set(held_out_ids):
            raise ValueError("Training identities must be unique and disjoint from holdouts")
        values = np.asarray([grids[name] for name in ids], dtype=float)
        if (values.ndim != 2 or not np.isfinite(values).all() or (values <= 0).any()
                or (np.diff(values, axis=1) < 0).any()):
            raise ValueError("Invalid per-array training quantiles")
        return cls(tuple(values.mean(axis=0)), ids)

    def probe_values(self, percentiles):
        values = np.asarray(percentiles, dtype=float)
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError("Invalid per-array percentile")
        return np.log2(np.interp(values, np.linspace(0, 1, len(self.target)), self.target))

    def transform(self, ranks, mapping):
        if not ranks.index.is_unique or not ranks.columns.is_unique:
            raise ValueError("Duplicate probe or sample identity")
        valid = mapping.loc[mapping.mapping_status.eq("mapped_unique")].copy()
        if (valid.probe_id.duplicated().any() or valid.gene_symbol.isna().any()
                or valid.gene_symbol.eq("").any()):
            raise ValueError("Invalid unique gene mapping")
        if not set(valid.probe_id).issubset(ranks.index):
            raise ValueError("Mapped manufacturer probes missing from input")
        data = pd.DataFrame(self.probe_values(ranks.loc[valid.probe_id]),
                            index=valid.gene_symbol, columns=ranks.columns)
        return data.groupby(level=0, sort=True).median().T

    def payload(self):
        return {"method": "training_quantile_median_pm_percentile_v1",
                "target": list(self.target), "training_ids": list(self.training_ids)}

    def digest(self):
        return hashlib.sha256(json.dumps(self.payload(), sort_keys=True).encode()).hexdigest()

    def save(self, path):
        path.write_text(json.dumps(self.payload(), sort_keys=True) + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path):
        data = json.loads(path.read_text(encoding="utf-8"))
        if data["method"] != "training_quantile_median_pm_percentile_v1":
            raise ValueError("Unknown reference method")
        return cls(tuple(data["target"]), tuple(data["training_ids"]))
