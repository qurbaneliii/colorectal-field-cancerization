import numpy as np
import pandas as pd
import pytest

from src.preprocessing.frozen_quantile import FrozenQuantile


def test_reference_train_only_and_single_array_invariance(tmp_path):
    grids = {"train": [1, 2, 4], "test": [100, 200, 400]}
    ref = FrozenQuantile.fit(grids, ["train"], ["test"])
    mapping = pd.DataFrame({"probe_id": ["p1", "p2"], "gene_symbol": ["A", "B"],
                            "mapping_status": ["mapped_unique"] * 2})
    ranks = pd.DataFrame({"test": [.1, .9], "other": [.5, .7]}, index=["p1", "p2"])
    alone = ref.transform(ranks[["test"]], mapping)
    batch = ref.transform(ranks, mapping).loc[["test"]]
    pd.testing.assert_frame_equal(alone, batch)
    grids["test"] = [1e5, 2e5, 3e5]
    assert FrozenQuantile.fit(grids, ["train"]).digest() == ref.digest()
    path = tmp_path / "ref.json"
    ref.save(path)
    assert FrozenQuantile.load(path).digest() == ref.digest()
    np.testing.assert_allclose(ref.probe_values([0, .5, 1]), [0, 1, 2])
    with pytest.raises(ValueError, match="disjoint"):
        FrozenQuantile.fit(grids, ["train"], ["train"])
    with pytest.raises(ValueError, match="missing"):
        ref.transform(ranks.loc[["p1"]], mapping)
    with pytest.raises(ValueError, match="percentile"):
        ref.probe_values([np.nan])
