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


def test_reference_and_mapping_validation():
    from dataclasses import FrozenInstanceError
    ref = FrozenQuantile((1., 2., 4.), ("a",))
    with pytest.raises(FrozenInstanceError):
        ref.target = (2., 3., 4.)
    with pytest.raises(ValueError, match="Invalid"):
        FrozenQuantile.fit({"a": [1., np.nan, 4.]}, ["a"])
    with pytest.raises(ValueError, match="Invalid"):
        FrozenQuantile((2., 1.), ("a",))
    with pytest.raises(ValueError, match="unique"):
        FrozenQuantile.fit({"a": [1., 2.]}, ["a", "a"])
    mapping = pd.DataFrame({"probe_id": ["p1", "p2"], "gene_symbol": ["B", "A"],
                            "mapping_status": ["mapped_unique"] * 2})
    ranks = pd.DataFrame({"b": [.1, .9], "a": [.3, .6]}, index=["p1", "p2"])
    x = ref.transform(ranks, mapping)
    assert list(x.columns) == ["A", "B"]
    pd.testing.assert_frame_equal(x.loc[["a", "b"]], ref.transform(ranks[["a", "b"]], mapping))


def test_numeric_model_roundtrip_and_probability_order(tmp_path):
    import json
    from scipy.special import expit
    from experiments.frozen_reference.run_comparison import fit_predict
    x = pd.DataFrame({"A": [-2., -1., 1., 2.], "B": [0., 1., 1., 0.]})
    path = tmp_path / "model.json"
    model, scores = fit_predict(x, np.array([False, False, True, True]), x, path)
    data = json.loads(path.read_text())
    reconstructed = expit(((x[data["features"]] - data["mean"]) / data["scale"])
                           @ np.asarray(data["coefficients"]) + data["intercept"])
    np.testing.assert_allclose(scores, reconstructed, atol=1e-12)
    assert list(model.classes_) == [False, True]
    assert scores[0] < .5 < scores[-1]


def test_paired_cluster_bootstrap_identical_scores_have_zero_delta():
    from experiments.frozen_reference.run_comparison import paired_bootstrap
    rows = [{"sample_id": f"s{i}", "patient_id": f"d{i//2}", "observed": bool(i % 2),
             "method": method, "score": .8 if i % 2 else .2}
            for i in range(8) for method in ["joint_rma", "frozen"]]
    for result in paired_bootstrap(pd.DataFrame(rows), .5, iterations=20):
        assert result["delta_lower"] == result["delta_upper"] == 0
        assert result["valid_draws"] == 20
