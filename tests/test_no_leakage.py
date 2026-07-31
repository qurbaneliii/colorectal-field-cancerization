from __future__ import annotations

import numpy as np
import yaml

from src.modeling.pipelines import (
    GroupCalibratedLinearSVC,
    TrainingVarianceQuantile,
    build_pipeline,
    parameter_grid,
)


def test_variance_filter_is_fit_on_training_only():
    train = np.array([[0.0, 0.0], [1.0, 0.1], [2.0, 0.2]])
    outer_test = np.array([[0.0, -1000.0], [0.0, 1000.0]])
    selector = TrainingVarianceQuantile(quantile=0.75).fit(train)
    expected = np.quantile(np.var(train, axis=0), 0.75)
    assert selector.threshold_ == expected
    assert selector.get_support().tolist() == [True, False]
    selector.transform(outer_test)
    assert selector.threshold_ == expected


def test_model_training_reproducible_with_fixed_seed(root):
    rng = np.random.default_rng(7)
    x = rng.normal(size=(40, 20))
    y = np.array(["a"] * 20 + ["b"] * 20)
    parameters = {
        "variance_quantile__quantile": 0.0,
        "univariate__k": 10,
        "model__C": 0.2,
        "model__l1_ratio": 0.5,
    }
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    first = build_pipeline("elastic_net", 99, config["modeling"]).set_params(**parameters).fit(x, y)
    second = build_pipeline("elastic_net", 99, config["modeling"]).set_params(**parameters).fit(x, y)
    np.testing.assert_array_equal(first.predict(x), second.predict(x))
    np.testing.assert_allclose(
        first.named_steps["model"].coef_, second.named_steps["model"].coef_
    )


def test_runtime_parameter_grids_match_configuration(root):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    modeling = config["modeling"]
    elastic = parameter_grid("elastic_net", modeling)
    assert elastic["variance_quantile__quantile"] == modeling["variance_quantiles"]
    assert elastic["univariate__k"] == modeling["max_features"]
    assert elastic["model__C"] == modeling["elastic_net"]["C"]
    assert elastic["model__l1_ratio"] == modeling["elastic_net"]["l1_ratio"]


def test_svm_calibration_preserves_groups():
    rng = np.random.default_rng(17)
    groups = np.repeat(np.arange(12), 2)
    y = np.tile(np.array(["adjacent_normal", "tumor"]), 12)
    x = rng.normal(size=(24, 8)) + (y == "tumor")[:, None]
    model = GroupCalibratedLinearSVC(C=0.2, calibration_splits=3, random_state=7)
    model.fit(x, y, groups=groups)
    probability = model.predict_proba(x)
    assert probability.shape == (24, 2)
    np.testing.assert_allclose(probability.sum(axis=1), 1.0)
