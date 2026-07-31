from __future__ import annotations

import numpy as np

from src.modeling.pipelines import TrainingVarianceQuantile, build_pipeline


def test_variance_filter_is_fit_on_training_only():
    train = np.array([[0.0, 0.0], [1.0, 0.1], [2.0, 0.2]])
    outer_test = np.array([[0.0, -1000.0], [0.0, 1000.0]])
    selector = TrainingVarianceQuantile(quantile=0.75).fit(train)
    expected = np.quantile(np.var(train, axis=0), 0.75)
    assert selector.threshold_ == expected
    assert selector.get_support().tolist() == [True, False]
    selector.transform(outer_test)
    assert selector.threshold_ == expected


def test_model_training_reproducible_with_fixed_seed():
    rng = np.random.default_rng(7)
    x = rng.normal(size=(40, 20))
    y = np.array(["a"] * 20 + ["b"] * 20)
    parameters = {
        "variance_quantile__quantile": 0.0,
        "univariate__k": 10,
        "model__C": 0.2,
        "model__l1_ratio": 0.5,
    }
    first = build_pipeline("elastic_net", 99).set_params(**parameters).fit(x, y)
    second = build_pipeline("elastic_net", 99).set_params(**parameters).fit(x, y)
    np.testing.assert_array_equal(first.predict(x), second.predict(x))
    np.testing.assert_allclose(
        first.named_steps["model"].coef_, second.named_steps["model"].coef_
    )
