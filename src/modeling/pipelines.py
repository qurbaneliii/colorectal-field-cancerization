from __future__ import annotations

import numpy as np
from joblib import Memory
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


class TrainingVarianceQuantile(BaseEstimator, TransformerMixin):
    """Select features by a variance quantile learned only in fit()."""

    def __init__(self, quantile: float = 0.5):
        self.quantile = quantile

    def fit(self, x, y=None):
        values = np.asarray(x)
        self.variances_ = np.var(values, axis=0)
        self.threshold_ = float(np.quantile(self.variances_, self.quantile))
        self.support_ = self.variances_ >= self.threshold_
        if not self.support_.any():
            self.support_[int(np.argmax(self.variances_))] = True
        self.n_features_in_ = values.shape[1]
        return self

    def transform(self, x):
        return np.asarray(x)[:, self.support_]

    def get_support(self, indices: bool = False):
        if indices:
            return np.flatnonzero(self.support_)
        return self.support_.copy()


def build_pipeline(model_name: str, seed: int) -> Pipeline:
    if model_name == "elastic_net":
        estimator = LogisticRegression(
            solver="saga",
            l1_ratio=0.5,
            max_iter=2000,
            class_weight="balanced",
            random_state=seed,
            tol=1e-3,
        )
    elif model_name == "linear_svm":
        estimator = LinearSVC(
            class_weight="balanced",
            dual="auto",
            max_iter=10000,
            random_state=seed,
        )
    elif model_name == "random_forest":
        estimator = RandomForestClassifier(
            n_estimators=200,
            class_weight="balanced",
            random_state=seed,
            n_jobs=1,
        )
    else:
        raise ValueError(f"Unknown model: {model_name}")
    memory = Memory(location="data/interim/sklearn_pipeline_cache", verbose=0)
    return Pipeline(
        [
            ("zero_variance", VarianceThreshold()),
            ("variance_quantile", TrainingVarianceQuantile()),
            ("univariate", SelectKBest(score_func=f_classif)),
            ("scale", StandardScaler()),
            ("model", estimator),
        ],
        memory=memory,
    )


def parameter_grid(model_name: str) -> dict[str, list[object]]:
    common = {
        "variance_quantile__quantile": [0.5],
        "univariate__k": [20, 50],
    }
    if model_name == "elastic_net":
        return {
            **common,
            "model__C": [0.05, 0.2],
            "model__l1_ratio": [0.2, 0.8],
        }
    if model_name == "linear_svm":
        return {
            **common,
            "univariate__k": [50],
            "model__C": [0.05, 0.2],
        }
    if model_name == "random_forest":
        return {
            **common,
            "univariate__k": [50],
            "model__max_features": ["sqrt"],
            "model__min_samples_leaf": [1, 3],
        }
    raise ValueError(model_name)


def selected_feature_names(pipeline: Pipeline, feature_names: list[str]) -> list[str]:
    names = np.asarray(feature_names, dtype=object)
    names = names[pipeline.named_steps["zero_variance"].get_support()]
    names = names[pipeline.named_steps["variance_quantile"].get_support()]
    names = names[pipeline.named_steps["univariate"].get_support()]
    return names.astype(str).tolist()
