from __future__ import annotations

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
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


class GroupCalibratedLinearSVC(BaseEstimator):
    """Linear SVM with calibration splits that preserve patient/donor groups."""

    def __init__(self, C: float = 1.0, calibration_splits: int = 3, random_state: int = 0):
        self.C = C
        self.calibration_splits = calibration_splits
        self.random_state = random_state

    def fit(self, x, y, groups=None):
        if groups is None:
            raise ValueError("GroupCalibratedLinearSVC requires patient/donor groups")
        values = np.asarray(x)
        labels = np.asarray(y)
        groups = np.asarray(groups)
        if len(groups) != len(labels):
            raise ValueError("Calibration groups must match the training rows")
        per_class_groups = [np.unique(groups[labels == label]).size for label in np.unique(labels)]
        n_splits = min(int(self.calibration_splits), min(per_class_groups))
        if n_splits < 2:
            raise ValueError("At least two patient/donor groups per class are required for calibration")
        splitter = StratifiedGroupKFold(
            n_splits=n_splits, shuffle=True, random_state=self.random_state
        )
        splits = list(splitter.split(values, labels, groups))
        for train, test in splits:
            overlap = set(groups[train]).intersection(groups[test])
            if overlap:
                raise AssertionError(f"Calibration group leakage detected: {sorted(overlap)[:5]}")
        base = LinearSVC(
            C=self.C,
            class_weight="balanced",
            dual="auto",
            max_iter=10000,
            random_state=self.random_state,
        )
        self.calibrated_ = CalibratedClassifierCV(estimator=base, cv=splits, method="sigmoid")
        self.calibrated_.fit(values, labels)
        self.classes_ = self.calibrated_.classes_
        self.n_features_in_ = values.shape[1]
        return self

    def predict(self, x):
        return self.calibrated_.predict(x)

    def predict_proba(self, x):
        return self.calibrated_.predict_proba(x)


def build_pipeline(model_name: str, seed: int, modeling_config: dict) -> Pipeline:
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
        estimator = GroupCalibratedLinearSVC(
            calibration_splits=int(modeling_config["linear_svm"]["calibration_splits"]),
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
    return Pipeline(
        [
            ("zero_variance", VarianceThreshold()),
            ("variance_quantile", TrainingVarianceQuantile()),
            ("univariate", SelectKBest(score_func=f_classif)),
            ("scale", StandardScaler()),
            ("model", estimator),
        ],
        memory=None,
    )


def parameter_grid(model_name: str, modeling_config: dict) -> dict[str, list[object]]:
    common = {
        "variance_quantile__quantile": list(modeling_config["variance_quantiles"]),
        "univariate__k": list(modeling_config["max_features"]),
    }
    if model_name == "elastic_net":
        return {
            **common,
            "model__C": list(modeling_config["elastic_net"]["C"]),
            "model__l1_ratio": list(modeling_config["elastic_net"]["l1_ratio"]),
        }
    if model_name == "linear_svm":
        return {
            **common,
            "model__C": list(modeling_config["linear_svm"]["C"]),
        }
    if model_name == "random_forest":
        return {
            **common,
            "model__n_estimators": list(modeling_config["random_forest"]["n_estimators"]),
            "model__max_features": list(modeling_config["random_forest"]["max_features"]),
            "model__min_samples_leaf": list(modeling_config["random_forest"]["min_samples_leaf"]),
        }
    raise ValueError(model_name)


def selected_feature_names(pipeline: Pipeline, feature_names: list[str]) -> list[str]:
    names = np.asarray(feature_names, dtype=object)
    names = names[pipeline.named_steps["zero_variance"].get_support()]
    names = names[pipeline.named_steps["variance_quantile"].get_support()]
    names = names[pipeline.named_steps["univariate"].get_support()]
    return names.astype(str).tolist()
