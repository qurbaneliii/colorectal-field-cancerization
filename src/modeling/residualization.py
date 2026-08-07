from __future__ import annotations

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class TrainingCompositionPCs:
    """Learn a limited composition-score representation from training rows only."""

    def __init__(self, n_components: int = 2):
        self.n_components = n_components

    def fit(self, composition_scores: np.ndarray):
        scores = np.asarray(composition_scores, dtype=float)
        if scores.ndim != 2:
            raise ValueError("Composition scores must be a two-dimensional matrix")
        if not 1 <= int(self.n_components) <= min(scores.shape):
            raise ValueError("Invalid number of composition principal components")
        self.scaler_ = StandardScaler().fit(scores)
        standardized = self.scaler_.transform(scores)
        self.pca_ = PCA(n_components=int(self.n_components)).fit(standardized)
        self.n_features_in_ = scores.shape[1]
        return self

    def transform(self, composition_scores: np.ndarray) -> np.ndarray:
        scores = np.asarray(composition_scores, dtype=float)
        if scores.ndim != 2 or scores.shape[1] != self.n_features_in_:
            raise ValueError("Composition-score columns do not match the training matrix")
        return self.pca_.transform(self.scaler_.transform(scores))

    def fit_transform(self, composition_scores: np.ndarray) -> np.ndarray:
        return self.fit(composition_scores).transform(composition_scores)


class TrainingCovariateResidualizer:
    """Remove linear covariate effects using coefficients learned on training rows only."""

    def fit(self, expression: np.ndarray, covariates: np.ndarray):
        expression = np.asarray(expression, dtype=float)
        covariates = np.asarray(covariates, dtype=float)
        self.covariate_mean_ = covariates.mean(axis=0)
        self.covariate_scale_ = covariates.std(axis=0, ddof=0)
        self.covariate_scale_[self.covariate_scale_ == 0] = 1.0
        standardized = (covariates - self.covariate_mean_) / self.covariate_scale_
        design = np.column_stack([np.ones(len(standardized)), standardized])
        self.coefficients_ = np.linalg.lstsq(design, expression, rcond=None)[0]
        self.n_features_in_ = expression.shape[1]
        return self

    def transform(self, expression: np.ndarray, covariates: np.ndarray) -> np.ndarray:
        expression = np.asarray(expression, dtype=float)
        covariates = np.asarray(covariates, dtype=float)
        standardized = (covariates - self.covariate_mean_) / self.covariate_scale_
        design = np.column_stack([np.ones(len(standardized)), standardized])
        covariate_effect = design[:, 1:] @ self.coefficients_[1:, :]
        return expression - covariate_effect

    def fit_transform(self, expression: np.ndarray, covariates: np.ndarray) -> np.ndarray:
        return self.fit(expression, covariates).transform(expression, covariates)
