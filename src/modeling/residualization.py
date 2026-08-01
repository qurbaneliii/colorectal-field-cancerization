from __future__ import annotations

import numpy as np


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
