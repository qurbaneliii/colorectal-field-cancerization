from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.provenance import expression_path, result_root
from src.modeling.evaluation import (
    classification_metrics,
    grouped_bootstrap_metrics,
    prediction_probabilities,
)
from src.modeling.pipelines import build_pipeline
from src.modeling.residualization import TrainingCovariateResidualizer
from src.modeling.splitters import stratified_group_splits
from scripts.run_compact_panel_analysis import (
    choose_panel,
    inner_panel_scores,
    panel_specs,
    rank_training_features,
)


def demographic_matrix(metadata: pd.DataFrame) -> np.ndarray:
    return np.column_stack(
        [
            metadata["age"].astype(float),
            metadata["sex"].eq("Male").astype(float),
            metadata["location"].eq("Right").astype(float),
        ]
    )


def match_age_sex(metadata: pd.DataFrame, caliper: float) -> pd.Index:
    healthy = metadata[metadata["tissue_class"].eq("healthy")].sort_values("geo_accession")
    adjacent = metadata[metadata["tissue_class"].eq("adjacent_normal")].copy()
    available = set(adjacent.index)
    retained: list[int] = []
    for index, row in healthy.iterrows():
        candidates = adjacent.loc[list(available)]
        candidates = candidates[
            candidates["sex"].eq(row["sex"])
            & candidates["age"].sub(float(row["age"])).abs().le(caliper)
        ].copy()
        if candidates.empty:
            continue
        candidates["distance"] = candidates["age"].sub(float(row["age"])).abs()
        chosen = int(candidates.sort_values(["distance", "geo_accession"]).index[0])
        retained.extend([int(index), chosen])
        available.remove(chosen)
    return pd.Index(retained)


def elastic_panel(parameters: dict[str, object], seed: int, modeling: dict) -> Pipeline:
    model = build_pipeline("elastic_net", seed, modeling)
    model.set_params(
        **{
            **parameters,
            "variance_quantile__quantile": 0.0,
            "univariate__k": "all",
            "memory": None,
        }
    )
    return model


def evaluate(
    scenario: str,
    expression: np.ndarray | None,
    metadata: pd.DataFrame,
    parameters: dict[str, object],
    modeling: dict,
    covariates: np.ndarray | None = None,
    demographic_only: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    labels = metadata["tissue_class"].to_numpy()
    groups = metadata["donor_or_patient_group"].astype(str).to_numpy()
    samples = metadata["geo_accession"].to_numpy()
    predictions: list[dict[str, object]] = []
    metrics: list[dict[str, object]] = []
    selections: list[dict[str, object]] = []
    for repeat, seed in enumerate(modeling["seeds"]):
        for fold, (train, test) in enumerate(
            stratified_group_splits(labels, groups, int(modeling["outer_splits"]), int(seed))
        ):
            if demographic_only:
                model = Pipeline(
                    [
                        ("scale", StandardScaler()),
                        (
                            "model",
                            LogisticRegression(
                                solver="saga",
                                C=float(parameters["model__C"]),
                                l1_ratio=float(parameters["model__l1_ratio"]),
                                class_weight="balanced",
                                max_iter=5000,
                                random_state=int(seed),
                            ),
                        ),
                    ]
                )
                train_x = demographic_matrix(metadata.iloc[train])
                test_x = demographic_matrix(metadata.iloc[test])
            else:
                if expression is None:
                    raise AssertionError("Expression is required for molecular scenarios")
                train_x = expression[train]
                test_x = expression[test]
                if covariates is not None:
                    residualizer = TrainingCovariateResidualizer()
                    train_x = residualizer.fit_transform(train_x, covariates[train])
                    test_x = residualizer.transform(test_x, covariates[test])
                model = elastic_panel(parameters, int(seed), modeling)
            model.fit(train_x, labels[train])
            predicted = model.predict(test_x)
            classes, probability = prediction_probabilities(model, test_x)
            measured = classification_metrics(labels[test], predicted, probability, classes)
            metrics.append(
                {
                    "scenario": scenario,
                    "repeat": repeat,
                    "outer_fold": fold,
                    "n_train": len(train),
                    "n_test": len(test),
                    **measured,
                }
            )
            if not demographic_only:
                coefficients = np.asarray(model.named_steps["model"].coef_).reshape(-1)
                for gene_index, coefficient in enumerate(coefficients):
                    selections.append(
                        {
                            "scenario": scenario,
                            "repeat": repeat,
                            "outer_fold": fold,
                            "gene_index": gene_index,
                            "coefficient": float(coefficient),
                            "selected_nonzero": bool(abs(coefficient) > 1e-12),
                        }
                    )
            for position, sample_index in enumerate(test):
                row = {
                    "scenario": scenario,
                    "repeat": repeat,
                    "outer_fold": fold,
                    "sample_id": samples[sample_index],
                    "patient_id": groups[sample_index],
                    "y_true": labels[sample_index],
                    "y_pred": predicted[position],
                }
                row.update(
                    {
                        f"probability_{label}": float(probability[position, class_index])
                        for class_index, label in enumerate(classes)
                    }
                )
                predictions.append(row)
    return pd.DataFrame(metrics), pd.DataFrame(predictions), pd.DataFrame(selections)


def average_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    probability_columns = sorted(column for column in frame if column.startswith("probability_"))
    averaged = frame.groupby("sample_id", as_index=False).agg(
        {
            "patient_id": "first",
            "y_true": "first",
            **{column: "mean" for column in probability_columns},
        }
    )
    classes = np.asarray([column.removeprefix("probability_") for column in probability_columns])
    averaged["y_pred"] = classes[np.argmax(averaged[probability_columns].to_numpy(), axis=1)]
    return averaged


def main() -> None:
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    provenance = config["modeling"]["expression_provenance"]
    tables = result_root(ROOT / paths["tables"], provenance)
    metrics_root = result_root(ROOT / paths["metrics"], provenance)
    reports = result_root(ROOT / paths["reports"], provenance)
    signature = pd.read_csv(tables / "task_b_final_signature.csv")
    genes = signature["gene_symbol"].astype(str).tolist()
    expression = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE44076", "gene", provenance)
    ).set_index("gene_symbol")
    metadata = pd.read_csv(ROOT / paths["metadata"] / "gse44076_samples.csv")
    metadata = metadata[
        metadata["inclusion_status"].eq("included")
        & metadata["tissue_class"].isin(["healthy", "adjacent_normal"])
    ].reset_index(drop=True)
    panel = expression.loc[genes, metadata["geo_accession"]].T.to_numpy(dtype=float)
    card = json.loads((ROOT / paths["models"] / "task_b_final_model_card.json").read_text())
    parameters = card["hyperparameters"]
    demographic = demographic_matrix(metadata)
    composition_frame = pd.read_csv(tables / "tissue_composition_scores.csv").set_index(
        "geo_accession"
    )
    score_columns = [
        column
        for column in composition_frame
        if column
        not in {
            "patient_id",
            "donor_or_patient_group",
            "tissue_class",
            "age",
            "sex",
            "location",
        }
    ]
    composition = composition_frame.loc[metadata["geo_accession"], score_columns].to_numpy()
    combined_covariates = np.column_stack([demographic, composition])
    matched_indices = match_age_sex(
        metadata, float(config["modeling"]["demographic_matching_age_caliper"])
    )
    qc_diagnostics = pd.read_csv(tables / "GSE44076_raw_cel_qc_diagnostics.csv")
    borderline_samples = set(
        qc_diagnostics.loc[
            qc_diagnostics["independent_failure_metrics"].eq(1), "geo_accession"
        ].astype(str)
    )
    qc_keep = ~metadata["geo_accession"].astype(str).isin(borderline_samples)

    scenarios = [
        ("primary_fixed_panel", panel, metadata, None, False),
        ("demographic_only", None, metadata, None, True),
        ("age_sex_location_residualized", panel, metadata, demographic, False),
        ("composition_and_demographic_residualized", panel, metadata, combined_covariates, False),
        (
            "age_sex_matched_fixed_panel",
            panel[matched_indices],
            metadata.loc[matched_indices].reset_index(drop=True),
            None,
            False,
        ),
        (
            "borderline_sample_excluded_fixed_panel",
            panel[qc_keep],
            metadata.loc[qc_keep].reset_index(drop=True),
            None,
            False,
        ),
    ]
    all_metrics: list[pd.DataFrame] = []
    all_predictions: list[pd.DataFrame] = []
    all_selections: list[pd.DataFrame] = []
    for name, values, frame, covariates, demographic_only in scenarios:
        metric, prediction, selection = evaluate(
            name,
            values,
            frame,
            parameters,
            config["modeling"],
            covariates=covariates,
            demographic_only=demographic_only,
        )
        all_metrics.append(metric)
        all_predictions.append(prediction)
        if not selection.empty:
            selection["gene_symbol"] = selection["gene_index"].map(dict(enumerate(genes)))
            all_selections.append(selection)
    fold_metrics = pd.concat(all_metrics, ignore_index=True)
    predictions = pd.concat(all_predictions, ignore_index=True)
    fold_metrics.to_csv(metrics_root / "task_b_confounding_fold_metrics.csv", index=False)
    predictions.to_csv(metrics_root / "task_b_confounding_oof_predictions.csv", index=False)

    summary_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    for scenario, frame in predictions.groupby("scenario"):
        averaged = average_predictions(frame)
        classes = np.asarray(sorted(averaged["y_true"].unique()))
        probabilities = averaged[[f"probability_{label}" for label in classes]].to_numpy()
        measured = classification_metrics(
            averaged["y_true"].to_numpy(),
            averaged["y_pred"].to_numpy(),
            probabilities,
            classes,
            include_calibration=True,
        )
        summary_rows.append(
            {
                "scenario": scenario,
                "n_samples": len(averaged),
                "n_patients": averaged["patient_id"].nunique(),
                "panel_genes_locked": scenario != "demographic_only",
                **measured,
            }
        )
        for row in grouped_bootstrap_metrics(
            averaged,
            ["f1_macro", "balanced_accuracy", "roc_auc", "brier_score"],
            int(config["modeling"]["bootstrap_iterations"]),
            int(config["project"]["random_seed"]),
        ):
            bootstrap_rows.append({"scenario": scenario, **row})
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(tables / "task_b_confounding_sensitivity.csv", index=False)
    pd.DataFrame(bootstrap_rows).to_csv(
        tables / "task_b_confounding_sensitivity_bootstrap_ci.csv", index=False
    )
    selection_frame = pd.concat(all_selections, ignore_index=True)
    selection_summary = (
        selection_frame.groupby(["scenario", "gene_symbol"], as_index=False)
        .agg(
            selected_folds=("selected_nonzero", "sum"),
            selection_frequency=("selected_nonzero", "mean"),
            median_coefficient=("coefficient", "median"),
        )
    )
    selection_summary["panel_was_locked_before_sensitivity"] = True
    selection_summary.to_csv(tables / "task_b_confounding_gene_sensitivity.csv", index=False)

    qc_metadata = metadata.loc[qc_keep].reset_index(drop=True)
    qc_expression = expression.loc[:, qc_metadata["geo_accession"]].T.to_numpy(dtype=float)
    qc_labels = qc_metadata["tissue_class"].to_numpy()
    qc_groups = qc_metadata["donor_or_patient_group"].astype(str).to_numpy()
    feature_names = expression.index.astype(str).tolist()
    qc_inner = inner_panel_scores(
        qc_expression,
        qc_labels,
        qc_groups,
        feature_names,
        parameters,
        int(config["project"]["random_seed"]) + 9100,
        config["modeling"],
    )
    qc_size, _ = choose_panel(qc_inner, config["modeling"])
    qc_ranked, _ = rank_training_features(
        qc_expression,
        qc_labels,
        feature_names,
        parameters,
        int(config["project"]["random_seed"]) + 9101,
        config["modeling"],
    )
    qc_genes = panel_specs(qc_ranked, config["modeling"]["compact_panels"])[qc_size]
    qc_gene_comparison = pd.DataFrame(
        {
            "gene_symbol": sorted(set(genes) | set(qc_genes)),
        }
    )
    qc_gene_comparison["primary_panel"] = qc_gene_comparison["gene_symbol"].isin(genes)
    qc_gene_comparison["qc_excluded_reselected_panel"] = qc_gene_comparison[
        "gene_symbol"
    ].isin(qc_genes)
    qc_gene_comparison["primary_panel_size"] = len(genes)
    qc_gene_comparison["qc_excluded_panel_size"] = len(qc_genes)
    qc_gene_comparison["selection_origin"] = (
        "full-development grouped inner CV; sensitivity excludes custom-QC borderline sample"
    )
    qc_gene_comparison.to_csv(
        tables / "qc_exclusion_model_gene_comparison.csv", index=False
    )

    primary = summary[summary["scenario"].eq("primary_fixed_panel")].iloc[0]
    report = f"""# Task B confounding sensitivity

Task B is the primary healthy-versus-adjacent endpoint. All molecular
sensitivity analyses use the five-gene panel locked by grouped inner CV. The
covariate residualizers are fitted on each outer-training fold and applied to
the held-out fold, preventing outcome or test-row leakage. The demographic-only
model uses age, sex, and left/right location. The matching analysis uses
one-to-one, without-replacement same-sex matching within a
{config['modeling']['demographic_matching_age_caliper']}-year age caliper.

- Locked panel: {', '.join(genes)}
- Primary fixed-panel aggregated OOF macro-F1: {primary['f1_macro']:.4f}
- Matched subset: {len(matched_indices)} samples ({len(matched_indices) // 2} per class)
- Custom-QC exclusion sensitivity: {len(qc_metadata)} samples; reselected panel
  ({qc_size} genes): {', '.join(qc_genes)}
- Independent Task B validation: unavailable
- Interpretation: internally validated association model; not clinically ready

## Aggregated repeated OOF estimates

{summary.to_markdown(index=False)}

Gene coefficients and selection frequencies for the confounding scenarios are
reported separately. Those scenarios keep the primary panel locked, so
coefficient shrinkage to zero is their predefined gene-level sensitivity
measure. The QC-only analysis additionally repeats full-development grouped
inner-CV panel selection after omitting the borderline array; it is reported as
a sensitivity panel and never replaces the primary model.
"""
    (reports / "task_b_confounding_sensitivity.md").write_text(report, encoding="utf-8")
    print(summary.to_json(orient="records", indent=2))


if __name__ == "__main__":
    main()
