from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.provenance import artifact_sha256
from src.modeling.external_result import external_report, external_summary, write_external_report


ARTIFACTS = {
    "metrics": "results/metrics/external_validation_metrics.csv",
    "predictions": "results/tables/supplementary_external_predictions.csv",
    "threshold_grid": "results/tables/task_c_primary_threshold_selection.csv",
    "full_card": "models/task_c_primary_full_signature_model_card.json",
    "transport_card": "models/task_c_cross_platform_transport_model_card.json",
    "signature_audit": "results/tables/external_common_signature_audit.csv",
}
INTERVALS = "results/metrics/external_validation_bootstrap_ci.csv"
REPORT = "reports/external_validation_report.json"


def read_components(root: Path) -> dict:
    return {
        name: (
            json.loads((root / path).read_text(encoding="utf-8"))
            if path.endswith(".json")
            else pd.read_csv(root / path)
        )
        for name, path in ARTIFACTS.items()
    }


def verify(root: Path) -> dict:
    components = read_components(root)
    expected = external_summary(**components)
    metadata = pd.read_csv(
        root / "data/metadata/gse41258_samples.csv", dtype={"donor_or_patient_group": str}
    ).set_index("geo_accession")
    predictions = components["predictions"]
    matched = metadata.loc[predictions["sample_id"]]
    if (
        not matched["inclusion_status"].eq("included").all()
        or matched["tissue_class"].tolist() != predictions["y_true"].tolist()
        or matched["donor_or_patient_group"].tolist() != predictions["patient_id"].astype(str).tolist()
    ):
        raise ValueError("External predictions disagree with eligible GEO tissue and patient metadata")
    active = json.loads((root / REPORT).read_text(encoding="utf-8"))
    evidence = active.pop("evidence")
    if active != expected:
        raise ValueError("Active external JSON differs from the canonical result artifacts")
    if (root / REPORT).with_suffix(".md").read_text(encoding="utf-8") != external_report(expected):
        raise ValueError("Active external Markdown differs from the canonical result artifacts")
    for table in ("table_6_external_validation.csv", "table_9_external_validation.csv"):
        pd.testing.assert_frame_equal(
            pd.read_csv(root / "results/tables" / table),
            components["metrics"],
            check_exact=False,
            rtol=1e-12,
            atol=1e-12,
        )
    for name in (*ARTIFACTS.values(), INTERVALS):
        if artifact_sha256(root / name) != evidence["sha256"][name]:
            raise ValueError(f"External evidence hash mismatch for {name}")
    # Generated evidence also binds the report to the actual model and matrix bytes.
    # These files are deliberately absent from a checkout of a historical snapshot.
    if evidence["mode"] == "model_execution":
        for name, digest in evidence["sha256"].items():
            if artifact_sha256(root / name) != digest:
                raise ValueError(f"External execution input or output hash mismatch for {name}")
    elif evidence["mode"] != "historical_prediction_replay":
        raise ValueError("Unknown external result evidence mode")
    return expected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reconcile-historical-snapshot", metavar="SOURCE_COMMIT")
    args = parser.parse_args()
    if args.reconcile_historical_snapshot:
        if (ROOT / REPORT).exists():
            raise FileExistsError("Archive the old report before creating a historical replay")
        digests = {name: artifact_sha256(ROOT / name) for name in (*ARTIFACTS.values(), INTERVALS)}
        for name, digest in digests.items():
            recorded = subprocess.check_output(
                ["git", "show", f"{args.reconcile_historical_snapshot}:{name}"], cwd=ROOT
            )
            if hashlib.sha256(recorded.replace(b"\r\n", b"\n")).hexdigest() != digest:
                raise ValueError(f"Historical evidence differs from the declared snapshot: {name}")
        summary = external_summary(**read_components(ROOT))
        summary["evidence"] = {
            "mode": "historical_prediction_replay",
            "snapshot_commit": args.reconcile_historical_snapshot,
            "sha256": digests,
        }
        write_external_report(ROOT / "reports", summary)
    result = verify(ROOT)
    print(f"External result verified at locked threshold {result['locked_threshold']:.2f}")


if __name__ == "__main__":
    main()
