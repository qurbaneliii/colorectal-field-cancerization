"""Screen explicit metadata in approved compute; synthetic mode is mandatory by default."""
import argparse
from datetime import date
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd

from src.data.independent_eligibility import screen
from src.data.provenance import file_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=["synthetic", "authorized-research"], required=True)
    parser.add_argument("--approval-record", type=Path)
    args = parser.parse_args()
    approval_hash = None
    if args.mode == "authorized-research":
        if args.approval_record is None:
            parser.error("Real execution requires an approval record in the approved environment")
        approval = json.loads(args.approval_record.read_text())
        fields = ["dataset", "approval_id", "authorized_investigator", "institution", "expiry",
                  "approved_compute", "ethics_record", "permitted_purpose"]
        if any(not str(approval.get(field, "")).strip() for field in fields):
            parser.error("Approval record lacks required fields; schema checking does not grant authorization")
        if date.fromisoformat(approval["expiry"]) < date.today():
            parser.error("Approval record has expired")
        approval_hash = file_sha256(args.approval_record)
    data = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    if args.mode == "synthetic" and not data.sample_id.str.startswith("SYNTHETIC_").all():
        parser.error("Synthetic mode requires explicitly prefixed sample identifiers")
    result = screen(data)
    args.output.mkdir(parents=True, exist_ok=False)
    result.to_csv(args.output / "eligibility.csv", index=False)
    result.loc[result.eligible].to_csv(args.output / "eligible_metadata.csv", index=False)
    report = {"mode": args.mode, "metadata_sha256": file_sha256(args.metadata),
              "approval_record_sha256": approval_hash, "primary_site": "sigmoid",
              "minimum_distance_cm": 10, "eligible_samples": int(result.eligible.sum()),
              "eligible_donors": int(result.loc[result.eligible].donor_id.nunique()),
              "exclusions": result.loc[~result.eligible].exclusion_reason.value_counts().to_dict(),
              "empirical_validation": "NOT RUN"}
    (args.output / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
