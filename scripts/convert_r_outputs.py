from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    converted = []
    for accession in ("GSE44076", "GSE41258"):
        for level in ("probe", "gene"):
            parquet = ROOT / "data/processed" / f"{accession}_{level}_expression.parquet"
            fallback = parquet.with_suffix(".csv.gz")
            if parquet.exists():
                continue
            if not fallback.exists():
                raise FileNotFoundError(
                    f"Missing both R parquet and CSV fallback for {accession} {level} expression"
                )
            frame = pd.read_csv(fallback)
            frame.to_parquet(parquet, index=False)
            converted.append(str(parquet.relative_to(ROOT)))
    print(f"Converted {len(converted)} R CSV fallback files to parquet")


if __name__ == "__main__":
    main()
