from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    converted = []
    for accession in ("GSE44076", "GSE41258"):
        for level in ("probe", "gene"):
            parquet = (
                ROOT
                / "data/processed"
                / f"{accession}_{level}_expression_raw_cel_rma.parquet"
            )
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
    primary = pd.read_parquet(
        ROOT / "data/processed/GSE44076_gene_expression_raw_cel_rma.parquet",
        columns=["gene_symbol"],
    )
    external = pd.read_parquet(
        ROOT / "data/processed/GSE41258_gene_expression_raw_cel_rma.parquet",
        columns=["gene_symbol"],
    )
    common = sorted(set(primary["gene_symbol"].astype(str)).intersection(
        external["gene_symbol"].astype(str)
    ))
    common_path = (
        ROOT / "data/processed/common_genes_GSE44076_GSE41258_raw_cel_rma.txt"
    )
    common_path.write_text("\n".join(common) + "\n", encoding="utf-8")
    print(
        f"Converted {len(converted)} R CSV fallback files to parquet; "
        f"saved {len(common)} common raw-CEL genes"
    )


if __name__ == "__main__":
    main()
