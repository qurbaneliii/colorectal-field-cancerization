suppressPackageStartupMessages({
  library(affy)
  library(hgu133acdf)
  library(hgu133a.db)
  library(arrayQualityMetrics)
})
source("R/utils.R")

metadata_path <- "data/metadata/gse41258_samples.csv"
if (!file.exists(metadata_path)) stop("Run scripts/run_data_audit.py first")
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
cel_dir <- "data/interim/cel/GSE41258"
cel_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
assert_platform(cel_files, "GPL96 / HG-U133A")
if (length(cel_files) != nrow(metadata)) {
  stop("CEL/metadata mismatch: ", length(cel_files), " versus ", nrow(metadata))
}

raw <- affy::ReadAffy(filenames = cel_files, cdfname = "hgu133acdf")
eset <- affy::rma(raw)
probe_expr <- Biobase::exprs(eset)
if (median(probe_expr, na.rm = TRUE) > 20) stop("RMA output is not on expected log2 scale")
gsm <- sub("_.*$", "", basename(sampleNames(eset)))
colnames(probe_expr) <- gsm
metadata <- metadata[match(colnames(probe_expr), metadata$geo_accession), ]
stopifnot(!anyNA(metadata$geo_accession))

mapping <- unambiguous_mapping(rownames(probe_expr), hgu133a.db)
write.csv(mapping, "data/metadata/GSE41258_probe_gene_mapping.csv", row.names = FALSE)
gene_expr <- aggregate_probes_by_gene(probe_expr, mapping)
write_expression(probe_expr, "data/processed/GSE41258_probe_expression.parquet", "probe_id")
write_expression(gene_expr, "data/processed/GSE41258_gene_expression.parquet")
write.csv(metadata, "data/processed/GSE41258_sample_metadata.csv", row.names = FALSE)
saveRDS(eset, "data/interim/GSE41258_rma_eset.rds", compress = "xz")
save_basic_qc(gene_expr, metadata, "GSE41258")

arrayQualityMetrics::arrayQualityMetrics(
  expressionset = eset, outdir = "results/figures/GSE41258_array_quality_metrics",
  force = TRUE, do.logtransform = FALSE
)
