suppressPackageStartupMessages({
  library(oligo)
  library(pd.hg.u219)
  library(hgu219.db)
  library(arrayQualityMetrics)
})
source("R/utils.R")

metadata <- read.csv("data/metadata/gse44076_samples.csv", stringsAsFactors = FALSE)
cel_files <- list.files("data/interim/cel/GSE44076", pattern = "\\.CEL\\.gz$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
platform <- validate_platform(cel_files, "GSE44076", "GPL13667", "HG-U219",
                              "pd.hg.u219 / hgu219.db", "oligo")
if (length(cel_files) != 246L || nrow(metadata) != 246L) stop("GSE44076 expected 246 CEL/metadata rows")
metadata <- validate_cel_metadata_alignment(cel_files, metadata)

raw <- oligo::read.celfiles(cel_files, verbose = FALSE)
if (ncol(raw) != 246L) stop("GSE44076 full CEL load count mismatch")
eset <- oligo::rma(raw)
probe_expr <- Biobase::exprs(eset)
if (nrow(probe_expr) != 49386L) stop("Unexpected HG-U219 core probe-set count: ", nrow(probe_expr))
if (!all(is.finite(probe_expr)) || median(probe_expr) > 20) stop("Invalid GSE44076 RMA scale/values")
colnames(probe_expr) <- extract_gsm(sampleNames(eset))
metadata <- metadata[match(colnames(probe_expr), toupper(metadata$geo_accession)), , drop = FALSE]
if (anyNA(metadata$geo_accession)) stop("GSE44076 post-RMA metadata alignment failed")

mapping <- complete_probe_mapping(rownames(probe_expr), hgu219.db)
write.csv(mapping, "data/metadata/GSE44076_probe_gene_mapping_raw_cel_rma.csv", row.names = FALSE)
gene_expr <- aggregate_probes_by_gene(probe_expr, mapping)
write_expression(probe_expr, "data/processed/GSE44076_probe_expression_raw_cel_rma.parquet",
                 "probe_id", "GSE44076")
write_expression(gene_expr, "data/processed/GSE44076_gene_expression_raw_cel_rma.parquet",
                 "gene_symbol", "GSE44076")
metadata$expression_provenance <- "raw_cel_rma"
write.csv(metadata, "data/processed/GSE44076_sample_metadata_raw_cel_rma.csv", row.names = FALSE)
saveRDS(eset, "data/interim/GSE44076_raw_cel_rma_eset.rds", compress = "xz")
qc <- save_qc_bundle(raw, gene_expr, metadata, "GSE44076")
update_sample_exclusion_log(metadata, qc, "GSE44076")
platform$normalized_probe_sets <- nrow(probe_expr)
platform$normalized_genes <- nrow(gene_expr)
write.csv(platform, "data/metadata/GSE44076_platform_validation.csv", row.names = FALSE)
if (identical(Sys.getenv("RUN_ARRAY_QUALITY_METRICS"), "1")) {
  arrayQualityMetrics::arrayQualityMetrics(
    eset, outdir = "results/figures/GSE44076_array_quality_metrics",
    force = TRUE, do.logtransform = FALSE
  )
}
