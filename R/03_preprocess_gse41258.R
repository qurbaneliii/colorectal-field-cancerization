suppressPackageStartupMessages({
  library(affy)
  library(hgu133acdf)
  library(hgu133a.db)
  library(affyPLM)
  library(arrayQualityMetrics)
})
source("R/utils.R")

metadata <- read.csv("data/metadata/gse41258_samples.csv", stringsAsFactors = FALSE)
cel_files <- list.files("data/interim/cel/GSE41258", pattern = "\\.CEL\\.gz$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
platform <- validate_platform(cel_files, "GSE41258", "GPL96", "HG-U133A",
                              "hgu133acdf / hgu133a.db", "affy")
if (length(cel_files) != 390L || nrow(metadata) != 390L) stop("GSE41258 expected 390 CEL/metadata rows")
metadata <- validate_cel_metadata_alignment(cel_files, metadata)

# ReadAffy detects the chip and CDF from the CEL header. No unverified cdfname is forced.
raw <- affy::ReadAffy(filenames = cel_files)
detected_annotation <- Biobase::annotation(raw)
detected_cdf <- affy::cdfName(raw)
if (!grepl("hgu133a", normalized_token(paste(detected_annotation, detected_cdf)), fixed = TRUE)) {
  stop("GSE41258 full CEL load did not resolve HG-U133A: ", detected_annotation, " / ", detected_cdf)
}
eset <- affy::rma(raw)
probe_expr <- Biobase::exprs(eset)
if (nrow(probe_expr) != 22283L) stop("Unexpected HG-U133A probe-set count: ", nrow(probe_expr))
if (!all(is.finite(probe_expr)) || median(probe_expr) > 20) stop("Invalid GSE41258 RMA scale/values")
colnames(probe_expr) <- extract_gsm(sampleNames(eset))
Biobase::sampleNames(eset) <- colnames(probe_expr)
metadata <- metadata[match(colnames(probe_expr), toupper(metadata$geo_accession)), , drop = FALSE]
if (anyNA(metadata$geo_accession)) stop("GSE41258 post-RMA metadata alignment failed")

mapping <- complete_probe_mapping(rownames(probe_expr), hgu133a.db)
write.csv(mapping, "data/metadata/GSE41258_probe_gene_mapping_raw_cel_rma.csv", row.names = FALSE)
gene_expr <- aggregate_probes_by_gene(probe_expr, mapping)
write_expression(probe_expr, "data/processed/GSE41258_probe_expression_raw_cel_rma.parquet",
                 "probe_id", "GSE41258")
write_expression(gene_expr, "data/processed/GSE41258_gene_expression_raw_cel_rma.parquet",
                 "gene_symbol", "GSE41258")
metadata$expression_provenance <- "raw_cel_rma"
write.csv(metadata, "data/processed/GSE41258_sample_metadata_raw_cel_rma.csv", row.names = FALSE)
saveRDS(eset, "data/interim/GSE41258_raw_cel_rma_eset.rds", compress = "xz")
qc <- save_qc_bundle(raw, gene_expr, metadata, "GSE41258")
update_sample_exclusion_log(metadata, qc, "GSE41258")
run_affyplm_nuse(raw, metadata, "GSE41258")
platform$detected_annotation <- detected_annotation
platform$detected_cdf <- detected_cdf
platform$normalized_probe_sets <- nrow(probe_expr)
platform$normalized_genes <- nrow(gene_expr)
write.csv(platform, "data/metadata/GSE41258_platform_validation.csv", row.names = FALSE)
if (!identical(Sys.getenv("SKIP_ARRAY_QUALITY_METRICS"), "1")) {
  run_array_quality_metrics(eset, metadata, "GSE41258")
}
