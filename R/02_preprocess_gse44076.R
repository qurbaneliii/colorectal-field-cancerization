suppressPackageStartupMessages({
  library(oligo)
  library(pd.hg.u219)
  library(hgu219.db)
  library(arrayQualityMetrics)
})
source("R/utils.R")

metadata_path <- "data/metadata/gse44076_samples.csv"
if (!file.exists(metadata_path)) stop("Run scripts/run_data_audit.py first")
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
cel_dir <- "data/interim/cel/GSE44076"
cel_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
assert_platform(cel_files, "GPL13667 / HG-U219")
if (length(cel_files) != nrow(metadata)) {
  stop("CEL/metadata mismatch: ", length(cel_files), " versus ", nrow(metadata))
}

raw <- oligo::read.celfiles(cel_files)
if (!inherits(annotation(raw), "character")) warning("Could not inspect array annotation")
eset <- oligo::rma(raw, target = "core")
probe_expr <- Biobase::exprs(eset)
if (median(probe_expr, na.rm = TRUE) > 20) stop("RMA output is not on expected log2 scale")

gsm <- sub("_.*$", "", basename(sampleNames(eset)))
colnames(probe_expr) <- gsm
metadata <- metadata[match(colnames(probe_expr), metadata$geo_accession), ]
stopifnot(!anyNA(metadata$geo_accession))

mapping <- unambiguous_mapping(rownames(probe_expr), hgu219.db)
write.csv(mapping, "data/metadata/GSE44076_probe_gene_mapping.csv", row.names = FALSE)
gene_expr <- aggregate_probes_by_gene(probe_expr, mapping)

write_expression(probe_expr, "data/processed/GSE44076_probe_expression.parquet", "probe_id")
write_expression(gene_expr, "data/processed/GSE44076_gene_expression.parquet")
write.csv(metadata, "data/processed/GSE44076_sample_metadata.csv", row.names = FALSE)
saveRDS(eset, "data/interim/GSE44076_rma_eset.rds", compress = "xz")
save_basic_qc(gene_expr, metadata, "GSE44076")

arrayQualityMetrics::arrayQualityMetrics(
  expressionset = eset, outdir = "results/figures/GSE44076_array_quality_metrics",
  force = TRUE, do.logtransform = FALSE
)
