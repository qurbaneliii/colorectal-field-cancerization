suppressPackageStartupMessages(library(limma))
source("R/utils.R")

read_expression <- function(path) {
  if (file.exists(path) && requireNamespace("arrow", quietly = TRUE)) {
    x <- as.data.frame(arrow::read_parquet(path))
  } else {
    fallback <- sub("\\.parquet$", ".csv.gz", path)
    if (!file.exists(fallback)) stop("Missing raw-CEL expression file: ", path)
    x <- data.table::fread(fallback, data.table = FALSE)
  }
  ids <- x[[1]]
  mat <- as.matrix(x[-1])
  rownames(mat) <- ids
  if (!all(is.finite(mat))) stop("Expression contains non-finite values")
  mat
}

format_table <- function(fit, coefficient, comparison, gene_to_entrez, provenance,
                         mapping_status = "mapped_unique") {
  if (!coefficient %in% colnames(fit$coefficients)) {
    stop("Requested coefficient '", coefficient, "' absent; available: ",
         paste(colnames(fit$coefficients), collapse = ", "))
  }
  tab <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "P",
                         adjust.method = "BH", confint = 0.95)
  standard_error <- fit$stdev.unscaled[, coefficient] * fit$sigma
  result <- data.frame(
    gene_symbol = rownames(tab),
    entrez_id = unname(gene_to_entrez[rownames(tab)]),
    comparison = comparison,
    log2_fold_change = tab$logFC,
    standard_error = standard_error[rownames(tab)],
    moderated_statistic = tab$t,
    p_value = tab$P.Value,
    adjusted_p_value = tab$adj.P.Val,
    average_expression = tab$AveExpr,
    confidence_interval_lower = tab$CI.L,
    confidence_interval_upper = tab$CI.R,
    direction = ifelse(tab$logFC >= 0, "up", "down"),
    mapping_status = mapping_status,
    analysis_provenance = provenance,
    stringsAsFactors = FALSE
  )
  result
}

expr <- read_expression("data/processed/GSE44076_gene_expression_raw_cel_rma.parquet")
meta <- read.csv("data/processed/GSE44076_sample_metadata_raw_cel_rma.csv",
                 stringsAsFactors = FALSE)
meta <- meta[match(colnames(expr), meta$geo_accession), , drop = FALSE]
if (anyNA(meta$geo_accession)) stop("Raw-CEL expression/metadata alignment failed")
mapping <- read.csv("data/metadata/GSE44076_probe_gene_mapping_raw_cel_rma.csv",
                    stringsAsFactors = FALSE)
valid_mapping <- mapping[mapping$mapping_status == "mapped_unique", ]
gene_to_entrez <- tapply(valid_mapping$entrez_id, valid_mapping$gene_symbol,
                         function(x) sort(unique(x[x != ""]))[[1]])
dir.create("results/differential_expression", recursive = TRUE, showWarnings = FALSE)

run_unpaired <- function(level_a, level_b, label) {
  keep <- meta$tissue_class %in% c(level_a, level_b)
  tissue <- factor(meta$tissue_class[keep], levels = c(level_b, level_a))
  design <- model.matrix(~ 0 + tissue)
  expected_columns <- paste0("tissue", c(level_b, level_a))
  if (!identical(colnames(design), expected_columns)) {
    stop(label, ": unexpected design columns: ", paste(colnames(design), collapse = ", "))
  }
  contrast_text <- paste0("`", expected_columns[[2]], "`-`", expected_columns[[1]], "`")
  contrast <- limma::makeContrasts(contrasts = contrast_text, levels = design)
  colnames(contrast) <- label
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr[, keep, drop = FALSE], design),
                                            contrast))
  format_table(fit, label, label, gene_to_entrez, "raw_cel_rma_limma")
}

adjacent_healthy <- run_unpaired("adjacent_normal", "healthy", "adjacent_vs_healthy")
tumor_healthy <- run_unpaired("tumor", "healthy", "tumor_vs_healthy")

paired <- meta$tissue_class %in% c("adjacent_normal", "tumor") & meta$pairing_status == "paired"
pair_meta <- droplevels(meta[paired, , drop = FALSE])
pair_expr <- expr[, paired, drop = FALSE]
pair_counts <- table(pair_meta$patient_id, pair_meta$tissue_class)
if (nrow(pair_counts) != 98L || any(pair_counts[, c("adjacent_normal", "tumor")] != 1L)) {
  stop("Paired limma cohort is not 98 complete one-to-one patient pairs")
}
pair_meta$patient_id <- factor(pair_meta$patient_id)
pair_meta$tissue_class <- factor(pair_meta$tissue_class,
                                 levels = c("adjacent_normal", "tumor"))
paired_design <- model.matrix(~ patient_id + tissue_class, data = pair_meta)
coefficient <- grep("^tissue_class", colnames(paired_design), value = TRUE)
if (!identical(coefficient, "tissue_classtumor")) {
  stop("Unexpected paired-design tissue coefficient: ", paste(coefficient, collapse = ", "))
}
if (qr(paired_design)$rank != ncol(paired_design)) stop("Paired fixed-effect design is rank deficient")
paired_fit <- limma::eBayes(limma::lmFit(pair_expr, paired_design))
tumor_adjacent <- format_table(paired_fit, coefficient, "tumor_vs_adjacent_paired",
                               gene_to_entrez, "raw_cel_rma_limma_patient_fixed_effect")

write.csv(adjacent_healthy,
          "results/differential_expression/raw_cel_adjacent_vs_healthy.csv", row.names = FALSE)
write.csv(tumor_healthy,
          "results/differential_expression/raw_cel_tumor_vs_healthy.csv", row.names = FALSE)
write.csv(tumor_adjacent,
          "results/differential_expression/raw_cel_tumor_vs_adjacent_paired.csv", row.names = FALSE)
all_de <- rbind(adjacent_healthy, tumor_healthy, tumor_adjacent)
write.csv(all_de, "results/tables/supplementary_full_raw_cel_de_results.csv", row.names = FALSE)

thresholds <- expand.grid(fdr = c(0.01, 0.05, 0.10), absolute_log2fc = c(0.25, 0.5, 1.0))
sensitivity <- do.call(rbind, lapply(split(all_de, all_de$comparison), function(frame) {
  do.call(rbind, lapply(seq_len(nrow(thresholds)), function(i) {
    data.frame(comparison = frame$comparison[[1]], thresholds[i, ],
               significant_genes = sum(frame$adjusted_p_value < thresholds$fdr[[i]] &
                                         abs(frame$log2_fold_change) >= thresholds$absolute_log2fc[[i]]))
  }))
}))
write.csv(sensitivity, "results/tables/raw_cel_de_threshold_sensitivity.csv", row.names = FALSE)

design_audit <- data.frame(
  analysis = c("adjacent_vs_healthy", "tumor_vs_healthy", "tumor_vs_adjacent_paired"),
  design = c("unpaired tissue indicator", "unpaired tissue indicator",
             "patient fixed effects plus tissue indicator"),
  samples = c(sum(meta$tissue_class %in% c("adjacent_normal", "healthy")),
              sum(meta$tissue_class %in% c("tumor", "healthy")), nrow(pair_meta)),
  patients_or_donors = c(length(unique(meta$donor_or_patient_group[
    meta$tissue_class %in% c("adjacent_normal", "healthy")])),
    length(unique(meta$donor_or_patient_group[meta$tissue_class %in% c("tumor", "healthy")])),
    length(unique(pair_meta$patient_id))),
  coefficient = c("tissueadjacent_normal-tissuehealthy", "tissuetumor-tissuehealthy",
                  coefficient), stringsAsFactors = FALSE
)
write.csv(design_audit, "results/tables/raw_cel_limma_design_audit.csv", row.names = FALSE)
