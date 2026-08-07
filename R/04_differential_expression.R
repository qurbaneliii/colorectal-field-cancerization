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

format_table <- function(fit, coefficient, model, contrast, gene_to_entrez, provenance,
                         mapping_status = "mapped_unique") {
  if (!coefficient %in% colnames(fit$coefficients)) {
    stop("Requested coefficient '", coefficient, "' absent; available: ",
         paste(colnames(fit$coefficients), collapse = ", "))
  }
  tab <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "P",
                         adjust.method = "BH", confint = 0.95)
  standard_error <- fit$stdev.unscaled[, coefficient] * fit$sigma
  data.frame(
    gene_symbol = rownames(tab),
    entrez_id = unname(gene_to_entrez[rownames(tab)]),
    model = model,
    contrast = contrast,
    comparison = contrast,
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
}

clean_covariates <- function(frame) {
  frame$age <- suppressWarnings(as.numeric(frame$age))
  frame$sex <- trimws(as.character(frame$sex))
  frame$location <- trimws(as.character(frame$location))
  frame$sex[frame$sex == ""] <- NA_character_
  frame$location[frame$location == ""] <- NA_character_
  frame
}

expr <- read_expression("data/processed/GSE44076_gene_expression_raw_cel_rma.parquet")
meta <- read.csv("data/processed/GSE44076_sample_metadata_raw_cel_rma.csv",
                 stringsAsFactors = FALSE)
meta <- clean_covariates(meta[match(colnames(expr), meta$geo_accession), , drop = FALSE])
if (anyNA(meta$geo_accession)) stop("Raw-CEL expression/metadata alignment failed")
mapping <- read.csv("data/metadata/GSE44076_probe_gene_mapping_raw_cel_rma.csv",
                    stringsAsFactors = FALSE)
valid_mapping <- mapping[mapping$mapping_status == "mapped_unique", ]
gene_to_entrez <- tapply(valid_mapping$entrez_id, valid_mapping$gene_symbol,
                         function(x) sort(unique(x[x != ""]))[[1]])
dir.create("results/differential_expression", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

design_audits <- list()

run_unpaired <- function(level_a, level_b, contrast, model, covariates = character(0),
                         excluded_samples = character(0)) {
  keep <- meta$tissue_class %in% c(level_a, level_b) &
    !meta$geo_accession %in% excluded_samples
  required <- c("tissue_class", covariates)
  complete <- stats::complete.cases(meta[keep, required, drop = FALSE])
  selected_indices <- which(keep)[complete]
  selected <- droplevels(meta[selected_indices, , drop = FALSE])
  selected$tissue <- stats::relevel(factor(selected$tissue_class), ref = level_b)
  if ("sex" %in% covariates) selected$sex <- factor(selected$sex)
  if ("location" %in% covariates) selected$location <- factor(selected$location)
  formula <- stats::reformulate(c("tissue", covariates))
  design <- stats::model.matrix(formula, data = selected)
  coefficient <- paste0("tissue", level_a)
  design_rank <- qr(design)$rank
  full_rank <- design_rank == ncol(design)
  minimum_category_n <- min(table(selected$tissue))
  if ("sex" %in% covariates) minimum_category_n <- min(minimum_category_n,
                                                        min(table(selected$sex)))
  if ("location" %in% covariates) minimum_category_n <- min(minimum_category_n,
                                                             min(table(selected$location)))
  design_audits[[paste(model, contrast, paste(excluded_samples, collapse = ";"), sep = "|")]] <<-
    data.frame(
      model = model,
      contrast = contrast,
      formula = paste(deparse(formula), collapse = ""),
      samples_available = sum(keep),
      samples_complete_case = nrow(selected),
      healthy_or_reference_n = sum(selected$tissue_class == level_b),
      comparison_n = sum(selected$tissue_class == level_a),
      age_complete = sum(!is.na(selected$age)),
      sex_complete = sum(!is.na(selected$sex)),
      location_complete = sum(!is.na(selected$location)),
      age_mean_reference = mean(selected$age[selected$tissue_class == level_b], na.rm = TRUE),
      age_mean_comparison = mean(selected$age[selected$tissue_class == level_a], na.rm = TRUE),
      sex_distribution_reference = paste(names(table(selected$sex[selected$tissue_class == level_b])),
                                         table(selected$sex[selected$tissue_class == level_b]),
                                         sep = ":", collapse = ";"),
      sex_distribution_comparison = paste(names(table(selected$sex[selected$tissue_class == level_a])),
                                          table(selected$sex[selected$tissue_class == level_a]),
                                          sep = ":", collapse = ";"),
      location_distribution_reference = paste(
        names(table(selected$location[selected$tissue_class == level_b])),
        table(selected$location[selected$tissue_class == level_b]), sep = ":", collapse = ";"),
      location_distribution_comparison = paste(
        names(table(selected$location[selected$tissue_class == level_a])),
        table(selected$location[selected$tissue_class == level_a]), sep = ":", collapse = ";"),
      design_rows = nrow(design),
      design_columns = ncol(design),
      design_rank = design_rank,
      full_rank = full_rank,
      minimum_factor_category_n = minimum_category_n,
      condition_number = kappa(design),
      coefficient = coefficient,
      excluded_samples = paste(excluded_samples, collapse = ";"),
      stringsAsFactors = FALSE
    )
  if (!full_rank) stop(model, "/", contrast, ": design matrix is rank deficient")
  if (minimum_category_n < 10L) {
    stop(model, "/", contrast, ": sparse factor category (n<10)")
  }
  if (!coefficient %in% colnames(design)) {
    stop(model, "/", contrast, ": expected coefficient absent: ", coefficient)
  }
  fit <- limma::eBayes(limma::lmFit(expr[, selected_indices, drop = FALSE], design))
  format_table(fit, coefficient, model, contrast, gene_to_entrez,
               paste0("raw_cel_rma_limma_", tolower(model)))
}

results <- list(
  adjacent_u0 = run_unpaired("adjacent_normal", "healthy", "adjacent_vs_healthy",
                             "U0_unadjusted"),
  adjacent_u1 = run_unpaired("adjacent_normal", "healthy", "adjacent_vs_healthy",
                             "U1_age_sex_adjusted", c("age", "sex")),
  adjacent_u2 = run_unpaired("adjacent_normal", "healthy", "adjacent_vs_healthy",
                             "U2_age_sex_location_adjusted", c("age", "sex", "location")),
  tumor_u0 = run_unpaired("tumor", "healthy", "tumor_vs_healthy", "U0_unadjusted"),
  tumor_u1 = run_unpaired("tumor", "healthy", "tumor_vs_healthy",
                          "U1_age_sex_adjusted", c("age", "sex")),
  tumor_u2 = run_unpaired("tumor", "healthy", "tumor_vs_healthy",
                          "U2_age_sex_location_adjusted", c("age", "sex", "location"))
)

paired <- meta$tissue_class %in% c("adjacent_normal", "tumor") &
  meta$pairing_status == "paired"
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
paired_coefficient <- grep("^tissue_class", colnames(paired_design), value = TRUE)
if (!identical(paired_coefficient, "tissue_classtumor")) {
  stop("Unexpected paired-design tissue coefficient: ", paste(paired_coefficient, collapse = ", "))
}
if (qr(paired_design)$rank != ncol(paired_design)) stop("Paired fixed-effect design is rank deficient")
paired_fit <- limma::eBayes(limma::lmFit(pair_expr, paired_design))
results$paired <- format_table(
  paired_fit, paired_coefficient, "P_patient_fixed_effect", "tumor_vs_adjacent",
  gene_to_entrez, "raw_cel_rma_limma_patient_fixed_effect"
)
design_audits$paired <- data.frame(
  model = "P_patient_fixed_effect", contrast = "tumor_vs_adjacent",
  formula = "~ patient_id + tissue_class", samples_available = nrow(pair_meta),
  samples_complete_case = nrow(pair_meta), healthy_or_reference_n = 98L,
  comparison_n = 98L, age_complete = sum(!is.na(pair_meta$age)),
  sex_complete = sum(!is.na(pair_meta$sex)),
  location_complete = sum(!is.na(pair_meta$location)),
  age_mean_reference = mean(pair_meta$age[pair_meta$tissue_class == "adjacent_normal"]),
  age_mean_comparison = mean(pair_meta$age[pair_meta$tissue_class == "tumor"]),
  sex_distribution_reference = "absorbed_by_patient_fixed_effects",
  sex_distribution_comparison = "absorbed_by_patient_fixed_effects",
  location_distribution_reference = "absorbed_by_patient_fixed_effects",
  location_distribution_comparison = "absorbed_by_patient_fixed_effects",
  design_rows = nrow(paired_design), design_columns = ncol(paired_design),
  design_rank = qr(paired_design)$rank, full_rank = TRUE,
  minimum_factor_category_n = 98L, condition_number = kappa(paired_design),
  coefficient = paired_coefficient, excluded_samples = "", stringsAsFactors = FALSE
)

output_paths <- c(
  adjacent_u0 = "results/differential_expression/raw_cel_adjacent_vs_healthy_unadjusted.csv",
  adjacent_u1 = "results/differential_expression/raw_cel_adjacent_vs_healthy_age_sex_adjusted.csv",
  adjacent_u2 = "results/differential_expression/raw_cel_adjacent_vs_healthy_age_sex_location_adjusted.csv",
  tumor_u0 = "results/differential_expression/raw_cel_tumor_vs_healthy_unadjusted.csv",
  tumor_u1 = "results/differential_expression/raw_cel_tumor_vs_healthy_age_sex_adjusted.csv",
  tumor_u2 = "results/differential_expression/raw_cel_tumor_vs_healthy_age_sex_location_adjusted.csv",
  paired = "results/differential_expression/raw_cel_tumor_vs_adjacent_patient_fixed_effect.csv"
)
for (name in names(output_paths)) write.csv(results[[name]], output_paths[[name]], row.names = FALSE)

# Backward-compatible aliases remain explicit about which model they contain.
write.csv(results$adjacent_u0,
          "results/differential_expression/raw_cel_adjacent_vs_healthy.csv", row.names = FALSE)
write.csv(results$tumor_u0,
          "results/differential_expression/raw_cel_tumor_vs_healthy.csv", row.names = FALSE)
write.csv(results$paired,
          "results/differential_expression/raw_cel_tumor_vs_adjacent_paired.csv", row.names = FALSE)

qc_path <- "results/tables/GSE44076_raw_cel_qc_diagnostics.csv"
if (file.exists(qc_path)) {
  qc <- read.csv(qc_path, stringsAsFactors = FALSE)
  borderline <- qc$geo_accession[qc$independent_failure_metrics == 1L]
  results$adjacent_u1_qc_excluded <- run_unpaired(
    "adjacent_normal", "healthy", "adjacent_vs_healthy",
    "U1_age_sex_adjusted_qc_sensitivity", c("age", "sex"), borderline
  )
  write.csv(results$adjacent_u1_qc_excluded,
            "results/differential_expression/raw_cel_adjacent_vs_healthy_age_sex_adjusted_qc_excluded.csv",
            row.names = FALSE)
}

all_de <- do.call(rbind, results)
write.csv(all_de, "results/tables/supplementary_full_raw_cel_de_results.csv", row.names = FALSE)

thresholds <- expand.grid(fdr = c(0.01, 0.05, 0.10), absolute_log2fc = c(0.25, 0.5, 1.0))
sensitivity <- do.call(rbind, lapply(split(all_de, list(all_de$contrast, all_de$model),
                                          drop = TRUE), function(frame) {
  do.call(rbind, lapply(seq_len(nrow(thresholds)), function(i) {
    data.frame(contrast = frame$contrast[[1]], model = frame$model[[1]], thresholds[i, ],
               significant_genes = sum(frame$adjusted_p_value < thresholds$fdr[[i]] &
                                         abs(frame$log2_fold_change) >=
                                           thresholds$absolute_log2fc[[i]]))
  }))
}))
write.csv(sensitivity, "results/tables/raw_cel_de_threshold_sensitivity.csv", row.names = FALSE)

design_audit <- do.call(rbind, design_audits)
rownames(design_audit) <- NULL
write.csv(design_audit, "results/tables/covariate_design_audit.csv", row.names = FALSE)
write.csv(design_audit, "results/tables/raw_cel_limma_design_audit.csv", row.names = FALSE)

report <- c(
  "# Covariate design audit", "",
  "Age, sex, and tumor location were taken only from deposited GSE44076 metadata.",
  "Stage was not used in healthy-versus-adjacent or tumor-versus-healthy models because",
  "it is structurally undefined for cancer-free donors. Age was modeled continuously;",
  "sex and location were modeled as factors. U2 was executed only after confirming",
  "a full-rank design and a minimum of ten observations in every factor category.", "",
  paste(capture.output(print(design_audit, row.names = FALSE)), collapse = "\n"), "",
  "The paired tumor-versus-adjacent model retains patient fixed effects, which absorb",
  "patient-level age, sex, and location."
)
writeLines(report, "reports/covariate_design_audit.md")
