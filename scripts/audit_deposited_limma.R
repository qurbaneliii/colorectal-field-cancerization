args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript scripts/audit_deposited_limma.R gene_matrix.csv output_directory")
}
stopifnot(as.character(packageVersion("limma")) == "3.64.3")

expression <- as.matrix(read.csv(args[[1]], row.names = 1L, check.names = FALSE))
metadata <- read.csv("data/metadata/gse44076_samples.csv", stringsAsFactors = FALSE)
metadata <- metadata[metadata$inclusion_status == "included", , drop = FALSE]
stopifnot(setequal(colnames(expression), metadata$geo_accession), all(is.finite(expression)))
reference <- read.csv("results/tables/supplementary_full_raw_cel_de_results.csv")
signature <- read.csv("results/tables/task_b_final_signature.csv")$gene_symbol
compare <- function(selected, formula, coefficient, model, comparison) {
  values <- expression[, selected$geo_accession, drop = FALSE]
  design <- model.matrix(formula, data = selected)
  stopifnot(nrow(design) == ncol(values), qr(design)$rank == ncol(design))
  fitted <- limma::eBayes(limma::lmFit(values, design))
  table <- limma::topTable(fitted, coef = coefficient, number = Inf, sort.by = "none")
  old <- reference[reference$model == model & reference$comparison == comparison, ]
  stopifnot(!anyDuplicated(old$gene_symbol), setequal(rownames(table), old$gene_symbol))
  old <- old[match(rownames(table), old$gene_symbol), ]
  original_hit <- old$adjusted_p_value < 0.05 & abs(old$log2_fold_change) >= 0.5
  deposited_hit <- table$adj.P.Val < 0.05 & abs(table$logFC) >= 0.5
  effect_threshold <- limma::treat(fitted, lfc = 0.5)
  treat_table <- limma::topTreat(
    effect_threshold, coef = coefficient, number = Inf, sort.by = "none"
  )
  stopifnot(identical(rownames(table), rownames(treat_table)))
  treat_hit <- treat_table$adj.P.Val < 0.05
  same_sign <- sign(old$log2_fold_change) == sign(table$logFC)
  matching <- match(intersect(signature, rownames(table)), rownames(table))
  list(
    summary = data.frame(
      model, comparison, samples = nrow(design), rank = qr(design)$rank, genes = nrow(table),
      pearson_log2fc = cor(old$log2_fold_change, table$logFC),
      spearman_log2fc = cor(old$log2_fold_change, table$logFC, method = "spearman"),
      sign_agreement = mean(same_sign), original_hits = sum(original_hit),
      deposited_hits = sum(deposited_hit), common_hits = sum(original_hit & deposited_hit),
      common_hits_same_sign = sum(original_hit & deposited_hit & same_sign),
      deposited_treat_hits = sum(treat_hit),
      deposited_provenance = "geo_deposited_matrix_raw_route_probe_annotation"
    ),
    panel = data.frame(
      model, comparison, gene_symbol = rownames(table)[matching],
      raw_archived_log2fc = old$log2_fold_change[matching],
      deposited_log2fc = table$logFC[matching], deposited_fdr = table$adj.P.Val[matching],
      deposited_treat_fdr = treat_table$adj.P.Val[matching]
    )
  )
}

field <- metadata[metadata$tissue_class %in% c("adjacent_normal", "healthy"), ]
field$tissue <- relevel(factor(field$tissue_class), ref = "healthy")
field$centered_age <- field$age - mean(field$age)
paired <- metadata[metadata$tissue_class %in% c("adjacent_normal", "tumor"), ]
paired$tissue <- relevel(factor(paired$tissue_class), ref = "adjacent_normal")
stopifnot(all(table(paired$patient_id, paired$tissue) == 1L))
results <- list(
  compare(field, ~tissue, "tissueadjacent_normal", "U0_unadjusted", "adjacent_vs_healthy"),
  compare(field, ~tissue + centered_age + sex, "tissueadjacent_normal",
          "U1_age_sex_adjusted", "adjacent_vs_healthy"),
  compare(field, ~tissue + centered_age + sex + location, "tissueadjacent_normal",
          "U2_age_sex_location_adjusted", "adjacent_vs_healthy"),
  compare(paired, ~factor(patient_id) + tissue, "tissuetumor",
          "P_patient_fixed_effect", "tumor_vs_adjacent")
)
dir.create(args[[2]], recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, lapply(results, `[[`, "summary")),
          file.path(args[[2]], "deposited_limma_concordance.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(results, `[[`, "panel")),
          file.path(args[[2]], "deposited_limma_panel_effects.csv"), row.names = FALSE)
