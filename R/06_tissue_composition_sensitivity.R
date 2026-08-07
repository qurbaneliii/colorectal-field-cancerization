suppressPackageStartupMessages({
  library(limma)
  library(MCPcounter)
  library(ggplot2)
  library(digest)
})
source("R/utils.R")

read_expression <- function(path) {
  if (file.exists(path) && requireNamespace("arrow", quietly = TRUE)) {
    x <- as.data.frame(arrow::read_parquet(path))
  } else {
    fallback <- sub("\\.parquet$", ".csv.gz", path)
    if (!file.exists(fallback)) stop("Missing expression file: ", path)
    x <- data.table::fread(fallback, data.table = FALSE)
  }
  ids <- x[[1]]
  mat <- as.matrix(x[-1])
  rownames(mat) <- ids
  mat
}

save_plot <- function(plot, stem, width = 12, height = 8) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, bg = "white")
}

signature_url <- paste0(
  "https://raw.githubusercontent.com/ebecht/MCPcounter/",
  "b6eac73e91c246fcff0bb1a5c68a816cd588fc48/Signatures/genes.txt"
)
signature_sha256 <- "408f6c5d02c8f9bd2f1599c598367881853520934f8ce398a4886a8b296922bb"
signature_path <- "data/reference/mcpcounter_genes_b6eac73e.tsv"
dir.create(dirname(signature_path), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(signature_path)) {
  download.file(signature_url, signature_path, mode = "wb", quiet = FALSE)
}
observed_sha256 <- digest::digest(file = signature_path, algo = "sha256")
if (!identical(observed_sha256, signature_sha256)) {
  stop("Pinned MCP-counter signature SHA-256 mismatch: ", observed_sha256)
}
marker_genes <- read.table(signature_path, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = "character")

expression <- read_expression("data/processed/GSE44076_gene_expression_raw_cel_rma.parquet")
metadata <- read.csv("data/processed/GSE44076_sample_metadata_raw_cel_rma.csv",
                     stringsAsFactors = FALSE)
metadata <- metadata[match(colnames(expression), metadata$geo_accession), , drop = FALSE]
if (anyNA(metadata$geo_accession)) stop("Composition metadata alignment failed")
metadata$age <- suppressWarnings(as.numeric(metadata$age))
metadata$sex <- factor(metadata$sex)
metadata$location <- factor(metadata$location)

estimates <- MCPcounter::MCPcounter.estimate(
  expression,
  featuresType = "HUGO_symbols",
  genes = marker_genes
)
if (identical(colnames(estimates), colnames(expression))) {
  score_matrix <- t(estimates)
} else if (identical(rownames(estimates), colnames(expression))) {
  score_matrix <- estimates
} else {
  stop("MCP-counter result does not align to GSE44076 samples")
}
scores <- data.frame(geo_accession = rownames(score_matrix), score_matrix,
                     check.names = FALSE, stringsAsFactors = FALSE)
scores <- merge(scores,
                metadata[c("geo_accession", "patient_id", "donor_or_patient_group",
                           "tissue_class", "age", "sex", "location")],
                by = "geo_accession", sort = FALSE)
scores <- scores[match(colnames(expression), scores$geo_accession), , drop = FALSE]
score_columns <- colnames(score_matrix)

group_rows <- list()
for (population in score_columns) {
  field <- scores$tissue_class %in% c("healthy", "adjacent_normal")
  frame <- droplevels(scores[field, , drop = FALSE])
  frame$tissue <- relevel(factor(frame$tissue_class), ref = "healthy")
  design <- model.matrix(~ tissue + age + sex, data = frame)
  if (qr(design)$rank != ncol(design)) stop(population, ": composition design is rank deficient")
  fit <- lm.fit(design, frame[[population]])
  coefficient <- "tissueadjacent_normal"
  position <- match(coefficient, colnames(design))
  residual_df <- nrow(design) - qr(design)$rank
  sigma <- sqrt(sum(fit$residuals^2) / residual_df)
  covariance <- chol2inv(qr.R(fit$qr)) * sigma^2
  estimate <- fit$coefficients[[position]]
  standard_error <- sqrt(covariance[position, position])
  statistic <- estimate / standard_error
  p_value <- 2 * pt(abs(statistic), df = residual_df, lower.tail = FALSE)
  group_rows[[population]] <- data.frame(
    population = population,
    method = "MCP-counter 1.2.0; age/sex-adjusted linear model",
    healthy_n = sum(frame$tissue_class == "healthy"),
    adjacent_normal_n = sum(frame$tissue_class == "adjacent_normal"),
    healthy_mean = mean(frame[[population]][frame$tissue_class == "healthy"]),
    adjacent_normal_mean = mean(frame[[population]][frame$tissue_class == "adjacent_normal"]),
    adjusted_difference = estimate,
    standard_error = standard_error,
    statistic = statistic,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}
group_comparison <- do.call(rbind, group_rows)
group_comparison$adjusted_p_value <- p.adjust(group_comparison$p_value, method = "BH")
rownames(group_comparison) <- NULL
write.csv(scores, "results/tables/tissue_composition_scores.csv", row.names = FALSE)
write.csv(group_comparison, "results/tables/tissue_composition_group_comparison.csv",
          row.names = FALSE)

long_scores <- reshape(
  scores[c("geo_accession", "tissue_class", score_columns)],
  varying = score_columns, v.names = "score", timevar = "population",
  times = score_columns, direction = "long"
)
long_scores$population <- factor(long_scores$population, levels = score_columns)
composition_plot <- ggplot(long_scores,
                           aes(x = tissue_class, y = score, fill = tissue_class)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.45) +
  facet_wrap(~ population, scales = "free_y", ncol = 3) +
  labs(title = "MCP-counter composition scores by GSE44076 tissue group",
       x = NULL, y = "MCP-counter score") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
save_plot(composition_plot, "results/figures/tissue_composition_by_group")

field <- metadata$tissue_class %in% c("healthy", "adjacent_normal") &
  complete.cases(metadata[c("age", "sex")])
field_meta <- droplevels(metadata[field, , drop = FALSE])
field_scores <- score_matrix[field_meta$geo_accession, , drop = FALSE]
score_pca <- prcomp(field_scores, center = TRUE, scale. = TRUE)
variance <- score_pca$sdev^2 / sum(score_pca$sdev^2)
pc_count <- 2L
pc_scores <- score_pca$x[, seq_len(pc_count), drop = FALSE]
colnames(pc_scores) <- paste0("composition_pc", seq_len(pc_count))
field_meta <- cbind(field_meta, pc_scores[field_meta$geo_accession, , drop = FALSE])
write.csv(data.frame(population = rownames(score_pca$rotation), score_pca$rotation,
                     check.names = FALSE),
          "results/tables/tissue_composition_pca_loadings.csv", row.names = FALSE)

field_meta$tissue <- relevel(factor(field_meta$tissue_class), ref = "healthy")
design <- model.matrix(~ tissue + age + sex + composition_pc1 + composition_pc2,
                       data = field_meta)
if (qr(design)$rank != ncol(design)) stop("Composition-adjusted limma design is rank deficient")
fit <- eBayes(lmFit(expression[, field_meta$geo_accession, drop = FALSE], design))
coefficient <- "tissueadjacent_normal"
tab <- topTable(fit, coef = coefficient, number = Inf, sort.by = "P",
                adjust.method = "BH", confint = 0.95)
standard_error <- fit$stdev.unscaled[, coefficient] * fit$sigma
mapping <- read.csv("data/metadata/GSE44076_probe_gene_mapping_raw_cel_rma.csv",
                    stringsAsFactors = FALSE)
valid_mapping <- mapping[mapping$mapping_status == "mapped_unique", ]
gene_to_entrez <- tapply(valid_mapping$entrez_id, valid_mapping$gene_symbol,
                         function(x) sort(unique(x[x != ""]))[[1]])
composition_de <- data.frame(
  gene_symbol = rownames(tab),
  entrez_id = unname(gene_to_entrez[rownames(tab)]),
  model = "U3_age_sex_composition_pc_adjusted",
  contrast = "adjacent_vs_healthy",
  comparison = "adjacent_vs_healthy",
  log2_fold_change = tab$logFC,
  standard_error = standard_error[rownames(tab)],
  moderated_statistic = tab$t,
  p_value = tab$P.Value,
  adjusted_p_value = tab$adj.P.Val,
  average_expression = tab$AveExpr,
  confidence_interval_lower = tab$CI.L,
  confidence_interval_upper = tab$CI.R,
  direction = ifelse(tab$logFC >= 0, "up", "down"),
  mapping_status = "mapped_unique",
  analysis_provenance = paste0(
    "raw_cel_rma_limma_MCPcounter_PC1_PC2; variance_explained=",
    paste(round(variance[seq_len(pc_count)], 4), collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(composition_de,
          "results/differential_expression/raw_cel_adjacent_vs_healthy_age_sex_composition_adjusted.csv",
          row.names = FALSE)
write.csv(composition_de, "results/tables/tissue_composition_sensitivity.csv",
          row.names = FALSE)

tier_path <- "results/tables/field_gene_evidence_tiers.csv"
tiers <- read.csv(tier_path, stringsAsFactors = FALSE)
composition_lookup <- composition_de[match(tiers$gene_symbol, composition_de$gene_symbol), ]
tiers$composition_adjusted_log2fc <- composition_lookup$log2_fold_change
tiers$composition_adjusted_fdr <- composition_lookup$adjusted_p_value
tiers$composition_robust <- tiers$composition_adjusted_fdr < 0.05 &
  abs(tiers$composition_adjusted_log2fc) >= 0.5 &
  sign(tiers$composition_adjusted_log2fc) == sign(tiers$u1_log2fc)
write.csv(tiers, tier_path, row.names = FALSE)
high <- tiers$evidence_tier == "Tier 1 - high-confidence"
high_table <- tiers[high, , drop = FALSE]
write.csv(high_table, "results/tables/final_high_confidence_field_signature.csv",
          row.names = FALSE)
write.csv(tiers[tiers$evidence_tier == "Tier 2 - provisional", , drop = FALSE],
          "results/tables/provisional_field_associated_genes.csv", row.names = FALSE)
write.csv(tiers[tiers$evidence_tier == "Tier 3 - exploratory", , drop = FALSE],
          "results/tables/exploratory_field_gene_universe.csv", row.names = FALSE)

stromal_markers <- unique(marker_genes[marker_genes[["Cell population"]] %in%
                                        c("Fibroblasts", "Endothelial cells"),
                                      "HUGO symbols"])
stromal_high <- high_table$gene_symbol %in% stromal_markers |
  grepl("^(COL[0-9]|ECM|MMP|TIMP|THBS|FN1$)", high_table$gene_symbol)
attenuated <- high_table$u1_log2fc != 0 &
  abs(high_table$composition_adjusted_log2fc) < 0.75 * abs(high_table$u1_log2fc)
report <- c(
  "# Tissue-composition sensitivity", "",
  "MCP-counter 1.2.0 was run on the independently RMA-normalized GSE44076",
  "gene-symbol matrix. The method's official HUGO marker table was pinned to",
  paste0("commit `b6eac73e` and SHA-256 `", signature_sha256, "`. MCP-counter"),
  "was developed and validated using transcriptomic data that included Affymetrix",
  "microarrays. Scores are comparable across samples within this cohort; they are",
  "not literal cell fractions.", "",
  paste(capture.output(print(group_comparison, row.names = FALSE)), collapse = "\n"), "",
  sprintf("PC1 and PC2 explained %.1f%% and %.1f%% of standardized score variance.",
          100 * variance[[1]], 100 * variance[[2]]),
  paste("The composition-adjusted sensitivity model included tissue, age, sex,",
        "and only these two prespecified composition PCs."),
  sprintf("Of %d high-confidence field genes, %d (%.1f%%) retained FDR < 0.05,",
          nrow(high_table), sum(high_table$composition_robust, na.rm = TRUE),
          100 * mean(high_table$composition_robust, na.rm = TRUE)),
  "|log2FC| >= 0.5, and direction after composition adjustment.",
  sprintf("Among %d stromal/ECM-flagged high-confidence genes, %d attenuated by at least 25%%.",
          sum(stromal_high), sum(stromal_high & attenuated, na.rm = TRUE)), "",
  paste("Attenuation does not show that the original association was false:",
        "composition can be part of the field microenvironment while confounding",
        "a strictly epithelial interpretation. Because score genes are derived from",
        "the same bulk transcriptome, this is an overadjustment-prone sensitivity",
        "analysis rather than a causal decomposition.")
)
writeLines(report, "reports/tissue_composition_sensitivity.md")
