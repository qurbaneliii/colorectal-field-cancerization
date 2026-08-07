suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
})

de_path <- "results/tables/supplementary_full_raw_cel_de_results.csv"
if (!file.exists(de_path)) stop("Run R/04_differential_expression.R first")
de <- read.csv(de_path, stringsAsFactors = FALSE)
if (!all(grepl("^raw_cel_rma_limma", de$analysis_provenance))) {
  stop("Enrichment requires raw-CEL limma inputs")
}
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

mapping <- unique(de[, c("gene_symbol", "entrez_id")])
mapping <- mapping[!is.na(mapping$entrez_id) & mapping$entrez_id != "", , drop = FALSE]
mapping$entrez_id <- as.character(mapping$entrez_id)
universe <- unique(mapping$entrez_id)

read_genes <- function(path, condition = NULL) {
  if (!file.exists(path)) return(character(0))
  frame <- read.csv(path, stringsAsFactors = FALSE)
  if (!is.null(condition)) frame <- condition(frame)
  unique(as.character(frame$gene_symbol))
}

high <- read.csv("results/tables/final_high_confidence_field_signature.csv",
                 stringsAsFactors = FALSE)
provisional <- read.csv("results/tables/provisional_field_associated_genes.csv",
                        stringsAsFactors = FALSE)
task_b <- read_genes("results/tables/task_b_final_signature.csv")
task_c <- read_genes("results/tables/final_compact_signature.csv")
stress_removed <- read_genes(
  "results/tables/high_confidence_field_signature_without_stress_genes.csv"
)

sets <- list(
  high_confidence_field_up = list(
    genes = high$gene_symbol[high$u1_log2fc > 0], comparison = "adjacent_vs_healthy_u1",
    direction = "up", tier = "high_confidence", minimum = 10L
  ),
  high_confidence_field_down = list(
    genes = high$gene_symbol[high$u1_log2fc < 0], comparison = "adjacent_vs_healthy_u1",
    direction = "down", tier = "high_confidence", minimum = 10L
  ),
  provisional_field_up = list(
    genes = provisional$gene_symbol[provisional$u1_log2fc > 0],
    comparison = "adjacent_vs_healthy_u1", direction = "up", tier = "provisional",
    minimum = 10L
  ),
  provisional_field_down = list(
    genes = provisional$gene_symbol[provisional$u1_log2fc < 0],
    comparison = "adjacent_vs_healthy_u1", direction = "down", tier = "provisional",
    minimum = 10L
  ),
  composition_robust_field = list(
    genes = high$gene_symbol[high$composition_robust %in% TRUE],
    comparison = "adjacent_vs_healthy_u3_composition_adjusted", direction = "both",
    tier = "composition_robust_high_confidence", minimum = 10L
  ),
  composition_sensitive_field = list(
    genes = high$gene_symbol[!(high$composition_robust %in% TRUE)],
    comparison = "adjacent_vs_healthy_u3_composition_adjusted", direction = "both",
    tier = "composition_sensitive_high_confidence", minimum = 10L
  ),
  high_confidence_without_stress = list(
    genes = stress_removed, comparison = "adjacent_vs_healthy_u1", direction = "both",
    tier = "high_confidence_stress_removed", minimum = 10L
  ),
  task_b_locked_panel = list(
    genes = task_b, comparison = "healthy_vs_adjacent_predictive", direction = "model",
    tier = "locked_model_panel", minimum = 2L
  ),
  task_c_locked_panel = list(
    genes = task_c, comparison = "tumor_vs_adjacent_predictive", direction = "model",
    tier = "locked_model_panel", minimum = 2L
  )
)

qc_de <- de[de$model == "U1_age_sex_adjusted_qc_sensitivity", , drop = FALSE]
if (nrow(qc_de)) {
  sets$qc_excluded_field_up <- list(
    genes = qc_de$gene_symbol[qc_de$adjusted_p_value < 0.05 &
                                qc_de$log2_fold_change >= 0.5],
    comparison = "adjacent_vs_healthy_qc_excluded", direction = "up",
    tier = "qc_sensitivity", minimum = 10L
  )
  sets$qc_excluded_field_down <- list(
    genes = qc_de$gene_symbol[qc_de$adjusted_p_value < 0.05 &
                                qc_de$log2_fold_change <= -0.5],
    comparison = "adjacent_vs_healthy_qc_excluded", direction = "down",
    tier = "qc_sensitivity", minimum = 10L
  )
}

save_plot <- function(plot, stem, width = 8, height = 6) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300,
         bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, bg = "white")
}

standardize <- function(frame, name, definition, database, ontology,
                        redundancy_pruned, analysis_type = "over_representation") {
  if (!nrow(frame)) return(NULL)
  required <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
                "qvalue", "geneID", "Count", "enrichmentScore", "NES", "setSize")
  for (column in setdiff(required, names(frame))) frame[[column]] <- NA
  data.frame(
    analysis_set = name,
    comparison = definition$comparison,
    direction = definition$direction,
    evidence_tier = definition$tier,
    analysis_type = analysis_type,
    database = database,
    ontology = ontology,
    ID = frame$ID,
    Description = frame$Description,
    GeneRatio = frame$GeneRatio,
    BgRatio = frame$BgRatio,
    p_value = frame$pvalue,
    adjusted_p_value = frame$p.adjust,
    q_value = frame$qvalue,
    gene_ids = frame$geneID,
    gene_count = frame$Count,
    enrichment_score = frame$enrichmentScore,
    normalized_enrichment_score = frame$NES,
    set_size = frame$setSize,
    redundancy_pruned = redundancy_pruned,
    interpretation = ifelse(
      definition$tier == "locked_model_panel",
      "exploratory small-panel over-representation; not mechanistic evidence",
      "association-level pathway result; not causal or cell-intrinsic evidence"
    ),
    analysis_provenance = "raw_cel_rma_tier_specific_enrichment",
    stringsAsFactors = FALSE
  )
}

all_results <- list()
set_audit <- list()
ranked_only <- identical(Sys.getenv("RANKED_ENRICHMENT_ONLY"), "1")
if (ranked_only) {
  existing <- read.csv("results/tables/supplementary_enrichment_results.csv",
                       stringsAsFactors = FALSE)
  existing <- existing[existing$analysis_type != "ranked_GSEA", , drop = FALSE]
  all_results$existing_thresholded <- existing
} else for (name in names(sets)) {
  definition <- sets[[name]]
  symbols <- unique(as.character(definition$genes))
  ids <- unique(mapping$entrez_id[mapping$gene_symbol %in% symbols])
  set_audit[[name]] <- data.frame(
    analysis_set = name, input_symbols = length(symbols), mapped_entrez = length(ids),
    minimum_required = definition$minimum,
    executed = length(ids) >= definition$minimum,
    stringsAsFactors = FALSE
  )
  if (length(ids) < definition$minimum) {
    warning(name, ": fewer than ", definition$minimum,
            " mapped genes; over-representation enrichment skipped")
    next
  }
  ego <- enrichGO(
    ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
    universe = universe, pAdjustMethod = "BH", readable = TRUE
  )
  ego_pruned <- tryCatch(
    simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(error) ego
  )
  reactome <- enrichPathway(
    ids, universe = universe, pAdjustMethod = "BH", readable = TRUE
  )
  go_frame <- as.data.frame(ego_pruned)
  reactome_frame <- as.data.frame(reactome)
  standardized_go <- standardize(
    go_frame, name, definition, "GO", "Biological Process", TRUE
  )
  standardized_reactome <- standardize(
    reactome_frame, name, definition, "Reactome", "Pathway", FALSE
  )
  if (!is.null(standardized_go)) {
    all_results[[paste0(name, "_go")]] <- standardized_go
    write.csv(standardized_go,
              file.path("results/enrichment", paste0(name, "_go_bp_pruned.csv")),
              row.names = FALSE)
    save_plot(
      dotplot(ego_pruned, showCategory = 15) + ggtitle(paste(name, "GO BP (pruned)")),
      file.path("results/figures", paste0(name, "_go_bp_dotplot"))
    )
  }
  if (!is.null(standardized_reactome)) {
    all_results[[paste0(name, "_reactome")]] <- standardized_reactome
    write.csv(standardized_reactome,
              file.path("results/enrichment", paste0(name, "_reactome.csv")),
              row.names = FALSE)
    save_plot(
      dotplot(reactome, showCategory = 15) + ggtitle(paste(name, "Reactome")),
      file.path("results/figures", paste0(name, "_reactome_dotplot"))
    )
  }
}

ranked_definitions <- list(
  adjacent_vs_healthy_u1 = list(
    comparison = "adjacent_vs_healthy_u1", direction = "ranked",
    tier = "full_rank_covariate_adjusted", source_comparison = "adjacent_vs_healthy",
    source_model = "U1_age_sex_adjusted"
  ),
  tumor_vs_adjacent_paired = list(
    comparison = "tumor_vs_adjacent_paired", direction = "ranked",
    tier = "full_rank_paired", source_comparison = "tumor_vs_adjacent",
    source_model = "P_patient_fixed_effect"
  )
)
for (comparison in names(ranked_definitions)) {
  definition <- ranked_definitions[[comparison]]
  frame <- de[de$comparison == definition$source_comparison &
                de$model == definition$source_model & !is.na(de$entrez_id) &
                de$entrez_id != "", , drop = FALSE]
  if (!nrow(frame)) next
  ranks <- frame$moderated_statistic
  names(ranks) <- frame$entrez_id
  collapsed <- tapply(ranks, names(ranks), function(x) x[which.max(abs(x))])
  ranks <- sort(as.numeric(collapsed), decreasing = TRUE)
  names(ranks) <- names(sort(collapsed, decreasing = TRUE))
  gsea <- gseGO(
    ranks, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
    pAdjustMethod = "BH", verbose = FALSE
  )
  gsea_frame <- as.data.frame(gsea)
  standardized <- standardize(
    gsea_frame, comparison, definition, "GO",
    "Biological Process", FALSE, "ranked_GSEA"
  )
  if (!is.null(standardized)) {
    all_results[[paste0(comparison, "_ranked")]] <- standardized
    write.csv(standardized,
              file.path("results/enrichment", paste0(comparison, "_ranked_go_bp.csv")),
              row.names = FALSE)
    save_plot(
      dotplot(gsea, showCategory = 15) + ggtitle(paste(comparison, "ranked GO BP")),
      file.path("results/figures", paste0(comparison, "_ranked_go_bp"))
    )
  }
}

audit <- if (ranked_only) {
  read.csv("results/tables/enrichment_analysis_set_audit.csv", stringsAsFactors = FALSE)
} else {
  do.call(rbind, set_audit)
}
if (!ranked_only) {
  write.csv(audit, "results/tables/enrichment_analysis_set_audit.csv", row.names = FALSE)
}
if (length(all_results)) {
  combined <- do.call(rbind, all_results)
  rownames(combined) <- NULL
  combined <- combined[order(combined$analysis_set, combined$adjusted_p_value), ]
  write.csv(combined, "results/tables/supplementary_enrichment_results.csv",
            row.names = FALSE)
  significant <- combined[is.finite(combined$adjusted_p_value) &
                            combined$adjusted_p_value < 0.05, , drop = FALSE]
  top <- do.call(rbind, lapply(split(significant, significant$analysis_set), head, 5L))
  if (is.null(top)) top <- significant
  write.csv(top, "results/tables/enrichment_top_terms_by_analysis_set.csv",
            row.names = FALSE)
} else {
  combined <- data.frame()
  top <- data.frame()
  write.csv(data.frame(note = "No enrichment result passed execution criteria"),
            "results/tables/supplementary_enrichment_results.csv", row.names = FALSE)
}

top_text <- if (nrow(top)) paste(capture.output(print(
  top[, c("analysis_set", "database", "Description", "adjusted_p_value")],
  row.names = FALSE)), collapse = "\n") else "No FDR-significant terms."
report <- paste0(
  "# Tier-specific enrichment report\n\n",
  "Enrichment was rerun separately for high-confidence up/down field genes, ",
  "provisional genes, composition-robust and composition-sensitive subsets, ",
  "stress-gene-removed genes, locked Task B/Task C panels, and full ranked ",
  "contrasts. GO Biological Process over-representation results were pruned ",
  "with semantic similarity (`simplify`, cutoff 0.7). Model-panel analyses ",
  "are explicitly exploratory because the sets contain only a few genes.\n\n",
  "All terms are associations and may reflect tissue composition; they are ",
  "not evidence of causal or epithelial-intrinsic mechanisms.\n\n",
  "## Analysis-set audit\n\n```\n",
  paste(capture.output(print(audit, row.names = FALSE)), collapse = "\n"),
  "\n```\n\n## Top FDR-significant terms\n\n```\n", top_text, "\n```\n"
)
writeLines(report, "reports/tier_specific_enrichment_report.md")
