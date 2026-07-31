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
if (!all(de$analysis_provenance %in% c("raw_cel_rma_limma",
                                      "raw_cel_rma_limma_patient_fixed_effect"))) {
  stop("Enrichment requires raw-CEL limma inputs")
}
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
universe <- unique(as.character(de$entrez_id[de$entrez_id != "" & !is.na(de$entrez_id)]))

sets <- list(
  adjacent_vs_healthy_up = subset(de, comparison == "adjacent_vs_healthy" &
                                       adjusted_p_value < 0.05 & log2_fold_change >= 0.5),
  adjacent_vs_healthy_down = subset(de, comparison == "adjacent_vs_healthy" &
                                         adjusted_p_value < 0.05 & log2_fold_change <= -0.5),
  tumor_vs_adjacent_up = subset(de, comparison == "tumor_vs_adjacent_paired" &
                                     adjusted_p_value < 0.05 & log2_fold_change >= 0.5),
  tumor_vs_adjacent_down = subset(de, comparison == "tumor_vs_adjacent_paired" &
                                       adjusted_p_value < 0.05 & log2_fold_change <= -0.5)
)

save_plot <- function(plot, stem, width = 8, height = 6) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300,
         bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, bg = "white")
}

all_results <- list()
for (name in names(sets)) {
  ids <- unique(as.character(sets[[name]]$entrez_id))
  ids <- ids[ids != "" & !is.na(ids)]
  if (length(ids) < 10L) {
    warning(name, ": fewer than 10 mapped genes; over-representation enrichment skipped")
    next
  }
  ego <- enrichGO(ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                  universe = universe, pAdjustMethod = "BH", readable = TRUE)
  reactome <- enrichPathway(ids, universe = universe, pAdjustMethod = "BH", readable = TRUE)
  go_frame <- as.data.frame(ego)
  reactome_frame <- as.data.frame(reactome)
  if (nrow(go_frame)) {
    go_frame$gene_set <- name
    go_frame$database <- "GO_BP"
    all_results[[paste0(name, "_go")]] <- go_frame
    write.csv(go_frame, file.path("results/enrichment", paste0(name, "_go_bp.csv")),
              row.names = FALSE)
    save_plot(dotplot(ego, showCategory = 20) + ggtitle(paste(name, "GO BP")),
              file.path("results/figures", paste0(name, "_go_bp_dotplot")))
  }
  if (nrow(reactome_frame)) {
    reactome_frame$gene_set <- name
    reactome_frame$database <- "Reactome"
    all_results[[paste0(name, "_reactome")]] <- reactome_frame
    write.csv(reactome_frame,
              file.path("results/enrichment", paste0(name, "_reactome.csv")), row.names = FALSE)
    save_plot(dotplot(reactome, showCategory = 20) + ggtitle(paste(name, "Reactome")),
              file.path("results/figures", paste0(name, "_reactome_dotplot")))
  }
}

for (comparison in unique(de$comparison)) {
  frame <- de[de$comparison == comparison & !is.na(de$entrez_id) & de$entrez_id != "", ]
  ranks <- frame$moderated_statistic
  names(ranks) <- frame$entrez_id
  collapsed <- tapply(ranks, names(ranks), function(x) x[which.max(abs(x))])
  ranks <- as.numeric(collapsed)
  names(ranks) <- names(collapsed)
  ranks <- sort(ranks, decreasing = TRUE)
  gsea <- gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                pAdjustMethod = "BH", verbose = FALSE)
  gsea_frame <- as.data.frame(gsea)
  if (nrow(gsea_frame)) {
    gsea_frame$gene_set <- comparison
    gsea_frame$database <- "GO_BP_ranked"
    all_results[[paste0(comparison, "_gsea")]] <- gsea_frame
    write.csv(gsea_frame,
              file.path("results/enrichment", paste0(comparison, "_ranked_go_bp.csv")),
              row.names = FALSE)
    save_plot(dotplot(gsea, showCategory = 20) + ggtitle(paste(comparison, "ranked GO BP")),
              file.path("results/figures", paste0(comparison, "_ranked_go_bp")))
  }
}

if (length(all_results)) {
  columns <- unique(unlist(lapply(all_results, names)))
  padded <- lapply(all_results, function(frame) {
    missing <- setdiff(columns, names(frame))
    for (column in missing) frame[[column]] <- NA
    frame[, columns]
  })
  write.csv(do.call(rbind, padded), "results/tables/supplementary_enrichment_results.csv",
            row.names = FALSE)
} else {
  write.csv(data.frame(note = "No enrichment result passed execution criteria"),
            "results/tables/supplementary_enrichment_results.csv", row.names = FALSE)
}
