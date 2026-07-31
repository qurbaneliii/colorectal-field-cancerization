suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

de_path <- "results/tables/supplementary_full_de_results.csv"
if (!file.exists(de_path)) stop("Run R/04_differential_expression.R first")
de <- read.csv(de_path, stringsAsFactors = FALSE)
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

universe <- unique(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = unique(de$gene_symbol), keytype = "SYMBOL",
  column = "ENTREZID", multiVals = "first"
))
universe <- universe[!is.na(universe)]

sets <- list(
  adjacent_vs_healthy_up = subset(de, comparison == "adjacent_normal_vs_healthy" &
                                       adj.P.Val < 0.05 & logFC >= 0.5)$gene_symbol,
  adjacent_vs_healthy_down = subset(de, comparison == "adjacent_normal_vs_healthy" &
                                         adj.P.Val < 0.05 & logFC <= -0.5)$gene_symbol,
  tumor_vs_adjacent = subset(de, comparison == "tumor_vs_adjacent_normal_paired" &
                                  adj.P.Val < 0.05 & abs(logFC) >= 0.5)$gene_symbol
)

for (name in names(sets)) {
  ids <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(sets[[name]]),
                               keytype = "SYMBOL", column = "ENTREZID",
                               multiVals = "first")
  ids <- unique(ids[!is.na(ids)])
  if (length(ids) < 10) {
    warning(name, ": fewer than 10 mapped genes; enrichment skipped")
    next
  }
  ego <- enrichGO(ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                  universe = universe, pAdjustMethod = "BH", readable = TRUE)
  write.csv(as.data.frame(ego), file.path("results/enrichment", paste0(name, "_go_bp.csv")),
            row.names = FALSE)
  pdf(file.path("results/figures", paste0(name, "_go_bp_dotplot.pdf")), 8, 6)
  print(dotplot(ego, showCategory = 20) + ggtitle(name))
  dev.off()
}
