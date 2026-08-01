required <- c(
  "affy", "affyio", "oligo", "pd.hg.u219", "hgu219.db", "hgu133acdf",
  "hgu133a.db", "AnnotationDbi", "org.Hs.eg.db", "Biobase", "limma",
  "arrayQualityMetrics", "affyPLM", "MCPcounter", "GEOquery", "matrixStats", "data.table", "ggplot2",
  "pheatmap", "clusterProfiler", "ReactomePA", "fgsea", "enrichplot", "renv"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))
if (getRversion() < "4.5.0" || getRversion() >= "4.6.0") stop("Expected R 4.5.x")
cat("R environment verification passed\n")
