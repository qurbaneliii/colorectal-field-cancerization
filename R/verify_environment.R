required <- c(
  "affy", "affyio", "oligo", "pd.hg.u219", "hgu219.db", "hgu133acdf",
  "hgu133a.db", "AnnotationDbi", "org.Hs.eg.db", "Biobase", "limma",
  "arrayQualityMetrics", "affyPLM", "MCPcounter", "GEOquery", "matrixStats", "data.table", "R.utils",
  "ggplot2", "svglite", "digest",
  "pheatmap", "clusterProfiler", "ReactomePA", "fgsea", "enrichplot", "renv"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))
lock <- renv::lockfile_read("renv.lock")
if (as.character(getRversion()) != lock$R$Version) {
  stop("Expected locked R version ", lock$R$Version, ", detected ", getRversion())
}
expected <- vapply(lock$Packages[required], function(record) record$Version, character(1))
observed <- vapply(required, function(package) packageDescription(package, fields = "Version"), character(1))
different <- observed != expected
if (any(different)) {
  stop("R package version mismatch: ", paste(
    sprintf("%s (locked %s, installed %s)", required[different],
            expected[different], observed[different]), collapse = ", "
  ))
}
cat("R and required packages match renv.lock\n")
