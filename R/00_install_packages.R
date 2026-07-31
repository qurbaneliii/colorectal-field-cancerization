options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("renv", "yaml", "data.table", "matrixStats", "ggplot2", "pheatmap")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c(
  "oligo", "pd.hg.u219", "hgu219.db", "affy", "hgu133acdf", "hgu133a.db",
  "AnnotationDbi", "org.Hs.eg.db", "Biobase", "limma", "arrayQualityMetrics",
  "GEOquery", "clusterProfiler", "ReactomePA"
)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

if (!requireNamespace("renv", quietly = TRUE)) stop("renv installation failed")
renv::snapshot(prompt = FALSE)
