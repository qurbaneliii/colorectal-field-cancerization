options(repos = c(CRAN = "https://cloud.r-project.org"))

if (getRversion() < "4.5.0" || getRversion() >= "4.6.0") {
  stop("This publication environment requires R 4.5.x; detected ", getRversion())
}
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!file.exists("renv/activate.R")) renv::init(bare = TRUE, restart = FALSE)
renv::load(project = getwd())
renv::settings$bioconductor.version("3.21")

if (!requireNamespace("BiocManager", quietly = TRUE)) renv::install("BiocManager")
if (as.character(BiocManager::version()) != "3.21") {
  BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
}

cran <- c(
  "yaml", "data.table", "matrixStats", "ggplot2", "pheatmap", "arrow",
  "msigdbr", "rmarkdown", "knitr", "jsonlite", "digest"
)
bioc <- c(
  "affy", "affyio", "affyPLM", "oligo", "pd.hg.u219", "hgu219.db",
  "hgu133acdf", "hgu133a.db", "AnnotationDbi", "org.Hs.eg.db", "Biobase",
  "limma", "arrayQualityMetrics", "GEOquery", "clusterProfiler", "ReactomePA",
  "fgsea", "enrichplot", "GO.db", "reactome.db"
)

missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) renv::install(missing_cran)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  for (package in missing_bioc) {
    BiocManager::install(package, ask = FALSE, update = FALSE, Ncpus = 1L)
  }
}

still_missing <- c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) stop("Required R packages remain unavailable: ",
                                paste(still_missing, collapse = ", "))
renv::snapshot(type = "all", prompt = FALSE)
message("Resolved R/Bioconductor environment and updated renv.lock")
