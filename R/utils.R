suppressPackageStartupMessages({
  library(Biobase)
  library(AnnotationDbi)
  library(matrixStats)
})

project_root <- function() {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

discover_files <- function(accession, pattern, root = "data/raw") {
  candidates <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE,
                           ignore.case = TRUE)
  candidates[grepl(accession, candidates, ignore.case = TRUE)]
}

assert_platform <- function(files, expected) {
  if (!length(files)) stop("No CEL files were discovered")
  message("Expected platform: ", expected, "; CEL files: ", length(files))
}

write_expression <- function(mat, path, id_column = "gene_symbol") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- data.frame(row_identifier = rownames(mat), mat, check.names = FALSE)
  names(out)[1] <- id_column
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(out, path)
  } else {
    fallback <- sub("\\.parquet$", ".csv.gz", path)
    data.table::fwrite(out, fallback)
    warning("R package 'arrow' is unavailable; wrote ", fallback,
            ". Run scripts/convert_r_outputs.py after installing pyarrow.")
  }
}

unambiguous_mapping <- function(probe_ids, annotation_package) {
  map <- AnnotationDbi::select(
    annotation_package, keys = unique(probe_ids),
    keytype = "PROBEID", columns = c("SYMBOL", "ENTREZID")
  )
  map <- map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), ]
  symbol_count <- aggregate(SYMBOL ~ PROBEID, map, function(x) length(unique(x)))
  valid <- symbol_count$PROBEID[symbol_count$SYMBOL == 1]
  map <- map[map$PROBEID %in% valid, ]
  map <- map[!duplicated(map[c("PROBEID", "SYMBOL")]), ]
  map
}

aggregate_probes_by_gene <- function(expr, mapping) {
  mapping <- mapping[match(rownames(expr), mapping$PROBEID, nomatch = 0), ]
  expr <- expr[rownames(expr) %in% mapping$PROBEID, , drop = FALSE]
  mapping <- mapping[match(rownames(expr), mapping$PROBEID), ]
  groups <- split(seq_len(nrow(expr)), mapping$SYMBOL)
  gene <- t(vapply(groups, function(i) matrixStats::colMedians(expr[i, , drop = FALSE]),
                   numeric(ncol(expr))))
  colnames(gene) <- colnames(expr)
  gene
}

save_basic_qc <- function(expr, metadata, prefix) {
  dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
  png(file.path("results/figures", paste0(prefix, "_normalized_boxplot.png")),
      width = 2400, height = 1600, res = 300)
  boxplot(expr, outline = FALSE, las = 2, cex.axis = 0.3,
          main = paste(prefix, "normalized expression"), ylab = "log2 expression")
  dev.off()

  pc <- prcomp(t(expr), scale. = FALSE)
  cls <- factor(metadata$tissue_class)
  png(file.path("results/figures", paste0(prefix, "_normalized_pca.png")),
      width = 2100, height = 1800, res = 300)
  plot(pc$x[, 1], pc$x[, 2], col = as.integer(cls), pch = 19,
       xlab = "PC1", ylab = "PC2", main = paste(prefix, "normalized PCA"))
  legend("topright", levels(cls), col = seq_along(levels(cls)), pch = 19)
  dev.off()
}
