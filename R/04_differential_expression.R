suppressPackageStartupMessages({
  library(limma)
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

expr <- read_expression("data/processed/GSE44076_gene_expression.parquet")
meta <- read.csv("data/processed/GSE44076_sample_metadata.csv", stringsAsFactors = FALSE)
meta <- meta[match(colnames(expr), meta$geo_accession), ]
stopifnot(!anyNA(meta$geo_accession))
dir.create("results/differential_expression", recursive = TRUE, showWarnings = FALSE)

run_unpaired <- function(level_a, level_b, label) {
  keep <- meta$tissue_class %in% c(level_a, level_b)
  y <- factor(meta$tissue_class[keep], levels = c(level_b, level_a))
  design <- model.matrix(~ y)
  fit <- eBayes(lmFit(expr[, keep, drop = FALSE], design))
  tab <- topTable(fit, coef = "y" %+% level_a, number = Inf, sort.by = "P")
  tab$gene_symbol <- rownames(tab)
  tab$comparison <- label
  tab$direction <- ifelse(tab$logFC > 0, "up", "down")
  tab
}

`%+%` <- function(a, b) paste0(a, b)

adj_healthy <- run_unpaired("adjacent_normal", "healthy", "adjacent_normal_vs_healthy")
tumor_healthy <- run_unpaired("tumor", "healthy", "tumor_vs_healthy")

paired <- meta$tissue_class %in% c("adjacent_normal", "tumor") & meta$pairing_status == "paired"
pair_meta <- meta[paired, ]
pair_expr <- expr[, paired, drop = FALSE]
pair_meta$patient_id <- factor(pair_meta$patient_id)
pair_meta$tissue_class <- factor(pair_meta$tissue_class,
                                 levels = c("adjacent_normal", "tumor"))
design <- model.matrix(~ patient_id + tissue_class, pair_meta)
fit <- eBayes(lmFit(pair_expr, design))
tumor_adjacent <- topTable(fit, coef = "tissue_classtumor", number = Inf, sort.by = "P")
tumor_adjacent$gene_symbol <- rownames(tumor_adjacent)
tumor_adjacent$comparison <- "tumor_vs_adjacent_normal_paired"
tumor_adjacent$direction <- ifelse(tumor_adjacent$logFC > 0, "up", "down")

all_de <- rbind(adj_healthy, tumor_healthy, tumor_adjacent)
all_de$annotation_status <- "mapped_unique_gene_symbol"
write.csv(adj_healthy, "results/differential_expression/adjacent_normal_vs_healthy.csv",
          row.names = FALSE)
write.csv(tumor_healthy, "results/differential_expression/tumor_vs_healthy.csv",
          row.names = FALSE)
write.csv(tumor_adjacent,
          "results/differential_expression/tumor_vs_adjacent_normal_paired.csv",
          row.names = FALSE)
write.csv(all_de, "results/tables/supplementary_full_de_results.csv", row.names = FALSE)
