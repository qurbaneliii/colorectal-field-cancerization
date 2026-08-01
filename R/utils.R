suppressPackageStartupMessages({
  library(Biobase)
  library(AnnotationDbi)
  library(matrixStats)
})

project_root <- function() normalizePath(getwd(), winslash = "/", mustWork = TRUE)

discover_files <- function(accession, pattern, root = "data/raw") {
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  candidates <- candidates[file.info(candidates)$isdir %in% FALSE]
  candidates[grepl(accession, basename(candidates), ignore.case = TRUE) &
               grepl(pattern, basename(candidates), ignore.case = TRUE, perl = TRUE)]
}

extract_gsm <- function(paths) {
  matches <- regexpr("GSM[0-9]+", basename(paths), ignore.case = TRUE)
  if (any(matches < 0L)) stop("CEL filenames without a GSM identifier: ",
                              paste(basename(paths[matches < 0L]), collapse = ", "))
  toupper(regmatches(basename(paths), matches))
}

normalized_token <- function(value) tolower(gsub("[^a-zA-Z0-9]", "", value))

validate_platform <- function(files, accession, expected_geo_platform, expected_chip,
                              annotation_package, reader) {
  if (!length(files)) stop(accession, ": no CEL files were discovered")
  if (anyDuplicated(tolower(basename(files)))) stop(accession, ": duplicate CEL filenames")
  expected_token <- normalized_token(expected_chip)
  if (reader == "oligo") {
    first <- oligo::read.celfiles(files[[1]], verbose = FALSE)
    detected_annotation <- as.character(Biobase::annotation(first))
    detected_cdf <- tryCatch(class(oligo::getCdfInfo(first))[[1]],
                             error = function(e) "not_applicable")
  } else if (reader == "affy") {
    first <- affy::ReadAffy(filenames = files[[1]])
    detected_annotation <- as.character(Biobase::annotation(first))
    detected_cdf <- as.character(affy::cdfName(first))
  } else {
    stop("Unknown CEL reader: ", reader)
  }
  observed <- normalized_token(paste(detected_annotation, detected_cdf))
  platform_pass <- grepl(expected_token, observed, fixed = TRUE)
  if (!platform_pass) {
    stop(accession, ": detected CEL annotation/CDF '", detected_annotation, " / ",
         detected_cdf, "' does not match expected chip ", expected_chip)
  }
  data.frame(
    accession = accession,
    expected_geo_platform = expected_geo_platform,
    expected_chip = expected_chip,
    detected_annotation = detected_annotation,
    detected_cdf = detected_cdf,
    annotation_package = annotation_package,
    reader = reader,
    cel_count = length(files),
    first_cel = basename(files[[1]]),
    platform_match = platform_pass,
    stringsAsFactors = FALSE
  )
}

validate_cel_metadata_alignment <- function(files, metadata) {
  gsm <- extract_gsm(files)
  if (anyDuplicated(gsm)) stop("Duplicate GSM identifiers parsed from CEL filenames")
  expected <- toupper(metadata$geo_accession)
  missing_cel <- setdiff(expected, gsm)
  unexpected_cel <- setdiff(gsm, expected)
  if (length(missing_cel) || length(unexpected_cel)) {
    stop("CEL/metadata GSM mismatch. Missing CEL: ", paste(head(missing_cel), collapse = ", "),
         "; unexpected CEL: ", paste(head(unexpected_cel), collapse = ", "))
  }
  metadata[match(gsm, expected), , drop = FALSE]
}

write_expression <- function(mat, path, id_column = "gene_symbol", accession,
                             provenance = "raw_cel_rma") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- data.frame(row_identifier = rownames(mat), mat, check.names = FALSE)
  names(out)[1] <- id_column
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(out, path)
    written_path <- path
  } else {
    written_path <- sub("\\.parquet$", ".csv.gz", path)
    data.table::fwrite(out, written_path)
    warning("R package 'arrow' is unavailable; wrote ", written_path)
  }
  sidecar <- data.frame(
    accession = accession,
    artifact = basename(written_path),
    provenance = provenance,
    rows = nrow(mat),
    samples = ncol(mat),
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
  write.csv(sidecar, paste0(written_path, ".provenance.csv"), row.names = FALSE)
}

complete_probe_mapping <- function(probe_ids, annotation_package) {
  raw <- AnnotationDbi::select(
    annotation_package, keys = unique(probe_ids), keytype = "PROBEID",
    columns = c("SYMBOL", "ENTREZID")
  )
  raw$SYMBOL[is.na(raw$SYMBOL)] <- ""
  raw$ENTREZID[is.na(raw$ENTREZID)] <- ""
  symbol_counts <- tapply(raw$SYMBOL[raw$SYMBOL != ""], raw$PROBEID[raw$SYMBOL != ""],
                          function(x) length(unique(x)))
  entrez_counts <- tapply(raw$ENTREZID[raw$ENTREZID != ""], raw$PROBEID[raw$ENTREZID != ""],
                          function(x) length(unique(x)))
  raw$symbol_count <- unname(symbol_counts[raw$PROBEID])
  raw$entrez_count <- unname(entrez_counts[raw$PROBEID])
  raw$symbol_count[is.na(raw$symbol_count)] <- 0L
  raw$entrez_count[is.na(raw$entrez_count)] <- 0L
  raw$mapping_status <- ifelse(
    raw$symbol_count == 0L | raw$entrez_count == 0L, "missing",
    ifelse(raw$symbol_count == 1L & raw$entrez_count == 1L,
           "mapped_unique", "ambiguous")
  )
  raw <- raw[!duplicated(raw[c("PROBEID", "SYMBOL", "ENTREZID")]), ]
  names(raw)[names(raw) == "PROBEID"] <- "probe_id"
  names(raw)[names(raw) == "SYMBOL"] <- "gene_symbol"
  names(raw)[names(raw) == "ENTREZID"] <- "entrez_id"
  raw[order(match(raw$probe_id, probe_ids), raw$gene_symbol), ]
}

aggregate_probes_by_gene <- function(expr, mapping) {
  valid <- mapping[mapping$mapping_status == "mapped_unique" &
                     mapping$gene_symbol != "" & mapping$entrez_id != "", ]
  valid <- valid[!duplicated(valid$probe_id), ]
  shared <- intersect(rownames(expr), valid$probe_id)
  if (!length(shared)) stop("No uniquely annotated probes overlap expression")
  expr <- expr[shared, , drop = FALSE]
  valid <- valid[match(shared, valid$probe_id), ]
  groups <- split(seq_len(nrow(expr)), valid$gene_symbol)
  gene <- t(vapply(groups, function(i) matrixStats::colMedians(expr[i, , drop = FALSE]),
                   numeric(ncol(expr))))
  colnames(gene) <- colnames(expr)
  if (anyDuplicated(rownames(gene))) stop("Gene aggregation did not create unique symbols")
  gene
}

sampled_intensity_matrix <- function(raw_object, max_probes = 50000L) {
  values <- if (inherits(raw_object, "FeatureSet")) {
    oligo::intensity(raw_object)
  } else {
    Biobase::exprs(raw_object)
  }
  if (nrow(values) > max_probes) {
    index <- unique(round(seq(1, nrow(values), length.out = max_probes)))
    values <- values[index, , drop = FALSE]
  }
  values
}

robust_z <- function(x) {
  center <- stats::median(x, na.rm = TRUE)
  spread <- stats::mad(x, center = center, constant = 1, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  0.67448975 * (x - center) / spread
}

save_qc_bundle <- function(raw_object, normalized_expr, metadata, prefix,
                           robust_z_threshold = 5) {
  figure_dir <- "results/figures"
  table_dir <- "results/tables"
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  raw_values <- sampled_intensity_matrix(raw_object)
  normalized_values <- normalized_expr
  stems <- c("pre_rma_boxplot", "post_rma_boxplot", "pre_rma_density_all_arrays",
             "post_rma_density_all_arrays", "rle", "ma_diagnostic", "pca",
             "sample_correlation", "hierarchical_clustering")
  devices <- list(
    png = function(path) png(path, width = 2700, height = 2100, res = 300),
    pdf = function(path) pdf(path, width = 9, height = 7),
    svg = function(path) svg(path, width = 9, height = 7)
  )
  draw <- function(stem, code) {
    for (extension in names(devices)) {
      devices[[extension]](file.path(figure_dir, paste0(prefix, "_", stem, ".", extension)))
      force(code())
      dev.off()
    }
  }
  draw(stems[[1]], function() boxplot(log2(raw_values + 1), outline = FALSE, las = 2,
                                      cex.axis = 0.25, ylab = "log2 raw intensity",
                                      main = paste(prefix, "pre-RMA distributions")))
  draw(stems[[2]], function() boxplot(normalized_values, outline = FALSE, las = 2,
                                      cex.axis = 0.25, ylab = "RMA log2 expression",
                                      main = paste(prefix, "post-RMA distributions")))
  tissue <- factor(metadata$tissue_class)
  tissue_colors <- setNames(grDevices::hcl.colors(length(levels(tissue)), "Dark 3"),
                            levels(tissue))
  density_bundle <- function(values) {
    lapply(seq_len(ncol(values)), function(i) stats::density(values[, i], na.rm = TRUE))
  }
  raw_density <- density_bundle(log2(raw_values + 1))
  normalized_density <- density_bundle(normalized_values)
  draw_density <- function(curves, xlab, title) {
    xlim <- range(vapply(curves, function(curve) range(curve$x), numeric(2)))
    ylim <- c(0, max(vapply(curves, function(curve) max(curve$y), numeric(1))))
    plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "Density", main = title)
    for (i in seq_along(curves)) {
      lines(curves[[i]]$x, curves[[i]]$y,
            col = grDevices::adjustcolor(tissue_colors[[as.character(tissue[[i]])]],
                                         alpha.f = 0.20),
            lwd = 0.55)
    }
    legend("topright", legend = levels(tissue), col = tissue_colors, lwd = 2,
           bty = "n", cex = 0.72, title = "Tissue")
    mtext(sprintf("Every array is shown (n=%d)", length(curves)), side = 3,
          line = 0.2, cex = 0.72)
  }
  draw(stems[[3]], function() draw_density(raw_density, "log2 raw intensity",
                                            paste(prefix, "pre-RMA density")))
  draw(stems[[4]], function() draw_density(normalized_density, "RMA log2 expression",
                                            paste(prefix, "post-RMA density")))
  rle <- sweep(normalized_values, 1, matrixStats::rowMedians(normalized_values), "-")
  draw(stems[[5]], function() boxplot(rle, outline = FALSE, las = 2, cex.axis = 0.25,
                                      ylab = "Relative log expression",
                                      main = paste(prefix, "RLE")))
  ref <- matrixStats::rowMedians(normalized_values)
  draw(stems[[6]], function() {
    sample_index <- seq_len(min(6L, ncol(normalized_values)))
    plot(NULL, xlim = range((normalized_values[, sample_index] + ref) / 2),
         ylim = range(normalized_values[, sample_index] - ref), xlab = "A",
         ylab = "M", main = paste(prefix, "MA diagnostics (first six arrays)"))
    for (i in sample_index) points((normalized_values[, i] + ref) / 2,
                                   normalized_values[, i] - ref, pch = 16,
                                   cex = 0.1, col = adjustcolor(i, alpha.f = 0.25))
    abline(h = 0, lty = 2)
  })
  variable <- order(matrixStats::rowVars(normalized_values), decreasing = TRUE)
  variable <- head(variable, min(5000L, length(variable)))
  pc <- prcomp(t(normalized_values[variable, , drop = FALSE]), scale. = FALSE)
  cls <- factor(metadata$tissue_class)
  draw(stems[[7]], function() {
    plot(pc$x[, 1], pc$x[, 2], col = as.integer(cls), pch = 19,
         xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pc)$importance[2, 1]),
         ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pc)$importance[2, 2]),
         main = paste(prefix, "RMA PCA"))
    legend("topright", levels(cls), col = seq_along(levels(cls)), pch = 19)
  })
  correlation <- cor(normalized_values[variable, , drop = FALSE], method = "pearson")
  draw(stems[[8]], function() heatmap(correlation, Rowv = NA, Colv = NA, labRow = NA,
                                      labCol = NA, main = paste(prefix, "sample correlation")))
  draw(stems[[9]], function() plot(hclust(as.dist(1 - correlation)), labels = FALSE,
                                   main = paste(prefix, "hierarchical clustering"),
                                   xlab = "Samples"))
  medians <- matrixStats::colMedians(normalized_values)
  iqrs <- matrixStats::colIQRs(normalized_values)
  diagnostics <- data.frame(
    accession = prefix,
    geo_accession = colnames(normalized_values),
    sample_median = medians,
    sample_iqr = iqrs,
    median_robust_z = robust_z(medians),
    iqr_robust_z = robust_z(iqrs),
    missing_values = colSums(is.na(normalized_values)),
    non_finite_values = colSums(!is.finite(normalized_values)),
    stringsAsFactors = FALSE
  )
  diagnostics$independent_failure_metrics <-
    as.integer(abs(diagnostics$median_robust_z) > robust_z_threshold) +
    as.integer(abs(diagnostics$iqr_robust_z) > robust_z_threshold) +
    as.integer(diagnostics$missing_values > 0 | diagnostics$non_finite_values > 0)
  diagnostics$technical_outlier_candidate <- diagnostics$independent_failure_metrics >= 2L
  diagnostics$automatic_exclusion <- FALSE
  write.csv(diagnostics, file.path(table_dir, paste0(prefix, "_raw_cel_qc_diagnostics.csv")),
            row.names = FALSE)
  write.csv(data.frame(geo_accession = colnames(normalized_values), PC1 = pc$x[, 1],
                       PC2 = pc$x[, 2], tissue_class = metadata$tissue_class,
                       patient_id = metadata$patient_id),
            file.path(table_dir, paste0(prefix, "_raw_cel_pca_scores.csv")), row.names = FALSE)
  diagnostics
}

run_array_quality_metrics <- function(eset, metadata, prefix) {
  if (!requireNamespace("arrayQualityMetrics", quietly = TRUE)) {
    stop("arrayQualityMetrics is required for the publication workflow")
  }
  outdir <- file.path("results/qc", paste0(prefix, "_array_quality_metrics"))
  dir.create(dirname(outdir), recursive = TRUE, showWarnings = FALSE)
  sample_ids <- extract_gsm(Biobase::sampleNames(eset))
  Biobase::sampleNames(eset) <- sample_ids
  metadata <- metadata[match(sample_ids, toupper(metadata$geo_accession)), ,
                       drop = FALSE]
  if (anyNA(metadata$geo_accession)) stop(prefix, ": AQM metadata alignment failed")
  rownames(metadata) <- Biobase::sampleNames(eset)
  Biobase::pData(eset) <- metadata
  aqm <- arrayQualityMetrics::arrayQualityMetrics(
    eset,
    outdir = outdir,
    force = TRUE,
    do.logtransform = FALSE,
    intgroup = "tissue_class",
    spatial = FALSE,
    reporttitle = paste(prefix, "normalized raw-CEL RMA quality report")
  )
  summary <- data.frame(
    accession = prefix,
    geo_accession = Biobase::sampleNames(eset),
    stringsAsFactors = FALSE
  )
  criterion_columns <- character(0)
  criterion_titles <- character(0)
  for (module in aqm$modules) {
    outliers <- methods::slot(module, "outliers")
    statistic <- methods::slot(outliers, "statistic")
    if (!length(statistic)) next
    module_id <- methods::slot(module, "id")
    module_title <- methods::slot(module, "title")
    column <- paste0("flag_", gsub("[^a-zA-Z0-9]+", "_", tolower(module_id)))
    summary[[column]] <- seq_len(nrow(summary)) %in% methods::slot(outliers, "which")
    criterion_columns <- c(criterion_columns, column)
    criterion_titles <- c(criterion_titles, paste0(column, "=", module_title))
  }
  if (length(criterion_columns)) {
    summary$aqm_flag_count <- rowSums(summary[, criterion_columns, drop = FALSE])
    summary$aqm_flagged_any <- summary$aqm_flag_count > 0
    summary$flagged_criteria <- apply(summary[, criterion_columns, drop = FALSE], 1,
                                      function(flags) paste(criterion_columns[flags],
                                                             collapse = ";"))
  } else {
    summary$aqm_flag_count <- 0L
    summary$aqm_flagged_any <- FALSE
    summary$flagged_criteria <- ""
  }
  summary$criterion_definitions <- paste(criterion_titles, collapse = ";")
  summary$exclusion_changed <- FALSE
  summary$exclusion_rule <- paste(
    "AQM flags are review signals; exclusion still requires at least two",
    "independent severe technical failures"
  )
  write.csv(summary,
            file.path("results/tables", paste0(prefix,
                                               "_array_quality_metrics_summary.csv")),
            row.names = FALSE)
  summary
}

run_affyplm_nuse <- function(raw_object, metadata, prefix) {
  if (!requireNamespace("affyPLM", quietly = TRUE)) {
    stop("affyPLM is required for GPL96 NUSE")
  }
  plm <- affyPLM::fitPLM(
    raw_object,
    output.param = list(weights = FALSE, residuals = FALSE,
                        varcov = "none", resid.SE = TRUE),
    verbosity.level = 1
  )
  dir.create("data/interim", recursive = TRUE, showWarnings = FALSE)
  saveRDS(plm, file.path("data/interim", paste0(prefix, "_affyplm_plmset.rds")),
          compress = FALSE)
  values <- affyPLM::NUSE(plm, type = "values")
  raw_stats <- affyPLM::NUSE(plm, type = "stats")
  if (ncol(raw_stats) == ncol(values)) {
    stats <- as.data.frame(t(raw_stats))
    colnames(stats) <- rownames(raw_stats)
  } else {
    stats <- as.data.frame(raw_stats)
  }
  if (nrow(stats) == ncol(values)) {
    stats$geo_accession <- extract_gsm(colnames(values))
  } else if (nrow(stats) == length(Biobase::sampleNames(raw_object))) {
    stats$geo_accession <- extract_gsm(Biobase::sampleNames(raw_object))
  } else {
    stop(prefix, ": unexpected NUSE statistics dimensions")
  }
  stats$accession <- prefix
  stats$method <- "affyPLM probe-level model; NUSE is valid for GPL96 AffyBatch only"
  write.csv(stats, file.path("results/tables", paste0(prefix, "_nuse_summary.csv")),
            row.names = FALSE)
  draw_nuse <- function() {
    boxplot(values, outline = FALSE, las = 2, cex.axis = 0.25,
            ylab = "Normalized unscaled standard error",
            main = paste(prefix, "affyPLM NUSE"))
    abline(h = 1, lty = 2)
  }
  for (extension in c("png", "pdf", "svg")) {
    path <- file.path("results/figures", paste0(prefix, "_nuse.", extension))
    if (extension == "png") png(path, width = 2700, height = 2100, res = 300)
    if (extension == "pdf") pdf(path, width = 9, height = 7)
    if (extension == "svg") svg(path, width = 9, height = 7)
    draw_nuse()
    dev.off()
  }
  invisible(stats)
}

update_sample_exclusion_log <- function(metadata, diagnostics, accession) {
  path <- "data/metadata/sample_exclusion_log.csv"
  candidates <- diagnostics$technical_outlier_candidate
  current <- data.frame(
    accession = accession,
    geo_accession = diagnostics$geo_accession,
    qc_candidate = candidates,
    exclusion_applied = FALSE,
    decision = ifelse(candidates, "borderline_include_primary_and_review_sensitivity",
                      "include_prespecified_technically_valid_cohort"),
    evidence = paste0("independent_failure_metrics=", diagnostics$independent_failure_metrics),
    rule = "No PCA-only exclusion; exclusion requires severe failure supported by at least two independent QC metrics",
    stringsAsFactors = FALSE
  )
  if (file.exists(path)) {
    previous <- read.csv(path, stringsAsFactors = FALSE)
    previous <- previous[previous$accession != accession, , drop = FALSE]
    current <- rbind(previous, current)
  }
  write.csv(current[order(current$accession, current$geo_accession), ], path, row.names = FALSE)
}
