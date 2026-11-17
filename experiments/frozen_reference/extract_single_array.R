# Single-array sufficient statistics; no outcomes or cross-sample fitting.
# Rscript --vanilla experiments/frozen_reference/extract_single_array.R CEL_ROOT OUT
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages({library(oligo); library(affy); library(preprocessCore)})
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
for (accession in c("GSE44076", "GSE41258")) {
  meta <- read.csv(file.path("data/metadata", paste0(tolower(accession), "_samples.csv")))
  meta <- meta[meta$inclusion_status == "included", ]
  files <- list.files(file.path(args[1], accession), recursive=TRUE, full.names=TRUE,
                      pattern="[Cc][Ee][Ll]\\.gz$")
  ids <- regmatches(basename(files), regexpr("GSM[0-9]+", basename(files)))
  stopifnot(!anyDuplicated(ids), all(meta$geo_accession %in% ids))
  dest <- file.path(args[2], accession)
  dir.create(dest, showWarnings=FALSE)
  probe_indices <- NULL
  for (gsm in meta$geo_accession) {
    output <- file.path(dest, paste0(gsm, ".rds"))
    # Resume only verified per-CEL and code-identical statistics.
    cel <- files[match(gsm, ids)]
    hash <- digest::digest(file=cel, algo="sha256")
    code_hash <- digest::digest(file="experiments/frozen_reference/extract_single_array.R", algo="sha256")
    if (file.exists(output)) {
      old <- readRDS(output)
      if (identical(old$cel_sha256, hash) && identical(old$code_sha256, code_hash)) next
      stop("Stale cache must be inspected: ", output)
    }
    if (is.null(probe_indices) && accession == "GSE44076") {
      raw <- oligo::read.celfiles(cel, verbose=FALSE)
      stopifnot(Biobase::annotation(raw) == "pd.hg.u219")
      info <- oligo::getProbeInfo(raw, field=c("fid", "fsetid"), sortBy="fid")
      values <- as.numeric(Biobase::exprs(raw)[info$fid, 1])
      stopifnot(isTRUE(all.equal(values, as.numeric(oligo::pm(raw)))))
      labels <- info$man_fsetid
      probe_indices <- info$fid
      rm(raw)
    } else if (is.null(probe_indices)) {
      raw <- affy::ReadAffy(filenames=cel)
      stopifnot(affy::cdfName(raw) == "HG-U133A")
      indices <- affy::indexProbes(raw, "pm")
      values <- as.numeric(Biobase::exprs(raw)[unlist(indices, use.names=FALSE), 1])
      labels <- rep(names(indices), lengths(indices))
      probe_indices <- unlist(indices, use.names=FALSE)
      rm(raw)
    }
    raw_cel <- affyio::read.celfile(cel, intensity.means.only=TRUE)
    values <- as.numeric(raw_cel$INTENSITY$MEAN[probe_indices])
    stopifnot(length(values) == length(labels), all(is.finite(values)))
    corrected <- as.numeric(preprocessCore::rma.background.correct(matrix(values, ncol=1)))
    stopifnot(all(is.finite(corrected)), all(corrected > 0))
    ranks <- (rank(corrected, ties.method="average") - 0.5) / length(corrected)
    summary <- tapply(ranks, labels, median)
    quantiles <- quantile(corrected, probs=seq(0, 1, length.out=4097), names=FALSE, type=7)
    result <- list(sample=gsm, cel_sha256=hash, code_sha256=code_hash,
                   pm_count=length(values), percentile=summary, quantiles=quantiles)
    saveRDS(result, output)
    write.csv(data.frame(probe_id=names(summary), percentile=as.numeric(summary)),
              file.path(dest, paste0(gsm, "_ranks.csv")), row.names=FALSE)
    writeBin(as.double(quantiles), file.path(dest, paste0(gsm, "_quantiles.bin")),
             size=8, endian="little")
    write.csv(data.frame(sample_id=gsm, cel_sha256=hash, pm_count=length(values),
                         code_sha256=code_hash), file.path(dest, paste0(gsm, "_input.csv")),
              row.names=FALSE)
    cat(accession, gsm, length(values), length(summary), "PASS\n")
    rm(raw_cel)
  }
}
writeLines(capture.output(sessionInfo()), file.path(args[2], "R_session.txt"))
