# Independent raw-PM check before summarization: single vs multicolumn correction.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages({library(oligo); library(affy); library(preprocessCore)})
rows <- list()
for (acc in c("GSE44076", "GSE41258")) {
  meta <- read.csv(file.path("data/metadata", paste0(tolower(acc), "_samples.csv")))
  meta <- meta[meta$inclusion_status == "included", ]
  # Deterministic first and last sample of every tissue class, not performance selected.
  selected <- unique(unlist(lapply(split(meta$geo_accession, meta$tissue_class),
                                  function(x) c(head(x, 1), tail(x, 1)))))
  files <- list.files(file.path(args[1], acc), full.names=TRUE, recursive=TRUE, pattern="[Cc][Ee][Ll]\\.gz$")
  ids <- regmatches(basename(files), regexpr("GSM[0-9]+", basename(files)))
  if (acc == "GSE44076") {
    raw <- oligo::read.celfiles(files[match(selected[1], ids)], verbose=FALSE)
    info <- oligo::getProbeInfo(raw, field=c("fid", "fsetid"), sortBy="fid")
    indices <- info$fid
  } else {
    raw <- affy::ReadAffy(filenames=files[match(selected[1], ids)])
    indices <- unlist(affy::indexProbes(raw, "pm"), use.names=FALSE)
  }
  values <- vapply(selected, function(gsm) {
    cel <- affyio::read.celfile(files[match(gsm, ids)], intensity.means.only=TRUE)
    cel$INTENSITY$MEAN[indices]
  }, numeric(length(indices)))
  batch <- preprocessCore::rma.background.correct(values, copy=TRUE)
  for (i in seq_along(selected)) {
    alone <- as.numeric(preprocessCore::rma.background.correct(matrix(values[,i], ncol=1)))
    difference <- max(abs(alone-batch[,i]))
    relative <- max(abs(alone-batch[,i])/pmax(abs(alone), 1e-12))
    rank_diff <- max(abs(rank(alone)-rank(batch[,i])))
    stopifnot(difference <= 1e-12, relative <= 1e-12, rank_diff == 0)
    rows[[length(rows)+1L]] <- data.frame(accession=acc, sample_id=selected[i],
        pm_count=length(indices), corrected_pm_max_abs=difference,
        corrected_pm_max_relative=relative, pm_rank_max_abs=rank_diff, status="PASS")
  }
}
write.csv(do.call(rbind, rows), args[2], row.names=FALSE)
