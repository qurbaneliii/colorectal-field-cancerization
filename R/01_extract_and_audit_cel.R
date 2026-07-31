source("R/utils.R")

for (accession in c("GSE44076", "GSE41258")) {
  tar_files <- discover_files(accession, paste0(accession, ".*_RAW\\.tar$"))
  if (length(tar_files) != 1L) {
    stop(accession, ": expected one raw TAR, observed ", length(tar_files))
  }
  destination <- file.path("data/interim/cel", accession)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  members <- utils::untar(tar_files[[1]], list = TRUE)
  if (!all(grepl("\\.CEL\\.gz$", members, ignore.case = TRUE))) {
    warning(accession, ": non-CEL members were found in raw TAR")
  }
  existing <- list.files(destination, recursive = TRUE)
  if (length(existing) != length(members)) {
    utils::untar(tar_files[[1]], exdir = destination)
  }
  write.csv(
    data.frame(accession = accession, archive = tar_files[[1]], member = members,
               stringsAsFactors = FALSE),
    file.path("data/metadata", paste0(accession, "_tar_members.csv")),
    row.names = FALSE
  )
}
