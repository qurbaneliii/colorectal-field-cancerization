source("R/utils.R")

expected_counts <- c(GSE44076 = 246L, GSE41258 = 390L)
audit_rows <- list()
for (accession in names(expected_counts)) {
  tar_files <- discover_files(accession, "_RAW\\.tar$")
  if (length(tar_files) != 1L) stop(accession, ": expected one raw TAR, observed ", length(tar_files))
  destination <- file.path("data/interim/cel", accession)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  members <- utils::untar(tar_files[[1]], list = TRUE)
  cel_members <- members[grepl("\\.CEL\\.gz$", members, ignore.case = TRUE)]
  if (length(cel_members) != expected_counts[[accession]]) {
    stop(accession, ": expected ", expected_counts[[accession]], " CEL members; observed ",
         length(cel_members))
  }
  if (anyDuplicated(tolower(basename(cel_members)))) stop(accession, ": duplicate member names")
  existing <- list.files(destination, pattern = "\\.CEL\\.gz$", recursive = TRUE,
                         full.names = TRUE, ignore.case = TRUE)
  if (length(existing) != length(cel_members)) utils::untar(tar_files[[1]], exdir = destination)
  extracted <- list.files(destination, pattern = "\\.CEL\\.gz$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)
  if (length(extracted) != length(cel_members)) stop(accession, ": extraction count mismatch")
  readable <- vapply(extracted, function(path) {
    connection <- gzfile(path, "rb")
    on.exit(close(connection), add = TRUE)
    length(readBin(connection, what = "raw", n = 64L)) == 64L
  }, logical(1))
  if (!all(readable)) stop(accession, ": unreadable compressed CEL files: ",
                           paste(basename(extracted[!readable]), collapse = ", "))
  metadata <- read.csv(file.path("data/metadata", paste0(tolower(accession), "_samples.csv")),
                       stringsAsFactors = FALSE)
  aligned <- validate_cel_metadata_alignment(extracted, metadata)
  member_table <- data.frame(
    accession = accession,
    archive = basename(tar_files[[1]]),
    member = cel_members,
    gsm = extract_gsm(cel_members),
    extracted_path = file.path("data/interim/cel", accession, cel_members),
    compressed_cel_readable = readable[match(tolower(basename(cel_members)),
                                              tolower(basename(extracted)))],
    metadata_match = extract_gsm(cel_members) %in% toupper(metadata$geo_accession),
    stringsAsFactors = FALSE
  )
  write.csv(member_table, file.path("data/metadata", paste0(accession, "_tar_members.csv")),
            row.names = FALSE)
  audit_rows[[accession]] <- data.frame(
    accession = accession, archive = basename(tar_files[[1]]),
    archive_members = length(members), cel_members = length(cel_members),
    unique_cel_names = length(unique(tolower(basename(cel_members)))),
    extracted_cels = length(extracted), readable_cels = sum(readable),
    metadata_matches = nrow(aligned), status = "PASS", stringsAsFactors = FALSE
  )
}
audit <- do.call(rbind, audit_rows)
write.csv(audit, "data/metadata/cel_archive_audit.csv", row.names = FALSE)
report <- c("# CEL archive audit", "", "Raw archives were read and extracted without modifying `data/raw/`.", "",
            paste(capture.output(print(audit, row.names = FALSE)), collapse = "\n"), "",
            "All CEL filenames were unique, readable through gzip, and matched one GEO sample identifier.")
writeLines(report, "reports/cel_archive_audit.md")
