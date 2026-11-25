# Download the original CEL archive and inspect its contents.

stopifnot(
  file.exists("ra-oa-synovial-expression.Rproj"),
  dir.exists("data_raw"),
  dir.exists("docs")
)

archive_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE55nnn/GSE55235/suppl/GSE55235_RAW.tar"
)

archive_path <- "data_raw/GSE55235_RAW.tar"

# Allow time for the download; restore the previous setting afterward.
download_archive <- function() {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = max(600, old_timeout))
  
  download.file(
    archive_url,
    destfile = archive_path,
    mode = "wb",
    method = "libcurl"
  )
}

if (!file.exists(archive_path)) {
  download_archive()
}

stopifnot(
  file.exists(archive_path),
  file.info(archive_path)$size > 0
)

# List contents without extracting them.
archive_contents <- untar(archive_path, list = TRUE)

cel_files <- archive_contents[
  grepl("\\.cel(\\.gz)?$", archive_contents, ignore.case = TRUE)
]

writeLines(
  archive_contents,
  "docs/GSE55235_raw_archive_contents.txt"
)

checksum <- unname(tools::md5sum(archive_path))

write.csv(
  data.frame(
    file = archive_path,
    source_url = archive_url,
    size_bytes = file.info(archive_path)$size,
    md5_downloaded_file = checksum
  ),
  "docs/GSE55235_raw_archive_record.csv",
  row.names = FALSE
)

cat("Archive size (MB):",
    round(file.info(archive_path)$size / 1e6, 2), "\n")

cat("CEL files in archive:", length(cel_files), "\n")
cat("\nFirst six CEL filenames:\n")
print(head(cel_files))

cat("\nArchive record saved:",
    file.exists("docs/GSE55235_raw_archive_record.csv"), "\n")

stopifnot(length(cel_files) == 30L)