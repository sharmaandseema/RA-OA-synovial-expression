# Extract CEL files and match them to the sample manifest.

archive_path <- "data_raw/GSE55235_RAW.tar"
cel_directory <- "data_raw/GSE55235_CEL"

stopifnot(file.exists(archive_path))
dir.create(cel_directory, recursive = TRUE, showWarnings = FALSE)

untar(archive_path, exdir = cel_directory)

cel_paths <- list.files(
  cel_directory,
  pattern = "\\.cel(\\.gz)?$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

manifest <- read.csv(
  "data_processed/GSE55235_sample_manifest.csv",
  stringsAsFactors = FALSE
)

# Extract the GEO sample accession from each filename.
cel_ids <- sub(
  "^(GSM[0-9]+)_.*$",
  "\\1",
  basename(cel_paths)
)

stopifnot(
  length(cel_paths) == 30L,
  !anyDuplicated(cel_ids),
  setequal(cel_ids, manifest$sample_id)
)

# Put files into the same order as the sample manifest.
matched_paths <- cel_paths[match(manifest$sample_id, cel_ids)]

file_manifest <- data.frame(
  sample_id = manifest$sample_id,
  condition = manifest$condition,
  cel_path = matched_paths,
  size_bytes = file.info(matched_paths)$size,
  stringsAsFactors = FALSE
)

stopifnot(
  all(file.exists(file_manifest$cel_path)),
  all(file_manifest$size_bytes > 0)
)

output_path <- "data_processed/GSE55235_cel_file_manifest.csv"

write.csv(file_manifest, output_path, row.names = FALSE)

print(
  data.frame(
    sample_id = file_manifest$sample_id,
    condition = file_manifest$condition,
    cel_file = basename(file_manifest$cel_path)
  ),
  row.names = FALSE
)

cat("\nCEL files matched:", nrow(file_manifest), "\n")
cat("File manifest saved:", file.exists(output_path), "\n")