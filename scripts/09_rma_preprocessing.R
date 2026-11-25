# Process all 30 CEL files together using RMA.

stopifnot(
  file.exists("ra-oa-synovial-expression.Rproj"),
  requireNamespace("affy", quietly = TRUE),
  requireNamespace("hgu133acdf", quietly = TRUE)
)

file_manifest <- read.csv(
  "data_processed/GSE55235_cel_file_manifest.csv",
  stringsAsFactors = FALSE
)

manifest <- read.csv(
  "data_processed/GSE55235_sample_manifest.csv",
  stringsAsFactors = FALSE
)

stopifnot(
  identical(file_manifest$sample_id, manifest$sample_id),
  nrow(file_manifest) == 30L,
  all(file.exists(file_manifest$cel_path))
)

# Read CEL files in the verified manifest order.
raw_arrays <- affy::ReadAffy(
  filenames = file_manifest$cel_path
)

cat("Array type:", Biobase::annotation(raw_arrays), "\n")
cat("Arrays read:", ncol(raw_arrays), "\n")

# Check that the imported names retain the expected sample order.
imported_ids <- sub(
  ".*(GSM[0-9]+).*",
  "\\1",
  Biobase::sampleNames(raw_arrays)
)

stopifnot(
  identical(imported_ids, manifest$sample_id),
  tolower(gsub(
    "[^A-Za-z0-9]", "",
    Biobase::annotation(raw_arrays)
  )) == "hgu133a"
)

Biobase::sampleNames(raw_arrays) <- manifest$sample_id

# RMA produces log2-scale expression values.
rma_eset <- affy::rma(raw_arrays)
rma_expression <- Biobase::exprs(rma_eset)

stopifnot(
  identical(colnames(rma_expression), manifest$sample_id),
  nrow(rma_expression) == 22283L,
  ncol(rma_expression) == 30L,
  all(is.finite(rma_expression))
)

# Attach the checked sample information.
sample_data <- manifest
rownames(sample_data) <- sample_data$sample_id

Biobase::pData(rma_eset) <- sample_data

output_path <- "data_processed/GSE55235_RMA_expression_set.rds"
saveRDS(rma_eset, output_path)

writeLines(
  c(
    "Input: 30 original GSE55235 CEL files.",
    "Method: affy::rma with default settings.",
    "All 30 arrays were processed together.",
    "Steps: background correction, quantile normalization,",
    "and median-polish probe-set summarization.",
    "Output expression values are on the log2 scale.",
    "Primary planned comparison: RA versus OA."
  ),
  "docs/RMA_preprocessing_notes.txt"
)

cat("\nRMA probe sets:", nrow(rma_expression), "\n")
cat("RMA samples:", ncol(rma_expression), "\n")

cat("\nRMA expression quantiles:\n")
print(quantile(
  rma_expression,
  probs = c(0, 0.25, 0.5, 0.75, 1)
))

cat("\nRMA expression set saved:",
    file.exists(output_path), "\n")