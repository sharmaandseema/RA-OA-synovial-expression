# Download and inspect the processed GEO expression data.

stopifnot(
  file.exists("ra-oa-synovial-expression.Rproj"),
  dir.exists("data_raw"),
  dir.exists("data_processed")
)

geo_data <- GEOquery::getGEO(
  "GSE55235",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = "data_raw"
)

# Check platforms before selecting the expression set.
platforms <- vapply(
  geo_data,
  Biobase::annotation,
  character(1)
)

print(platforms)

selected <- which(platforms == "GPL96")
stopifnot(length(selected) == 1L)

eset <- geo_data[[selected]]
expression_matrix <- Biobase::exprs(eset)
sample_metadata <- Biobase::pData(eset)

stopifnot(
  identical(
    colnames(expression_matrix),
    rownames(sample_metadata)
  )
)

# Save the imported data without transforming it.
expression_path <- "data_processed/GSE55235_GPL96_expression_set.rds"
metadata_path <- "data_processed/GSE55235_sample_metadata_original.csv"

saveRDS(eset, expression_path)

write.csv(
  sample_metadata,
  metadata_path,
  row.names = TRUE
)

cat("\nProbe sets:", nrow(expression_matrix), "\n")
cat("Samples:", ncol(expression_matrix), "\n")

cat("\nSample titles:\n")
print(
  data.frame(
    sample_id = rownames(sample_metadata),
    title = sample_metadata$title
  ),
  row.names = FALSE
)

cat("\nExpression-value quantiles:\n")
print(
  quantile(
    expression_matrix,
    probs = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE
  )
)

cat("\nExpression set saved:", file.exists(expression_path), "\n")
cat("Sample metadata saved:", file.exists(metadata_path), "\n")