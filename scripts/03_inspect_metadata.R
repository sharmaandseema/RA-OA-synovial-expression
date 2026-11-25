# Inspect processing notes and sample characteristics.

eset <- readRDS(
  "data_processed/GSE55235_GPL96_expression_set.rds"
)

metadata <- Biobase::pData(eset)

fields <- grep(
  "characteristics|processing|protocol|source_name",
  names(metadata),
  value = TRUE,
  ignore.case = TRUE
)

report <- capture.output({
  cat("Expression dimensions: probe sets x samples\n")
  print(dim(Biobase::exprs(eset)))
  
  cat("\nMetadata fields available:\n")
  print(names(metadata))
  
  cat("\nProcessing notes and sample characteristics:\n")
  
  for (field in fields) {
    cat("\n---", field, "---\n")
    values <- unique(as.character(metadata[[field]]))
    values <- values[!is.na(values) & nzchar(values)]
    cat(values, sep = "\n")
    cat("\n")
  }
})

output_path <- "docs/GSE55235_metadata_inspection.txt"
writeLines(report, output_path)

cat(report, sep = "\n")
cat("\nMetadata inspection saved:", file.exists(output_path), "\n")