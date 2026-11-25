# Create and validate the sample manifest.

eset <- readRDS(
  "data_processed/GSE55235_GPL96_expression_set.rds"
)

metadata <- Biobase::pData(eset)

disease <- trimws(as.character(metadata[["disease state:ch1"]]))

group_mapping <- c(
  "healthy control" = "Healthy",
  "osteoarthritis" = "OA",
  "synovial tissue isolated from osteoarthritic joint" = "OA",
  "rheumatoid arthritis" = "RA"
)

condition <- unname(group_mapping[disease])

# Stop if any disease description is unrecognized.
stopifnot(!anyNA(condition))

manifest <- data.frame(
  sample_id = rownames(metadata),
  title = as.character(metadata$title),
  source_disease_description = disease,
  condition = condition,
  primary_comparison = condition %in% c("OA", "RA"),
  stringsAsFactors = FALSE
)

# Independently check the group against the sample title.
title_matches <- (
  manifest$condition == "Healthy" &
    grepl("healthy joint", manifest$title, ignore.case = TRUE)
) | (
  manifest$condition == "OA" &
    grepl("osteoarthritic joint", manifest$title, ignore.case = TRUE)
) | (
  manifest$condition == "RA" &
    grepl("rheumatoid arthritis joint", manifest$title, ignore.case = TRUE)
)

stopifnot(
  all(title_matches),
  !anyDuplicated(manifest$sample_id),
  identical(manifest$sample_id, colnames(Biobase::exprs(eset))),
  sum(condition == "Healthy") == 10L,
  sum(condition == "OA") == 10L,
  sum(condition == "RA") == 10L
)

output_path <- "data_processed/GSE55235_sample_manifest.csv"

write.csv(manifest, output_path, row.names = FALSE)

print(
  manifest[, c("sample_id", "condition", "primary_comparison")],
  row.names = FALSE
)

cat("\nSamples per group:\n")
print(table(manifest$condition))

cat(
  "\nSamples in the primary RA versus OA comparison:",
  sum(manifest$primary_comparison),
  "\n"
)

cat("Sample manifest saved:", file.exists(output_path), "\n")