# scripts/12_annotate_limma_results.R
# Add gene annotations while preserving one row per probe set.

# 1. Install missing annotation packages
required_packages <- c("AnnotationDbi", "hgu133a.db")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  BiocManager::install(
    missing_packages,
    ask = FALSE,
    update = FALSE
  )
}

library(AnnotationDbi)
library(hgu133a.db)

dir.create("results", showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)

# 2. Read the saved limma results
de_results <- read.csv(
  "results/RA_vs_OA_limma_all_probes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  all(c("probe_id", "logFC", "adj.P.Val") %in% names(de_results)),
  !anyNA(de_results$probe_id),
  !anyDuplicated(de_results$probe_id)
)

# 3. Retrieve all mappings for recognized probe IDs
annotation_db <- hgu133a.db::hgu133a.db

valid_probe_ids <- intersect(
  de_results$probe_id,
  AnnotationDbi::keys(annotation_db, keytype = "PROBEID")
)

stopifnot(length(valid_probe_ids) > 0L)

annotation_long <- AnnotationDbi::select(
  annotation_db,
  keys = valid_probe_ids,
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "PROBEID"
)

# Save the full mapping table, including multiple mappings.
write.csv(
  annotation_long,
  "results/HG_U133A_probe_gene_mapping.csv",
  row.names = FALSE
)

# 4. Add annotations without duplicating result rows
mapping_rows <- split(
  seq_len(nrow(annotation_long)),
  annotation_long$PROBEID
)

clean_values <- function(values) {
  values <- as.character(values)
  unique(values[!is.na(values) & nzchar(values)])
}

collapse_mapping <- function(probe, column) {
  rows <- mapping_rows[[probe]]
  
  if (is.null(rows)) {
    return(NA_character_)
  }
  
  values <- clean_values(annotation_long[[column]][rows])
  
  if (length(values) == 0L) {
    return(NA_character_)
  }
  
  paste(values, collapse = "; ")
}

annotated_results <- de_results

annotated_results$gene_symbol <- vapply(
  de_results$probe_id,
  collapse_mapping,
  character(1),
  column = "SYMBOL"
)

annotated_results$gene_name <- vapply(
  de_results$probe_id,
  collapse_mapping,
  character(1),
  column = "GENENAME"
)

annotated_results$entrez_id <- vapply(
  de_results$probe_id,
  collapse_mapping,
  character(1),
  column = "ENTREZID"
)

annotated_results$n_entrez_mappings <- vapply(
  de_results$probe_id,
  function(probe) {
    rows <- mapping_rows[[probe]]
    
    if (is.null(rows)) {
      return(0L)
    }
    
    length(clean_values(annotation_long$ENTREZID[rows]))
  },
  integer(1)
)

# Put annotations before the statistical results.
annotation_columns <- c(
  "probe_id", "gene_symbol", "gene_name",
  "entrez_id", "n_entrez_mappings"
)

annotated_results <- annotated_results[
  , c(
    annotation_columns,
    setdiff(names(de_results), "probe_id")
  )
]

stopifnot(
  nrow(annotated_results) == nrow(de_results),
  identical(annotated_results$probe_id, de_results$probe_id)
)

# 5. Save annotated tables
significant <- (
  !is.na(annotated_results$adj.P.Val) &
    annotated_results$adj.P.Val < 0.05
)

large_change <- abs(annotated_results$logFC) >= 1

write.csv(
  annotated_results,
  "results/RA_vs_OA_limma_annotated_all_probes.csv",
  row.names = FALSE
)

write.csv(
  annotated_results[significant, , drop = FALSE],
  "results/RA_vs_OA_limma_annotated_FDR05.csv",
  row.names = FALSE
)

write.csv(
  annotated_results[
    significant & large_change, , drop = FALSE
  ],
  "results/RA_vs_OA_limma_annotated_FDR05_logFC1.csv",
  row.names = FALSE
)

# Record annotation versions and session information.
writeLines(
  c(
    paste("Annotation date:", Sys.Date()),
    paste(
      "hgu133a.db version:",
      packageVersion("hgu133a.db")
    ),
    "",
    capture.output(sessionInfo())
  ),
  "docs/session_info_annotation.txt"
)

# 6. Print summary
cat("\n--- Annotation summary ---\n")
cat("Total probe sets:", nrow(annotated_results), "\n")
cat(
  "Probe sets with gene symbols:",
  sum(!is.na(annotated_results$gene_symbol)), "\n"
)
cat(
  "Probe sets mapping to multiple Entrez genes:",
  sum(annotated_results$n_entrez_mappings > 1L), "\n"
)

cat("\n--- Differential expression summary ---\n")
cat("FDR < 0.05:", sum(significant), "\n")
cat(
  "Higher in RA:",
  sum(significant & annotated_results$logFC > 0), "\n"
)
cat(
  "Lower in RA:",
  sum(significant & annotated_results$logFC < 0), "\n"
)
cat(
  "FDR < 0.05 and |log2FC| >= 1:",
  sum(significant & large_change), "\n"
)

cat("\nTop 10 annotated probe sets:\n")
print(
  head(
    annotated_results[
      , c(
        "probe_id", "gene_symbol",
        "n_entrez_mappings", "logFC", "adj.P.Val"
      )
    ],
    10
  ),
  row.names = FALSE
)

cat(
  "\nAnnotated results saved:",
  file.exists(
    "results/RA_vs_OA_limma_annotated_all_probes.csv"
  ),
  "\n"
)