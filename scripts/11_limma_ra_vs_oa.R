# scripts/11_limma_ra_vs_oa.R
# Primary comparison: RA versus OA
# Input: RMA-normalized log2 expression
# Output: probe-level differential expression results

library(Biobase)
library(limma)

dir.create("results", showWarnings = FALSE)
dir.create("data_processed", showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)

# 1. Read the saved RMA expression set
eset <- readRDS(
  "data_processed/GSE55235_RMA_expression_set.rds"
)

x <- exprs(eset)
metadata <- pData(eset)

stopifnot(
  identical(colnames(x), rownames(metadata)),
  "condition" %in% colnames(metadata),
  all(is.finite(x))
)

# 2. Select OA and RA samples
selected <- as.character(metadata$condition) %in% c("OA", "RA")

x_ra_oa <- x[, selected, drop = FALSE]
metadata_ra_oa <- metadata[selected, , drop = FALSE]

metadata_ra_oa$condition <- factor(
  as.character(metadata_ra_oa$condition),
  levels = c("OA", "RA")
)

stopifnot(
  identical(colnames(x_ra_oa), rownames(metadata_ra_oa)),
  sum(metadata_ra_oa$condition == "OA") == 10L,
  sum(metadata_ra_oa$condition == "RA") == 10L
)

cat("Samples included:\n")
print(table(metadata_ra_oa$condition))

cat("\nProbe sets included:", nrow(x_ra_oa), "\n")

# 3. Construct the design matrix
# Separate mean expression estimates for OA and RA
design <- model.matrix(
  ~ 0 + condition,
  data = metadata_ra_oa
)

colnames(design) <- levels(metadata_ra_oa$condition)
rownames(design) <- colnames(x_ra_oa)

stopifnot(
  identical(rownames(design), colnames(x_ra_oa)),
  qr(design)$rank == ncol(design)
)

# 4. Define RA minus OA
contrast_matrix <- makeContrasts(
  RA_vs_OA = RA - OA,
  levels = design
)

# RMA expression is already log2-normalized.
# Do not log-transform again.
fit <- lmFit(x_ra_oa, design)
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ra_oa <- eBayes(fit_contrast)

# 5. Extract all probe-level results
de_results <- topTable(
  fit_ra_oa,
  coef = "RA_vs_OA",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

de_results <- data.frame(
  probe_id = rownames(de_results),
  de_results,
  row.names = NULL,
  check.names = FALSE
)

# FDR-significant results
significant <- (
  !is.na(de_results$adj.P.Val) &
    de_results$adj.P.Val < 0.05
)

de_fdr05 <- de_results[significant, , drop = FALSE]

# Additional descriptive subset:
# FDR < 0.05 and absolute log2 fold change >= 1
# This corresponds to an estimated >= 2-fold difference.
de_fdr05_fc2 <- de_results[
  significant & abs(de_results$logFC) >= 1,
  ,
  drop = FALSE
]

# 6. Save results and analysis inputs
write.csv(
  de_results,
  "results/RA_vs_OA_limma_all_probes.csv",
  row.names = FALSE
)

write.csv(
  de_fdr05,
  "results/RA_vs_OA_limma_FDR05.csv",
  row.names = FALSE
)

write.csv(
  de_fdr05_fc2,
  "results/RA_vs_OA_limma_FDR05_logFC1.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    sample_id = rownames(metadata_ra_oa),
    metadata_ra_oa,
    row.names = NULL,
    check.names = FALSE
  ),
  "results/RA_vs_OA_sample_metadata.csv",
  row.names = FALSE
)

write.csv(
  design,
  "results/RA_vs_OA_design_matrix.csv",
  row.names = TRUE
)

saveRDS(
  list(
    fit = fit_ra_oa,
    design = design,
    contrast = contrast_matrix,
    metadata = metadata_ra_oa,
    results = de_results
  ),
  "data_processed/RA_vs_OA_limma_analysis.rds"
)

writeLines(
  capture.output(sessionInfo()),
  "docs/session_info_limma.txt"
)

# 7. Print a concise summary
cat("\n--- RA versus OA results ---\n")
cat("Probe sets tested:", nrow(de_results), "\n")
cat("FDR < 0.05:", nrow(de_fdr05), "\n")
cat(
  "Higher in RA (FDR < 0.05):",
  sum(de_fdr05$logFC > 0), "\n"
)
cat(
  "Lower in RA (FDR < 0.05):",
  sum(de_fdr05$logFC < 0), "\n"
)
cat(
  "FDR < 0.05 and |log2FC| >= 1:",
  nrow(de_fdr05_fc2), "\n"
)

cat("\nTop 10 probe sets:\n")
print(
  head(
    de_results[
      , c("probe_id", "logFC", "AveExpr", "P.Value", "adj.P.Val")
    ],
    10
  ),
  row.names = FALSE
)

cat(
  "\nFull results saved:",
  file.exists("results/RA_vs_OA_limma_all_probes.csv"),
  "\n"
)