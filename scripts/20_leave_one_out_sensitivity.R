# scripts/20_leave_one_out_sensitivity.R
# Internal sensitivity analysis for eight selected genes.
# RMA normalization is kept fixed; limma is refitted each time.
# This is not independent validation or a reason to remove samples.

library(Biobase)
library(limma)

dir.create("results/sensitivity", recursive = TRUE, showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)

# 1. Read saved inputs
eset <- readRDS(
  "data_processed/GSE55235_RMA_expression_set.rds"
)

chosen <- read.csv(
  "results/RA_vs_OA_selected_gene_probe_statistics.csv",
  stringsAsFactors = FALSE
)

x <- exprs(eset)
metadata <- pData(eset)

stopifnot(
  identical(colnames(x), rownames(metadata)),
  "condition" %in% names(metadata)
)

selected <- as.character(metadata$condition) %in% c("OA", "RA")
metadata <- metadata[selected, , drop = FALSE]
x <- x[, selected, drop = FALSE]

metadata$condition <- factor(
  as.character(metadata$condition),
  levels = c("OA", "RA")
)

stopifnot(
  ncol(x) == 20L,
  sum(metadata$condition == "OA") == 10L,
  sum(metadata$condition == "RA") == 10L,
  identical(colnames(x), rownames(metadata)),
  all(chosen$probe_id %in% rownames(x))
)

# 2. Refit after omitting each sample
run_results <- vector("list", ncol(x))

for (i in seq_len(ncol(x))) {
  
  remaining <- seq_len(ncol(x)) != i
  metadata_run <- metadata[remaining, , drop = FALSE]
  
  design_run <- model.matrix(
    ~ 0 + condition,
    data = metadata_run
  )
  
  colnames(design_run) <- c("OA", "RA")
  rownames(design_run) <- rownames(metadata_run)
  
  contrast_run <- makeContrasts(
    RA_vs_OA = RA - OA,
    levels = design_run
  )
  
  fit_run <- lmFit(x[, remaining, drop = FALSE], design_run)
  fit_run <- contrasts.fit(fit_run, contrast_run)
  fit_run <- eBayes(fit_run)
  
  # Adjust across ALL tested probes in this refit.
  adjusted_p <- p.adjust(
    fit_run$p.value[, "RA_vs_OA"],
    method = "BH"
  )
  
  probe_rows <- match(
    chosen$probe_id,
    rownames(fit_run$coefficients)
  )
  
  run_results[[i]] <- data.frame(
    omitted_sample = colnames(x)[i],
    omitted_condition = as.character(metadata$condition[i]),
    gene_symbol = chosen$gene_symbol,
    probe_id = chosen$probe_id,
    logFC = fit_run$coefficients[probe_rows, "RA_vs_OA"],
    adj.P.Val = adjusted_p[probe_rows],
    stringsAsFactors = FALSE
  )
  
  cat("Completed", i, "of", ncol(x), "refits\n")
}

sensitivity_results <- do.call(rbind, run_results)

stopifnot(
  nrow(sensitivity_results) == nrow(chosen) * ncol(x),
  all(is.finite(sensitivity_results$logFC)),
  all(is.finite(sensitivity_results$adj.P.Val))
)

# 3. Summarize each selected gene
gene_summaries <- lapply(seq_len(nrow(chosen)), function(j) {
  
  gene_runs <- sensitivity_results[
    sensitivity_results$probe_id == chosen$probe_id[j],
    ,
    drop = FALSE
  ]
  
  stopifnot(nrow(gene_runs) == 20L)
  
  data.frame(
    gene_symbol = chosen$gene_symbol[j],
    probe_id = chosen$probe_id[j],
    full_data_logFC = chosen$logFC[j],
    full_data_FDR = chosen$adj.P.Val[j],
    minimum_LOO_logFC = min(gene_runs$logFC),
    maximum_LOO_logFC = max(gene_runs$logFC),
    same_direction_all_runs = all(
      sign(gene_runs$logFC) == sign(chosen$logFC[j])
    ),
    runs_FDR05 = sum(gene_runs$adj.P.Val < 0.05),
    runs_FDR05_logFC1 = sum(
      gene_runs$adj.P.Val < 0.05 &
        abs(gene_runs$logFC) >= 1
    ),
    total_runs = nrow(gene_runs)
  )
})

sensitivity_summary <- do.call(rbind, gene_summaries)

# 4. Save results
write.csv(
  sensitivity_results,
  "results/sensitivity/selected_genes_leave_one_out_all_runs.csv",
  row.names = FALSE
)

summary_path <- paste0(
  "results/sensitivity/",
  "selected_genes_leave_one_out_summary.csv"
)

write.csv(
  sensitivity_summary,
  summary_path,
  row.names = FALSE
)

writeLines(
  c(
    "Leave-one-out sensitivity analysis of eight selected genes.",
    "Twenty limma refits, omitting one OA or RA sample each time.",
    "Original RMA normalization retained.",
    "BH adjustment recalculated across all probes within each refit.",
    "Genes selected after reviewing the full-data results.",
    "This analysis is not independent validation.",
    "",
    capture.output(sessionInfo())
  ),
  "docs/session_info_leave_one_out.txt"
)

cat("\n--- Leave-one-out sensitivity summary ---\n")
print(sensitivity_summary, row.names = FALSE, digits = 4)

cat("\nSummary saved:", file.exists(summary_path), "\n")