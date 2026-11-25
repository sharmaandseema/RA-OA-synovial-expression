# scripts/21_write_analysis_summary.R
# Generate a Markdown report from saved analysis outputs.

dir.create("docs", showWarnings = FALSE)

# 1. Read saved results
de <- read.csv(
  "results/RA_vs_OA_limma_annotated_all_probes.csv",
  stringsAsFactors = FALSE
)

selected_genes <- read.csv(
  "results/enrichment/RA_vs_OA_selected_genes.csv",
  stringsAsFactors = FALSE
)

background <- read.csv(
  "results/enrichment/RA_vs_OA_background_genes.csv"
)

go_summary <- read.csv(
  "results/enrichment/GO_BP_analysis_summary.csv",
  stringsAsFactors = FALSE
)

go_higher <- read.csv(
  "results/enrichment/GO_BP_higher_in_RA_terms_with_gene_symbols.csv",
  stringsAsFactors = FALSE
)

go_lower <- read.csv(
  "results/enrichment/GO_BP_lower_in_RA_terms_with_gene_symbols.csv",
  stringsAsFactors = FALSE
)

sensitivity <- read.csv(
  "results/sensitivity/selected_genes_leave_one_out_summary.csv",
  stringsAsFactors = FALSE
)

# 2. Helper for Markdown tables
md_table <- function(tab) {
  
  clean_text <- function(values) {
    values <- as.character(values)
    values[is.na(values)] <- ""
    values <- gsub("|", "/", values, fixed = TRUE)
    gsub("[\r\n]+", " ", values)
  }
  
  tab[] <- lapply(tab, clean_text)
  
  c(
    paste0("| ", paste(names(tab), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(tab)), collapse = " | "), " |"),
    apply(
      tab,
      1,
      function(row) {
        paste0("| ", paste(row, collapse = " | "), " |")
      }
    )
  )
}

# 3. Calculate counts from saved tables
fdr05 <- !is.na(de$adj.P.Val) & de$adj.P.Val < 0.05
selected_probes <- fdr05 & abs(de$logFC) >= 1

probe_summary <- data.frame(
  Measure = c(
    "Probe sets tested",
    "Probe sets with FDR < 0.05",
    "Probe sets with FDR < 0.05 and abs(log2FC) >= 1",
    "Selected probe sets higher in RA",
    "Selected probe sets lower in RA",
    "Background genes for enrichment",
    "Selected representative-probe genes for enrichment"
  ),
  Count = c(
    nrow(de),
    sum(fdr05),
    sum(selected_probes),
    sum(selected_probes & de$logFC > 0),
    sum(selected_probes & de$logFC < 0),
    nrow(background),
    nrow(selected_genes)
  )
)

# 4. Prepare concise report tables
prepare_go_table <- function(tab) {
  
  tab <- tab[order(tab$p.adjust, tab$pvalue), , drop = FALSE]
  tab <- head(tab, 10)
  
  data.frame(
    GO_ID = tab$ID,
    Biological_process = tab$Description,
    Gene_ratio = tab$GeneRatio,
    FDR = formatC(tab$p.adjust, format = "e", digits = 2),
    Genes = tab$Count,
    check.names = FALSE
  )
}

sensitivity_table <- data.frame(
  Gene = sensitivity$gene_symbol,
  Full_log2FC = round(sensitivity$full_data_logFC, 3),
  Minimum_LOO_log2FC = round(sensitivity$minimum_LOO_logFC, 3),
  Maximum_LOO_log2FC = round(sensitivity$maximum_LOO_logFC, 3),
  Same_direction_all_runs = sensitivity$same_direction_all_runs,
  Runs_passing_both_thresholds = paste0(
    sensitivity$runs_FDR05_logFC1,
    "/",
    sensitivity$total_runs
  )
)

# 5. Assemble report
report <- c(
  "# RA versus OA synovial gene-expression analysis",
  "",
  paste("Report generated:", Sys.Date()),
  "",
  "## Study and comparison",
  "",
  paste(
    "Dataset: GSE55235, platform GPL96 (Affymetrix HG-U133A).",
    "The dataset contains 10 Healthy, 10 OA and 10 RA samples.",
    "The primary comparison uses the 10 RA and 10 OA samples",
    "as independent groups."
  ),
  "",
  "## Analysis methods",
  "",
  paste(
    "Expression was normalized using RMA.",
    "PCA was used for exploratory sample assessment.",
    "Limma fitted a group-only model with the RA-minus-OA contrast",
    "to RMA log2 expression.",
    "Positive log2 fold changes indicate higher expression in RA.",
    "Benjamini-Hochberg adjustment was applied across tested probes."
  ),
  "",
  paste(
    "Probe annotations were obtained from hgu133a.db.",
    "The enrichment analysis retained probes mapping to one Entrez gene.",
    "One representative probe per gene was selected by highest",
    "average expression across the RA and OA samples,",
    "independently of statistical significance.",
    "Selection then required the original probe-level FDR < 0.05",
    "and absolute log2 fold change >= 1."
  ),
  "",
  paste(
    "GO Biological Process over-representation analysis was performed",
    "separately for genes higher and lower in RA.",
    "The background comprised all eligible representative-probe genes;",
    "only GO-annotated genes contribute to the enrichment calculation.",
    "Term sizes were restricted to 10-500 background genes.",
    "BH correction was applied separately within each direction."
  ),
  "",
  "## Differential expression summary",
  "",
  md_table(probe_summary),
  "",
  "Probe counts and gene counts are distinct quantities.",
  "",
  "## GO enrichment summary",
  "",
  md_table(go_summary),
  "",
  "### Genes higher in RA: top 10 significant terms",
  "",
  md_table(prepare_go_table(go_higher)),
  "",
  "### Genes lower in RA: top 10 significant terms",
  "",
  md_table(prepare_go_table(go_lower)),
  "",
  paste(
    "GO terms overlap and must not be interpreted as independent findings.",
    "Enrichment indicates over-representation of annotated genes;",
    "it does not establish activation or inhibition of a process."
  ),
  "",
  "## Selected-gene sensitivity analysis",
  "",
  paste(
    "Eight genes were selected after reviewing the analysis.",
    "Limma was refitted 20 times, omitting one sample per refit.",
    "BH adjustment was recalculated across all probes in each refit.",
    "Original RMA normalization was held fixed."
  ),
  "",
  md_table(sensitivity_table),
  "",
  "These overlapping refits provide internal sensitivity assessment, not independent validation.",
  "",
  "## Figures",
  "",
  "### PCA",
  "",
  "![PCA](../figures/RMA_PCA_all_samples.png)",
  "",
  "### Volcano plot",
  "",
  "![Volcano plot](../figures/RA_vs_OA_volcano.png)",
  "",
  "### Selected-gene heatmap",
  "",
  "![Heatmap](../figures/RA_vs_OA_top30_heatmap.png)",
  "",
  "### GO terms: genes higher in RA",
  "",
  "![Higher-in-RA enrichment](../figures/GO_BP_higher_in_RA_dotplot.png)",
  "",
  "### GO terms: genes lower in RA",
  "",
  "![Lower-in-RA enrichment](../figures/GO_BP_lower_in_RA_dotplot.png)",
  "",
  "### Individual-sample gene expression",
  "",
  "![Selected genes](../figures/RA_vs_OA_selected_gene_expression.png)",
  "",
  "## Interpretation limits",
  "",
  "- This is an exploratory analysis of one dataset.",
  "- The fitted model includes disease group only; potential covariate confounding remains.",
  "- Bulk expression cannot by itself distinguish changes in cell composition from changes within cells.",
  "- Gene and heatmap selections used the same samples being displayed.",
  "- Probe mapping and representative-probe selection affect the enrichment inputs.",
  "- Leave-one-out sensitivity was assessed for eight selected genes, not the full result set.",
  "- No independent cohort or experimental validation was performed.",
  "",
  "## Reproducibility",
  "",
  paste(
    "Numbered scripts are stored in scripts/.",
    "Results, figures, annotation mappings and session information",
    "are saved in the project directories.",
    "See results/enrichment/selection_rules.txt for gene-selection rules."
  )
)

# 6. Save report
report_path <- "docs/analysis_summary.md"
writeLines(report, report_path, useBytes = TRUE)

cat("\n--- Analysis report ---\n")
cat("Report saved:", file.exists(report_path), "\n")
cat("Path:", report_path, "\n")

cat("\nKey counts:\n")
print(probe_summary, row.names = FALSE)