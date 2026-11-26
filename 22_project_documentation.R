# scripts/22_project_documentation.R
# Check main deliverables and write the project README.

dir.create("docs", showWarnings = FALSE)

# 1. Check required analysis outputs
required_files <- c(
  "data_processed/GSE55235_RMA_expression_set.rds",
  "data_processed/RA_vs_OA_limma_analysis.rds",
  "results/RA_vs_OA_limma_annotated_all_probes.csv",
  "results/enrichment/RA_vs_OA_representative_probes.csv",
  "results/enrichment/RA_vs_OA_background_genes.csv",
  "results/enrichment/RA_vs_OA_selected_genes.csv",
  "results/enrichment/GO_BP_analysis_summary.csv",
  "results/enrichment/GO_BP_higher_in_RA_term_gene_evidence.csv",
  "results/enrichment/GO_BP_lower_in_RA_term_gene_evidence.csv",
  "results/sensitivity/selected_genes_leave_one_out_summary.csv",
  "figures/RMA_PCA_all_samples.png",
  "figures/RA_vs_OA_volcano.png",
  "figures/RA_vs_OA_top30_heatmap.png",
  "figures/GO_BP_higher_in_RA_dotplot.png",
  "figures/GO_BP_lower_in_RA_dotplot.png",
  "figures/RA_vs_OA_selected_gene_expression.png",
  "docs/analysis_summary.md"
)

file_audit <- data.frame(
  path = required_files,
  exists = file.exists(required_files),
  size_bytes = as.numeric(file.info(required_files)$size),
  stringsAsFactors = FALSE
)

file_audit$passed <- (
  file_audit$exists &
    !is.na(file_audit$size_bytes) &
    file_audit$size_bytes > 0
)

write.csv(
  file_audit,
  "docs/project_file_audit.csv",
  row.names = FALSE
)

if (!all(file_audit$passed)) {
  print(file_audit[!file_audit$passed, ], row.names = FALSE)
  stop("Some required outputs are missing or empty. Send the table above.")
}

# 2. Read final counts
de <- read.csv(
  "results/RA_vs_OA_limma_annotated_all_probes.csv"
)

background <- read.csv(
  "results/enrichment/RA_vs_OA_background_genes.csv"
)

selected_genes <- read.csv(
  "results/enrichment/RA_vs_OA_selected_genes.csv"
)

fdr05 <- !is.na(de$adj.P.Val) & de$adj.P.Val < 0.05
selected <- fdr05 & abs(de$logFC) >= 1

# 3. List actual project scripts
script_files <- sort(
  list.files(
    "scripts",
    pattern = "\\.[Rr]$",
    full.names = FALSE
  )
)

script_links <- paste0(
  "- [", script_files, "](scripts/", script_files, ")"
)

# 4. Prepare README
readme <- c(
  "# RA versus OA synovial gene expression reanalysis",
  "",
  "Exploratory analysis of GSE55235 using the GPL96 HG-U133A platform.",
  "",
  "## Original study and project scope",
  "",
  paste(
    "This project is an exploratory reanalysis of GSE55235,",
    "originally published by Woetzel et al. (2014).",
    "The original study used MAS5 preprocessing and rule-based",
    "classifiers across multiple centres. This project applies",
    "RMA preprocessing, limma differential expression and GO",
    "enrichment to the Berlin cohort."
  ),
  "",
  "Original publication:",
  "Woetzel D, Huber R, Kupfer P, et al.",
  "Identification of rheumatoid arthritis and osteoarthritis",
  "patients by transcriptome-based rule set generation.",
  "Arthritis Research & Therapy. 2014;16:R84.",
  "https://doi.org/10.1186/ar4526",
  "",
  "## Study design",
  "",
  "- Dataset: 10 Healthy, 10 OA and 10 RA samples.",
  "- Primary comparison: 10 RA versus 10 OA samples, treated as independent groups.",
  "- Positive log2 fold changes indicate higher expression in RA.",
  "",
  "## Main results",
  "",
  paste0("- Probe sets tested: ", nrow(de), "."),
  paste0("- Probe sets with BH FDR < 0.05: ", sum(fdr05), "."),
  paste0(
    "- Probe sets with FDR < 0.05 and |log2FC| >= 1: ",
    sum(selected), "."
  ),
  paste0(
    "- Selected probe sets higher/lower in RA: ",
    sum(selected & de$logFC > 0), "/",
    sum(selected & de$logFC < 0), "."
  ),
  paste0(
    "- Enrichment inputs: ", nrow(selected_genes),
    " selected genes; ", nrow(background), " background genes."
  ),
  "",
  "See the [analysis report](docs/analysis_summary.md) for GO results,",
  "selected-gene sensitivity results, figures and limitations.",
  "",
  "## Workflow",
  "",
  "1. Inspect sample metadata and prepare the sample manifest.",
  "2. Process raw arrays using RMA and assess exploratory QC.",
  "3. Fit the group-only RA-minus-OA limma comparison.",
  "4. Annotate probes and create differential-expression figures.",
  "5. Select representative probes and prepare enrichment inputs.",
  "6. Run GO Biological Process enrichment separately by direction.",
  "7. Inspect contributing genes and individual-sample expression.",
  "8. Assess leave-one-out sensitivity for eight selected genes.",
  "9. Generate the analysis report and project documentation.",
  "",
  "## Running the scripts",
  "",
  "Open the existing RStudio project so the project root is the working directory.",
  "Run the numbered analysis scripts in order using Source.",
  "Scripts use relative paths and saved outputs from earlier steps.",
  "",
  "`16a_fix_enrichment_packages.R` is a troubleshooting script.",
  "Use it only if needed to repair dependencies before running enrichment.",
  "It is not a routine analysis step.",
  "",
  "The script list below is an inventory, not an automated execution sequence.",
  "",
  script_links,
  "",
  "## Project folders",
  "",
  "| Folder | Contents |",
  "| --- | --- |",
  "| scripts/ | Analysis and documentation scripts |",
  "| data_processed/ | Saved expression objects and model objects |",
  "| results/ | Differential-expression tables and plotting inputs |",
  "| results/enrichment/ | Gene lists, background and GO results |",
  "| results/sensitivity/ | Selected-gene leave-one-out results |",
  "| figures/ | Saved analysis figures |",
  "| docs/ | Report, session information and file audit |",
  "",
  "## Interpretation and reproducibility",
  "",
  "- Probe-level BH-adjusted results are retained after representative-probe selection.",
  "- GO enrichment uses the eligible measured-gene background.",
  "- GO terms overlap; their counts do not represent independent findings.",
  "- Enrichment does not establish process activation or inhibition.",
  "- Heatmap and example genes were selected using these same samples.",
  "- Leave-one-out refits retain the original RMA normalization.",
  "- Sensitivity results concern eight selected genes, not every tested probe.",
  "- The model does not adjust for additional covariates.",
  "- Bulk expression may reflect cell-composition differences.",
  "- Independent validation has not been performed.",
  "",
  "Package versions are recorded in the session-information files under docs/.",
  "The file audit checks existence and nonzero size, not scientific correctness.",
  "",
  paste("Documentation generated:", Sys.Date())
)

# Preserve an existing README before replacing it.
backup_path <- "docs/README_before_documentation.md"

if (file.exists("README.md") && !file.exists(backup_path)) {
  backup_ok <- file.copy("README.md", backup_path)
  
  if (!backup_ok) {
    stop("Could not back up the existing README.")
  }
}

writeLines(readme, "README.md", useBytes = TRUE)

# 5. Print summary
cat("\n--- Project documentation summary ---\n")
cat("Required files checked:", nrow(file_audit), "\n")
cat("Missing or empty files:", sum(!file_audit$passed), "\n")
cat("R scripts found:", length(script_files), "\n")
cat("README saved:", file.exists("README.md"), "\n")
cat(
  "File audit saved:",
  file.exists("docs/project_file_audit.csv"),
  "\n"
)
