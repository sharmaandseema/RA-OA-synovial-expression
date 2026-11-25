# RA versus OA synovial gene expression reanalysis

## Original study and project scope

This project is an exploratory reanalysis of GSE55235,
originally published by Woetzel et al. (2014).

The original study used MAS5 preprocessing and rule-based
classifiers across multiple centres. This project applies
RMA preprocessing, limma differential expression, and GO
enrichment to the Berlin cohort.

Original publication:
Woetzel D, Huber R, Kupfer P, et al.
Identification of rheumatoid arthritis and osteoarthritis
patients by transcriptome-based rule set generation.
Arthritis Research & Therapy. 2014;16:R84.
https://doi.org/10.1186/ar4526

## Study design

- Dataset: 10 Healthy, 10 OA and 10 RA samples.
- Primary comparison: 10 RA versus 10 OA samples, treated as independent groups.
- Positive log2 fold changes indicate higher expression in RA.

## Main results

- Probe sets tested: 22283.
- Probe sets with BH FDR < 0.05: 1922.
- Probe sets with FDR < 0.05 and |log2FC| >= 1: 350.
- Selected probe sets higher/lower in RA: 197/153.
- Enrichment inputs: 205 selected genes; 12650 background genes.

See the [analysis report](docs/analysis_summary.md) for GO results,
selected-gene sensitivity results, figures and limitations.

## Workflow

1. Inspect sample metadata and prepare the sample manifest.
2. Process raw arrays using RMA and assess exploratory QC.
3. Fit the group-only RA-minus-OA limma comparison.
4. Annotate probes and create differential-expression figures.
5. Select representative probes and prepare enrichment inputs.
6. Run GO Biological Process enrichment separately by direction.
7. Inspect contributing genes and individual-sample expression.
8. Assess leave-one-out sensitivity for eight selected genes.
9. Generate the analysis report and project documentation.

## Running the scripts

Open the existing RStudio project so the project root is the working directory.
Run the numbered analysis scripts in order using Source.
Scripts use relative paths and saved outputs from earlier steps.

`16a_fix_enrichment_packages.R` is a troubleshooting script.
Use it only if needed to repair dependencies before running enrichment.
It is not a routine analysis step.

The script list below is an inventory, not an automated execution sequence.

- [00_setup.R](scripts/00_setup.R)
- [01_prepare_packages.R](scripts/01_prepare_packages.R)
- [02_download_geo.R](scripts/02_download_geo.R)
- [03_inspect_metadata.R](scripts/03_inspect_metadata.R)
- [04_sample_manifest.R](scripts/04_sample_manifest.R)
- [05_expression_distribution_qc.R](scripts/05_expression_distribution_qc.R)
- [06_download_raw_arrays.R](scripts/06_download_raw_arrays.R)
- [07_extract_and_match_cel.R](scripts/07_extract_and_match_cel.R)
- [08_prepare_array_packages.R](scripts/08_prepare_array_packages.R)
- [09_rma_preprocessing.R](scripts/09_rma_preprocessing.R)
- [10_rma_pca_qc.R](scripts/10_rma_pca_qc.R)
- [11_limma_ra_vs_oa.R](scripts/11_limma_ra_vs_oa.R)
- [12_annotate_limma_results.R](scripts/12_annotate_limma_results.R)
- [13_volcano_plot.R](scripts/13_volcano_plot.R)
- [14_heatmap_top_genes.R](scripts/14_heatmap_top_genes.R)
- [15_prepare_enrichment_inputs.R](scripts/15_prepare_enrichment_inputs.R)
- [16_go_bp_enrichment.R](scripts/16_go_bp_enrichment.R)
- [16a_fix_enrichment_packages.R](scripts/16a_fix_enrichment_packages.R)
- [17_go_bp_dotplots.R](scripts/17_go_bp_dotplots.R)-
- [18_go_term_gene_tables.R](scripts/18_go_term_gene_tables.R)
- [19_selected_gene_expression.R](scripts/19_selected_gene_expression.R)
- [20_leave_one_out_sensitivity.R](scripts/20_leave_one_out_sensitivity.R)
- [21_write_analysis_summary.R](scripts/21_write_analysis_summary.R)
- [22_project_documentation.R](scripts/22_project_documentation.R)

## Project folders

| Folder | Contents |
| --- | --- |
| scripts/ | Analysis and documentation scripts |
| data_processed/ | Saved expression objects and model objects |
| results/ | Differential-expression tables and plotting inputs |
| results/enrichment/ | Gene lists, background and GO results |
| results/sensitivity/ | Selected-gene leave-one-out results |
| figures/ | Saved analysis figures |
| docs/ | Report, session information and file audit |

## Interpretation and reproducibility

- Probe-level BH-adjusted results are retained after representative-probe selection.
- GO enrichment uses the eligible measured-gene background.
- GO terms overlap; their counts do not represent independent findings.
- Enrichment does not establish process activation or inhibition.
- Heatmap and example genes were selected using these same samples.
- Leave-one-out refits retain the original RMA normalization.
- Sensitivity results concern eight selected genes, not every tested probe.
- The model does not adjust for additional covariates.
- Bulk expression may reflect cell-composition differences.
- Independent validation has not been performed.

Package versions are recorded in the session-information files under docs/.
The file audit checks existence and nonzero size, not scientific correctness.

Documentation generated: 2026-11-25
