# scripts/15_prepare_enrichment_inputs.R
# Prepare unique Entrez gene lists for enrichment analysis.
# Statistical results retain their original probe-level BH adjustment.

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

# 1. Read all annotated results
de <- read.csv(
  "results/RA_vs_OA_limma_annotated_all_probes.csv",
  stringsAsFactors = FALSE,
  colClasses = c(entrez_id = "character")
)

stopifnot(
  all(
    c(
      "probe_id", "gene_symbol", "entrez_id",
      "n_entrez_mappings", "AveExpr", "logFC", "adj.P.Val"
    ) %in% names(de)
  )
)

# 2. Keep valid, unambiguous Entrez mappings
eligible <- (
  !is.na(de$n_entrez_mappings) &
    de$n_entrez_mappings == 1L &
    !is.na(de$entrez_id) &
    nzchar(de$entrez_id) &
    is.finite(de$AveExpr) &
    is.finite(de$logFC) &
    is.finite(de$adj.P.Val) &
    de$adj.P.Val >= 0 &
    de$adj.P.Val <= 1 &
    !grepl("^AFFX", de$probe_id)
)

candidates <- de[eligible, , drop = FALSE]

stopifnot(nrow(candidates) > 0L)

# 3. Choose one representative probe per gene
# Highest average expression across the RA and OA samples.
# Probe ID provides a reproducible tie-break.
candidates <- candidates[
  order(-candidates$AveExpr, candidates$probe_id),
  ,
  drop = FALSE
]

representatives <- candidates[
  !duplicated(candidates$entrez_id),
  ,
  drop = FALSE
]

rownames(representatives) <- NULL

stopifnot(!anyDuplicated(representatives$entrez_id))

# 4. Define background and selected gene sets
# Background includes all eligible represented genes,
# including those that are not differentially expressed.
universe <- representatives$entrez_id

selected <- (
  representatives$adj.P.Val < 0.05 &
    abs(representatives$logFC) >= 1
)

genes_selected <- representatives[
  selected, , drop = FALSE
]

genes_higher_ra <- representatives[
  selected & representatives$logFC > 0,
  ,
  drop = FALSE
]

genes_lower_ra <- representatives[
  selected & representatives$logFC < 0,
  ,
  drop = FALSE
]

stopifnot(
  all(genes_selected$entrez_id %in% universe),
  length(
    intersect(
      genes_higher_ra$entrez_id,
      genes_lower_ra$entrez_id
    )
  ) == 0L,
  nrow(genes_selected) ==
    nrow(genes_higher_ra) + nrow(genes_lower_ra)
)

# 5. Save inputs with annotations and probe statistics
write.csv(
  representatives,
  "results/enrichment/RA_vs_OA_representative_probes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(entrez_id = universe),
  "results/enrichment/RA_vs_OA_background_genes.csv",
  row.names = FALSE
)

write.csv(
  genes_selected,
  "results/enrichment/RA_vs_OA_selected_genes.csv",
  row.names = FALSE
)

write.csv(
  genes_higher_ra,
  "results/enrichment/RA_vs_OA_genes_higher_in_RA.csv",
  row.names = FALSE
)

write.csv(
  genes_lower_ra,
  "results/enrichment/RA_vs_OA_genes_lower_in_RA.csv",
  row.names = FALSE
)

# 6. Save the selection rules
writeLines(
  c(
    "Analysis: RA versus OA",
    "Gene identifier: Entrez ID",
    "Only probes mapping to exactly one Entrez gene retained.",
    "AFFX control probes excluded.",
    paste(
      "Representative probe: highest AveExpr across RA/OA;",
      "ties resolved by probe ID."
    ),
    "Selection: original probe-level BH FDR < 0.05 and |log2FC| >= 1.",
    "BH-adjusted values were not recalculated after probe selection.",
    "Background: all eligible genes with a representative probe.",
    "Higher and lower expression lists prepared separately."
  ),
  "results/enrichment/selection_rules.txt"
)

# 7. Print summary
cat("\n--- Enrichment input summary ---\n")
cat("Original probe sets:", nrow(de), "\n")
cat("Eligible probe sets:", nrow(candidates), "\n")
cat("Unique background genes:", length(universe), "\n")
cat("Selected genes:", nrow(genes_selected), "\n")
cat("Higher in RA:", nrow(genes_higher_ra), "\n")
cat("Lower in RA:", nrow(genes_lower_ra), "\n")

cat(
  "\nBackground file saved:",
  file.exists(
    "results/enrichment/RA_vs_OA_background_genes.csv"
  ),
  "\n"
)

cat("\nFirst 10 selected genes:\n")
print(
  head(
    genes_selected[
      order(genes_selected$adj.P.Val),
      c("gene_symbol", "entrez_id", "probe_id", "logFC", "adj.P.Val")
    ],
    10
  ),
  row.names = FALSE
)