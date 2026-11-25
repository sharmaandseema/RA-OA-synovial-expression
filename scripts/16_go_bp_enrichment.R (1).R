# scripts/16_go_bp_enrichment.R
# GO Biological Process over-representation analysis.
# Higher and lower expression lists are analysed separately.

library(clusterProfiler)
library(org.Hs.eg.db)

dir.create(
  "results/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create("docs", showWarnings = FALSE)

# 1. Read gene lists
read_gene_ids <- function(path) {
  tab <- read.csv(
    path,
    stringsAsFactors = FALSE,
    colClasses = c(entrez_id = "character")
  )
  
  stopifnot("entrez_id" %in% names(tab))
  
  ids <- tab$entrez_id
  unique(ids[!is.na(ids) & nzchar(ids)])
}

background <- read_gene_ids(
  "results/enrichment/RA_vs_OA_background_genes.csv"
)

higher_ra <- read_gene_ids(
  "results/enrichment/RA_vs_OA_genes_higher_in_RA.csv"
)

lower_ra <- read_gene_ids(
  "results/enrichment/RA_vs_OA_genes_lower_in_RA.csv"
)

stopifnot(
  length(background) > 0L,
  all(higher_ra %in% background),
  all(lower_ra %in% background),
  length(intersect(higher_ra, lower_ra)) == 0L
)

cat("\n--- Input gene counts ---\n")
cat("Background genes:", length(background), "\n")
cat("Higher in RA:", length(higher_ra), "\n")
cat("Lower in RA:", length(lower_ra), "\n")

# 2. Define a reusable enrichment function
run_go_bp <- function(genes, label) {
  
  cat("\n--- GO BP:", label, "---\n")
  
  prefix <- file.path(
    "results/enrichment",
    paste0("GO_BP_", label)
  )
  
  go_result <- NULL
  
  # Empty table structure also allows empty results to be saved.
  result_table <- data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer()
  )
  
  if (length(genes) > 0L) {
    
    go_result <- clusterProfiler::enrichGO(
      gene = genes,
      universe = background,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE
    )
    
    if (!is.null(go_result)) {
      # Preserve the complete tested-term table before filtering.
      result_table <- go_result@result
    }
  }
  
  if (nrow(result_table) > 0L) {
    result_table <- result_table[
      order(result_table$p.adjust, result_table$pvalue),
      ,
      drop = FALSE
    ]
  }
  
  # Apply the reporting threshold explicitly.
  significant_terms <- result_table[
    !is.na(result_table$p.adjust) &
      result_table$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  write.csv(
    result_table,
    paste0(prefix, "_all_tested_terms.csv"),
    row.names = FALSE
  )
  
  write.csv(
    significant_terms,
    paste0(prefix, "_FDR05.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    go_result,
    paste0(prefix, "_enrichResult.rds")
  )
  
  cat("Input genes:", length(genes), "\n")
  cat("GO terms tested:", nrow(result_table), "\n")
  cat("GO terms with FDR < 0.05:", nrow(significant_terms), "\n")
  
  if (nrow(significant_terms) > 0L) {
    cat("\nTop 10 significant terms:\n")
    
    print(
      head(
        significant_terms[
          , c("ID", "Description", "GeneRatio", "p.adjust", "Count")
        ],
        10
      ),
      row.names = FALSE
    )
  } else {
    cat("No GO BP terms meet FDR < 0.05.\n")
  }
  
  data.frame(
    direction = label,
    input_genes = length(genes),
    terms_tested = nrow(result_table),
    terms_FDR05 = nrow(significant_terms)
  )
}

# 3. Run each direction separately
summary_higher <- run_go_bp(higher_ra, "higher_in_RA")
summary_lower <- run_go_bp(lower_ra, "lower_in_RA")

go_summary <- rbind(summary_higher, summary_lower)

write.csv(
  go_summary,
  "results/enrichment/GO_BP_analysis_summary.csv",
  row.names = FALSE
)

# 4. Record analysis settings and package versions
writeLines(
  c(
    "Analysis: GO Biological Process over-representation",
    "Input identifier: Entrez ID",
    "Directions analysed separately: higher and lower in RA",
    "Background: saved representative-probe background genes",
    "Only GO-annotated background genes contribute to enrichment.",
    "GO term size limits: 10 to 500 background genes",
    "Multiple testing: BH adjustment separately within each direction",
    "Reported significance: adjusted P value < 0.05",
    "",
    capture.output(sessionInfo())
  ),
  "docs/session_info_GO_BP.txt"
)

cat("\n--- GO BP analysis summary ---\n")
print(go_summary, row.names = FALSE)

cat(
  "\nSummary saved:",
  file.exists(
    "results/enrichment/GO_BP_analysis_summary.csv"
  ),
  "\n"
)