# scripts/18_go_term_gene_tables.R
# Link enriched GO terms to their contributing genes.
# No additional statistical tests are performed.

output_dir <- "results/enrichment"

# 1. Read representative probes
representatives <- read.csv(
  file.path(output_dir, "RA_vs_OA_representative_probes.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(entrez_id = "character")
)

stopifnot(!anyDuplicated(representatives$entrez_id))

# 2. Process one direction at a time
export_term_genes <- function(direction) {
  
  terms <- read.csv(
    file.path(
      output_dir,
      paste0("GO_BP_", direction, "_FDR05.csv")
    ),
    stringsAsFactors = FALSE
  )
  
  stopifnot(nrow(terms) > 0L)
  
  terms <- terms[
    order(terms$p.adjust, terms$pvalue, terms$ID),
    ,
    drop = FALSE
  ]
  
  # geneID contains slash-separated Entrez IDs.
  members <- lapply(
    strsplit(as.character(terms$geneID), "/", fixed = TRUE),
    unique
  )
  
  stopifnot(all(lengths(members) == terms$Count))
  
  term_rows <- rep(
    seq_len(nrow(terms)),
    lengths(members)
  )
  
  gene_ids <- unlist(members, use.names = FALSE)
  
  gene_rows <- match(
    gene_ids,
    representatives$entrez_id
  )
  
  stopifnot(!anyNA(gene_rows))
  
  # One row per GO-term/gene pair.
  evidence <- data.frame(
    direction = direction,
    GO_ID = terms$ID[term_rows],
    GO_description = terms$Description[term_rows],
    GO_adjusted_P = terms$p.adjust[term_rows],
    entrez_id = gene_ids,
    gene_symbol = representatives$gene_symbol[gene_rows],
    probe_id = representatives$probe_id[gene_rows],
    logFC_RA_minus_OA = representatives$logFC[gene_rows],
    probe_adjusted_P = representatives$adj.P.Val[gene_rows],
    stringsAsFactors = FALSE
  )
  
  # Check that the contributing genes match the intended direction.
  stopifnot(
    all(evidence$probe_adjusted_P < 0.05),
    all(abs(evidence$logFC_RA_minus_OA) >= 1)
  )
  
  if (direction == "higher_in_RA") {
    stopifnot(all(evidence$logFC_RA_minus_OA > 0))
  } else {
    stopifnot(all(evidence$logFC_RA_minus_OA < 0))
  }
  
  write.csv(
    evidence,
    file.path(
      output_dir,
      paste0("GO_BP_", direction, "_term_gene_evidence.csv")
    ),
    row.names = FALSE
  )
  
  # 3. Add gene symbols to the term-level table
  terms$overlapping_gene_symbols <- vapply(
    members,
    function(ids) {
      rows <- match(ids, representatives$entrez_id)
      symbols <- representatives$gene_symbol[rows]
      
      missing_symbol <- is.na(symbols) | !nzchar(symbols)
      symbols[missing_symbol] <- paste0(
        "Entrez:", ids[missing_symbol]
      )
      
      paste(sort(unique(symbols)), collapse = "; ")
    },
    character(1)
  )
  
  write.csv(
    terms,
    file.path(
      output_dir,
      paste0("GO_BP_", direction, "_terms_with_gene_symbols.csv")
    ),
    row.names = FALSE
  )
  
  # 4. Count how many significant terms contain each gene
  # This measures annotation recurrence, not biological importance.
  gene_frequency <- as.data.frame(
    table(evidence$entrez_id),
    stringsAsFactors = FALSE
  )
  
  names(gene_frequency) <- c(
    "entrez_id", "number_of_significant_GO_terms"
  )
  
  gene_frequency$entrez_id <- as.character(
    gene_frequency$entrez_id
  )
  
  rows <- match(
    gene_frequency$entrez_id,
    representatives$entrez_id
  )
  
  gene_frequency$gene_symbol <- representatives$gene_symbol[rows]
  gene_frequency$logFC <- representatives$logFC[rows]
  
  gene_frequency <- gene_frequency[
    order(
      -gene_frequency$number_of_significant_GO_terms,
      gene_frequency$entrez_id
    ),
    ,
    drop = FALSE
  ]
  
  write.csv(
    gene_frequency,
    file.path(
      output_dir,
      paste0("GO_BP_", direction, "_gene_recurrence.csv")
    ),
    row.names = FALSE
  )
  
  cat("\n---", direction, "---\n")
  cat("Significant GO terms:", nrow(terms), "\n")
  cat("Term-gene pairs:", nrow(evidence), "\n")
  cat("Distinct contributing genes:", nrow(gene_frequency), "\n")
  
  cat("\nGenes appearing in the most significant GO terms:\n")
  print(
    head(
      gene_frequency[
        , c(
          "gene_symbol",
          "number_of_significant_GO_terms",
          "logFC"
        )
      ],
      10
    ),
    row.names = FALSE
  )
  
  invisible(evidence)
}

# 5. Export both directions
evidence_higher <- export_term_genes("higher_in_RA")
evidence_lower <- export_term_genes("lower_in_RA")

cat("\nGO term-gene tables saved successfully.\n")