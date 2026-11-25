# scripts/14_heatmap_top_genes.R
# Exploratory heatmap of selected RA-versus-OA results.
# Gene selection uses these same samples; this is not validation.

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(Biobase)
library(pheatmap)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# 1. Load saved data
eset <- readRDS(
  "data_processed/GSE55235_RMA_expression_set.rds"
)

de <- read.csv(
  "results/RA_vs_OA_limma_annotated_all_probes.csv",
  stringsAsFactors = FALSE,
  colClasses = c(entrez_id = "character")
)

x <- exprs(eset)
metadata <- pData(eset)

stopifnot(
  identical(colnames(x), rownames(metadata)),
  "condition" %in% colnames(metadata)
)

# 2. Select OA and RA samples, displayed in group order
selected <- as.character(metadata$condition) %in% c("OA", "RA")
metadata <- metadata[selected, , drop = FALSE]

metadata$condition <- factor(
  as.character(metadata$condition),
  levels = c("OA", "RA")
)

metadata <- metadata[
  order(metadata$condition, rownames(metadata)),
  ,
  drop = FALSE
]

stopifnot(
  sum(metadata$condition == "OA") == 10L,
  sum(metadata$condition == "RA") == 10L
)

x <- x[, rownames(metadata), drop = FALSE]

# 3. Select significant probes with one Entrez gene mapping
eligible <- (
  !is.na(de$adj.P.Val) &
    de$adj.P.Val < 0.05 &
    abs(de$logFC) >= 1 &
    !is.na(de$n_entrez_mappings) &
    de$n_entrez_mappings == 1L &
    !is.na(de$gene_symbol) &
    nzchar(de$gene_symbol) &
    !is.na(de$entrez_id) &
    nzchar(de$entrez_id) &
    de$probe_id %in% rownames(x)
)

candidates <- de[eligible, , drop = FALSE]

# Exclude constant rows before selecting genes.
probe_sd <- apply(
  x[candidates$probe_id, , drop = FALSE],
  1,
  sd
)

candidates <- candidates[
  is.finite(probe_sd) & probe_sd > 0,
  ,
  drop = FALSE
]

# Choose the most statistically significant probe per gene.
candidates <- candidates[
  order(candidates$adj.P.Val, candidates$P.Value, candidates$probe_id),
  ,
  drop = FALSE
]

candidates <- candidates[
  !duplicated(candidates$entrez_id),
  ,
  drop = FALSE
]

stopifnot(nrow(candidates) >= 2L)

top_genes <- head(candidates, 30)

# 4. Calculate row z-scores across the 20 samples
heatmap_expression <- x[top_genes$probe_id, , drop = FALSE]
heatmap_z <- t(scale(t(heatmap_expression)))

stopifnot(
  all(is.finite(heatmap_z)),
  identical(colnames(heatmap_z), rownames(metadata))
)

# Include probe IDs so each displayed row is traceable.
row_labels <- paste0(
  top_genes$gene_symbol,
  " [", top_genes$probe_id, "]"
)

annotation_col <- data.frame(
  Condition = metadata$condition,
  row.names = rownames(metadata)
)

# 5. Build the heatmap
# Columns remain in OA-then-RA order; rows are clustered.
color_limit <- max(abs(heatmap_z))

heatmap_plot <- pheatmap(
  heatmap_z,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = list(
    Condition = c(OA = "#0072B2", RA = "#D55E00")
  ),
  labels_row = row_labels,
  show_colnames = FALSE,
  gaps_col = sum(metadata$condition == "OA"),
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(101),
  breaks = seq(
    -color_limit, color_limit, length.out = 102
  ),
  border_color = NA,
  fontsize_row = 8,
  main = paste0(
    "RA versus OA: ", nrow(top_genes),
    " selected genes\nRow z-scores; samples ordered OA then RA"
  ),
  silent = TRUE
)

# 6. Save PNG and PDF
draw_heatmap <- function() {
  grid::grid.newpage()
  grid::grid.draw(heatmap_plot$gtable)
}

png_path <- "figures/RA_vs_OA_top30_heatmap.png"
pdf_path <- "figures/RA_vs_OA_top30_heatmap.pdf"

png(png_path, width = 3300, height = 3000, res = 300)
tryCatch(draw_heatmap(), finally = dev.off())

pdf(pdf_path, width = 11, height = 10)
tryCatch(draw_heatmap(), finally = dev.off())

# 7. Save selected rows and plotted values
write.csv(
  top_genes,
  "results/RA_vs_OA_heatmap_selected_genes.csv",
  row.names = FALSE
)

write.csv(
  heatmap_z,
  "results/RA_vs_OA_heatmap_row_zscores.csv",
  row.names = TRUE
)

write.csv(
  data.frame(
    sample_id = rownames(metadata),
    condition = metadata$condition
  ),
  "results/RA_vs_OA_heatmap_sample_order.csv",
  row.names = FALSE
)

# Display in RStudio
draw_heatmap()

cat("\n--- Heatmap summary ---\n")
cat("Genes displayed:", nrow(top_genes), "\n")
cat("Samples displayed:", ncol(heatmap_z), "\n")
cat("OA samples:", sum(metadata$condition == "OA"), "\n")
cat("RA samples:", sum(metadata$condition == "RA"), "\n")
cat("PNG saved:", file.exists(png_path), "\n")
cat("PDF saved:", file.exists(pdf_path), "\n")