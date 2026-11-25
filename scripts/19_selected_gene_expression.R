# scripts/19_selected_gene_expression.R
# Inspect individual-sample expression for eight selected genes.
# These genes were selected after reviewing the analysis.
# Plots are descriptive; no additional significance tests are added.

library(Biobase)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# 1. Load expression and representative probes
eset <- readRDS(
  "data_processed/GSE55235_RMA_expression_set.rds"
)

representatives <- read.csv(
  "results/enrichment/RA_vs_OA_representative_probes.csv",
  stringsAsFactors = FALSE,
  colClasses = c(entrez_id = "character")
)

x <- exprs(eset)
metadata <- pData(eset)

stopifnot(
  identical(colnames(x), rownames(metadata)),
  "condition" %in% names(metadata)
)

# 2. Select OA and RA samples
selected <- as.character(metadata$condition) %in% c("OA", "RA")
metadata <- metadata[selected, , drop = FALSE]

metadata$condition <- factor(
  as.character(metadata$condition),
  levels = c("OA", "RA")
)

x <- x[, rownames(metadata), drop = FALSE]

stopifnot(
  sum(metadata$condition == "OA") == 10L,
  sum(metadata$condition == "RA") == 10L
)

# 3. Explicit, exploratory gene selection
genes_to_plot <- c(
  "PTPRC", "CXCL13", "CCL5", "XBP1",
  "SFRP1", "APOE", "MGP", "ANKH"
)

# If a symbol occurs for multiple Entrez IDs, prefer the
# representative probe with highest average expression.
representatives <- representatives[
  order(-representatives$AveExpr, representatives$probe_id),
  ,
  drop = FALSE
]

gene_rows <- match(genes_to_plot, representatives$gene_symbol)

if (anyNA(gene_rows)) {
  stop(
    paste(
      "Missing selected genes:",
      paste(genes_to_plot[is.na(gene_rows)], collapse = ", ")
    )
  )
}

chosen <- representatives[gene_rows, , drop = FALSE]

stopifnot(
  all(chosen$probe_id %in% rownames(x)),
  !anyDuplicated(chosen$probe_id)
)

# 4. Prepare individual-sample expression
expression_selected <- x[
  chosen$probe_id,
  ,
  drop = FALSE
]

n_genes <- nrow(expression_selected)
n_samples <- ncol(expression_selected)

plot_data <- data.frame(
  gene_symbol = rep(chosen$gene_symbol, times = n_samples),
  probe_id = rep(chosen$probe_id, times = n_samples),
  sample_id = rep(colnames(expression_selected), each = n_genes),
  condition = rep(metadata$condition, each = n_genes),
  log2_expression = as.vector(expression_selected),
  stringsAsFactors = FALSE
)

plot_data$condition <- factor(
  plot_data$condition,
  levels = c("OA", "RA")
)

plot_data$gene_symbol <- factor(
  plot_data$gene_symbol,
  levels = genes_to_plot
)

stopifnot(
  nrow(plot_data) == 160L,
  all(is.finite(plot_data$log2_expression))
)

# 5. Draw expression distributions
expression_plot <- ggplot(
  plot_data,
  aes(
    x = condition,
    y = log2_expression,
    fill = condition
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.45
  ) +
  geom_point(
    aes(colour = condition),
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 42
    ),
    size = 1.8,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ gene_symbol,
    ncol = 4,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(OA = "#0072B2", RA = "#D55E00")
  ) +
  scale_colour_manual(
    values = c(OA = "#0072B2", RA = "#D55E00")
  ) +
  labs(
    title = "GSE55235: selected gene expression",
    subtitle = "RMA-normalized expression; each point represents one sample",
    x = NULL,
    y = "Log2 expression",
    caption = paste(
      "Exploratory genes selected after reviewing differential expression",
      "and enrichment results. Panels use separate y-axis scales."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

# 6. Save figures and supporting tables
png_path <- "figures/RA_vs_OA_selected_gene_expression.png"
pdf_path <- "figures/RA_vs_OA_selected_gene_expression.pdf"

ggsave(
  png_path,
  plot = expression_plot,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  pdf_path,
  plot = expression_plot,
  width = 12,
  height = 7,
  bg = "white"
)

write.csv(
  plot_data,
  "results/RA_vs_OA_selected_gene_expression_values.csv",
  row.names = FALSE
)

write.csv(
  chosen,
  "results/RA_vs_OA_selected_gene_probe_statistics.csv",
  row.names = FALSE
)

print(expression_plot)

cat("\n--- Selected gene expression summary ---\n")
cat("Genes plotted:", n_genes, "\n")
cat("Samples per gene:", n_samples, "\n")
cat("PNG saved:", file.exists(png_path), "\n")
cat("PDF saved:", file.exists(pdf_path), "\n")

cat("\nSelected representative probe statistics:\n")
print(
  chosen[, c("gene_symbol", "probe_id", "logFC", "adj.P.Val")],
  row.names = FALSE
)