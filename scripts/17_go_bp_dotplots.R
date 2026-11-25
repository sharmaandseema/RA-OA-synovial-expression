# scripts/17_go_bp_dotplots.R
# Plot significant GO BP enrichment results.

library(ggplot2)

dir.create("figures", showWarnings = FALSE)

# 1. Read significant results
higher <- read.csv(
  "results/enrichment/GO_BP_higher_in_RA_FDR05.csv",
  stringsAsFactors = FALSE
)

lower <- read.csv(
  "results/enrichment/GO_BP_lower_in_RA_FDR05.csv",
  stringsAsFactors = FALSE
)

# 2. Prepare the terms to display
prepare_terms <- function(tab, number = 15L) {
  
  required <- c(
    "ID", "Description", "GeneRatio", "p.adjust", "Count"
  )
  
  stopifnot(
    all(required %in% names(tab)),
    nrow(tab) > 0L,
    all(is.finite(tab$p.adjust)),
    all(tab$p.adjust < 0.05)
  )
  
  tab <- tab[
    order(tab$p.adjust, tab$pvalue, tab$ID),
    ,
    drop = FALSE
  ]
  
  tab <- head(tab, number)
  
  tab$gene_ratio_numeric <- vapply(
    strsplit(as.character(tab$GeneRatio), "/", fixed = TRUE),
    function(parts) {
      values <- as.numeric(parts)
      values[1] / values[2]
    },
    numeric(1)
  )
  
  stopifnot(
    all(is.finite(tab$gene_ratio_numeric)),
    all(tab$gene_ratio_numeric > 0),
    all(tab$gene_ratio_numeric <= 1)
  )
  
  tab$minus_log10_fdr <- -log10(
    pmax(tab$p.adjust, .Machine$double.xmin)
  )
  
  # Most significant term appears at the top.
  tab$term_id <- factor(tab$ID, levels = rev(tab$ID))
  
  tab
}

higher_plot_data <- prepare_terms(higher, 15L)
lower_plot_data <- prepare_terms(lower, 15L)

# Shared scales make the two figures easier to compare.
all_plot_data <- rbind(higher_plot_data, lower_plot_data)

max_ratio <- min(
  1,
  max(all_plot_data$gene_ratio_numeric) * 1.15
)

max_count <- max(all_plot_data$Count)
max_score <- max(all_plot_data$minus_log10_fdr)

# 3. Define plotting and saving function
make_dotplot <- function(tab, direction, title) {
  
  term_labels <- setNames(
    vapply(
      tab$Description,
      function(label) {
        paste(strwrap(label, width = 42), collapse = "\n")
      },
      character(1)
    ),
    tab$ID
  )
  
  plot_object <- ggplot(
    tab,
    aes(
      x = gene_ratio_numeric,
      y = term_id,
      size = Count,
      colour = minus_log10_fdr
    )
  ) +
    geom_point() +
    scale_y_discrete(labels = term_labels) +
    scale_x_continuous(
      limits = c(0, max_ratio),
      labels = function(values) {
        paste0(round(values * 100), "%")
      }
    ) +
    scale_size_area(
      name = "Overlapping genes",
      max_size = 9,
      limits = c(0, max_count)
    ) +
    scale_colour_gradient(
      name = "-log10(FDR)",
      low = "#4575B4",
      high = "#D73027",
      limits = c(0, max_score)
    ) +
    labs(
      title = title,
      subtitle = paste(
        nrow(tab),
        "terms shown; ranked by BH-adjusted P value"
      ),
      x = "Gene ratio",
      y = NULL,
      caption = paste(
        "Gene ratio: overlapping genes / GO-annotated input genes.",
        "BH adjustment performed separately within each direction.",
        sep = "\n"
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0),
      legend.position = "right"
    )
  
  prefix <- file.path(
    "figures",
    paste0("GO_BP_", direction, "_dotplot")
  )
  
  plot_height <- max(6, nrow(tab) * 0.55 + 1.6)
  
  ggsave(
    paste0(prefix, ".png"),
    plot = plot_object,
    width = 11,
    height = plot_height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    paste0(prefix, ".pdf"),
    plot = plot_object,
    width = 11,
    height = plot_height,
    units = "in",
    bg = "white"
  )
  
  write.csv(
    tab,
    file.path(
      "results/enrichment",
      paste0("GO_BP_", direction, "_plotted_terms.csv")
    ),
    row.names = FALSE
  )
  
  print(plot_object)
  
  cat("\n---", direction, "---\n")
  cat("Terms plotted:", nrow(tab), "\n")
  cat("PNG saved:", file.exists(paste0(prefix, ".png")), "\n")
  cat("PDF saved:", file.exists(paste0(prefix, ".pdf")), "\n")
  
  invisible(plot_object)
}

# 4. Generate both plots
go_plot_higher <- make_dotplot(
  higher_plot_data,
  "higher_in_RA",
  "GO Biological Process: genes higher in RA"
)

go_plot_lower <- make_dotplot(
  lower_plot_data,
  "lower_in_RA",
  "GO Biological Process: genes lower in RA"
)