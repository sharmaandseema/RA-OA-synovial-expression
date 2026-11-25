# scripts/13_volcano_plot.R
# Volcano plot: RA versus OA, at probe-set level

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# 1. Read annotated results
volcano_data <- read.csv(
  "results/RA_vs_OA_limma_annotated_all_probes.csv",
  stringsAsFactors = FALSE
)

stopifnot(
  all(
    c("probe_id", "logFC", "adj.P.Val") %in%
      names(volcano_data)
  )
)

# 2. Set display thresholds
fdr_threshold <- 0.05
logfc_threshold <- 1

valid <- (
  is.finite(volcano_data$logFC) &
    is.finite(volcano_data$adj.P.Val) &
    volcano_data$adj.P.Val >= 0 &
    volcano_data$adj.P.Val <= 1
)

cat("Probe sets excluded for invalid values:", sum(!valid), "\n")

volcano_data <- volcano_data[valid, , drop = FALSE]
stopifnot(nrow(volcano_data) > 0L)

# Protect against log10(0).
volcano_data$minus_log10_fdr <- -log10(
  pmax(volcano_data$adj.P.Val, .Machine$double.xmin)
)

volcano_data$category <- "Below display thresholds"

higher_ra <- (
  volcano_data$adj.P.Val < fdr_threshold &
    volcano_data$logFC >= logfc_threshold
)

lower_ra <- (
  volcano_data$adj.P.Val < fdr_threshold &
    volcano_data$logFC <= -logfc_threshold
)

volcano_data$category[higher_ra] <- "Higher in RA"
volcano_data$category[lower_ra] <- "Lower in RA"

plot_colors <- c(
  "Below display thresholds" = "grey75",
  "Higher in RA" = "#D55E00",
  "Lower in RA" = "#0072B2"
)

n_higher <- sum(higher_ra)
n_lower <- sum(lower_ra)

# 3. Define plotting function
draw_volcano <- function() {
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Draw highlighted points last.
  point_order <- order(
    volcano_data$category != "Below display thresholds"
  )
  
  x_limit <- max(abs(volcano_data$logFC)) * 1.08
  y_limit <- max(volcano_data$minus_log10_fdr) * 1.12
  
  plot(
    volcano_data$logFC[point_order],
    volcano_data$minus_log10_fdr[point_order],
    col = unname(
      plot_colors[volcano_data$category[point_order]]
    ),
    pch = 16,
    cex = 0.65,
    xlim = c(-x_limit, x_limit),
    ylim = c(0, y_limit),
    xlab = "Log2 fold change (RA minus OA)",
    ylab = expression(-log[10]("BH-adjusted P value")),
    main = "GSE55235: RA versus OA",
    las = 1
  )
  
  abline(
    v = c(-logfc_threshold, logfc_threshold),
    h = -log10(fdr_threshold),
    lty = 2,
    col = "grey40"
  )
  
  legend(
    "topright",
    legend = c(
      paste0("Higher in RA (", n_higher, ")"),
      paste0("Lower in RA (", n_lower, ")"),
      "Below display thresholds"
    ),
    col = unname(
      plot_colors[
        c(
          "Higher in RA",
          "Lower in RA",
          "Below display thresholds"
        )
      ]
    ),
    pch = 16,
    bty = "n",
    cex = 0.85
  )
  
  mtext(
    "Highlighted: FDR < 0.05 and |log2 fold change| >= 1",
    side = 3,
    line = 0.35,
    cex = 0.8
  )
}

# 4. Save PNG and PDF
png_path <- "figures/RA_vs_OA_volcano.png"
pdf_path <- "figures/RA_vs_OA_volcano.pdf"

png(
  png_path,
  width = 2400,
  height = 2000,
  res = 300
)
tryCatch(draw_volcano(), finally = dev.off())

pdf(pdf_path, width = 8, height = 6.7)
tryCatch(draw_volcano(), finally = dev.off())

# 5. Save the plotted data
write.csv(
  volcano_data,
  "results/RA_vs_OA_volcano_data.csv",
  row.names = FALSE
)

# 6. Show the plot in RStudio
draw_volcano()

cat("\n--- Volcano plot summary ---\n")
cat("Probe sets plotted:", nrow(volcano_data), "\n")
cat("Higher in RA:", n_higher, "\n")
cat("Lower in RA:", n_lower, "\n")
cat("Total highlighted:", n_higher + n_lower, "\n")
cat("PNG saved:", file.exists(png_path), "\n")
cat("PDF saved:", file.exists(pdf_path), "\n")