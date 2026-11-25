# Inspect expression distributions without changing the imported data.

eset <- readRDS(
  "data_processed/GSE55235_GPL96_expression_set.rds"
)

manifest <- read.csv(
  "data_processed/GSE55235_sample_manifest.csv",
  stringsAsFactors = FALSE
)

x <- Biobase::exprs(eset)

stopifnot(
  identical(colnames(x), manifest$sample_id),
  all(is.finite(x)),
  all(x >= 0)
)

# Summarize each sample on the original scale.
q <- t(apply(
  x, 2, quantile,
  probs = c(0, 0.25, 0.5, 0.75, 1),
  names = FALSE
))

qc <- data.frame(
  sample_id = manifest$sample_id,
  condition = manifest$condition,
  minimum = q[, 1],
  percentile_25 = q[, 2],
  median = q[, 3],
  percentile_75 = q[, 4],
  maximum = q[, 5],
  zero_percent = 100 * colMeans(x == 0),
  row.names = NULL
)

table_path <- "results/expression_distribution_qc.csv"
write.csv(qc, table_path, row.names = FALSE)

# Transformation used only to make the distributions readable.
display_values <- log2(x + 1)

group_colours <- c(
  Healthy = "#999999",
  OA = "#0072B2",
  RA = "#D55E00"
)

draw_distribution <- function() {
  old_par <- par(mar = c(8, 5, 4, 1))
  on.exit(par(old_par))
  
  boxplot(
    display_values,
    names = sub("GSM", "", manifest$sample_id),
    col = unname(group_colours[manifest$condition]),
    outline = FALSE,
    las = 2,
    cex.axis = 0.7,
    ylab = "log2(expression value + 1)",
    main = "GSE55235: sample expression distributions"
  )
  
  legend(
    "topright",
    legend = names(group_colours),
    fill = group_colours,
    bty = "n"
  )
}

figure_path <- "figures/expression_distribution_qc.png"

png(figure_path, width = 1800, height = 1100, res = 150)
draw_distribution()
dev.off()

draw_distribution()

print(qc, row.names = FALSE, digits = 4)

cat("\nQC table saved:", file.exists(table_path), "\n")
cat("QC figure saved:", file.exists(figure_path), "\n")