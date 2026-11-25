# Inspect sample variation after RMA preprocessing.

eset <- readRDS(
  "data_processed/GSE55235_RMA_expression_set.rds"
)

x <- Biobase::exprs(eset)
metadata <- Biobase::pData(eset)

stopifnot(
  identical(colnames(x), rownames(metadata)),
  all(is.finite(x))
)

# Use all probe sets with non-zero variance.
probe_variance <- apply(x, 1, var)
keep <- is.finite(probe_variance) & probe_variance > 0

# Centre each probe set without scaling to unit variance.
pca <- prcomp(
  t(x[keep, , drop = FALSE]),
  center = TRUE,
  scale. = FALSE
)

variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)

scores <- data.frame(
  sample_id = rownames(pca$x),
  condition = metadata$condition,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  row.names = NULL
)

write.csv(
  scores,
  "results/RMA_PCA_sample_scores.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    component = paste0("PC", seq_along(variance_percent)),
    variance_percent = variance_percent
  ),
  "results/RMA_PCA_variance.csv",
  row.names = FALSE
)

group_colours <- c(
  Healthy = "#999999",
  OA = "#0072B2",
  RA = "#D55E00"
)

# Short labels retain the GEO accession suffix.
labels <- sub("^GSM133", "", scores$sample_id)

draw_pca <- function() {
  plot(
    scores$PC1,
    scores$PC2,
    col = unname(group_colours[scores$condition]),
    pch = 19,
    xlab = sprintf("PC1 (%.1f%%)", variance_percent[1]),
    ylab = sprintf("PC2 (%.1f%%)", variance_percent[2]),
    main = "GSE55235: PCA after RMA",
    sub = "All 30 samples; labels show GEO accession suffixes",
    ylim = extendrange(scores$PC2, f = 0.15)
  )
  
  text(
    scores$PC1,
    scores$PC2,
    labels = labels,
    pos = 3,
    cex = 0.65
  )
  
  legend(
    "topright",
    legend = names(group_colours),
    col = group_colours,
    pch = 19,
    bty = "n"
  )
}

figure_path <- "figures/RMA_PCA_all_samples.png"

png(figure_path, width = 1500, height = 1100, res = 150)
draw_pca()
dev.off()

draw_pca()

cat("Probe sets used:", sum(keep), "\n")
cat("\nVariance explained by the first five PCs (%):\n")
print(round(variance_percent[1:5], 2))

cat("\nPCA figure saved:", file.exists(figure_path), "\n")
