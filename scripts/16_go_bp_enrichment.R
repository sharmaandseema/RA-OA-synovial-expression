# scripts/16a_fix_enrichment_packages.R
# Retry with fresh repository indexes.

options(
  timeout = 600,
  download.file.method = "libcurl",
  BioC_mirror = "https://bioconductor.org"
)

# Remove only temporary repository-index caches.
repo_cache <- list.files(
  tempdir(),
  pattern = "^repos_.*\\.rds$",
  full.names = TRUE
)

cat("Temporary repository caches found:", length(repo_cache), "\n")

if (length(repo_cache) > 0L) {
  unlink(repo_cache)
}

# Restore the standard repositories for this Bioconductor version.
options(
  repos = BiocManager::repositories()
)

cat("\nRepositories:\n")
print(getOption("repos"))

# Install and check dependencies in order.
packages <- c("DOSE", "enrichplot", "clusterProfiler")

for (pkg in packages) {
  
  cat("\n--- Checking", pkg, "---\n")
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE,
      type = "binary"
    )
  }
  
  library(pkg, character.only = TRUE)
  cat("Successfully loaded:", pkg, "\n")
}

cat("\nAll enrichment packages loaded successfully.\n")
