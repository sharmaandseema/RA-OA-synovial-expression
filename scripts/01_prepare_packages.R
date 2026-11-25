
# 01_prepare_packages.R
# Install missing packages and record their versions.

stopifnot(
  file.exists("ra-oa-synovial-expression.Rproj"),
  dir.exists("docs")
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

packages <- c("GEOquery", "Biobase", "limma")

available <- vapply(
  packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (any(!available)) {
  BiocManager::install(
    packages[!available],
    ask = FALSE,
    update = FALSE
  )
}

available <- vapply(
  packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

versions <- vapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}, character(1))

package_checks <- data.frame(
  package = packages,
  available = unname(available),
  version = unname(versions)
)

write.csv(
  package_checks,
  "docs/package_versions_setup.csv",
  row.names = FALSE
)

print(package_checks, row.names = FALSE)

stopifnot(all(available))