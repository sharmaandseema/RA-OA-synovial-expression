# Prepare packages for Affymetrix U133A array processing.

stopifnot(
  file.exists("ra-oa-synovial-expression.Rproj"),
  dir.exists("docs")
)

packages <- c("affy", "hgu133acdf")

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

versions <- vapply(
  packages,
  function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  },
  character(1)
)

package_checks <- data.frame(
  package = packages,
  available = unname(available),
  version = unname(versions)
)

write.csv(
  package_checks,
  "docs/array_package_versions.csv",
  row.names = FALSE
)

print(package_checks, row.names = FALSE)
stopifnot(all(available))