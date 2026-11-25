# 00_setup.R
# Set up the RA–OA synovial gene-expression project.

stopifnot(
  file.exists("ra-oa-synovial-expression.Rproj")
)

folders <- c(
  "scripts",
  "data_raw",
  "data_processed",
  "results",
  "figures",
  "docs"
)

for (folder in folders) {
  dir.create(folder, showWarnings = FALSE)
}

writeLines(
  capture.output(sessionInfo()),
  "docs/session_info_setup.txt"
)

print(data.frame(
  check = c(
    "Setup script saved",
    "Session record saved"
  ),
  passed = c(
    file.exists("scripts/00_setup.R"),
    file.exists("docs/session_info_setup.txt")
  )
))