# Install the R packages required by the Gümüşler MPM reproducibility workflow.
# Run once from the repository root:
#   source("00_install_packages.R")

required_packages <- c(
  "sp",
  "gstat",
  "openxlsx",
  "raster",
  "terra",
  "sf",
  "ranger",
  "xgboost",
  "e1071",
  "pROC",
  "ggplot2",
  "scales"
)

installed <- rownames(installed.packages())
missing <- setdiff(required_packages, installed)

if (length(missing) == 0L) {
  message("All required R packages are already installed.")
} else {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
}

message("\nPackage check:")
for (pkg in required_packages) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (ok) as.character(utils::packageVersion(pkg)) else "NOT INSTALLED"
  message(sprintf("  %-12s %s", pkg, ver))
}

message("\nBefore archiving the final release, save sessionInfo() to docs/sessionInfo.txt.")
