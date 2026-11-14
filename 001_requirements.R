# Install and load packages.
packages <- c("glmnet", "Matrix", "knitr", "rmarkdown")
missing_packages <- setdiff(packages, rownames(installed.packages()))
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
source("002_requirements_lite.R")
