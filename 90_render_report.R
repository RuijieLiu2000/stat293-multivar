# Knit retained results only.
library(rmarkdown)
rmarkdown::render(
  "FinalReport.Rmd",
  output_file = "FinalReport.pdf",
  clean = TRUE,
  quiet = TRUE
)
