# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT 6: SUMMARY TABLE OF DATASET CHARACTERISTICS (TABLE 1 IN PAPER)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loads the processed empirical datasets from `data/individual/` and computes:
#   - Response range [min, max]
#   - Sample size (n) and number of numerical predictors (p)
#   - Global percentage of missing predictor cells (m_global)
#   - Number (and percentage) of rows containing at least one missing value (n_miss)
#
# Output: Formatted LaTeX table code corresponding to Table 1 (`tab:datasets`)
# in the manuscript.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)
library(purrr)
library(tibble)
library(knitr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. LOAD DATASETS USED IN THE PAPER
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dataset_names <- c(
  "Crime",
  "Garment",
  "Auto",
  "Automobile",
  "Ozone"
)

datasets <- lapply(dataset_names, function(name) {
  readRDS(file.path("data/individual", paste0(name, ".rds")))
})
names(datasets) <- dataset_names

# Source metadata for table annotation
source_info <- tibble::tribble(
  ~internal_name, ~source,   ~id,
  "Crime",        "UCI",     "183",
  "Garment",      "UCI",     "597",
  "Auto",         "UCI",     "9",
  "Automobile",   "UCI",     "10",
  "Ozone",        "mlbench", "-"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. COMPUTE SUMMARY METRICS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tab1 <- map_df(names(datasets), function(name) {
  d   <- datasets[[name]]
  src <- source_info %>% filter(internal_name == name)
  
  n_rows     <- nrow(d$X)
  n_cols     <- ncol(d$X)
  total_miss <- sum(is.na(as.matrix(d$X)))
  rows_miss  <- sum(apply(is.na(d$X), 1, any))
  
  data.frame(
    Dataset          = d$name,
    Source           = src$source,
    ID               = src$id,
    `Response range` = paste0(
      "$[", round(min(d$y, na.rm = TRUE), 2), ", ",
      round(max(d$y, na.rm = TRUE), 2), "]$"
    ),
    n                = n_rows,
    p                = n_cols,
    m_global         = round(100 * total_miss / (n_rows * n_cols), 2),
    n_miss           = sprintf("%d (%.2f\\%%)", rows_miss, 100 * rows_miss / n_rows),
    stringsAsFactors = FALSE
  )
}) %>%
  arrange(desc(n))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. PRINT LATEX TABLE (TABLE 1 IN PAPER)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kable(
  tab1,
  format   = "latex",
  escape   = FALSE,
  booktabs = TRUE,
  linesep  = "",
  col.names = c(
    "Dataset",
    "Source",
    "UCI ID",
    "Response range",
    "$n$",
    "$p$",
    "$m_{global}(\\%)$",
    "$n_{miss}$"
  ),
  caption  = paste(
    "Summary of the datasets used in the empirical evaluation:",
    "$n$ and $p$ denote the number of observations and the number of numerical predictors, respectively;",
    "$m_{global}$ is the percentage of missing predictor cells among all $n\\times p$ predictor values;",
    "$n_{miss}$ is the number of observations containing at least one missing predictor value."
  ),
  label    = "tab:datasets"
)
