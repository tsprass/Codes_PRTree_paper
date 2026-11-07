# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT 4: DOWNLOAD AND PREPARE SELECTED DATASETS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script prepares the datasets listed below for analysis.
# Input: UCI ML Repository and the mlbench package
# Output: one RDS file per dataset in data/individual/

# Load functions
source("03 - functions.R")

# Define the selected datasets manually below.
cat("Downloading selected UCI datasets...\n")

# Fetch datasets with missing values
datasets <- list(
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # large n
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Communities and Crime (1994, 127); "?" denotes missing values.
  Crime = fetch_safe(id = 183),
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # medium n
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Productivity Prediction of Garment Employees (1197, 14)
  Garment = fetch_safe(id = 597),
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # small n
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Auto MPG (398, 7)
  Auto = fetch_safe(id = 9),
  # Automobile (205, 25)
  Automobile = fetch_safe(id = 10)
)

# Remove NULLs (datasets that failed to download) and report which ones failed
failed <- names(datasets)[sapply(datasets, is.null)]
if (length(failed) > 0) {
  warning("Failed to download: ", paste(failed, collapse = ", "), immediate. = TRUE)
}
datasets <- datasets[!sapply(datasets, is.null)]
cat("Successfully downloaded", length(datasets), "datasets\n")

# Crime: drop the 5 non-predictive attributes documented by UCI
# (state, county, community, communityname and fold), keeping the 122 predictors
# and the goal (ViolentCrimesPerPop).
if ("Crime" %in% names(datasets)) {
  nm <- colnames(datasets$Crime$data)
  nm <- nm[! nm %in% c("state", "county", "community", "communityname", "fold")]
  datasets$Crime$data <- datasets$Crime$data[, nm]
}

# Add R datasets
cat("Adding R datasets...\n")

# Add Ozone --- small (587 13)
# V4 (daily maximum one-hour-average ozone) is the response. V1-V3 (month,
# day of month and day of week) are factors and are dropped by prepare_data().
data(Ozone, package = "mlbench")
datasets$Ozone <- list(
  data = as.data.frame(Ozone),
  response = "V4",
  name = "Ozone",
  n = nrow(Ozone),
  p = ncol(Ozone) - 1
)

# Prepare each dataset
cat("\nPreparing datasets for analysis...\n")
datasets <- purrr::map(datasets, prepare_data)

# Create directory for individual datasets
if (!dir.exists("data/individual")) dir.create("data/individual", recursive = TRUE)

# Save each dataset separately
for (name in names(datasets)) {
  saveRDS(datasets[[name]], paste0("data/individual/", name, ".rds"))
  cat("Saved:", name, "\n")
}

cat("All datasets saved individually in 'data/individual/'\n")

# Quick summary
summary_df <- purrr::map_df(names(datasets), function(name) {
  d <- datasets[[name]]
  data.frame(
    Dataset = name,
    Response = d$response,
    n = d$n,
    p = d$p,
    has_missing = d$has_natural_missing,
    missing_pct = round(100 * sum(is.na(d$X)) / (d$n * d$p), 2),
    stringsAsFactors = FALSE
  )
}) %>% dplyr::arrange(dplyr::desc(n))

print(summary_df)
