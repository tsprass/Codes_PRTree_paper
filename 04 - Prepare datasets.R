# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT 4: DOWNLOAD AND PREPARE SELECTED DATASETS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script downloads the datasets selected in 02-Filter.R and prepares them
# Input:  regression_metadata from 02-Filter.R
# Output: data/processed_datasets_ready.rds

# Load functions
source("03 - functions.R")

# create the directory to save the data
if (!dir.exists("data")) dir.create("data")

# Load filtered dataset list (you need to save this from 02-Filter.R)
# For now, we'll define them manually based on what you selected
cat("🌐 Downloading selected UCI datasets...\n")

# # Fetch datasets
# datasets <- list(
#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # very large n 
#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Superconductivty Data (21263  81)
#   Super = fetch_safe(id = 464),
#   # Bike Sharing (17389 13)
#   Bike = fetch_safe(id = 275),
#   
#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # large n 
#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Combined Cycle Power Plant (9568 4)
#   Power = fetch_safe(id = 294),
#   # Air Quality (9358 15)
#   AirQuality = fetch_safe(id = 360),
#   
#   # ~~~~~~~~~~~~~~~~~~~~~~~~
#   # mediun n
#   # ~~~~~~~~~~~~~~~~~~~~~~~~
#   # Abalone (4177 8)
#   Abalone = fetch_safe(id = 1),           
#   # Airfoil Self-Noise (1503 5)
#   Airfoil = fetch_safe(id = 291),  
#   # Concrete Compressive Strength (1030 8)
#   Concrete = fetch_safe(id = 165),
#   
#   # ~~~~~~~~~~~~~~~~~~~~~~~~
#   # small n
#   # ~~~~~~~~~~~~~~~~~~~~~~~~  
#   # Energy Efficiency (768 8)
#   Energy = fetch_safe(id = 242),
#   # Auto MPG (398 7)
#   Auto = fetch_safe(id = 9)
# )

# Fetch datasets with missing values
datasets <- list(
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # large n
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Communities and Crime (1994, 127), missing tem "?"
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
  Automobile = fetch_safe(id = 10),
  # Bone marrow transplant: children (187, 36)
  BoneMarrow = fetch_safe(id = 565),
  # Tennis Major Tournament Match Statistics (127, 42)
  Tennis = fetch_safe(id = 300)
)

# Remove NULLs
datasets <- datasets[!sapply(datasets, is.null)]
cat("✅ Successfully downloaded", length(datasets), "datasets\n")

nm <- colnames(datasets$Crime$data)
nm <- nm[! nm %in% c("state", "county", "community", "communityname", "fold")]
datasets$Crime$data <- datasets$Crime$data[, nm]

# Add R datasets
cat("📦 Adding R datasets...\n")

# Add Ozone --- small (587 13)
data(Ozone, package = "mlbench")
datasets$Ozone <- list(
  data = as.data.frame(Ozone),
  response = "V4",
  name = "Ozone",
  n = nrow(Ozone),
  p = ncol(Ozone) - 1
)

# # Add Hitters --- small (322 20)
# data(Hitters, package = "ISLR")
# datasets$Hitters <- list(
#   data = Hitters,
#   response = "Salary",
#   name = "Hitters",
#   n = nrow(Hitters),
#   p = ncol(Hitters) - 1
# )
# 
# # Add Boston --- small (506  14)
# data(Boston, package = "MASS")
# datasets$Boston <- list(
#   data = Boston,
#   response = "medv",
#   name = "Boston Housing",
#   n = nrow(Boston),
#   p = ncol(Boston) - 1
# )

# Prepare each dataset
cat("\n🔧 Preparing datasets for analysis...\n")
datasets <- purrr::map(datasets, prepare_data)

# Create directory for individual datasets
if (!dir.exists("data/individual")) dir.create("data/individual", recursive = TRUE)

# Save each dataset separately
for (name in names(datasets)) {
  saveRDS(datasets[[name]], paste0("data/individual/", name, ".rds"))
  cat("✅ Saved:", name, "\n")
}

cat("🎉 All datasets saved individually in 'data/individual/'\n")

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
