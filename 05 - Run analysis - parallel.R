res_dir <- "results"
rng_kind <- "default"

#res_dir <- "resultsRNG2"
#rng_kind <- "L'Ecuyer-CMRG"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN PRTree vs CART EXPERIMENTS - full experiment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PACKAGES AND FILES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(parallel)
library(pbapply)
library(dplyr)
library(purrr)

# Load required functions
source("03 - functions.R")

# Set the directories for individual data and results 
wd <- getwd()
data_dir <- "data/individual"

if (!dir.exists(res_dir)) dir.create(res_dir)
individual_files <- list.files(res_dir, pattern = "\\.rds$")
processed_datasets <- gsub("\\.rds$", "", individual_files)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMMON CONFIGURATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seed = 1234            # (for CART and PRTree)
set.seed(seed + 1000)  # global experiment

# List of configurations (need update for each dataset)
config <- PRTree::pr_tree_control_cv(
  method = "kfold",     # k-fold cross-validation
  n_rep = 15,                # 10 folds
  stratify = FALSE,          # Stratify by missing pattern
  max_terminal_nodes = 100,  # Maximum terminal nodes (large to allow full growth)
  cp = 0.001,                # Complexity parameter
  by_node = TRUE,
  tiny_sigma = 0,            # Near-zero sigma candidate
  #prop_valid = 0,            # 10% for validation (from training)
  verbose = FALSE,            # Suppress messages
  update_final = FALSE
)

# add CART control parameters
config$maxdepth <- config$max_terminal_nodes - 1
config$xval <- 10

prop_min <- 0.01
dist_grid <- tribble(
  ~dist,    ~df, ~shape, ~sdlog,
  # Gaussian baseline
  "norm",   NA,  NA,     NA,
  # Student-t
  # df controls tail heaviness
  # lower df => heavier tails
  "t",      3,   NA,     NA,
  "t",      5,   NA,     NA,
  # Lognormal
  # sdlog controls skewness/spread
  "lnorm",  NA,  NA,     0.5,
  "lnorm",  NA,  NA,     1.0,
  "lnorm",  NA,  NA,     2.0,
  # Gamma
  # shape controls asymmetry
  "gamma",  NA,  0.5,    NA,
  #"gamma",  NA,  1.0,    NA,
  "gamma",  NA,  2.0,    NA,
  "gamma",  NA,  5.0,    NA
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATASETS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# List of datasets sorted by n (largest to smallest)
dataset_names <- c(
  ##------ no missing -------
  # "Super",        # 21,263 rows, 81 features ← SEQUENCIAL
  # "Bike",         # 17,379 rows, 15 features
  # "Power",        # 9,568 rows, 4 features
  # "AirQuality",   # 9,357 rows, 12 features
  # "Abalone",      # 4,177 rows, 7 features
  # "Airfoil",      # 1,503 rows, 5 features
  # "Concrete",     # 1,030 rows, 8 features
  # "Energy",       # 768 rows, 9 features
  # "Boston",       # 506 rows, 13 features
  # "Hitters",      # 263 rows, 16 features
  ##------ with natural missing -------
  "Automobile",   # 205 rows, 17 features
  "Ozone",        # 361 rows, 9 features  
  "Auto",          # 398 rows, 7 features  
  "Garment",      # 1,197 rows, 14 features  
  "Crime"         # 1,994 rows, 122 features
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECKPOINT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datasets_to_process <- setdiff(dataset_names, processed_datasets)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REPOR THE CURRENT STATUS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n", strrep("=", 70), "\n")
cat("RUNNING EXPERIMENTS\n")
cat(strrep("=", 70), "\n")

cat("Datasets total:", length(dataset_names), "\n")
cat("Already processed:", length(processed_datasets), "\n")
cat("Remaining:", length(datasets_to_process), "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WORKER FUNCTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_single_dataset <- function(name, config, data_dir, prop_min, seed,
                               dist_grid, res_dir, wd) {
  
  log_file <- file.path(res_dir, paste0(name, ".log"))
  zz <- file(log_file, open = "wt")
  sink(zz)
  sink(zz, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(zz)
  }, add = TRUE)
  
  tryCatch({
    cat("\n", strrep("=", 70), "\n")
    cat("DATASET:", name, "\n")
    cat("START :", as.character(Sys.time()), "\n")
    cat(strrep("=", 70), "\n\n")
    
    data <- readRDS(file.path(data_dir, paste0(name, ".rds")))
    cfg <- config
    n_train <- (cfg$n_rep - 1) * nrow(data$X) / cfg$n_rep
    cfg$prop_valid <- 1 / (cfg$n_rep - 1)
    cfg$n_min <- nmin(prop = prop_min, n_train = n_train)
    cfg$minsplit <- cfg$n_min
    cfg$n_candidates <- as.integer(max(round(ncol(data$X)/2), 3))
    
    cat("n            =", nrow(data$X), "\n")
    cat("p            =", ncol(data$X), "\n")
    cat("minsplit     =", cfg$minsplit, "\n")
    cat("n_candidates =", cfg$n_candidates, "\n\n")
    
    cat("[", Sys.time(), "] Generating folds...\n")
    set.seed(seed)
    cfg$fold_idx <- get_folds(config = cfg, seed = seed, X = data$X)
    
    cat("[", Sys.time(), "] Running experiment...\n")
    result <- run_dataset_experiment(
      data = data,
      config = cfg,
      verbose = TRUE,
      seed = seed,
      dist_grid = dist_grid
    )
    
    cat("[", Sys.time(), "] Saving results...\n")
    saveRDS(result, file.path(res_dir, paste0(name, ".rds")))
    
    rm(result)
    gc(verbose = FALSE)
    
    cat("\n", strrep("=", 70), "\n")
    cat("FINISHED:", name, "\n")
    cat("END     :", as.character(Sys.time()), "\n")
    cat(strrep("=", 70), "\n")
    
    return(name)
    
  }, error = function(e) {
    
    cat("\n", strrep("=", 70), "\n")
    cat("ERROR IN:", name, "\n")
    cat("TIME    :", as.character(Sys.time()), "\n")
    cat("MESSAGE :", e$message, "\n")
    cat(strrep("=", 70), "\n")
    
    return(NULL)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PARALLEL EXECUTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (length(datasets_to_process) > 0) { 
  # Print the current status
  cat("\n", strrep("=", 70), "\n")
  cat("RUNNING IN PARALLEL\n")
  cat(strrep("=", 70), "\n")
  
  # Strictly force 4 workers
  n_cores <- min(4, length(datasets_to_process))
  cat("Using", n_cores, "cores\n")
  
  # Create a native Windows process cluster (bypasses local network/firewall socket limitations)
  cl <- makeCluster(n_cores)
  
  # Explicitly export necessary global variables into the isolated node spaces
  clusterExport(
    cl,
    varlist = c(
      "datasets_to_process",
      "config",
      "data_dir",
      "prop_min",
      "seed",
      "dist_grid",
      "res_dir",
      "run_single_dataset",
      "wd",
      "rng_kind"
    )
  )
  
  # Initialize background workers
  clusterEvalQ(cl, {
    RNGkind(rng_kind)
    library(parallel)
    library(PRTree)
    library(dplyr)
    library(purrr)
    source("03 - functions.R")
    setwd(wd)
    NULL
  })
  # Run the native sifting array loop
  
  pboptions(type = "timer")
  
  out <- parallel::parLapplyLB(
    cl,
    datasets_to_process,
    run_single_dataset,
    config = config,
    data_dir = data_dir,
    prop_min = prop_min,
    seed = seed,
    dist_grid = dist_grid,
    res_dir = res_dir,
    wd = wd
  )
  
  # Securely close server channels and stop hardware threads
  stopCluster(cl)
  cat("\nParallel execution complete.\n")
} 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n", strrep("=", 70), "\n")
cat("LOADING RESULTS\n")
cat(strrep("=", 70), "\n")

result_files <- file.path(res_dir, paste0(dataset_names, ".rds"))
existing_files <- result_files[file.exists(result_files)]

if(length(existing_files) > 0) {
  all_results_raw <- list_rbind(lapply(existing_files, readRDS))
  saveRDS(all_results_raw, paste0(res_dir, "/prtree_cart_raw.rds"))
  write.csv(all_results_raw, paste0(res_dir, "/prtree_cart_raw.csv"), row.names = FALSE)
  cat("\nFinal results saved.\n")
} else {
  cat("\nNo output data targets located for aggregation.\n")
}

cat("All processing complete.\n")
