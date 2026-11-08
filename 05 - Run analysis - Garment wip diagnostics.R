# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT 5b: GARMENT WIP DIAGNOSTICS (DIAGNOSTIC COMPLEMENT TO SCRIPT 5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the Garment dataset, 'wip' (work in progress) is the only predictor with
# missing values (506 missing out of 1197 rows). Script 5 shows that PRTree
# did not outperform CART on Garment.
#
# This script investigates why by refitting PRTree and CART across 6 conditions
# under the exact same 15 cross-validation folds, seeds, and tree-growing settings:
#   1. without wip, defaults (baseline without the missing feature)
#   2. with wip, defaults
#   3. with wip, prop_x = p_min = 0 (removes minimum probability split constraint)
#   4. with wip, prop_x = p_min = 0, sigma = 0 (disables smoothing completely)
#   5. with wip, prop_x = p_min = 0, grid up to 4 sd (wider sigma candidate range)
#   6. with wip, prop_x = p_min = 0, sigma = 0 for discrete predictors (<= 12 values)
#
# Output (in results_garment_wip/):
#   - Garment_wip_diagnostics.rds : Full list with CART and PRTree per-fold RMSE and sigma
#   - Garment_wip_diagnostics.csv : Per-fold RMSE comparisons in tabular format
#   - Summary metrics printed to console
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. PACKAGES AND SETTINGS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(parallel)
library(PRTree)
library(dplyr)
library(purrr)

source("03 - functions.R")

seed    <- 1234
n_cores <- 4
res_dir <- "results_garment_wip"
out_file <- file.path(res_dir, "Garment_wip_diagnostics.rds")
csv_file  <- file.path(res_dir, "Garment_wip_diagnostics.csv")

if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. BASE CONFIGURATION (MATCHING SCRIPT 5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
config <- PRTree::pr_tree_control_cv(
  method             = "kfold",
  n_rep              = 15,
  stratify           = FALSE,
  max_terminal_nodes = 100,
  cp                 = 0.001,
  by_node            = TRUE,
  tiny_sigma         = 0,
  verbose            = FALSE,
  update_final       = FALSE
)

# CART parameters
config$maxdepth <- config$max_terminal_nodes - 1
config$xval     <- 10

prop_min <- 0.01

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. PREPARE DATASETS: WITH WIP AND WITHOUT WIP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dataset 1: Full Garment (includes wip)
data_wip      <- readRDS("data/individual/Garment.rds")
data_wip$name <- "Garment"
data_wip$n    <- nrow(data_wip$X)
data_wip$p    <- ncol(data_wip$X)

# Dataset 2: Garment without wip
data_now      <- data_wip
data_now$X    <- data_wip$X[, colnames(data_wip$X) != "wip", drop = FALSE]
data_now$name <- "Garment without wip"
data_now$p    <- ncol(data_now$X)

stopifnot(!anyNA(data_now$X))

# Configure cross-validation parameters
cfg            <- config
n_train        <- (cfg$n_rep - 1) * data_wip$n / cfg$n_rep
cfg$prop_valid <- 1 / (cfg$n_rep - 1)
cfg$n_min      <- nmin(prop = prop_min, n_train = n_train)
cfg$minsplit   <- cfg$n_min
cfg$n_candidates <- as.integer(max(round(data_wip$p / 2), 3))

# Cross-validation folds (shared across all conditions and CART)
set.seed(seed)
cfg$fold_idx   <- get_folds(config = cfg, seed = seed, X = data_wip$X)

# PRTree standard options
cfg$fill_type  <- 2
cfg$proxy_crit <- "both"
cfg$dist       <- "norm"
cfg$grid_size  <- 8

# Custom sigma grid: smooth only continuous predictors (discrete <= 12 unique values get sigma = 0)
discrete <- sapply(as.data.frame(data_wip$X), function(x) length(unique(x[!is.na(x)]))) <= 12
grid_cont <- unclass(PRTree::expand_sigma_grid(
  X          = data_wip$X,
  grid_size  = 8,
  tiny_sigma = 0,
  verbose    = FALSE
))
attr(grid_cont, "has_tiny") <- NULL
grid_cont[, discrete]       <- 0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. EXPERIMENTAL CONDITIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
conditions <- list(
  list(
    label  = "without wip, defaults",
    wip    = FALSE,
    prop_x = 0.1
  ),
  list(
    label  = "with wip, defaults",
    wip    = TRUE,
    prop_x = 0.1
  ),
  list(
    label  = "with wip, prop_x = p_min = 0",
    wip    = TRUE,
    prop_x = 0
  ),
  list(
    label  = "with wip, prop_x = p_min = 0, sigma = 0",
    wip    = TRUE,
    prop_x = 0,
    args   = list(grid_size = 0)
  ),
  list(
    label  = "with wip, prop_x = p_min = 0, grid up to 4 sd",
    wip    = TRUE,
    prop_x = 0,
    args   = list(max_mult = 4)
  ),
  list(
    label  = "with wip, prop_x = p_min = 0, sigma = 0 for discrete predictors",
    wip    = TRUE,
    prop_x = 0,
    args   = list(sigma_grid = grid_cont)
  )
)

# Worker function for one condition
run_condition <- function(cond, data_wip, data_now, cfg, seed) {
  data   <- if (cond$wip) data_wip else data_now
  config <- cfg

  config$prop_x <- cond$prop_x
  config$p_min  <- if (cond$prop_x == 0) 0 else 0.05

  for (a in names(cond$args)) {
    config[[a]] <- cond$args[[a]]
  }

  set.seed(seed)
  cv <- PRTree::pr_tree_cv(data$y, data$X, control_cv = config, verbose = FALSE)

  list(
    label = cond$label,
    wip   = cond$wip,
    rmse  = cv$rmse_by_rep$test,
    sigma = cv$sigma_matrix
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. EXECUTION (PARALLEL OVER CONDITIONS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (file.exists(out_file)) {
  cat("\nResults found in", out_file, "- skipping execution.\n")
  out <- readRDS(out_file)
} else {
  n_cl <- min(n_cores, length(conditions))
  cat("\nRunning CART and", length(conditions), "PRTree conditions across", n_cl, "cores...\n")
  t0 <- Sys.time()

  # Baseline CART (with and without wip)
  cart <- list(
    wip   = run_cart_cv(data_wip, config = cfg, seed = seed),
    nowip = run_cart_cv(data_now, config = cfg, seed = seed)
  )

  # Parallel execution of PRTree conditions
  cl <- makeCluster(n_cl)
  wd <- getwd()
  clusterExport(cl, "wd")
  invisible(clusterEvalQ(cl, {
    setwd(wd)
    library(PRTree)
    library(dplyr)
    source("03 - functions.R")
    NULL
  }))

  fits <- parLapplyLB(
    cl,
    conditions,
    run_condition,
    data_wip = data_wip,
    data_now = data_now,
    cfg      = cfg,
    seed     = seed
  )
  stopCluster(cl)

  out <- list(cart = cart, fits = fits)
  saveRDS(out, out_file)
  cat("Execution finished in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes.\n")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. COMPARISONS AND EXPORT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build summary table across all folds:
#   - diff   : mean(RMSE_PRTree - RMSE_CART)
#   - se     : standard error of the mean difference
#   - Delta  : median((RMSE_PRTree - RMSE_CART) / RMSE_CART)
#   - win    : proportion of folds where PRTree beat CART (diff < 0)
#   - sigma0 : number of folds where all selected sigmas were zero
tab <- bind_rows(lapply(out$fits, function(f) {
  cart_rmse <- if (f$wip) out$cart$wip else out$cart$nowip
  d         <- f$rmse - cart_rmse
  
  data.frame(
    condition    = f$label,
    rmse         = mean(f$rmse),
    rmse_cart    = mean(cart_rmse),
    diff         = mean(d),
    se           = sd(d) / sqrt(length(d)),
    Delta_median = median(d / cart_rmse),
    win          = mean(d < 0),
    sigma0       = sum(rowSums(f$sigma) == 0),
    stringsAsFactors = FALSE
  )
}))

# Save detailed per-fold CSV table
fold_df <- bind_rows(lapply(out$fits, function(f) {
  cart_rmse <- if (f$wip) out$cart$wip else out$cart$nowip
  data.frame(
    condition = f$label,
    fold      = seq_along(f$rmse),
    rmse      = f$rmse,
    rmse_cart = cart_rmse,
    stringsAsFactors = FALSE
  )
}))
write.csv(fold_df, csv_file, row.names = FALSE)

# Print console summary
cat("\n", strrep("=", 80), "\n")
cat("GARMENT WIP DIAGNOSTICS: SUMMARY RESULTS\n")
cat(strrep("=", 80), "\n\n")

cat("CART (mean RMSE over 15 folds):\n")
cat("  - With wip    :", round(mean(out$cart$wip), 4), "(sd =", round(sd(out$cart$wip), 4), ")\n")
cat("  - Without wip :", round(mean(out$cart$nowip), 4), "(sd =", round(sd(out$cart$nowip), 4), ")\n\n")

print(tab, digits = 3, row.names = FALSE)
cat("\nResults saved to:\n  -", out_file, "\n  -", csv_file, "\n")
