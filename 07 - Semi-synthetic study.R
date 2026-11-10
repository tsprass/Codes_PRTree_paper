# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT 7: MASTER SEMI-SYNTHETIC SIMULATION STUDY
#
# Evaluates PRTree with internal strategies (fill_type 0, 1, 2) versus
# impute-then-PRTree and CART across benchmark datasets with complete covariates
# (Auto MPG, Ozone, Concrete) under three distinct missingness scenarios:
#
#   Scenario 1: Main Paper - Section 6
#     - Missingness: One-sided (upper-tail) censorship in the strongest predictors
#     - Target set: ceil(p/2) predictors most correlated with y (|cor(X_j, y)|)
#     - MNAR: highest values removed (sensor saturation thresholds)
#     - MAR : rows with highest guide values flagged; strong predictors deleted with pi0 = 0.5
#     - MCAR: random cells in strong predictors deleted independently
#     - Replications: R = 100 (5-fold cross-validation)
#     - Output: results_semisynthetic_upper/
#
#   Scenario 2: Exploratory Sensitivity - Strong Predictors
#     - Missingness: Two-sided censorship in the strongest predictors
#     - MNAR: lowest and highest values removed (half each)
#     - MAR : extreme values (|Z|) in guides trigger deletion
#     - Replications: R = 10 (5-fold cross-validation)
#     - Output: results_semisynthetic_strong/
#
#   Scenario 3: Exploratory Baseline - All Predictors
#     - Missingness: Induced indiscriminately across all predictors (p variables)
#     - MNAR: lowest and highest values removed across all columns
#     - MAR : random partition of predictors into guides and targets
#     - Replications: R = 10 (5-fold cross-validation)
#     - Output: results_semisynthetic_mech/
#
# Output artifacts per scenario:
#   <dataset>.rds              raw task results (checkpointed; skipped if already done)
#   semisynthetic_boxplots.pdf Figure 4 in paper (relative RMSE boxplots)
#   tab_semisynthetic.tex      Table 3 in paper (median & IQR of relative RMSE)
#   semisynthetic_heatmaps.pdf Figure 5 in paper (pairwise win rate heatmaps of PR-f2)
#   tab_semisynthetic_time.tex Table 4 in paper (median runtime & ratio to PR-f2)
#   semisynthetic_time.pdf     Figure 6 in paper (median runtime per fold heatmaps)
#
# Run from the folder "Applic".
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("07 - functions_semisynthetic.R")

# ------------------------------------------------------------------------------
# User Configuration
# ------------------------------------------------------------------------------
n_cores     <- 4       # Number of CPU cores for parallel execution
clean_start <- FALSE   # If TRUE, removes previous .rds files and restarts from scratch.
                       # If FALSE, uses fault-tolerant checkpoints (skips finished tasks).

# Select which scenario(s) to execute:
#   1 = Upper-tail in strong predictors (Main Paper - Section 6, R = 100)
#   2 = Two-tail in strong predictors   (Exploratory sensitivity, R = 10)
#   3 = Two-tail across all predictors  (Exploratory baseline,    R = 10)
# Examples: 1:3 (runs all), 1 (runs only main paper), c(2, 3) (runs exploratory)
run_scenarios <- 1:3

# ------------------------------------------------------------------------------
# Scenarios Definition
# ------------------------------------------------------------------------------
scenarios <- list(
  list(
    id        = 1,
    name      = "1. Upper-tail in strong predictors (Main Paper - Section 6)",
    res_dir   = "results_semisynthetic_upper",
    version   = "rev1-upper",
    inject_fn = inject_upper_tail_strong,
    n_reps    = 100
  ),
  list(
    id        = 2,
    name      = "2. Two-tail in strong predictors (Exploratory sensitivity)",
    res_dir   = "results_semisynthetic_strong",
    version   = "rev1-strong",
    inject_fn = inject_two_tail_strong,
    n_reps    = 10
  ),
  list(
    id        = 3,
    name      = "3. Two-tail across all predictors (Exploratory baseline)",
    res_dir   = "results_semisynthetic_mech",
    version   = "rev1-mech",
    inject_fn = inject_two_tail_all,
    n_reps    = 10
  )
)

# Filter scenarios according to user selection
active_scenarios <- Filter(function(s) s$id %in% run_scenarios, scenarios)

total_t0 <- Sys.time()

# ------------------------------------------------------------------------------
# Execution Loop
# ------------------------------------------------------------------------------
for (sc in active_scenarios) {
  cat("\n", strrep("=", 75), "\n", sep = "")
  cat(sprintf("STARTING SCENARIO: %s\n", sc$name))
  cat(sprintf("Output directory: %s | Replications: %d | Cores: %d\n",
              sc$res_dir, sc$n_reps, n_cores))
  cat(strrep("=", 75), "\n\n", sep = "")

  if (clean_start && dir.exists(sc$res_dir)) {
    old_files <- list.files(sc$res_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(old_files) > 0) {
      cat(sprintf("Clean start enabled: removing %d previous .rds files from %s...\n",
                  length(old_files), sc$res_dir))
      file.remove(old_files)
    }
  }

  design <- list(
    datasets   = c("Auto", "Ozone", "Concrete"),
    mechanisms = c("MCAR", "MAR", "MNAR"),
    rates      = c(0.10, 0.25, 0.40),
    R          = sc$n_reps,
    K          = 5,
    pi0        = 0.5,
    k_knn      = 5,
    m_mice     = 5,
    prop_min   = 0.01,
    dist       = "norm",
    proxy      = "both",
    seed       = 2026
  )

  t0_sc <- Sys.time()

  # 1. Run parallel simulation with checkpointing
  run_semisynthetic_study(
    design    = design,
    inject_fn = sc$inject_fn,
    res_dir   = sc$res_dir,
    version   = sc$version,
    n_cores   = n_cores
  )

  # 2. Generate LaTeX tables, PDF figures and terminal summaries
  generate_semisynthetic_summaries(
    design  = design,
    res_dir = sc$res_dir,
    version = sc$version
  )

  cat(sprintf("\n[COMPLETED] %s in %.1f minutes.\n",
              sc$name, as.numeric(difftime(Sys.time(), t0_sc, units = "mins"))))
}

cat("\n", strrep("=", 75), "\n", sep = "")
cat(sprintf("ALL REQUESTED SCENARIOS FINISHED in %.1f minutes!\n",
            as.numeric(difftime(Sys.time(), total_t0, units = "mins"))))
cat(strrep("=", 75), "\n", sep = "")
