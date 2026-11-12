# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SHARED FUNCTIONS FOR SEMI-SYNTHETIC SIMULATION STUDIES
#
# Contains data loaders, missingness injection mechanisms, imputation pipelines,
# PRTree and CART model fitting helpers, parallel task execution with checkpointing,
# and LaTeX table and figure generation.
# Sourced by:
#   - 07 - Semi-synthetic study.R (Main Paper - Section 6)
#   - 07 - Semi-synthetic study - exploratory all predictors.R
#   - 07 - Semi-synthetic study - exploratory strong predictors.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

# Loads the specified benchmark dataset and restricts to complete cases.
# This ensures that synthetic missingness is injected into a clean, 100% complete
# ground-truth baseline, allowing exact computation of relative loss and benchmarks.
load_complete <- function(name) {
  # Search both the standard individual data directory and complete-cases directory
  f1 <- file.path("data/individual", paste0(name, ".rds"))
  f2 <- file.path("data/individual_complete", paste0(name, ".rds"))
  d <- readRDS(if (file.exists(f1)) f1 else f2)

  # Filter complete cases only (no pre-existing missing values in covariates)
  cc <- stats::complete.cases(d$X)
  X <- d$X[cc, , drop = FALSE]

  # Ensure column names are syntactically valid and unique for formula interfaces (y ~ .)
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  list(name = name, X = X, y = d$y[cc])
}

# ==============================================================================
# 2. MISSINGNESS INJECTION MECHANISMS
# ==============================================================================

# Helper: identify the ceil(p/2) predictors most correlated with the response y.
# Ranking is based on absolute Pearson correlation |cor(X_j, y)|.
# Note: X is strictly complete at this stage, as predictors are selected prior to
# missingness injection. Setting cy[is.na(cy)] <- 0 provides a safety guard
# against constant / zero-variance columns where sd(X_j) == 0 (division by zero).
get_strong_predictors <- function(X, y) {
  cy <- abs(suppressWarnings(stats::cor(X, y)))[, 1]
  cy[is.na(cy)] <- 0  # Zero-variance columns receive zero correlation
  sort(order(cy, decreasing = TRUE)[seq_len(ceiling(ncol(X) / 2))])
}

# Mechanism A (Main Paper - Section 6): Upper-tail deletion in strongest predictors
#
# Missing values are injected ONLY into the set A of ceil(p/2) strongest predictors.
# The remaining floor(p/2) predictors (set O) remain 100% fully observed.
#
# - MCAR: Each cell in A is deleted independently with probability q.
#         By probability theory, the probability of a row remaining complete in A is
#         (1 - q)^na. Setting (1 - q)^na = 1 - r yields q = 1 - (1 - r)^(1 / na),
#         guaranteeing that on average exactly a fraction r of rows will be incomplete.
#
# - MAR : Deletion in strong predictors depends only on observed predictors (set O).
#         Predictors in O are standardized (Z = scale(X)). For each row, the guide
#         signal is s_i = max_{k in O}(Z_ik), measuring how high the observed guides are.
#         The r*n rows with the highest guide values are flagged. In each flagged row,
#         each strong predictor is deleted with probability pi0 = 0.5 (repeated until
#         at least one predictor is deleted to ensure the row becomes incomplete).
#
# - MNAR (Upper-tail): Models sensor saturation / upper-detection limits. In each
#         strong predictor, the highest values are removed. A bisection search (40 steps)
#         calibrates the exact cutoff/fraction q so the proportion of incomplete rows matches r.
inject_upper_tail_strong <- function(X, y, mech, r, pi0) {
  n <- nrow(X)
  p <- ncol(X)
  A <- get_strong_predictors(X, y)   # Indices of strongest predictors (target set)
  O <- setdiff(seq_len(p), A)        # Indices of remaining predictors (guide set)
  na <- length(A)
  Z <- scale(X)
  M <- matrix(FALSE, n, p)           # Missingness mask (TRUE = deleted)

  if (mech == "MCAR") {
    # Analytic cell deletion probability yielding target incomplete row rate r
    q <- 1 - (1 - r)^(1 / na)
    M[, A] <- matrix(runif(n * na) < q, n, na)

  } else if (mech == "MAR") {
    # Highest values in the observed guides trigger missingness in target set A
    s <- apply(Z[, O, drop = FALSE], 1, max)
    flag <- rank(-s, ties.method = "random") <= round(r * n)  # Top r*n extreme rows
    for (i in which(flag)) {
      repeat {
        d <- runif(na) < pi0
        if (any(d)) break
      }
      M[i, A[d]] <- TRUE
    }

  } else {
    # Upper-tail MNAR: delete the k highest values in each predictor in A
    rk <- apply(X[, A, drop = FALSE], 2, rank, ties.method = "random")
    cells <- function(q) {
      tot <- round(q * n * na)                                # Total cells to delete
      k <- rep(tot %/% na, na) + (seq_len(na) <= tot %% na)   # Cells per predictor
      sweep(rk, 2, n - k, ">")                                # Delete k highest values
    }
    inc <- function(q) mean(rowSums(cells(q)) > 0)

    # Bisection search to calibrate cell proportion q for target row rate r
    lo <- 0
    hi <- 1
    for (it in 1:40) {
      mid <- (lo + hi) / 2
      if (inc(mid) < r) lo <- mid else hi <- mid
    }
    M[, A] <- cells(hi)
  }

  Xm <- X
  Xm[M] <- NA
  Xm
}

# Mechanism B (Exploratory): Two-tail deletion in strongest predictors
#
# Same target set A (strongest predictors) and guide set O as Mechanism A, but:
# - MAR : The guide signal is s_i = max_{k in O} |Z_ik| (extremeness in either tail).
# - MNAR: Both lowest and highest values are removed (two-sided censorship, half each).
inject_two_tail_strong <- function(X, y, mech, r, pi0) {
  n <- nrow(X)
  p <- ncol(X)
  A <- get_strong_predictors(X, y)
  O <- setdiff(seq_len(p), A)
  na <- length(A)
  Z <- scale(X)
  M <- matrix(FALSE, n, p)

  if (mech == "MCAR") {
    q <- 1 - (1 - r)^(1 / na)
    M[, A] <- matrix(runif(n * na) < q, n, na)

  } else if (mech == "MAR") {
    # Extremeness in either tail (|Z_ik|) of observed guides flags rows
    s <- apply(abs(Z[, O, drop = FALSE]), 1, max)
    flag <- rank(-s, ties.method = "random") <= round(r * n)
    for (i in which(flag)) {
      repeat {
        d <- runif(na) < pi0
        if (any(d)) break
      }
      M[i, A[d]] <- TRUE
    }

  } else {
    # Two-tail MNAR: delete lowest and highest values in each strong predictor
    rk <- apply(X[, A, drop = FALSE], 2, rank, ties.method = "random")
    cells <- function(q) {
      tot <- round(q * n * na)
      k <- rep(tot %/% na, na) + (seq_len(na) <= tot %% na)
      lo <- floor(k / 2)
      hi <- k - lo
      sweep(rk, 2, lo, "<=") | sweep(rk, 2, n - hi, ">")       # Half lowest, half highest
    }
    inc <- function(q) mean(rowSums(cells(q)) > 0)
    lo <- 0
    hi <- 1
    for (it in 1:40) {
      mid <- (lo + hi) / 2
      if (inc(mid) < r) lo <- mid else hi <- mid
    }
    lo <- 0
    hi <- 1
    for (it in 1:40) {
      mid <- (lo + hi) / 2
      if (inc(mid) < r) lo <- mid else hi <- mid
    }
    M[, A] <- cells(hi)
  }

  Xm <- X
  Xm[M] <- NA
  Xm
}

# Mechanism C (Exploratory): Two-tail deletion across all predictors
#
# Missingness is induced indiscriminately across ALL predictors (p variables):
# - MCAR: Random cells deleted across all p columns with q = 1 - (1 - r)^(1 / p).
# - MAR : A random subset of floor(p/2) predictors is sampled to serve as guides O.
#         The remaining predictors form the target set A.
# - MNAR: Two-tail deletion (lowest and highest values, half each) across all p columns.
inject_two_tail_all <- function(X, y, mech, r, pi0) {
  n <- nrow(X)
  p <- ncol(X)
  Z <- scale(X)

  if (mech == "MCAR") {
    q <- 1 - (1 - r)^(1 / p)
    M <- matrix(runif(n * p) < q, n, p)

  } else if (mech == "MAR") {
    # Randomly partition predictors into guides (O) and targets (A)
    O <- sort(sample.int(p, floor(p / 2)))
    A <- setdiff(seq_len(p), O)
    s <- apply(abs(Z[, O, drop = FALSE]), 1, max)
    flag <- rank(-s, ties.method = "random") <= round(r * n)
    M <- matrix(FALSE, n, p)
    for (i in which(flag)) {
      repeat {
        d <- runif(length(A)) < pi0
        if (any(d)) break
      }
      M[i, A[d]] <- TRUE
    }

  } else {
    # Two-tail MNAR across all predictors
    rk <- apply(X, 2, rank, ties.method = "random")
    cells <- function(q) {
      tot <- round(q * n * p)
      k <- rep(tot %/% p, p) + (seq_len(p) <= tot %% p)
      lo <- floor(k / 2)
      hi <- k - lo
      sweep(rk, 2, lo, "<=") | sweep(rk, 2, n - hi, ">")
    }
    inc <- function(q) mean(rowSums(cells(q)) > 0)
    lo <- 0
    hi <- 1
    for (it in 1:40) {
      mid <- (lo + hi) / 2
      if (inc(mid) < r) lo <- mid else hi <- mid
    }
    M <- cells(hi)
  }

  Xm <- X
  Xm[M] <- NA
  Xm
}

# ==============================================================================
# 3. IMPUTATION METHODS
#
# Crucial methodological rules:
#   1. Training fold only: Imputation parameters (means, donors, random forests,
#      MICE chained models) are estimated strictly on Xtr.
#   2. The test fold Xte is imputed using the models/statistics learned on Xtr,
#      strictly preventing data leakage.
#   3. The response y is NEVER used in the imputation pipelines, matching
#      unsupervised realistic deployment settings.
# ==============================================================================

# Mean imputation (baseline / simple benchmark):
# Replaces missing values with the column-wise mean computed from training cases.
imp_mean <- function(Xtr, Xte) {
  mu <- colMeans(Xtr, na.rm = TRUE)
  f <- function(X) {
    for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mu[j]
    X
  }
  list(list(train = f(Xtr), test = f(Xte)))
}

# k-Nearest Neighbors (kNN) imputation:
# Imputes missing cells using the average of the k nearest donors from Xtr.
# Distance metric: standardized Euclidean distance over all mutually observed features.
imp_knn <- function(Xtr, Xte, k = 5) {
  mu <- colMeans(Xtr, na.rm = TRUE)
  s <- apply(Xtr, 2, sd, na.rm = TRUE)
  s[s == 0] <- 1  # Guard against division by zero in constant features
  Ztr <- sweep(sweep(Xtr, 2, mu), 2, s, "/")

  f <- function(X, self) {
    Z <- sweep(sweep(X, 2, mu), 2, s, "/")
    for (i in which(!stats::complete.cases(X))) {
      obs <- which(!is.na(Z[i, ]))  # Features observed in the recipient row

      # Standardized Euclidean distance to all training donors
      d <- sqrt(rowMeans((Ztr[, obs, drop = FALSE] -
                            matrix(Z[i, obs], nrow(Ztr), length(obs), byrow = TRUE))^2,
                         na.rm = TRUE))
      if (self) d[i] <- Inf  # In the training fold, a row cannot serve as its own donor

      # Impute each missing cell with the average of the top k nearest donors
      for (j in which(is.na(X[i, ]))) {
        ok <- which(!is.na(Xtr[, j]) & is.finite(d))
        nn <- ok[order(d[ok])][seq_len(min(k, length(ok)))]
        X[i, j] <- mean(Xtr[nn, j])
      }
    }
    X
  }
  list(list(train = f(Xtr, TRUE), test = f(Xte, FALSE)))
}

# missForest: Non-parametric iterative Random Forest imputation
# 1. Fits missForest::missForest on training fold Xtr.
# 2. To impute the test fold Xte without leakage, trains separate randomForest models
#    on Xtr for each incomplete feature, and iteratively predicts missing cells in Xte
#    until convergence (relative squared change <= 1e-8) or maxiter iterations.
imp_missforest <- function(Xtr, Xte, seed = 1, ntree = 100, maxiter = 5) {
  set.seed(seed)
  # Fit missForest on the training fold
  tr <- as.matrix(missForest::missForest(as.data.frame(Xtr), ntree = ntree, verbose = FALSE)$ximp)

  # Impute test fold using random forests trained on the imputed training data
  te <- Xte
  Mte <- is.na(Xte)
  cols <- which(colSums(Mte) > 0)
  if (length(cols) > 0) {
    mu <- colMeans(tr)
    for (j in cols) te[Mte[, j], j] <- mu[j]  # Initial mean imputation for test

    # Train random forest models for each missing column
    rfs <- lapply(cols, function(j)
      randomForest::randomForest(x = tr[, -j, drop = FALSE], y = tr[, j], ntree = ntree))
    names(rfs) <- as.character(cols)

    # Iteratively update test predictions until convergence
    for (it in seq_len(maxiter)) {
      old <- te
      for (j in cols) te[Mte[, j], j] <- predict(rfs[[as.character(j)]], te[Mte[, j], -j, drop = FALSE])
      if (sum((te - old)^2) <= 1e-8 * sum(te^2)) break
    }
  }
  list(list(train = tr, test = te))
}

# MICE: Multivariate Imputation by Chained Equations
# Uses Predictive Mean Matching (PMM) with m imputations.
# To prevent test-data leakage, Xtr and Xte are stacked but test rows are flagged
# with ignore = ign: MICE fits imputation models strictly on training rows, but
# draws imputations for both training and test rows.
imp_mice <- function(Xtr, Xte, m = 5, seed = 1) {
  Xall <- as.data.frame(rbind(Xtr, Xte))
  ign <- c(rep(FALSE, nrow(Xtr)), rep(TRUE, nrow(Xte)))
  imp <- mice::mice(Xall, m = m, method = "pmm", ignore = ign, printFlag = FALSE, seed = seed)
  lapply(seq_len(m), function(k) {
    out <- as.matrix(mice::complete(imp, k))
    list(train = out[!ign, , drop = FALSE], test = out[ign, , drop = FALSE])
  })
}

# ==============================================================================
# 4. MODEL FITTING AND TIMING HELPERS
# ==============================================================================

# Minimum leaf/split size function: dynamically scales with sample size (>= 3)
n_min_fun <- function(n_train, prop_min = 0.01) as.integer(max(3, floor(prop_min * n_train)))

# Hyperparameter control builder for PRTree (matching Section 5 and 6 of paper)
# - max_terminal_nodes: 100
# - cp: complexity parameter 0.001
# - prop_hold: internal validation holdout fraction ~ 1 / (K - 1)
# - fill_type: 0 (substitute with proxy), 1 (distribute weights), 2 (impute via distribution)
pr_ctrl <- function(n_train, p, fill, K, proxy = "both", dist = "norm", prop_min = 0.01) {
  PRTree::pr_tree_control(
    max_terminal_nodes = 100, cp = 0.001, by_node = TRUE, tiny_sigma = 0,
    n_min = n_min_fun(n_train, prop_min),
    n_candidates = as.integer(max(round(p / 2), 3)),
    prop_hold = 1 / (K - 1),
    fill_type = fill, proxy_crit = proxy, dist = dist
  )
}

# Helper: extracts flat numeric vector from PRTree prediction object
pred_vec <- function(p) {
  if (is.list(p)) p <- if (!is.null(p$yhat)) p$yhat else p[[1]]
  as.numeric(p)
}

# Fits PRTree with specified internal fill strategy and returns test predictions
fit_pr <- function(Xtr, ytr, Xte, fill, K, seed, proxy = "both", dist = "norm", prop_min = 0.01) {
  set.seed(seed)
  fit <- PRTree::pr_tree(y = ytr, X = Xtr, control = pr_ctrl(nrow(Xtr), ncol(Xtr), fill, K, proxy, dist, prop_min))
  pred_vec(predict(fit, Xte))
}

# Fits standard CART (rpart) with surrogate splits, 10-fold CV pruning, and returns test predictions
fit_cart <- function(Xtr, ytr, Xte, seed, prop_min = 0.01) {
  ctrl <- rpart::rpart.control(cp = 0.001, xval = 10, maxdepth = 30,
                               minsplit = n_min_fun(nrow(Xtr), prop_min))
  set.seed(seed)
  fit <- rpart::rpart(y ~ ., data = data.frame(Xtr, y = ytr), method = "anova", control = ctrl)
  # Prune to optimal complexity parameter via minimum cross-validation error
  fit <- rpart::prune(fit, cp = fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"])
  as.numeric(predict(fit, newdata = as.data.frame(Xte)))
}

# Root Mean Squared Error (RMSE) evaluation metric
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

# Timing wrapper: executes an expression and returns list(value = result, time = elapsed_seconds)
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- expr
  list(value = v, time = proc.time()[["elapsed"]] - t0)
}

fmt_time <- function(sec) {
  if (!is.finite(sec)) return("--:--:--")
  sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, round(sec %% 60))
}

progress <- function(done, total, t0) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  eta <- if (done > 0) el / done * (total - done) else NA
  w <- 30
  k <- round(w * done / total)
  cat(sprintf("\r  [%s%s] %d/%d tasks | elapsed %s | ETA %s   ",
              strrep("=", k), strrep(" ", w - k), done, total, fmt_time(el), fmt_time(eta)))
  flush.console()
}

# ==============================================================================
# 5. CORE SIMULATION TASK WORKER
#
# Executes a single simulation task:
#   1. Injects missingness according to the designated mechanism and rate.
#   2. Randomly splits the data into K cross-validation folds.
#   3. For each fold k = 1, ..., K:
#      - Evaluates 18 model/imputation combinations:
#          * Complete-data baselines: PR-full, CART-full
#          * PRTree internal strategies: PR-f0, PR-f1, PR-f2
#          * CART with surrogate splits: CART-sur
#          * Single imputations + PRTree/CART: mean, kNN, missForest
#          * Multiple imputations (MICE):
#              - MICE+PR / MICE+CART: average of m=5 fitted trees
#              - MICE1+PR / MICE1+CART: tree fitted to single imputation m=1
#              - MICEavg+PR / MICEavg+CART: tree fitted to average imputed matrix
#      - Records out-of-sample RMSE, imputation time, and model fitting time.
# ==============================================================================
run_task_semisynthetic <- function(task, design, inject_fn, version) {
  d <- load_complete(task$dataset)

  # Standardized reproducible seed formula across all studies:
  # seed = base_seed + rep*1000 + mech_index*100 + round(rate*100)
  seed <- design$seed + task$rep * 1000 + match(task$mech, design$mechanisms) * 100 + round(task$rate * 100)
  set.seed(seed)

  # 1. Inject missing values into covariates
  Xm <- inject_fn(d$X, d$y, task$mech, task$rate, design$pi0)

  # 2. Partition data into K cross-validation folds
  folds <- sample(rep(seq_len(design$K), length.out = nrow(Xm)))
  res <- list()

  # Helper to record results per fold, tracking RMSE and runtimes separately
  add <- function(fold, method, yte, yhat, time_fit, time_imp = 0) {
    res[[length(res) + 1]] <<- data.frame(
      id = task$id, dataset = task$dataset, mech = task$mech, rate = task$rate,
      rep = task$rep, fold = fold, method = method, rmse = rmse(yte, yhat),
      time_imp = time_imp, time_fit = time_fit,
      pct_inc = mean(!stats::complete.cases(Xm)), version = version,
      stringsAsFactors = FALSE
    )
  }

  # 3. Cross-validation loop over K folds
  for (k in seq_len(design$K)) {
    tr <- folds != k
    te <- !tr
    ytr <- d$y[tr]
    yte <- d$y[te]
    s <- seed + k
    Xtr <- Xm[tr, , drop = FALSE]
    Xte <- Xm[te, , drop = FALSE]

    # --------------------------------------------------------------------------
    # A) Complete-data references (fitted to ground-truth X without missing values)
    # --------------------------------------------------------------------------
    f <- timed(fit_pr(d$X[tr, ], ytr, d$X[te, ], 2, design$K, s, design$proxy, design$dist, design$prop_min))
    add(k, "PR-full", yte, f$value, f$time)
    f <- timed(fit_cart(d$X[tr, ], ytr, d$X[te, ], s, design$prop_min))
    add(k, "CART-full", yte, f$value, f$time)

    # --------------------------------------------------------------------------
    # B) Internal PRTree strategies (f0, f1, f2) and CART with surrogate splits
    # --------------------------------------------------------------------------
    for (ft in 0:2) {
      f <- timed(fit_pr(Xtr, ytr, Xte, ft, design$K, s, design$proxy, design$dist, design$prop_min))
      add(k, paste0("PR-f", ft), yte, f$value, f$time)
    }
    f <- timed(fit_cart(Xtr, ytr, Xte, s, design$prop_min))
    add(k, "CART-sur", yte, f$value, f$time)

    # --------------------------------------------------------------------------
    # C) Impute-then-fit pipelines: single imputations (mean, kNN, missForest)
    # --------------------------------------------------------------------------
    for (nm in c("mean", "kNN", "missForest")) {
      im <- timed(switch(nm,
                         mean       = imp_mean(Xtr, Xte),
                         kNN        = imp_knn(Xtr, Xte, design$k_knn),
                         missForest = imp_missforest(Xtr, Xte, s))[[1]])
      D <- im$value

      # Fit PRTree (fill_type = 2) on the imputed data
      f <- timed(fit_pr(D$train, ytr, D$test, 2, design$K, s, design$proxy, design$dist, design$prop_min))
      add(k, paste0(nm, "+PR"), yte, f$value, f$time, im$time)

      # Fit CART on the imputed data
      f <- timed(fit_cart(D$train, ytr, D$test, s, design$prop_min))
      add(k, paste0(nm, "+CART"), yte, f$value, f$time, im$time)
    }

    # --------------------------------------------------------------------------
    # D) Multiple Imputation (MICE):
    #    - MICE+PR / MICE+CART: ensemble of m=5 trees fitted to 5 imputed datasets
    #    - MICE1+PR / MICE1+CART: single tree fitted to the 1st imputation
    #    - MICEavg+PR / MICEavg+CART: single tree fitted to mean of 5 imputed datasets
    # --------------------------------------------------------------------------
    im5 <- timed(imp_mice(Xtr, Xte, design$m_mice, s))
    L <- im5$value
    im1 <- timed(imp_mice(Xtr, Xte, 1, s))$time

    F_pr <- lapply(L, function(D) {
      timed(fit_pr(D$train, ytr, D$test, 2, design$K, s, design$proxy, design$dist, design$prop_min))
    })
    F_cart <- lapply(L, function(D) {
      timed(fit_cart(D$train, ytr, D$test, s, design$prop_min))
    })
    P_pr <- sapply(F_pr, `[[`, "value")
    T_pr <- sapply(F_pr, `[[`, "time")
    P_cart <- sapply(F_cart, `[[`, "value")
    T_cart <- sapply(F_cart, `[[`, "time")

    # Ensemble prediction averaging predictions across m trees
    add(k, "MICE+PR",    yte, rowMeans(P_pr),   sum(T_pr),   im5$time)
    add(k, "MICE+CART",  yte, rowMeans(P_cart), sum(T_cart), im5$time)

    # Single MICE imputation (m = 1)
    add(k, "MICE1+PR",   yte, P_pr[, 1],        T_pr[1],     im1)
    add(k, "MICE1+CART", yte, P_cart[, 1],      T_cart[1],   im1)

    # Average of imputed datasets (mean matrix)
    Davg <- list(train = Reduce(`+`, lapply(L, `[[`, "train")) / length(L),
                 test  = Reduce(`+`, lapply(L, `[[`, "test"))  / length(L))
    f <- timed(fit_pr(Davg$train, ytr, Davg$test, 2, design$K, s, design$proxy, design$dist, design$prop_min))
    add(k, "MICEavg+PR", yte, f$value, f$time, im5$time)
    f <- timed(fit_cart(Davg$train, ytr, Davg$test, s, design$prop_min))
    add(k, "MICEavg+CART", yte, f$value, f$time, im5$time)
  }
  dplyr::bind_rows(res)
}

# ==============================================================================
# 6. SIMULATION STUDY RUNNER (PARALLELIZED & CHECKPOINTED)
#
# Coordinates parallel execution of simulation tasks over a socket cluster.
# Features:
#   - Fault-tolerant checkpointing: saves results to <dataset>.rds after each batch.
#   - Restartability: skips tasks whose IDs already exist in the target .rds file.
#   - Clean error trapping with progress bar and elapsed/ETA estimation.
# ==============================================================================
run_semisynthetic_study <- function(design, inject_fn, res_dir, version, n_cores = 4) {
  if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
  res_file <- function(ds) file.path(res_dir, paste0(ds, ".rds"))

  # Check version compatibility to prevent mixing incompatible simulation outputs
  for (ds in design$datasets) {
    f <- res_file(ds)
    if (file.exists(f) && !identical(unique(readRDS(f)$version), version)) {
      stop(f, " was produced by an older version of the study. Delete ", res_dir, " to rerun.")
    }
  }

  # Expand simulation grid over replications, missing rates, mechanisms, and datasets
  grid <- expand.grid(rep = seq_len(design$R), rate = design$rates,
                      mech = design$mechanisms, dataset = design$datasets,
                      stringsAsFactors = FALSE)
  tasks <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    list(dataset = g$dataset, mech = g$mech, rate = g$rate, rep = g$rep,
         id = sprintf("%s_%s_r%02d_rep%02d", g$dataset, g$mech, round(100 * g$rate), g$rep))
  })

  # Initialize parallel cluster and export required objects
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("design", "version", "timed", "load_complete", "get_strong_predictors",
                      "inject_fn", "imp_mean", "imp_knn", "imp_missforest", "imp_mice",
                      "n_min_fun", "pr_ctrl", "pred_vec", "fit_pr", "fit_cart", "rmse",
                      "run_task_semisynthetic"), envir = environment())
  invisible(clusterEvalQ(cl, {
    library(PRTree)
    library(dplyr)
    NULL
  }))
  errors <- character(0)

  # Process dataset by dataset
  for (ds in design$datasets) {
    f <- res_file(ds)
    res_ds <- if (file.exists(f)) readRDS(f) else NULL
    done_ids <- unique(res_ds$id)
    todo <- Filter(function(t) t$dataset == ds && !(t$id %in% done_ids), tasks)
    n_tot <- sum(vapply(tasks, function(t) t$dataset == ds, logical(1)))
    cat(sprintf("\n%s: %d of %d tasks already done\n", ds, n_tot - length(todo), n_tot))
    if (length(todo) == 0) next

    t0 <- Sys.time()
    n_done <- 0
    progress(n_done, length(todo), t0)

    # Process tasks in batches equal to n_cores for periodic saving (checkpoints)
    batches <- split(todo, ceiling(seq_along(todo) / n_cores))
    for (b in batches) {
      out <- parLapply(cl, b, function(t) {
        tryCatch(
          run_task_semisynthetic(t, design, inject_fn, version),
          error = function(e) paste("ERROR", t$id, conditionMessage(e))
        )
      })
      ok <- vapply(out, is.data.frame, logical(1))
      errors <- c(errors, unlist(out[!ok]))
      if (any(ok)) {
        res_ds <- dplyr::bind_rows(res_ds, out[ok])
        saveRDS(res_ds, f)  # Save incremental checkpoint
      }
      n_done <- n_done + length(b)
      progress(n_done, length(todo), t0)
    }
    cat("\n")
  }
  stopCluster(cl)
  if (length(errors) > 0) {
    print(errors)
    stop("Some tasks failed. Check errors above.")
  }
}

# ==============================================================================
# 7. SUMMARY, FIGURES AND LATEX TABLES GENERATOR
#
# Generates all publication-quality outputs and prints numbers quoted in the text:
#   - semisynthetic_boxplots.pdf: Figure 4 in the paper (relative RMSE boxplots)
#   - tab_semisynthetic.tex     : Table 3 in the paper (median & IQR of relative RMSE)
#   - tab_semisynthetic_ci.tex  : Table 3 with 95% replication-level bootstrap CIs of the medians
#   - tab_semisynthetic_winci.csv: win rate of PR-f2 against each competitor with 95% bootstrap CIs
#   - semisynthetic_heatmaps.pdf: Figure 5 in the paper (PR-f2 pairwise win rate heatmaps)
#   - tab_semisynthetic_time.tex: Table 4 in the paper (runtime per fold & ratio to PR-f2)
#   - semisynthetic_time.pdf    : Figure 6 in the paper (runtime per fold heatmaps)
#   - Terminal summaries        : Pairwise ranking, complete-data sanity check,
#                                 cost-gain ratios, and numbers cited in Section 6.
# ==============================================================================

# Replication-level cluster bootstrap. Each replication draws its own missingness mask and
# fold partition, so replications are independent whereas the folds of one replication are
# not. The bootstrap therefore resamples whole replications (keeping all their folds) and
# returns the percentile interval of stat() over n_boot resamples.
boot_rep_ci <- function(x, rep, stat, n_boot = 2000, level = 0.95) {
  idx <- split(x, factor(rep, unique(rep)))
  n_rep <- length(idx)
  draws <- matrix(sample.int(n_rep, n_rep * n_boot, replace = TRUE), nrow = n_rep, ncol = n_boot)

  b <- apply(draws, 2, function(i) {
    stat(unlist(idx[i], use.names = FALSE))
  })

  alpha <- (1 - level) / 2
  c(
    est = stat(x),
    lo  = unname(quantile(b, alpha)),
    hi  = unname(quantile(b, 1 - alpha))
  )
}

generate_semisynthetic_summaries <- function(design, res_dir, version, n_boot = 2000, seed_boot = 2026) {
  res_file <- function(ds) file.path(res_dir, paste0(ds, ".rds"))
  labels <- c(Auto = "Auto MPG", Ozone = "Ozone", Concrete = "Concrete")

  # Load and combine results across benchmark datasets
  d <- dplyr::bind_rows(lapply(res_file(design$datasets), readRDS)) %>%
    mutate(mech = factor(mech, design$mechanisms),
           dataset = factor(labels[dataset], labels[design$datasets]),
           rate_lab = factor(paste0(100 * rate, "%"), paste0(100 * design$rates, "%")))

  # Reshape to wide format for fold-by-fold pairwise comparisons
  w <- d %>% select(dataset, mech, rate, rep, fold, method, rmse) %>%
    pivot_wider(names_from = method, values_from = rmse)

  # ----------------------------------------------------------------------------
  # Sanity Check: CART-full vs CART-sur
  # Evaluates how often CART without missingness produces identical predictions
  # to CART with surrogate splits on data with missingness (expected ~ 0%).
  # ----------------------------------------------------------------------------
  cat("\nSanity check - share of folds where CART-full and CART-sur coincide (should be ~0):\n")
  print(w %>% group_by(dataset) %>% summarise(identical = mean(abs(`CART-full` - `CART-sur`) < 1e-12)))

  # Short label mapping for compact tables and heatmaps
  short <- c("CART-full" = "CART-full", "PR-full" = "full", "PR-f0" = "f0", "PR-f1" = "f1", "PR-f2" = "f2",
             "mean+PR" = "mean", "kNN+PR" = "kNN", "missForest+PR" = "mF",
             "MICE+PR" = "MICE", "MICE1+PR" = "MICE1", "MICEavg+PR" = "MICEavg",
             "CART-sur" = "CART-sur", "MICE+CART" = "MICE+CART")

  cand <- c("PR-f0", "PR-f1", "PR-f2", "CART-sur",
            paste0(c("mean", "kNN", "missForest", "MICE", "MICE1", "MICEavg"), "+PR"),
            paste0(c("mean", "kNN", "missForest", "MICE", "MICE1", "MICEavg"), "+CART"))
  xx <- d %>% filter(method %in% cand) %>% select(dataset, mech, rate, rep, fold, method, rmse)

  # ----------------------------------------------------------------------------
  # Overall Pairwise Win Rate Ranking
  # Calculates the fraction of folds in which method A achieved lower RMSE than method B
  # ----------------------------------------------------------------------------
  ranking <- xx %>%
    inner_join(xx, by = c("dataset", "mech", "rate", "rep", "fold"), relationship = "many-to-many") %>%
    filter(method.x != method.y) %>%
    group_by(method = method.x) %>% summarise(win = mean(rmse.x < rmse.y), .groups = "drop") %>%
    arrange(desc(win))
  cat("\nRanking of the procedures (overall pairwise win rate):\n")
  print(ranking, n = 30)

  # Select the top 6 procedures (ensuring PR-f2 is included) for figures and tables
  best <- head(ranking$method, 6)
  if (!"PR-f2" %in% best) best <- c(head(best, 5), "PR-f2")
  best <- ranking$method[ranking$method %in% best]
  nice <- function(m) {
    sub("\\+PR$", " + PRTree", sub("^PR-f2$", "fill_type = 2", sub("\\+CART$", " + CART", m)))
  }

  # ----------------------------------------------------------------------------
  # Relative RMSE Loss: E = RMSE / RMSE(PRTree, complete data) - 1
  # Represents the proportional inflation in test prediction error due to missingness
  # ----------------------------------------------------------------------------
  loss <- d %>% group_by(dataset, mech, rate, rep, fold) %>%
    mutate(E = rmse / rmse[method == "PR-full"] - 1) %>% ungroup()

  pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  bx <- loss %>% filter(method %in% best) %>% mutate(method = factor(nice(method), nice(best)))

  # Publication color palette (Okabe-Ito / ColorBrewer)
  all_cols <- c("PR-f2" = "#CC0000", "MICE+PR" = "#6A3D9A", "missForest+PR" = "#009E73",
                "kNN+PR" = "#0072B2", "MICEavg+PR" = "#E69F00", "MICE1+PR" = "#CC79A7",
                "mean+PR" = "#999999", "MICE+CART" = "#56B4E9", "missForest+CART" = "#A6D96A",
                "kNN+CART" = "#9ECAE1", "MICEavg+CART" = "#FDD49E", "MICE1+CART" = "#F1B6DA",
                "mean+CART" = "#D9D9D9", "CART-sur" = "#444444", "PR-f1" = "#F4A582", "PR-f0" = "#B2182B")
  cols <- setNames(all_cols[best], nice(best))
  lim <- c("Auto MPG" = 1.25, "Ozone" = Inf, "Concrete" = Inf)
  n_m <- nlevels(bx$method)
  w_d <- 0.85

  # ----------------------------------------------------------------------------
  # Figure 4: Boxplots of Relative RMSE Loss across missing rates & mechanisms
  # ----------------------------------------------------------------------------
  plot_box <- function(dd, first, last) {
    ymax <- max(dd$E[dd$E < lim[[as.character(dd$dataset[1])]]])
    above <- dd %>% group_by(dataset, mech, rate_lab, method) %>%
      summarise(n_out = sum(E > ymax), .groups = "drop") %>% filter(n_out > 0) %>%
      mutate(y = ymax + 0.05, label = paste0(n_out, "\u25B2"),
             xpos = as.numeric(rate_lab) + (as.numeric(method) - (n_m + 1) / 2) * w_d / n_m)
    ggplot(dd, aes(rate_lab, E, fill = method)) +
      geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
      geom_boxplot(alpha = 0.85, outlier.size = 0.4, linewidth = 0.3,
                   position = position_dodge(width = w_d), width = 0.75) +
      geom_text(data = above, aes(x = xpos, y = y, label = label), inherit.aes = FALSE, colour = "red", size = 3) +
      scale_fill_manual(values = cols, name = NULL) +
      facet_grid(dataset ~ mech) +
      coord_cartesian(ylim = c(min(dd$E), 1.1 * ymax)) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.text.x = if (first) element_text(face = "bold") else element_blank(),
            axis.text.x = if (last) element_text() else element_blank(),
            axis.ticks.x = if (last) element_line() else element_blank()) +
      guides(fill = guide_legend(nrow = 1))
  }

  pb <- ggpubr::ggarrange(
    plotlist = lapply(seq_along(levels(bx$dataset)), function(i)
      plot_box(bx[bx$dataset == levels(bx$dataset)[i], ], i == 1, i == nlevels(bx$dataset))),
    ncol = 1, align = "v", common.legend = TRUE, legend = "bottom")
  pb <- ggpubr::annotate_figure(pb, left = "RMSE / RMSE(PRTree, complete data) - 1",
                                bottom = "Proportion of incomplete rows r")
  ggsave(file.path(res_dir, "semisynthetic_boxplots.pdf"), pb, width = 11, height = 8, device = pdf_dev)

  # ----------------------------------------------------------------------------
  # Table 3: LaTeX Table of Median (IQR) Relative RMSE Loss
  # Bold denotes the best (lowest median loss) method in each scenario row.
  # ----------------------------------------------------------------------------
  tab_cols <- c(best[!grepl("CART", best)], best[grepl("CART", best)])
  tb <- loss %>% filter(method %in% best) %>%
    group_by(dataset, mech, rate, method) %>%
    summarise(med = median(E), iqr = IQR(E), inc = 100 * mean(pct_inc), .groups = "drop") %>%
    group_by(dataset, mech, rate) %>% mutate(inc = mean(inc), best_m = med == min(med)) %>% ungroup() %>%
    mutate(med = round(med, 2) + 0, iqr = round(iqr, 2) + 0,
           cell = ifelse(best_m, sprintf("\\textbf{%.2f} (%.2f)", med, iqr), sprintf("%.2f (%.2f)", med, iqr)),
           method = factor(method, tab_cols)) %>%
    select(dataset, mech, rate, inc, method, cell) %>%
    pivot_wider(names_from = method, values_from = cell) %>% arrange(dataset, mech, rate)

  lines <- c()
  for (ds in levels(tb$dataset)) {
    sub <- tb %>% filter(dataset == ds)
    if (nrow(sub) == 0) {
      next
    }
    lines <- c(lines, sprintf("\\multicolumn{%d}{l}{\\textit{%s}}\\\\", length(tab_cols) + 3, ds))
    for (i in seq_len(nrow(sub))) {
      mech <- if (i == 1 || sub$mech[i] != sub$mech[i - 1]) as.character(sub$mech[i]) else ""
      lines <- c(lines, paste0(paste(c(mech, round(100 * sub$rate[i]), sprintf("%.0f", sub$inc[i]),
                                       unlist(sub[i, tab_cols])), collapse = " & "), "\\\\"))
    }
    lines <- c(lines, "\\midrule")
  }
  lines[length(lines)] <- "\\bottomrule"

  hdr <- c("PR-f2" = "\\texttt{f2}", "PR-f1" = "\\texttt{f1}", "PR-f0" = "\\texttt{f0}",
           "mean+PR" = "mean", "kNN+PR" = "$k$NN", "missForest+PR" = "mF", "MICE+PR" = "MICE",
           "MICE1+PR" = "MICE$_1$", "MICEavg+PR" = "MICE$_{\\mathrm{avg}}$",
           "CART-sur" = "sur.", "mean+CART" = "mean", "kNN+CART" = "$k$NN", "missForest+CART" = "mF",
           "MICE+CART" = "MICE", "MICE1+CART" = "MICE$_1$", "MICEavg+CART" = "MICE$_{\\mathrm{avg}}$")
  is_cart <- grepl("CART", tab_cols)
  grp <- c(sprintf("\\multicolumn{%d}{c}{PRTree}", sum(!is_cart)),
           if (any(is_cart)) sprintf("\\multicolumn{%d}{c}{CART}", sum(is_cart)))
  rules <- c(sprintf("\\cmidrule(lr){4-%d}", 3 + sum(!is_cart)),
             if (any(is_cart)) sprintf("\\cmidrule(lr){%d-%d}", 4 + sum(!is_cart), 3 + length(tab_cols)))
  writeLines(c(sprintf("\\begin{tabular}{llr%s}", strrep("c", length(tab_cols))), "\\toprule",
               paste0(" & & & ", paste(grp, collapse = " & "), "\\\\"),
               paste(rules, collapse = " "),
               paste0("Mech. & $r$ & Inc. & ", paste(hdr[tab_cols], collapse = " & "), "\\\\"),
               "\\midrule", lines, "\\end{tabular}"),
             file.path(res_dir, "tab_semisynthetic.tex"))

  # ----------------------------------------------------------------------------
  # Table 3 with uncertainty: 95% replication-level bootstrap CIs of the median loss
  # Same layout as tab_semisynthetic.tex but with three decimals and \best{} instead of \textbf{};
  # each scenario row is followed by a row with the percentile interval [lo, hi]
  # of the median (replications resampled as a whole).
  # ----------------------------------------------------------------------------
  set.seed(seed_boot)

  ci <- loss %>%
    filter(method %in% best) %>%
    group_by(dataset, mech, rate, method) %>%
    group_modify(~ as.data.frame(as.list(boot_rep_ci(.x$E, .x$rep, median, n_boot)))) %>%
    ungroup() %>%
    mutate(
      method = factor(method, tab_cols),
      ci     = sprintf("{\\scriptsize [%.3f, %.3f]}", round(lo, 3) + 0, round(hi, 3) + 0)
    ) %>%
    select(dataset, mech, rate, method, ci) %>%
    pivot_wider(names_from = method, values_from = ci)

  # Median (IQR) with three decimals. The best method (lowest median) is wrapped in \best{},
  # which must be defined in the paper preamble (e.g. \newcommand{\best}[1]{\textbf{#1}}),
  # so that its emphasis can be changed independently of any other markup in the manuscript.
  tb3 <- loss %>%
    filter(method %in% best) %>%
    group_by(dataset, mech, rate, method) %>%
    summarise(
      med = median(E),
      iqr = IQR(E),
      .groups = "drop"
    ) %>%
    group_by(dataset, mech, rate) %>%
    mutate(best_m = (med == min(med))) %>%
    ungroup() %>%
    mutate(
      cell = ifelse(
        best_m,
        sprintf("\\best{%.3f (%.3f)}", round(med, 3) + 0, round(iqr, 3) + 0),
        sprintf("%.3f (%.3f)", round(med, 3) + 0, round(iqr, 3) + 0)
      ),
      method = factor(method, tab_cols)
    ) %>%
    select(dataset, mech, rate, method, cell) %>%
    pivot_wider(names_from = method, values_from = cell)

  lines_ci <- c()
  for (ds in levels(tb$dataset)) {
    sub <- tb %>% filter(dataset == ds)
    if (nrow(sub) == 0) {
      next
    }

    lines_ci <- c(lines_ci, sprintf("\\multicolumn{%d}{l}{\\textit{%s}}\\\\", length(tab_cols) + 3, ds))

    for (i in seq_len(nrow(sub))) {
      mech <- if (i == 1 || sub$mech[i] != sub$mech[i - 1]) {
        as.character(sub$mech[i])
      } else {
        ""
      }

      ci_i  <- ci %>% filter(dataset == ds, mech == sub$mech[i], rate == sub$rate[i])
      tb3_i <- tb3 %>% filter(dataset == ds, mech == sub$mech[i], rate == sub$rate[i])

      row_vals_est <- c(
        mech,
        round(100 * sub$rate[i]),
        sprintf("%.0f", sub$inc[i]),
        unlist(tb3_i[1, tab_cols])
      )
      row_vals_ci <- c(
        "", "", "",
        unlist(ci_i[1, tab_cols])
      )

      lines_ci <- c(
        lines_ci,
        paste0(paste(row_vals_est, collapse = " & "), "\\\\"),
        paste0(paste(row_vals_ci, collapse = " & "), "\\\\")
      )
    }

    lines_ci <- c(lines_ci, "\\midrule")
  }
  lines_ci[length(lines_ci)] <- "\\bottomrule"

  writeLines(
    c(
      sprintf("\\begin{tabular}{llr%s}", strrep("c", length(tab_cols))),
      "\\toprule",
      paste0(" & & & ", paste(grp, collapse = " & "), "\\\\"),
      paste(rules, collapse = " "),
      paste0("Mech. & $r$ & Inc. & ", paste(hdr[tab_cols], collapse = " & "), "\\\\"),
      "\\midrule",
      lines_ci,
      "\\end{tabular}"
    ),
    file.path(res_dir, "tab_semisynthetic_ci.tex")
  )

  # ----------------------------------------------------------------------------
  # Figure 5: Heatmaps of PRTree (fill_type = 2) Pairwise Win Rate
  # Displays the percentage of folds where PR-f2 achieved lower RMSE than competitor.
  # ----------------------------------------------------------------------------
  grid_heatmaps <- function(dat, make_one, legend_title) {
    hm <- lapply(levels(dat$dataset), function(ds) make_one(dplyr::filter(dat, dataset == ds)) +
                   labs(title = ds, x = NULL, y = NULL) + theme_bw(base_size = 14) +
                   theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 14)))
    leg <- ggpubr::get_legend(hm[[1]] + guides(fill = guide_colorbar(
      direction = "horizontal", title.position = "top", title = legend_title,
      barwidth = unit(14, "cm"), barheight = unit(1.5, "cm"))))
    hm <- lapply(hm, function(p) p + theme(legend.position = "none"))
    ggpubr::ggarrange(plotlist = c(hm, list(ggpubr::as_ggplot(leg))), ncol = 2,
                      nrow = ceiling((length(hm) + 1) / 2))
  }

  comp <- c("PR-f0", "PR-f1", "CART-sur", "mean+PR", "kNN+PR", "missForest+PR",
            "MICE+PR", "MICE1+PR", "MICEavg+PR", "MICE+CART")
  hw <- w %>% pivot_longer(all_of(comp), names_to = "competitor", values_to = "rmse_c") %>%
    group_by(dataset, mech, rate, competitor) %>%
    summarise(win = mean(`PR-f2` < rmse_c), .groups = "drop") %>%
    mutate(competitor = factor(short[competitor], rev(short[comp])),
           col = factor(paste0(mech, "\n", 100 * rate, "%"),
                        paste0(rep(design$mechanisms, each = length(design$rates)), "\n",
                               100 * design$rates, "%")))

  # 95% replication-level bootstrap CIs of the win rate of PR-f2 (all rates pooled)
  win_ci <- w %>%
    pivot_longer(all_of(comp), names_to = "competitor", values_to = "rmse_c") %>%
    mutate(win = as.numeric(`PR-f2` < rmse_c)) %>%
    group_by(dataset, mech, competitor) %>%
    group_modify(~ as.data.frame(as.list(boot_rep_ci(.x$win, paste(.x$rate, .x$rep), mean, n_boot)))) %>%
    ungroup() %>%
    mutate(across(c(est, lo, hi), ~ round(., 3)))

  write.csv(win_ci, file.path(res_dir, "tab_semisynthetic_winci.csv"), row.names = FALSE)

  cat("\nWin rate of PR-f2 with 95% replication-level bootstrap CI (rates pooled):\n")
  print(
    win_ci %>% filter(competitor %in% c("MICE+PR", "MICE+CART", "CART-sur", "PR-f0", "PR-f1")),
    n = 100
  )

  ph <- grid_heatmaps(hw, function(x) ggplot(x, aes(col, competitor, fill = win)) +
                        geom_tile(colour = "white") +
                        geom_text(aes(label = sprintf("%.0f%%", 100 * win)), size = 4.5) +
                        scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2",
                                             midpoint = 0.5, limits = c(0, 1)),
                      "f2 win rate")
  ggsave(file.path(res_dir, "semisynthetic_heatmaps.pdf"), ph, width = 15, height = 10, device = pdf_dev)

  # ---- numbers quoted in the text ----------------------------------------------
  cat("\nDataset sizes:\n")
  print(dplyr::bind_rows(lapply(design$datasets, function(n) { x <- load_complete(n)
    data.frame(dataset = n, n = nrow(x$X), p = ncol(x$X)) })))
  cat("\nComplete-data references: mean RMSE and share of folds with PR-full < CART-full:\n")
  print(w %>% group_by(dataset) %>%
          summarise(PR_full = mean(`PR-full`), CART_full = mean(`CART-full`),
                    win = mean(`PR-full` < `CART-full`)))
  cat("\nMedian RMSE/RMSE(PR-full) - 1 by dataset, mechanism and rate:\n")
  print(d %>% group_by(dataset, mech, rate, rep, fold) %>%
          mutate(E = rmse / rmse[method == "PR-full"] - 1) %>% ungroup() %>%
          filter(method %in% c(cand, "CART-full")) %>%
          group_by(dataset, mech, rate, method) %>% summarise(med = round(median(E), 3), .groups = "drop") %>%
          pivot_wider(names_from = method, values_from = med), n = 30, width = Inf)
  cat("\nWin rate of PR-f2 against each competitor, by dataset and mechanism:\n")
  print(hw %>% group_by(dataset, mech, competitor) %>% summarise(win = round(mean(win), 2), .groups = "drop") %>%
          pivot_wider(names_from = competitor, values_from = win), n = 30, width = Inf)
  cat("\nWin rate of PR-f2 against each competitor, by dataset, mechanism and rate:\n")
  print(hw %>% mutate(win = round(win, 2)) %>% select(-col) %>%
          pivot_wider(names_from = competitor, values_from = win), n = 100, width = Inf)
  cat("\nPRTree vs CART after the same imputation (share of folds with lower RMSE):\n")
  print(w %>% group_by(dataset) %>%
          summarise(mean = mean(`mean+PR` < `mean+CART`), kNN = mean(`kNN+PR` < `kNN+CART`),
                    mF = mean(`missForest+PR` < `missForest+CART`), MICE = mean(`MICE+PR` < `MICE+CART`),
                    MICE1 = mean(`MICE1+PR` < `MICE1+CART`), MICEavg = mean(`MICEavg+PR` < `MICEavg+CART`)) %>%
          mutate(across(where(is.numeric), ~ round(., 2))), width = Inf)
  cat("\nShare of folds in which MICE+PR (5 trees) beats PR-full (complete data):\n")
  print(w %>% group_by(dataset, mech) %>% summarise(share = round(mean(`MICE+PR` < `PR-full`), 2), .groups = "drop"))

  # Computational cost table and heatmap
  tm <- d %>% mutate(total = time_imp + time_fit) %>%
    group_by(dataset, method) %>%
    summarise(imp = median(time_imp), fit = median(time_fit), med = median(total),
              q1 = quantile(total, 0.25), q3 = quantile(total, 0.75), .groups = "drop")
  ref_f2 <- tm %>% filter(method == "PR-f2") %>% select(dataset, f2 = med)
  tm <- tm %>% left_join(ref_f2, by = "dataset") %>% mutate(ratio_f2 = med / f2)

  cat("\nComputational cost per fold (seconds): median imputation, fitting and total time,\n",
      "quartiles of the total, and ratio of the median total to that of PR-f2:\n")
  print(tm %>% mutate(across(c(imp, fit, med, q1, q3), ~ round(., 2)), ratio_f2 = round(ratio_f2, 2)) %>%
          select(-f2) %>% arrange(dataset, med), n = 100, width = Inf)

  time_cols <- c("PR-f2", "PR-f1", "PR-f0", "mean+PR", "kNN+PR", "missForest+PR", "MICE1+PR",
                 "MICEavg+PR", "MICE+PR", "CART-sur", "MICE+CART")
  tt <- tm %>% filter(method %in% time_cols) %>%
    mutate(cell = sprintf("%.2f (%.2f)", med, ratio_f2), method = factor(method, time_cols)) %>%
    select(dataset, method, cell) %>% pivot_wider(names_from = dataset, values_from = cell) %>%
    arrange(method)
  tlab <- c("PR-f2" = "\\texttt{fill\\_type} = 2", "PR-f1" = "\\texttt{fill\\_type} = 1",
            "PR-f0" = "\\texttt{fill\\_type} = 0", "mean+PR" = "mean + PRTree", "kNN+PR" = "$k$NN + PRTree",
            "missForest+PR" = "missForest + PRTree", "MICE1+PR" = "MICE$_1$ + PRTree",
            "MICEavg+PR" = "MICE$_{\\mathrm{avg}}$ + PRTree", "MICE+PR" = "MICE + PRTree",
            "CART-sur" = "CART (surrogate splits)", "MICE+CART" = "MICE + CART")
  writeLines(c(sprintf("\\begin{tabular}{l%s}", strrep("c", ncol(tt) - 1)), "\\toprule",
               paste0("Procedure & ", paste(names(tt)[-1], collapse = " & "), "\\\\"), "\\midrule",
               apply(tt, 1, function(r) paste0(tlab[r[["method"]]], " & ", paste(r[-1], collapse = " & "), "\\\\")),
               "\\bottomrule", "\\end{tabular}"),
             file.path(res_dir, "tab_semisynthetic_time.tex"))

  time_lab <- c("MICE+PR" = "MICE", "MICEavg+PR" = "MICEavg", "MICE1+PR" = "MICE1", "kNN+PR" = "kNN",
                "missForest+PR" = "mF", "mean+PR" = "mean", "PR-f2" = "f2", "PR-f1" = "f1", "PR-f0" = "f0",
                "CART-sur" = "CART-sur", "MICE+CART" = "MICE+CART", "MICEavg+CART" = "MICEavg+CART",
                "MICE1+CART" = "MICE1+CART", "kNN+CART" = "kNN+CART", "missForest+CART" = "mF+CART",
                "mean+CART" = "mean+CART")
  ht <- d %>% filter(method %in% best) %>%
    mutate(total = time_imp + time_fit) %>%
    group_by(dataset, mech, rate, method) %>% summarise(t = median(total), .groups = "drop") %>%
    mutate(method = factor(time_lab[method], rev(time_lab[best])),
           col = factor(paste0(mech, "\n", 100 * rate, "%"),
                        paste0(rep(design$mechanisms, each = length(design$rates)), "\n",
                               100 * design$rates, "%")))
  t_lim <- range(pmax(ht$t, 0.005))
  pt <- grid_heatmaps(ht, function(x) ggplot(x, aes(col, method, fill = pmax(t, 0.005))) +
                        geom_tile(colour = "white") +
                        geom_text(aes(label = ifelse(t < 10, sprintf("%.2f", t), sprintf("%.1f", t))), size = 4.5) +
                        scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                                             trans = "log10", midpoint = 0, limits = t_lim),
                      "Median time per fold (s)")
  ggsave(file.path(res_dir, "semisynthetic_time.pdf"), pt, width = 15, height = 7, device = pdf_dev)

  # ---- accuracy versus cost, relative to PR-f2 ----------------------------------
  #   gain G = 1 - RMSE / RMSE(PR-f2) (fold by fold; positive = more accurate than PR-f2)
  #   cost   = median total time per fold / median total time of PR-f2 (by dataset)
  wb <- d %>% filter(method %in% best) %>% mutate(total = time_imp + time_fit) %>%
    select(dataset, mech, rate, rep, fold, method, rmse, total)
  cb <- wb %>% filter(method != "PR-f2") %>%
    left_join(wb %>% filter(method == "PR-f2") %>%
                select(dataset, mech, rate, rep, fold, rmse_f2 = rmse),
              by = c("dataset", "mech", "rate", "rep", "fold")) %>%
    mutate(G = 1 - rmse / rmse_f2)
  cost <- wb %>% group_by(dataset, method) %>% summarise(t = median(total), .groups = "drop") %>%
    group_by(dataset) %>% mutate(cost = t / t[method == "PR-f2"]) %>% ungroup() %>%
    select(dataset, method, cost)
  cat("\nAccuracy versus cost relative to PR-f2: median gain by mechanism and cost ratio:\n")
  print(cb %>% group_by(dataset, mech, method) %>% summarise(G = round(median(G), 3), .groups = "drop") %>%
          pivot_wider(names_from = mech, values_from = G) %>%
          left_join(cost %>% mutate(cost = round(cost, 2)), by = c("dataset", "method")) %>%
          arrange(dataset, method), n = 50, width = Inf)
  cat("\nMedian gain relative to PR-f2 under MNAR, by rate:\n")
  print(cb %>% filter(mech == "MNAR") %>% group_by(dataset, rate, method) %>%
          summarise(G = round(median(G), 3), .groups = "drop") %>%
          pivot_wider(names_from = method, values_from = G), n = 50, width = Inf)

  cat(sprintf("\nOutputs generated successfully in %s:\n", res_dir))
  cat("  - semisynthetic_boxplots.pdf\n  - semisynthetic_heatmaps.pdf\n  - semisynthetic_time.pdf\n")
  cat("  - tab_semisynthetic.tex\n  - tab_semisynthetic_ci.tex\n  - tab_semisynthetic_winci.csv\n  - tab_semisynthetic_time.tex\n")
}
