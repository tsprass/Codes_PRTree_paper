# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT 9: SEMI-SYNTHETIC STUDY, ONE-SIDED MNAR IN THE STRONGEST PREDICTORS
# (sensitivity analysis; scripts 7 and 8 are left untouched)
#
# Same as script 8 (missing values only in the ceil(p/2) predictors most correlated with y),
# but MAR and MNAR are ONE-SIDED (upper tail):
#   MNAR: in each predictor of A the HIGHEST values are deleted, as if the measuring device
#         saturated above a ceiling and recorded nothing. The number of deleted values per
#         predictor is calibrated so that the proportion of incomplete rows is r.
#   MAR : the r*n rows with the HIGHEST max_k z_ik over the always observed predictors O
#         (no absolute value) are flagged, and in each flagged row every predictor in A is
#         deleted with probability pi0 (draw repeated if nothing is deleted), as if the
#         device failed whenever a recorded variable was high.
#   MCAR: as in script 8 (random cells of A).
# The labels "MCAR", "MAR" and "MNAR" are kept so that the summaries still run.
#
# Same procedures, folds, rates and outputs as script 7. The only difference is WHERE
# the missing values are placed: the ceil(p/2) predictors most correlated with the
# response (|cor(X_j, y)|, computed once on the complete data) form the set A and are the
# only ones that can lose values; the remaining floor(p/2) predictors (set O) are always
# observed. In the MAR mechanism O is the set of guides. The response only chooses A and O;
# it never decides which cells are deleted.
#
# Internal strategies (fill_type 0/1/2) versus impute-then-PRTree and
# impute-then-CART under MCAR / MAR / MNAR missingness injected into fully
# observed predictors (complete cases of Auto MPG and Ozone; Concrete).
#
# Procedures (per fold):
#   PR-full, CART-full            : fitted to the complete data (reference)
#   PR-f0, PR-f1, PR-f2           : PRTree with internal strategies
#   CART-sur                      : CART with surrogate splits
#   X+PR, X+CART, X in {mean, kNN, missForest}: impute, then fit one tree
#   MICE+PR, MICE+CART            : m = 5 imputations, predictions averaged over 5 trees
#   MICE1+PR, MICE1+CART          : one tree fitted to the first MICE imputation
#   MICEavg+PR, MICEavg+CART      : one tree fitted to the average of the 5 imputations
# Imputation models use the TRAINING fold only and never use the response.
#
# Output (folder results_semisynthetic_upper):
#   <dataset>.rds            raw results (one file per dataset, checkpointed)
#   semisynthetic_boxplots.pdf  boxplots of RMSE/RMSE(PR-full) - 1 for the 6 best procedures
#   semisynthetic_heatmaps.pdf  win rate of PR-f2 against each competitor, as in Section 5
#   tab_semisynthetic.tex       median (IQR) of RMSE/RMSE(PR-full) - 1 for the 6 best procedures
#   tab_semisynthetic_time.tex  median time per fold (s) and ratio to fill_type = 2
#   semisynthetic_time.pdf      heatmap of the median time per fold for the 6 best procedures
#   numbers quoted in the text are printed at the end
#
# Run from the folder "Applic". START FROM AN EMPTY results_semisynthetic_upper FOLDER.
# The script can be interrupted and restarted: finished tasks are skipped.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(parallel)
library(dplyr)

library(tidyr)
library(ggplot2)
library(ggpubr)     # combines the panels of the boxplot figures

res_dir <- "results_semisynthetic_upper"
if (!dir.exists(res_dir)) dir.create(res_dir)
version <- "rev1-upper" # stamp stored in the results files

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DESIGN (edit here to reduce cost)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
design <- list(
  datasets   = c("Auto", "Ozone", "Concrete"),
  mechanisms = c("MCAR", "MAR", "MNAR"),  # MAR and MNAR one-sided (upper tail)
  rates      = c(0.10, 0.25, 0.40),  # proportion of incomplete rows
  R          = 10,        # replications (new mask + new folds each)
  K          = 5,         # folds
  pi0        = 0.5,       # MAR: probability of deleting each non-guide predictor in a flagged row
  k_knn      = 5,         # neighbours for kNN imputation
  m_mice     = 5,         # number of multiple imputations
  prop_min   = 0.01,      # n_min = 1% of training size (as in Section 5)
  dist       = "norm",    # PRTree smoothing density
  proxy      = "both",    # PRTree proxy criterion
  seed       = 2026
)
n_cores <- 4

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA: complete predictor matrices
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_complete <- function(name) {
  f1 <- file.path("data/individual", paste0(name, ".rds"))
  f2 <- file.path("data/individual_complete", paste0(name, ".rds"))
  d <- readRDS(if (file.exists(f1)) f1 else f2)
  cc <- stats::complete.cases(d$X)          # complete cases only
  X <- d$X[cc, , drop = FALSE]
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  list(name = name, X = X, y = d$y[cc])
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MISSINGNESS INJECTION (only in the strong predictors A; O = complement, always observed)
#   r = target proportion of INCOMPLETE ROWS (rows with at least one missing predictor);
#   the response is never involved in the deletion.
#   MCAR: predictors in A; every cell is deleted with the same probability
#         q = 1 - (1 - r)^(1/|A|), independently of the data.
#   MAR : the r*n rows with the largest max_k z_ik (one-sided, no absolute value), k in O
#         (ties broken at random), are flagged and, in each flagged row, every predictor in A is deleted with
#         probability pi0 (the draw is repeated if nothing is deleted). Depends only on
#         observed values.
#   MNAR (one-sided): predictors in A; in each one the highest values are deleted (ranks
#         with random tie-breaking, so that discrete predictors are handled), about q*n
#         cells per predictor, with q found by bisection so that the proportion of
#         incomplete rows is r.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inject_missing <- function(X, mech, r, pi0, A) {
  n <- nrow(X); p <- ncol(X)
  O <- setdiff(seq_len(p), A)
  na <- length(A)
  Z <- scale(X)
  M <- matrix(FALSE, n, p)
  if (mech == "MCAR") {
    q <- 1 - (1 - r)^(1 / na)
    M[, A] <- matrix(runif(n * na) < q, n, na)
  } else if (mech == "MAR") {
    s <- apply(Z[, O, drop = FALSE], 1, max)                  # how high the guides are
    flag <- rank(-s, ties.method = "random") <= round(r * n)  # the r*n most extreme rows
    for (i in which(flag)) {
      repeat { d <- runif(na) < pi0; if (any(d)) break }
      M[i, A[d]] <- TRUE
    }
  } else {
    rk <- apply(X[, A, drop = FALSE], 2, rank, ties.method = "random")
    cells <- function(q) {
      tot <- round(q * n * na)                      # total number of cells to delete
      k <- rep(tot %/% na, na) + (seq_len(na) <= tot %% na)
      sweep(rk, 2, n - k, ">")                       # the k highest values of each predictor
    }
    inc <- function(q) mean(rowSums(cells(q)) > 0)
    lo <- 0; hi <- 1
    for (it in 1:40) { mid <- (lo + hi) / 2; if (inc(mid) < r) lo <- mid else hi <- mid }
    M[, A] <- cells(hi)
  }
  Xm <- X; Xm[M] <- NA
  Xm
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPUTATION (fitted on the training fold only; response never used)
#   each function returns list(train = , test = ) or, for MICE, a list of m
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
imp_mean <- function(Xtr, Xte) {
  mu <- colMeans(Xtr, na.rm = TRUE)
  f <- function(X) { for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mu[j]; X }
  list(list(train = f(Xtr), test = f(Xte)))
}

imp_knn <- function(Xtr, Xte, k = 5) {
  mu <- colMeans(Xtr, na.rm = TRUE); s <- apply(Xtr, 2, sd, na.rm = TRUE); s[s == 0] <- 1
  Ztr <- sweep(sweep(Xtr, 2, mu), 2, s, "/")
  f <- function(X, self) {
    Z <- sweep(sweep(X, 2, mu), 2, s, "/")
    for (i in which(!stats::complete.cases(X))) {
      obs <- which(!is.na(Z[i, ]))
      d <- sqrt(rowMeans((Ztr[, obs, drop = FALSE] -
                            matrix(Z[i, obs], nrow(Ztr), length(obs), byrow = TRUE))^2,
                         na.rm = TRUE))
      if (self) d[i] <- Inf                       # a row is not its own donor
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

imp_missforest <- function(Xtr, Xte, seed = 1, ntree = 100, maxiter = 5) {
  # Training fold only: missForest is fitted to the training predictors.
  # missForest has no predict step, so the test predictors are imputed with
  # random forests (one per incomplete predictor) fitted to the imputed training
  # data, iterating from a training-mean start (same idea as missForest).
  # The test fold is never used to fit any imputation model.
  set.seed(seed)
  tr <- as.matrix(missForest::missForest(as.data.frame(Xtr), ntree = ntree,
                                         verbose = FALSE)$ximp)
  te <- Xte; Mte <- is.na(Xte)
  cols <- which(colSums(Mte) > 0)
  if (length(cols) > 0) {
    mu <- colMeans(tr)
    for (j in cols) te[Mte[, j], j] <- mu[j]
    rfs <- lapply(cols, function(j)
      randomForest::randomForest(x = tr[, -j, drop = FALSE], y = tr[, j], ntree = ntree))
    names(rfs) <- as.character(cols)
    for (it in seq_len(maxiter)) {
      old <- te
      for (j in cols) te[Mte[, j], j] <-
        predict(rfs[[as.character(j)]], te[Mte[, j], -j, drop = FALSE])
      if (sum((te - old)^2) <= 1e-8 * sum(te^2)) break
    }
  }
  list(list(train = tr, test = te))
}

imp_mice <- function(Xtr, Xte, m = 5, seed = 1) {
  Xall <- as.data.frame(rbind(Xtr, Xte))
  ign <- c(rep(FALSE, nrow(Xtr)), rep(TRUE, nrow(Xte)))   # models use training rows only
  imp <- mice::mice(Xall, m = m, method = "pmm", ignore = ign,
                    printFlag = FALSE, seed = seed)
  lapply(seq_len(m), function(k) {
    out <- as.matrix(mice::complete(imp, k))
    list(train = out[!ign, , drop = FALSE], test = out[ign, , drop = FALSE])
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODELS (same settings as Section 5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Minimum node size: 1% of the training size, but at least 3
# (with minsplit = 1, rpart sets minbucket = 0 and returns a root-only tree).
n_min_fun <- function(n_train) as.integer(max(3, floor(design$prop_min * n_train)))

pr_ctrl <- function(n_train, p, fill, K) {
  PRTree::pr_tree_control(
    max_terminal_nodes = 100, cp = 0.001, by_node = TRUE, tiny_sigma = 0,
    n_min = n_min_fun(n_train),
    n_candidates = as.integer(max(round(p / 2), 3)),
    prop_hold = 1 / (K - 1),     # validation sample ~ one fold (as in Section 5)
    fill_type = fill, proxy_crit = design$proxy, dist = design$dist
  )
}

pred_vec <- function(p) {        # CHECK: format returned by predict.prtree
  if (is.list(p)) p <- if (!is.null(p$yhat)) p$yhat else p[[1]]
  as.numeric(p)
}

fit_pr <- function(Xtr, ytr, Xte, fill, K, seed) {
  set.seed(seed)
  fit <- PRTree::pr_tree(y = ytr, X = Xtr, control = pr_ctrl(nrow(Xtr), ncol(Xtr), fill, K))
  pred_vec(predict(fit, Xte))
}

fit_cart <- function(Xtr, ytr, Xte, seed) {
  ctrl <- rpart::rpart.control(cp = 0.001, xval = 10, maxdepth = 30,
                               minsplit = n_min_fun(nrow(Xtr)))
  set.seed(seed)
  fit <- rpart::rpart(y ~ ., data = data.frame(Xtr, y = ytr), method = "anova", control = ctrl)
  fit <- rpart::prune(fit, cp = fit$cptable[which.min(fit$cptable[, "xerror"]), "CP"])
  as.numeric(predict(fit, newdata = as.data.frame(Xte)))
}

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ONE TASK = dataset x mechanism x rate x replication
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# timing helper: evaluates expr and returns list(value, elapsed seconds)
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- expr
  list(value = v, time = proc.time()[["elapsed"]] - t0)
}

run_task <- function(task) {
  d <- load_complete(task$dataset)
  seed <- design$seed + task$rep * 1000 + match(task$mech, c("MNAR", "MAR", "MCAR")) * 100 +  # fixed order keeps the seeds of saved results
   
    round(task$rate * 100)
  # A = the ceil(p/2) predictors most correlated with y (same for every replication)
  cy <- abs(suppressWarnings(cor(d$X, d$y)))[, 1]; cy[is.na(cy)] <- 0
  A <- sort(order(cy, decreasing = TRUE)[seq_len(ceiling(ncol(d$X) / 2))])
  set.seed(seed)
  Xm <- inject_missing(d$X, task$mech, task$rate, design$pi0, A)
  folds <- sample(rep(seq_len(design$K), length.out = nrow(Xm)))
  res <- list()
  # time_imp: imputation time (s); time_fit: time to fit the tree(s) and predict (s)
  add <- function(fold, method, yte, yhat, time_fit, time_imp = 0) {
    res[[length(res) + 1]] <<- data.frame(
      id = task$id, dataset = task$dataset, mech = task$mech, rate = task$rate,
      rep = task$rep, fold = fold, method = method, rmse = rmse(yte, yhat),
      time_imp = time_imp, time_fit = time_fit,
      pct_inc = mean(!stats::complete.cases(Xm)), version = version)
  }
  for (k in seq_len(design$K)) {
    tr <- folds != k; te <- !tr
    ytr <- d$y[tr]; yte <- d$y[te]; s <- seed + k
    Xtr <- Xm[tr, , drop = FALSE]; Xte <- Xm[te, , drop = FALSE]
    # references on complete data
    f <- timed(fit_pr(d$X[tr, ], ytr, d$X[te, ], 2, design$K, s));  add(k, "PR-full",   yte, f$value, f$time)
    f <- timed(fit_cart(d$X[tr, ], ytr, d$X[te, ], s));              add(k, "CART-full", yte, f$value, f$time)
    # internal strategies and CART with surrogate splits (no imputation)
    for (ft in 0:2) {
      f <- timed(fit_pr(Xtr, ytr, Xte, ft, design$K, s)); add(k, paste0("PR-f", ft), yte, f$value, f$time)
    }
    f <- timed(fit_cart(Xtr, ytr, Xte, s)); add(k, "CART-sur", yte, f$value, f$time)
    # single imputations: one tree each
    for (nm in c("mean", "kNN", "missForest")) {
      im <- timed(switch(nm,
                         mean       = imp_mean(Xtr, Xte),
                         kNN        = imp_knn(Xtr, Xte, design$k_knn),
                         missForest = imp_missforest(Xtr, Xte, s))[[1]])
      D <- im$value
      f <- timed(fit_pr(D$train, ytr, D$test, 2, design$K, s))
      add(k, paste0(nm, "+PR"), yte, f$value, f$time, im$time)
      f <- timed(fit_cart(D$train, ytr, D$test, s))
      add(k, paste0(nm, "+CART"), yte, f$value, f$time, im$time)
    }
    # MICE: m imputations (used by MICE, MICE1 and MICEavg)
    im5 <- timed(imp_mice(Xtr, Xte, design$m_mice, s)); L <- im5$value
    im1 <- timed(imp_mice(Xtr, Xte, 1, s))$time          # cost of a single imputation (timing only)
    F_pr   <- lapply(L, function(D) timed(fit_pr(D$train, ytr, D$test, 2, design$K, s)))
    F_cart <- lapply(L, function(D) timed(fit_cart(D$train, ytr, D$test, s)))
    P_pr   <- sapply(F_pr, `[[`, "value");   T_pr   <- sapply(F_pr, `[[`, "time")
    P_cart <- sapply(F_cart, `[[`, "value"); T_cart <- sapply(F_cart, `[[`, "time")
    add(k, "MICE+PR",    yte, rowMeans(P_pr),   sum(T_pr),   im5$time)   # m trees (standard use)
    add(k, "MICE+CART",  yte, rowMeans(P_cart), sum(T_cart), im5$time)
    add(k, "MICE1+PR",   yte, P_pr[, 1],        T_pr[1],     im1)        # one tree, one imputation
    add(k, "MICE1+CART", yte, P_cart[, 1],      T_cart[1],   im1)
    Davg <- list(train = Reduce(`+`, lapply(L, `[[`, "train")) / length(L),
                 test  = Reduce(`+`, lapply(L, `[[`, "test"))  / length(L))
    f <- timed(fit_pr(Davg$train, ytr, Davg$test, 2, design$K, s))
    add(k, "MICEavg+PR", yte, f$value, f$time, im5$time)
    f <- timed(fit_cart(Davg$train, ytr, Davg$test, s))
    add(k, "MICEavg+CART", yte, f$value, f$time, im5$time)
  }
  dplyr::bind_rows(res)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TASKS, PROGRESS BAR AND PARALLEL EXECUTION (one results file per dataset)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_id <- function(dataset, mech, rate, rep)
  sprintf("%s_%s_r%02d_rep%02d", dataset, mech, round(100 * rate), rep)
res_file <- function(ds) file.path(res_dir, paste0(ds, ".rds"))

grid <- expand.grid(rep = seq_len(design$R), rate = design$rates,
                    mech = design$mechanisms, dataset = design$datasets,
                    stringsAsFactors = FALSE)
tasks <- lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  list(dataset = g$dataset, mech = g$mech, rate = g$rate, rep = g$rep,
       id = make_id(g$dataset, g$mech, g$rate, g$rep))
})

fmt_time <- function(sec) {
  if (!is.finite(sec)) return("--:--:--")
  sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, round(sec %% 60))
}
progress <- function(done, total, t0) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  eta <- if (done > 0) el / done * (total - done) else NA
  w <- 30; k <- round(w * done / total)
  cat(sprintf("\r  [%s%s] %d/%d tasks | elapsed %s | ETA %s   ",
              strrep("=", k), strrep(" ", w - k), done, total, fmt_time(el), fmt_time(eta)))
  flush.console()
}

# Quick test before the full run (uncomment): one task, sequential
# library(PRTree); t0 <- Sys.time(); print(run_task(tasks[[1]])); Sys.time() - t0

# refuse to mix with results from older versions of the study
for (ds in design$datasets) {
  f <- res_file(ds)
  if (file.exists(f) && !identical(unique(readRDS(f)$version), version))
    stop(f, " was produced by an older version of the study. ",
         "Delete the folder ", res_dir, " and run again.")
}

cl <- makeCluster(n_cores)
clusterExport(cl, c("design", "version", "timed", "load_complete", "inject_missing",
                    "imp_mean", "imp_knn", "imp_missforest", "imp_mice", "n_min_fun",
                    "pr_ctrl", "pred_vec", "fit_pr", "fit_cart", "rmse", "run_task"))
invisible(clusterEvalQ(cl, { library(PRTree); library(dplyr); NULL }))
errors <- character(0)

for (ds in design$datasets) {
  f <- res_file(ds)
  res_ds <- if (file.exists(f)) readRDS(f) else NULL
  done_ids <- unique(res_ds$id)
  todo <- Filter(function(t) t$dataset == ds && !(t$id %in% done_ids), tasks)
  n_tot <- sum(vapply(tasks, function(t) t$dataset == ds, logical(1)))
  cat(sprintf("\n%s: %d of %d tasks already done\n", ds, n_tot - length(todo), n_tot))
  if (length(todo) == 0) next
  t0 <- Sys.time(); n_done <- 0
  progress(n_done, length(todo), t0)
  batches <- split(todo, ceiling(seq_along(todo) / n_cores))   # one task per core
  for (b in batches) {
    out <- parLapply(cl, b, function(t) tryCatch(run_task(t), error = function(e)
      paste("ERROR", t$id, conditionMessage(e))))
    ok <- vapply(out, is.data.frame, logical(1))
    errors <- c(errors, unlist(out[!ok]))
    if (any(ok)) {
      res_ds <- dplyr::bind_rows(res_ds, out[ok])
      saveRDS(res_ds, f)                                        # checkpoint
    }
    n_done <- n_done + length(b)
    progress(n_done, length(todo), t0)
  }
  cat("\n")
}
stopCluster(cl)
if (length(errors) > 0) { print(errors); stop("Some tasks failed (see above). Run the script again.") }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARIES (figures and numbers for Section 6), in the style of Section 5
#   Relative errors are computed fold by fold with respect to PRTree fitted to
#   the complete data (PR-full) and to PRTree with fill_type = 2 (PR-f2).
#   Win rate of PR-f2 = share of folds in which PR-f2 has a lower RMSE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
labels <- c(Auto = "Auto MPG", Ozone = "Ozone", Concrete = "Concrete")
d <- dplyr::bind_rows(lapply(res_file(design$datasets), readRDS)) %>%
  mutate(mech = factor(mech, design$mechanisms),
         dataset = factor(labels[dataset], labels[design$datasets]),
         rate_lab = factor(paste0(100 * rate, "%"), paste0(100 * design$rates, "%")))
w <- d %>% select(dataset, mech, rate, rep, fold, method, rmse) %>%
  pivot_wider(names_from = method, values_from = rmse)

cat("\nSanity check - share of folds where CART-full and CART-sur coincide (should be ~0):\n")
print(w %>% group_by(dataset) %>% summarise(identical = mean(abs(`CART-full` - `CART-sur`) < 1e-12)))

# short labels used in the figures (defined in the text of Section 6)
short <- c("CART-full" = "CART-full", "PR-full" = "full", "PR-f0" = "f0", "PR-f1" = "f1", "PR-f2" = "f2",
           "mean+PR" = "mean", "kNN+PR" = "kNN", "missForest+PR" = "mF",
           "MICE+PR" = "MICE", "MICE1+PR" = "MICE1", "MICEavg+PR" = "MICEavg",
           "CART-sur" = "CART-sur", "MICE+CART" = "MICE+CART")

# ---- Ranking of the procedures by win rate -----------------------------------
#   For each procedure: share of pairwise comparisons (same fold) in which it has
#   a lower RMSE than another procedure, over all datasets, mechanisms and rates.
#   The references fitted to the complete data (PR-full, CART-full) are excluded.
cand <- c("PR-f0", "PR-f1", "PR-f2", "CART-sur",
          paste0(c("mean", "kNN", "missForest", "MICE", "MICE1", "MICEavg"), "+PR"),
          paste0(c("mean", "kNN", "missForest", "MICE", "MICE1", "MICEavg"), "+CART"))
xx <- d %>% filter(method %in% cand) %>% select(dataset, mech, rate, rep, fold, method, rmse)
ranking <- xx %>%
  inner_join(xx, by = c("dataset", "mech", "rate", "rep", "fold"), relationship = "many-to-many") %>%
  filter(method.x != method.y) %>%
  group_by(method = method.x) %>% summarise(win = mean(rmse.x < rmse.y), .groups = "drop") %>%
  arrange(desc(win))
cat("\nRanking of the procedures (overall pairwise win rate):\n"); print(ranking, n = 30)
best <- head(ranking$method, 6)
if (!"PR-f2" %in% best) best <- c(head(best, 5), "PR-f2")      # fill_type = 2 is always shown
best <- ranking$method[ranking$method %in% best]                # keep the ranking order
nice <- function(m) sub("\\+PR$", " + PRTree", sub("^PR-f2$", "fill_type = 2", sub("\\+CART$", " + CART", m)))

# loss with respect to the complete data, fold by fold: RMSE / RMSE(PR-full) - 1
# (same scale as the relative error Delta of Section 5)
loss <- d %>% group_by(dataset, mech, rate, rep, fold) %>%
  mutate(E = rmse / rmse[method == "PR-full"] - 1) %>% ungroup()

# ---- Figure: boxplots of the loss of the 6 best procedures --------------------
pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
bx <- loss %>% filter(method %in% best) %>%
  mutate(method = factor(nice(method), nice(best)))
all_cols <- c("PR-f2" = "#CC0000", "MICE+PR" = "#6A3D9A", "missForest+PR" = "#009E73",
              "kNN+PR" = "#0072B2", "MICEavg+PR" = "#E69F00", "MICE1+PR" = "#CC79A7",
              "mean+PR" = "#999999", "MICE+CART" = "#56B4E9", "missForest+CART" = "#A6D96A",
              "kNN+CART" = "#9ECAE1", "MICEavg+CART" = "#FDD49E", "MICE1+CART" = "#F1B6DA",
              "mean+CART" = "#D9D9D9", "CART-sur" = "#444444", "PR-f1" = "#F4A582", "PR-f0" = "#B2182B")
cols <- setNames(all_cols[best], nice(best))
pb <- ggplot(bx, aes(rate_lab, E, fill = method)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.4, linewidth = 0.3,
               position = position_dodge(width = 0.85), width = 0.75) +
  scale_fill_manual(values = cols, name = NULL) +
  facet_grid(dataset ~ mech, scales = "free_y") +
  labs(x = "Proportion of incomplete rows r", y = "RMSE / RMSE(PRTree, complete data) - 1") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = guide_legend(nrow = 1))
ggsave(file.path(res_dir, "semisynthetic_boxplots.pdf"), pb, width = 11, height = 8, device = pdf_dev)

# ---- Table: median (IQR) of the loss of the 6 best procedures -----------------
tab_cols <- c(best[!grepl("CART", best)], best[grepl("CART", best)])
tb <- loss %>% filter(method %in% best) %>%
  group_by(dataset, mech, rate, method) %>%
  summarise(med = median(E), iqr = IQR(E), inc = 100 * mean(pct_inc), .groups = "drop") %>%
  group_by(dataset, mech, rate) %>% mutate(inc = mean(inc), best_m = med == min(med)) %>% ungroup() %>%
  mutate(med = round(med, 2) + 0, iqr = round(iqr, 2) + 0,                 # "+ 0" avoids "-0.00"
         cell = ifelse(best_m, sprintf("\\textbf{%.2f} (%.2f)", med, iqr), sprintf("%.2f (%.2f)", med, iqr)),
         method = factor(method, tab_cols)) %>%
  select(dataset, mech, rate, inc, method, cell) %>%
  pivot_wider(names_from = method, values_from = cell) %>% arrange(dataset, mech, rate)
lines <- c()
for (ds in levels(tb$dataset)) {
  sub <- tb %>% filter(dataset == ds); if (nrow(sub) == 0) next
  lines <- c(lines, sprintf("\\multicolumn{%d}{l}{\\textit{%s}}\\\\", length(tab_cols) + 3, ds))
  for (i in seq_len(nrow(sub))) {
    mech <- if (i == 1 || sub$mech[i] != sub$mech[i - 1]) as.character(sub$mech[i]) else ""
    lines <- c(lines, paste0(paste(c(mech, round(100 * sub$rate[i]), sprintf("%.0f", sub$inc[i]),
                                     unlist(sub[i, tab_cols])), collapse = " & "), "\\\\"))
  }
  lines <- c(lines, "\\midrule")
}
lines[length(lines)] <- "\\bottomrule"
# short column headers; the procedures are grouped under "PRTree" and "CART"
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

# arranges one heatmap per dataset in a 2 x 2 grid with the legend in the free slot
# (same layout as the heatmaps of Section 5)
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

# ---- Figure: win rate of PR-f2 against each competitor ----------------------
comp <- c("PR-f0", "PR-f1", "CART-sur", "mean+PR", "kNN+PR", "missForest+PR",
          "MICE+PR", "MICE1+PR", "MICEavg+PR", "MICE+CART")
hw <- w %>% pivot_longer(all_of(comp), names_to = "competitor", values_to = "rmse_c") %>%
  group_by(dataset, mech, rate, competitor) %>%
  summarise(win = mean(`PR-f2` < rmse_c), .groups = "drop") %>%
  mutate(competitor = factor(short[competitor], rev(short[comp])),
         col = factor(paste0(mech, "\n", 100 * rate, "%"),
                      paste0(rep(design$mechanisms, each = length(design$rates)), "\n",
                             100 * design$rates, "%")))
ph <- grid_heatmaps(hw, function(x) ggplot(x, aes(col, competitor, fill = win)) +
                      geom_tile(colour = "white") +
                      geom_text(aes(label = sprintf("%.0f%%", 100 * win)), size = 4.5) +
                      scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2",
                                           midpoint = 0.5, limits = c(0, 1)),
                    "f2 win rate")
ggsave(file.path(res_dir, "semisynthetic_heatmaps.pdf"), ph, width = 15, height = 10, device = pdf_dev)

# ---- numbers quoted in the text ----------------------------------------------
cat("\nDataset sizes:\n")
print(bind_rows(lapply(design$datasets, function(n) { x <- load_complete(n)
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

# ---- computational cost --------------------------------------------------------
#   total time per fold = imputation + fitting/prediction (seconds).
#   Tasks run in parallel (n_cores at a time), so absolute times are inflated by
#   competition for the machine; comparisons between procedures remain valid.
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
             apply(tt, 1, function(r) paste0(tlab[r[["method"]]], " & ",
                                             paste(r[-1], collapse = " & "), "\\\\")),
             "\\bottomrule", "\\end{tabular}"),
           file.path(res_dir, "tab_semisynthetic_time.tex"))

# ---- Figure: median time per fold of the 6 best procedures -------------------
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
                      geom_text(aes(label = ifelse(t < 10, sprintf("%.2f", t), sprintf("%.1f", t))),
                                size = 4.5) +
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
