# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HELPER FUNCTIONS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)

# Safely fetch datasets from UCI repository ----
fetch_safe <- function(id, response_name = NULL) {
  tryCatch({
    data <- ucimlrepo::fetch_ucirepo(id = id)

    # Get original data
    df <- data$data$original

    # If no original, try features + targets
    if (is.null(df) && !is.null(data$data$features)) {
      df <- data$data$features
      if (!is.null(data$data$targets)) {
        df <- bind_cols(df, data$data$targets)
      }
    }

    if (is.null(df)) {
      return(NULL)
    }

    # Convert UCI missing markers
    df[df == "?"] <- NA

    # Recover numeric columns
    df[] <- lapply(df, function(x) {
      type.convert(x, as.is = TRUE)
    })

    # Determine response variable
    if (is.null(response_name)) {

      # Try to find target from metadata
      if (!is.null(data$variables)) {
        target <- data$variables$name[data$variables$role == "Target"]
        if (length(target) > 0) response_name <- target[1]
      }

      # If still NULL, flag the response as not found
      if (is.null(response_name)) {
        response_name <- "Not found"
      }
    }

    list(
      data = df,
      response = response_name,
      name = data$metadata$name,
      n = nrow(df),
      p = ncol(df) - 1,
      id = id
    )

  },
  error = function(e) {
    warning("Failed to fetch UCI dataset id = ", id, ": ", conditionMessage(e),
            call. = FALSE, immediate. = TRUE)
    NULL
  }
  )
}


# Prepare data for analysis ----
# 1. Remove rows with NA in response
# 2. Separate X and y
# 3. Remove categorical variables (no conversion)
# 4. Check for natural missing in predictors
# 5. Convert to matrix format
prepare_data <- function(d) {
  # Stop early if the response variable is not a column of the data
  if (!d$response %in% names(d$data)) {
    stop("Response variable '", d$response, "' not found in dataset '", d$name, "'")
  }

  # Remove rows with NA in response
  df <- d$data[!is.na(d$data[[d$response]]), ]

  # Separate X and y
  y <- df[[d$response]]
  X <- df[, !names(df) %in% d$response, drop = FALSE] %>%
    select(where(is.numeric))

  # Check for missing values in predictors
  has_natural_missing <- any(is.na(X))

  # Convert to matrix format
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.vector(y)
  storage.mode(y) <- "double"

  if (any(is.na(y))) {
    warning(paste0("fail to process y for ", d$name, "\n"), immediate. = TRUE)
  }

  list(
    X = X,
    y = y,
    name = d$name,
    response = d$response,
    n = nrow(X),
    p = ncol(X),
    has_natural_missing = has_natural_missing
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CROSS-VALIDATION FUNCTIONS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute the minimum leaf size
nmin <- function(prop = 0.01, n_train = 0, fixed = FALSE, fixed.size = 5){
  if (fixed) return(as.integer(fixed.size))
  as.integer(max(1, floor(prop * n_train)))
}

# `config` is the object created by PRTree::pr_tree_control_cv() (it accepts any
# argument of PRTree::pr_tree_control(), plus method, n_rep and stratify).
# The scripts and functions below add or use these fields:
#   - n_rep, stratify : number of folds and stratification (read by get_folds())
#   - fold_idx        : precomputed fold indices; created by get_folds() if NULL
#   - maxdepth, xval, minsplit, cp, ... : matching rpart.control() arguments,
#                       copied into the CART control by run_cart_cv()
#   - fill_type, proxy_crit, dist, df, shape, sdlog : overwritten by
#                       run_dataset_experiment() for each PRTree configuration

get_folds <- function(config, seed, X){
  # set seed to ensure the same results as in PRTree
  set.seed(seed)
  if (config$stratify) {
    is_NA <- apply(X, 1, anyNA)
    fold_idx <- PRTree:::.sample.stratified_folds(
      n_obs = nrow(X), n_folds = config$n_rep, is_NA = is_NA
    )
  } else {
    fold_idx <- sample(rep(1:config$n_rep, length.out = nrow(X)))
  }
  return(fold_idx)
}


# Run CART cross-validation ----
#  - returns VECTOR of RMSEs
run_cart_cv <- function(data, config, seed) {

  fold_idx <- config$fold_idx
  n_folds <- config$n_rep
  test_rmse <- numeric(n_folds)

  # Create CART control with defaults, then update with whatever is in config
  cart_control <- rpart::rpart.control()
  config_params <- intersect(names(cart_control), names(config))
  if (length(config_params) > 0) {
    cart_control[config_params] <- config[config_params]
  }
  cart_control$maxdepth <- min(30, cart_control$maxdepth)

  # Run CV
  for (fold in 1:n_folds) {
    test_idx <- which(fold_idx == fold)
    train_idx <- which(fold_idx != fold)

    X_train <- data$X[train_idx, , drop = FALSE]
    y_train <- data$y[train_idx]
    X_test <- data$X[test_idx, , drop = FALSE]
    y_test <- data$y[test_idx]

    train_df <- as.data.frame(X_train)
    train_df$y <- y_train

    # set seed again to ensure the same results
    # when running again (CART uses internal CV)
    set.seed(seed + fold)
    fit <- rpart::rpart(
      y ~ .,
      data = train_df,
      method = "anova",
      control = cart_control
    )
    # prune the tree
    fit <- rpart::prune(fit, cp = fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])

    test_df <- as.data.frame(X_test)
    pred <- predict(fit, newdata = test_df)

    test_rmse[fold] <- sqrt(mean((y_test - pred)^2))
  }

  return(test_rmse)
}

# Run PRTree cross-validation  ----
#  - returns VECTOR of RMSEs
run_pr_cv <- function(data, config) {
  cv <- PRTree::pr_tree_cv(data$y, data$X, control_cv = config, verbose = FALSE)
  return(cv$rmse_by_rep$test)
}

# Run complete experiment for one dataset ----
#  - returns RAW values
run_dataset_experiment <- function(
    data, config, verbose = FALSE, seed = 1234,
    fill_options = c(0,1,2),
    proxy_options = c("mean", "var", "both"),
    dist_grid = NULL) {

  # print the current status
  if (verbose) {
    cat("\n", strrep("=", 60), "\n")
    cat("Processing:", data$name, "\n")
    cat("  n =", data$n, ", p =", data$p, "\n")
    cat(strrep("=", 60), "\n")
  }

  # create a list to return the results
  results_list <- list()

  # create a dummy grid so the code will not break
  if (is.null(dist_grid)) {
    dist_grid <- tibble::tribble(
      ~dist,    ~df, ~shape, ~sdlog,
      "norm",   NA,  NA,     NA
    )
  }

  if (is.null(config$fold_idx)) {
    config$fold_idx <- get_folds(config = config, seed = seed, X = data$X)
  }

  # run CART
  if (verbose) cat("[", Sys.time(), "] Running CART...\n")
  cart_values <- run_cart_cv(data, config = config, seed = seed)
  for (rep in seq_along(cart_values)) {
    results_list[[paste0("CART_", rep)]] <- data.frame(
      method = "CART",
      fill = "--",
      proxy = "--",
      dist = "--",
      par = 0,
      rep = rep,
      rmse = cart_values[rep],
      stringsAsFactors = FALSE
    )
  }

  # run PRTree
  # loop1: fill_type
  for (fill in fill_options) {
    config$fill_type <- fill

    # loop2: proxy criterion
    for (proxy in proxy_options) {
      config$proxy_crit <- proxy

      # loop3: distribution
      for (i in seq_len(nrow(dist_grid))) {
        dist_cfg <- dist_grid[i, ]
        config$dist <- dist_cfg$dist
        config$df <- if (!is.na(dist_cfg$df)) dist_cfg$df else NULL
        config$shape <- if (!is.na(dist_cfg$shape)) dist_cfg$shape else NULL
        config$sdlog <- if (!is.na(dist_cfg$sdlog)) dist_cfg$sdlog else NULL

        if (verbose) {
          cat(paste0("[", Sys.time(), "] Running PRTree (fill = ",
          fill, ", proxy = ", proxy, ", dist = ", config$dist, ")...\n")
          )
        }
        set.seed(seed)
        pr_values <- run_pr_cv(data = data, config = config)

        for (rep in seq_along(pr_values)) {
          par <- switch(config$dist,
            "norm" = 0,
            "t" = config$df,
            "lnorm" = config$sdlog,
            "gamma" = config$shape,
            stop("Unknown distribution: ", config$dist)
          )
          key <- paste0("PRTree_", fill, "_", proxy, "_", config$dist, "_", par, "_",  rep)
          results_list[[key]] <- data.frame(
            method = "PRTree",
            fill = as.character(fill),
            proxy = proxy,
            dist = config$dist,
            par = par,
            rep = rep,
            rmse = pr_values[rep],
            stringsAsFactors = FALSE
          )
        }
      } # loop3
    } # loop2
  } # loop1

  bind_rows(results_list) %>%
    mutate(dataset = data$name, .before = 1)
}
