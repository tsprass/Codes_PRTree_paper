library(dplyr)
library(ggplot2)
library(forcats)
library(readr)
library(purrr)
library(knitr)


source("03 - functions.R")

# # Metadata: run 02 - Filter to get statistics for metadata 
# 
# # Processed data: 
# data_file <- "data/processed_datasets.rds"
# cat("📁 Loading pre-processed data from", data_file, "\n")
# datasets <- readRDS(data_file)
# 
# # Quick summary (datasets) ----
# summary_df <- purrr::map_df(names(datasets), function(name) {
#   d <- datasets[[name]]
#   data.frame(
#     Dataset = name,
#     Response = d$response,
#     n = d$n,
#     p = d$p,
#     has_missing = d$has_natural_missing,
#     missing_pct = round(100 * sum(is.na(d$X)) / (d$n * d$p), 2),
#     stringsAsFactors = FALSE
#   )
# }) %>% arrange(desc(n))
# 
# print(summary_df)



# ============================================================
# LOAD FINAL DATASETS
# ============================================================

individual_files <- list.files("data/individual", pattern = "\\.rds$")
processed_datasets <- gsub("\\.rds$", "", individual_files)

datasets <- lapply(
  processed_datasets, 
  function(name) readRDS(paste0("data/individual/", name, ".rds")))
names(datasets) <- processed_datasets

# ============================================================
# KEEP ONLY DATASETS USED IN THE PAPER
# ============================================================
keep <- c(
  "Crime",
  "Garment",
  "Auto",
  "Automobile",
  #"BoneMarrow",
  #"Tennis",
  "Ozone"
)
datasets <- datasets[keep]

# ============================================================
# SOURCE INFORMATION
# ============================================================
source_info <- tibble::tribble(
  ~internal_name, ~source,     ~id,
  "Crime",         "UCI",       "183",
  "Garment",      "UCI",       "597",
  "Auto",         "UCI",       "9",
  "Automobile",   "UCI",       "10",
 # "BoneMarrow",   "UCI",       "565",
#  "Tennis",       "UCI",       "300",
  "Ozone",        "mlbench",   "-"
)

# ============================================================
# BUILD TABLE
# ============================================================

tab1 <- map_df(names(datasets), function(name) {
  
  d <- datasets[[name]]
  src <- source_info %>% filter(internal_name == name)
  
  data.frame(
    Dataset = d$name,
    Source = src$source,
    ID = src$id,
    Response = d$response,
    `Response range` = paste0(
      "$[",
      round(min(d$y, na.rm = TRUE), 2),
      ", ",
      round(max(d$y, na.rm = TRUE), 2),
      "]$"
    ),
    n = nrow(d$X),
    p = ncol(d$X),
    `missing (global)` = round(
      100 * mean(is.na(as.matrix(d$X))),
      2
    ),
    `n rows` = sum(
        apply(is.na(d$X), 1, any)
      ),
    stringsAsFactors = FALSE
  )
})

# ============================================================
# ORDER
# ============================================================

tab1 <- tab1 %>% arrange(desc(n))

# ============================================================
# PRINT LATEX TABLE
# ============================================================

kable(
  tab1,
  format = "latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = "",
  caption = "Summary of the datasets used in the empirical evaluation.",
  label = "tab:datasets"
)


tab1 <- map_df(names(datasets), function(name) {
  
  d <- datasets[[name]]
  src <- source_info %>% filter(internal_name == name)
  
  data.frame(
    Dataset = d$name,
    Source = src$source,
    ID = src$id,
    `Response range` = paste0(
      "$[",
      round(min(d$y, na.rm = TRUE), 2),
      ", ",
      round(max(d$y, na.rm = TRUE), 2),
      "]$"
    ),
    n = nrow(d$X),
    p = ncol(d$X),
    m_global = round(
      100 * mean(is.na(as.matrix(d$X))),
      2
    ),
    n_miss = sum(
      apply(is.na(d$X), 1, any)
    ),
    stringsAsFactors = FALSE
  )
}) %>%
  arrange(desc(n))

kable(
  tab1,
  format = "latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = "",
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
  caption = paste(
    "Summary of the datasets used in the empirical evaluation:",
    "$n$ and $p$ denote the dataset dimensions;",
    "$m_{global}$ is the percentage of missing predictor cells among all $n\times p$ predictor values;",
    "$n_{miss}$ is the number of observations containing at least one missing predictor value."
  ),
  label = "tab:datasets"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PAPER FIGURE - BOXPLOTS BY DATASET
# PRTree vs CART under correlated missingness
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ---------------------------------------------------------------------------
# INPUT / OUTPUT
# ---------------------------------------------------------------------------
result_files <- "results/prtree_cart_raw.rds"
output_dir   <- "paper_outputs"

if (!dir.exists(output_dir)) dir.create(output_dir)

# ---------------------------------------------------------------------------
# LOAD RESULTS
# ---------------------------------------------------------------------------
results <- readRDS(result_files)

# ---------------------------------------------------------------------------
# CLEAN / STANDARDIZE
# ---------------------------------------------------------------------------
results <- results |>
  mutate(
    rho = as.numeric(rho),
    fill = as.character(fill),
    method_label = case_when(
      method == "CART" ~ "CART",
      method == "PRTree" & fill == "0" ~ "PRTree (fill = 0)",
      method == "PRTree" & fill == "1" ~ "PRTree (fill = 1)",
      method == "PRTree" & fill == "2" ~ "PRTree (fill = 2)",
      TRUE ~ method
    ),
    method_label = factor(
      method_label,
      levels = c("CART", "PRTree (fill = 0)", "PRTree (fill = 1)", "PRTree (fill = 2)")
    ),
    rho_f = factor(rho, levels = sort(unique(rho)))
  )

# ---------------------------------------------------------------------------
# NORMALIZATION
# NRMSE_d(m, rho, k) = RMSE_d(m, rho, k) / min_{m',rho',k'} RMSE_d(m',rho',k')
# ---------------------------------------------------------------------------
dataset_best <- results |>
  group_by(dataset) |>
  summarise(
    best_rmse = min(rmse, na.rm = TRUE),
    .groups = "drop"
  )

results_norm <- results |>
  left_join(dataset_best, by = "dataset") |>
  mutate(
    nrmse = rmse / best_rmse
  )

# ---------------------------------------------------------------------------
# VISUAL SETTINGS
# ---------------------------------------------------------------------------
method_colors <- c(
  "CART" = "#000000",
  "PRTree (fill = 0)" = "#D55E00",
  "PRTree (fill = 1)" = "#0072B2",
  "PRTree (fill = 2)" = "#CC79A7"
)

method_shapes <- c(
  "CART" = 16,
  "PRTree (fill = 0)" = 17,
  "PRTree (fill = 1)" = 15,
  "PRTree (fill = 2)" = 18
)

# ---------------------------------------------------------------------------
# MAIN FIGURE: BOXPLOTS BY DATASET
# ---------------------------------------------------------------------------
fig_box_dataset <- ggplot(
  results_norm,
  aes(
    x = rho_f,
    y = nrmse,
    fill = method_label,
    color = method_label
  )
) +
  #geom_hline(yintercept = 1, linetype = 2, linewidth = 0.45, color = "gray45") +
  geom_boxplot(
    width = 0.72,
    alpha = 0.22,
    outlier.shape = NA,
    position = position_dodge(width = 0.82),
    linewidth = 0.50
  ) +
  geom_point(
    aes(shape = method_label),
    size = 1.2,
    alpha = 0.35,
    position = position_jitterdodge(
      jitter.width = 0.08,
      jitter.height = 0,
      dodge.width = 0.82
    )
  ) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(
    x = expression(rho),
    y = "Normalized RMSE",
    fill = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

fig_box_dataset

# ---------------------------------------------------------------------------
# SAVE DATA USED IN THE FIGURE
# ---------------------------------------------------------------------------
#write_csv(results_norm, file.path(output_dir, "results_normalized_fold_level.csv"))

cat("\nDone. Figure saved to:", file.path(output_dir, "fig_boxplot_by_dataset.pdf"), "\n")
