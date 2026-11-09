# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT 6b: RELATIVE ERROR BOXPLOTS AND WIN-RATE HEATMAPS (FIGURES 1 & 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loads consolidated experimental results from `results/prtree_cart_raw.rds`
# and computes the relative error of PRTree against CART:
#
#   delta = (RMSE_PRTree - RMSE_CART) / RMSE_CART
#
# Generates and saves:
#   1. delta_boxplots.pdf : Boxplots of relative error by fill_type for each dataset (Figure 1 in paper)
#   2. heatmaps.pdf       : Grid of 2D heatmaps showing PRTree win rates (delta < 0) (Figure 2 in paper)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. LOAD AND PREPARE RESULTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results <- readRDS("results/prtree_cart_raw.rds")
results$fill <- as.character(results$fill)

# Map long metadata names to standardized short names
dataset_map <- c(
  "Productivity Prediction of Garment Employees" = "Garment",
  "Communities and Crime"                        = "Crime"
)

results <- results %>%
  mutate(dataset = recode(dataset, !!!dataset_map))

# Baseline CART reference per dataset and fold repetition
cart_ref <- results %>%
  filter(method == "CART") %>%
  select(dataset, rep, cart_rmse = rmse)

# PRTree relative error (delta) and formatted distribution labels
plot_data <- results %>%
  filter(method == "PRTree") %>%
  left_join(cart_ref, by = c("dataset", "rep")) %>%
  mutate(
    delta = rmse / cart_rmse - 1,
    dist_label = case_when(
      dist == "norm"  ~ "N",
      dist == "t"     ~ paste0("t(", par, ")"),
      dist == "lnorm" ~ paste0("LN(", par, ")"),
      dist == "gamma" ~ paste0("G(", par, ")")
    ),
    dist_order = case_when(
      dist == "norm"  ~ 1,
      dist == "t"     ~ 2,
      dist == "lnorm" ~ 3,
      dist == "gamma" ~ 4
    )
  )

# Order distribution factor levels
dist_levels <- plot_data %>%
  arrange(dist_order, par) %>%
  distinct(dist_label) %>%
  pull(dist_label)

plot_data$dist_label <- factor(plot_data$dist_label, levels = dist_levels)

# Consistent fill_type color scheme (Okabe-Ito friendly)
fill_cols <- c(
  "0" = "#0072B2",  # blue
  "1" = "#009E73",  # bluish green
  "2" = "#D55E00"   # vermillion
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. FIGURE 1: RELATIVE ERROR BOXPLOTS BY fill_type (delta_boxplots.pdf)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ymax <- max(plot_data$delta[plot_data$delta <= 0.6])

# Annotate points trimmed above ymax
above_lim <- plot_data %>%
  group_by(dataset, fill) %>%
  summarise(n_out = sum(delta > ymax), .groups = "drop") %>%
  filter(n_out > 0) %>%
  mutate(
    y     = ymax - 0.05,
    label = paste0(n_out, "\u25B2"),
    xpos  = c("0" = 1.15, "1" = 2.15, "2" = 3.15)[fill]
  )

p_box <- ggplot(plot_data, aes(x = fill, y = delta, fill = fill)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5, linewidth = 0.3) +
  geom_text(
    data        = above_lim,
    aes(x = xpos, y = y, label = label),
    inherit.aes = FALSE,
    size        = 5,
    col         = "red"
  ) +
  geom_hline(
    yintercept = seq(-0.6, 0.6, 0.1),
    linetype   = "dashed",
    color      = "grey40",
    linewidth  = 0.3
  ) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  scale_fill_manual(values = fill_cols, name = "fill_type") +
  facet_grid(~ dataset, scales = "free_y") +
  coord_cartesian(ylim = c(min(plot_data$delta), ymax)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    legend.position   = "none",
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold")
  )

ggsave(
  "delta_boxplots.pdf",
  plot   = p_box,
  device = cairo_pdf,
  width  = 15,
  height = 5
)
cat("Saved: delta_boxplots.pdf\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. FIGURE 2: WIN-RATE HEATMAPS (heatmaps.pdf)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute win rate per configuration: proportion of folds where delta < 0
heat_tbl <- plot_data %>%
  mutate(
    row_lab = paste0(
      "f", fill, "-",
      c("mean" = "M", "var" = "V", "both" = "B")[proxy]
    ),
    col_lab = dist_label,
    win     = delta < 0
  ) %>%
  group_by(dataset, row_lab, col_lab) %>%
  summarise(
    win_rate  = mean(win),
    med_delta = median(delta),
    .groups   = "drop"
  )

# Build individual heatmap panels
datasets_ordered <- unique(heat_tbl$dataset)
heatmaps <- lapply(datasets_ordered, function(ds) {
  dat <- heat_tbl %>% filter(dataset == ds)

  ggplot(dat, aes(x = col_lab, y = row_lab, fill = win_rate)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = sprintf("%.0f%%", 100 * win_rate)),
      size = 5
    ) +
    scale_fill_gradient2(
      low      = "#D55E00",
      mid      = "white",
      high     = "#0072B2",
      midpoint = 0.5,
      limits   = c(0, 1),
      name     = "PRTree\nwin rate"
    ) +
    labs(title = ds, x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 1),
      plot.title  = element_text(face = "bold", size = 14)
    )
})

# Extract common horizontal colorbar legend
p_leg <- heatmaps[[1]] +
  guides(
    fill = guide_colorbar(
      direction      = "horizontal",
      title.position = "top",
      barwidth       = unit(14, "cm"),
      barheight      = unit(1.5, "cm")
    )
  )
shared_legend <- cowplot::get_legend(p_leg)

# Remove legends from individual plots for grid arrangement
heatmaps_noleg <- lapply(heatmaps, function(p) {
  p + theme(legend.position = "none")
})

# Arrange panels into 2 columns x 3 rows grid (with legend in 6th slot)
h_grid <- ggarrange(
  heatmaps_noleg[[1]],
  heatmaps_noleg[[2]],
  heatmaps_noleg[[3]],
  heatmaps_noleg[[4]],
  heatmaps_noleg[[5]],
  ggpubr::as_ggplot(shared_legend),
  ncol = 2,
  nrow = 3
)

ggsave(
  "heatmaps.pdf",
  plot   = h_grid,
  device = cairo_pdf,
  width  = 15,
  height = 15
)
cat("Saved: heatmaps.pdf\n")
