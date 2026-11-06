library(dplyr)
library(purrr)
library(ggplot2)

# ============================================================
# READ ALL RESULTS
# ============================================================

# files <- list.files(
#   "results",
#   pattern = "\\.rds$",
#   full.names = TRUE
# )
# results <- map_dfr(files, readRDS)

results <- readRDS("results/prtree_cart_raw.rds")
results$fill <- as.character(results$fill)
dataset_map <- c(
  "Productivity Prediction of Garment Employees" = "Garment",
  "Communities and Crime" = "Crime"
)

results <- results %>%
  mutate(
    dataset = recode(dataset, !!!dataset_map)
  )


# ============================================================
# CART REFERENCE + RELATIVE ERROR
# ============================================================

cart_ref <- results %>%
  filter(method == "CART") %>%
  select(dataset, rep, cart_rmse = rmse)

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

# ============================================================
# ORDER DISTRIBUTIONS
# ============================================================

dist_levels <- plot_data %>%
  arrange(dist_order, par) %>%
  distinct(dist_label) %>%
  pull(dist_label)

plot_data$dist_label <- factor(
  plot_data$dist_label,
  levels = dist_levels
)

# ============================================================
# COLORS = fill_type
# ============================================================

fill_cols <- c(
  "0" = "#0072B2",  # blue
  "1" = "#009E73",  # bluish green
  "2" = "#D55E00"   # vermillion
)


# ============================================================
# BOXPLOT 
# ============================================================

ymax <- max(plot_data$delta[plot_data$delta < 1.1])

fill_offset <- c("0" = -0.27,"1" = 0.00,"2" = 0.27)

above_lim <- plot_data %>%
  group_by(dataset, dist_label, fill) %>%
  summarise(n_out = sum(delta > ymax), .groups = "drop") %>%
  filter(n_out > 0) %>%
  mutate(
    y = ymax + 0.05,
    label = paste0(n_out, "\u25B2"),
    xpos = as.numeric(factor(dist_label, levels = levels(plot_data$dist_label))) +
      fill_offset[fill]
  )

ggplot(
  plot_data,
  aes(x = dist_label, y = delta, fill = fill)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5, linewidth = 0.3) +
  geom_text(
    data = above_lim,
    aes(x = xpos, y = y, label = label, color = fill),
    inherit.aes = FALSE, size = 4
  ) +
  scale_fill_manual(values = fill_cols, name = "fill_type") +
  scale_color_manual(values = fill_cols,guide = "none") +
  facet_grid(dataset ~ ., scales = "free_x") +
  coord_cartesian(ylim = c(min(plot_data$delta),1.1 * ymax)) +
  labs(x = NULL,y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 1),
    strip.text = element_text(face = "bold")
  )


# ============================================================
# BOXPLOT BY fill_type
# ============================================================

ymax <- max(plot_data$delta[plot_data$delta <= 0.6])

above_lim <- plot_data %>%
  group_by(dataset, fill) %>%
  summarise(n_out = sum(delta > ymax), .groups = "drop") %>%
  filter(n_out > 0) %>%
  mutate(y = ymax - 0.05, 
         label = paste0(n_out, "\u25B2"),
         xpos = c("0" = 1.15, "1" = 2.15, "2" = 3.15)[fill])

p <- ggplot(
  plot_data,
  aes(x = fill, y = delta, fill = fill)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5, linewidth = 0.3) +
  geom_text(
    data = above_lim,
    aes(x = xpos, y = y, label = label),
    inherit.aes = FALSE, size = 5, col = "red") + 
  geom_hline(
    yintercept = seq(-0.6, 0.6, 0.1),
    linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  scale_fill_manual(values = fill_cols,name = "fill_type") +
  facet_grid( ~ dataset, scales = "free_y") +
  coord_cartesian(ylim = c(min(plot_data$delta), ymax))  +
  labs(x = NULL,y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
p

ggsave(
  "delta_boxplots.pdf",
  plot = p,
  device = cairo_pdf,
  width = 15,
  height = 5
)

# ============================================================
# WIN RATE TABLE
# ============================================================

heat_tbl <- plot_data %>%
  mutate(
    # row_lab = paste0(
    #   "fill=",
    #   fill,
    #   "\nproxy=",
    #   proxy
    # ),
    row_lab = paste0(
      "f", fill, "-",
      c(
        "mean" = "M",
        "var" = "V",
        "both" = "B"
      )[proxy]
    ),
    col_lab = dist_label,
    win = delta < 0
  ) %>%
  group_by(
    dataset,
    row_lab,
    col_lab
  ) %>%
  summarise(
    win_rate = mean(win),
    med_delta = median(delta),
    .groups = "drop"
  )

# ============================================================
# WIN RATE HEATMAPS
# ============================================================

heatmaps <- lapply(
  unique(heat_tbl$dataset),
  function(ds) {
    
    dat <- heat_tbl %>%
      filter(dataset == ds)
    
    ggplot(
      dat,
      aes(
        col_lab,
        row_lab,
        fill = win_rate
      )
    ) +
      geom_tile(color = "white") +
      geom_text(
        aes(
          label = sprintf(
            "%.0f%%",
            100 * win_rate
          )
        ),
        size = 5
      ) +
      scale_fill_gradient2(
        low = "#D55E00",
        mid = "white",
        high = "#0072B2",
        midpoint = 0.5,
        limits = c(0, 1),
        name = "PRTree\nwin rate"
      ) +
      labs(
        title = ds,
        x = NULL,
        y = NULL
      ) +
      theme_bw(base_size = 16) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 0,
          hjust = 1
        ),
        plot.title = element_text(
          face = "bold",
          size = 14
        )
      )
  }
)

heatmaps[[1]]
heatmaps[[2]]
heatmaps[[3]]
heatmaps[[4]]
heatmaps[[5]]


h <- ggplot(
  heat_tbl,
  aes(
    col_lab,
    row_lab,
    fill = win_rate
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = sprintf("%.0f%%", 100 * win_rate)
    ),
    size = 4.5
  ) +
  scale_fill_gradient2(
    low = "#D55E00",
    mid = "white",
    high = "#0072B2",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "PRTree\nwin rate"
  ) +
  facet_wrap(
    dataset ~ .
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold"
    )
  )
h

library(cowplot)

p_leg <- heatmaps[[1]] +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(14, "cm"),
      barheight = unit(1.5, "cm")
    )
  ) 

leg <- get_legend(p_leg)

heatmaps_noleg <- lapply(
  heatmaps,
  function(p) p + theme(legend.position = "none")
)

h <- ggarrange(
  heatmaps_noleg[[1]],
  heatmaps_noleg[[2]],
  heatmaps_noleg[[3]],
  heatmaps_noleg[[4]],
  heatmaps_noleg[[5]],
  ggpubr::as_ggplot(leg),
  ncol = 2,
  nrow = 3
)

ggsave(
  "heatmaps.pdf",
  plot = h,
  device = cairo_pdf,
  width = 15,
  height = 15
)
