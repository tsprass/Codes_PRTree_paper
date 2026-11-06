library(dplyr)
library(ggplot2)
library(stringr)

# ============================================================
# LOAD + FILTER
# ============================================================

results_base <- readRDS("results/prtree_cart_raw.rds")
#rem <- c("Tennis Major Tournament Match Statistics", "Bone marrow transplant: children")
results <- results_base #%>% 
  # filter(!dataset %in% rem,
  #    !(dist == "lnorm" & (par == 0.1 | par == 1)),
  #    !(dist == "gamma" & (par == 5 | par == 1)),
  #    !(dist == "t" & par == 1))


results$fill <- as.character(results$fill)
results$rmse <- results$rmse

# ============================================================
# CART REFERENCE + DELTA
# ============================================================

cart_ref <- results %>% filter(method == "CART") %>%
  select(dataset, rep, cart_rmse = rmse)

results <- results %>% left_join(cart_ref, by = c("dataset", "rep")) %>%
  mutate(delta = rmse / cart_rmse - 1, 
         prtree_win = delta < 0)


worst <- results %>% filter(delta > 5)

# ============================================================
# LABELS + FACTORS
# ============================================================

results <- results %>%
  mutate(
    method_group = case_when(
      method == "CART" ~ "CART", fill == "0" ~ "fill_type = 0", 
      fill == "1" ~ "fill_type = 1", fill == "2" ~ "fill_type = 2"
    ),
    dist_label = case_when(
      dist == "norm"  ~ "norm",
      dist == "t"     ~ paste0("t (df = ", par, ")"),
      dist == "lnorm" ~ paste0("lnorm\n(sd = ", par, ")"),
      dist == "gamma" ~ paste0("gamma\n(shape = ", par, ")")
    ),
    dist_order = case_when(
      dist == "norm"  ~ 1, dist == "t"     ~ 2,
      dist == "lnorm" ~ 3, dist == "gamma" ~ 4
    ),
    proxy_order = case_when(
      proxy == "mean" ~ 1, proxy == "var"  ~ 2, proxy == "both" ~ 3, 
      TRUE ~ 99
    ),
    config_key = ifelse(
      method == "CART", "CART", paste0(method_group, "__", dist_label, "__", proxy)),
    config_label = ifelse(method == "CART", "CART", dist_label)
  )

plot_data <- results %>% filter(method != "CART") %>%
  arrange(fill, dist_order, par, proxy_order)
config_levels <- plot_data %>% distinct(config_key) %>% pull(config_key)
plot_data$config_key <- factor(plot_data$config_key, levels = config_levels)
label_tbl <- plot_data %>% distinct(config_key, config_label)
label_map <- setNames(label_tbl$config_label, label_tbl$config_key)

# ============================================================
# COLORS
# ============================================================

make_pal <- function(keys, cols) setNames(colorRampPalette(cols)(length(keys)),keys)

fill0_keys <- config_levels[str_detect(config_levels, "fill_type = 0")]
fill1_keys <- config_levels[str_detect(config_levels, "fill_type = 1")]
fill2_keys <- config_levels[str_detect(config_levels, "fill_type = 2")]
color_map <- c(
  make_pal(fill0_keys, c("#c6dbef", "#08306b")),
  make_pal(fill1_keys, c("#c7e9c0", "#00441b")),
  make_pal(fill2_keys, c("#fcbba1", "#67000d"))
)

dist_levels <- plot_data %>% arrange(dist_order, par) %>%
  distinct(dist_label) %>% pull(dist_label)
plot_data$dist_label <- factor(plot_data$dist_label, levels = dist_levels)

# ============================================================
# DELTA BOXPLOTS
# ============================================================

delta_boxplots <- lapply(
  unique(plot_data$dataset),
  function(ds) {
    dat <- plot_data %>% filter(dataset == ds, delta < 5.95)
    ggplot(dat, aes(x = dist_label, y = delta, fill = config_key)) +
      geom_boxplot(alpha = 0.95, linewidth = 0.25, outlier.size = 2) +
      geom_hline(yintercept = 0, linewidth = 0.4, color = "red") +
      scale_fill_manual(values = color_map) +
      facet_grid(method_group ~ proxy, scales = "free_x", space = "free_x") +
      labs(title = ds, x = NULL, 
           y = expression(frac(RMSE[PRTree] - RMSE[CART], RMSE[CART]))) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
            plot.title = element_text(face = "bold", size = 14 ),
            panel.grid.minor = element_blank())
  }
)

delta_boxplots[[1]]
delta_boxplots[[2]]
delta_boxplots[[3]]
delta_boxplots[[4]]
delta_boxplots[[5]]


dat <- plot_data %>% filter(delta < 1)
ggplot(dat, aes(x = dist_label, y = delta, fill = config_key)) +
  geom_boxplot(alpha = 0.95, linewidth = 0.25, outlier.size = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "red") +
  scale_fill_manual(values = color_map) +
  facet_grid(dataset ~ method_group, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        plot.title = element_text(face = "bold", size = 14 ),
        panel.grid.minor = element_blank())

ggplot(plot_data, aes(x = dist_label, y = delta, fill = config_key)) +
  geom_boxplot(alpha = 0.95, linewidth = 0.25, outlier.size = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "red") +
  scale_fill_manual(values = color_map) +
  facet_grid(dataset ~ method_group, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-0.4, 0.55)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

# ============================================================
# HEATMAP TABLE
# ============================================================

heat_tbl <- plot_data %>%
  mutate(row_lab = paste0(gsub("fill_type = ", "fill=", method_group), "\nproxy=", proxy),
    col_lab = case_when(
      dist == "norm"  ~ "Gaussian",
      dist == "t"     ~ paste0("t(", par, ")"),
      dist == "lnorm" ~ paste0("LN(", par, ")"),
      dist == "gamma" ~ paste0("G(", par, ")")
    )
  ) %>%
  group_by(dataset, row_lab, col_lab) %>%
  summarise(
    med_delta = median(delta, na.rm = TRUE),
    win_rate = mean(prtree_win, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# HEATMAPS
# ============================================================

heatmaps <- lapply(
  unique(heat_tbl$dataset),
  function(ds) {
    dat <- heat_tbl %>% filter(dataset == ds)
    ggplot(dat, aes(col_lab, row_lab, fill = med_delta)) +
      geom_tile(color = "white") +
      geom_text(aes(label = paste0(sprintf("%+.1f", 100 * med_delta), "%\n", 
                                   round(100 * win_rate),"%")), size = 3) +
      scale_fill_gradient2(low = "#0072B2",
                           mid = "white",
                           high = "#D55E00",
        midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish,
        name = "Median\nrelative\nchange") +
      labs(title = ds,x = NULL,y = NULL) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", size = 14))
  }
)

heatmaps[[1]]
heatmaps[[2]]
heatmaps[[3]]
heatmaps[[4]]
heatmaps[[5]]

# ============================================================
# BEST PRTREE CONFIG PER DATASET + REP
# ============================================================

best_prtree <- plot_data %>%
  group_by(dataset, rep) %>%
  slice_min(delta, n = 1, with_ties = FALSE) %>%
  ungroup()

compare_tbl <- best_prtree %>%
  transmute(dataset,rep,best_config = paste0("fill=", fill," | ", proxy," | ", dist_label),delta)

# ============================================================
# SINGLE BOXPLOT
# ============================================================

ggplot(compare_tbl, aes(dataset, delta, fill = dataset)) +
  geom_boxplot(alpha = 0.9, linewidth = 0.3, outlier.size = 1) +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.5) + 
  labs(x = NULL,y = expression(frac(RMSE[Best~PRTree] - RMSE[CART],RMSE[CART]))) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20,hjust = 1),
    panel.grid.minor = element_blank()
  )


heatmaps <- lapply(unique(heat_tbl$dataset), function(ds) {
  dat <- heat_tbl %>% filter(dataset == ds)
  ggplot(dat, aes(col_lab, row_lab, fill = win_rate)) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = sprintf("%.0f%%",100*win_rate)),
      size = 3
    ) +
    scale_fill_gradient2(low = "#D55E00",
                           mid = "white",
                           high = "#0072B2",
      midpoint = 0.5,limits = c(0, 1),name = "PRTree\nwin rate"
    ) +
    labs(title = ds,x = NULL,y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 14)
    )
})


heatmaps[[1]]
heatmaps[[2]]
heatmaps[[3]]
heatmaps[[4]]
heatmaps[[5]]


library(ggpubr)

ggarrange(
  heatmaps[[1]],
  heatmaps[[2]],
  heatmaps[[3]],
  heatmaps[[4]],
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

