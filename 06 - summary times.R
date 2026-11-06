library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(ggplot2)

# ============================================================
# LOG FILES
# ============================================================

files <- list.files("results",  pattern = "\\.log$", full.names = TRUE)

# ============================================================
# PARSE LOG
# ============================================================

parse_log <- function(file){
  
  txt <- readLines(file)
  dataset <- str_remove(basename(file), "\\.log$")
  start_line <- txt[str_detect(txt, "START")][1]
  end_line <- txt[str_detect(txt, "END")][1]
  start_time <- ymd_hms(
    str_extract(start_line,
      "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"))
  end_time <- ymd_hms(
    str_extract(end_line,
      "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"))
  cfg_lines <- txt[str_detect(txt, "Running PRTree")]
  cfg_times <- ymd_hms(
    str_extract(cfg_lines,
      "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"))
  cart_total <- as.numeric(difftime(cfg_times[1], start_time, units = "secs"))
  elapsed_sec <- c(
    diff(as.numeric(cfg_times)), 
    as.numeric(difftime(end_time,tail(cfg_times, 1),units = "secs")))
  fill <- str_extract(cfg_lines,"(?<=fill = )\\d")
  proxy <- str_extract(cfg_lines,"(?<=proxy = )[a-z]+")
  dist <- str_extract(cfg_lines,"(?<=dist = )[a-z]+")
  cfg_tbl <- tibble(dataset, fill, proxy, dist, elapsed_sec)
  
  summary_tbl <- tibble(
    dataset = dataset,
    cart_total_sec = cart_total,
    prtree_total_sec = sum(elapsed_sec),
    prtree_min_sec = min(elapsed_sec),
    prtree_max_sec = max(elapsed_sec)
  )
  
  list(cfg = cfg_tbl, summary = summary_tbl)
}

# ============================================================
# READ ALL LOGS
# ============================================================

logs <- map(files, parse_log)
cfg_tbl <- map_dfr(logs, "cfg")
time_summary <- map_dfr(logs, "summary") %>%
  mutate(
    cart_total_sec = round(cart_total_sec, 2),
    prtree_total_min = round(prtree_total_sec / 60, 2),
    prtree_min_sec = round(prtree_min_sec, 2),
    prtree_max_sec = round(prtree_max_sec, 2)
  )

# ============================================================
# ORIGINAL TABLE
# ============================================================

tbla <- time_summary %>%
  transmute(
    Dataset = dataset,
    `CART total (s)` = cart_total_sec,
    `PRTree total (min)` = prtree_total_min,
    `PRTree min (s)` = prtree_min_sec,
    `PRTree max (s)` = prtree_max_sec
  )

knitr::kable(
  tbla,
  format = "latex",
  booktabs = TRUE,
  digits = 2,
  caption = "Computational cost of CART and PRTree."
)

# ============================================================
# HEATMAP TABLE
# ============================================================
time_heat <- cfg_tbl %>%
  group_by(dataset, fill, proxy, dist) %>%
  mutate(idx = row_number()) %>%
  ungroup() %>%
  mutate(
    row_lab = paste0("f", fill, "-",
      case_when(
        proxy == "both" ~ "B",
        proxy == "mean" ~ "M",
        proxy == "var"  ~ "V"
      )),
    col_lab = case_when(
      dist == "norm" ~ "N",
      dist == "t" & idx == 1 ~ "t(3)",
      dist == "t" & idx == 2 ~ "t(5)",
      dist == "lnorm" & idx == 1 ~ "LN(.5)",
      dist == "lnorm" & idx == 2 ~ "LN(1)",
      dist == "lnorm" & idx == 3 ~ "LN(2)",
      dist == "gamma" & idx == 1 ~ "G(.5)",
      dist == "gamma" & idx == 2 ~ "G(2)",
      dist == "gamma" & idx == 3 ~ "G(5)"
    )) %>%
  group_by(dataset, row_lab, col_lab) %>%
  summarise(mean_time = mean(elapsed_sec) / 60, .groups = "drop")

dist_levels <- c(
  "N",
  "t(3)",
  "t(5)",
  "LN(.5)",
  "LN(1)",
  "LN(2)",
  "G(.5)",
  "G(2)",
  "G(5)"
)

time_heat$col_lab <- factor(
  time_heat$col_lab,
  levels = dist_levels
)

global_min <- min(time_heat$mean_time)
global_max <- max(time_heat$mean_time)


time_heatmaps <- lapply(
  unique(time_heat$dataset),
  function(ds){
    
    dat <- time_heat %>%
      filter(dataset == ds)
    
    ggplot(
      dat,
      aes(
        col_lab,
        row_lab,
        fill = mean_time
      )
    ) +
      geom_tile(color = "white") +
      geom_text(
        aes(
          label = sprintf("%.1f", mean_time)
        ),
        size = 4.5
      ) +
      scale_fill_gradient2(
        low = "#0072B2",
        mid = "white",
        high = "#D55E00",
        limits = c(0, 15),
        midpoint = 7,
        oob = scales::squish,
        name = "Time\n(min)"
      ) +
      labs(
        title = ds,
        x = NULL,
        y = NULL
      ) +
      theme_bw(base_size = 15) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5
        ),
        plot.title = element_text(face = "bold")
      )
  }
)

p_leg <- time_heatmaps[[1]] +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(14, "cm"),
      barheight = unit(1.5, "cm")
    )
  )

leg <- get_legend(p_leg)

time_heatmaps_noleg <- lapply(
  time_heatmaps, 
  function(p) p + theme(legend.position = "none")
)

h_time_final <- ggarrange(
  time_heatmaps_noleg[[1]],
  time_heatmaps_noleg[[2]],
  time_heatmaps_noleg[[3]],
  time_heatmaps_noleg[[4]],
  time_heatmaps_noleg[[5]],
  ggpubr::as_ggplot(leg),
  ncol = 2,
  nrow = 3
)


ggsave(
  "heatmapsTime.pdf",
  plot = h_time_final,
  device = cairo_pdf,
  width = 15,
  height = 15
)
