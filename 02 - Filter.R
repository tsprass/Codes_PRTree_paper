# ============================================================================
# VISUAL DATASET EXPLORER
# ============================================================================

#-----------------------------------------------------------------------------
# Load and filter the metadata
#-----------------------------------------------------------------------------

# Load the metadata of all UCI datasets (saved in the previous step)
data <- readRDS("data/uci_all_metadata.rds")
head(data)
View(data)

# Regression datasets that contain missing values, sorted by number of instances
with_missing <- data %>%
  dplyr::filter(
    !is.na(tasks),
    grepl("Regression", tasks, ignore.case = TRUE),  # grep direto (base R)
    !is.na(instances),
    has_missing == "Yes"
  ) %>%
  dplyr::arrange(dplyr::desc(instances))

with_missing$id

# Regression datasets with 300-30000 instances and 3-99 features (final selection)
regression_metadata <- data %>%
  dplyr::filter(
    !is.na(tasks),
    grepl("Regression", tasks, ignore.case = TRUE),  # grep direto (base R)
    !is.na(instances),
    instances > 300 & instances < 30000,
    features > 2 & features < 100,
  ) %>%
  dplyr::arrange(dplyr::desc(instances))

View(regression_metadata)
View(regression_metadata %>% dplyr::arrange(id))

#-----------------------------------------------------------------------------
# 1. SUMMARY STATISTICS
#-----------------------------------------------------------------------------
summary_stats <- data.frame(
  Metric = c("Total datasets", "Min instances", "Max instances", "Mean instances",
             "Min features", "Max features", "Mean features",
             "Has missing", "No missing", "Has target name"),
  Value = c(
    nrow(regression_metadata),
    min(regression_metadata$instances, na.rm = TRUE),
    max(regression_metadata$instances, na.rm = TRUE),
    round(mean(regression_metadata$instances, na.rm = TRUE), 0),
    min(regression_metadata$features, na.rm = TRUE),
    max(regression_metadata$features, na.rm = TRUE),
    round(mean(regression_metadata$features, na.rm = TRUE), 1),
    sum(regression_metadata$has_missing == "Yes", na.rm = TRUE),
    sum(regression_metadata$has_missing == "No", na.rm = TRUE),
    sum(!is.na(regression_metadata$target_names))
  )
)

print(summary_stats, row.names = FALSE)

#-----------------------------------------------------------------------------
# 2. DISTRIBUTION PLOTS
#-----------------------------------------------------------------------------

# Distribution of instances (log scale)
p1 <- regression_metadata %>%
  ggplot2::ggplot(ggplot2::aes(x = instances)) +
  ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(title = "Dataset Size Distribution", x = "Instances (log10)", y = "Count") +
  ggplot2::theme_minimal()

# Distribution of features
p2 <- regression_metadata %>%
  ggplot2::ggplot(ggplot2::aes(x = features)) +
  ggplot2::geom_histogram(bins = 30, fill = "coral", color = "black") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(title = "Feature Count Distribution", x = "Features (log10)", y = "Count") +
  ggplot2::theme_minimal()

# Missing values
p3 <- regression_metadata %>%
  dplyr::count(has_missing) %>%
  ggplot2::ggplot(ggplot2::aes(x = has_missing, y = n, fill = has_missing)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(title = "Missing Values", x = "Has Missing?", y = "Count") +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values = c("Yes" = "coral", "No" = "steelblue")) +
  ggplot2::theme(legend.position = "none")

# Target availability
p4 <- regression_metadata %>%
  dplyr::mutate(has_target = ifelse(is.na(target_names), "No", "Yes")) %>%
  dplyr::count(has_target) %>%
  ggplot2::ggplot(ggplot2::aes(x = has_target, y = n, fill = has_target)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(title = "Target Variable", x = "Has Target Name?", y = "Count") +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values = c("Yes" = "steelblue", "No" = "coral")) +
  ggplot2::theme(legend.position = "none")

# Arrange plots
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
