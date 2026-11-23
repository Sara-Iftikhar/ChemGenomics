library(tidyverse)
library(patchwork)

# ============================================================
# 1. LOAD DATA
# ============================================================

base <- "/Users/daneshm/Documents/ChemicalGenomic/revision"
out <- file.path(base, "ML/ML_results_aggregated")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

agg <- read_csv(file.path(out, "all_models_mean_sd_five_runs.csv"), show_col_types = FALSE)
feat <- read_csv(file.path(base, "feature_count/feature_counts_aggregated.csv"), show_col_types = FALSE)
groups <- read_csv(file.path(dirname(base), "condition_group.csv"), show_col_types = FALSE)

clean <- function(x) gsub("[^A-Za-z0-9]", "", x)

cols <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

# ============================================================
# 2. CONDITION CATEGORIES
# ============================================================

categories <- groups %>%
  mutate(target = clean(File_name)) %>%
  filter(!is.na(Chemical), !is.na(Category)) %>%
  distinct(target, .keep_all = TRUE) %>%
  select(target, Category)

# ============================================================
# 3. BEST MODEL AND OVERFITTING GAPS
# ============================================================

agg <- agg %>% mutate(target = clean(target))

best <- agg %>%
  filter(split == "test", is.finite(rho_mean)) %>%
  group_by(target) %>%
  slice_max(rho_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(target, model, mode, feature_groups)

gaps <- agg %>%
  inner_join(best, by = c("target", "model", "mode", "feature_groups")) %>%
  filter(split %in% c("train_final", "validation_cv", "test")) %>%
  group_by(target, split) %>%
  summarise(rho = mean(rho_mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = split, values_from = rho) %>%
  mutate(
    train_test_gap = train_final - test,
    validation_test_gap = validation_cv - test
  )

# ============================================================
# 4. COMBINE FEATURE COUNTS AND GAPS
# ============================================================

plot_df <- gaps %>%
  left_join(
    feat %>%
      mutate(target = clean(target)) %>%
      select(target, pangenome_mean, pangenome_sd,
             snp_mean, snp_sd, unitig_mean, unitig_sd),
    by = "target"
  ) %>%
  left_join(categories, by = "target") %>%
  filter(!is.na(Category)) %>%
  mutate(Category = factor(Category, levels = names(cols))) %>%
  arrange(Category, target) %>%
  mutate(Condition_number = row_number())

cat("Total conditions:", nrow(plot_df), "\n")
print(plot_df %>% count(Category), n = Inf)

write_csv(plot_df, file.path(out, "FigureS10_All_Data.csv"))
write_csv(
  plot_df %>% select(Condition_number, target, Category),
  file.path(out, "FigureS10_Condition_Index.csv")
)

# ============================================================
# 5. COMMON THEME
# ============================================================

theme_panel <- theme_bw(base_size = 10, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, face = "bold",
                                colour = "black", margin = margin(r = 8)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey50", linewidth = 0.5),
    plot.margin = margin(3, 5, 3, 5),
    legend.position = "none"
  )

xscale <- scale_x_continuous(
  limits = c(0.5, max(plot_df$Condition_number) + 0.5),
  expand = c(0, 0)
)

# ============================================================
# 6. TRAIN–TEST AND VALIDATION–TEST GAPS
# ============================================================

plot_gap <- function(variable, title) {
  ggplot(plot_df, aes(Condition_number, .data[[variable]], fill = Category)) +
    geom_col(width = 0.85) +
    geom_hline(yintercept = 0, linetype = "dotted",
               colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = cols) +
    xscale +
    labs(title = title, y = expression(Delta*rho)) +
    theme_panel
}

# ============================================================
# 7. FEATURE COUNTS WITH SD
# ============================================================

plot_features <- function(mean_col, sd_col, title) {
  
  d <- plot_df %>%
    mutate(
      mean_value = .data[[mean_col]],
      sd_value = .data[[sd_col]],
      ymin = pmax(mean_value - sd_value, 0),
      ymax = mean_value + sd_value,
      ymin_plot = pmax(ymin, 1)
    ) %>%
    filter(is.finite(mean_value), is.finite(sd_value), mean_value > 0)
  
  ggplot(d, aes(Condition_number, mean_value, fill = Category)) +
    geom_col(width = 0.85) +
    geom_errorbar(
      aes(ymin = ymin_plot, ymax = ymax),
      width = 0.2, linewidth = 0.25, colour = "grey30"
    ) +
    scale_fill_manual(values = cols) +
    scale_y_log10(
      limits = c(1, NA),
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      labels = scales::label_comma()
    ) +
    xscale +
    labs(title = title, y = "Feature count") +
    theme_panel
}

# ============================================================
# 8. GENERATE FIVE PANELS
# ============================================================

p1 <- plot_gap("train_test_gap", "Train–test gap")
p2 <- plot_gap("validation_test_gap", "Validation–test gap")
p3 <- plot_features("pangenome_mean", "pangenome_sd", "Pangenome")
p4 <- plot_features("snp_mean", "snp_sd", "SNP")
p5 <- plot_features("unitig_mean", "unitig_sd", "Unitig") +
  scale_x_continuous(
    limits = c(0.5, max(plot_df$Condition_number) + 0.5),
    breaks = seq(5, max(plot_df$Condition_number), 5),
    expand = c(0, 0)
  ) +
  labs(x = "Condition") +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 8, colour = "black",
                               angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_line(linewidth = 0.4)
  )

# ============================================================
# 9. COMBINE AND SAVE
# ============================================================

final_plot <- (p1 / p2 / p3 / p4 / p5) +
  plot_layout(heights = rep(1, 5))

print(final_plot)

ggsave(
  file.path(out, "FigureS10_Overfitting_FeatureCounts_SD.pdf"),
  final_plot, width = 15, height = 7, device = cairo_pdf
)

ggsave(
  file.path(out, "FigureS10_Overfitting_FeatureCounts_SD.png"),
  final_plot, width = 15, height = 7, dpi = 600, bg = "white"
)

cat("Completed: five-panel figure, mean ± SD, condition indices.\n")