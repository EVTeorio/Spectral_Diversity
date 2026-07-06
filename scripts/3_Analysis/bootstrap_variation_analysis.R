user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "scales", "knitr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(knitr)

setwd("C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity")

analysis_dir <- "reports/analysis"
figure_dir <- "reports/figures/bootstrap_variation"
table_dir <- "reports/tables/bootstrap_variation"

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

scale_levels <- c("10m", "20m", "50m")
scale_labels <- c("10m" = "10 m", "20m" = "20 m", "50m" = "50 m")
N_BOOT_EXPECTED <- 70
BOOT_SAMPLE_PIXELS <- 5000
PAIR_SAMPLE <- 10000
EXACT_PAIR_THRESHOLD <- 2000000

read_scale_data <- function(scale) {
  summary_path <- file.path("Quad_Values", paste0(scale, "_SA_entropy_smooth_masked_5nm_summary.csv"))
  wide_path <- file.path("Quad_Values", paste0(scale, "_SA_entropy_smooth_masked_5nm_boot_wide.csv"))

  if (!file.exists(summary_path)) {
    stop("Missing summary file: ", summary_path)
  }
  if (!file.exists(wide_path)) {
    stop("Missing bootstrap wide file: ", wide_path)
  }

  summary_df <- read_csv(summary_path, show_col_types = FALSE) |>
    mutate(
      scale = scale,
      quad_id = as.character(quad_id)
    )

  wide_df <- read_csv(wide_path, show_col_types = FALSE) |>
    mutate(
      scale = scale,
      quad_id = as.character(quad_id)
    )

  boot_cols <- grep("^boot_", names(wide_df), value = TRUE)
  if (length(boot_cols) != N_BOOT_EXPECTED) {
    stop("Expected ", N_BOOT_EXPECTED, " bootstrap columns for ", scale, " but found ", length(boot_cols))
  }

  boot_long <- wide_df |>
    pivot_longer(
      cols = all_of(boot_cols),
      names_to = "boot_iter",
      values_to = "boot_entropy"
    ) |>
    mutate(
      boot_iter = as.integer(sub("^boot_", "", boot_iter))
    )

  list(summary = summary_df, boot_long = boot_long)
}

scale_data <- lapply(scale_levels, read_scale_data)
summary_all <- bind_rows(lapply(scale_data, `[[`, "summary")) |>
  mutate(
    scale = factor(scale, levels = scale_levels, labels = scale_labels),
    boot_se = boot_sd / sqrt(N_BOOT_EXPECTED),
    valid_entropy = is.finite(spectral_entropy)
  )

boot_long_all <- bind_rows(lapply(scale_data, `[[`, "boot_long")) |>
  mutate(
    scale = factor(scale, levels = scale_levels, labels = scale_labels)
  )

bootstrap_quads <- summary_all |>
  filter(method == "bootstrap_mean", is.finite(boot_sd), is.finite(boot_cv))

calc_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || stats::sd(x) == 0) {
    return(NA_real_)
  }

  mean((x - mean(x))^3) / stats::sd(x)^3
}

replicate_diagnostics <- boot_long_all |>
  filter(is.finite(boot_entropy)) |>
  group_by(scale, quad_id) |>
  summarise(
    n_boot = n(),
    boot_mean_from_reps = mean(boot_entropy),
    boot_sd_from_reps = sd(boot_entropy),
    boot_cv_from_reps = boot_sd_from_reps / boot_mean_from_reps,
    boot_q25 = quantile(boot_entropy, 0.25),
    boot_q75 = quantile(boot_entropy, 0.75),
    boot_iqr = IQR(boot_entropy),
    boot_min = min(boot_entropy),
    boot_max = max(boot_entropy),
    boot_range = boot_max - boot_min,
    low_fence = boot_q25 - 1.5 * boot_iqr,
    high_fence = boot_q75 + 1.5 * boot_iqr,
    outlier_reps = sum(boot_entropy < low_fence | boot_entropy > high_fence),
    .groups = "drop"
  )

analysis_df <- bootstrap_quads |>
  left_join(
    replicate_diagnostics |>
      select(scale, quad_id, n_boot, boot_range, outlier_reps),
    by = c("scale", "quad_id")
  ) |>
  mutate(
    n_boot = if_else(is.na(n_boot), N_BOOT_EXPECTED, n_boot),
    boot_se = boot_sd / sqrt(n_boot),
    boot_ci_t = qt(0.975, df = pmax(n_boot - 1, 1)),
    boot_ci_half_width = boot_ci_t * boot_se,
    boot_ci_low = spectral_entropy - boot_ci_half_width,
    boot_ci_high = spectral_entropy + boot_ci_half_width,
    boot_ci_width = boot_ci_high - boot_ci_low,
    boot_ci_half_width_pct = boot_ci_half_width / spectral_entropy
  )

scale_summary <- summary_all |>
  group_by(scale) |>
  summarise(
    summary_rows = n(),
    valid_entropy_n = sum(valid_entropy),
    missing_entropy_n = sum(!valid_entropy),
    bootstrap_quads_n = sum(method == "bootstrap_mean"),
    exact_quads_n = sum(method == "exact_all_pixels"),
    insufficient_quads_n = sum(method == "insufficient_pixels"),
    between_sd = sd(spectral_entropy, na.rm = TRUE),
    between_var = var(spectral_entropy, na.rm = TRUE),
    between_iqr = IQR(spectral_entropy, na.rm = TRUE),
    mean_entropy = mean(spectral_entropy, na.rm = TRUE),
    median_entropy = median(spectral_entropy, na.rm = TRUE),
    entropy_skewness = calc_skewness(spectral_entropy),
    .groups = "drop"
  ) |>
  left_join(
    analysis_df |>
      group_by(scale) |>
      summarise(
        within_sd_mean = mean(boot_sd, na.rm = TRUE),
        within_sd_median = median(boot_sd, na.rm = TRUE),
        within_sd_p90 = quantile(boot_sd, 0.90, na.rm = TRUE),
        within_sd_max = max(boot_sd, na.rm = TRUE),
        within_cv_mean = mean(boot_cv, na.rm = TRUE),
        within_cv_median = median(boot_cv, na.rm = TRUE),
        within_cv_p90 = quantile(boot_cv, 0.90, na.rm = TRUE),
        within_cv_max = max(boot_cv, na.rm = TRUE),
        boot_se_mean = mean(boot_se, na.rm = TRUE),
        boot_se_median = median(boot_se, na.rm = TRUE),
        boot_se_p90 = quantile(boot_se, 0.90, na.rm = TRUE),
        boot_ci_half_width_mean = mean(boot_ci_half_width, na.rm = TRUE),
        boot_ci_half_width_median = median(boot_ci_half_width, na.rm = TRUE),
        boot_ci_half_width_p90 = quantile(boot_ci_half_width, 0.90, na.rm = TRUE),
        boot_ci_half_width_max = max(boot_ci_half_width, na.rm = TRUE),
        boot_ci_half_width_pct_median = median(boot_ci_half_width_pct, na.rm = TRUE),
        boot_ci_half_width_pct_p90 = quantile(boot_ci_half_width_pct, 0.90, na.rm = TRUE),
        mean_within_var = mean(boot_sd^2, na.rm = TRUE),
        median_outlier_reps = median(outlier_reps, na.rm = TRUE),
        quads_cv_gt_005 = sum(boot_cv > 0.05, na.rm = TRUE),
        quads_cv_gt_010 = sum(boot_cv > 0.10, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "scale"
  ) |>
  mutate(
    sd_ratio_mean_within_to_between = within_sd_mean / between_sd,
    sd_ratio_median_within_to_between = within_sd_median / between_sd,
    se_ratio_mean_to_between = boot_se_mean / between_sd,
    within_to_between_var_ratio = mean_within_var / between_var,
    icc_like_between_fraction = between_var / (between_var + mean_within_var),
    pct_cv_gt_005 = quads_cv_gt_005 / bootstrap_quads_n,
    pct_cv_gt_010 = quads_cv_gt_010 / bootstrap_quads_n,
    stability_class = case_when(
      se_ratio_mean_to_between <= 0.05 & within_to_between_var_ratio <= 0.10 ~ "Strong",
      se_ratio_mean_to_between <= 0.10 & within_to_between_var_ratio <= 0.25 ~ "Acceptable",
      TRUE ~ "Review"
    )
  )

top_unstable <- analysis_df |>
  arrange(scale, desc(boot_cv), desc(boot_sd)) |>
  group_by(scale) |>
  slice_head(n = 15) |>
  ungroup()

top_ci <- analysis_df |>
  arrange(scale, desc(boot_ci_half_width), desc(boot_sd)) |>
  group_by(scale) |>
  slice_head(n = 15) |>
  ungroup()

top_unstable_long <- boot_long_all |>
  semi_join(top_unstable, by = c("scale", "quad_id")) |>
  filter(is.finite(boot_entropy)) |>
  left_join(
    top_unstable |>
      select(scale, quad_id, boot_cv, boot_sd),
    by = c("scale", "quad_id")
  ) |>
  mutate(
    quad_id_ordered = reorder(quad_id, boot_cv)
  )

write_csv(scale_summary, file.path(table_dir, "bootstrap_variation_scale_summary.csv"))
write_csv(analysis_df, file.path(table_dir, "bootstrap_variation_quadrat_diagnostics.csv"))
write_csv(top_unstable, file.path(table_dir, "bootstrap_variation_top_unstable_quadrats.csv"))
write_csv(top_ci, file.path(table_dir, "bootstrap_variation_widest_ci_quadrats.csv"))

theme_report <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

variation_long <- scale_summary |>
  select(scale, between_sd, within_sd_mean, within_sd_median, boot_se_mean, boot_ci_half_width_mean) |>
  pivot_longer(
    cols = c(between_sd, within_sd_mean, within_sd_median, boot_se_mean, boot_ci_half_width_mean),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      between_sd = "Between-quadrat SD",
      within_sd_mean = "Mean within-quadrat bootstrap SD",
      within_sd_median = "Median within-quadrat bootstrap SD",
      boot_se_mean = "Mean SE of bootstrap mean",
      boot_ci_half_width_mean = "Mean 95% CI half-width"
    )
  )

p1 <- ggplot(variation_long, aes(x = scale, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "Bootstrap uncertainty compared with between-quadrat variation",
    subtitle = "Most SA entropy values are means of 70 bootstrap iterations; exact all-pair entropy is limited to small quadrats.",
    x = NULL,
    y = "SA entropy mean units"
  ) +
  theme_report
ggsave(file.path(figure_dir, "within_vs_between_sd.png"), p1, width = 11, height = 6.2, dpi = 300, bg = "white")

p2 <- ggplot(analysis_df, aes(x = boot_cv)) +
  geom_histogram(bins = 45, fill = "#3C8DAD", color = "white") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#9B2226") +
  geom_vline(xintercept = 0.10, linetype = "dotted", color = "#5F0F40") +
  facet_wrap(~ scale, scales = "free_y") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Distribution of within-quadrat bootstrap CV",
    x = "Bootstrap coefficient of variation",
    y = "Quadrat count"
  ) +
  theme_report
ggsave(file.path(figure_dir, "bootstrap_cv_distribution.png"), p2, width = 9, height = 5.5, dpi = 300, bg = "white")

p3 <- ggplot(analysis_df, aes(x = spectral_entropy, y = boot_sd, color = scale)) +
  geom_point(alpha = 0.55, size = 1.7) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  scale_color_brewer(palette = "Dark2", name = "Scale") +
  labs(
    title = "Bootstrap SD by SA entropy mean",
    subtitle = "Bootstrap rows summarize 70 pixel-subsampling iterations, not exhaustive all-retained-pixel pair calculations.",
    x = "SA entropy mean",
    y = "Within-quadrat bootstrap SD"
  ) +
  theme_report
ggsave(file.path(figure_dir, "bootstrap_sd_vs_entropy.png"), p3, width = 9, height = 5.5, dpi = 300, bg = "white")

p4 <- ggplot(analysis_df, aes(x = n_pixels, y = boot_cv, color = scale)) +
  geom_point(alpha = 0.55, size = 1.7) +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Dark2", name = "Scale") +
  labs(
    title = "Bootstrap CV versus retained sunlit pixel count",
    x = "Valid sunlit pixels after masking (log scale)",
    y = "Bootstrap coefficient of variation"
  ) +
  theme_report
ggsave(file.path(figure_dir, "pixel_count_vs_bootstrap_cv.png"), p4, width = 9, height = 5.5, dpi = 300, bg = "white")

ratio_long <- scale_summary |>
  select(scale, within_to_between_var_ratio, icc_like_between_fraction) |>
  pivot_longer(
    cols = c(within_to_between_var_ratio, icc_like_between_fraction),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      within_to_between_var_ratio = "Mean within variance / between variance",
      icc_like_between_fraction = "Between fraction of total variance"
    )
  )

p5 <- ggplot(ratio_long, aes(x = scale, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Paired", name = NULL) +
  labs(
    title = "Variance partition diagnostic",
    x = NULL,
    y = "Ratio"
  ) +
  theme_report
ggsave(file.path(figure_dir, "variance_partition_diagnostic.png"), p5, width = 9, height = 5.5, dpi = 300, bg = "white")

p6 <- ggplot(top_unstable_long, aes(x = quad_id_ordered, y = boot_entropy)) +
  geom_boxplot(outlier.alpha = 0.35, fill = "#DDB892", color = "#4A4E69") +
  coord_flip() +
  facet_wrap(~ scale, scales = "free_y") +
  labs(
    title = "Bootstrap replicate spread for the most variable quadrats",
    x = "Quadrat ID, ordered by bootstrap CV",
    y = "Bootstrap SA entropy value"
  ) +
  theme_report
ggsave(file.path(figure_dir, "top_unstable_quadrat_bootstrap_distributions.png"), p6, width = 10, height = 8, dpi = 300, bg = "white")

uncertainty_long <- analysis_df |>
  transmute(
    scale,
    `Bootstrap CV (%)` = boot_cv * 100,
    `Bootstrap SD` = boot_sd,
    `SE of mean` = boot_se,
    `95% CI half-width` = boot_ci_half_width
  ) |>
  pivot_longer(-scale, names_to = "metric", values_to = "value")

p8 <- ggplot(uncertainty_long, aes(x = scale, y = value, fill = scale)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.62) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_fill_brewer(palette = "Dark2", name = "Scale") +
  labs(
    title = "Bootstrap variation and mean uncertainty by scale",
    subtitle = "SE and 95 percent CI summarize confidence in the 70-iteration mean; SD and CV show raw replicate spread.",
    x = "Quadrat scale",
    y = NULL
  ) +
  theme_report +
  theme(legend.position = "none")
ggsave(file.path(figure_dir, "bootstrap_uncertainty_by_scale.png"), p8, width = 10, height = 7.2, dpi = 300, bg = "white")

p9 <- ggplot(analysis_df, aes(x = boot_ci_half_width)) +
  geom_histogram(bins = 45, fill = "#4F9D69", color = "white") +
  facet_wrap(~ scale, scales = "free_y") +
  labs(
    title = "Distribution of 95 percent CI half-widths for the bootstrap mean",
    subtitle = "Narrower intervals mean higher confidence in the reported 70-iteration mean.",
    x = "95 percent CI half-width",
    y = "Quadrat count"
  ) +
  theme_report
ggsave(file.path(figure_dir, "bootstrap_ci_width_distribution.png"), p9, width = 9, height = 5.5, dpi = 300, bg = "white")

top_ci_plot <- top_ci |>
  mutate(quad_id_ordered = reorder(quad_id, boot_ci_half_width))

p10 <- ggplot(top_ci_plot, aes(x = quad_id_ordered, y = spectral_entropy)) +
  geom_errorbar(aes(ymin = boot_ci_low, ymax = boot_ci_high), width = 0.22, color = "#7A3E48") +
  geom_point(size = 1.8, color = "#1B4D89") +
  coord_flip() +
  facet_wrap(~ scale, scales = "free_y") +
  labs(
    title = "Widest 95 percent CIs around the bootstrap mean",
    subtitle = "Intervals are boot_mean +/- t * boot_sd / sqrt(70), shown for the highest-uncertainty quadrats in each scale.",
    x = "Quadrat ID, ordered by CI half-width",
    y = "SA entropy mean with 95 percent CI"
  ) +
  theme_report
ggsave(file.path(figure_dir, "bootstrap_mean_ci_widest_quadrats.png"), p10, width = 10, height = 8, dpi = 300, bg = "white")

entropy_reference <- summary_all |>
  filter(valid_entropy) |>
  group_by(scale) |>
  summarise(
    mean_entropy = mean(spectral_entropy),
    median_entropy = median(spectral_entropy),
    .groups = "drop"
  )

p7 <- ggplot(summary_all |> filter(valid_entropy), aes(x = spectral_entropy)) +
  geom_histogram(bins = 35, fill = "#6A994E", color = "white", alpha = 0.85) +
  geom_vline(
    data = entropy_reference,
    aes(xintercept = mean_entropy),
    color = "#9B2226",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    data = entropy_reference,
    aes(xintercept = median_entropy),
    color = "#005F73",
    linetype = "solid",
    linewidth = 0.8
  ) +
  facet_wrap(~ scale, scales = "free_y") +
  labs(
    title = "Distribution of SA entropy means by scale",
    subtitle = "Dashed red = mean; solid blue = median. Most values are means of 70 bootstrap iterations.",
    x = "SA entropy mean",
    y = "Quadrat count"
  ) +
  theme_report
ggsave(file.path(figure_dir, "spectral_entropy_histograms_by_scale.png"), p7, width = 9, height = 5.8, dpi = 300, bg = "white")

format_number <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

scale_summary_report <- scale_summary |>
  transmute(
    Scale = as.character(scale),
    `Valid SA values` = valid_entropy_n,
    `Bootstrap quads` = bootstrap_quads_n,
    `Exact quads` = exact_quads_n,
    `Scale mean entropy` = format_number(mean_entropy),
    `Median entropy` = format_number(median_entropy),
    `Skewness` = format_number(entropy_skewness),
    `Between SD` = format_number(between_sd),
    `Median within SD` = format_number(within_sd_median),
    `Mean within SD` = format_number(within_sd_mean),
    `Mean boot mean SE` = format_number(boot_se_mean),
    `Median 95% CI half-width` = format_number(boot_ci_half_width_median),
    `P90 95% CI half-width` = format_number(boot_ci_half_width_p90),
    `Within/Between variance` = percent(within_to_between_var_ratio, accuracy = 0.1),
    `Between fraction` = percent(icc_like_between_fraction, accuracy = 0.1),
    `Median CV` = percent(within_cv_median, accuracy = 0.1),
    `P90 CV` = percent(within_cv_p90, accuracy = 0.1),
    `Class` = stability_class
  )

top_unstable_report <- top_unstable |>
  group_by(scale) |>
  slice_head(n = 8) |>
  ungroup() |>
  transmute(
    Scale = as.character(scale),
    Quadrat = quad_id,
    `Pixels` = n_pixels,
    `Entropy mean` = format_number(spectral_entropy),
    `Boot SD` = format_number(boot_sd),
    `95% CI low` = format_number(boot_ci_low),
    `95% CI high` = format_number(boot_ci_high),
    `Boot CV` = percent(boot_cv, accuracy = 0.1),
    `Outlier reps` = outlier_reps
  )

top_ci_report <- top_ci |>
  group_by(scale) |>
  slice_head(n = 8) |>
  ungroup() |>
  transmute(
    Scale = as.character(scale),
    Quadrat = quad_id,
    `Pixels` = n_pixels,
    `Entropy mean` = format_number(spectral_entropy),
    `95% CI low` = format_number(boot_ci_low),
    `95% CI high` = format_number(boot_ci_high),
    `95% CI half-width` = format_number(boot_ci_half_width),
    `Half-width % of mean` = percent(boot_ci_half_width_pct, accuracy = 0.1)
  )

overall_max_cv <- max(analysis_df$boot_cv, na.rm = TRUE)
overall_p95_cv <- quantile(analysis_df$boot_cv, 0.95, na.rm = TRUE)
worst_quad <- analysis_df |>
  arrange(desc(boot_cv), desc(boot_sd)) |>
  slice(1)

validity_statement <- if (all(scale_summary$stability_class %in% c("Strong", "Acceptable"))) {
  paste(
    "The bootstrap results support using the per-quadrat spectral entropy means.",
    "For all scales, the uncertainty of the reported bootstrap mean is small relative to between-quadrat variation."
  )
} else {
  paste(
    "The bootstrap means are usable with quality-control flags, but raw within-quadrat bootstrap variation is not uniformly small enough to ignore.",
    "A subset of quadrats, especially at 10 m and 20 m, should be flagged or included in sensitivity checks before downstream modeling."
  )
}

report_path <- file.path(analysis_dir, "20260618_bootstrap_variation_analysis.md")

report_lines <- c(
  "# Bootstrap Variation Analysis for Spectral Heterogeneity",
  "",
  "Date: 2026-06-18",
  "",
  "Updated: 2026-07-03",
  "",
  "## Purpose",
  "",
  "This report evaluates whether variation among bootstrap iterations within quadrats is small enough relative to variation among quadrats to use the SA entropy mean values in downstream biodiversity analyses.",
  "",
  "The key distinction is that `boot_sd` measures how much spectral entropy changes across repeated pixel subsamples within the same quadrat, while `boot_sd / sqrt(70)` is the standard error of the reported `boot_mean`.",
  "",
  "## Inputs",
  "",
  "- `Quad_Values/10m_SA_entropy_smooth_masked_5nm_summary.csv` and matching bootstrap files",
  "- `Quad_Values/20m_SA_entropy_smooth_masked_5nm_summary.csv` and matching bootstrap files",
  "- `Quad_Values/50m_SA_entropy_smooth_masked_5nm_summary.csv` and matching bootstrap files",
  "",
  "All inputs were generated from current smoothed 5 nm spectra using the shadow mask threshold of `0.0305476` at `563 nm`.",
  "",
  paste0(
    "For `method = \"bootstrap_mean\"` rows, the reported SA value is the mean of ",
    N_BOOT_EXPECTED,
    " bootstrap iterations. Each iteration samples up to ",
    comma(BOOT_SAMPLE_PIXELS),
    " retained pixels without replacement; when that sample still has too many possible pairs, the entropy is calculated from ",
    comma(PAIR_SAMPLE),
    " sampled spectral-angle pairs. Rows with `method = \"exact_all_pixels\"` are the small subset below ",
    comma(EXACT_PAIR_THRESHOLD),
    " possible retained-pixel pairs, roughly no more than 2,000 retained pixels, and use all retained pixel pairs."
  ),
  "",
  "## Main Conclusion",
  "",
  validity_statement,
  "",
  "Decision: use the reported `spectral_entropy` / `boot_mean` values as the main heterogeneity means, but do not treat bootstrap rows as exhaustive all-pair entropy. Carry `boot_sd`, `boot_cv`, `boot_se`, and the 95 percent CI fields forward as quality-control fields, and run sensitivity checks that flag or exclude high-CV or wide-CI quadrats.",
  "",
  paste0(
    "Across all bootstrap-mean quadrats, the 95th percentile bootstrap CV was ",
    percent(overall_p95_cv, accuracy = 0.1),
    ", and the maximum bootstrap CV was ",
    percent(overall_max_cv, accuracy = 0.1),
    " in quadrat `",
    worst_quad$quad_id,
    "` at ",
    as.character(worst_quad$scale),
    "."
  ),
  "",
  "The most important validation metric is the ratio of mean within-quadrat bootstrap variance to between-quadrat variance. Lower values mean the bootstrap noise is small relative to the spatial signal among quadrats.",
  "",
  paste0(
    "The median 95 percent CI half-width around the bootstrap mean was ",
    paste(
      paste0(
        as.character(scale_summary$scale),
        ": ",
        format_number(scale_summary$boot_ci_half_width_median),
        " (",
        percent(scale_summary$boot_ci_half_width_pct_median, accuracy = 0.1),
        " of the mean)"
      ),
      collapse = "; "
    ),
    "."
  ),
  "",
  paste0(
    "Using a 5% bootstrap-CV flag, ",
    scale_summary$quads_cv_gt_005[scale_summary$scale == "10 m"],
    " of ",
    scale_summary$bootstrap_quads_n[scale_summary$scale == "10 m"],
    " bootstrap-mean 10 m quadrats, ",
    scale_summary$quads_cv_gt_005[scale_summary$scale == "20 m"],
    " of ",
    scale_summary$bootstrap_quads_n[scale_summary$scale == "20 m"],
    " 20 m quadrats, and ",
    scale_summary$quads_cv_gt_005[scale_summary$scale == "50 m"],
    " of ",
    scale_summary$bootstrap_quads_n[scale_summary$scale == "50 m"],
    " 50 m quadrats would be flagged."
  ),
  "",
  "## Scale Summary",
  "",
  paste(capture.output(knitr::kable(scale_summary_report, format = "markdown", align = "r")), collapse = "\n"),
  "",
  "## Figures",
  "",
  "### Within Versus Between SD",
  "",
  "![Within versus between SD](../figures/bootstrap_variation/within_vs_between_sd.png)",
  "",
  "### Bootstrap Variation And Mean Uncertainty By Scale",
  "",
  "This figure separates raw replicate spread (`boot_sd`, `boot_cv`) from uncertainty in the 70-iteration mean (`boot_se`, 95 percent CI half-width).",
  "",
  "![Bootstrap variation and mean uncertainty by scale](../figures/bootstrap_variation/bootstrap_uncertainty_by_scale.png)",
  "",
  "### CI Width Distribution",
  "",
  "![CI width distribution](../figures/bootstrap_variation/bootstrap_ci_width_distribution.png)",
  "",
  "### Widest CI Quadrats",
  "",
  "![Widest CIs around bootstrap mean](../figures/bootstrap_variation/bootstrap_mean_ci_widest_quadrats.png)",
  "",
  "### Bootstrap CV Distribution",
  "",
  "Dashed and dotted reference lines mark 5% and 10% bootstrap CV.",
  "",
  "![Bootstrap CV distribution](../figures/bootstrap_variation/bootstrap_cv_distribution.png)",
  "",
  "### Bootstrap SD Versus SA Entropy Mean",
  "",
  "![Bootstrap SD versus SA entropy mean](../figures/bootstrap_variation/bootstrap_sd_vs_entropy.png)",
  "",
  "### SA Entropy Mean Histograms",
  "",
  "These histograms show the distribution and skewness of the final SA entropy means at each scale. The dashed red line is the scale mean and the solid blue line is the scale median. Most quadrat values are means of 70 bootstrap iterations.",
  "",
  "![SA entropy mean histograms by scale](../figures/bootstrap_variation/spectral_entropy_histograms_by_scale.png)",
  "",
  "### Pixel Count Versus Bootstrap CV",
  "",
  "![Pixel count versus bootstrap CV](../figures/bootstrap_variation/pixel_count_vs_bootstrap_cv.png)",
  "",
  "### Variance Partition Diagnostic",
  "",
  "![Variance partition diagnostic](../figures/bootstrap_variation/variance_partition_diagnostic.png)",
  "",
  "### Most Variable Quadrat Replicate Distributions",
  "",
  "![Top unstable quadrat bootstrap distributions](../figures/bootstrap_variation/top_unstable_quadrat_bootstrap_distributions.png)",
  "",
  "## Most Variable Quadrats",
  "",
  paste(capture.output(knitr::kable(top_unstable_report, format = "markdown", align = "r")), collapse = "\n"),
  "",
  "## Widest Bootstrap-Mean Confidence Intervals",
  "",
  paste(capture.output(knitr::kable(top_ci_report, format = "markdown", align = "r")), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "- Median bootstrap CV is the best quick read on typical within-quadrat stability.",
  "- The mean bootstrap SD can be pulled upward by a small number of quadrats with occasional low or high bootstrap replicates.",
  "- The SE and 95 percent CI around the bootstrap mean are much smaller than the raw bootstrap SD because each reported bootstrap value averages 70 iterations.",
  "- The raw within-quadrat bootstrap variance is large enough in some quadrats that it should be represented in downstream model diagnostics.",
  "- Exact all-pixel values were used where computationally reasonable; most quadrats use `bootstrap_mean` because all-pixel pair counts were too large.",
  "- For bootstrap-mean quadrats, `n_pixels` and `n_pairs` describe the retained pixel pool and possible pair count, not the number of pixel pairs actually evaluated in the entropy calculation.",
  "- Some bootstrap variation may reflect the sampled-pair approximation used for computationally large bootstrap replicates; high-CV quadrats are the best candidates for reruns with larger pair samples if more precision is needed.",
  "",
  "## Recommendation",
  "",
  "Use `spectral_entropy` from the per-scale summary CSVs as the primary spectral heterogeneity mean, but carry `boot_sd`, `boot_cv`, `boot_se`, CI fields, and `method` into downstream analysis tables as quality-control fields. Only `method = \"exact_all_pixels\"` rows represent exhaustive all-retained-pixel pairwise entropy. In model sensitivity checks, consider excluding or flagging quadrats with bootstrap CV above 5%, wide 95 percent CIs, or many outlying bootstrap replicates.",
  "",
  "## Output Tables",
  "",
  "- `reports/tables/bootstrap_variation/bootstrap_variation_scale_summary.csv`",
  "- `reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv`",
  "- `reports/tables/bootstrap_variation/bootstrap_variation_top_unstable_quadrats.csv`",
  "- `reports/tables/bootstrap_variation/bootstrap_variation_widest_ci_quadrats.csv`"
)

writeLines(report_lines, report_path)

cat("Bootstrap variation analysis complete\n")
cat("Report:", report_path, "\n")
cat("Figures:", figure_dir, "\n")
cat("Tables:", table_dir, "\n")
