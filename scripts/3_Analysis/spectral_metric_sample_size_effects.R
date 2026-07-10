user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(unique(c(user_library, .libPaths())))
}

required_packages <- c("terra", "dplyr", "tidyr", "readr", "ggplot2", "scales", "knitr", "alphahull")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(knitr)
})

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

Sys.setenv(RUN_SPECTRAL_HET_WORKFLOW = "false")
source("scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R")

REPORT_PATH <- "reports/analysis/20260704_pca_metric_sample_size_effects.md"
SA_DESIGN_PATH <- "reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_design.csv"
PCA_BASIS <- Sys.getenv("PCA_BASIS", "regular_PCA")

if (!PCA_BASIS %in% c("regular_PCA", "standardized_PCA")) {
  stop("PCA_BASIS must be either regular_PCA or standardized_PCA", call. = FALSE)
}

PCA_BASIS_LABEL <- if (PCA_BASIS == "standardized_PCA") "Standardized PCA" else "Regular PCA"
PCA_BASIS_RDS <- if (PCA_BASIS == "standardized_PCA") STANDARDIZED_PCA_RDS else PCA_RDS

if (PCA_BASIS == "standardized_PCA") {
  REPORT_PATH <- "reports/analysis/20260707_standardized_pca_metric_sample_size_effects.md"
}

N_BOOT_EXPERIMENT <- 50
EXPERIMENT_SEED <- 20260703

metric_catalog <- data.frame(
  base_metric_id = c("pca_mean_distance", "spectral_rao_q", "alpha_hull_area"),
  metric_id = c("pca_mean_distance", "spectral_rao_q", "alpha_hull_area"),
  metric_label = c("PCA Mean Distance", "Spectral Rao's Q", "Alpha-Hull Area"),
  metric_column = c("pca_euclidean_mean", "rao_q_pca", "alpha_hull_area"),
  method_column = c("metric_method", "rao_q_distance", "alpha_hull_method"),
  n_points_column = c("n_pixels", "n_pixels", "alpha_hull_n_points"),
  y_label = c(
    "Mean Euclidean distance in PCA space",
    "Rao's Q in PCA space",
    "Alpha-hull area in PC1-PC2 space"
  ),
  stringsAsFactors = FALSE
)

if (PCA_BASIS == "standardized_PCA") {
  metric_catalog$metric_id <- paste0("standardized_PCA_", metric_catalog$metric_id)
  metric_catalog$metric_label <- paste("Standardized", metric_catalog$metric_label)
  metric_catalog$y_label <- sub("PCA", "standardized PCA", metric_catalog$y_label, fixed = TRUE)
}

requested_metrics <- Sys.getenv("SAMPLE_SIZE_METRICS", "")
if (nzchar(requested_metrics)) {
  requested_metrics <- trimws(strsplit(requested_metrics, ",", fixed = TRUE)[[1]])
  requested_metric_ids <- unique(c(
    requested_metrics[requested_metrics %in% metric_catalog$metric_id],
    metric_catalog$metric_id[metric_catalog$base_metric_id %in% requested_metrics]
  ))
  unknown_metrics <- setdiff(requested_metrics, c(metric_catalog$metric_id, metric_catalog$base_metric_id))
  if (length(unknown_metrics) > 0) {
    stop("Unknown SAMPLE_SIZE_METRICS value(s): ", paste(unknown_metrics, collapse = ", "), call. = FALSE)
  }
  metric_catalog <- metric_catalog |>
    filter(.data$metric_id %in% requested_metric_ids)
}

dir.create(dirname(REPORT_PATH), recursive = TRUE, showWarnings = FALSE)

format_pct <- function(x, accuracy = 0.1) {
  scales::percent(x, accuracy = accuracy)
}

format_number <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

seed_from_values <- function(...) {
  text <- paste(..., collapse = "_")
  EXPERIMENT_SEED + sum(utf8ToInt(text))
}

find_raster_file <- function(scale, quad_id) {
  scale_spec_dir <- c(
    "10m" = "Quad_Spectra/10m_smooth_5nm",
    "20m" = "Quad_Spectra/20m_smooth_5nm",
    "50m" = "Quad_Spectra/50m_smooth_5nm"
  )
  file <- file.path(scale_spec_dir[[scale]], quad_id)
  if (!file.exists(file)) {
    stop("Missing raster file for ", scale, " quadrat ", quad_id, ": ", file, call. = FALSE)
  }
  file
}

calculate_selected_metrics <- function(scores, quad_id) {
  n_pixels <- nrow(scores)
  run_alpha <- "alpha_hull_area" %in% metric_catalog$base_metric_id

  if (n_pixels < MIN_VALID_PIXELS) {
    return(data.frame(
      n_pixels = n_pixels,
      pca_euclidean_mean = NA_real_,
      metric_method = "insufficient_pixels",
      rao_q_pca = NA_real_,
      rao_q_distance = "squared_euclidean_pc1_pc3",
      alpha_hull_area = NA_real_,
      alpha_hull_method = "insufficient_pixels",
      alpha_hull_n_points = n_pixels,
      stringsAsFactors = FALSE
    ))
  }

  center <- colMeans(scores)
  centered <- sweep(scores, 2, center, "-")
  squared_radius <- rowSums(centered^2)
  distances <- sqrt(squared_radius)
  alpha_summary <- if (run_alpha) {
    alpha_hull_area_2d(scores[, 1:2, drop = FALSE], quad_id)
  } else {
    data.frame(
      alpha_hull_area = NA_real_,
      alpha_hull_method = "not_requested",
      alpha_hull_n_points = n_pixels,
      stringsAsFactors = FALSE
    )
  }

  data.frame(
    n_pixels = n_pixels,
    pca_euclidean_mean = mean(distances),
    metric_method = "all_sampled_pixels",
    rao_q_pca = 2 * mean(squared_radius),
    rao_q_distance = "squared_euclidean_pc1_pc3",
    alpha_hull_area = alpha_summary$alpha_hull_area,
    alpha_hull_method = alpha_summary$alpha_hull_method,
    alpha_hull_n_points = alpha_summary$alpha_hull_n_points,
    stringsAsFactors = FALSE
  )
}

long_metric_rows <- function(metric_values, metadata) {
  bind_rows(lapply(seq_len(nrow(metric_catalog)), function(i) {
    metric <- metric_catalog[i, ]
    data.frame(
      metadata,
      metric_id = metric$metric_id,
      metric_label = metric$metric_label,
      metric_value = metric_values[[metric$metric_column]],
      metric_method = metric_values[[metric$method_column]],
      metric_n_points = metric_values[[metric$n_points_column]],
      stringsAsFactors = FALSE
    )
  }))
}

run_condition <- function(scores, scale, quad_id, condition) {
  n_pixels <- nrow(scores)
  sample_size <- condition$sample_size
  set.seed(seed_from_values(scale, quad_id, condition$sample_rule))
  full_pixel_values <- NULL

  if (sample_size == n_pixels) {
    full_pixel_values <- calculate_selected_metrics(scores, quad_id)
  }

  bind_rows(lapply(seq_len(N_BOOT_EXPERIMENT), function(iteration) {
    if (sample_size == n_pixels) {
      metric_values <- full_pixel_values
    } else {
      sample_index <- sample.int(n_pixels, sample_size, replace = FALSE)
      metric_values <- calculate_selected_metrics(scores[sample_index, , drop = FALSE], quad_id)
    }

    metadata <- data.frame(
      scale = scale,
      quad_id = quad_id,
      n_pixels = n_pixels,
      sample_rule = condition$sample_rule,
      rule_type = condition$rule_type,
      rule_value = condition$rule_value,
      rule_label = condition$rule_label,
      rule_axis_label = condition$rule_axis_label,
      sample_size = sample_size,
      sample_fraction = sample_size / n_pixels,
      bootstrap_iter = iteration,
      stringsAsFactors = FALSE
    )

    long_metric_rows(metric_values, metadata)
  }))
}

make_summary_table <- function(bootstrap_long, design_table) {
  fixed_reference <- bootstrap_long |>
    filter(sample_rule == "fixed_4000") |>
    group_by(metric_id, scale, quad_id) |>
    summarise(fixed_4000_mean = mean(metric_value, na.rm = TRUE), .groups = "drop")

  bootstrap_long |>
    group_by(metric_id, metric_label, scale, quad_id, sample_rule, sample_size, sample_fraction) |>
    summarise(
      n_boot = sum(is.finite(metric_value)),
      metric_mean = ifelse(n_boot > 0, mean(metric_value, na.rm = TRUE), NA_real_),
      metric_median = ifelse(n_boot > 0, median(metric_value, na.rm = TRUE), NA_real_),
      metric_sd = ifelse(n_boot > 1, sd(metric_value, na.rm = TRUE), NA_real_),
      metric_min = ifelse(n_boot > 0, min(metric_value, na.rm = TRUE), NA_real_),
      metric_max = ifelse(n_boot > 0, max(metric_value, na.rm = TRUE), NA_real_),
      metric_method = paste(sort(unique(metric_method)), collapse = "; "),
      metric_n_points_min = ifelse(any(is.finite(metric_n_points)), min(metric_n_points, na.rm = TRUE), NA_real_),
      metric_n_points_max = ifelse(any(is.finite(metric_n_points)), max(metric_n_points, na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) |>
    mutate(
      metric_cv = ifelse(is.finite(metric_mean) & metric_mean != 0, metric_sd / abs(metric_mean), NA_real_),
      metric_se = metric_sd / sqrt(n_boot),
      ci_low = metric_mean - stats::qt(0.975, pmax(n_boot - 1, 1)) * metric_se,
      ci_high = metric_mean + stats::qt(0.975, pmax(n_boot - 1, 1)) * metric_se
    ) |>
    left_join(fixed_reference, by = c("metric_id", "scale", "quad_id")) |>
    mutate(delta_from_fixed_4000 = metric_mean - fixed_4000_mean) |>
    left_join(
      design_table |>
        select(
          scale, quad_id, sample_rule, rule_type, rule_value, rule_label,
          rule_axis_label, sample_label, plot_sample_label, source_n_pixels,
          source_method, selection_note
        ),
      by = c("scale", "quad_id", "sample_rule")
    ) |>
    arrange(metric_id, factor(scale, levels = c("10m", "20m", "50m")), quad_id, sample_fraction, sample_rule)
}

write_metric_tables <- function(metric_id, design_table, bootstrap_long, summary_table) {
  output_dir <- file.path("reports/tables/sample_size_effects", metric_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  metric_boot <- bootstrap_long |>
    filter(metric_id == !!metric_id) |>
    select(-metric_id, -metric_label)
  metric_summary <- summary_table |>
    filter(metric_id == !!metric_id) |>
    select(-metric_id, -metric_label)

  write_csv(design_table, file.path(output_dir, paste0(metric_id, "_sample_size_design.csv")))
  write_csv(metric_boot, file.path(output_dir, paste0(metric_id, "_sample_size_boot_long.csv")))
  write_csv(metric_summary, file.path(output_dir, paste0(metric_id, "_sample_size_summary.csv")))
}

make_metric_figures <- function(metric_id, bootstrap_long, summary_table, design_table) {
  metric <- metric_catalog |> filter(.data$metric_id == !!metric_id)
  metric_boot <- bootstrap_long |> filter(.data$metric_id == !!metric_id)
  metric_summary <- summary_table |> filter(.data$metric_id == !!metric_id)

  figure_dir <- file.path("reports/figures/sample_size_effects", metric_id)
  split_dir <- file.path(figure_dir, "distributions_by_sample_size")
  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)

  if (!any(is.finite(metric_boot$metric_value)) || !any(is.finite(metric_summary$metric_mean))) {
    writeLines(
      c(
        paste(metric$metric_label, "figures not generated."),
        "No finite metric values were available for this PCA basis and metric."
      ),
      file.path(figure_dir, paste0(metric_id, "_figure_skip_note.txt"))
    )
    return(invisible(FALSE))
  }

  line_color <- "#3b6f6d"
  point_color <- "#273947"
  box_fill <- "#cfe3dc"
  box_outline <- "#34515e"

  theme_report <- theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )

  metric_summary <- metric_summary |>
    mutate(
      scale = factor(scale, levels = c("10m", "20m", "50m")),
      quad_axis_label = paste0(quad_id, "\nn = ", comma(source_n_pixels), " px"),
      rule_axis_label = factor(rule_axis_label, levels = unique(design_table$rule_axis_label))
    )

  metric_boot <- metric_boot |>
    mutate(
      scale = factor(scale, levels = c("10m", "20m", "50m")),
      quad_axis_label = paste0(quad_id, "\nn = ", comma(source_n_pixels), " px"),
      rule_axis_label = factor(rule_axis_label, levels = unique(design_table$rule_axis_label))
    )

  p_mean_overview <- ggplot(metric_summary, aes(x = sample_fraction, y = metric_mean, group = quad_id)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.002, alpha = 0.25, color = line_color) +
    geom_line(alpha = 0.35, linewidth = 0.45, color = line_color) +
    geom_point(alpha = 0.75, size = 1.6, color = point_color) +
    facet_wrap(~ scale, scales = "free_x") +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste(metric$metric_label, "mean by sample size"),
      subtitle = "Compact overview; scale-specific figures show each quadrat separately.",
      x = "Sampled fraction of retained pixels",
      y = metric$y_label
    ) +
    theme_report
  ggsave(file.path(figure_dir, paste0(metric_id, "_mean_by_sample_size.png")), p_mean_overview, width = 12, height = 6, dpi = 320, bg = "white")

  p_cv_overview <- ggplot(metric_summary, aes(x = sample_fraction, y = metric_cv, group = quad_id)) +
    geom_line(alpha = 0.35, linewidth = 0.45, color = line_color) +
    geom_point(alpha = 0.75, size = 1.6, color = point_color) +
    facet_wrap(~ scale, scales = "free_x") +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      title = paste(metric$metric_label, "bootstrap CV by sample size"),
      subtitle = "Lower values indicate more stable bootstrap estimates for a quadrat.",
      x = "Sampled fraction of retained pixels",
      y = "Bootstrap coefficient of variation"
    ) +
    theme_report
  ggsave(file.path(figure_dir, paste0(metric_id, "_cv_by_sample_size.png")), p_cv_overview, width = 12, height = 6, dpi = 320, bg = "white")

  p_delta_overview <- ggplot(metric_summary, aes(x = sample_fraction, y = delta_from_fixed_4000, group = quad_id)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_line(alpha = 0.35, linewidth = 0.45, color = line_color) +
    geom_point(alpha = 0.75, size = 1.6, color = point_color) +
    facet_wrap(~ scale, scales = "free_x") +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Difference from the fixed 4,000-pixel rule",
      subtitle = "Compact overview; scale-specific figures show each quadrat separately.",
      x = "Sampled fraction of retained pixels",
      y = paste("Mean difference in", metric$metric_label)
    ) +
    theme_report
  ggsave(file.path(figure_dir, paste0(metric_id, "_delta_from_fixed_4000.png")), p_delta_overview, width = 12, height = 6, dpi = 320, bg = "white")

  p_reps_overview <- ggplot(metric_boot, aes(x = rule_axis_label, y = metric_value)) +
    geom_boxplot(outlier.alpha = 0.12, linewidth = 0.25, fill = box_fill, color = box_outline) +
    facet_wrap(~ scale, scales = "free_y") +
    scale_x_discrete(drop = TRUE) +
    labs(
      title = paste(metric$metric_label, "bootstrap replicate distributions by sample-size rule"),
      subtitle = "Compact overview across all selected quadrats; scale-specific figures separate quadrats.",
      x = "Sample-size rule",
      y = metric$y_label
    ) +
    theme_report +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
  ggsave(file.path(figure_dir, paste0(metric_id, "_replicate_distributions.png")), p_reps_overview, width = 12, height = 6, dpi = 320, bg = "white")

  for (scale_name in c("10m", "20m", "50m")) {
    scale_summary <- metric_summary |> filter(scale == scale_name)

    p_mean_scale <- ggplot(scale_summary, aes(x = sample_fraction, y = metric_mean, group = quad_id)) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.002, alpha = 0.7, color = line_color) +
      geom_line(alpha = 0.55, linewidth = 0.5, color = line_color) +
      geom_point(alpha = 0.8, size = 1.7, color = point_color) +
      facet_wrap(~ quad_axis_label, scales = "free_y", ncol = 4) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        title = paste(scale_name, metric$metric_label, "mean by sample size"),
        x = "Sampled fraction of retained pixels",
        y = metric$y_label
      ) +
      theme_report +
      theme(strip.text = element_text(size = 8))
    ggsave(file.path(figure_dir, paste0(metric_id, "_mean_by_sample_size_", scale_name, ".png")), p_mean_scale, width = 15, height = 10, dpi = 320, bg = "white")

    p_cv_scale <- ggplot(scale_summary, aes(x = sample_fraction, y = metric_cv, group = quad_id)) +
      geom_line(alpha = 0.55, linewidth = 0.5, color = line_color) +
      geom_point(alpha = 0.8, size = 1.7, color = point_color) +
      facet_wrap(~ quad_axis_label, scales = "free_y", ncol = 4) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
      labs(
        title = paste(scale_name, metric$metric_label, "bootstrap CV by sample size"),
        x = "Sampled fraction of retained pixels",
        y = "Bootstrap coefficient of variation"
      ) +
      theme_report +
      theme(strip.text = element_text(size = 8))
    ggsave(file.path(figure_dir, paste0(metric_id, "_cv_by_sample_size_", scale_name, ".png")), p_cv_scale, width = 15, height = 10, dpi = 320, bg = "white")

    p_delta_scale <- ggplot(scale_summary, aes(x = sample_fraction, y = delta_from_fixed_4000, group = quad_id)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_line(alpha = 0.55, linewidth = 0.5, color = line_color) +
      geom_point(alpha = 0.8, size = 1.7, color = point_color) +
      facet_wrap(~ quad_axis_label, scales = "free_y", ncol = 4) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        title = paste(scale_name, metric$metric_label, "difference from fixed 4,000 pixels"),
        x = "Sampled fraction of retained pixels",
        y = paste("Mean difference in", metric$metric_label)
      ) +
      theme_report +
      theme(strip.text = element_text(size = 8))
    ggsave(file.path(figure_dir, paste0(metric_id, "_delta_from_fixed_4000_", scale_name, ".png")), p_delta_scale, width = 15, height = 10, dpi = 320, bg = "white")

    rules_for_scale <- design_table |>
      filter(scale == scale_name) |>
      distinct(sample_rule, rule_axis_label) |>
      mutate(rule_order = match(sample_rule, unique(design_table$sample_rule))) |>
      arrange(rule_order)

    for (rule_i in seq_len(nrow(rules_for_scale))) {
      rule_row <- rules_for_scale[rule_i, ]
      split_data <- metric_boot |>
        filter(scale == scale_name, sample_rule == rule_row$sample_rule) |>
        mutate(quad_axis_label = factor(quad_axis_label, levels = unique(quad_axis_label)))

      p_split <- ggplot(split_data, aes(x = quad_axis_label, y = metric_value)) +
        geom_boxplot(outlier.alpha = 0.14, linewidth = 0.25, fill = box_fill, color = box_outline) +
        labs(
          title = paste(scale_name, metric$metric_label, "bootstrap distributions:", rule_row$rule_axis_label),
          subtitle = "Each box is one quadrat; labels include retained pixel counts.",
          x = "Quadrat",
          y = metric$y_label
        ) +
        theme_report +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7))

      ggsave(
        file.path(split_dir, paste0(metric_id, "_replicate_distribution_", scale_name, "_", rule_row$sample_rule, ".png")),
        p_split,
        width = 18,
        height = 6.5,
        dpi = 320,
        bg = "white"
      )
    }
  }
}

make_metric_mean_figures <- function(metric_id, summary_table, design_table) {
  metric <- metric_catalog |> filter(.data$metric_id == !!metric_id)
  metric_summary <- summary_table |> filter(.data$metric_id == !!metric_id)

  figure_dir <- file.path("reports/figures/sample_size_effects", metric_id)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  if (!any(is.finite(metric_summary$metric_mean))) {
    writeLines(
      c(
        paste(metric$metric_label, "mean figures not generated."),
        "No finite metric means were available for this PCA basis and metric."
      ),
      file.path(figure_dir, paste0(metric_id, "_mean_figure_skip_note.txt"))
    )
    return(invisible(FALSE))
  }

  line_color <- "#3b6f6d"
  point_color <- "#273947"

  theme_report <- theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )

  metric_summary <- metric_summary |>
    mutate(
      scale = factor(scale, levels = c("10m", "20m", "50m")),
      quad_axis_label = paste0(quad_id, "\nn = ", comma(source_n_pixels), " px"),
      rule_axis_label = factor(rule_axis_label, levels = unique(design_table$rule_axis_label))
    )

  p_mean_overview <- ggplot(metric_summary, aes(x = sample_fraction, y = metric_mean, group = quad_id)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.002, alpha = 0.25, color = line_color) +
    geom_line(alpha = 0.35, linewidth = 0.45, color = line_color) +
    geom_point(alpha = 0.75, size = 1.6, color = point_color) +
    facet_wrap(~ scale, scales = "free_x") +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste(metric$metric_label, "mean by sample size"),
      subtitle = "Error bars show 95% confidence intervals around the 50-iteration mean.",
      x = "Sampled fraction of retained pixels",
      y = metric$y_label
    ) +
    theme_report
  ggsave(file.path(figure_dir, paste0(metric_id, "_mean_by_sample_size.png")), p_mean_overview, width = 12, height = 6, dpi = 320, bg = "white")

  for (scale_name in c("10m", "20m", "50m")) {
    scale_summary <- metric_summary |> filter(scale == scale_name)

    p_mean_scale <- ggplot(scale_summary, aes(x = sample_fraction, y = metric_mean, group = quad_id)) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.002, alpha = 0.7, color = line_color) +
      geom_line(alpha = 0.55, linewidth = 0.5, color = line_color) +
      geom_point(alpha = 0.8, size = 1.7, color = point_color) +
      facet_wrap(~ quad_axis_label, scales = "free_y", ncol = 4) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      labs(
        title = paste(scale_name, metric$metric_label, "mean by sample size"),
        subtitle = "Error bars show 95% confidence intervals around the 50-iteration mean.",
        x = "Sampled fraction of retained pixels",
        y = metric$y_label
      ) +
      theme_report +
      theme(strip.text = element_text(size = 8))
    ggsave(file.path(figure_dir, paste0(metric_id, "_mean_by_sample_size_", scale_name, ".png")), p_mean_scale, width = 15, height = 10, dpi = 320, bg = "white")
  }
}

write_report <- function(design_table, summary_table) {
  selection_report <- design_table |>
    distinct(scale, quad_id, source_n_pixels, selection_note) |>
    transmute(
      Scale = scale,
      Quadrat = quad_id,
      `Retained pixels` = comma(source_n_pixels),
      `Selection note` = selection_note
    )

  design_report <- design_table |>
    transmute(
      Scale = scale,
      Quadrat = quad_id,
      `Retained pixels` = comma(source_n_pixels),
      `Sample rule` = sample_label,
      `Sample pixels` = comma(sample_size),
      `Sample percent` = format_pct(sample_fraction)
    )

  summary_report <- summary_table |>
    transmute(
      Metric = metric_label,
      Scale = as.character(scale),
      Quadrat = quad_id,
      `Sample rule` = sample_label,
      `Mean value` = format_number(metric_mean),
      `Boot SD` = format_number(metric_sd),
      `Boot CV` = format_pct(metric_cv),
      `95% CI low` = format_number(ci_low),
      `95% CI high` = format_number(ci_high),
      `Delta from fixed 4,000` = format_number(delta_from_fixed_4000),
      `Calculation method` = metric_method
    )

  metric_sections <- unlist(lapply(seq_len(nrow(metric_catalog)), function(i) {
    metric <- metric_catalog[i, ]
    scale_lines <- unlist(lapply(c("10m", "20m", "50m"), function(scale_name) {
      rules_for_scale <- design_table |>
        filter(scale == scale_name) |>
        distinct(sample_rule, rule_axis_label) |>
        mutate(rule_order = match(sample_rule, unique(design_table$sample_rule))) |>
        arrange(rule_order)

      c(
        paste0("### ", metric$metric_label, ": ", scale_name, " Scale Figures"),
        "",
        paste0(
          "![", metric$metric_label, " ", scale_name, " mean](../figures/sample_size_effects/",
          metric$metric_id, "/", metric$metric_id, "_mean_by_sample_size_", scale_name, ".png)"
        ),
        "",
        paste0(
          "![", metric$metric_label, " ", scale_name, " CV](../figures/sample_size_effects/",
          metric$metric_id, "/", metric$metric_id, "_cv_by_sample_size_", scale_name, ".png)"
        ),
        "",
        paste0(
          "![", metric$metric_label, " ", scale_name, " delta](../figures/sample_size_effects/",
          metric$metric_id, "/", metric$metric_id, "_delta_from_fixed_4000_", scale_name, ".png)"
        ),
        "",
        paste0("### ", metric$metric_label, ": ", scale_name, " Distribution Charts Split By Sample Size"),
        "",
        unlist(lapply(seq_len(nrow(rules_for_scale)), function(rule_i) {
          rule_row <- rules_for_scale[rule_i, ]
          paste0(
            "- [", scale_name, " ", rule_row$rule_axis_label,
            "](../figures/sample_size_effects/", metric$metric_id,
            "/distributions_by_sample_size/", metric$metric_id,
            "_replicate_distribution_", scale_name, "_", rule_row$sample_rule, ".png)"
          )
        })),
        ""
      )
    }))

    c(
      paste0("## ", metric$metric_label),
      "",
      paste0(
        "![", metric$metric_label, " mean](../figures/sample_size_effects/",
        metric$metric_id, "/", metric$metric_id, "_mean_by_sample_size.png)"
      ),
      "",
      paste0(
        "![", metric$metric_label, " CV](../figures/sample_size_effects/",
        metric$metric_id, "/", metric$metric_id, "_cv_by_sample_size.png)"
      ),
      "",
      paste0(
        "![", metric$metric_label, " delta](../figures/sample_size_effects/",
        metric$metric_id, "/", metric$metric_id, "_delta_from_fixed_4000.png)"
      ),
      "",
      paste0(
        "![", metric$metric_label, " replicate distributions](../figures/sample_size_effects/",
        metric$metric_id, "/", metric$metric_id, "_replicate_distributions.png)"
      ),
      "",
      scale_lines,
      "",
      "Output tables:",
      "",
      paste0("- `reports/tables/sample_size_effects/", metric$metric_id, "/", metric$metric_id, "_sample_size_design.csv`"),
      paste0("- `reports/tables/sample_size_effects/", metric$metric_id, "/", metric$metric_id, "_sample_size_boot_long.csv`"),
      paste0("- `reports/tables/sample_size_effects/", metric$metric_id, "/", metric$metric_id, "_sample_size_summary.csv`"),
      ""
    )
  }))

  report_lines <- c(
    paste0("# ", PCA_BASIS_LABEL, "-Derived Spectral Metric Sample-Size Effects"),
    "",
    if (PCA_BASIS == "standardized_PCA") "Date: 2026-07-07" else "Date: 2026-07-04",
    "",
    "## Purpose",
    "",
    paste0(
      "Extend the sample-size sensitivity analysis from spectral angle entropy to ",
      if (PCA_BASIS == "standardized_PCA") "vector-normalized standardized PCA" else "regular PCA",
      " mean distance, spectral Rao's Q, and alpha-hull area."
    ),
    "",
    "## Design",
    "",
    paste0("- Bootstrap iterations per quadrat x sample rule: ", N_BOOT_EXPERIMENT),
    paste0("- PCA basis: ", PCA_BASIS_LABEL, " (`", PCA_BASIS_RDS, "`)."),
    "- Quadrat and sample-size rules are reused from `reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_design.csv` so the non-SA metrics are evaluated on the same quadrats and retained-pixel sample sizes.",
    "- Each replicate samples retained pixels without replacement. If a rule resolves to 100% of retained pixels, the metric is calculated once from the full retained-pixel set and repeated across the 50 output rows, so full-pixel conditions have zero artificial bootstrap variation.",
    "- PCA mean distance and spectral Rao's Q use all sampled pixels in PCA space. Alpha-hull area uses all sampled PC1-PC2 points unless the existing alpha-hull helper has to remove duplicate points or returns an internal failure method.",
    "",
    "## Selected Quadrats",
    "",
    paste(capture.output(knitr::kable(selection_report, format = "markdown")), collapse = "\n"),
    "",
    "## Sample Rules",
    "",
    paste(capture.output(knitr::kable(design_report, format = "markdown")), collapse = "\n"),
    "",
    "## Summary Results",
    "",
    paste(capture.output(knitr::kable(summary_report, format = "markdown")), collapse = "\n"),
    "",
    "## Figures",
    "",
    "Mean-by-sample-size figures include 95 percent confidence interval error bars around the 50-iteration mean. Compact overview distribution plots are retained for each metric. Longer distribution plots are split by sample-size rule, with quadrats as boxes and retained pixel counts in the x-axis labels.",
    "",
    metric_sections
  )

  writeLines(report_lines, REPORT_PATH)
}

if (!file.exists(SA_DESIGN_PATH)) {
  stop("Missing SA design table: ", SA_DESIGN_PATH, call. = FALSE)
}

design_table <- read_csv(SA_DESIGN_PATH, show_col_types = FALSE) |>
  mutate(
    scale = factor(scale, levels = c("10m", "20m", "50m")),
    quad_id = as.character(quad_id),
    sample_rule = as.character(sample_rule)
  ) |>
  arrange(scale, quad_id, sample_fraction, sample_rule)

quad_metadata <- design_table |>
  distinct(scale, quad_id, source_n_pixels, selection_note, raster_file) |>
  mutate(scale = as.character(scale)) |>
  arrange(factor(scale, levels = c("10m", "20m", "50m")), quad_id)

if (identical(tolower(Sys.getenv("REPORT_ONLY", "false")), "true")) {
  summary_table <- bind_rows(lapply(seq_len(nrow(metric_catalog)), function(i) {
    metric <- metric_catalog[i, ]
    summary_path <- file.path(
      "reports/tables/sample_size_effects",
      metric$metric_id,
      paste0(metric$metric_id, "_sample_size_summary.csv")
    )
    if (!file.exists(summary_path)) {
      stop("Missing metric summary table for report-only mode: ", summary_path, call. = FALSE)
    }
    read_csv(summary_path, show_col_types = FALSE) |>
      mutate(
        metric_id = metric$metric_id,
        metric_label = metric$metric_label
      )
  }))

  write_report(design_table, summary_table)
  cat("Combined report updated:", REPORT_PATH, "\n")
  quit(save = "no")
}

if (identical(tolower(Sys.getenv("MEAN_FIGURES_ONLY", "false")), "true")) {
  summary_table <- bind_rows(lapply(seq_len(nrow(metric_catalog)), function(i) {
    metric <- metric_catalog[i, ]
    summary_path <- file.path(
      "reports/tables/sample_size_effects",
      metric$metric_id,
      paste0(metric$metric_id, "_sample_size_summary.csv")
    )
    if (!file.exists(summary_path)) {
      stop("Missing metric summary table for mean-figures-only mode: ", summary_path, call. = FALSE)
    }
    read_csv(summary_path, show_col_types = FALSE) |>
      mutate(
        metric_id = metric$metric_id,
        metric_label = metric$metric_label
      )
  }))

  for (metric_id in metric_catalog$metric_id) {
    make_metric_mean_figures(metric_id, summary_table, design_table)
  }

  cat("Mean figures updated from existing summary tables\n")
  cat("Metric figures root: reports/figures/sample_size_effects\n")
  quit(save = "no")
}

pca_object <- readRDS(PCA_BASIS_RDS)

bootstrap_long <- bind_rows(lapply(seq_len(nrow(quad_metadata)), function(i) {
  quad <- quad_metadata[i, ]
  message("Processing ", quad$scale, " ", quad$quad_id, " (", i, "/", nrow(quad_metadata), ")")

  raster_file <- find_raster_file(quad$scale, quad$quad_id)
  spectra <- read_masked_spectra(raster_file)
  scores <- project_scores(spectra, pca_object)

  conditions <- design_table |>
    filter(scale == quad$scale, quad_id == quad$quad_id) |>
    arrange(sample_fraction, sample_rule)

  bind_rows(lapply(seq_len(nrow(conditions)), function(condition_i) {
    run_condition(scores, quad$scale, quad$quad_id, conditions[condition_i, ])
  }))
}))

bootstrap_long <- bootstrap_long |>
  mutate(scale = factor(scale, levels = c("10m", "20m", "50m"))) |>
  left_join(
    design_table |>
      select(
        scale, quad_id, sample_rule, selection_note, source_n_pixels,
        source_method, sample_label, plot_sample_label
      ),
    by = c("scale", "quad_id", "sample_rule")
  ) |>
  arrange(metric_id, scale, quad_id, sample_fraction, sample_rule, bootstrap_iter)

summary_table <- make_summary_table(bootstrap_long, design_table)

for (metric_id in metric_catalog$metric_id) {
  write_metric_tables(metric_id, design_table, bootstrap_long, summary_table)
  make_metric_figures(metric_id, bootstrap_long, summary_table, design_table)
}

if (!identical(tolower(Sys.getenv("WRITE_SAMPLE_SIZE_REPORT", "true")), "false")) {
  write_report(design_table, summary_table)
}

cat("PCA-derived spectral metric sample-size effects complete\n")
cat("Report:", REPORT_PATH, "\n")
cat("Metric tables root: reports/tables/sample_size_effects\n")
cat("Metric figures root: reports/figures/sample_size_effects\n")
