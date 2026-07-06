user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(unique(c(user_library, .libPaths())))
}

required_packages <- c("terra", "dplyr", "tidyr", "readr", "ggplot2", "scales", "knitr")
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

Sys.setenv(RUN_SA_ENTROPY_WORKFLOW = "false")
source("scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R")

VARIATION_TYPE <- "sa_entropy"
OUTPUT_DIR <- file.path("reports/tables/sample_size_effects", VARIATION_TYPE)
FIGURE_DIR <- file.path("reports/figures/sample_size_effects", VARIATION_TYPE)
REPORT_PATH <- "reports/analysis/20260703_sa_entropy_sample_size_effects.md"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURE_DIR, "distributions_by_sample_size"), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REPORT_PATH), recursive = TRUE, showWarnings = FALSE)

N_BOOT_EXPERIMENT <- 50
MAX_SAMPLE_PIXELS <- 5000
PAIR_SAMPLE_EXPERIMENT <- 10000
EXPERIMENT_SEED <- 20260703

target_quads_per_scale <- c("10m" = 32L, "20m" = 16L, "50m" = 8L)

seed_quads <- data.frame(
  scale = c("10m", "10m", "20m", "20m", "50m", "50m"),
  quad_id = c("800_a", "112_b", "219", "302", "sub50_36", "sub50_60"),
  selection_note = c(
    "Original 10 m replacement: closest bootstrap-mean quadrat at or below 4,000 retained pixels",
    "Original random bootstrap-mean quadrat selected with seed 20260703",
    "Original random bootstrap-mean quadrat selected with seed 20260703",
    "Original random bootstrap-mean quadrat selected with seed 20260703",
    "Original random bootstrap-mean quadrat selected with seed 20260703",
    "Original random bootstrap-mean quadrat selected with seed 20260703"
  ),
  stringsAsFactors = FALSE
)

sample_rules <- data.frame(
  sample_rule = c(
    "pct_1", "pct_2", "pct_3", "fixed_1250",
    "fixed_2000", "fixed_3000", "fixed_4000", "fixed_6000", "fixed_8000"
  ),
  rule_type = c(
    "percent", "percent", "percent", "fixed",
    "fixed", "fixed", "fixed", "fixed", "fixed"
  ),
  rule_value = c(0.01, 0.02, 0.03, 1250, 2000, 3000, 4000, 6000, 8000),
  rule_label = c(
    "1% of pixels", "2% of pixels", "3% of pixels", "Fixed 1,250",
    "Fixed 2,000", "Fixed 3,000", "Fixed 4,000", "Fixed 6,000", "Fixed 8,000"
  ),
  rule_axis_label = c(
    "1%", "2%", "3%", "1,250",
    "2,000", "3,000", "4,000", "6,000", "8,000"
  ),
  applicable_scales = c(
    "all", "all", "all", "all",
    "10m,20m", "10m,20m", "all", "50m", "50m"
  ),
  stringsAsFactors = FALSE
)

scale_spec_dir <- c(
  "10m" = "Quad_Spectra/10m_smooth_5nm",
  "20m" = "Quad_Spectra/20m_smooth_5nm",
  "50m" = "Quad_Spectra/50m_smooth_5nm"
)

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

rule_applies_to_scale <- function(scale_name, applicable_scales) {
  applicable_scales == "all" || scale_name %in% strsplit(applicable_scales, ",", fixed = TRUE)[[1]]
}

read_sa_summary <- function(scale_name) {
  read_csv(
    file.path("Quad_Values", paste0(scale_name, "_SA_entropy_smooth_masked_5nm_summary.csv")),
    show_col_types = FALSE
  ) |>
    mutate(scale = scale_name, quad_id = as.character(quad_id))
}

select_additional_quads <- function(scale_name, scale_summary) {
  scale_seed_quads <- seed_quads |>
    filter(scale == scale_name)
  n_additional <- target_quads_per_scale[[scale_name]] - nrow(scale_seed_quads)

  if (n_additional < 0) {
    stop("More seed quadrats than requested for ", scale_name, call. = FALSE)
  }

  candidates <- scale_summary |>
    filter(
      method == "bootstrap_mean",
      is.finite(spectral_entropy),
      !quad_id %in% scale_seed_quads$quad_id
    ) |>
    arrange(quad_id)

  if (nrow(candidates) < n_additional) {
    stop("Not enough candidate quadrats for ", scale_name, call. = FALSE)
  }

  set.seed(seed_from_values("select", scale_name))
  sampled <- candidates[sample(seq_len(nrow(candidates)), n_additional), , drop = FALSE]

  bind_rows(
    scale_seed_quads,
    sampled |>
      transmute(
        scale,
        quad_id,
        selection_note = paste0(
          "Additional random bootstrap-mean quadrat selected with seed ",
          EXPERIMENT_SEED
        )
      )
  ) |>
    mutate(selection_order = row_number())
}

find_raster_file <- function(scale, quad_id) {
  file <- file.path(scale_spec_dir[[scale]], quad_id)
  if (!file.exists(file)) {
    stop("Missing raster file for ", scale, " quadrat ", quad_id, ": ", file, call. = FALSE)
  }
  file
}

read_quadrat_spectra_norm <- function(file) {
  raster <- terra::rast(file)
  raster_masked <- apply_shadow_mask(raster)
  normalize_spectra(terra::values(raster_masked, mat = TRUE))
}

resolve_sample_size <- function(n_pixels, rule_type, rule_value) {
  if (rule_type == "fixed") {
    return(min(as.integer(rule_value), n_pixels))
  }

  min(
    MAX_SAMPLE_PIXELS,
    max(1L, as.integer(round(n_pixels * rule_value)))
  )
}

entropy_for_sample <- function(sampled_norm) {
  possible_pairs <- pair_count(nrow(sampled_norm))
  if (possible_pairs <= MAX_BOOT_EXACT_PAIRS) {
    return(list(
      entropy = spectral_entropy_from_norm(sampled_norm),
      pair_method = "all_pairs_within_sample",
      pair_count_used = possible_pairs
    ))
  }

  list(
    entropy = spectral_entropy_from_sampled_pairs(
      sampled_norm,
      n_pairs = PAIR_SAMPLE_EXPERIMENT,
      n_bins = N_BINS
    ),
    pair_method = "sampled_pairs_within_sample",
    pair_count_used = PAIR_SAMPLE_EXPERIMENT
  )
}

entropy_for_full_pixel_condition <- function(x_norm) {
  possible_pairs <- pair_count(nrow(x_norm))
  if (possible_pairs <= MAX_STORED_PAIRS) {
    return(list(
      entropy = spectral_entropy_from_norm(x_norm),
      pair_method = "all_pairs_full_pixel_condition",
      pair_count_used = possible_pairs
    ))
  }

  set.seed(seed_from_values("full_pixel_pair_sample", nrow(x_norm)))
  list(
    entropy = spectral_entropy_from_sampled_pairs(
      x_norm,
      n_pairs = PAIR_SAMPLE_EXPERIMENT,
      n_bins = N_BINS
    ),
    pair_method = "deterministic_sampled_pairs_full_pixel_condition",
    pair_count_used = PAIR_SAMPLE_EXPERIMENT
  )
}

run_condition <- function(x_norm, scale, quad_id, condition) {
  n_pixels <- nrow(x_norm)
  sample_size <- condition$sample_size
  set.seed(seed_from_values(scale, quad_id, condition$sample_rule))
  full_pixel_result <- NULL

  if (sample_size == n_pixels) {
    full_pixel_result <- entropy_for_full_pixel_condition(x_norm)
  }

  bind_rows(lapply(seq_len(N_BOOT_EXPERIMENT), function(iteration) {
    if (sample_size == n_pixels) {
      result <- full_pixel_result
    } else {
      sample_index <- sample.int(n_pixels, sample_size, replace = FALSE)
      result <- entropy_for_sample(x_norm[sample_index, , drop = FALSE])
    }

    data.frame(
      scale = scale,
      quad_id = quad_id,
      n_pixels = n_pixels,
      sample_rule = condition$sample_rule,
      rule_type = condition$rule_type,
      rule_value = condition$rule_value,
      sample_size = sample_size,
      sample_fraction = sample_size / n_pixels,
      bootstrap_iter = iteration,
      spectral_entropy = result$entropy,
      pair_method = result$pair_method,
      pair_count_used = result$pair_count_used,
      stringsAsFactors = FALSE
    )
  }))
}

sa_summaries <- bind_rows(lapply(names(target_quads_per_scale), read_sa_summary))

selected_quads <- bind_rows(lapply(names(target_quads_per_scale), function(scale_name) {
  select_additional_quads(
    scale_name,
    sa_summaries |> filter(scale == scale_name)
  )
})) |>
  arrange(factor(scale, levels = names(target_quads_per_scale)), selection_order) |>
  select(scale, quad_id, selection_note)

quad_metadata <- selected_quads |>
  rowwise() |>
  mutate(
    raster_file = find_raster_file(scale, quad_id),
    summary_csv = file.path("Quad_Values", paste0(scale, "_SA_entropy_smooth_masked_5nm_summary.csv"))
  ) |>
  ungroup() |>
  left_join(
    sa_summaries |>
      select(scale, quad_id, source_n_pixels = n_pixels, source_method = method, source_entropy = spectral_entropy),
    by = c("scale", "quad_id")
  ) |>
  arrange(factor(scale, levels = names(target_quads_per_scale)), quad_id)

design_table <- tidyr::crossing(
  quad_metadata,
  sample_rules
) |>
  rowwise() |>
  filter(rule_applies_to_scale(scale, applicable_scales)) |>
  ungroup() |>
  rowwise() |>
  mutate(
    sample_size = resolve_sample_size(source_n_pixels, rule_type, rule_value),
    sample_fraction = sample_size / source_n_pixels,
    sample_label = if_else(
      rule_type == "fixed",
      paste0(rule_label, " (", format_pct(sample_fraction), ")"),
      paste0(rule_label, " = ", comma(sample_size), " pixels (", format_pct(sample_fraction), ")")
    ),
    plot_sample_label = if_else(
      rule_type == "fixed",
      paste0(comma(sample_size), "\n", format_pct(sample_fraction)),
      paste0(gsub(" of pixels", "", rule_label), "\n", comma(sample_size), " px")
    )
  ) |>
  ungroup() |>
  arrange(scale, quad_id, sample_size)

write_csv(design_table, file.path(OUTPUT_DIR, "sa_entropy_sample_size_design.csv"))

bootstrap_results <- list()
result_index <- 1

for (quad_i in seq_len(nrow(quad_metadata))) {
  quad <- quad_metadata[quad_i, ]
  cat("Reading", quad$scale, quad$quad_id, "\n")
  spectra_norm <- read_quadrat_spectra_norm(quad$raster_file)
  actual_n_pixels <- nrow(spectra_norm)

  conditions <- design_table |>
    filter(scale == quad$scale, quad_id == quad$quad_id) |>
    mutate(
      n_pixels = actual_n_pixels,
      sample_size = pmin(sample_size, actual_n_pixels),
      sample_fraction = sample_size / actual_n_pixels
    )

  for (condition_i in seq_len(nrow(conditions))) {
    condition <- conditions[condition_i, ]
    cat(
      "  ",
      condition$sample_rule,
      "sample_size=",
      condition$sample_size,
      "fraction=",
      format_pct(condition$sample_fraction),
      "\n"
    )
    bootstrap_results[[result_index]] <- run_condition(
      spectra_norm,
      quad$scale,
      quad$quad_id,
      condition
    )
    result_index <- result_index + 1
  }

  rm(spectra_norm)
  gc(verbose = FALSE)
}

bootstrap_long <- bind_rows(bootstrap_results) |>
  left_join(
    design_table |>
      select(
        scale, quad_id, sample_rule, rule_label, rule_axis_label, sample_label, plot_sample_label,
        selection_note, source_method, source_entropy
      ),
    by = c("scale", "quad_id", "sample_rule")
  ) |>
  mutate(
    scale = factor(scale, levels = c("10m", "20m", "50m")),
    sample_rule = factor(sample_rule, levels = sample_rules$sample_rule),
    rule_axis_label = factor(rule_axis_label, levels = sample_rules$rule_axis_label),
    plot_sample_label = factor(plot_sample_label, levels = unique(design_table$plot_sample_label)),
    scale_quad = paste0(scale, " / ", quad_id, "\nn = ", comma(n_pixels), " px"),
    quad_axis_label = paste0(quad_id, "\nn = ", comma(n_pixels), " px")
  )

summary_table <- bootstrap_long |>
  group_by(
    scale, quad_id, n_pixels, sample_rule, rule_type, rule_value,
    rule_label, rule_axis_label, sample_label, sample_size, sample_fraction,
    source_method, source_entropy
  ) |>
  summarise(
    n_boot = sum(is.finite(spectral_entropy)),
    entropy_mean = mean(spectral_entropy, na.rm = TRUE),
    entropy_sd = sd(spectral_entropy, na.rm = TRUE),
    entropy_cv = entropy_sd / entropy_mean,
    entropy_se = entropy_sd / sqrt(n_boot),
    ci_t = qt(0.975, df = pmax(n_boot - 1, 1)),
    ci_half_width = ci_t * entropy_se,
    ci_low = entropy_mean - ci_half_width,
    ci_high = entropy_mean + ci_half_width,
    entropy_min = min(spectral_entropy, na.rm = TRUE),
    entropy_max = max(spectral_entropy, na.rm = TRUE),
    pair_method = paste(unique(pair_method), collapse = "; "),
    .groups = "drop"
  ) |>
  group_by(scale, quad_id) |>
  mutate(
    fixed_4000_reference = entropy_mean[sample_rule == "fixed_4000"][1],
    delta_from_fixed_4000 = entropy_mean - fixed_4000_reference
  ) |>
  ungroup() |>
  mutate(
    scale_quad = paste0(scale, " / ", quad_id, "\nn = ", comma(n_pixels), " px"),
    quad_axis_label = paste0(quad_id, "\nn = ", comma(n_pixels), " px")
  )

write_csv(bootstrap_long, file.path(OUTPUT_DIR, "sa_entropy_sample_size_boot_long.csv"))
write_csv(summary_table, file.path(OUTPUT_DIR, "sa_entropy_sample_size_summary.csv"))

theme_report <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", color = "#263238"),
    plot.subtitle = element_text(color = "#546E7A"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

point_color <- "#1B4D89"
line_color <- "#587FAE"
box_fill <- "#7BA7C7"
box_outline <- "#263238"

scale_plot_dims <- data.frame(
  scale = factor(c("10m", "20m", "50m"), levels = c("10m", "20m", "50m")),
  ncol = c(4L, 4L, 4L),
  width = c(19, 19, 16),
  height = c(22, 13, 8)
)

save_scale_figures <- function(scale_name) {
  summary_scale <- summary_table |> filter(scale == scale_name)
  long_scale <- bootstrap_long |> filter(scale == scale_name)
  dims <- scale_plot_dims |> filter(scale == scale_name)

  p_mean_scale <- ggplot(summary_scale, aes(x = sample_size, y = entropy_mean)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 55, alpha = 0.7, color = line_color) +
    geom_point(size = 2.1, color = point_color) +
    geom_line(aes(group = quad_id), linewidth = 0.65, alpha = 0.75, color = line_color) +
    facet_wrap(~ scale_quad, scales = "free_x", ncol = dims$ncol) +
    scale_x_continuous(labels = comma) +
    labs(
      title = paste(scale_name, "SA entropy mean by bootstrap sample size"),
      subtitle = "Points are means of 50 bootstrap iterations; bars are 95 percent CIs around each mean.",
      x = "Sampled retained pixels per iteration",
      y = "SA entropy mean"
    ) +
    theme_report
  ggsave(
    file.path(FIGURE_DIR, paste0("sa_entropy_mean_by_sample_size_", scale_name, ".png")),
    p_mean_scale,
    width = dims$width,
    height = dims$height,
    dpi = 320,
    bg = "white"
  )

  p_cv_scale <- ggplot(summary_scale, aes(x = sample_size, y = entropy_cv)) +
    geom_point(size = 2.1, color = point_color) +
    geom_line(aes(group = quad_id), linewidth = 0.65, alpha = 0.75, color = line_color) +
    facet_wrap(~ scale_quad, scales = "free_x", ncol = dims$ncol) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      title = paste(scale_name, "bootstrap CV by sample size"),
      subtitle = "Lower CV means less replicate-to-replicate variation for the same quadrat and rule.",
      x = "Sampled retained pixels per iteration",
      y = "Bootstrap CV"
    ) +
    theme_report
  ggsave(
    file.path(FIGURE_DIR, paste0("sa_entropy_cv_by_sample_size_", scale_name, ".png")),
    p_cv_scale,
    width = dims$width,
    height = dims$height,
    dpi = 320,
    bg = "white"
  )

  p_delta_scale <- ggplot(summary_scale, aes(x = sample_fraction, y = delta_from_fixed_4000)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(size = 2.1, color = point_color) +
    geom_line(aes(group = quad_id), linewidth = 0.65, alpha = 0.75, color = line_color) +
    facet_wrap(~ scale_quad, scales = "free_x", ncol = dims$ncol) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste(scale_name, "difference from the fixed 4,000-pixel rule"),
      subtitle = "Each point is the sample-size rule mean minus that quadrat's fixed-4,000 mean.",
      x = "Sampled fraction of retained pixels",
      y = "Mean difference in SA entropy"
    ) +
    theme_report
  ggsave(
    file.path(FIGURE_DIR, paste0("sa_entropy_delta_from_fixed_4000_", scale_name, ".png")),
    p_delta_scale,
    width = dims$width,
    height = dims$height,
    dpi = 320,
    bg = "white"
  )

  distribution_rules <- summary_scale |>
    distinct(sample_rule, rule_axis_label) |>
    mutate(sample_rule_order = match(as.character(sample_rule), sample_rules$sample_rule)) |>
    arrange(sample_rule_order)

  for (rule_i in seq_len(nrow(distribution_rules))) {
    rule_row <- distribution_rules[rule_i, ]
    rule_long <- long_scale |>
      filter(sample_rule == rule_row$sample_rule) |>
      mutate(quad_axis_label = factor(quad_axis_label, levels = unique(quad_axis_label)))

    p_rule_dist <- ggplot(rule_long, aes(x = quad_axis_label, y = spectral_entropy)) +
      geom_boxplot(outlier.alpha = 0.35, linewidth = 0.35, fill = box_fill, color = box_outline) +
      labs(
        title = paste(scale_name, "SA entropy replicate distribution:", as.character(rule_row$rule_axis_label)),
        subtitle = "Each box is one quadrat; x-axis labels include retained pixel counts.",
        x = "Quadrat",
        y = "Bootstrap SA entropy value"
      ) +
      theme_report +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none"
      )

    ggsave(
      file.path(
        FIGURE_DIR,
        "distributions_by_sample_size",
        paste0("sa_entropy_replicate_distribution_", scale_name, "_", rule_row$sample_rule, ".png")
      ),
      p_rule_dist,
      width = ifelse(scale_name == "10m", 20, ifelse(scale_name == "20m", 15, 11)),
      height = 6,
      dpi = 320,
      bg = "white"
    )
  }
}

invisible(lapply(levels(summary_table$scale), save_scale_figures))

p_mean_overview <- ggplot(summary_table, aes(x = sample_size, y = entropy_mean, group = quad_id)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 55, alpha = 0.25, color = line_color) +
  geom_line(alpha = 0.35, linewidth = 0.45, color = line_color) +
  geom_point(alpha = 0.75, size = 1.6, color = point_color) +
  facet_wrap(~ scale, scales = "free_x") +
  scale_x_continuous(labels = comma) +
  labs(
    title = "SA entropy mean by bootstrap sample size",
    subtitle = "Compact overview; scale-specific figures show each quadrat separately.",
    x = "Sampled retained pixels per iteration",
    y = "SA entropy mean"
  ) +
  theme_report
ggsave(file.path(FIGURE_DIR, "sa_entropy_mean_by_sample_size.png"), p_mean_overview, width = 12, height = 6, dpi = 320, bg = "white")

p_cv_overview <- ggplot(summary_table, aes(x = sample_size, y = entropy_cv, group = quad_id)) +
  geom_line(alpha = 0.35, linewidth = 0.45, color = line_color) +
  geom_point(alpha = 0.75, size = 1.6, color = point_color) +
  facet_wrap(~ scale, scales = "free_x") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Bootstrap CV by sample size",
    subtitle = "Compact overview; scale-specific figures show each quadrat separately.",
    x = "Sampled retained pixels per iteration",
    y = "Bootstrap CV"
  ) +
  theme_report
ggsave(file.path(FIGURE_DIR, "sa_entropy_cv_by_sample_size.png"), p_cv_overview, width = 12, height = 6, dpi = 320, bg = "white")

p_delta_overview <- ggplot(summary_table, aes(x = sample_fraction, y = delta_from_fixed_4000, group = quad_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_line(alpha = 0.35, linewidth = 0.45, color = line_color) +
  geom_point(alpha = 0.75, size = 1.6, color = point_color) +
  facet_wrap(~ scale, scales = "free_x") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Difference from the fixed 4,000-pixel rule",
    subtitle = "Compact overview; scale-specific figures show each quadrat separately.",
    x = "Sampled fraction of retained pixels",
    y = "Mean difference in SA entropy"
  ) +
  theme_report
ggsave(file.path(FIGURE_DIR, "sa_entropy_delta_from_fixed_4000.png"), p_delta_overview, width = 12, height = 6, dpi = 320, bg = "white")

p_reps_overview <- ggplot(bootstrap_long, aes(x = rule_axis_label, y = spectral_entropy)) +
  geom_boxplot(outlier.alpha = 0.12, linewidth = 0.25, fill = box_fill, color = box_outline) +
  facet_wrap(~ scale, scales = "free_y") +
  scale_x_discrete(drop = TRUE) +
  labs(
    title = "Bootstrap replicate distributions by sample-size rule",
    subtitle = "Compact overview across all selected quadrats; scale-specific figures separate quadrats.",
    x = "Sample-size rule",
    y = "Bootstrap SA entropy value"
  ) +
  theme_report +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8), legend.position = "none")
ggsave(file.path(FIGURE_DIR, "sa_entropy_replicate_distributions.png"), p_reps_overview, width = 12, height = 6, dpi = 320, bg = "white")

design_report <- design_table |>
  transmute(
    Scale = scale,
    Quadrat = quad_id,
    `Retained pixels` = comma(source_n_pixels),
    `Sample rule` = sample_label,
    `Sample pixels` = comma(sample_size),
    `Sample percent` = format_pct(sample_fraction),
    `Source method` = source_method
  )

selection_report <- quad_metadata |>
  transmute(
    Scale = scale,
    Quadrat = quad_id,
    `Retained pixels` = comma(source_n_pixels),
    `Source SA mean` = format_number(source_entropy),
    `Selection note` = selection_note
  )

summary_report <- summary_table |>
  transmute(
    Scale = as.character(scale),
    Quadrat = quad_id,
    `Sample rule` = sample_label,
    `Mean entropy` = format_number(entropy_mean),
    `Boot SD` = format_number(entropy_sd),
    `Boot CV` = format_pct(entropy_cv),
    `95% CI low` = format_number(ci_low),
    `95% CI high` = format_number(ci_high),
    `Delta from fixed 4,000` = format_number(delta_from_fixed_4000),
    `Pair method` = pair_method
  )

scale_figure_lines <- unlist(lapply(c("10m", "20m", "50m"), function(scale_name) {
  c(
    paste0("### ", scale_name, " SA Entropy Mean By Sample Size"),
    "",
    paste0(
      "![", scale_name, " SA entropy mean by sample size]",
      "(../figures/sample_size_effects/sa_entropy/sa_entropy_mean_by_sample_size_",
      scale_name,
      ".png)"
    ),
    "",
    paste0("### ", scale_name, " Bootstrap CV By Sample Size"),
    "",
    paste0(
      "![", scale_name, " bootstrap CV by sample size]",
      "(../figures/sample_size_effects/sa_entropy/sa_entropy_cv_by_sample_size_",
      scale_name,
      ".png)"
    ),
    "",
    paste0("### ", scale_name, " Difference From Fixed 4,000 Pixels"),
    "",
    paste0(
      "![", scale_name, " difference from fixed 4,000 pixels]",
      "(../figures/sample_size_effects/sa_entropy/sa_entropy_delta_from_fixed_4000_",
      scale_name,
      ".png)"
    ),
    "",
    ""
  )
}))

split_distribution_lines <- unlist(lapply(c("10m", "20m", "50m"), function(scale_name) {
  rules_for_scale <- design_table |>
    filter(scale == scale_name) |>
    distinct(sample_rule, rule_axis_label) |>
    mutate(sample_rule_order = match(as.character(sample_rule), sample_rules$sample_rule)) |>
    arrange(sample_rule_order)

  c(
    paste0("### ", scale_name, " Distribution Charts Split By Sample Size"),
    "",
    unlist(lapply(seq_len(nrow(rules_for_scale)), function(rule_i) {
      rule_row <- rules_for_scale[rule_i, ]
      c(
        paste0(
          "- [", scale_name, " ", as.character(rule_row$rule_axis_label),
          "](../figures/sample_size_effects/sa_entropy/distributions_by_sample_size/sa_entropy_replicate_distribution_",
          scale_name,
          "_",
          rule_row$sample_rule,
          ".png)"
        )
      )
    })),
    ""
  )
}))

report_lines <- c(
  "# SA Entropy Sample-Size Effects",
  "",
  "Date: 2026-07-04",
  "",
  "## Purpose",
  "",
  "Compare how spectral angle entropy responds to different bootstrap pixel sample sizes across 10 m, 20 m, and 50 m quadrats.",
  "",
  "## Design",
  "",
  paste0(
    "- Bootstrap iterations per quadrat x sample rule: ",
    N_BOOT_EXPERIMENT
  ),
  paste0(
    "- Quadrat sample: ",
    target_quads_per_scale[["10m"]],
    " 10 m quadrats, ",
    target_quads_per_scale[["20m"]],
    " 20 m quadrats, and ",
    target_quads_per_scale[["50m"]],
    " 50 m quadrats. The original six quadrats were retained and additional quadrats were sampled reproducibly from `method = \"bootstrap_mean\"` rows."
  ),
  paste0(
    "- Fixed-count rules: 1,250 and 4,000 pixels at all scales; 2,000 and 3,000 pixels added for 10 m and 20 m; 6,000 and 8,000 pixels added for 50 m. Fixed-count rules are capped by available retained pixels."
  ),
  paste0(
    "- Percent rules: 1%, 2%, and 3% of retained pixels, capped at ",
    comma(MAX_SAMPLE_PIXELS),
    " pixels."
  ),
  paste0(
    "- Large sampled-pixel iterations use ",
    comma(PAIR_SAMPLE_EXPERIMENT),
    " sampled spectral-angle pairs when all within-sample pairs exceed `MAX_BOOT_EXACT_PAIRS`."
  ),
  "- When a rule resolves to 100% of retained pixels, the entropy is calculated once from the full retained-pixel set and repeated across the 50 output rows. This removes artificial bootstrap variation for full-pixel conditions.",
  "",
  "## 100% Pixel Check",
  "",
  "The earlier pilot output showed non-zero variation for the `800_a` fixed-4,000 rule even though that rule used 100% of the retained pixels. The source was not changing pixels; it was the sampled-pair fallback inside each replicate. This version treats full-pixel conditions as deterministic full-pixel conditions, so the 100% row has zero bootstrap SD and zero CI width.",
  "",
  "The 10 m quadrat `800_a` replaced `1304_a` because it is the closest bootstrap-mean 10 m quadrat at or below 4,000 retained pixels.",
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
  "### Compact Overview: SA Entropy Mean By Sample Size",
  "",
  "![SA entropy mean by sample size](../figures/sample_size_effects/sa_entropy/sa_entropy_mean_by_sample_size.png)",
  "",
  "### Compact Overview: Bootstrap CV By Sample Size",
  "",
  "![Bootstrap CV by sample size](../figures/sample_size_effects/sa_entropy/sa_entropy_cv_by_sample_size.png)",
  "",
  "### Compact Overview: Difference From Fixed 4,000 Pixels",
  "",
  "![Difference from fixed 4,000 pixels](../figures/sample_size_effects/sa_entropy/sa_entropy_delta_from_fixed_4000.png)",
  "",
  "### Compact Overview: Bootstrap Replicate Distributions",
  "",
  "![Bootstrap replicate distributions](../figures/sample_size_effects/sa_entropy/sa_entropy_replicate_distributions.png)",
  "",
  scale_figure_lines,
  "",
  "## Distribution Charts Split By Sample Size",
  "",
  "The compact distribution overview above is retained. The following long charts split replicate distributions by sample-size rule and show quadrats as boxes with retained pixel counts in the x-axis labels.",
  "",
  split_distribution_lines,
  "",
  "## Output Tables",
  "",
  "- `reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_design.csv`",
  "- `reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_boot_long.csv`",
  "- `reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_summary.csv`"
)

writeLines(report_lines, REPORT_PATH)

cat("SA entropy sample-size effects complete\n")
cat("Report:", REPORT_PATH, "\n")
cat("Tables:", OUTPUT_DIR, "\n")
cat("Figures:", FIGURE_DIR, "\n")
