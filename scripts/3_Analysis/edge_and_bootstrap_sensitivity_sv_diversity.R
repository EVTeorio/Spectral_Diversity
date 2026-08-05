PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"

INPUT_DATA <- file.path(
  PROJECT_DIR,
  "reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv"
)
BOOTSTRAP_DATA <- file.path(
  PROJECT_DIR,
  "reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv"
)

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/edge_bootstrap_sensitivity")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/edge_bootstrap_sensitivity")
REPORT_PATH <- file.path(PROJECT_DIR, "reports/analysis/20260725_edge_bootstrap_sensitivity_sv_diversity.md")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REPORT_PATH), recursive = TRUE, showWarnings = FALSE)

DIVERSITY_MEASURES <- c("phy_faith", "phy_rao", "phy_afaith", "sp_shannon")
SV_MEASURES <- c(
  "spec_sa",
  "spec_spca_mean",
  "spec_spca_med",
  "spec_spca_sd",
  "spec_spca_rao",
  "spec_spca_alpha",
  "spec_spca_convex",
  "spec_spca_hull3d_v",
  "spec_spca_hull3d_a"
)
SCALES <- c("10m", "20m")
N_RESAMPLE <- 2000
ENV_MEASURES <- c("env_elev", "env_tri5", "env_tri11", "env_tri21")

DISPLAY_NAMES <- c(
  spec_sa = "SA entropy",
  spec_spca_mean = "Std PCA mean distance",
  spec_spca_med = "Std PCA median distance",
  spec_spca_sd = "Std PCA distance SD",
  spec_spca_rao = "Std PCA Rao's Q",
  spec_spca_alpha = "Std PCA alpha-hull area",
  spec_spca_convex = "Std PCA convex-hull area",
  spec_spca_hull3d_v = "Std PCA 3D hull volume",
  spec_spca_hull3d_a = "Std PCA 3D hull area",
  phy_faith = "Faith's PD",
  phy_rao = "Phylogenetic Rao's Q",
  phy_afaith = "Abundance-weighted Faith's PD",
  sp_shannon = "Shannon diversity",
  env_elev = "Elevation",
  env_tri5 = "TRI 5x5",
  env_tri11 = "TRI 11x11",
  env_tri21 = "TRI 21x21"
)

display_name <- function(x) unname(DISPLAY_NAMES[x])

fmt_num <- function(x, digits = 3) {
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(as.numeric(x[ok]), digits = digits, format = "f")
  out
}

fmt_p <- function(x) {
  out <- rep("NA", length(x))
  ok <- !is.na(x)
  out[ok] <- ifelse(
    as.numeric(x[ok]) < 0.001,
    formatC(as.numeric(x[ok]), digits = 2, format = "e"),
    formatC(as.numeric(x[ok]), digits = 3, format = "f")
  )
  out
}

markdown_table <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  header <- paste0("| ", paste(names(data), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |")
  rows <- apply(data, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(header, separator, rows)
}

relative_path <- function(paths) {
  normalized_project <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = FALSE)
  normalized_paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", normalized_project), "/"), "", normalized_paths)
}

safe_cor_lm <- function(data, scale_name, sv_measure, diversity_measure, analysis_set) {
  pair_data <- data[
    data$scale == scale_name &
      !is.na(data[[sv_measure]]) &
      !is.na(data[[diversity_measure]]),
  ]

  if (analysis_set == "non_edge") {
    pair_data <- pair_data[!pair_data$edge_flag, ]
  } else if (analysis_set == "edge_only") {
    pair_data <- pair_data[pair_data$edge_flag, ]
  }

  result <- data.frame(
    scale = scale_name,
    analysis_set = analysis_set,
    sv_measure = sv_measure,
    sv_label = display_name(sv_measure),
    diversity_measure = diversity_measure,
    diversity_label = display_name(diversity_measure),
    n = nrow(pair_data),
    edge_n = sum(pair_data$edge_flag),
    pearson_r = NA_real_,
    r_squared = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    spearman_r = NA_real_,
    spearman_p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(pair_data) < 3 ||
      length(unique(pair_data[[sv_measure]])) < 2 ||
      length(unique(pair_data[[diversity_measure]])) < 2) {
    return(result)
  }

  model_data <- data.frame(
    response = pair_data[[sv_measure]],
    predictor = pair_data[[diversity_measure]]
  )
  fit <- lm(response ~ predictor, data = model_data)
  fit_summary <- summary(fit)
  f_values <- fit_summary$fstatistic
  pearson <- cor.test(model_data$predictor, model_data$response, method = "pearson")
  spearman <- suppressWarnings(cor.test(
    model_data$predictor,
    model_data$response,
    method = "spearman",
    exact = FALSE
  ))

  result$pearson_r <- unname(pearson$estimate)
  result$r_squared <- result$pearson_r^2
  result$slope <- unname(coef(fit)[["predictor"]])
  result$intercept <- unname(coef(fit)[["(Intercept)"]])
  result$f_statistic <- unname(f_values[["value"]])
  result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
  result$spearman_r <- unname(spearman$estimate)
  result$spearman_p_value <- spearman$p.value
  result
}

run_edge_correlations <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (sv_measure in SV_MEASURES) {
      for (diversity_measure in DIVERSITY_MEASURES) {
        for (analysis_set in c("all_quads", "non_edge", "edge_only")) {
          rows[[index]] <- safe_cor_lm(data, scale_name, sv_measure, diversity_measure, analysis_set)
          index <- index + 1
        }
      }
    }
  }
  results <- do.call(rbind, rows)
  all_results <- results[results$analysis_set == "all_quads", ]
  non_edge_results <- results[results$analysis_set == "non_edge", ]
  edge_only_results <- results[results$analysis_set == "edge_only", ]
  comparison <- merge(
    all_results,
    non_edge_results,
    by = c("scale", "sv_measure", "sv_label", "diversity_measure", "diversity_label"),
    suffixes = c("_all", "_non_edge"),
    sort = FALSE
  )
  comparison <- merge(
    comparison,
    edge_only_results,
    by = c("scale", "sv_measure", "sv_label", "diversity_measure", "diversity_label"),
    sort = FALSE
  )
  names(comparison) <- sub("^analysis_set$", "analysis_set_edge_only", names(comparison))
  names(comparison) <- sub("^n$", "n_edge_only", names(comparison))
  names(comparison) <- sub("^edge_n$", "edge_n_edge_only", names(comparison))
  names(comparison) <- sub("^pearson_r$", "pearson_r_edge_only", names(comparison))
  names(comparison) <- sub("^r_squared$", "r_squared_edge_only", names(comparison))
  names(comparison) <- sub("^slope$", "slope_edge_only", names(comparison))
  names(comparison) <- sub("^intercept$", "intercept_edge_only", names(comparison))
  names(comparison) <- sub("^f_statistic$", "f_statistic_edge_only", names(comparison))
  names(comparison) <- sub("^f_p_value$", "f_p_value_edge_only", names(comparison))
  names(comparison) <- sub("^spearman_r$", "spearman_r_edge_only", names(comparison))
  names(comparison) <- sub("^spearman_p_value$", "spearman_p_value_edge_only", names(comparison))
  comparison$delta_r_non_edge_minus_all <- comparison$pearson_r_non_edge - comparison$pearson_r_all
  comparison$delta_r2_non_edge_minus_all <- comparison$r_squared_non_edge - comparison$r_squared_all
  comparison$delta_r_edge_only_minus_non_edge <- comparison$pearson_r_edge_only - comparison$pearson_r_non_edge
  comparison$delta_r_edge_only_minus_all <- comparison$pearson_r_edge_only - comparison$pearson_r_all
  list(results = results, comparison = comparison)
}

save_edge_scatter_plot <- function(data, scale_name, sv_measure) {
  plot_data <- data[data$scale == scale_name & !is.na(data[[sv_measure]]), ]
  file_path <- file.path(FIGURE_DIR, paste0(scale_name, "_", sv_measure, "_edge_highlight_scatter.png"))

  grDevices::png(file_path, width = 2600, height = 2100, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  par(mfrow = c(2, 2), mar = c(4.4, 4.7, 3.2, 1), oma = c(0, 0, 2.1, 0))
  for (diversity_measure in DIVERSITY_MEASURES) {
    pair_data <- plot_data[!is.na(plot_data[[diversity_measure]]), ]
    x <- pair_data[[diversity_measure]]
    y <- pair_data[[sv_measure]]
    edge <- pair_data$edge_flag
    plot(
      x[!edge],
      y[!edge],
      pch = 16,
      col = grDevices::adjustcolor("#2F6F73", alpha.f = 0.45),
      xlab = display_name(diversity_measure),
      ylab = display_name(sv_measure),
      main = display_name(diversity_measure),
      xlim = range(x, na.rm = TRUE),
      ylim = range(y, na.rm = TRUE)
    )
    points(
      x[edge],
      y[edge],
      pch = 21,
      bg = grDevices::adjustcolor("#D55E00", alpha.f = 0.65),
      col = "#8A3B00",
      cex = 0.85,
      lwd = 0.7
    )
    if (nrow(pair_data) >= 3) {
      abline(lm(y ~ x), col = "#222222", lwd = 2)
    }
    non_edge_data <- pair_data[!pair_data$edge_flag, ]
    if (nrow(non_edge_data) >= 3) {
      abline(lm(non_edge_data[[sv_measure]] ~ non_edge_data[[diversity_measure]]), col = "#2F6F73", lwd = 2, lty = 2)
    }
    edge_data <- pair_data[pair_data$edge_flag, ]
    if (nrow(edge_data) >= 3) {
      abline(lm(edge_data[[sv_measure]] ~ edge_data[[diversity_measure]]), col = "#D55E00", lwd = 2, lty = 3)
    }
    legend(
      "topleft",
      legend = c("non-edge", "edge", "all-quads fit", "non-edge fit", "edge-only fit"),
      pch = c(16, 21, NA, NA, NA),
      pt.bg = c(NA, "#D55E00", NA, NA, NA),
      col = c("#2F6F73", "#8A3B00", "#222222", "#2F6F73", "#D55E00"),
      lty = c(NA, NA, 1, 2, 3),
      lwd = c(NA, NA, 2, 2, 2),
      bty = "n",
      cex = 0.68
    )
  }
  mtext(paste(scale_name, display_name(sv_measure), "edge-quadrat sensitivity"), outer = TRUE, cex = 1.15, font = 2)
  file_path
}

save_all_edge_scatter_plots <- function(data) {
  files <- character()
  for (scale_name in SCALES) {
    for (sv_measure in SV_MEASURES) {
      files <- c(files, save_edge_scatter_plot(data, scale_name, sv_measure))
    }
  }
  files
}

read_bootstrap_data <- function() {
  boot <- read.csv(BOOTSTRAP_DATA, stringsAsFactors = FALSE, check.names = FALSE)
  boot$scale <- gsub(" ", "", boot$scale)
  boot$quad_id <- as.character(boot$quad_id)
  keep_cols <- c(
    "quad_id", "scale", "boot_sd", "boot_cv", "boot_ci_width",
    "boot_ci_half_width_pct", "n_pixels", "n_boot", "valid_entropy"
  )
  boot[, keep_cols]
}

join_bootstrap_fields <- function(data, boot) {
  data$row_order <- seq_len(nrow(data))
  merged <- merge(data, boot, by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
  merged <- merged[order(merged$row_order), ]
  merged$row_order <- NULL

  merged$boot_cv_class_5_10 <- NA_character_
  has_cv <- !is.na(merged$boot_cv)
  merged$boot_cv_class_5_10[has_cv & merged$boot_cv <= 0.05] <- "cv_le_5pct"
  merged$boot_cv_class_5_10[has_cv & merged$boot_cv > 0.05 & merged$boot_cv <= 0.10] <- "cv_5_to_10pct"
  merged$boot_cv_class_5_10[has_cv & merged$boot_cv > 0.10] <- "cv_gt_10pct"

  merged$boot_sd_tertile <- NA_character_
  for (scale_name in SCALES) {
    idx <- which(merged$scale == scale_name & !is.na(merged$boot_sd))
    if (length(idx) >= 3) {
      breaks <- unique(as.numeric(quantile(merged$boot_sd[idx], probs = c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE)))
      if (length(breaks) == 4) {
        merged$boot_sd_tertile[idx] <- as.character(cut(
          merged$boot_sd[idx],
          breaks = breaks,
          include.lowest = TRUE,
          labels = c("low_boot_sd", "mid_boot_sd", "high_boot_sd")
        ))
      }
    }
  }
  merged
}

run_bootstrap_sensitivity <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (diversity_measure in DIVERSITY_MEASURES) {
      base_filter <- data$scale == scale_name &
        !data$edge_flag &
        !is.na(data$spec_sa) &
        !is.na(data[[diversity_measure]]) &
        !is.na(data$boot_cv)

      filters <- list(
        non_edge_all_bootstrap = base_filter,
        boot_cv_le_5pct = base_filter & data$boot_cv <= 0.05,
        boot_cv_le_10pct = base_filter & data$boot_cv <= 0.10,
        boot_cv_gt_5pct = base_filter & data$boot_cv > 0.05,
        boot_cv_gt_10pct = base_filter & data$boot_cv > 0.10
      )

      for (filter_name in names(filters)) {
        subset_data <- data[filters[[filter_name]], ]
        result <- data.frame(
          scale = scale_name,
          bootstrap_subset = filter_name,
          diversity_measure = diversity_measure,
          diversity_label = display_name(diversity_measure),
          n = nrow(subset_data),
          mean_boot_sd = mean(subset_data$boot_sd, na.rm = TRUE),
          mean_boot_cv = mean(subset_data$boot_cv, na.rm = TRUE),
          pearson_r = NA_real_,
          r_squared = NA_real_,
          f_p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        if (nrow(subset_data) >= 3 &&
            length(unique(subset_data$spec_sa)) >= 2 &&
            length(unique(subset_data[[diversity_measure]])) >= 2) {
          fit <- lm(subset_data$spec_sa ~ subset_data[[diversity_measure]])
          f_values <- summary(fit)$fstatistic
          pearson <- cor.test(subset_data[[diversity_measure]], subset_data$spec_sa)
          result$pearson_r <- unname(pearson$estimate)
          result$r_squared <- result$pearson_r^2
          result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
        }
        rows[[index]] <- result
        index <- index + 1
      }

      for (tertile in c("low_boot_sd", "mid_boot_sd", "high_boot_sd")) {
        subset_data <- data[base_filter & data$boot_sd_tertile == tertile, ]
        result <- data.frame(
          scale = scale_name,
          bootstrap_subset = tertile,
          diversity_measure = diversity_measure,
          diversity_label = display_name(diversity_measure),
          n = nrow(subset_data),
          mean_boot_sd = mean(subset_data$boot_sd, na.rm = TRUE),
          mean_boot_cv = mean(subset_data$boot_cv, na.rm = TRUE),
          pearson_r = NA_real_,
          r_squared = NA_real_,
          f_p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        if (nrow(subset_data) >= 3 &&
            length(unique(subset_data$spec_sa)) >= 2 &&
            length(unique(subset_data[[diversity_measure]])) >= 2) {
          fit <- lm(subset_data$spec_sa ~ subset_data[[diversity_measure]])
          f_values <- summary(fit)$fstatistic
          pearson <- cor.test(subset_data[[diversity_measure]], subset_data$spec_sa)
          result$pearson_r <- unname(pearson$estimate)
          result$r_squared <- result$pearson_r^2
          result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
        }
        rows[[index]] <- result
        index <- index + 1
      }
    }
  }
  do.call(rbind, rows)
}

save_bootstrap_scatter_plot <- function(data, scale_name) {
  plot_data <- data[
    data$scale == scale_name &
      !data$edge_flag &
      !is.na(data$spec_sa) &
      !is.na(data$boot_cv),
  ]
  file_path <- file.path(FIGURE_DIR, paste0(scale_name, "_sa_entropy_bootstrap_cv_sensitivity.png"))

  grDevices::png(file_path, width = 2600, height = 2100, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  point_cols <- c(
    cv_le_5pct = grDevices::adjustcolor("#2F6F73", alpha.f = 0.55),
    cv_5_to_10pct = grDevices::adjustcolor("#E69F00", alpha.f = 0.70),
    cv_gt_10pct = grDevices::adjustcolor("#C44E52", alpha.f = 0.75)
  )

  par(mfrow = c(2, 2), mar = c(4.4, 4.7, 3.2, 1), oma = c(0, 0, 2.1, 0))
  for (diversity_measure in DIVERSITY_MEASURES) {
    pair_data <- plot_data[!is.na(plot_data[[diversity_measure]]), ]
    x <- pair_data[[diversity_measure]]
    y <- pair_data$spec_sa
    classes <- pair_data$boot_cv_class_5_10
    plot(
      x,
      y,
      pch = 16,
      col = point_cols[classes],
      xlab = display_name(diversity_measure),
      ylab = display_name("spec_sa"),
      main = display_name(diversity_measure)
    )
    if (nrow(pair_data) >= 3) {
      abline(lm(y ~ x), col = "#222222", lwd = 2)
    }
    stable_data <- pair_data[pair_data$boot_cv <= 0.05, ]
    if (nrow(stable_data) >= 3) {
      abline(lm(stable_data$spec_sa ~ stable_data[[diversity_measure]]), col = "#2F6F73", lwd = 2, lty = 2)
    }
    legend(
      "topleft",
      legend = c("CV <= 5%", "5% < CV <= 10%", "CV > 10%", "all fit", "CV <= 5% fit"),
      col = c("#2F6F73", "#E69F00", "#C44E52", "#222222", "#2F6F73"),
      pch = c(16, 16, 16, NA, NA),
      lty = c(NA, NA, NA, 1, 2),
      lwd = c(NA, NA, NA, 2, 2),
      bty = "n",
      cex = 0.72
    )
  }
  mtext(paste(scale_name, "SA entropy bootstrap-CV sensitivity"), outer = TRUE, cex = 1.15, font = 2)
  file_path
}

save_bootstrap_plots <- function(data) {
  c(save_bootstrap_scatter_plot(data, "10m"), save_bootstrap_scatter_plot(data, "20m"))
}

summarize_edge_nonedge_ranges <- function(data) {
  variables <- c(DIVERSITY_MEASURES, SV_MEASURES)
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (variable in variables) {
      for (group_name in c("non_edge", "edge_only")) {
        group_data <- data[data$scale == scale_name, ]
        group_data <- if (group_name == "edge_only") {
          group_data[group_data$edge_flag, ]
        } else {
          group_data[!group_data$edge_flag, ]
        }
        values <- group_data[[variable]]
        values <- values[!is.na(values)]
        rows[[index]] <- data.frame(
          scale = scale_name,
          variable = variable,
          variable_label = display_name(variable),
          edge_group = group_name,
          n = length(values),
          mean = mean(values, na.rm = TRUE),
          sd = sd(values, na.rm = TRUE),
          min = min(values, na.rm = TRUE),
          q25 = as.numeric(quantile(values, 0.25, na.rm = TRUE)),
          median = median(values, na.rm = TRUE),
          q75 = as.numeric(quantile(values, 0.75, na.rm = TRUE)),
          max = max(values, na.rm = TRUE),
          range = diff(range(values, na.rm = TRUE)),
          iqr = IQR(values, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        index <- index + 1
      }
    }
  }

  summary <- do.call(rbind, rows)
  non_edge <- summary[summary$edge_group == "non_edge", ]
  edge_only <- summary[summary$edge_group == "edge_only", ]
  comparison <- merge(
    edge_only,
    non_edge,
    by = c("scale", "variable", "variable_label"),
    suffixes = c("_edge_only", "_non_edge"),
    sort = FALSE
  )
  comparison$range_ratio_edge_to_non_edge <- comparison$range_edge_only / comparison$range_non_edge
  comparison$sd_ratio_edge_to_non_edge <- comparison$sd_edge_only / comparison$sd_non_edge
  comparison$iqr_ratio_edge_to_non_edge <- comparison$iqr_edge_only / comparison$iqr_non_edge
  list(summary = summary, comparison = comparison)
}

run_equal_n_resampling <- function(data, edge_comparison) {
  set.seed(20260725)
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (sv_measure in SV_MEASURES) {
      for (diversity_measure in DIVERSITY_MEASURES) {
        pair_data <- data[
          data$scale == scale_name &
            !is.na(data[[sv_measure]]) &
            !is.na(data[[diversity_measure]]),
        ]
        edge_data <- pair_data[pair_data$edge_flag, ]
        non_edge_data <- pair_data[!pair_data$edge_flag, ]
        n_edge <- nrow(edge_data)
        edge_r <- edge_comparison$pearson_r_edge_only[
          edge_comparison$scale == scale_name &
            edge_comparison$sv_measure == sv_measure &
            edge_comparison$diversity_measure == diversity_measure
        ]

        result <- data.frame(
          scale = scale_name,
          sv_measure = sv_measure,
          sv_label = display_name(sv_measure),
          diversity_measure = diversity_measure,
          diversity_label = display_name(diversity_measure),
          n_edge = n_edge,
          n_non_edge_available = nrow(non_edge_data),
          edge_r = edge_r,
          resample_reps = N_RESAMPLE,
          resampled_non_edge_mean_r = NA_real_,
          resampled_non_edge_sd_r = NA_real_,
          resampled_non_edge_q025_r = NA_real_,
          resampled_non_edge_median_r = NA_real_,
          resampled_non_edge_q975_r = NA_real_,
          edge_r_percentile_in_resampled_non_edge = NA_real_,
          abs_edge_r_tail_probability = NA_real_,
          same_direction_tail_probability = NA_real_,
          stringsAsFactors = FALSE
        )

        if (n_edge >= 3 &&
            nrow(non_edge_data) >= n_edge &&
            !is.na(edge_r) &&
            length(unique(non_edge_data[[sv_measure]])) >= 2 &&
            length(unique(non_edge_data[[diversity_measure]])) >= 2) {
          sampled_r <- replicate(N_RESAMPLE, {
            selected <- sample(seq_len(nrow(non_edge_data)), n_edge, replace = FALSE)
            sample_data <- non_edge_data[selected, ]
            if (length(unique(sample_data[[sv_measure]])) < 2 ||
                length(unique(sample_data[[diversity_measure]])) < 2) {
              return(NA_real_)
            }
            cor(sample_data[[diversity_measure]], sample_data[[sv_measure]], use = "complete.obs")
          })
          sampled_r <- sampled_r[!is.na(sampled_r)]
          result$resample_reps <- length(sampled_r)
          result$resampled_non_edge_mean_r <- mean(sampled_r)
          result$resampled_non_edge_sd_r <- sd(sampled_r)
          result$resampled_non_edge_q025_r <- as.numeric(quantile(sampled_r, 0.025))
          result$resampled_non_edge_median_r <- median(sampled_r)
          result$resampled_non_edge_q975_r <- as.numeric(quantile(sampled_r, 0.975))
          result$edge_r_percentile_in_resampled_non_edge <- mean(sampled_r <= edge_r)
          result$abs_edge_r_tail_probability <- mean(abs(sampled_r) >= abs(edge_r))
          result$same_direction_tail_probability <- if (edge_r >= 0) {
            mean(sampled_r >= edge_r)
          } else {
            mean(sampled_r <= edge_r)
          }
        }
        rows[[index]] <- result
        index <- index + 1
      }
    }
  }
  do.call(rbind, rows)
}

save_resampling_summary_plot <- function(resampling_results) {
  primary <- resampling_results[
    resampling_results$sv_measure %in% c("spec_sa", "spec_spca_mean"),
  ]
  file_path <- file.path(FIGURE_DIR, "equal_n_non_edge_resampling_primary_measures.png")

  short_diversity_names <- c(
    phy_faith = "Faith PD",
    phy_rao = "Phylo Rao Q",
    phy_afaith = "A-weighted Faith PD",
    sp_shannon = "Shannon"
  )

  grDevices::png(file_path, width = 3400, height = 2300, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  par(mfrow = c(2, 2), mar = c(5.2, 7.2, 3.2, 1.2), oma = c(0, 0, 2.1, 0))
  for (scale_name in SCALES) {
    for (sv_measure in c("spec_sa", "spec_spca_mean")) {
      plot_data <- primary[primary$scale == scale_name & primary$sv_measure == sv_measure, ]
      y_positions <- seq_len(nrow(plot_data))
      xlim <- range(c(
        plot_data$resampled_non_edge_q025_r,
        plot_data$resampled_non_edge_q975_r,
        plot_data$edge_r
      ), na.rm = TRUE)
      plot(
        NA,
        xlim = xlim,
        ylim = c(0.5, length(y_positions) + 0.5),
        yaxt = "n",
        xlab = "Pearson r",
        ylab = "",
        main = paste(scale_name, display_name(sv_measure)),
        bty = "n"
      )
      axis(2, at = y_positions, labels = unname(short_diversity_names[plot_data$diversity_measure]), las = 2, cex.axis = 0.76)
      abline(v = 0, col = "gray75", lty = 3)
      segments(
        plot_data$resampled_non_edge_q025_r,
        y_positions,
        plot_data$resampled_non_edge_q975_r,
        y_positions,
        col = "#2F6F73",
        lwd = 4
      )
      points(plot_data$resampled_non_edge_mean_r, y_positions, pch = 16, col = "#2F6F73", cex = 0.9)
      points(plot_data$edge_r, y_positions, pch = 21, bg = "#D55E00", col = "#8A3B00", cex = 1.1)
      legend(
        "topright",
        legend = c("non-edge equal-n 95% interval", "non-edge equal-n mean", "edge-only r"),
        col = c("#2F6F73", "#2F6F73", "#8A3B00"),
        pt.bg = c(NA, NA, "#D55E00"),
        lwd = c(4, NA, NA),
        pch = c(NA, 16, 21),
        bty = "n",
        cex = 0.72
      )
    }
  }
  mtext("Equal-n resampling: edge-only r versus non-edge subsets", outer = TRUE, cex = 1.15, font = 2)
  file_path
}

summarize_environment_edge_nonedge <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (env_measure in ENV_MEASURES) {
      for (group_name in c("non_edge", "edge_only")) {
        group_data <- data[data$scale == scale_name, ]
        group_data <- if (group_name == "edge_only") {
          group_data[group_data$edge_flag, ]
        } else {
          group_data[!group_data$edge_flag, ]
        }
        values <- group_data[[env_measure]]
        values <- values[!is.na(values)]
        rows[[index]] <- data.frame(
          scale = scale_name,
          env_measure = env_measure,
          env_label = display_name(env_measure),
          edge_group = group_name,
          n = length(values),
          mean = mean(values, na.rm = TRUE),
          sd = sd(values, na.rm = TRUE),
          median = median(values, na.rm = TRUE),
          min = min(values, na.rm = TRUE),
          max = max(values, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        index <- index + 1
      }
    }
  }

  summary <- do.call(rbind, rows)
  edge_only <- summary[summary$edge_group == "edge_only", ]
  non_edge <- summary[summary$edge_group == "non_edge", ]
  comparison <- merge(
    edge_only,
    non_edge,
    by = c("scale", "env_measure", "env_label"),
    suffixes = c("_edge_only", "_non_edge"),
    sort = FALSE
  )
  comparison$mean_difference_edge_minus_non_edge <- comparison$mean_edge_only - comparison$mean_non_edge
  comparison$sd_ratio_edge_to_non_edge <- comparison$sd_edge_only / comparison$sd_non_edge
  comparison
}

run_environment_correlations <- function(data) {
  variables <- c(DIVERSITY_MEASURES, "spec_sa", "spec_spca_mean")
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (edge_group in c("all_quads", "non_edge", "edge_only")) {
      group_data <- data[data$scale == scale_name, ]
      if (edge_group == "non_edge") {
        group_data <- group_data[!group_data$edge_flag, ]
      } else if (edge_group == "edge_only") {
        group_data <- group_data[group_data$edge_flag, ]
      }
      for (env_measure in ENV_MEASURES) {
        for (variable in variables) {
          pair_data <- group_data[!is.na(group_data[[env_measure]]) & !is.na(group_data[[variable]]), ]
          result <- data.frame(
            scale = scale_name,
            edge_group = edge_group,
            env_measure = env_measure,
            env_label = display_name(env_measure),
            variable = variable,
            variable_label = display_name(variable),
            n = nrow(pair_data),
            pearson_r = NA_real_,
            r_squared = NA_real_,
            p_value = NA_real_,
            stringsAsFactors = FALSE
          )
          if (nrow(pair_data) >= 3 &&
              length(unique(pair_data[[env_measure]])) >= 2 &&
              length(unique(pair_data[[variable]])) >= 2) {
            test <- cor.test(pair_data[[env_measure]], pair_data[[variable]], method = "pearson")
            result$pearson_r <- unname(test$estimate)
            result$r_squared <- result$pearson_r^2
            result$p_value <- test$p.value
          }
          rows[[index]] <- result
          index <- index + 1
        }
      }
    }
  }
  do.call(rbind, rows)
}

run_environment_adjusted_models <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (edge_group in c("all_quads", "non_edge", "edge_only")) {
      group_data <- data[data$scale == scale_name, ]
      if (edge_group == "non_edge") {
        group_data <- group_data[!group_data$edge_flag, ]
      } else if (edge_group == "edge_only") {
        group_data <- group_data[group_data$edge_flag, ]
      }
      for (sv_measure in c("spec_sa", "spec_spca_mean")) {
        for (diversity_measure in DIVERSITY_MEASURES) {
          model_data <- group_data[
            complete.cases(group_data[, c(sv_measure, diversity_measure, ENV_MEASURES)]),
          ]
          result <- data.frame(
            scale = scale_name,
            edge_group = edge_group,
            sv_measure = sv_measure,
            sv_label = display_name(sv_measure),
            diversity_measure = diversity_measure,
            diversity_label = display_name(diversity_measure),
            n = nrow(model_data),
            simple_r = NA_real_,
            simple_r2 = NA_real_,
            env_only_r2 = NA_real_,
            adjusted_r2 = NA_real_,
            diversity_incremental_r2_after_env = NA_real_,
            diversity_slope_adjusted = NA_real_,
            diversity_p_adjusted = NA_real_,
            env_elev_p = NA_real_,
            env_tri5_p = NA_real_,
            env_tri11_p = NA_real_,
            env_tri21_p = NA_real_,
            stringsAsFactors = FALSE
          )
          if (nrow(model_data) >= 8 &&
              length(unique(model_data[[sv_measure]])) >= 2 &&
              length(unique(model_data[[diversity_measure]])) >= 2) {
            simple_fit <- lm(model_data[[sv_measure]] ~ model_data[[diversity_measure]])
            env_fit <- lm(model_data[[sv_measure]] ~ env_elev + env_tri5 + env_tri11 + env_tri21, data = model_data)
            adjusted_fit <- lm(
              model_data[[sv_measure]] ~ model_data[[diversity_measure]] + env_elev + env_tri5 + env_tri11 + env_tri21,
              data = model_data
            )
            simple_summary <- summary(simple_fit)
            env_summary <- summary(env_fit)
            adjusted_summary <- summary(adjusted_fit)
            result$simple_r <- cor(model_data[[diversity_measure]], model_data[[sv_measure]])
            result$simple_r2 <- simple_summary$r.squared
            result$env_only_r2 <- env_summary$r.squared
            result$adjusted_r2 <- adjusted_summary$r.squared
            result$diversity_incremental_r2_after_env <- adjusted_summary$r.squared - env_summary$r.squared
            coef_table <- coef(adjusted_summary)
            result$diversity_slope_adjusted <- coef_table[2, "Estimate"]
            result$diversity_p_adjusted <- coef_table[2, "Pr(>|t|)"]
            for (env_measure in ENV_MEASURES) {
              if (env_measure %in% rownames(coef_table)) {
                result[[paste0(env_measure, "_p")]] <- coef_table[env_measure, "Pr(>|t|)"]
              }
            }
          }
          rows[[index]] <- result
          index <- index + 1
        }
      }
    }
  }
  do.call(rbind, rows)
}

save_environment_boxplot <- function(data) {
  file_path <- file.path(FIGURE_DIR, "edge_nonedge_environment_boxplots.png")
  grDevices::png(file_path, width = 2600, height = 2100, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  par(mfrow = c(2, 4), mar = c(5, 4.2, 3.1, 1), oma = c(0, 0, 2.1, 0))
  for (scale_name in SCALES) {
    for (env_measure in ENV_MEASURES) {
      plot_data <- data[data$scale == scale_name & !is.na(data[[env_measure]]), ]
      boxplot(
        plot_data[[env_measure]] ~ ifelse(plot_data$edge_flag, "edge", "non-edge"),
        col = c("#9FC4C6", "#E5A36F"),
        border = "#333333",
        ylab = display_name(env_measure),
        xlab = "",
        main = paste(scale_name, display_name(env_measure))
      )
    }
  }
  mtext("Environmental distributions by edge status", outer = TRUE, cex = 1.15, font = 2)
  file_path
}

write_report <- function(edge_comparison, bootstrap_results, range_comparison, resampling_results, env_comparison, env_correlations, env_adjusted_models, edge_plot_files, bootstrap_plot_files, resampling_plot_files, env_plot_files) {
  strongest_edge <- edge_comparison[
    order(abs(edge_comparison$delta_r_non_edge_minus_all), decreasing = TRUE),
  ][seq_len(min(10, nrow(edge_comparison))), ]
  strongest_edge_table <- data.frame(
    Scale = strongest_edge$scale,
    `SV measure` = strongest_edge$sv_label,
    `Diversity measure` = strongest_edge$diversity_label,
    `n all` = strongest_edge$n_all,
    `edge n all` = strongest_edge$edge_n_all,
    `r all` = fmt_num(strongest_edge$pearson_r_all, 3),
    `r non-edge` = fmt_num(strongest_edge$pearson_r_non_edge, 3),
    `delta r` = fmt_num(strongest_edge$delta_r_non_edge_minus_all, 3),
    check.names = FALSE
  )

  primary_sv <- edge_comparison[
    edge_comparison$sv_measure %in% c("spec_sa", "spec_spca_mean"),
  ]
  primary_sv_table <- data.frame(
    Scale = primary_sv$scale,
    `SV measure` = primary_sv$sv_label,
    `Diversity measure` = primary_sv$diversity_label,
    `r all` = fmt_num(primary_sv$pearson_r_all, 3),
    `R2 all` = fmt_num(primary_sv$r_squared_all, 3),
    `r non-edge` = fmt_num(primary_sv$pearson_r_non_edge, 3),
    `R2 non-edge` = fmt_num(primary_sv$r_squared_non_edge, 3),
    `r edge-only` = fmt_num(primary_sv$pearson_r_edge_only, 3),
    `R2 edge-only` = fmt_num(primary_sv$r_squared_edge_only, 3),
    `delta r` = fmt_num(primary_sv$delta_r_non_edge_minus_all, 3),
    `edge-only minus non-edge r` = fmt_num(primary_sv$delta_r_edge_only_minus_non_edge, 3),
    check.names = FALSE
  )

  edge_only_ranked <- edge_comparison[
    order(abs(edge_comparison$pearson_r_edge_only), decreasing = TRUE),
  ][seq_len(min(12, nrow(edge_comparison))), ]
  edge_only_ranked_table <- data.frame(
    Scale = edge_only_ranked$scale,
    `SV measure` = edge_only_ranked$sv_label,
    `Diversity measure` = edge_only_ranked$diversity_label,
    `n edge-only` = edge_only_ranked$n_edge_only,
    `r edge-only` = fmt_num(edge_only_ranked$pearson_r_edge_only, 3),
    `R2 edge-only` = fmt_num(edge_only_ranked$r_squared_edge_only, 3),
    `F p-value edge-only` = fmt_p(edge_only_ranked$f_p_value_edge_only),
    `r non-edge` = fmt_num(edge_only_ranked$pearson_r_non_edge, 3),
    check.names = FALSE
  )

  primary_variables <- range_comparison[
    range_comparison$variable %in% c(DIVERSITY_MEASURES, "spec_sa", "spec_spca_mean"),
  ]
  range_table <- data.frame(
    Scale = primary_variables$scale,
    Variable = primary_variables$variable_label,
    `n edge` = primary_variables$n_edge_only,
    `edge range` = fmt_num(primary_variables$range_edge_only, 3),
    `non-edge range` = fmt_num(primary_variables$range_non_edge, 3),
    `edge/non-edge range` = fmt_num(primary_variables$range_ratio_edge_to_non_edge, 3),
    `edge SD` = fmt_num(primary_variables$sd_edge_only, 3),
    `non-edge SD` = fmt_num(primary_variables$sd_non_edge, 3),
    `edge/non-edge SD` = fmt_num(primary_variables$sd_ratio_edge_to_non_edge, 3),
    check.names = FALSE
  )

  primary_resampling <- resampling_results[
    resampling_results$sv_measure %in% c("spec_sa", "spec_spca_mean"),
  ]
  resampling_table <- data.frame(
    Scale = primary_resampling$scale,
    `SV measure` = primary_resampling$sv_label,
    `Diversity measure` = primary_resampling$diversity_label,
    `n edge` = primary_resampling$n_edge,
    `edge r` = fmt_num(primary_resampling$edge_r, 3),
    `non-edge equal-n mean r` = fmt_num(primary_resampling$resampled_non_edge_mean_r, 3),
    `non-edge equal-n 2.5% r` = fmt_num(primary_resampling$resampled_non_edge_q025_r, 3),
    `non-edge equal-n 97.5% r` = fmt_num(primary_resampling$resampled_non_edge_q975_r, 3),
    `abs-tail prob` = fmt_num(primary_resampling$abs_edge_r_tail_probability, 3),
    `same-direction tail prob` = fmt_num(primary_resampling$same_direction_tail_probability, 3),
    check.names = FALSE
  )

  env_distribution_table <- data.frame(
    Scale = env_comparison$scale,
    Variable = env_comparison$env_label,
    `edge mean` = fmt_num(env_comparison$mean_edge_only, 3),
    `non-edge mean` = fmt_num(env_comparison$mean_non_edge, 3),
    `edge minus non-edge` = fmt_num(env_comparison$mean_difference_edge_minus_non_edge, 3),
    `edge/non-edge SD` = fmt_num(env_comparison$sd_ratio_edge_to_non_edge, 3),
    check.names = FALSE
  )

  env_sv_correlations <- env_correlations[
    env_correlations$variable %in% c("spec_sa", "spec_spca_mean"),
  ]
  env_sv_correlations <- env_sv_correlations[
    order(abs(env_sv_correlations$pearson_r), decreasing = TRUE),
  ][seq_len(min(16, nrow(env_sv_correlations))), ]
  env_correlation_table <- data.frame(
    Scale = env_sv_correlations$scale,
    Group = env_sv_correlations$edge_group,
    Environment = env_sv_correlations$env_label,
    Variable = env_sv_correlations$variable_label,
    n = env_sv_correlations$n,
    r = fmt_num(env_sv_correlations$pearson_r, 3),
    R2 = fmt_num(env_sv_correlations$r_squared, 3),
    `p-value` = fmt_p(env_sv_correlations$p_value),
    check.names = FALSE
  )

  primary_adjusted <- env_adjusted_models[
    env_adjusted_models$edge_group %in% c("non_edge", "edge_only"),
  ]
  adjusted_table <- data.frame(
    Scale = primary_adjusted$scale,
    Group = primary_adjusted$edge_group,
    `SV measure` = primary_adjusted$sv_label,
    `Diversity measure` = primary_adjusted$diversity_label,
    n = primary_adjusted$n,
    `simple R2` = fmt_num(primary_adjusted$simple_r2, 3),
    `env-only R2` = fmt_num(primary_adjusted$env_only_r2, 3),
    `adjusted R2` = fmt_num(primary_adjusted$adjusted_r2, 3),
    `diversity incremental R2` = fmt_num(primary_adjusted$diversity_incremental_r2_after_env, 3),
    `diversity p adjusted` = fmt_p(primary_adjusted$diversity_p_adjusted),
    check.names = FALSE
  )

  bootstrap_summary <- bootstrap_results[
    bootstrap_results$bootstrap_subset %in% c("non_edge_all_bootstrap", "boot_cv_le_5pct", "boot_cv_gt_5pct", "low_boot_sd", "high_boot_sd"),
  ]
  bootstrap_summary <- data.frame(
    Scale = bootstrap_summary$scale,
    Subset = bootstrap_summary$bootstrap_subset,
    `Diversity measure` = bootstrap_summary$diversity_label,
    n = bootstrap_summary$n,
    `mean boot SD` = fmt_num(bootstrap_summary$mean_boot_sd, 4),
    `mean boot CV` = fmt_num(bootstrap_summary$mean_boot_cv, 4),
    r = fmt_num(bootstrap_summary$pearson_r, 3),
    R2 = fmt_num(bootstrap_summary$r_squared, 3),
    `F p-value` = fmt_p(bootstrap_summary$f_p_value),
    check.names = FALSE
  )

  lines <- c(
    "# Edge Quadrat And Bootstrap Sensitivity For SV-Diversity Correlations",
    "",
    "Date: 2026-07-25",
    "",
    "## Purpose",
    "",
    "This analysis compares how edge quadrats affect correlations between diversity measures and spectral variation at the 10 m and 20 m scales. It also checks whether bootstrap variation in spectral angle entropy changes the apparent strength of biodiversity-SV correlations.",
    "",
    "## Inputs",
    "",
    "- `reports/tables/multiscale_spectral_biodiversity/sv_diversity_analysis_dataset.csv`",
    "- `reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv`",
    "",
    "## Measures",
    "",
    "- Diversity measures: `phy_faith`, `phy_rao`, `phy_afaith`, and `sp_shannon`.",
    "- Spectral variation measures: `spec_sa` and all standardized PCA spectral variation fields in the current analysis dataset.",
    "",
    "## Edge Quadrat Comparison",
    "",
    "Each pairwise correlation was run three ways within 10 m and 20 m quadrats: all quadrats, non-edge quadrats only, and edge quadrats only. Scatter plots show all quadrats, highlight edge quadrats, and include all-quadrat, non-edge, and edge-only regression lines.",
    "",
    "### Primary Measure Summary",
    "",
    markdown_table(primary_sv_table),
    "",
    "### Largest Edge-Removal Changes",
    "",
    markdown_table(strongest_edge_table),
    "",
    "### Strongest Edge-Only Correlations",
    "",
    markdown_table(edge_only_ranked_table),
    "",
    "## Edge Versus Non-Edge Range Check",
    "",
    "This check compares the spread of diversity and spectral variables in edge and non-edge quadrats. Smaller edge-only ranges or standard deviations can make the scatterplots look less diffuse even when the underlying relationship is not fundamentally stronger.",
    "",
    markdown_table(range_table),
    "",
    "## Equal-Sample-Size Non-Edge Resampling",
    "",
    paste0("For each SV-diversity pair, non-edge quadrats were randomly downsampled to the matching edge-only sample size for ", N_RESAMPLE, " iterations. The edge-only correlation was then compared with the distribution of same-sized non-edge correlations."),
    "",
    markdown_table(resampling_table),
    "",
    "## Environmental Context",
    "",
    "This section checks whether elevation or topographic roughness differs between edge and non-edge quadrats, correlates with spectral variation, or reduces the diversity-SV relationship when included in the same model. Environmental adjusted models use `env_elev`, `env_tri5`, `env_tri11`, and `env_tri21` together.",
    "",
    "### Edge Versus Non-Edge Environmental Distributions",
    "",
    markdown_table(env_distribution_table),
    "",
    "### Strongest Environment-SV Correlations",
    "",
    markdown_table(env_correlation_table),
    "",
    "### Diversity Models Before And After Environmental Adjustment",
    "",
    markdown_table(adjusted_table),
    "",
    "## Bootstrap Variation Sensitivity",
    "",
    "Bootstrap sensitivity was evaluated for `spec_sa`, because the bootstrap diagnostics table describes uncertainty in spectral angle entropy. The analysis compares non-edge correlations across all bootstrap-valid quadrats, low-CV subsets, high-CV subsets, and scale-specific bootstrap-SD tertiles.",
    "",
    markdown_table(bootstrap_summary),
    "",
    "## Output Tables",
    "",
    "- `reports/tables/edge_bootstrap_sensitivity/edge_inclusion_correlation_results.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/edge_inclusion_correlation_comparison.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/edge_nonedge_range_comparison.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/equal_n_non_edge_resampling_results.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/edge_nonedge_environment_comparison.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/environment_correlations.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/environment_adjusted_diversity_models.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/bootstrap_variation_correlation_sensitivity.csv`",
    "- `reports/tables/edge_bootstrap_sensitivity/analysis_dataset_with_bootstrap_fields.csv`",
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(c(edge_plot_files, bootstrap_plot_files, resampling_plot_files, env_plot_files)), "`"),
    "",
    "## Interpretation Notes",
    "",
    "- A positive `delta r` means the correlation became stronger or more positive after edge quadrats were removed.",
    "- A negative `delta r` means the correlation became weaker or more negative after edge quadrats were removed.",
    "- `edge-only minus non-edge r` compares the edge-quadrat correlation directly against the non-edge correlation.",
    "- In the equal-n resampling table, small tail probabilities mean the edge-only correlation is unusual relative to randomly selected non-edge subsets of the same size.",
    "- Environmental adjusted models are screening diagnostics, not final spatial inference; the TRI variables are related topographic summaries and should be interpreted cautiously when fitted together.",
    "- Bootstrap sensitivity applies directly to SA entropy only; standardized PCA metrics do not have bootstrap SD/CV fields in the current diagnostics table."
  )
  writeLines(lines, REPORT_PATH)
}

run_analysis <- function() {
  data <- read.csv(INPUT_DATA, stringsAsFactors = FALSE, check.names = FALSE)
  data$quad_id <- as.character(data$quad_id)
  data$edge_flag <- as.logical(data$edge_flag)

  missing_columns <- setdiff(c(DIVERSITY_MEASURES, SV_MEASURES, ENV_MEASURES, "edge_flag", "scale"), names(data))
  if (length(missing_columns) > 0) {
    stop("Missing expected columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  data <- data[data$scale %in% SCALES, ]
  edge_outputs <- run_edge_correlations(data)
  edge_plot_files <- save_all_edge_scatter_plots(data)
  range_outputs <- summarize_edge_nonedge_ranges(data)
  resampling_results <- run_equal_n_resampling(data, edge_outputs$comparison)
  resampling_plot_files <- save_resampling_summary_plot(resampling_results)
  env_comparison <- summarize_environment_edge_nonedge(data)
  env_correlations <- run_environment_correlations(data)
  env_adjusted_models <- run_environment_adjusted_models(data)
  env_plot_files <- save_environment_boxplot(data)

  boot <- read_bootstrap_data()
  data_with_boot <- join_bootstrap_fields(data, boot)
  bootstrap_results <- run_bootstrap_sensitivity(data_with_boot)
  bootstrap_plot_files <- save_bootstrap_plots(data_with_boot)

  write.csv(edge_outputs$results, file.path(TABLE_DIR, "edge_inclusion_correlation_results.csv"), row.names = FALSE)
  write.csv(edge_outputs$comparison, file.path(TABLE_DIR, "edge_inclusion_correlation_comparison.csv"), row.names = FALSE)
  write.csv(range_outputs$summary, file.path(TABLE_DIR, "edge_nonedge_range_summary.csv"), row.names = FALSE)
  write.csv(range_outputs$comparison, file.path(TABLE_DIR, "edge_nonedge_range_comparison.csv"), row.names = FALSE)
  write.csv(resampling_results, file.path(TABLE_DIR, "equal_n_non_edge_resampling_results.csv"), row.names = FALSE)
  write.csv(env_comparison, file.path(TABLE_DIR, "edge_nonedge_environment_comparison.csv"), row.names = FALSE)
  write.csv(env_correlations, file.path(TABLE_DIR, "environment_correlations.csv"), row.names = FALSE)
  write.csv(env_adjusted_models, file.path(TABLE_DIR, "environment_adjusted_diversity_models.csv"), row.names = FALSE)
  write.csv(bootstrap_results, file.path(TABLE_DIR, "bootstrap_variation_correlation_sensitivity.csv"), row.names = FALSE)
  write.csv(data_with_boot, file.path(TABLE_DIR, "analysis_dataset_with_bootstrap_fields.csv"), row.names = FALSE)

  write_report(
    edge_outputs$comparison,
    bootstrap_results,
    range_outputs$comparison,
    resampling_results,
    env_comparison,
    env_correlations,
    env_adjusted_models,
    edge_plot_files,
    bootstrap_plot_files,
    resampling_plot_files,
    env_plot_files
  )

  message("Edge and bootstrap sensitivity analysis complete.")
  invisible(list(
    edge_results = edge_outputs$results,
    edge_comparison = edge_outputs$comparison,
    range_comparison = range_outputs$comparison,
    resampling_results = resampling_results,
    env_comparison = env_comparison,
    env_correlations = env_correlations,
    env_adjusted_models = env_adjusted_models,
    bootstrap_results = bootstrap_results
  ))
}

if (sys.nframe() == 0) {
  run_analysis()
}
