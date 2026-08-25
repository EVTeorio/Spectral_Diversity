PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-08-19"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/spectral_biodiversity_all_metrics")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/spectral_biodiversity_all_metrics")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260819_spectral_biodiversity_all_metric_scatter_analysis.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260819_spectral_biodiversity_all_metric_scatter_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260819_spectral_biodiversity_all_metric_scatter_analysis_validation.md")

SCALES <- c("10m", "20m", "50m")
SCALE_LABELS <- c("10m" = "10 m", "20m" = "20 m", "50m" = "50 m")

SPECTRAL_METRICS <- c(
  "spec_sa",
  "spec_alpha",
  "spec_rao_q",
  "spec_pca_mean",
  "spec_spca_alpha",
  "spec_spca_rao",
  "spec_spca_mean"
)

BIODIVERSITY_METRICS <- c(
  "phy_faith",
  "phy_rao",
  "phy_afaith",
  "sp_rich",
  "sp_shannon",
  "sp_simpson",
  "sp_even"
)

DISPLAY_NAMES <- c(
  spec_sa = "Spectral angle entropy",
  spec_alpha = "Raw PCA alpha hull",
  spec_rao_q = "Raw PCA Rao's Q",
  spec_pca_mean = "Raw PCA mean distance",
  spec_spca_alpha = "Std PCA alpha hull",
  spec_spca_rao = "Std PCA Rao's Q",
  spec_spca_mean = "Std PCA mean distance",
  phy_faith = "Faith's PD",
  phy_rao = "Phylogenetic Rao's Q",
  phy_afaith = "Abundance-weighted Faith's PD",
  sp_rich = "Species richness",
  sp_shannon = "Shannon diversity",
  sp_simpson = "Simpson diversity",
  sp_even = "Species evenness"
)

SHORT_NAMES <- c(
  spec_sa = "sa_entropy",
  spec_alpha = "raw_pca_alpha_hull",
  spec_rao_q = "raw_pca_rao_q",
  spec_pca_mean = "raw_pca_mean_distance",
  spec_spca_alpha = "std_pca_alpha_hull",
  spec_spca_rao = "std_pca_rao_q",
  spec_spca_mean = "std_pca_mean_distance"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

display_name <- function(metric) unname(DISPLAY_NAMES[metric])

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

read_scale_table <- function(scale_name) {
  file_path <- file.path(PROJECT_DIR, paste0("quadrat_analysis_", scale_name, ".csv"))
  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  data$quad_id <- as.character(data$quad_id)
  data$scale <- as.character(data$scale)
  data
}

load_analysis_data <- function() {
  data <- do.call(rbind, lapply(SCALES, read_scale_table))
  required <- c("quad_id", "scale", SPECTRAL_METRICS, BIODIVERSITY_METRICS)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data[, required, drop = FALSE]
}

relationship_row <- function(data, scale_name, spectral_metric, biodiversity_metric) {
  scale_data <- data[data$scale == scale_name, , drop = FALSE]
  complete <- !is.na(scale_data[[spectral_metric]]) & !is.na(scale_data[[biodiversity_metric]])
  row <- data.frame(
    scale = scale_name,
    spectral_metric = spectral_metric,
    spectral_label = display_name(spectral_metric),
    biodiversity_metric = biodiversity_metric,
    biodiversity_label = display_name(biodiversity_metric),
    n = sum(complete),
    pearson_r = NA_real_,
    r_squared = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    stringsAsFactors = FALSE
  )
  if (sum(complete) >= 3) {
    model_data <- scale_data[complete, c(biodiversity_metric, spectral_metric)]
    names(model_data) <- c("x", "y")
    fit <- lm(y ~ x, data = model_data)
    fit_summary <- summary(fit)
    f_values <- fit_summary$fstatistic
    pearson <- cor.test(model_data$x, model_data$y, method = "pearson")
    row$pearson_r <- unname(pearson$estimate)
    row$r_squared <- fit_summary$r.squared
    row$f_statistic <- unname(f_values[["value"]])
    row$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
    row$slope <- unname(coef(fit)[["x"]])
    row$intercept <- unname(coef(fit)[["(Intercept)"]])
  }
  row
}

summarize_relationships <- function(data) {
  rows <- list()
  index <- 1
  for (spectral_metric in SPECTRAL_METRICS) {
    for (biodiversity_metric in BIODIVERSITY_METRICS) {
      for (scale_name in SCALES) {
        rows[[index]] <- relationship_row(data, scale_name, spectral_metric, biodiversity_metric)
        index <- index + 1
      }
    }
  }
  do.call(rbind, rows)
}

plot_panel <- function(data, relationships, scale_name, spectral_metric, biodiversity_metric) {
  scale_data <- data[data$scale == scale_name, , drop = FALSE]
  complete <- !is.na(scale_data[[spectral_metric]]) & !is.na(scale_data[[biodiversity_metric]])
  point_col <- if (grepl("^phy_", biodiversity_metric)) "#315C8A" else "#3D7B45"
  plot(
    scale_data[[biodiversity_metric]],
    scale_data[[spectral_metric]],
    pch = 16,
    cex = 0.50,
    col = adjustcolor(point_col, alpha.f = 0.58),
    xlab = display_name(biodiversity_metric),
    ylab = display_name(spectral_metric),
    main = SCALE_LABELS[[scale_name]]
  )

  row <- relationships[
    relationships$scale == scale_name &
      relationships$spectral_metric == spectral_metric &
      relationships$biodiversity_metric == biodiversity_metric,
  ]
  if (nrow(row) == 1 && !is.na(row$pearson_r) && sum(complete) >= 3) {
    fit <- lm(scale_data[[spectral_metric]][complete] ~ scale_data[[biodiversity_metric]][complete])
    abline(fit, col = "#111111", lwd = 2)
    legend(
      "topleft",
      legend = c(
        paste0("n=", row$n),
        paste0("r=", fmt_num(row$pearson_r, 3)),
        paste0("R2=", fmt_num(row$r_squared, 3)),
        paste0("p=", fmt_p(row$f_p_value))
      ),
      bty = "n",
      cex = 0.72
    )
  }
}

save_spectral_metric_figures <- function(data, relationships) {
  files <- character(length(SPECTRAL_METRICS))
  for (i in seq_along(SPECTRAL_METRICS)) {
    spectral_metric <- SPECTRAL_METRICS[i]
    file_path <- file.path(
      FIGURE_DIR,
      paste0("01_", SHORT_NAMES[[spectral_metric]], "_vs_all_biodiversity_metrics_by_scale.png")
    )
    png(file_path, width = 3600, height = 6200, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.5, 2.3, 0.8), oma = c(0, 0, 2.8, 0))
    for (biodiversity_metric in BIODIVERSITY_METRICS) {
      for (scale_name in SCALES) {
        plot_panel(data, relationships, scale_name, spectral_metric, biodiversity_metric)
      }
    }
    mtext(
      paste0(display_name(spectral_metric), " vs biodiversity metrics across scales"),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    par(old_par)
    dev.off()
    files[i] <- file_path
  }
  files
}

write_reports <- function(data, relationships, figure_files) {
  top <- relationships[order(-abs(relationships$pearson_r)), ]
  top <- head(top, 20)
  top_table <- data.frame(
    Scale = top$scale,
    Spectral = top$spectral_label,
    Biodiversity = top$biodiversity_label,
    n = top$n,
    r = fmt_num(top$pearson_r, 3),
    R2 = fmt_num(top$r_squared, 3),
    `p-value` = fmt_p(top$f_p_value),
    check.names = FALSE
  )

  missing_summary <- do.call(rbind, lapply(SCALES, function(scale_name) {
    scale_data <- data[data$scale == scale_name, , drop = FALSE]
    data.frame(
      scale = scale_name,
      rows = nrow(scale_data),
      min_complete_n = min(relationships$n[relationships$scale == scale_name]),
      max_complete_n = max(relationships$n[relationships$scale == scale_name]),
      stringsAsFactors = FALSE
    )
  }))

  analysis_lines <- c(
    "# Spectral Biodiversity All-Metric Scatter Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis creates one scatterplot figure for each of the seven focal spectral heterogeneity measures. In each figure, rows are the seven biodiversity metrics and columns are the 10 m, 20 m, and 50 m quadrat scales.",
    "",
    "## Data",
    "",
    "- Input tables are `quadrat_analysis_10m.csv`, `quadrat_analysis_20m.csv`, and `quadrat_analysis_50m.csv`.",
    "- Each panel uses all complete quadrat records for that spectral-biodiversity metric pair; no edge filtering is applied.",
    "- Panels report sample size, Pearson r, linear-model R2, and the simple-regression F-test p-value.",
    "",
    "## Strongest Absolute Correlations",
    "",
    markdown_table(top_table),
    "",
    "## Complete Case Ranges",
    "",
    markdown_table(missing_summary),
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Tables",
    "",
    "- `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_dataset.csv`",
    "- `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_relationships.csv`"
  )
  writeLines(analysis_lines, ANALYSIS_REPORT)

  task_lines <- c(
    "# Task Report: Spectral Biodiversity All-Metric Scatter Figures",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Objective",
    "",
    "Create seven scale-column scatterplot figures, one for each spectral heterogeneity measure, comparing that metric with all seven biodiversity metrics.",
    "",
    "## Outputs",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "- `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_dataset.csv`",
    "- `reports/tables/spectral_biodiversity_all_metrics/spectral_biodiversity_all_metric_relationships.csv`"
  )
  writeLines(task_lines, TASK_REPORT)

  validation_lines <- c(
    "# Validation: Spectral Biodiversity All-Metric Scatter Figures",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Checks",
    "",
    paste0("- Dataset rows: ", nrow(data)),
    paste0("- Relationship rows: ", nrow(relationships)),
    paste0("- Expected relationship rows: ", length(SPECTRAL_METRICS) * length(BIODIVERSITY_METRICS) * length(SCALES)),
    paste0("- Missing Pearson r values: ", sum(is.na(relationships$pearson_r))),
    paste0("- Output figures present: ", sum(file.exists(figure_files)), " of ", length(figure_files)),
    paste0("- Output tables present: ", sum(file.exists(file.path(TABLE_DIR, c(
      "spectral_biodiversity_all_metric_dataset.csv",
      "spectral_biodiversity_all_metric_relationships.csv"
    )))), " of 2"),
    "",
    "## Result",
    "",
    "The requested all-metric spectral-biodiversity scatterplot figures and companion tables were produced."
  )
  writeLines(validation_lines, VALIDATION_REPORT)
}

run_spectral_biodiversity_all_metric_scatter_analysis <- function() {
  data <- load_analysis_data()
  relationships <- summarize_relationships(data)
  figure_files <- save_spectral_metric_figures(data, relationships)

  write.csv(data, file.path(TABLE_DIR, "spectral_biodiversity_all_metric_dataset.csv"), row.names = FALSE)
  write.csv(relationships, file.path(TABLE_DIR, "spectral_biodiversity_all_metric_relationships.csv"), row.names = FALSE)
  write_reports(data, relationships, figure_files)

  message("Spectral biodiversity all-metric scatter analysis complete.")
  invisible(list(data = data, relationships = relationships, figure_files = figure_files))
}

if (sys.nframe() == 0) {
  run_spectral_biodiversity_all_metric_scatter_analysis()
}
