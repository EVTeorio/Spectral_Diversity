PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-08-17"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/bootstrap_metric_relationships")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/bootstrap_metric_relationships")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260817_bootstrap_metric_stat_relationship_analysis.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260817_bootstrap_metric_stat_relationship_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260817_bootstrap_metric_stat_relationship_analysis_validation.md")

SCALES <- c("10m", "20m", "50m")
SCALE_LABELS <- c("10m" = "10 m", "20m" = "20 m", "50m" = "50 m")
STATS <- c("metric_sd", "metric_cv", "metric_se")
STAT_LABELS <- c(
  metric_sd = "Bootstrap SD",
  metric_cv = "Bootstrap CV",
  metric_se = "Bootstrap SE"
)

SAMPLE_SIZE_SUMMARIES <- data.frame(
  metric_id = c(
    "sa_entropy",
    "pca_mean_distance",
    "spectral_rao_q",
    "alpha_hull_area",
    "standardized_PCA_pca_mean_distance",
    "standardized_PCA_spectral_rao_q",
    "standardized_PCA_alpha_hull_area"
  ),
  metric_label = c(
    "Spectral angle entropy",
    "Raw PCA mean distance",
    "Raw PCA Rao's Q",
    "Raw PCA alpha hull",
    "Std PCA mean distance",
    "Std PCA Rao's Q",
    "Std PCA alpha hull"
  ),
  path = file.path(
    PROJECT_DIR,
    "reports/tables/sample_size_effects",
    c(
      "sa_entropy/sa_entropy_sample_size_summary.csv",
      "pca_mean_distance/pca_mean_distance_sample_size_summary.csv",
      "spectral_rao_q/spectral_rao_q_sample_size_summary.csv",
      "alpha_hull_area/alpha_hull_area_sample_size_summary.csv",
      "standardized_PCA_pca_mean_distance/standardized_PCA_pca_mean_distance_sample_size_summary.csv",
      "standardized_PCA_spectral_rao_q/standardized_PCA_spectral_rao_q_sample_size_summary.csv",
      "standardized_PCA_alpha_hull_area/standardized_PCA_alpha_hull_area_sample_size_summary.csv"
    )
  ),
  stringsAsFactors = FALSE
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

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

relationship_row <- function(data, x_col, y_col) {
  complete <- !is.na(data[[x_col]]) & !is.na(data[[y_col]])
  result <- data.frame(
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
    model_data <- data[complete, c(x_col, y_col)]
    names(model_data) <- c("x", "y")
    fit <- lm(y ~ x, data = model_data)
    fit_summary <- summary(fit)
    f_values <- fit_summary$fstatistic
    pearson <- cor.test(model_data$x, model_data$y, method = "pearson")
    result$pearson_r <- unname(pearson$estimate)
    result$r_squared <- fit_summary$r.squared
    result$f_statistic <- unname(f_values[["value"]])
    result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
    result$slope <- unname(coef(fit)[["x"]])
    result$intercept <- unname(coef(fit)[["(Intercept)"]])
  }
  result
}

read_final_sa_bootstrap <- function() {
  path <- file.path(PROJECT_DIR, "reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv")
  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  data$scale <- gsub(" ", "", data$scale)
  data.frame(
    source = "final_sa_entropy",
    scale = data$scale,
    quad_id = as.character(data$quad_id),
    metric_id = "sa_entropy",
    metric_label = "Spectral angle entropy",
    sample_rule = "final_70_iteration_bootstrap",
    n_boot = data$n_boot,
    metric_value = data$spectral_entropy,
    metric_sd = data$boot_sd,
    metric_cv = data$boot_cv,
    metric_se = data$boot_se,
    stringsAsFactors = FALSE
  )
}

read_fixed4000_sample_bootstrap <- function() {
  rows <- vector("list", nrow(SAMPLE_SIZE_SUMMARIES))
  for (i in seq_len(nrow(SAMPLE_SIZE_SUMMARIES))) {
    config <- SAMPLE_SIZE_SUMMARIES[i, ]
    data <- read.csv(config$path, stringsAsFactors = FALSE, check.names = FALSE)
    data <- data[data$sample_rule == "fixed_4000", , drop = FALSE]
    if (config$metric_id == "sa_entropy") {
      metric_value <- data$entropy_mean
      metric_sd <- data$entropy_sd
      metric_cv <- data$entropy_cv
      metric_se <- data$entropy_se
    } else {
      metric_value <- data$metric_mean
      metric_sd <- data$metric_sd
      metric_cv <- data$metric_cv
      metric_se <- data$metric_se
    }
    rows[[i]] <- data.frame(
      source = "sample_size_fixed4000",
      scale = data$scale,
      quad_id = as.character(data$quad_id),
      metric_id = config$metric_id,
      metric_label = config$metric_label,
      sample_rule = data$sample_rule,
      n_boot = data$n_boot,
      metric_value = metric_value,
      metric_sd = metric_sd,
      metric_cv = metric_cv,
      metric_se = metric_se,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_relationships <- function(data) {
  rows <- list()
  index <- 1
  for (source_name in unique(data$source)) {
    source_data <- data[data$source == source_name, , drop = FALSE]
    for (scale_name in SCALES) {
      scale_data <- source_data[source_data$scale == scale_name, , drop = FALSE]
      for (metric_id in unique(scale_data$metric_id)) {
        metric_data <- scale_data[scale_data$metric_id == metric_id, , drop = FALSE]
        metric_label <- unique(metric_data$metric_label)[1]
        for (stat in STATS) {
          rel <- relationship_row(metric_data, "metric_value", stat)
          rows[[index]] <- cbind(
            data.frame(
              source = source_name,
              scale = scale_name,
              metric_id = metric_id,
              metric_label = metric_label,
              statistic = stat,
              statistic_label = STAT_LABELS[[stat]],
              stringsAsFactors = FALSE
            ),
            rel
          )
          index <- index + 1
        }
      }
    }
  }
  do.call(rbind, rows)
}

plot_panel <- function(panel_data, metric_label, stat, relationship_table, source_name, scale_name) {
  plot(
    panel_data$metric_value,
    panel_data[[stat]],
    pch = 16,
    cex = 0.58,
    col = adjustcolor("#315C8A", alpha.f = 0.62),
    xlab = metric_label,
    ylab = STAT_LABELS[[stat]],
    main = paste(metric_label, "vs", STAT_LABELS[[stat]])
  )
  row <- relationship_table[
    relationship_table$source == source_name &
      relationship_table$scale == scale_name &
      relationship_table$metric_label == metric_label &
      relationship_table$statistic == stat,
  ]
  complete <- !is.na(panel_data$metric_value) & !is.na(panel_data[[stat]])
  if (nrow(row) == 1 && !is.na(row$pearson_r) && sum(complete) >= 3) {
    fit <- lm(panel_data[[stat]][complete] ~ panel_data$metric_value[complete])
    abline(fit, col = "#111111", lwd = 2)
    legend(
      "topleft",
      legend = c(
        paste0("n=", row$n),
        paste0("r=", fmt_num(row$pearson_r, 2)),
        paste0("R2=", fmt_num(row$r_squared, 2)),
        paste0("p=", fmt_p(row$f_p_value))
      ),
      bty = "n",
      cex = 0.72
    )
  }
}

save_final_sa_figures <- function(final_sa, relationships) {
  files <- character(length(SCALES))
  for (i in seq_along(SCALES)) {
    scale_name <- SCALES[i]
    file_path <- file.path(FIGURE_DIR, paste0("01_final_sa_entropy_vs_bootstrap_stats_", scale_name, ".png"))
    plot_data <- final_sa[final_sa$scale == scale_name, , drop = FALSE]
    png(file_path, width = 2700, height = 900, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
    for (stat in STATS) {
      plot_panel(plot_data, "Spectral angle entropy", stat, relationships, "final_sa_entropy", scale_name)
    }
    mtext(paste0(scale_name, " final SA entropy vs 70-iteration bootstrap statistics"), outer = TRUE, cex = 1.1, font = 2)
    par(old_par)
    dev.off()
    files[i] <- file_path
  }
  files
}

save_fixed4000_figures <- function(fixed4000, relationships) {
  files <- character(length(SCALES))
  metric_order <- SAMPLE_SIZE_SUMMARIES$metric_label
  for (i in seq_along(SCALES)) {
    scale_name <- SCALES[i]
    file_path <- file.path(FIGURE_DIR, paste0("02_fixed4000_metrics_vs_bootstrap_stats_", scale_name, ".png"))
    png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
    for (metric_label in metric_order) {
      plot_data <- fixed4000[fixed4000$scale == scale_name & fixed4000$metric_label == metric_label, , drop = FALSE]
      for (stat in STATS) {
        plot_panel(plot_data, metric_label, stat, relationships, "sample_size_fixed4000", scale_name)
      }
    }
    mtext(paste0(scale_name, " fixed-4,000 sample-size bootstrap means vs statistics"), outer = TRUE, cex = 1.1, font = 2)
    par(old_par)
    dev.off()
    files[i] <- file_path
  }
  files
}

write_reports <- function(final_sa, fixed4000, relationships, figure_files) {
  strongest <- relationships[order(relationships$source, relationships$scale, -abs(relationships$pearson_r)), ]
  strongest <- do.call(rbind, lapply(split(strongest, paste(strongest$source, strongest$scale)), function(x) x[seq_len(min(3, nrow(x))), ]))
  strongest_table <- data.frame(
    Source = strongest$source,
    Scale = strongest$scale,
    Metric = strongest$metric_label,
    Statistic = strongest$statistic_label,
    n = strongest$n,
    r = fmt_num(strongest$pearson_r, 3),
    R2 = fmt_num(strongest$r_squared, 3),
    `p-value` = fmt_p(strongest$f_p_value),
    check.names = FALSE
  )

  lines <- c(
    "# Bootstrap Metric Statistic Relationship Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis plots spectral metric values against bootstrap variability statistics: SD, CV, and SE.",
    "",
    "## Data Sources",
    "",
    "- Final-pipeline SA entropy bootstrap diagnostics use `reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv` and summarize the 70-iteration SA entropy bootstrap.",
    "- All-metric bootstrap-stat panels use the fixed-4,000-pixel condition from the sample-size sensitivity branch. This is intentionally labeled separately from the final pipeline because the current visible main pipeline does not contain per-iteration PCA metric CSVs.",
    "",
    "## Strongest Relationships",
    "",
    markdown_table(strongest_table),
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Output Tables",
    "",
    "- `reports/tables/bootstrap_metric_relationships/final_sa_entropy_bootstrap_metric_stats.csv`",
    "- `reports/tables/bootstrap_metric_relationships/fixed4000_sample_size_bootstrap_metric_stats.csv`",
    "- `reports/tables/bootstrap_metric_relationships/bootstrap_metric_stat_relationships.csv`"
  )
  writeLines(lines, ANALYSIS_REPORT)

  task_lines <- c(
    "# Task Report: Bootstrap Metric Statistic Relationship Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Objective",
    "",
    "Create scale-specific scatterplot figures showing bootstrapped spectral metric values against metric SD, metric CV, and metric SE.",
    "",
    "## Outputs Created",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "- `reports/tables/bootstrap_metric_relationships/final_sa_entropy_bootstrap_metric_stats.csv`",
    "- `reports/tables/bootstrap_metric_relationships/fixed4000_sample_size_bootstrap_metric_stats.csv`",
    "- `reports/tables/bootstrap_metric_relationships/bootstrap_metric_stat_relationships.csv`",
    "",
    "## Notes",
    "",
    "- Final-pipeline bootstrap-stat figures are available for SA entropy.",
    "- All-seven-metric figures use the fixed-4,000-pixel condition from the sample-size branch and should be interpreted as bootstrap stability diagnostics, not final-value uncertainty for PCA metrics."
  )
  writeLines(task_lines, TASK_REPORT)

  validation_lines <- c(
    "# Validation: Bootstrap Metric Statistic Relationship Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Checks",
    "",
    paste0("- Final SA entropy rows: ", nrow(final_sa)),
    paste0("- Fixed-4,000 all-metric rows: ", nrow(fixed4000)),
    paste0("- Relationship summary rows: ", nrow(relationships)),
    paste0("- Missing relationship Pearson r values: ", sum(is.na(relationships$pearson_r))),
    paste0("- Output figures present: ", sum(file.exists(figure_files)), " of ", length(figure_files)),
    paste0("- Output tables present: ", sum(file.exists(file.path(TABLE_DIR, c(
      "final_sa_entropy_bootstrap_metric_stats.csv",
      "fixed4000_sample_size_bootstrap_metric_stats.csv",
      "bootstrap_metric_stat_relationships.csv"
    )))), " of 3"),
    "",
    "## Result",
    "",
    "The requested bootstrap metric-statistic scatterplot figures and companion tables were produced."
  )
  writeLines(validation_lines, VALIDATION_REPORT)
}

run_bootstrap_metric_stat_relationship_analysis <- function() {
  final_sa <- read_final_sa_bootstrap()
  fixed4000 <- read_fixed4000_sample_bootstrap()
  all_data <- rbind(final_sa, fixed4000)
  relationships <- summarize_relationships(all_data)
  figure_files <- c(
    save_final_sa_figures(final_sa, relationships),
    save_fixed4000_figures(fixed4000, relationships)
  )

  write.csv(final_sa, file.path(TABLE_DIR, "final_sa_entropy_bootstrap_metric_stats.csv"), row.names = FALSE)
  write.csv(fixed4000, file.path(TABLE_DIR, "fixed4000_sample_size_bootstrap_metric_stats.csv"), row.names = FALSE)
  write.csv(relationships, file.path(TABLE_DIR, "bootstrap_metric_stat_relationships.csv"), row.names = FALSE)
  write_reports(final_sa, fixed4000, relationships, figure_files)

  message("Bootstrap metric-stat relationship analysis complete.")
  invisible(list(
    final_sa = final_sa,
    fixed4000 = fixed4000,
    relationships = relationships,
    figure_files = figure_files
  ))
}

if (sys.nframe() == 0) {
  run_bootstrap_metric_stat_relationship_analysis()
}
