PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-08-31"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/spectral_axis_power2_shannon_phylogenetic")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/spectral_axis_power2_shannon_phylogenetic")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260831_spectral_axis_power2_shannon_phylogenetic_side_investigation.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260831_spectral_axis_power2_shannon_phylogenetic_side_investigation_validation.md")

SCALES <- c("10m", "20m", "50m")
SCALE_LABELS <- c("10m" = "10 m", "20m" = "20 m", "50m" = "50 m")

SPECTRAL_METRICS <- c(
  "spec_sa",
  "spec_alpha",
  "spec_rao_q",
  "spec_pca_mean",
  "spec_convex",
  "spec_hull3d_v",
  "spec_hull3d_a",
  "spec_spca_alpha",
  "spec_spca_rao",
  "spec_spca_mean",
  "spec_spca_convex",
  "spec_spca_hull3d_v",
  "spec_spca_hull3d_a"
)

BIODIVERSITY_METRICS <- c(
  "phy_faith",
  "phy_rao",
  "phy_afaith",
  "sp_shannon"
)

DISPLAY_NAMES <- c(
  spec_sa = "Spectral angle entropy",
  spec_alpha = "Raw PCA alpha hull",
  spec_rao_q = "Raw PCA Rao's Q",
  spec_pca_mean = "Raw PCA mean distance",
  spec_convex = "Raw PCA convex hull",
  spec_hull3d_v = "Raw PCA 3D hull volume",
  spec_hull3d_a = "Raw PCA 3D hull area",
  spec_spca_alpha = "Std PCA alpha hull",
  spec_spca_rao = "Std PCA Rao's Q",
  spec_spca_mean = "Std PCA mean distance",
  spec_spca_convex = "Std PCA convex hull",
  spec_spca_hull3d_v = "Std PCA 3D hull volume",
  spec_spca_hull3d_a = "Std PCA 3D hull area",
  phy_faith = "Faith's PD",
  phy_rao = "Phylogenetic Rao's Q",
  phy_afaith = "Abundance-weighted Faith's PD",
  sp_shannon = "Shannon diversity"
)

SHORT_NAMES <- c(
  spec_sa = "sa_entropy",
  spec_alpha = "raw_pca_alpha_hull",
  spec_rao_q = "raw_pca_rao_q",
  spec_pca_mean = "raw_pca_mean_distance",
  spec_convex = "raw_pca_convex_hull",
  spec_hull3d_v = "raw_pca_3d_hull_volume",
  spec_hull3d_a = "raw_pca_3d_hull_area",
  spec_spca_alpha = "std_pca_alpha_hull",
  spec_spca_rao = "std_pca_rao_q",
  spec_spca_mean = "std_pca_mean_distance",
  spec_spca_convex = "std_pca_convex_hull",
  spec_spca_hull3d_v = "std_pca_3d_hull_volume",
  spec_spca_hull3d_a = "std_pca_3d_hull_area"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

display_name <- function(metric) {
  ifelse(metric %in% names(DISPLAY_NAMES), unname(DISPLAY_NAMES[metric]), metric)
}

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
  escaped_project <- gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", normalized_project)
  sub(paste0("^", escaped_project, "/"), "", normalized_paths)
}

power2_values <- function(values) {
  as.numeric(values)^2
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

safe_cor <- function(x, y, method = "pearson") {
  complete <- is.finite(x) & is.finite(y)
  if (sum(complete) < 3 || stats::sd(x[complete]) == 0 || stats::sd(y[complete]) == 0) {
    return(NA_real_)
  }
  suppressWarnings(unname(stats::cor(x[complete], y[complete], method = method)))
}

relationship_row <- function(data, scale_name, spectral_metric, biodiversity_metric) {
  scale_data <- data[data$scale == scale_name, , drop = FALSE]
  x <- as.numeric(scale_data[[biodiversity_metric]])
  y <- power2_values(scale_data[[spectral_metric]])
  complete <- is.finite(x) & is.finite(y)

  row <- data.frame(
    scenario = "spectral_power2_biodiversity_identity",
    scale = scale_name,
    spectral_metric = spectral_metric,
    spectral_label = display_name(spectral_metric),
    spectral_transform = "power_2",
    biodiversity_metric = biodiversity_metric,
    biodiversity_label = display_name(biodiversity_metric),
    biodiversity_transform = "identity",
    model_formula = "spectral_metric^2 ~ biodiversity_metric",
    x_axis = "biodiversity_identity",
    y_axis = "spectral_power_2",
    n = sum(complete),
    pearson_r = NA_real_,
    spearman_r = NA_real_,
    r_squared = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    stringsAsFactors = FALSE
  )

  if (sum(complete) >= 4 && stats::sd(x[complete]) > 0 && stats::sd(y[complete]) > 0) {
    fit <- stats::lm(y[complete] ~ x[complete])
    fit_summary <- summary(fit)
    f_values <- fit_summary$fstatistic
    row$pearson_r <- safe_cor(x, y, "pearson")
    row$spearman_r <- safe_cor(x, y, "spearman")
    row$r_squared <- fit_summary$r.squared
    row$f_statistic <- unname(f_values[["value"]])
    row$f_p_value <- stats::pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
    row$slope <- unname(stats::coef(fit)[[2]])
    row$intercept <- unname(stats::coef(fit)[[1]])
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
  out <- do.call(rbind, rows)
  out$spectral_group <- ifelse(
    out$spectral_metric == "spec_sa",
    "spectral_angle_entropy",
    ifelse(grepl("^spec_spca", out$spectral_metric), "standardized_pca", "raw_pca")
  )
  out$biodiversity_group <- ifelse(grepl("^phy_", out$biodiversity_metric), "phylogenetic", "shannon")
  out$abs_pearson_r <- abs(out$pearson_r)
  out$spectral_order <- match(out$spectral_metric, SPECTRAL_METRICS)
  out$biodiversity_order <- match(out$biodiversity_metric, BIODIVERSITY_METRICS)
  out$scale_order <- match(out$scale, SCALES)
  out <- out[order(out$spectral_order, out$biodiversity_order, out$scale_order), ]
  out[, setdiff(names(out), c("spectral_order", "biodiversity_order", "scale_order")), drop = FALSE]
}

axis_label <- function(metric, transformed = FALSE) {
  label <- display_name(metric)
  if (transformed) {
    paste0(label, "^2")
  } else {
    label
  }
}

plot_panel <- function(data, relationships, scale_name, spectral_metric, biodiversity_metric) {
  scale_data <- data[data$scale == scale_name, , drop = FALSE]
  x <- as.numeric(scale_data[[biodiversity_metric]])
  y <- power2_values(scale_data[[spectral_metric]])
  complete <- is.finite(x) & is.finite(y)

  point_col <- if (grepl("^phy_", biodiversity_metric)) "#315C8A" else "#3D7B45"
  plot(
    x,
    y,
    pch = 16,
    cex = 0.52,
    col = grDevices::adjustcolor(point_col, alpha.f = 0.58),
    xlab = axis_label(biodiversity_metric, transformed = FALSE),
    ylab = axis_label(spectral_metric, transformed = TRUE),
    main = SCALE_LABELS[[scale_name]]
  )

  row <- relationships[
    relationships$scale == scale_name &
      relationships$spectral_metric == spectral_metric &
      relationships$biodiversity_metric == biodiversity_metric,
  ]

  if (nrow(row) == 1 && !is.na(row$pearson_r) && sum(complete) >= 4) {
    fit <- stats::lm(y[complete] ~ x[complete])
    graphics::abline(fit, col = "#111111", lwd = 2)
    graphics::legend(
      "topleft",
      legend = c(
        paste0("n=", row$n),
        paste0("r=", fmt_num(row$pearson_r, 3)),
        paste0("R2=", fmt_num(row$r_squared, 3)),
        paste0("p=", fmt_p(row$f_p_value))
      ),
      bty = "o",
      bg = grDevices::adjustcolor("white", alpha.f = 0.88),
      box.col = NA,
      cex = 0.68
    )
  }
}

save_spectral_metric_figures <- function(data, relationships) {
  files <- character(length(SPECTRAL_METRICS))
  for (i in seq_along(SPECTRAL_METRICS)) {
    spectral_metric <- SPECTRAL_METRICS[i]
    file_path <- file.path(
      FIGURE_DIR,
      paste0("spectral_power2_", SHORT_NAMES[[spectral_metric]], "_vs_shannon_phylogenetic_by_scale.png")
    )

    grDevices::png(file_path, width = 3600, height = 3600, res = 300)
    old_par <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow = c(length(BIODIVERSITY_METRICS), length(SCALES)), mar = c(4.8, 5.0, 2.5, 0.8), oma = c(0, 0, 3.0, 0))
    for (biodiversity_metric in BIODIVERSITY_METRICS) {
      for (scale_name in SCALES) {
        plot_panel(data, relationships, scale_name, spectral_metric, biodiversity_metric)
      }
    }
    graphics::mtext(
      paste0(axis_label(spectral_metric, transformed = TRUE), " vs Shannon and phylogenetic diversity"),
      outer = TRUE,
      cex = 1.15,
      font = 2
    )
    graphics::par(old_par)
    grDevices::dev.off()
    files[i] <- file_path
  }
  files
}

write_reports <- function(data, relationships, figure_files) {
  top <- relationships[order(-relationships$abs_pearson_r, relationships$scale), ]
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

  complete_summary <- do.call(rbind, lapply(SCALES, function(scale_name) {
    scale_relationships <- relationships[relationships$scale == scale_name, , drop = FALSE]
    scale_data <- data[data$scale == scale_name, , drop = FALSE]
    data.frame(
      scale = scale_name,
      rows = nrow(scale_data),
      min_complete_n = min(scale_relationships$n),
      max_complete_n = max(scale_relationships$n),
      stringsAsFactors = FALSE
    )
  }))

  analysis_lines <- c(
    "# Side Investigation: Spectral-Axis Power-2 Transform",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This side investigation repeats the spectral-biodiversity scatter comparison with a power-2 transformation applied only to the spectral metric axis. Biodiversity metrics remain on their original scale.",
    "",
    "## Scope",
    "",
    "- Spectral metrics: spectral angle entropy plus raw and standardized PCA alpha hull, Rao's Q, mean distance, convex hull, 3D hull volume, and 3D hull area.",
    "- Biodiversity metrics: Faith's PD, phylogenetic Rao's Q, abundance-weighted Faith's PD, and Shannon diversity.",
    "- Spatial scales: 10 m, 20 m, and 50 m.",
    "- Each panel reports sample size, Pearson r, linear-model R2, and the simple-regression F-test p-value.",
    "",
    "## Strongest Absolute Pearson Correlations",
    "",
    markdown_table(top_table),
    "",
    "## Complete Case Ranges",
    "",
    markdown_table(complete_summary),
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Table",
    "",
    "- `reports/tables/spectral_axis_power2_shannon_phylogenetic/spectral_axis_power2_shannon_phylogenetic_summary.csv`"
  )
  writeLines(analysis_lines, ANALYSIS_REPORT)

  validation_lines <- c(
    "# Validation: Spectral-Axis Power-2 Shannon/Phylogenetic Side Investigation",
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
    paste0("- Output table present: ", file.exists(file.path(TABLE_DIR, "spectral_axis_power2_shannon_phylogenetic_summary.csv"))),
    "",
    "## Result",
    "",
    "The side-investigation scatter figures and clean relationship summary CSV were produced."
  )
  writeLines(validation_lines, VALIDATION_REPORT)
}

run_spectral_axis_power2_side_investigation <- function() {
  data <- load_analysis_data()
  relationships <- summarize_relationships(data)
  summary_file <- file.path(TABLE_DIR, "spectral_axis_power2_shannon_phylogenetic_summary.csv")
  utils::write.csv(relationships, summary_file, row.names = FALSE)
  figure_files <- save_spectral_metric_figures(data, relationships)
  write_reports(data, relationships, figure_files)

  stray_plot <- file.path(PROJECT_DIR, "Rplots.pdf")
  if (file.exists(stray_plot)) {
    unlink(stray_plot)
  }

  message("Spectral-axis power-2 Shannon/phylogenetic side investigation complete.")
  invisible(list(
    data = data,
    relationships = relationships,
    summary_file = summary_file,
    figure_files = figure_files
  ))
}

if (sys.nframe() == 0) {
  run_spectral_axis_power2_side_investigation()
}
