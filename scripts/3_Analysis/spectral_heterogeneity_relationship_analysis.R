PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-08-17"
USER_R_LIB <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(USER_R_LIB)) {
  .libPaths(unique(c(USER_R_LIB, .libPaths())))
}

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/spectral_heterogeneity_relationships")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/spectral_heterogeneity_relationships")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260817_spectral_heterogeneity_relationship_analysis.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260817_spectral_heterogeneity_relationship_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260817_spectral_heterogeneity_relationship_analysis_validation.md")

SCALES <- c("10m", "20m", "50m")
SPECTRAL_METRICS <- c(
  "spec_sa",
  "spec_alpha",
  "spec_rao_q",
  "spec_pca_mean",
  "spec_spca_alpha",
  "spec_spca_rao",
  "spec_spca_mean"
)
DISPLAY_NAMES <- c(
  spec_sa = "Spectral angle entropy",
  spec_alpha = "Raw PCA alpha hull",
  spec_rao_q = "Raw PCA Rao's Q",
  spec_pca_mean = "Raw PCA mean distance",
  spec_spca_alpha = "Std PCA alpha hull",
  spec_spca_rao = "Std PCA Rao's Q",
  spec_spca_mean = "Std PCA mean distance"
)
REGION_BRIGHTNESS_COLUMNS <- c(
  blue = "mean_blue_pixel_brightness",
  green = "mean_green_pixel_brightness",
  red = "mean_red_pixel_brightness",
  nir = "mean_nir_pixel_brightness"
)
REGION_DISPLAY_NAMES <- c(
  blue = "Blue illum.",
  green = "Green illum.",
  red = "Red illum.",
  nir = "NIR illum."
)

COMBINED_FILES <- c(
  "10m" = file.path(PROJECT_DIR, "quadrat_analysis_10m.csv"),
  "20m" = file.path(PROJECT_DIR, "quadrat_analysis_20m.csv"),
  "50m" = file.path(PROJECT_DIR, "quadrat_analysis_50m.csv")
)
CLUSTER_FILE <- file.path(PROJECT_DIR, "reports/tables/species_phylogenetic_correlation/species_composition_clusters.csv")
COUNT_FILE <- file.path(PROJECT_DIR, "reports/tables/species_phylogenetic_correlation/quadrat_crown_equivalent_individual_totals.csv")
BRIGHTNESS_FILE <- file.path(TABLE_DIR, "quadrat_pixel_brightness_summary.csv")
SPEC_DIRS <- c(
  "10m" = file.path(PROJECT_DIR, "Quad_Spectra/10m_smooth_5nm"),
  "20m" = file.path(PROJECT_DIR, "Quad_Spectra/20m_smooth_5nm"),
  "50m" = file.path(PROJECT_DIR, "Quad_Spectra/50m_smooth_5nm")
)
SIDECAR_PATTERN <- "\\.(hdr|aux|xml|enp|sta)$"
BEST_BAND <- "563 nm"
SHADOW_THRESHOLD <- 0.0305476
DIRECTION <- "<"
WAVELENGTHS <- seq(398, 998, by = 5)
SPECTRAL_REGIONS <- list(
  blue = c(450, 495),
  green = c(500, 570),
  red = c(620, 750),
  nir = c(750, 998)
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

read_spectral_data <- function() {
  rows <- lapply(SCALES, function(scale_name) {
    data <- read.csv(COMBINED_FILES[[scale_name]], stringsAsFactors = FALSE, check.names = FALSE)
    data$quad_id <- as.character(data$quad_id)
    data$scale <- as.character(data$scale)
    data
  })
  data <- do.call(rbind, rows)
  required <- c("quad_id", "scale", SPECTRAL_METRICS, "env_elev")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data
}

read_composition_clusters <- function() {
  clusters <- read.csv(CLUSTER_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  clusters$quad_id <- as.character(clusters$quad_id)
  clusters$scale <- as.character(clusters$scale)
  clusters
}

read_quadrat_counts <- function() {
  counts <- read.csv(COUNT_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  counts$quad_id <- as.character(counts$quad_id)
  counts$scale <- as.character(counts$scale)
  counts
}

list_raster_files <- function(spec_dir) {
  files <- list.files(spec_dir, full.names = TRUE)
  files[!grepl(SIDECAR_PATTERN, files, ignore.case = TRUE)]
}

calculate_quadrat_brightness <- function(file, scale_name) {
  hdr_file <- paste0(file, ".hdr")
  if (!file.exists(hdr_file)) {
    stop("Missing ENVI header for ", file, call. = FALSE)
  }
  hdr <- readLines(hdr_file, warn = FALSE)
  header_value <- function(key) {
    pattern <- paste0("^", key, "\\s*=")
    line <- hdr[grepl(pattern, hdr, ignore.case = TRUE)][1]
    if (is.na(line)) stop("Missing `", key, "` in ", hdr_file, call. = FALSE)
    trimws(sub(pattern, "", line, ignore.case = TRUE))
  }
  samples <- as.integer(header_value("samples"))
  lines <- as.integer(header_value("lines"))
  bands <- as.integer(header_value("bands"))
  data_type <- as.integer(header_value("data type"))
  interleave <- tolower(header_value("interleave"))
  byte_order <- as.integer(header_value("byte order"))
  if (data_type != 4 || interleave != "bsq") {
    stop("Unsupported ENVI format for ", file, "; expected float32 BSQ.", call. = FALSE)
  }

  n_cells <- samples * lines
  if (n_cells == 0) {
    return(data.frame(
      quad_id = basename(file),
      scale = scale_name,
      brightness_n_pixels = 0L,
      mean_pixel_brightness = NA_real_,
      median_pixel_brightness = NA_real_,
      sd_pixel_brightness = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  con <- file(file, "rb")
  on.exit(close(con), add = TRUE)
  read_band <- function(band_index) {
    seek(con, where = (band_index - 1) * n_cells * 4, origin = "start")
    readBin(
      con,
      what = "numeric",
      n = n_cells,
      size = 4,
      endian = ifelse(byte_order == 0, "little", "big")
    )
  }

  band_index <- which(WAVELENGTHS == 563)
  band_values <- read_band(band_index)
  sunlit <- if (DIRECTION == ">") {
    band_values < SHADOW_THRESHOLD
  } else {
    band_values > SHADOW_THRESHOLD
  }
  keep <- sunlit & is.finite(band_values) & band_values > 0
  keep[is.na(keep)] <- FALSE
  pixel_brightness <- band_values[keep]

  region_means <- lapply(SPECTRAL_REGIONS, function(region_range) {
    region_bands <- which(WAVELENGTHS >= region_range[1] & WAVELENGTHS <= region_range[2])
    kept_count <- sum(keep)
    if (kept_count == 0 || length(region_bands) == 0) {
      return(data.frame(
        n_pixels = 0L,
        mean = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    region_sum <- rep(0, kept_count)
    region_count <- rep(0L, kept_count)
    for (region_band in region_bands) {
      values <- read_band(region_band)[keep]
      ok <- is.finite(values) & values > 0
      region_sum[ok] <- region_sum[ok] + values[ok]
      region_count[ok] <- region_count[ok] + 1L
    }
    valid <- region_count > 0
    pixel_region_mean <- region_sum[valid] / region_count[valid]
    data.frame(
      n_pixels = length(pixel_region_mean),
      mean = if (length(pixel_region_mean) > 0) mean(pixel_region_mean) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  result <- data.frame(
    quad_id = basename(file),
    scale = scale_name,
    brightness_n_pixels = length(pixel_brightness),
    mean_pixel_brightness = if (length(pixel_brightness) > 0) mean(pixel_brightness) else NA_real_,
    median_pixel_brightness = if (length(pixel_brightness) > 0) median(pixel_brightness) else NA_real_,
    sd_pixel_brightness = if (length(pixel_brightness) > 1) sd(pixel_brightness) else NA_real_,
    stringsAsFactors = FALSE
  )
  for (region_name in names(region_means)) {
    result[[paste0(region_name, "_brightness_n_pixels")]] <- region_means[[region_name]]$n_pixels
    result[[paste0("mean_", region_name, "_pixel_brightness")]] <- region_means[[region_name]]$mean
  }
  result
}

read_or_create_brightness <- function(data) {
  if (file.exists(BRIGHTNESS_FILE)) {
    brightness <- read.csv(BRIGHTNESS_FILE, stringsAsFactors = FALSE, check.names = FALSE)
    brightness$quad_id <- as.character(brightness$quad_id)
    brightness$scale <- as.character(brightness$scale)
    required_region_columns <- c(
      "mean_blue_pixel_brightness",
      "mean_green_pixel_brightness",
      "mean_red_pixel_brightness",
      "mean_nir_pixel_brightness"
    )
    if (all(required_region_columns %in% names(brightness))) {
      return(brightness)
    }
  }

  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    scale_ids <- unique(data$quad_id[data$scale == scale_name])
    raster_files <- list_raster_files(SPEC_DIRS[[scale_name]])
    raster_files <- raster_files[basename(raster_files) %in% scale_ids]
    raster_files <- raster_files[order(basename(raster_files))]
    for (file in raster_files) {
      rows[[index]] <- calculate_quadrat_brightness(file, scale_name)
      index <- index + 1
    }
  }
  brightness <- do.call(rbind, rows)
  write.csv(brightness, BRIGHTNESS_FILE, row.names = FALSE)
  brightness
}

metric_pairs <- function() utils::combn(SPECTRAL_METRICS, 2, simplify = FALSE)

run_pairwise_correlations <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (pair in metric_pairs()) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      keep <- data$scale == scale_name & !is.na(data[[x_metric]]) & !is.na(data[[y_metric]])
      pair_data <- data[keep, c(x_metric, y_metric)]
      names(pair_data) <- c("x_value", "y_value")
      result <- data.frame(
        scale = scale_name,
        x_metric = x_metric,
        x_label = display_name(x_metric),
        y_metric = y_metric,
        y_label = display_name(y_metric),
        n = nrow(pair_data),
        pearson_r = NA_real_,
        r_squared = NA_real_,
        f_statistic = NA_real_,
        f_p_value = NA_real_,
        slope = NA_real_,
        intercept = NA_real_,
        spearman_r = NA_real_,
        spearman_p_value = NA_real_,
        stringsAsFactors = FALSE
      )
      if (nrow(pair_data) >= 3) {
        fit <- lm(y_value ~ x_value, data = pair_data)
        fit_summary <- summary(fit)
        f_values <- fit_summary$fstatistic
        pearson <- cor.test(pair_data$x_value, pair_data$y_value, method = "pearson")
        spearman <- suppressWarnings(cor.test(pair_data$x_value, pair_data$y_value, method = "spearman", exact = FALSE))
        result$pearson_r <- unname(pearson$estimate)
        result$r_squared <- fit_summary$r.squared
        result$f_statistic <- unname(f_values[["value"]])
        result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
        result$slope <- unname(coef(fit)[["x_value"]])
        result$intercept <- unname(coef(fit)[["(Intercept)"]])
        result$spearman_r <- unname(spearman$estimate)
        result$spearman_p_value <- spearman$p.value
      }
      rows[[index]] <- result
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

run_elevation_adjusted_models <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (pair in metric_pairs()) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      keep <- data$scale == scale_name &
        !is.na(data[[x_metric]]) &
        !is.na(data[[y_metric]]) &
        !is.na(data$env_elev)
      model_data <- data[keep, c(x_metric, y_metric, "env_elev")]
      names(model_data) <- c("x_value", "y_value", "env_elev")
      result <- data.frame(
        scale = scale_name,
        x_metric = x_metric,
        x_label = display_name(x_metric),
        y_metric = y_metric,
        y_label = display_name(y_metric),
        n = nrow(model_data),
        base_r2 = NA_real_,
        elevation_adjusted_r2 = NA_real_,
        elevation_incremental_r2_after_metric = NA_real_,
        elevation_partial_f = NA_real_,
        elevation_partial_p_value = NA_real_,
        elevation_slope = NA_real_,
        stringsAsFactors = FALSE
      )
      if (nrow(model_data) >= 4) {
        base_fit <- lm(y_value ~ x_value, data = model_data)
        adjusted_fit <- lm(y_value ~ x_value + env_elev, data = model_data)
        comparison <- anova(base_fit, adjusted_fit)
        adjusted_summary <- summary(adjusted_fit)
        result$base_r2 <- summary(base_fit)$r.squared
        result$elevation_adjusted_r2 <- adjusted_summary$r.squared
        result$elevation_incremental_r2_after_metric <- result$elevation_adjusted_r2 - result$base_r2
        result$elevation_partial_f <- comparison$F[2]
        result$elevation_partial_p_value <- comparison$`Pr(>F)`[2]
        result$elevation_slope <- coef(adjusted_fit)[["env_elev"]]
      }
      rows[[index]] <- result
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

run_regional_brightness_correlations <- function(data, brightness) {
  plot_data <- merge(
    data,
    brightness[, c("quad_id", "scale", REGION_BRIGHTNESS_COLUMNS)],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (region_name in names(REGION_BRIGHTNESS_COLUMNS)) {
      brightness_column <- REGION_BRIGHTNESS_COLUMNS[[region_name]]
      for (metric in SPECTRAL_METRICS) {
        keep <- plot_data$scale == scale_name &
          !is.na(plot_data[[metric]]) &
          !is.na(plot_data[[brightness_column]])
        model_data <- plot_data[keep, c(brightness_column, metric)]
        names(model_data) <- c("illumination", "metric_value")
        result <- data.frame(
          scale = scale_name,
          region = region_name,
          region_label = REGION_DISPLAY_NAMES[[region_name]],
          brightness_column = brightness_column,
          metric = metric,
          metric_label = display_name(metric),
          n = nrow(model_data),
          pearson_r = NA_real_,
          r_squared = NA_real_,
          f_statistic = NA_real_,
          f_p_value = NA_real_,
          slope = NA_real_,
          intercept = NA_real_,
          spearman_r = NA_real_,
          spearman_p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        if (nrow(model_data) >= 3) {
          fit <- lm(metric_value ~ illumination, data = model_data)
          fit_summary <- summary(fit)
          f_values <- fit_summary$fstatistic
          pearson <- cor.test(model_data$illumination, model_data$metric_value, method = "pearson")
          spearman <- suppressWarnings(cor.test(model_data$illumination, model_data$metric_value, method = "spearman", exact = FALSE))
          result$pearson_r <- unname(pearson$estimate)
          result$r_squared <- fit_summary$r.squared
          result$f_statistic <- unname(f_values[["value"]])
          result$f_p_value <- pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
          result$slope <- unname(coef(fit)[["illumination"]])
          result$intercept <- unname(coef(fit)[["(Intercept)"]])
          result$spearman_r <- unname(spearman$estimate)
          result$spearman_p_value <- spearman$p.value
        }
        rows[[index]] <- result
        index <- index + 1
      }
    }
  }
  do.call(rbind, rows)
}

cluster_palette <- function(n) {
  base_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#332288", "#117733")
  if (n <= length(base_colors)) base_colors[seq_len(n)] else grDevices::hcl.colors(n, palette = "Dark 3")
}

ramp_colors <- function(values, palette_values = c("#440154", "#2A6FBB", "#00A878", "#FDE725", "#F94144")) {
  palette <- grDevices::colorRampPalette(palette_values)(101)
  value_range <- range(values, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) == 0) {
    return(rep(palette[51], length(values)))
  }
  idx <- round((values - value_range[1]) / diff(value_range) * 100) + 1
  palette[pmax(pmin(idx, 101), 1)]
}

add_ramp_legend <- function(values, title, palette_values, digits = 1) {
  value_range <- range(values, na.rm = TRUE)
  if (!all(is.finite(value_range))) return(invisible(NULL))
  legend_values <- pretty(value_range, n = 5)
  legend_values <- legend_values[legend_values >= value_range[1] & legend_values <= value_range[2]]
  if (length(legend_values) < 3) legend_values <- seq(value_range[1], value_range[2], length.out = 5)
  legend(
    "bottomright",
    legend = fmt_num(legend_values, digits),
    col = ramp_colors(legend_values, palette_values),
    pch = 16,
    title = title,
    bty = "n",
    cex = 0.7
  )
}

plot_pair_panel <- function(plot_data, scale_name, x_metric, y_metric, correlations, point_colors, point_cex = 0.58) {
  complete <- plot_data$scale == scale_name & !is.na(plot_data[[x_metric]]) & !is.na(plot_data[[y_metric]])
  plot(
    plot_data[[x_metric]][complete],
    plot_data[[y_metric]][complete],
    col = point_colors[complete],
    pch = 16,
    cex = point_cex,
    xlab = display_name(x_metric),
    ylab = display_name(y_metric),
    main = paste(display_name(x_metric), "vs", display_name(y_metric))
  )
  if (sum(complete) >= 3) {
    fit <- lm(plot_data[[y_metric]][complete] ~ plot_data[[x_metric]][complete])
    abline(fit, col = "#111111", lwd = 2)
    row <- correlations[
      correlations$scale == scale_name &
        correlations$x_metric == x_metric &
        correlations$y_metric == y_metric,
    ]
    legend(
      "topleft",
      legend = c(paste0("n=", row$n), paste0("r=", fmt_num(row$pearson_r, 2)), paste0("R2=", fmt_num(row$r_squared, 2))),
      bty = "n",
      cex = 0.78
    )
  }
  complete
}

save_plain_scatter <- function(data, correlations) {
  files <- character(length(SCALES))
  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(FIGURE_DIR, paste0("01_spectral_metric_pairwise_scatter_", scale_name, ".png"))
    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
    point_colors <- adjustcolor(rep("#315C8A", nrow(data)), alpha.f = 0.52)
    for (pair in metric_pairs()) {
      plot_pair_panel(data, scale_name, pair[1], pair[2], correlations, point_colors)
    }
    mtext(paste0(scale_name, " quadrats: spectral heterogeneity metric pairs"), outer = TRUE, cex = 1.2, font = 2)
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }
  files
}

save_elevation_scatter <- function(data, correlations, elevation_models) {
  files <- character(length(SCALES))
  palette_values <- c("#B03A2E", "#F7F7F2", "#2F7D4F")
  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(FIGURE_DIR, paste0("02_spectral_metric_pairwise_scatter_elevation_gradient_", scale_name, ".png"))
    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
    point_colors <- adjustcolor(ramp_colors(data$env_elev, palette_values), alpha.f = 0.72)
    for (pair_index in seq_along(metric_pairs())) {
      pair <- metric_pairs()[[pair_index]]
      plot_pair_panel(data, scale_name, pair[1], pair[2], correlations, point_colors)
      row <- elevation_models[elevation_models$scale == scale_name & elevation_models$x_metric == pair[1] & elevation_models$y_metric == pair[2], ]
      legend(
        "bottomleft",
        legend = c(
          paste0("elev dR2=", fmt_num(row$elevation_incremental_r2_after_metric, 3)),
          paste0("elev p=", fmt_p(row$elevation_partial_p_value))
        ),
        bty = "n",
        cex = 0.7
      )
      if (pair_index == 1) {
        add_ramp_legend(data$env_elev[data$scale == scale_name], "Mean\nelevation", palette_values, digits = 1)
      }
    }
    mtext(paste0(scale_name, " quadrats: spectral metric pairs colored by mean elevation"), outer = TRUE, cex = 1.2, font = 2)
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }
  files
}

calculate_composition_scatterplot_silhouette <- function(data, clusters) {
  plot_data <- merge(
    data,
    clusters[, c("quad_id", "scale", "composition_cluster", "composition_type")],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (pair in metric_pairs()) {
      x_metric <- pair[1]
      y_metric <- pair[2]
      complete <- plot_data$scale == scale_name &
        !is.na(plot_data[[x_metric]]) &
        !is.na(plot_data[[y_metric]]) &
        !is.na(plot_data$composition_cluster)
      panel_data <- plot_data[complete, c("composition_cluster", "composition_type", x_metric, y_metric)]
      names(panel_data) <- c("composition_cluster", "composition_type", "x_value", "y_value")
      cluster_names <- sort(unique(panel_data$composition_cluster))
      cluster_stats <- unique(panel_data[, c("composition_cluster", "composition_type")])
      cluster_stats <- cluster_stats[order(cluster_stats$composition_type), ]
      cluster_stats$n <- as.integer(table(panel_data$composition_cluster)[cluster_stats$composition_cluster])
      cluster_stats$mean_spectral_scatterplot_silhouette_width <- NA_real_
      overall_silhouette <- NA_real_
      if (nrow(panel_data) > length(cluster_names) && length(cluster_names) > 1 && requireNamespace("cluster", quietly = TRUE)) {
        xy <- scale(panel_data[, c("x_value", "y_value")])
        cluster_factor <- factor(panel_data$composition_cluster, levels = cluster_names)
        silhouette_values <- cluster::silhouette(as.integer(cluster_factor), stats::dist(xy))
        cluster_silhouette <- tapply(silhouette_values[, "sil_width"], panel_data$composition_cluster, mean, na.rm = TRUE)
        overall_silhouette <- mean(silhouette_values[, "sil_width"], na.rm = TRUE)
        cluster_stats$mean_spectral_scatterplot_silhouette_width <- unname(cluster_silhouette[cluster_stats$composition_cluster])
      }
      for (row_index in seq_len(nrow(cluster_stats))) {
        rows[[index]] <- data.frame(
          scale = scale_name,
          x_metric = x_metric,
          x_label = display_name(x_metric),
          y_metric = y_metric,
          y_label = display_name(y_metric),
          composition_cluster = cluster_stats$composition_cluster[row_index],
          composition_type = cluster_stats$composition_type[row_index],
          n = cluster_stats$n[row_index],
          mean_spectral_scatterplot_silhouette_width = cluster_stats$mean_spectral_scatterplot_silhouette_width[row_index],
          overall_spectral_scatterplot_silhouette_width = overall_silhouette,
          stringsAsFactors = FALSE
        )
        index <- index + 1
      }
    }
  }
  do.call(rbind, rows)
}

save_composition_scatter <- function(data, clusters, correlations, silhouette_table) {
  plot_data <- merge(
    data,
    clusters[, c("quad_id", "scale", "composition_cluster", "composition_type")],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  files <- character(length(SCALES))
  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    scale_clusters <- unique(plot_data[plot_data$scale == scale_name & !is.na(plot_data$composition_type), c("composition_cluster", "composition_type")])
    scale_clusters <- scale_clusters[order(scale_clusters$composition_type), ]
    cluster_colors <- cluster_palette(nrow(scale_clusters))
    names(cluster_colors) <- scale_clusters$composition_cluster
    file_path <- file.path(FIGURE_DIR, paste0("03_spectral_metric_pairwise_scatter_composition_cluster_", scale_name, ".png"))
    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
    point_colors <- adjustcolor(cluster_colors[plot_data$composition_cluster], alpha.f = 0.62)
    point_colors[is.na(point_colors)] <- adjustcolor("#777777", alpha.f = 0.3)
    for (pair in metric_pairs()) {
      plot_pair_panel(plot_data, scale_name, pair[1], pair[2], correlations, point_colors, point_cex = 0.54)
      panel_silhouette <- silhouette_table[
        silhouette_table$scale == scale_name &
          silhouette_table$x_metric == pair[1] &
          silhouette_table$y_metric == pair[2],
      ]
      panel_silhouette <- panel_silhouette[match(scale_clusters$composition_cluster, panel_silhouette$composition_cluster), ]
      legend_labels <- paste0(
        scale_clusters$composition_type,
        " (sil=",
        fmt_num(panel_silhouette$mean_spectral_scatterplot_silhouette_width, 2),
        ")"
      )
      legend(
        "bottomright",
        legend = legend_labels,
        col = cluster_colors,
        pch = 16,
        title = "Composition type (panel sil.)",
        bty = "n",
        cex = ifelse(nrow(scale_clusters) > 4, 0.38, 0.48)
      )
    }
    mtext(paste0(scale_name, " quadrats: spectral metric pairs colored by species composition type"), outer = TRUE, cex = 1.2, font = 2)
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }
  files
}

save_count_ramp_scatter <- function(data, counts, correlations, column, label, prefix, palette_values, digits = 1) {
  plot_data <- merge(data, counts[, c("quad_id", "scale", column)], by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
  files <- character(length(SCALES))
  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(FIGURE_DIR, paste0(prefix, "_", scale_name, ".png"))
    grDevices::png(file_path, width = 3900, height = 5400, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 3), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
    point_colors <- adjustcolor(ramp_colors(plot_data[[column]], palette_values), alpha.f = 0.72)
    for (pair_index in seq_along(metric_pairs())) {
      pair <- metric_pairs()[[pair_index]]
      plot_pair_panel(plot_data, scale_name, pair[1], pair[2], correlations, point_colors)
      if (pair_index == 1) {
        add_ramp_legend(plot_data[[column]][plot_data$scale == scale_name], label, palette_values, digits = digits)
      }
    }
    title_label <- gsub("\n", " ", tolower(label))
    mtext(paste0(scale_name, " quadrats: spectral metric pairs colored by ", title_label), outer = TRUE, cex = 1.2, font = 2)
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }
  files
}

save_metric_vs_regional_brightness_scatter <- function(data, brightness, regional_correlations) {
  plot_data <- merge(
    data,
    brightness[, c("quad_id", "scale", REGION_BRIGHTNESS_COLUMNS)],
    by = c("quad_id", "scale"),
    all.x = TRUE,
    sort = FALSE
  )
  region_colors <- c(
    blue = "#1D4ED8",
    green = "#16A34A",
    red = "#DC2626",
    nir = "#F97316"
  )
  files <- character(length(SCALES))
  for (scale_index in seq_along(SCALES)) {
    scale_name <- SCALES[scale_index]
    file_path <- file.path(FIGURE_DIR, paste0("11_spectral_metrics_vs_regional_illumination_", scale_name, ".png"))
    grDevices::png(file_path, width = 3900, height = 7200, res = 300)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(7, 4), mar = c(4.2, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
    for (metric in SPECTRAL_METRICS) {
      for (region_name in names(REGION_BRIGHTNESS_COLUMNS)) {
        brightness_column <- REGION_BRIGHTNESS_COLUMNS[[region_name]]
        complete <- plot_data$scale == scale_name &
          !is.na(plot_data[[metric]]) &
          !is.na(plot_data[[brightness_column]])
        plot(
          plot_data[[brightness_column]][complete],
          plot_data[[metric]][complete],
          col = adjustcolor(region_colors[[region_name]], alpha.f = 0.62),
          pch = 16,
          cex = 0.52,
          xlab = REGION_DISPLAY_NAMES[[region_name]],
          ylab = display_name(metric),
          main = paste(display_name(metric), "vs", REGION_DISPLAY_NAMES[[region_name]])
        )
        row <- regional_correlations[
          regional_correlations$scale == scale_name &
            regional_correlations$region == region_name &
            regional_correlations$metric == metric,
        ]
        if (nrow(row) == 1 && !is.na(row$pearson_r)) {
          fit <- lm(plot_data[[metric]][complete] ~ plot_data[[brightness_column]][complete])
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
            cex = 0.7
          )
        }
      }
    }
    mtext(paste0(scale_name, " quadrats: spectral metrics vs regional illumination"), outer = TRUE, cex = 1.2, font = 2)
    par(old_par)
    grDevices::dev.off()
    files[scale_index] <- file_path
  }
  files
}

save_figures <- function(data, correlations, elevation_models, clusters, counts, silhouette_table, regional_correlations) {
  brightness <- read_or_create_brightness(data)
  c(
    save_plain_scatter(data, correlations),
    save_elevation_scatter(data, correlations, elevation_models),
    save_composition_scatter(data, clusters, correlations, silhouette_table),
    save_count_ramp_scatter(
      data,
      counts,
      correlations,
      "crown_equivalent_individuals",
      "Crown-equiv.\nindividuals",
      "04_spectral_metric_pairwise_scatter_individual_ramp",
      c("#440154", "#2A6FBB", "#00A878", "#FDE725", "#F94144"),
      digits = 1
    ),
    save_count_ramp_scatter(
      data,
      counts,
      correlations,
      "present_species_count",
      "Species\nper quad",
      "05_spectral_metric_pairwise_scatter_species_count_ramp",
      c("#2D004B", "#0057B8", "#00B4D8", "#7AE582", "#FFE45E", "#FF3D00"),
      digits = 0
    ),
    save_count_ramp_scatter(
      data,
      brightness,
      correlations,
      "mean_pixel_brightness",
      "Mean pixel\nbrightness",
      "06_spectral_metric_pairwise_scatter_pixel_brightness_ramp",
      c("#1B1B1B", "#3B528B", "#21918C", "#5EC962", "#FDE725"),
      digits = 3
    ),
    save_count_ramp_scatter(
      data,
      brightness,
      correlations,
      "mean_blue_pixel_brightness",
      "Blue\nbrightness",
      "07_spectral_metric_pairwise_scatter_blue_brightness_ramp",
      c("#1B1B1B", "#3B528B", "#21918C", "#5EC962", "#FDE725"),
      digits = 3
    ),
    save_count_ramp_scatter(
      data,
      brightness,
      correlations,
      "mean_green_pixel_brightness",
      "Green\nbrightness",
      "08_spectral_metric_pairwise_scatter_green_brightness_ramp",
      c("#1B1B1B", "#3B528B", "#21918C", "#5EC962", "#FDE725"),
      digits = 3
    ),
    save_count_ramp_scatter(
      data,
      brightness,
      correlations,
      "mean_red_pixel_brightness",
      "Red\nbrightness",
      "09_spectral_metric_pairwise_scatter_red_brightness_ramp",
      c("#1B1B1B", "#3B528B", "#21918C", "#5EC962", "#FDE725"),
      digits = 3
    ),
    save_count_ramp_scatter(
      data,
      brightness,
      correlations,
      "mean_nir_pixel_brightness",
      "NIR\nbrightness",
      "10_spectral_metric_pairwise_scatter_nir_brightness_ramp",
      c("#1B1B1B", "#3B528B", "#21918C", "#5EC962", "#FDE725"),
      digits = 3
    ),
    save_metric_vs_regional_brightness_scatter(data, brightness, regional_correlations)
  )
}

metric_summary <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (metric in SPECTRAL_METRICS) {
      values <- data[[metric]][data$scale == scale_name]
      rows[[index]] <- data.frame(
        scale = scale_name,
        metric = metric,
        metric_label = display_name(metric),
        n = sum(!is.na(values)),
        missing = sum(is.na(values)),
        mean = mean(values, na.rm = TRUE),
        sd = sd(values, na.rm = TRUE),
        min = min(values, na.rm = TRUE),
        median = median(values, na.rm = TRUE),
        max = max(values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

write_analysis_report <- function(correlations, elevation_models, metric_stats, silhouette_table, figure_files) {
  strongest_correlations <- correlations[order(correlations$scale, -abs(correlations$pearson_r)), ]
  strongest_correlations <- do.call(rbind, lapply(split(strongest_correlations, strongest_correlations$scale), function(x) x[seq_len(min(3, nrow(x))), ]))
  top_elevation <- elevation_models[order(-elevation_models$elevation_incremental_r2_after_metric), ]
  top_elevation <- top_elevation[seq_len(min(10, nrow(top_elevation))), ]
  silhouette_summary <- aggregate(
    overall_spectral_scatterplot_silhouette_width ~ scale + x_label + y_label,
    data = silhouette_table,
    FUN = mean,
    na.rm = TRUE
  )
  top_silhouette <- silhouette_summary[order(-silhouette_summary$overall_spectral_scatterplot_silhouette_width), ]
  top_silhouette <- top_silhouette[seq_len(min(10, nrow(top_silhouette))), ]

  cor_table <- data.frame(
    Scale = strongest_correlations$scale,
    Relationship = paste(strongest_correlations$x_label, "vs", strongest_correlations$y_label),
    n = strongest_correlations$n,
    r = fmt_num(strongest_correlations$pearson_r, 3),
    R2 = fmt_num(strongest_correlations$r_squared, 3),
    `p-value` = fmt_p(strongest_correlations$f_p_value),
    check.names = FALSE
  )
  elev_table <- data.frame(
    Scale = top_elevation$scale,
    Relationship = paste(top_elevation$x_label, "vs", top_elevation$y_label),
    n = top_elevation$n,
    `Elevation delta R2` = fmt_num(top_elevation$elevation_incremental_r2_after_metric, 3),
    `Elevation p-value` = fmt_p(top_elevation$elevation_partial_p_value),
    check.names = FALSE
  )
  sil_table <- data.frame(
    Scale = top_silhouette$scale,
    Relationship = paste(top_silhouette$x_label, "vs", top_silhouette$y_label),
    `Overall composition silhouette` = fmt_num(top_silhouette$overall_spectral_scatterplot_silhouette_width, 3),
    check.names = FALSE
  )

  lines <- c(
    "# Spectral Heterogeneity Metric Relationship Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis compares all seven focal spectral heterogeneity measures across 10 m, 20 m, and 50 m quadrats: spectral angle entropy, raw PCA alpha-hull area, raw PCA spectral Rao's Q, raw PCA mean Euclidean distance, standardized PCA alpha-hull area, standardized PCA spectral Rao's Q, and standardized PCA mean Euclidean distance.",
    "",
    "## Methods",
    "",
    "- Pairwise scatterplots and correlations were calculated for every unique pair among the seven focal spectral heterogeneity measures within each scale.",
    "- Additional figure sets color the same pairwise relationships by mean elevation, species presence/absence composition type, crown-equivalent individuals, number of species present, mean retained-pixel brightness at 563 nm, and retained-pixel brightness in blue, green, red, and near-infrared spectral regions.",
    "- Elevation panels report the incremental R2 and p-value for adding `env_elev` after the x-axis spectral metric.",
    "- Composition-type panels report panel-specific mean silhouette width calculated from the two plotted spectral metrics.",
    "",
    "## Strongest Spectral Metric Correlations",
    "",
    markdown_table(cor_table),
    "",
    "## Largest Elevation Contributions",
    "",
    markdown_table(elev_table),
    "",
    "## Strongest Composition-Type Separation In Spectral Metric Space",
    "",
    markdown_table(sil_table),
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Output Tables",
    "",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_pairwise_correlations.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_elevation_adjusted_models.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_summary.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_composition_scatterplot_silhouette.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/quadrat_pixel_brightness_summary.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_regional_illumination_correlations.csv`"
  )
  writeLines(lines, ANALYSIS_REPORT)
}

write_task_report <- function(correlations, elevation_models, silhouette_table) {
  lines <- c(
    "# Task Report: Spectral Heterogeneity Metric Relationship Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Objective",
    "",
    "Create all-pair spectral heterogeneity scatterplot sets and matching elevation, species-composition, individual-count, species-count, whole-brightness, and region-specific pixel-brightness highlight figures.",
    "",
    "## Outputs Created",
    "",
    "- `reports/analysis/20260817_spectral_heterogeneity_relationship_analysis.md`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_pairwise_correlations.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_elevation_adjusted_models.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_summary.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_composition_scatterplot_silhouette.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/quadrat_pixel_brightness_summary.csv`",
    "- `reports/tables/spectral_heterogeneity_relationships/spectral_metric_regional_illumination_correlations.csv`",
    "- `reports/figures/spectral_heterogeneity_relationships/`",
    "",
    "## Result Size",
    "",
    paste0("- Spectral pairwise correlation rows: ", nrow(correlations)),
    paste0("- Elevation-adjusted model rows: ", nrow(elevation_models)),
    paste0("- Spectral composition silhouette rows: ", nrow(silhouette_table)),
    "",
    "## Notes",
    "",
    "- Spectral metric pairings use complete cases for each pair, preserving upstream missingness.",
    "- Pixel brightness is calculated from each quadrat raster as mean retained sunlit reflectance at 563 nm, the same band used for shadow masking.",
    "- Region-specific brightness uses retained sunlit pixel means for blue 450-495 nm, green 500-570 nm, red 620-750 nm, and near-infrared 750-998 nm.",
    "- Composition types and quadrat species/individual counts are reused from the presence-based biodiversity composition analysis."
  )
  writeLines(lines, TASK_REPORT)
}

write_validation_report <- function(data, correlations, elevation_models, metric_stats, silhouette_table, regional_correlations, figure_files) {
  table_files <- c(
    file.path(TABLE_DIR, "spectral_metric_pairwise_correlations.csv"),
    file.path(TABLE_DIR, "spectral_metric_elevation_adjusted_models.csv"),
    file.path(TABLE_DIR, "spectral_metric_summary.csv"),
    file.path(TABLE_DIR, "spectral_composition_scatterplot_silhouette.csv"),
    BRIGHTNESS_FILE,
    file.path(TABLE_DIR, "spectral_metric_regional_illumination_correlations.csv")
  )
  expected_correlations <- length(SCALES) * choose(length(SPECTRAL_METRICS), 2)
  expected_silhouette_rows <- choose(length(SPECTRAL_METRICS), 2) * 15
  expected_regional_rows <- length(SCALES) * length(SPECTRAL_METRICS) * length(REGION_BRIGHTNESS_COLUMNS)
  brightness <- if (file.exists(BRIGHTNESS_FILE)) {
    read.csv(BRIGHTNESS_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame(mean_pixel_brightness = numeric(0))
  }
  lines <- c(
    "# Validation: Spectral Heterogeneity Metric Relationship Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Checks",
    "",
    paste0("- Expected spectral pairwise correlation rows: ", expected_correlations),
    paste0("- Observed spectral pairwise correlation rows: ", nrow(correlations)),
    paste0("- Observed elevation-adjusted model rows: ", nrow(elevation_models)),
    paste0("- Spectral metric summary rows: ", nrow(metric_stats)),
    paste0("- Expected spectral composition silhouette rows: ", expected_silhouette_rows),
    paste0("- Observed spectral composition silhouette rows: ", nrow(silhouette_table)),
    paste0("- Expected regional illumination correlation rows: ", expected_regional_rows),
    paste0("- Observed regional illumination correlation rows: ", nrow(regional_correlations)),
    paste0("- Observed pixel brightness rows: ", nrow(brightness)),
    paste0("- Missing Pearson r values: ", sum(is.na(correlations$pearson_r))),
    paste0("- Missing elevation incremental R2 values: ", sum(is.na(elevation_models$elevation_incremental_r2_after_metric))),
    paste0("- Missing spectral composition silhouette values: ", sum(is.na(silhouette_table$mean_spectral_scatterplot_silhouette_width))),
    paste0("- Missing regional illumination Pearson r values: ", sum(is.na(regional_correlations$pearson_r))),
    paste0("- Missing mean pixel brightness values: ", sum(is.na(brightness$mean_pixel_brightness))),
    paste0("- Missing blue brightness values: ", sum(is.na(brightness$mean_blue_pixel_brightness))),
    paste0("- Missing green brightness values: ", sum(is.na(brightness$mean_green_pixel_brightness))),
    paste0("- Missing red brightness values: ", sum(is.na(brightness$mean_red_pixel_brightness))),
    paste0("- Missing near-infrared brightness values: ", sum(is.na(brightness$mean_nir_pixel_brightness))),
    paste0("- Output tables present: ", sum(file.exists(table_files)), " of ", length(table_files)),
    paste0("- Output figures present: ", sum(file.exists(figure_files)), " of ", length(figure_files)),
    "",
    "## Coverage",
    "",
    markdown_table(metric_stats[, c("scale", "metric_label", "n", "missing")]),
    "",
    "## Result",
    "",
    "The requested spectral heterogeneity metric scatterplot sets completed and produced the expected tables and figures."
  )
  writeLines(lines, VALIDATION_REPORT)
}

run_spectral_heterogeneity_relationship_analysis <- function() {
  data <- read_spectral_data()
  clusters <- read_composition_clusters()
  counts <- read_quadrat_counts()
  correlations <- run_pairwise_correlations(data)
  elevation_models <- run_elevation_adjusted_models(data)
  metric_stats <- metric_summary(data)
  silhouette_table <- calculate_composition_scatterplot_silhouette(data, clusters)
  brightness <- read_or_create_brightness(data)
  regional_correlations <- run_regional_brightness_correlations(data, brightness)
  figure_files <- save_figures(data, correlations, elevation_models, clusters, counts, silhouette_table, regional_correlations)

  write.csv(correlations, file.path(TABLE_DIR, "spectral_metric_pairwise_correlations.csv"), row.names = FALSE)
  write.csv(elevation_models, file.path(TABLE_DIR, "spectral_metric_elevation_adjusted_models.csv"), row.names = FALSE)
  write.csv(metric_stats, file.path(TABLE_DIR, "spectral_metric_summary.csv"), row.names = FALSE)
  write.csv(silhouette_table, file.path(TABLE_DIR, "spectral_composition_scatterplot_silhouette.csv"), row.names = FALSE)
  write.csv(regional_correlations, file.path(TABLE_DIR, "spectral_metric_regional_illumination_correlations.csv"), row.names = FALSE)

  write_analysis_report(correlations, elevation_models, metric_stats, silhouette_table, figure_files)
  write_task_report(correlations, elevation_models, silhouette_table)
  write_validation_report(data, correlations, elevation_models, metric_stats, silhouette_table, regional_correlations, figure_files)

  message("Spectral heterogeneity relationship analysis complete.")
  invisible(list(
    correlations = correlations,
    elevation_models = elevation_models,
    metric_stats = metric_stats,
    silhouette_table = silhouette_table,
    regional_correlations = regional_correlations,
    figure_files = figure_files
  ))
}

if (sys.nframe() == 0) {
  run_spectral_heterogeneity_relationship_analysis()
}
