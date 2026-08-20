PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
ANALYSIS_DATE <- "2026-08-19"

TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/final_research_direction")
FIGURE_DIR <- file.path(PROJECT_DIR, "reports/figures/final_research_direction")
ANALYSIS_REPORT <- file.path(PROJECT_DIR, "reports/analysis/20260819_final_research_direction_analysis.md")
TASK_REPORT <- file.path(PROJECT_DIR, "reports/tasks/20260819_final_research_direction_analysis.md")
VALIDATION_REPORT <- file.path(PROJECT_DIR, "reports/validation/20260819_final_research_direction_analysis_validation.md")

SCALES <- c("10m", "20m", "50m")
PCA_METRICS <- c(
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
NON_PCA_METRICS <- c("spec_sa")
SPECTRAL_METRICS <- c(PCA_METRICS, NON_PCA_METRICS)
BIODIVERSITY_METRICS <- c(
  "phy_faith",
  "phy_rao",
  "phy_afaith",
  "sp_rich",
  "sp_shannon",
  "sp_simpson",
  "sp_even"
)
PRIORITY_SPATIAL_SPECTRAL <- c("spec_spca_alpha", "spec_spca_mean", "spec_spca_rao")
PRIORITY_SPATIAL_BIODIVERSITY <- c("phy_rao", "phy_afaith", "sp_shannon")

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
  sp_rich = "Species richness",
  sp_shannon = "Shannon diversity",
  sp_simpson = "Simpson diversity",
  sp_even = "Species evenness",
  env_elev = "Mean elevation",
  env_tri5 = "TRI 5x5",
  env_tri11 = "TRI 11x11",
  env_tri21 = "TRI 21x21",
  crown_equivalent_individuals = "Crown-equivalent individuals",
  present_species_count = "Present species count",
  mean_pixel_brightness = "563 nm brightness",
  mean_blue_pixel_brightness = "Blue brightness",
  mean_green_pixel_brightness = "Green brightness",
  mean_red_pixel_brightness = "Red brightness",
  mean_nir_pixel_brightness = "NIR brightness"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
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
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", normalized_project), "/"), "", normalized_paths)
}

read_scale_table <- function(scale_name) {
  file_path <- file.path(PROJECT_DIR, paste0("quadrat_analysis_", scale_name, ".csv"))
  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  data$quad_id <- as.character(data$quad_id)
  data$scale <- as.character(data$scale)
  data
}

load_combined_data <- function() {
  data <- do.call(rbind, lapply(SCALES, read_scale_table))
  required <- c("quad_id", "scale", "center_x", "center_y", SPECTRAL_METRICS, BIODIVERSITY_METRICS, "env_elev")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data
}

safe_cor <- function(x, y) {
  complete <- is.finite(x) & is.finite(y)
  if (sum(complete) < 3 || stats::sd(x[complete]) == 0 || stats::sd(y[complete]) == 0) {
    return(NA_real_)
  }
  unname(stats::cor(x[complete], y[complete]))
}

safe_spearman <- function(x, y) {
  complete <- is.finite(x) & is.finite(y)
  if (sum(complete) < 3 || stats::sd(x[complete]) == 0 || stats::sd(y[complete]) == 0) {
    return(NA_real_)
  }
  suppressWarnings(unname(stats::cor(x[complete], y[complete], method = "spearman")))
}

fit_relationship <- function(data, scale_name, x_metric, y_metric, x_transform = "identity", y_transform = "identity") {
  scale_data <- data[data$scale == scale_name, , drop = FALSE]
  x <- transform_values(scale_data[[x_metric]], x_transform)
  y <- transform_values(scale_data[[y_metric]], y_transform)
  complete <- is.finite(x) & is.finite(y)
  row <- data.frame(
    scale = scale_name,
    x_metric = x_metric,
    x_label = display_name(x_metric),
    y_metric = y_metric,
    y_label = display_name(y_metric),
    x_transform = x_transform,
    y_transform = y_transform,
    n = sum(complete),
    pearson_r = NA_real_,
    spearman_r = NA_real_,
    r_squared = NA_real_,
    f_statistic = NA_real_,
    f_p_value = NA_real_,
    slope = NA_real_,
    intercept = NA_real_,
    residual_shapiro_p = NA_real_,
    residual_abs_cor_fitted = NA_real_,
    stringsAsFactors = FALSE
  )
  if (sum(complete) >= 4 && stats::sd(x[complete]) > 0 && stats::sd(y[complete]) > 0) {
    fit <- stats::lm(y[complete] ~ x[complete])
    fit_summary <- summary(fit)
    f_values <- fit_summary$fstatistic
    residuals <- stats::residuals(fit)
    fitted <- stats::fitted(fit)
    row$pearson_r <- safe_cor(x, y)
    row$spearman_r <- safe_spearman(x, y)
    row$r_squared <- fit_summary$r.squared
    row$f_statistic <- unname(f_values[["value"]])
    row$f_p_value <- stats::pf(f_values[["value"]], f_values[["numdf"]], f_values[["dendf"]], lower.tail = FALSE)
    row$slope <- unname(stats::coef(fit)[[2]])
    row$intercept <- unname(stats::coef(fit)[[1]])
    if (length(residuals) >= 3 && length(residuals) <= 5000) {
      row$residual_shapiro_p <- tryCatch(stats::shapiro.test(residuals)$p.value, error = function(e) NA_real_)
    }
    row$residual_abs_cor_fitted <- safe_cor(abs(residuals), fitted)
  }
  row
}

transform_values <- function(x, transform_name) {
  x <- as.numeric(x)
  if (transform_name == "identity") {
    return(x)
  }
  if (transform_name == "log") {
    if (any(x <= 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(log(x))
  }
  if (transform_name == "log1p") {
    if (any(x < 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(log1p(x))
  }
  if (transform_name == "sqrt") {
    if (any(x < 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(sqrt(x))
  }
  if (transform_name == "inverse") {
    if (any(x == 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(1 / x)
  }
  if (transform_name == "power_0.25") {
    if (any(x < 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(x^0.25)
  }
  if (transform_name == "power_0.5") {
    if (any(x < 0, na.rm = TRUE)) return(rep(NA_real_, length(x)))
    return(x^0.5)
  }
  if (transform_name == "power_2") {
    return(x^2)
  }
  if (transform_name == "boxcox") {
    return(boxcox_transform(x))
  }
  if (transform_name == "yeojohnson") {
    return(yeojohnson_transform(x))
  }
  stop("Unknown transform: ", transform_name, call. = FALSE)
}

boxcox_transform <- function(x) {
  if (any(x <= 0, na.rm = TRUE)) {
    return(rep(NA_real_, length(x)))
  }
  lambda <- estimate_boxcox_lambda(x)
  if (is.na(lambda)) {
    return(rep(NA_real_, length(x)))
  }
  if (abs(lambda) < 1e-8) {
    log(x)
  } else {
    ((x^lambda) - 1) / lambda
  }
}

estimate_boxcox_lambda <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || any(x <= 0) || stats::sd(x) == 0) {
    return(NA_real_)
  }
  lambdas <- seq(-2, 2, by = 0.1)
  ll <- vapply(lambdas, function(lambda) {
    z <- if (abs(lambda) < 1e-8) log(x) else ((x^lambda) - 1) / lambda
    if (!all(is.finite(z)) || stats::var(z) <= 0) return(-Inf)
    -length(x) / 2 * log(stats::var(z)) + (lambda - 1) * sum(log(x))
  }, numeric(1))
  lambdas[which.max(ll)]
}

yeojohnson_transform <- function(x) {
  lambda <- estimate_yeojohnson_lambda(x)
  if (is.na(lambda)) {
    return(rep(NA_real_, length(x)))
  }
  yj_apply(x, lambda)
}

yj_apply <- function(x, lambda) {
  out <- rep(NA_real_, length(x))
  pos <- is.finite(x) & x >= 0
  neg <- is.finite(x) & x < 0
  if (abs(lambda) < 1e-8) {
    out[pos] <- log1p(x[pos])
  } else {
    out[pos] <- (((x[pos] + 1)^lambda) - 1) / lambda
  }
  if (abs(lambda - 2) < 1e-8) {
    out[neg] <- -log1p(-x[neg])
  } else {
    out[neg] <- -((((1 - x[neg])^(2 - lambda)) - 1) / (2 - lambda))
  }
  out
}

estimate_yeojohnson_lambda <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || stats::sd(x) == 0) {
    return(NA_real_)
  }
  lambdas <- seq(-2, 2, by = 0.1)
  ll <- vapply(lambdas, function(lambda) {
    z <- yj_apply(x, lambda)
    if (!all(is.finite(z)) || stats::var(z) <= 0) return(-Inf)
    jacobian <- sum(ifelse(x >= 0, log1p(x), -log1p(-x)))
    -length(x) / 2 * log(stats::var(z)) + (lambda - 1) * jacobian
  }, numeric(1))
  lambdas[which.max(ll)]
}

make_all_spectral_biodiversity_relationships <- function(data) {
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (spectral_metric in SPECTRAL_METRICS) {
      for (biodiversity_metric in BIODIVERSITY_METRICS) {
        rows[[index]] <- fit_relationship(
          data,
          scale_name,
          x_metric = biodiversity_metric,
          y_metric = spectral_metric
        )
        index <- index + 1
      }
    }
  }
  out <- do.call(rbind, rows)
  out$spectral_group <- ifelse(out$y_metric %in% NON_PCA_METRICS, "non_pca", ifelse(grepl("^spec_spca", out$y_metric), "standardized_pca", "raw_pca"))
  out$biodiversity_group <- ifelse(grepl("^phy_", out$x_metric), "phylogenetic", "species")
  out
}

make_pairwise_metric_relationships <- function(data, metrics, group_name) {
  pairs <- utils::combn(metrics, 2)
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (i in seq_len(ncol(pairs))) {
      rows[[index]] <- fit_relationship(
        data,
        scale_name,
        x_metric = pairs[1, i],
        y_metric = pairs[2, i]
      )
      rows[[index]]$metric_group <- group_name
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

make_transformation_relationships <- function(data) {
  transforms <- c("identity", "log", "log1p", "sqrt", "power_0.25", "power_0.5", "power_2", "boxcox", "yeojohnson")
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (spectral_metric in PCA_METRICS) {
      for (biodiversity_metric in BIODIVERSITY_METRICS) {
        for (x_transform in transforms) {
          for (y_transform in transforms) {
            rows[[index]] <- fit_relationship(
              data,
              scale_name,
              x_metric = biodiversity_metric,
              y_metric = spectral_metric,
              x_transform = x_transform,
              y_transform = y_transform
            )
            index <- index + 1
          }
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$abs_r <- abs(out$pearson_r)
  out
}

make_transformation_summary <- function(transform_results) {
  identity_rows <- transform_results[
    transform_results$x_transform == "identity" & transform_results$y_transform == "identity",
    c("scale", "x_metric", "y_metric", "n", "pearson_r", "spearman_r", "r_squared", "residual_abs_cor_fitted")
  ]
  names(identity_rows) <- c(
    "scale", "x_metric", "y_metric", "n",
    "identity_pearson_r", "identity_spearman_r", "identity_r_squared",
    "identity_residual_abs_cor_fitted"
  )
  best_rows <- do.call(rbind, lapply(split(transform_results, list(transform_results$scale, transform_results$x_metric, transform_results$y_metric), drop = TRUE), function(group) {
    group <- group[!is.na(group$abs_r), , drop = FALSE]
    if (nrow(group) == 0) return(NULL)
    group[order(-group$abs_r, -group$r_squared, group$residual_abs_cor_fitted), ][1, ]
  }))
  out <- merge(
    best_rows,
    identity_rows,
    by = c("scale", "x_metric", "y_metric", "n"),
    all.x = TRUE,
    sort = FALSE
  )
  out$delta_abs_r <- abs(out$pearson_r) - abs(out$identity_pearson_r)
  out$delta_r_squared <- out$r_squared - out$identity_r_squared
  out$best_transform_pair <- paste0(out$x_transform, " x / ", out$y_transform, " y")
  out[order(-out$delta_abs_r, -out$delta_r_squared), ]
}

read_optional_context <- function(data) {
  count_file <- file.path(PROJECT_DIR, "reports/tables/species_phylogenetic_correlation/quadrat_crown_equivalent_individual_totals.csv")
  cluster_file <- file.path(PROJECT_DIR, "reports/tables/species_phylogenetic_correlation/species_composition_clusters.csv")
  brightness_file <- file.path(PROJECT_DIR, "reports/tables/spectral_heterogeneity_relationships/quadrat_pixel_brightness_summary.csv")
  if (file.exists(count_file)) {
    counts <- read.csv(count_file, stringsAsFactors = FALSE, check.names = FALSE)
    counts$quad_id <- as.character(counts$quad_id)
    counts$scale <- as.character(counts$scale)
    data <- merge(data, counts[, c("quad_id", "scale", "crown_equivalent_individuals", "present_species_count")], by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
  }
  if (file.exists(cluster_file)) {
    clusters <- read.csv(cluster_file, stringsAsFactors = FALSE, check.names = FALSE)
    clusters$quad_id <- as.character(clusters$quad_id)
    clusters$scale <- as.character(clusters$scale)
    keep_cols <- intersect(c("quad_id", "scale", "composition_cluster", "composition_type"), names(clusters))
    data <- merge(data, clusters[, keep_cols], by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
  }
  if (file.exists(brightness_file)) {
    brightness <- read.csv(brightness_file, stringsAsFactors = FALSE, check.names = FALSE)
    brightness$quad_id <- as.character(brightness$quad_id)
    brightness$scale <- as.character(brightness$scale)
    keep_cols <- intersect(
      c("quad_id", "scale", "brightness_n_pixels", "mean_pixel_brightness", "mean_blue_pixel_brightness", "mean_green_pixel_brightness", "mean_red_pixel_brightness", "mean_nir_pixel_brightness"),
      names(brightness)
    )
    data <- merge(data, brightness[, keep_cols], by = c("quad_id", "scale"), all.x = TRUE, sort = FALSE)
  }
  data
}

make_driver_relationships <- function(data) {
  spectral_drivers <- c(
    "env_elev", "env_tri5", "env_tri11", "env_tri21",
    "brightness_n_pixels", "mean_pixel_brightness",
    "mean_blue_pixel_brightness", "mean_green_pixel_brightness",
    "mean_red_pixel_brightness", "mean_nir_pixel_brightness"
  )
  biodiversity_drivers <- c("present_species_count", "crown_equivalent_individuals", "sp_even", "sp_simpson")
  spectral_drivers <- intersect(spectral_drivers, names(data))
  biodiversity_drivers <- intersect(biodiversity_drivers, names(data))
  rows <- list()
  index <- 1
  for (scale_name in SCALES) {
    for (metric in SPECTRAL_METRICS) {
      for (driver in spectral_drivers) {
        rows[[index]] <- fit_relationship(data, scale_name, x_metric = driver, y_metric = metric)
        rows[[index]]$target_group <- "spectral"
        index <- index + 1
      }
    }
    for (metric in BIODIVERSITY_METRICS) {
      for (driver in biodiversity_drivers) {
        if (identical(metric, driver)) next
        rows[[index]] <- fit_relationship(data, scale_name, x_metric = driver, y_metric = metric)
        rows[[index]]$target_group <- "biodiversity"
        index <- index + 1
      }
    }
  }
  out <- do.call(rbind, rows)
  out$driver_metric <- out$x_metric
  out$target_metric <- out$y_metric
  out
}

make_pca_method_summary <- function() {
  rows <- list()
  for (scale_name in SCALES) {
    file_path <- file.path(PROJECT_DIR, "Quad_Values", paste0(scale_name, "_spectral_heterogeneity_smooth_masked_5nm_summary.csv"))
    data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    method_fields <- c(
      raw_pca_metric = "metric_method",
      standardized_pca_metric = "standardized_PCA_metric_method",
      raw_alpha_hull = "alpha_hull_method",
      standardized_alpha_hull = "standardized_PCA_alpha_hull_method",
      raw_3d_hull = "pca_hull_3d_method",
      standardized_3d_hull = "standardized_PCA_pca_hull_3d_method"
    )
    scale_rows <- list()
    for (field_name in names(method_fields)) {
      field <- method_fields[[field_name]]
      counts <- table(data[[field]], useNA = "ifany")
      scale_rows[[field_name]] <- data.frame(
        scale = scale_name,
        field_group = field_name,
        field = field,
        method = names(counts),
        n = as.integer(counts),
        stringsAsFactors = FALSE
      )
    }
    rows[[scale_name]] <- do.call(rbind, scale_rows)
  }
  do.call(rbind, rows)
}

moran_i <- function(values, x_coord, y_coord, k = 8, n_perm = 199, seed = 20260819) {
  complete <- is.finite(values) & is.finite(x_coord) & is.finite(y_coord)
  values <- values[complete]
  x_coord <- x_coord[complete]
  y_coord <- y_coord[complete]
  n <- length(values)
  if (n <= k + 2 || stats::sd(values) == 0) {
    return(data.frame(n = n, k = k, moran_i = NA_real_, permutation_p = NA_real_))
  }
  coords <- cbind(x_coord, y_coord)
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  neighbors <- t(apply(d, 1, function(row) order(row)[seq_len(k)]))
  edge_i <- rep(seq_len(n), each = k)
  edge_j <- as.vector(t(neighbors))
  edge_key <- paste(pmin(edge_i, edge_j), pmax(edge_i, edge_j), sep = "_")
  keep_edge <- !duplicated(edge_key)
  edge_i <- edge_i[keep_edge]
  edge_j <- edge_j[keep_edge]
  w_sum <- length(edge_i) * 2
  centered <- values - mean(values)
  observed <- n / w_sum * (2 * sum(centered[edge_i] * centered[edge_j])) / sum(centered^2)
  set.seed(seed)
  permuted <- replicate(n_perm, {
    z <- sample(centered)
    n / w_sum * (2 * sum(z[edge_i] * z[edge_j])) / sum(z^2)
  })
  p_value <- (sum(abs(permuted) >= abs(observed)) + 1) / (n_perm + 1)
  data.frame(n = n, k = k, moran_i = observed, permutation_p = p_value)
}

semivariogram_summary <- function(values, x_coord, y_coord, n_bins = 10, max_pairs = 200000, seed = 20260819) {
  complete <- is.finite(values) & is.finite(x_coord) & is.finite(y_coord)
  values <- values[complete]
  x_coord <- x_coord[complete]
  y_coord <- y_coord[complete]
  n <- length(values)
  if (n < 5 || stats::sd(values) == 0) {
    return(data.frame())
  }
  pair_count <- n * (n - 1) / 2
  if (pair_count > max_pairs) {
    set.seed(seed)
    pair_index <- cbind(sample.int(n, max_pairs, replace = TRUE), sample.int(n, max_pairs, replace = TRUE))
    pair_index <- pair_index[pair_index[, 1] != pair_index[, 2], , drop = FALSE]
  } else {
    pair_index <- utils::combn(n, 2)
    pair_index <- t(pair_index)
  }
  distances <- sqrt((x_coord[pair_index[, 1]] - x_coord[pair_index[, 2]])^2 + (y_coord[pair_index[, 1]] - y_coord[pair_index[, 2]])^2)
  semivariance <- 0.5 * (values[pair_index[, 1]] - values[pair_index[, 2]])^2
  breaks <- unique(stats::quantile(distances, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(breaks) < 3) {
    return(data.frame())
  }
  bins <- cut(distances, breaks = breaks, include.lowest = TRUE)
  aggregate(
    data.frame(distance = distances, semivariance = semivariance),
    by = list(distance_bin = bins),
    FUN = function(z) c(mean = mean(z), median = stats::median(z), n = length(z))
  )
}

make_spatial_diagnostics <- function(data) {
  rows <- list()
  variogram_rows <- list()
  index <- 1
  vindex <- 1
  for (scale_name in SCALES) {
    scale_data <- data[data$scale == scale_name, , drop = FALSE]
    variables <- unique(c(PRIORITY_SPATIAL_SPECTRAL, PRIORITY_SPATIAL_BIODIVERSITY))
    for (variable in variables) {
      result <- moran_i(scale_data[[variable]], scale_data$center_x, scale_data$center_y)
      result$scale <- scale_name
      result$diagnostic_type <- "variable"
      result$variable <- variable
      result$variable_label <- display_name(variable)
      result$model <- NA_character_
      rows[[index]] <- result
      index <- index + 1
    }
    for (spectral_metric in PRIORITY_SPATIAL_SPECTRAL) {
      for (biodiversity_metric in PRIORITY_SPATIAL_BIODIVERSITY) {
        complete <- is.finite(scale_data[[spectral_metric]]) & is.finite(scale_data[[biodiversity_metric]])
        if (sum(complete) >= 4) {
          fit <- lm(scale_data[[spectral_metric]][complete] ~ scale_data[[biodiversity_metric]][complete])
          residuals <- rep(NA_real_, nrow(scale_data))
          residuals[which(complete)] <- stats::residuals(fit)
          model_name <- paste0(spectral_metric, " ~ ", biodiversity_metric)
          result <- moran_i(residuals, scale_data$center_x, scale_data$center_y)
          result$scale <- scale_name
          result$diagnostic_type <- "residual"
          result$variable <- NA_character_
          result$variable_label <- NA_character_
          result$model <- model_name
          rows[[index]] <- result
          index <- index + 1
          if (spectral_metric %in% c("spec_spca_alpha", "spec_spca_mean") && biodiversity_metric %in% c("phy_rao", "phy_afaith")) {
            vg <- semivariogram_summary(residuals, scale_data$center_x, scale_data$center_y)
            if (nrow(vg) > 0) {
              flat <- data.frame(
                scale = scale_name,
                model = model_name,
                distance_bin = vg$distance_bin,
                distance_mean = vg$distance[, "mean"],
                distance_median = vg$distance[, "median"],
                pair_n = vg$distance[, "n"],
                semivariance_mean = vg$semivariance[, "mean"],
                semivariance_median = vg$semivariance[, "median"],
                stringsAsFactors = FALSE
              )
              variogram_rows[[vindex]] <- flat
              vindex <- vindex + 1
            }
          }
        }
      }
    }
  }
  list(
    moran = do.call(rbind, rows),
    variogram = if (length(variogram_rows) == 0) data.frame() else do.call(rbind, variogram_rows)
  )
}

plot_transformation_improvement <- function(summary_data) {
  plot_data <- summary_data[summary_data$y_metric %in% c("spec_spca_alpha", "spec_spca_mean", "spec_spca_rao") &
                              summary_data$x_metric %in% c("phy_rao", "phy_afaith", "sp_shannon"), ]
  file_path <- file.path(FIGURE_DIR, "01_transformation_delta_abs_r_priority_relationships.png")
  png(file_path, width = 3000, height = 1800, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(8, 4.5, 2.5, 1))
  labels <- paste(plot_data$scale, display_name(plot_data$y_metric), "vs", display_name(plot_data$x_metric), sep = "\n")
  bar_cols <- ifelse(plot_data$delta_abs_r > 0.05, "#1B6E6A", "#8AA7A9")
  barplot(
    plot_data$delta_abs_r,
    names.arg = labels,
    las = 2,
    cex.names = 0.48,
    col = bar_cols,
    border = NA,
    ylab = "Best transform gain in |Pearson r|",
    main = "Transformation gains for priority relationships"
  )
  abline(h = 0, lwd = 1)
  dev.off()
  file_path
}

plot_spatial_moran <- function(moran_data) {
  plot_data <- moran_data[!is.na(moran_data$moran_i), ]
  file_path <- file.path(FIGURE_DIR, "02_spatial_moran_i_priority_variables_and_residuals.png")
  png(file_path, width = 3200, height = 2200, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(8, 4.5, 3, 1))
  for (diag_type in c("variable", "residual")) {
    subset_data <- plot_data[plot_data$diagnostic_type == diag_type, ]
    subset_data <- subset_data[order(subset_data$scale, subset_data$moran_i), ]
    labels <- if (diag_type == "variable") {
      paste(subset_data$scale, subset_data$variable_label, sep = "\n")
    } else {
      paste(subset_data$scale, subset_data$model, sep = "\n")
    }
    cols <- ifelse(subset_data$permutation_p < 0.05, "#8B3A3A", "#7B8EA0")
    barplot(
      subset_data$moran_i,
      names.arg = labels,
      las = 2,
      cex.names = ifelse(diag_type == "variable", 0.55, 0.38),
      col = cols,
      border = NA,
      ylab = "Moran's I",
      main = ifelse(diag_type == "variable", "Variables", "Model residuals")
    )
    abline(h = 0, lwd = 1)
  }
  dev.off()
  file_path
}

plot_metric_heatmap <- function(relationship_data, metrics, file_name, title) {
  plot_data <- relationship_data[relationship_data$x_metric %in% metrics & relationship_data$y_metric %in% metrics, ]
  file_path <- file.path(FIGURE_DIR, file_name)
  png(file_path, width = 3000, height = 1800, res = 300)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 3), mar = c(8, 8, 3, 1))
  for (scale_name in SCALES) {
    mat <- matrix(NA_real_, nrow = length(metrics), ncol = length(metrics), dimnames = list(display_name(metrics), display_name(metrics)))
    diag(mat) <- 1
    for (i in seq_len(nrow(plot_data))) {
      row <- plot_data[i, ]
      if (row$scale == scale_name) {
        mat[display_name(row$y_metric), display_name(row$x_metric)] <- row$pearson_r
        mat[display_name(row$x_metric), display_name(row$y_metric)] <- row$pearson_r
      }
    }
    image(
      seq_len(ncol(mat)),
      seq_len(nrow(mat)),
      t(mat[nrow(mat):1, ]),
      zlim = c(-1, 1),
      col = colorRampPalette(c("#7A2E2E", "#F7F7F2", "#245E8A"))(101),
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = paste(title, scale_name)
    )
    axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.55)
    axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.55)
    box()
  }
  dev.off()
  file_path
}

write_outputs <- function(outputs) {
  write.csv(outputs$all_spectral_biodiversity, file.path(TABLE_DIR, "all_spectral_biodiversity_relationships.csv"), row.names = FALSE)
  write.csv(outputs$spectral_pairwise, file.path(TABLE_DIR, "spectral_metric_pairwise_relationships.csv"), row.names = FALSE)
  write.csv(outputs$biodiversity_pairwise, file.path(TABLE_DIR, "biodiversity_metric_pairwise_relationships.csv"), row.names = FALSE)
  write.csv(outputs$transformation_results, file.path(TABLE_DIR, "transformation_spectral_biodiversity_relationships.csv"), row.names = FALSE)
  write.csv(outputs$transformation_summary, file.path(TABLE_DIR, "transformation_best_relationship_summary.csv"), row.names = FALSE)
  write.csv(outputs$driver_relationships, file.path(TABLE_DIR, "metric_driver_relationships.csv"), row.names = FALSE)
  write.csv(outputs$pca_method_summary, file.path(TABLE_DIR, "pca_metric_computation_method_summary.csv"), row.names = FALSE)
  write.csv(outputs$spatial$moran, file.path(TABLE_DIR, "spatial_moran_diagnostics.csv"), row.names = FALSE)
  write.csv(outputs$spatial$variogram, file.path(TABLE_DIR, "spatial_variogram_summary.csv"), row.names = FALSE)
}

make_report <- function(outputs, figure_files) {
  all_rel <- outputs$all_spectral_biodiversity
  top_pca <- all_rel[all_rel$y_metric %in% PCA_METRICS, ]
  top_pca <- top_pca[order(-abs(top_pca$pearson_r)), ]
  top_pca <- head(top_pca, 20)
  transform_top <- outputs$transformation_summary[order(-outputs$transformation_summary$delta_abs_r), ]
  transform_top <- head(transform_top, 15)
  driver_top <- outputs$driver_relationships[order(-abs(outputs$driver_relationships$pearson_r)), ]
  driver_top <- head(driver_top, 20)
  spatial_sig <- outputs$spatial$moran[order(outputs$spatial$moran$permutation_p, -abs(outputs$spatial$moran$moran_i)), ]
  spatial_sig <- head(spatial_sig, 20)
  pca_method_summary <- outputs$pca_method_summary
  pca_method_summary <- pca_method_summary[pca_method_summary$method != "NA", ]

  top_pca_table <- data.frame(
    Scale = top_pca$scale,
    Spectral = top_pca$y_label,
    Biodiversity = top_pca$x_label,
    n = top_pca$n,
    r = fmt_num(top_pca$pearson_r),
    R2 = fmt_num(top_pca$r_squared),
    p = fmt_p(top_pca$f_p_value),
    stringsAsFactors = FALSE
  )
  transform_table <- data.frame(
    Scale = transform_top$scale,
    Spectral = transform_top$y_label,
    Biodiversity = transform_top$x_label,
    Best = transform_top$best_transform_pair,
    Identity_r = fmt_num(transform_top$identity_pearson_r),
    Best_r = fmt_num(transform_top$pearson_r),
    Delta_abs_r = fmt_num(transform_top$delta_abs_r),
    stringsAsFactors = FALSE
  )
  driver_table <- data.frame(
    Scale = driver_top$scale,
    Target = driver_top$y_label,
    Driver = driver_top$x_label,
    Group = driver_top$target_group,
    n = driver_top$n,
    r = fmt_num(driver_top$pearson_r),
    R2 = fmt_num(driver_top$r_squared),
    p = fmt_p(driver_top$f_p_value),
    stringsAsFactors = FALSE
  )
  spatial_table <- data.frame(
    Scale = spatial_sig$scale,
    Type = spatial_sig$diagnostic_type,
    Variable_or_model = ifelse(is.na(spatial_sig$model), spatial_sig$variable_label, spatial_sig$model),
    n = spatial_sig$n,
    Moran_I = fmt_num(spatial_sig$moran_i),
    p = fmt_p(spatial_sig$permutation_p),
    stringsAsFactors = FALSE
  )
  pca_method_table <- data.frame(
    Scale = pca_method_summary$scale,
    Field = pca_method_summary$field_group,
    Method = pca_method_summary$method,
    n = pca_method_summary$n,
    stringsAsFactors = FALSE
  )

  lines <- c(
    "# Final Research Direction Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Purpose",
    "",
    "This analysis completes the current `Checklist.md` analysis-needed items that can be completed from existing 10 m, 20 m, and 50 m outputs. It compiles PCA-centered spectral-biodiversity relationships, compares spectral and biodiversity metrics, tests transformations, summarizes metric drivers, and runs spatial autocorrelation diagnostics for priority variables and residuals.",
    "",
    "The discarded 100 m optimal-scale task is not run because no 100 m products are currently part of the validated workflow.",
    "",
    "## Top PCA Spectral-Biodiversity Relationships",
    "",
    markdown_table(top_pca_table),
    "",
    "## Largest Transformation Gains",
    "",
    markdown_table(transform_table),
    "",
    "## Strongest Metric Driver Relationships",
    "",
    markdown_table(driver_table),
    "",
    "## Strongest Spatial Autocorrelation Diagnostics",
    "",
    markdown_table(spatial_table),
    "",
    "## PCA Metric Computation Method Check",
    "",
    markdown_table(pca_method_table),
    "",
    "## Figures",
    "",
    paste0("- `", relative_path(figure_files), "`"),
    "",
    "## Output Tables",
    "",
    paste0("- `", relative_path(list.files(TABLE_DIR, pattern = "\\.csv$", full.names = TRUE)), "`"),
    "",
    "## Interpretation Notes",
    "",
    "- Vector-normalized PCA metrics are reported as the primary PCA evidence layer, with raw PCA metrics retained as brightness-sensitive comparisons.",
    "- Transformation gains should be interpreted as diagnostics, not as automatic replacements for untransformed metrics. A transform is useful only if it improves correlation and residual behavior while remaining interpretable.",
    "- Moran's I diagnostics use eight-nearest-neighbor spatial weights and 199 permutations. Significant residual autocorrelation means coefficient p-values from ordinary least squares should remain cautious.",
    "- Driver relationships are screening summaries. They identify plausible influences on metric values but do not by themselves establish causation."
  )
  writeLines(lines, ANALYSIS_REPORT)
}

make_task_report <- function() {
  lines <- c(
    "# Task Report: Final Research Direction Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Objective",
    "",
    "Complete the analysis-needed items in `Checklist.md` that can be addressed from current validated 10 m, 20 m, and 50 m outputs.",
    "",
    "## Outputs Created Or Replaced",
    "",
    paste0("- `", relative_path(ANALYSIS_REPORT), "`"),
    paste0("- `", relative_path(TABLE_DIR), "/`"),
    paste0("- `", relative_path(FIGURE_DIR), "/`"),
    "",
    "## Notes",
    "",
    "- The standalone optimal-scale task remains discarded because the checklist now references 100 m and no current 100 m workflow exists.",
    "- Existing heavy raster-derived outputs were reused where possible to avoid redoing upstream processing.",
    "- The analysis script overwrites the final-research-direction tables and figures on rerun."
  )
  writeLines(lines, TASK_REPORT)
}

make_validation_report <- function(outputs, figure_files) {
  expected_tables <- c(
    "all_spectral_biodiversity_relationships.csv",
    "spectral_metric_pairwise_relationships.csv",
    "biodiversity_metric_pairwise_relationships.csv",
    "transformation_spectral_biodiversity_relationships.csv",
    "transformation_best_relationship_summary.csv",
    "metric_driver_relationships.csv",
    "pca_metric_computation_method_summary.csv",
    "spatial_moran_diagnostics.csv",
    "spatial_variogram_summary.csv"
  )
  table_paths <- file.path(TABLE_DIR, expected_tables)
  lines <- c(
    "# Validation: Final Research Direction Analysis",
    "",
    paste0("Date: ", ANALYSIS_DATE),
    "",
    "## Checks",
    "",
    paste0("- All spectral-biodiversity rows: ", nrow(outputs$all_spectral_biodiversity)),
    paste0("- Spectral pairwise rows: ", nrow(outputs$spectral_pairwise)),
    paste0("- Biodiversity pairwise rows: ", nrow(outputs$biodiversity_pairwise)),
    paste0("- Transformation rows: ", nrow(outputs$transformation_results)),
    paste0("- Transformation summary rows: ", nrow(outputs$transformation_summary)),
    paste0("- Driver relationship rows: ", nrow(outputs$driver_relationships)),
    paste0("- PCA computation method rows: ", nrow(outputs$pca_method_summary)),
    paste0("- Moran diagnostic rows: ", nrow(outputs$spatial$moran)),
    paste0("- Variogram rows: ", nrow(outputs$spatial$variogram)),
    paste0("- Expected tables present: ", sum(file.exists(table_paths)), " of ", length(table_paths)),
    paste0("- Figures present: ", sum(file.exists(figure_files)), " of ", length(figure_files)),
    "",
    "## Result",
    "",
    ifelse(all(file.exists(table_paths)) && all(file.exists(figure_files)), "The final research direction analyses completed and produced the expected tables and figures.", "Some expected outputs are missing; inspect the script log and output folders.")
  )
  writeLines(lines, VALIDATION_REPORT)
}

main <- function() {
  data <- load_combined_data()
  data <- read_optional_context(data)
  all_spectral_biodiversity <- make_all_spectral_biodiversity_relationships(data)
  spectral_pairwise <- make_pairwise_metric_relationships(data, SPECTRAL_METRICS, "spectral")
  biodiversity_pairwise <- make_pairwise_metric_relationships(data, BIODIVERSITY_METRICS, "biodiversity")
  transformation_results <- make_transformation_relationships(data)
  transformation_summary <- make_transformation_summary(transformation_results)
  driver_relationships <- make_driver_relationships(data)
  pca_method_summary <- make_pca_method_summary()
  spatial <- make_spatial_diagnostics(data)

  outputs <- list(
    all_spectral_biodiversity = all_spectral_biodiversity,
    spectral_pairwise = spectral_pairwise,
    biodiversity_pairwise = biodiversity_pairwise,
    transformation_results = transformation_results,
    transformation_summary = transformation_summary,
    driver_relationships = driver_relationships,
    pca_method_summary = pca_method_summary,
    spatial = spatial
  )
  write_outputs(outputs)
  figure_files <- c(
    plot_transformation_improvement(transformation_summary),
    plot_spatial_moran(spatial$moran),
    plot_metric_heatmap(spectral_pairwise, SPECTRAL_METRICS, "03_spectral_metric_correlation_heatmap.png", "Spectral metrics"),
    plot_metric_heatmap(biodiversity_pairwise, BIODIVERSITY_METRICS, "04_biodiversity_metric_correlation_heatmap.png", "Biodiversity metrics")
  )
  make_report(outputs, figure_files)
  make_task_report()
  make_validation_report(outputs, figure_files)
  stray_plot <- file.path(PROJECT_DIR, "Rplots.pdf")
  if (file.exists(stray_plot)) {
    unlink(stray_plot)
  }
}

main()
