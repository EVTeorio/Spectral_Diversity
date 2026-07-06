user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

# Calculates the current all-scale spectral heterogeneity products from
# smoothed 5 nm quadrat spectra. The workflow adds three metrics to the
# existing spectral-angle entropy product:
#   1. mean Euclidean distance from the quadrat centroid in global PCA space
#   2. alpha-hull area in global PC1-PC2 space
#   3. Rao's Q using equal pixel weights and squared Euclidean distance in
#      global PC1-PC3 space

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
SIDECAR_PATTERN <- "\\.(hdr|aux|xml|enp|sta)$"
OUTPUT_DIR <- "Quad_Values"
OUTPUT_SHP_DIR <- "Quad_Values/Spectral_diversitySHPs"
FIGURE_DIR <- "reports/figures/spectral_heterogeneity"
PROGRESS_LOG <- "logs/spectral_heterogeneity_all_metrics_progress.log"
PCA_RDS <- file.path(OUTPUT_SHP_DIR, "global_pca_smooth_masked_5nm.rds")
PCA_VARIANCE_CSV <- file.path(OUTPUT_SHP_DIR, "global_pca_smooth_masked_5nm_variance_explained.csv")
PCA_VARIANCE_PNG <- file.path(FIGURE_DIR, "global_pca_smooth_masked_5nm_variance_explained.png")

# Current shadow mask carried forward from the validated SA entropy workflow.
BEST_BAND <- "563 nm"
SHADOW_THRESHOLD <- 0.0305476
DIRECTION <- "<"

PCA_SAMPLE_PER_RASTER <- as.integer(Sys.getenv("PCA_SAMPLE_PER_RASTER", "500"))
PCA_MAX_ROWS <- as.integer(Sys.getenv("PCA_MAX_ROWS", "800000"))
PCA_N_COMPONENTS <- as.integer(Sys.getenv("PCA_N_COMPONENTS", "3"))
MIN_VALID_PIXELS <- as.integer(Sys.getenv("PCA_MIN_VALID_PIXELS", "10"))

ALPHA_HULL_ALPHA <- as.numeric(Sys.getenv("ALPHA_HULL_ALPHA", "1"))
ALPHA_MAX_POINTS <- as.integer(Sys.getenv("ALPHA_MAX_POINTS", "10000"))
ALPHA_FALLBACK_POINTS <- as.integer(Sys.getenv("ALPHA_FALLBACK_POINTS", "5000"))
ALPHA_TIMEOUT <- as.numeric(Sys.getenv("ALPHA_TIMEOUT", "10"))
HULL_3D_MAX_POINTS <- as.integer(Sys.getenv("HULL_3D_MAX_POINTS", "10000"))

RANDOM_SEED <- as.integer(Sys.getenv("SPECTRAL_HET_RANDOM_SEED", "42"))
PCA_EXCLUSION_POLICY_ID <- "manual_atmospheric_cloud_exclusions_20260622"

EXCLUDED_20M_ATMOSPHERIC <- as.character(c(
  1424, 1423, 1422, 1420, 1421, 1419, 1418, 1414,
  1521, 1522, 1523, 1524, 1520, 1519,
  1624, 1622, 1623, 1621, 1620,
  1724, 1723, 1722, 1721,
  1824, 1823, 1822,
  1923, 1924, 1922, 1921,
  1322, 1321, 1319, 1320, 1318, 1317, 1316, 1315, 1314, 1313,
  1221, 1220, 1219, 1216, 1215, 1213, 1214, 1212, 1211,
  1120, 1119, 1115, 1114, 1113, 1112, 1111, 1110,
  1014, 1013, 1010, 1009,
  909, 908, 24
))

EXCLUDED_50M_ATMOSPHERIC <- c(
  "sub50_80", "sub50_79", "sub50_71", "sub50_70", "sub50_62", "sub50_53"
)

SCALE_CONFIG <- list(
  list(
    scale = "10m",
    spec_dir = "Quad_Spectra/10m_smooth_5nm",
    shp_path = "Quad_Scale_SHPs/PR_10m.shp",
    join_field = "sub_id",
    n_cores = 8
  ),
  list(
    scale = "20m",
    spec_dir = "Quad_Spectra/20m_smooth_5nm",
    shp_path = "Quad_Scale_SHPs/PR_20m.shp",
    join_field = "Name",
    n_cores = 4
  ),
  list(
    scale = "50m",
    spec_dir = "Quad_Spectra/50m_smooth_5nm",
    shp_path = "Quad_Scale_SHPs/PR_50m.shp",
    join_field = "Name",
    n_cores = 2
  )
)

log_progress <- function(message) {
  text <- paste(Sys.time(), message)
  cat(text, "\n")
  cat(text, "\n", file = PROGRESS_LOG, append = TRUE)
}

require_runtime_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Required R package is not installed or visible on .libPaths(): ",
      package,
      call. = FALSE
    )
  }
}

list_raster_files <- function(spec_dir) {
  files <- list.files(spec_dir, full.names = TRUE)
  files[!grepl(SIDECAR_PATTERN, files, ignore.case = TRUE)]
}

get_quad_id <- function(file, scale) {
  quad_id <- basename(file)

  if (scale == "20m") {
    return(regmatches(quad_id, regexpr("^\\d+", quad_id)))
  }

  quad_id
}

excluded_quad_ids_for_scale <- function(scale) {
  if (scale == "20m") {
    return(EXCLUDED_20M_ATMOSPHERIC)
  }

  if (scale == "10m") {
    return(as.vector(outer(
      EXCLUDED_20M_ATMOSPHERIC,
      c("_a", "_b", "_c", "_d"),
      paste0
    )))
  }

  if (scale == "50m") {
    return(EXCLUDED_50M_ATMOSPHERIC)
  }

  character()
}

is_manual_excluded <- function(quad_id, scale) {
  as.character(quad_id) %in% excluded_quad_ids_for_scale(scale)
}

seed_from_id <- function(id, offset = 0) {
  RANDOM_SEED + offset + sum(utf8ToInt(as.character(id)))
}

apply_shadow_mask <- function(raster) {
  band_mask <- raster[[BEST_BAND]]
  sunlit_mask <- if (DIRECTION == ">") {
    band_mask < SHADOW_THRESHOLD
  } else {
    band_mask > SHADOW_THRESHOLD
  }

  terra::mask(raster, sunlit_mask, maskvalues = 0)
}

clean_spectra <- function(x) {
  if (nrow(x) == 0) {
    return(x)
  }

  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]
  x
}

read_masked_spectra <- function(file) {
  raster <- terra::rast(file)
  raster_masked <- apply_shadow_mask(raster)
  clean_spectra(terra::values(raster_masked, mat = TRUE))
}

sample_spectra_for_pca <- function(file, scale, sample_n = PCA_SAMPLE_PER_RASTER) {
  quad_id <- get_quad_id(file, scale)
  x <- read_masked_spectra(file)

  if (nrow(x) == 0) {
    return(x)
  }

  if (nrow(x) > sample_n) {
    set.seed(seed_from_id(paste0(scale, "_", quad_id), offset = 1000))
    x <- x[sample.int(nrow(x), sample_n), , drop = FALSE]
  }

  x
}

build_or_load_global_pca <- function(scale_config = SCALE_CONFIG) {
  require_runtime_package("terra")

  rebuild <- identical(tolower(Sys.getenv("PCA_REBUILD", "false")), "true")
  if (file.exists(PCA_RDS) && !rebuild) {
    pca_object <- readRDS(PCA_RDS)
    if (identical(pca_object$exclusion_policy_id, PCA_EXCLUSION_POLICY_ID)) {
      log_progress(paste("Loading existing global PCA:", PCA_RDS))
      return(pca_object)
    }
    log_progress(paste(
      "Existing global PCA does not match exclusion policy",
      PCA_EXCLUSION_POLICY_ID,
      "and will be rebuilt"
    ))
  }

  raster_index <- do.call(rbind, lapply(scale_config, function(config) {
    files <- list_raster_files(config$spec_dir)
    data.frame(
      scale = config$scale,
      file = files,
      quad_id = vapply(files, get_quad_id, character(1), scale = config$scale),
      stringsAsFactors = FALSE
    )
  }))
  raster_index$manual_excluded <- mapply(
    is_manual_excluded,
    raster_index$quad_id,
    raster_index$scale
  )
  excluded_summary <- stats::aggregate(
    manual_excluded ~ scale,
    raster_index,
    sum
  )
  names(excluded_summary)[2] <- "excluded_rasters"
  log_progress(paste(
    "Manual PCA exclusions by scale:",
    paste(paste(excluded_summary$scale, excluded_summary$excluded_rasters, sep = "="), collapse = ", ")
  ))
  raster_index <- raster_index[!raster_index$manual_excluded, , drop = FALSE]

  effective_sample_per_raster <- PCA_SAMPLE_PER_RASTER
  if (is.finite(PCA_MAX_ROWS) && PCA_MAX_ROWS > 0 && nrow(raster_index) * effective_sample_per_raster > PCA_MAX_ROWS) {
    effective_sample_per_raster <- max(1, floor(PCA_MAX_ROWS / nrow(raster_index)))
  }

  log_progress(paste(
    "Building global PCA sample from",
    nrow(raster_index),
    "rasters with up to",
    effective_sample_per_raster,
    "pixels per eligible raster"
  ))

  x_list <- vector("list", nrow(raster_index))
  for (i in seq_len(nrow(raster_index))) {
    if (i %% 100 == 0) {
      log_progress(paste("PCA sampling raster", i, "of", nrow(raster_index)))
    }
    x_list[[i]] <- sample_spectra_for_pca(
      raster_index$file[i],
      raster_index$scale[i],
      sample_n = effective_sample_per_raster
    )
  }

  x_global <- do.call(rbind, x_list)

  log_progress(paste(
    "Running prcomp on",
    nrow(x_global),
    "pixels and",
    ncol(x_global),
    "bands"
  ))

  x_scaled <- scale(x_global)
  pca <- stats::prcomp(x_scaled, center = FALSE, scale. = FALSE)
  center <- attr(x_scaled, "scaled:center")
  scale_values <- attr(x_scaled, "scaled:scale")
  scale_values[!is.finite(scale_values) | scale_values == 0] <- 1

  pca_object <- list(
    pca = pca,
    center = center,
    scale = scale_values,
    n_components = PCA_N_COMPONENTS,
    n_sample_rows = nrow(x_global),
    n_bands = ncol(x_global),
    requested_sample_per_raster = PCA_SAMPLE_PER_RASTER,
    effective_sample_per_raster = effective_sample_per_raster,
    max_rows = PCA_MAX_ROWS,
    eligible_rasters = nrow(raster_index),
    excluded_rasters_by_scale = excluded_summary,
    exclusion_policy_id = PCA_EXCLUSION_POLICY_ID,
    excluded_20m_atmospheric = EXCLUDED_20M_ATMOSPHERIC,
    excluded_50m_atmospheric = EXCLUDED_50M_ATMOSPHERIC,
    best_band = BEST_BAND,
    shadow_threshold = SHADOW_THRESHOLD,
    direction = DIRECTION,
    sampling_note = "Uniform per-eligible-quadrat sample within each scale after manual atmospheric/cloud exclusions; final metrics use all retained pixels except excluded quadrats and hull metrics above configured point caps.",
    created_at = Sys.time()
  )

  saveRDS(pca_object, PCA_RDS)
  log_progress(paste("Saved global PCA:", PCA_RDS))
  pca_object
}

write_pca_diagnostics <- function(pca_object) {
  variance <- pca_object$pca$sdev^2
  variance_df <- data.frame(
    pc_axis = seq_along(variance),
    variance = variance,
    pct_variance = variance / sum(variance),
    cumulative_pct_variance = cumsum(variance) / sum(variance),
    stringsAsFactors = FALSE
  )

  write.csv(variance_df, PCA_VARIANCE_CSV, row.names = FALSE)

  plot_n <- min(30, nrow(variance_df))
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  png(PCA_VARIANCE_PNG, width = 1800, height = 1100, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  bar_positions <- barplot(
    height = variance_df$pct_variance[seq_len(plot_n)] * 100,
    names.arg = paste0("PC", seq_len(plot_n)),
    las = 2,
    ylim = c(0, 100),
    col = "#4C78A8",
    border = NA,
    main = "Global PCA Variance Explained",
    ylab = "Variance explained (%)",
    xlab = "Principal component axis"
  )
  lines(
    x = bar_positions,
    y = variance_df$cumulative_pct_variance[seq_len(plot_n)] * 100,
    type = "b",
    pch = 16,
    col = "#D55E00",
    lwd = 2
  )
  legend(
    "topright",
    legend = c("Axis variance", "Cumulative variance"),
    fill = c("#4C78A8", NA),
    border = NA,
    lty = c(NA, 1),
    pch = c(NA, 16),
    col = c("#4C78A8", "#D55E00"),
    bty = "n"
  )

  invisible(variance_df)
}

project_scores <- function(x, pca_object) {
  x_scaled <- sweep(x, 2, pca_object$center, "-")
  x_scaled <- sweep(x_scaled, 2, pca_object$scale, "/")
  x_scaled %*% pca_object$pca$rotation[, seq_len(PCA_N_COMPONENTS), drop = FALSE]
}

convex_hull_area_2d <- function(points) {
  points <- unique(points)
  if (nrow(points) < 3) {
    return(NA_real_)
  }

  hull_index <- grDevices::chull(points[, 1], points[, 2])
  hull <- points[c(hull_index, hull_index[1]), , drop = FALSE]
  abs(sum(hull[-1, 1] * hull[-nrow(hull), 2] - hull[-nrow(hull), 1] * hull[-1, 2])) / 2
}

alpha_hull_area_2d <- function(points, quad_id) {
  points <- unique(points)
  n_available <- nrow(points)

  if (n_available < 4) {
    return(data.frame(
      alpha_hull_area = NA_real_,
      alpha_hull_method = "insufficient_points",
      alpha_hull_n_points = n_available,
      alpha_hull_alpha = ALPHA_HULL_ALPHA,
      pca_convex_hull_area = convex_hull_area_2d(points),
      stringsAsFactors = FALSE
    ))
  }

  method <- "all_pixels"
  points_used <- points
  if (n_available > ALPHA_MAX_POINTS) {
    set.seed(seed_from_id(quad_id, offset = 2000))
    points_used <- points[sample.int(n_available, ALPHA_MAX_POINTS), , drop = FALSE]
    method <- "sampled_pixels"
  }

  area <- tryCatch({
    hull <- alphahull::ahull(points_used[, 1], points_used[, 2], alpha = ALPHA_HULL_ALPHA)
    alphahull::areaahull(hull, timeout = ALPHA_TIMEOUT)
  }, error = function(e) {
    NA_real_
  })

  if (!is.finite(area)) {
    fallback_n <- min(ALPHA_FALLBACK_POINTS, n_available)
    if (fallback_n >= 4) {
      set.seed(seed_from_id(quad_id, offset = 2500))
      points_used <- points[sample.int(n_available, fallback_n), , drop = FALSE]
      area <- tryCatch({
        hull <- alphahull::ahull(points_used[, 1], points_used[, 2], alpha = ALPHA_HULL_ALPHA)
        alphahull::areaahull(hull, timeout = ALPHA_TIMEOUT)
      }, error = function(e) {
        NA_real_
      })
      method <- "fallback_sampled_pixels"
    }
  }

  if (!is.finite(area)) {
    method <- paste(method, "alpha_failed", sep = "_")
  }

  data.frame(
    alpha_hull_area = area,
    alpha_hull_method = method,
    alpha_hull_n_points = nrow(points_used),
    alpha_hull_alpha = ALPHA_HULL_ALPHA,
    pca_convex_hull_area = convex_hull_area_2d(points_used),
    stringsAsFactors = FALSE
  )
}

pca_hull_volume_3d <- function(points, quad_id) {
  points <- unique(points)
  n_available <- nrow(points)

  if (n_available < 4) {
    return(data.frame(
      pca_hull_volume_3d = NA_real_,
      pca_hull_area_3d = NA_real_,
      pca_hull_3d_method = "insufficient_points",
      pca_hull_3d_n_points = n_available,
      stringsAsFactors = FALSE
    ))
  }

  method <- "all_pixels"
  if (n_available > HULL_3D_MAX_POINTS) {
    set.seed(seed_from_id(quad_id, offset = 3000))
    points <- points[sample.int(n_available, HULL_3D_MAX_POINTS), , drop = FALSE]
    method <- "sampled_pixels"
  }

  hull <- tryCatch({
    geometry::convhulln(points, options = "FA")
  }, error = function(e) {
    NULL
  })

  if (is.null(hull)) {
    return(data.frame(
      pca_hull_volume_3d = NA_real_,
      pca_hull_area_3d = NA_real_,
      pca_hull_3d_method = paste(method, "hull_failed", sep = "_"),
      pca_hull_3d_n_points = nrow(points),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    pca_hull_volume_3d = hull$vol,
    pca_hull_area_3d = hull$area,
    pca_hull_3d_method = method,
    pca_hull_3d_n_points = nrow(points),
    stringsAsFactors = FALSE
  )
}

calculate_metrics_from_scores <- function(scores, quad_id) {
  n_pixels <- nrow(scores)

  if (n_pixels < MIN_VALID_PIXELS) {
    return(data.frame(
      quad_id = quad_id,
      n_pixels = n_pixels,
      manual_excluded = FALSE,
      metric_method = "insufficient_pixels",
      pca_euclidean_mean = NA_real_,
      pca_euclidean_median = NA_real_,
      pca_euclidean_sd = NA_real_,
      rao_q_pca = NA_real_,
      rao_q_distance = "squared_euclidean_pc1_pc3",
      alpha_hull_area = NA_real_,
      alpha_hull_method = "insufficient_pixels",
      alpha_hull_n_points = n_pixels,
      alpha_hull_alpha = ALPHA_HULL_ALPHA,
      pca_convex_hull_area = NA_real_,
      pca_hull_volume_3d = NA_real_,
      pca_hull_area_3d = NA_real_,
      pca_hull_3d_method = "insufficient_pixels",
      pca_hull_3d_n_points = n_pixels,
      stringsAsFactors = FALSE
    ))
  }

  center <- colMeans(scores)
  centered <- sweep(scores, 2, center, "-")
  squared_radius <- rowSums(centered^2)
  distances <- sqrt(squared_radius)
  alpha_summary <- alpha_hull_area_2d(scores[, 1:2, drop = FALSE], quad_id)
  hull_3d_summary <- pca_hull_volume_3d(scores[, 1:3, drop = FALSE], quad_id)

  cbind(
    data.frame(
      quad_id = quad_id,
      n_pixels = n_pixels,
      manual_excluded = FALSE,
      metric_method = "all_pixels",
      pca_euclidean_mean = mean(distances),
      pca_euclidean_median = stats::median(distances),
      pca_euclidean_sd = stats::sd(distances),
      rao_q_pca = 2 * mean(squared_radius),
      rao_q_distance = "squared_euclidean_pc1_pc3",
      stringsAsFactors = FALSE
    ),
    alpha_summary,
    hull_3d_summary
  )
}

manual_excluded_metrics <- function(quad_id) {
  data.frame(
    quad_id = quad_id,
    n_pixels = NA_integer_,
    manual_excluded = TRUE,
    metric_method = "manual_excluded",
    pca_euclidean_mean = NA_real_,
    pca_euclidean_median = NA_real_,
    pca_euclidean_sd = NA_real_,
    rao_q_pca = NA_real_,
    rao_q_distance = "squared_euclidean_pc1_pc3",
    alpha_hull_area = NA_real_,
    alpha_hull_method = "manual_excluded",
    alpha_hull_n_points = NA_integer_,
    alpha_hull_alpha = ALPHA_HULL_ALPHA,
    pca_convex_hull_area = NA_real_,
    pca_hull_volume_3d = NA_real_,
    pca_hull_area_3d = NA_real_,
    pca_hull_3d_method = "manual_excluded",
    pca_hull_3d_n_points = NA_integer_,
    stringsAsFactors = FALSE
  )
}

process_raster_metrics <- function(file, scale, pca_object) {
  quad_id <- get_quad_id(file, scale)

  if (is_manual_excluded(quad_id, scale)) {
    return(manual_excluded_metrics(quad_id))
  }

  x <- read_masked_spectra(file)

  if (nrow(x) == 0) {
    return(calculate_metrics_from_scores(matrix(numeric(), ncol = PCA_N_COMPONENTS), quad_id))
  }

  scores <- project_scores(x, pca_object)
  calculate_metrics_from_scores(scores, quad_id)
}

prefix_sa_columns <- function(sa_df) {
  names(sa_df) <- ifelse(names(sa_df) == "quad_id", "quad_id", paste0("sa_", names(sa_df)))
  names(sa_df)[names(sa_df) == "sa_spectral_entropy"] <- "sa_entropy"
  sa_df
}

write_metric_csvs <- function(scale, metrics_df, combined_df) {
  pca_csv <- file.path(
    OUTPUT_DIR,
    paste0(scale, "_PCA_metrics_smooth_masked_5nm_summary.csv")
  )
  combined_csv_root <- file.path(
    OUTPUT_DIR,
    paste0(scale, "_spectral_heterogeneity_smooth_masked_5nm_summary.csv")
  )
  combined_csv_shp_dir <- file.path(
    OUTPUT_SHP_DIR,
    paste0("spectral_heterogeneity_", scale, "_smooth_masked_5nm_summary.csv")
  )

  write.csv(metrics_df, pca_csv, row.names = FALSE)
  write.csv(combined_df, combined_csv_root, row.names = FALSE)
  write.csv(combined_df, combined_csv_shp_dir, row.names = FALSE)

  list(
    pca_csv = pca_csv,
    combined_csv_root = combined_csv_root,
    combined_csv_shp_dir = combined_csv_shp_dir
  )
}

attach_combined_to_quads <- function(shp_path, join_field, combined_df, out_shp) {
  quads <- terra::vect(shp_path)
  join_values <- as.character(quads[[join_field]][, 1])
  result_index <- match(join_values, as.character(combined_df$quad_id))

  quads$quad_id <- join_values
  quads$sa_ent <- combined_df$sa_entropy[result_index]
  quads$pca_eucl <- combined_df$pca_euclidean_mean[result_index]
  quads$alpha_h <- combined_df$alpha_hull_area[result_index]
  quads$hull3d_v <- combined_df$pca_hull_volume_3d[result_index]
  quads$hull3d_a <- combined_df$pca_hull_area_3d[result_index]
  quads$rao_q <- combined_df$rao_q_pca[result_index]
  quads$pix_n <- combined_df$n_pixels[result_index]
  quads$pca_excl <- combined_df$manual_excluded[result_index]
  quads$sa_bcv <- combined_df$sa_boot_cv[result_index]
  quads$sa_meth <- combined_df$sa_method[result_index]
  quads$ah_meth <- combined_df$alpha_hull_method[result_index]
  quads$h3d_meth <- combined_df$pca_hull_3d_method[result_index]

  terra::writeVector(quads, out_shp, overwrite = TRUE)
}

process_scale <- function(config, pca_object) {
  require_runtime_package("terra")
  require_runtime_package("alphahull")
  require_runtime_package("geometry")

  scale <- config$scale
  raster_files <- list_raster_files(config$spec_dir)
  quad_ids <- vapply(raster_files, get_quad_id, character(1), scale = scale)
  n_manual_excluded <- sum(vapply(quad_ids, is_manual_excluded, logical(1), scale = scale))
  n_cores <- min(config$n_cores, max(1, parallel::detectCores() - 3))

  log_progress(paste(
    "Starting PCA/Rao/alpha metrics for",
    scale,
    "using",
    length(raster_files),
    "rasters,",
    n_manual_excluded,
    "manual exclusions, and",
    n_cores,
    "cores"
  ))

  cl <- parallel::makeCluster(n_cores, type = "SOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    library(terra)
    library(alphahull)
    library(geometry)
  })
  parallel::clusterExport(
    cl,
    c(
      "BEST_BAND", "SHADOW_THRESHOLD", "DIRECTION",
      "PCA_N_COMPONENTS", "MIN_VALID_PIXELS", "ALPHA_HULL_ALPHA",
      "ALPHA_MAX_POINTS", "ALPHA_FALLBACK_POINTS", "ALPHA_TIMEOUT",
      "HULL_3D_MAX_POINTS", "RANDOM_SEED",
      "EXCLUDED_20M_ATMOSPHERIC", "EXCLUDED_50M_ATMOSPHERIC",
      "get_quad_id", "excluded_quad_ids_for_scale", "is_manual_excluded",
      "seed_from_id", "apply_shadow_mask", "clean_spectra",
      "read_masked_spectra", "project_scores", "convex_hull_area_2d",
      "alpha_hull_area_2d", "pca_hull_volume_3d", "calculate_metrics_from_scores",
      "manual_excluded_metrics", "process_raster_metrics", "pca_object"
    ),
    envir = environment()
  )

  result_list <- parallel::parLapply(
    cl,
    raster_files,
    function(file) process_raster_metrics(file, scale, pca_object)
  )

  metrics_df <- do.call(rbind, result_list)
  metrics_df <- metrics_df[order(metrics_df$quad_id), , drop = FALSE]

  sa_csv <- file.path(
    OUTPUT_DIR,
    paste0(scale, "_SA_entropy_smooth_masked_5nm_summary.csv")
  )
  if (!file.exists(sa_csv)) {
    stop("Missing SA entropy summary needed for combined output: ", sa_csv)
  }

  sa_df <- read.csv(sa_csv, stringsAsFactors = FALSE)
  sa_df$quad_id <- as.character(sa_df$quad_id)
  sa_df <- prefix_sa_columns(sa_df)

  combined_df <- merge(metrics_df, sa_df, by = "quad_id", all = TRUE)
  combined_df$scale <- scale
  combined_df <- combined_df[order(combined_df$quad_id), , drop = FALSE]

  csv_outputs <- write_metric_csvs(scale, metrics_df, combined_df)
  out_shp <- file.path(
    OUTPUT_SHP_DIR,
    paste0("spectral_heterogeneity_", scale, "_smooth_masked_5nm.shp")
  )
  attach_combined_to_quads(config$shp_path, config$join_field, combined_df, out_shp)

  log_progress(paste(
    "Finished",
    scale,
    "outputs:",
    csv_outputs$combined_csv_shp_dir,
    out_shp
  ))

  c(csv_outputs, list(out_shp = out_shp))
}

run_spectral_heterogeneity_all_metrics <- function(scale_config = SCALE_CONFIG) {
  require_runtime_package("terra")
  require_runtime_package("alphahull")
  require_runtime_package("geometry")

  setwd(PROJECT_DIR)
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUTPUT_SHP_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(PROGRESS_LOG), recursive = TRUE, showWarnings = FALSE)

  pca_object <- build_or_load_global_pca(scale_config)
  write_pca_diagnostics(pca_object)
  outputs <- lapply(scale_config, process_scale, pca_object = pca_object)

  if (requireNamespace("beepr", quietly = TRUE)) {
    beepr::beep(3)
  }

  outputs
}

if (!identical(Sys.getenv("RUN_SPECTRAL_HET_WORKFLOW"), "false")) {
  run_spectral_heterogeneity_all_metrics()
}
