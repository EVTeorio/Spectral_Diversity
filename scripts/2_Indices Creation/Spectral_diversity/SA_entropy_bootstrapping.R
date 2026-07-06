user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

# Spectral angle entropy from the current smoothed and 5 nm resampled
# quadrat spectra. Exact all-retained-pixel entropy is used only when pair
# counts are small enough; larger quadrats fall back to bootstrap mean entropy.

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
SIDECAR_PATTERN <- "\\.(hdr|aux|xml|enp|sta)$"
OUTPUT_DIR <- "Quad_Values"
OUTPUT_SHP_DIR <- "Quad_Values/Spectral_diversitySHPs"
PROGRESS_LOG <- "logs/sa_entropy_bootstrapping_progress.log"

N_BOOT <- as.integer(Sys.getenv("SA_N_BOOT", "70"))
N_SAMPLE <- as.integer(Sys.getenv("SA_N_SAMPLE", "5000"))
N_PAIR_SAMPLE <- as.integer(Sys.getenv("SA_N_PAIR_SAMPLE", "10000"))
N_BINS <- as.integer(Sys.getenv("SA_N_BINS", "40"))
MAX_EXACT_PAIRS <- as.numeric(Sys.getenv("SA_MAX_EXACT_PAIRS", "2000000"))
MAX_BOOT_EXACT_PAIRS <- as.numeric(Sys.getenv("SA_MAX_BOOT_EXACT_PAIRS", "250000"))
MAX_STORED_PAIRS <- as.numeric(Sys.getenv("SA_MAX_STORED_PAIRS", "20000000"))
ANGLE_BLOCK_SIZE <- as.integer(Sys.getenv("SA_ANGLE_BLOCK_SIZE", "250"))
RANDOM_SEED <- as.integer(Sys.getenv("SA_RANDOM_SEED", "42"))

# Current shadow mask from the prior spectral-angle entropy workflow.
# With DIRECTION set to "<", pixels with 563 nm reflectance > 0.0305476
# are retained as sunlit pixels.
BEST_BAND <- "563 nm"
SHADOW_THRESHOLD <- 0.0305476
DIRECTION <- "<"

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

optional_beep <- function(sound = 1) {
  if (requireNamespace("beepr", quietly = TRUE)) {
    beepr::beep(sound)
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

pair_count <- function(n_pixels) {
  n_pixels * (n_pixels - 1) / 2
}

normalize_spectra <- function(x) {
  if (nrow(x) == 0) {
    return(x)
  }

  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]

  if (nrow(x) == 0) {
    return(x)
  }

  norms <- sqrt(rowSums(x^2))
  x <- x[is.finite(norms) & norms > 0, , drop = FALSE]
  norms <- norms[is.finite(norms) & norms > 0]

  x / norms
}

pairwise_angles <- function(x_norm, block_size = ANGLE_BLOCK_SIZE) {
  n_pixels <- nrow(x_norm)
  n_pairs <- pair_count(n_pixels)

  if (n_pixels < 3) {
    return(numeric())
  }

  if (n_pairs > MAX_STORED_PAIRS) {
    stop("Pair count exceeds MAX_STORED_PAIRS: ", n_pairs)
  }

  angles <- numeric(n_pairs)
  write_index <- 1

  for (block_start in seq(1, n_pixels - 1, by = block_size)) {
    block_end <- min(block_start + block_size - 1, n_pixels - 1)
    dots <- x_norm[block_start:block_end, , drop = FALSE] %*% t(x_norm)

    for (row_index in seq_len(nrow(dots))) {
      source_index <- block_start + row_index - 1
      values <- dots[row_index, (source_index + 1):n_pixels]
      values <- pmin(pmax(as.numeric(values), -1), 1)
      next_index <- write_index + length(values) - 1
      angles[write_index:next_index] <- acos(values)
      write_index <- next_index + 1
    }
  }

  angles
}

shannon_entropy <- function(x, n_bins = N_BINS) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  h <- hist(x, breaks = n_bins, plot = FALSE)
  probs <- h$counts / sum(h$counts)
  probs <- probs[probs > 0]
  -sum(probs * log(probs))
}

spectral_entropy_from_norm <- function(x_norm, n_bins = N_BINS) {
  if (nrow(x_norm) < 3) {
    return(NA_real_)
  }

  angles <- pairwise_angles(x_norm)
  shannon_entropy(angles, n_bins = n_bins)
}

sample_pairwise_angles <- function(x_norm, n_pairs = N_PAIR_SAMPLE) {
  n_pixels <- nrow(x_norm)

  if (n_pixels < 3) {
    return(numeric())
  }

  first_index <- sample.int(n_pixels, n_pairs, replace = TRUE)
  second_index <- sample.int(n_pixels, n_pairs, replace = TRUE)
  same_index <- first_index == second_index

  while (any(same_index)) {
    second_index[same_index] <- sample.int(n_pixels, sum(same_index), replace = TRUE)
    same_index <- first_index == second_index
  }

  dots <- rowSums(x_norm[first_index, , drop = FALSE] * x_norm[second_index, , drop = FALSE])
  acos(pmin(pmax(dots, -1), 1))
}

spectral_entropy_from_sampled_pairs <- function(x_norm, n_pairs = N_PAIR_SAMPLE, n_bins = N_BINS) {
  angles <- sample_pairwise_angles(x_norm, n_pairs = n_pairs)
  shannon_entropy(angles, n_bins = n_bins)
}

seed_from_quad_id <- function(quad_id) {
  RANDOM_SEED + sum(utf8ToInt(as.character(quad_id)))
}

bootstrap_entropy <- function(x_norm, quad_id) {
  n_pixels <- nrow(x_norm)

  if (n_pixels < 10) {
    return(rep(NA_real_, N_BOOT))
  }

  set.seed(seed_from_quad_id(quad_id))

  vapply(seq_len(N_BOOT), function(iteration) {
    sample_size <- min(N_SAMPLE, n_pixels)
    sample_index <- sample(n_pixels, sample_size, replace = FALSE)
    sampled_norm <- x_norm[sample_index, , drop = FALSE]

    if (pair_count(nrow(sampled_norm)) <= MAX_BOOT_EXACT_PAIRS) {
      return(spectral_entropy_from_norm(sampled_norm))
    }

    spectral_entropy_from_sampled_pairs(sampled_norm)
  }, numeric(1))
}

summarize_bootstrap <- function(values) {
  if (all(is.na(values))) {
    return(data.frame(
      boot_mean = NA_real_,
      boot_sd = NA_real_,
      boot_cv = NA_real_,
      boot_median = NA_real_,
      boot_min = NA_real_,
      boot_max = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  boot_mean <- mean(values, na.rm = TRUE)
  boot_sd <- sd(values, na.rm = TRUE)

  data.frame(
    boot_mean = boot_mean,
    boot_sd = boot_sd,
    boot_cv = boot_sd / boot_mean,
    boot_median = median(values, na.rm = TRUE),
    boot_min = min(values, na.rm = TRUE),
    boot_max = max(values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
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

process_raster <- function(file, scale) {
  raster <- terra::rast(file)
  quad_id <- get_quad_id(file, scale)

  raster_masked <- apply_shadow_mask(raster)
  spectra <- terra::values(raster_masked, mat = TRUE)
  spectra_norm <- normalize_spectra(spectra)
  n_pixels <- nrow(spectra_norm)
  n_pairs <- pair_count(n_pixels)

  if (n_pixels < 10) {
    return(list(
      summary = data.frame(
        quad_id = quad_id,
        n_pixels = n_pixels,
        n_pairs = n_pairs,
        method = "insufficient_pixels",
        spectral_entropy = NA_real_,
        exact_entropy = NA_real_,
        boot_mean = NA_real_,
        boot_sd = NA_real_,
        boot_cv = NA_real_,
        boot_median = NA_real_,
        boot_min = NA_real_,
        boot_max = NA_real_,
        stringsAsFactors = FALSE
      ),
      boot = data.frame(
        quad_id = quad_id,
        n_pixels = n_pixels,
        boot_iter = seq_len(N_BOOT),
        spectral_entropy = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  if (n_pairs <= MAX_EXACT_PAIRS) {
    exact_entropy <- spectral_entropy_from_norm(spectra_norm)
    return(list(
      summary = data.frame(
        quad_id = quad_id,
        n_pixels = n_pixels,
        n_pairs = n_pairs,
        method = "exact_all_pixels",
        spectral_entropy = exact_entropy,
        exact_entropy = exact_entropy,
        boot_mean = NA_real_,
        boot_sd = NA_real_,
        boot_cv = NA_real_,
        boot_median = NA_real_,
        boot_min = NA_real_,
        boot_max = NA_real_,
        stringsAsFactors = FALSE
      ),
      boot = data.frame(
        quad_id = quad_id,
        n_pixels = n_pixels,
        boot_iter = seq_len(N_BOOT),
        spectral_entropy = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  boot_values <- bootstrap_entropy(spectra_norm, quad_id)
  boot_summary <- summarize_bootstrap(boot_values)

  list(
    summary = cbind(
      data.frame(
        quad_id = quad_id,
        n_pixels = n_pixels,
        n_pairs = n_pairs,
        method = "bootstrap_mean",
        spectral_entropy = boot_summary$boot_mean,
        exact_entropy = NA_real_,
        stringsAsFactors = FALSE
      ),
      boot_summary
    ),
    boot = data.frame(
      quad_id = quad_id,
      n_pixels = n_pixels,
      boot_iter = seq_len(N_BOOT),
      spectral_entropy = boot_values,
      stringsAsFactors = FALSE
    )
  )
}

attach_results_to_quads <- function(shp_path, join_field, summary_df, out_shp) {
  quads <- terra::vect(shp_path)
  join_values <- as.character(quads[[join_field]][, 1])
  result_index <- match(join_values, as.character(summary_df$quad_id))

  quads$spec_ent <- summary_df$spectral_entropy[result_index]
  quads$method <- summary_df$method[result_index]
  quads$pix_n <- summary_df$n_pixels[result_index]
  quads$pair_n <- summary_df$n_pairs[result_index]
  quads$exact_ent <- summary_df$exact_entropy[result_index]
  quads$boot_mean <- summary_df$boot_mean[result_index]
  quads$boot_sd <- summary_df$boot_sd[result_index]
  quads$boot_cv <- summary_df$boot_cv[result_index]

  terra::writeVector(quads, out_shp, overwrite = TRUE)
}

process_scale <- function(config) {
  require_runtime_package("terra")

  scale <- config$scale
  raster_files <- list_raster_files(config$spec_dir)
  n_cores <- min(config$n_cores, max(1, parallel::detectCores() - 3))

  log_progress(paste(
    "Starting spectral angle entropy for",
    scale,
    "using",
    length(raster_files),
    "rasters and",
    n_cores,
    "cores"
  ))

  cl <- parallel::makeCluster(n_cores, type = "SOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    library(terra)
  })
  parallel::clusterExport(
    cl,
    c(
      "BEST_BAND", "SHADOW_THRESHOLD", "DIRECTION",
      "N_BOOT", "N_SAMPLE", "N_PAIR_SAMPLE", "N_BINS",
      "MAX_EXACT_PAIRS", "MAX_BOOT_EXACT_PAIRS", "MAX_STORED_PAIRS", "ANGLE_BLOCK_SIZE",
      "RANDOM_SEED", "get_quad_id", "pair_count", "normalize_spectra",
      "pairwise_angles", "shannon_entropy", "spectral_entropy_from_norm",
      "sample_pairwise_angles", "spectral_entropy_from_sampled_pairs",
      "seed_from_quad_id", "bootstrap_entropy", "summarize_bootstrap",
      "apply_shadow_mask", "process_raster"
    ),
    envir = environment()
  )

  result_list <- parallel::parLapply(
    cl,
    raster_files,
    function(file) process_raster(file, scale)
  )

  summary_df <- do.call(rbind, lapply(result_list, `[[`, "summary"))
  summary_df <- summary_df[order(summary_df$quad_id), , drop = FALSE]

  boot_df_long <- do.call(rbind, lapply(result_list, `[[`, "boot"))
  boot_df_long <- boot_df_long[order(boot_df_long$quad_id, boot_df_long$boot_iter), , drop = FALSE]

  boot_df_wide <- stats::reshape(
    boot_df_long,
    idvar = c("quad_id", "n_pixels"),
    timevar = "boot_iter",
    direction = "wide"
  )
  names(boot_df_wide) <- sub("^spectral_entropy\\.", "boot_", names(boot_df_wide))
  boot_df_wide <- boot_df_wide[order(boot_df_wide$quad_id), , drop = FALSE]

  summary_csv <- file.path(
    OUTPUT_DIR,
    paste0(scale, "_SA_entropy_smooth_masked_5nm_summary.csv")
  )
  boot_long_csv <- file.path(
    OUTPUT_DIR,
    paste0(scale, "_SA_entropy_smooth_masked_5nm_boot_long.csv")
  )
  boot_wide_csv <- file.path(
    OUTPUT_DIR,
    paste0(scale, "_SA_entropy_smooth_masked_5nm_boot_wide.csv")
  )
  out_shp <- file.path(
    OUTPUT_SHP_DIR,
    paste0("SA_entropy_", scale, "_smooth_masked_5nm.shp")
  )

  write.csv(summary_df, summary_csv, row.names = FALSE)
  write.csv(boot_df_long, boot_long_csv, row.names = FALSE)
  write.csv(boot_df_wide, boot_wide_csv, row.names = FALSE)
  attach_results_to_quads(config$shp_path, config$join_field, summary_df, out_shp)

  log_progress(paste(
    "Finished",
    scale,
    "outputs:",
    summary_csv,
    boot_long_csv,
    boot_wide_csv,
    out_shp
  ))

  list(
    summary_csv = summary_csv,
    boot_long_csv = boot_long_csv,
    boot_wide_csv = boot_wide_csv,
    out_shp = out_shp
  )
}

run_sa_entropy_workflow <- function(scale_config = SCALE_CONFIG) {
  require_runtime_package("terra")

  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }
  if (!dir.exists(OUTPUT_SHP_DIR)) {
    dir.create(OUTPUT_SHP_DIR, recursive = TRUE)
  }

  setwd(PROJECT_DIR)
  outputs <- lapply(scale_config, process_scale)
  optional_beep(3)
  outputs
}

if (!identical(Sys.getenv("RUN_SA_ENTROPY_WORKFLOW"), "false")) {
  run_sa_entropy_workflow()
}
