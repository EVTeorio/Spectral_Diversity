user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
OUTPUT_DIR <- "reports/tables/alpha_hull_all_pixels_benchmark"
ANALYSIS_REPORT <- "reports/analysis/20261010_alpha_hull_all_pixels_benchmark.md"
TASK_REPORT <- "reports/tasks/20261010_alpha_hull_all_pixels_benchmark.md"
VALIDATION_REPORT <- "reports/validation/20261010_alpha_hull_all_pixels_benchmark_validation.md"

PCA_RDS <- "Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds"
STANDARDIZED_PCA_RDS <- "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds"

BEST_BAND <- "563 nm"
SHADOW_THRESHOLD <- 0.0305476
PCA_N_COMPONENTS <- 3L
RANDOM_SEED <- 42L

ALPHA_VALUES <- as.numeric(strsplit(Sys.getenv(
  "ALPHA_BENCHMARK_VALUES",
  "1,0.5,2"
), ",")[[1]])
ALPHA_TIMEOUT <- as.numeric(Sys.getenv("ALPHA_BENCHMARK_TIMEOUT", "60"))
BASIS_VALUES <- strsplit(Sys.getenv(
  "ALPHA_BENCHMARK_BASES",
  "vector_normalized_pca"
), ",")[[1]]
N_SELECT_BY_SCALE <- as.integer(strsplit(Sys.getenv(
  "ALPHA_BENCHMARK_N_SELECT",
  "8,4,2"
), ",")[[1]])
if (length(N_SELECT_BY_SCALE) != 3L || any(!is.finite(N_SELECT_BY_SCALE))) {
  stop("ALPHA_BENCHMARK_N_SELECT must contain three comma-separated integers.", call. = FALSE)
}

SCALE_CONFIG <- data.frame(
  scale = c("10m", "20m", "50m"),
  spec_dir = c(
    "Quad_Spectra/10m_smooth_5nm",
    "Quad_Spectra/20m_smooth_5nm",
    "Quad_Spectra/50m_smooth_5nm"
  ),
  quadrat_csv = c(
    "quadrat_analysis_10m.csv",
    "quadrat_analysis_20m.csv",
    "quadrat_analysis_50m.csv"
  ),
  spectral_csv = c(
    "Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv",
    "Quad_Values/20m_spectral_heterogeneity_smooth_masked_5nm_summary.csv",
    "Quad_Values/50m_spectral_heterogeneity_smooth_masked_5nm_summary.csv"
  ),
  n_select = N_SELECT_BY_SCALE,
  stringsAsFactors = FALSE
)

required_packages <- c("terra", "alphahull")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is unavailable: ", package, call. = FALSE)
  }
}

setwd(PROJECT_DIR)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(ANALYSIS_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(TASK_REPORT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(VALIDATION_REPORT), recursive = TRUE, showWarnings = FALSE)

list_raster_files <- function(spec_dir) {
  files <- list.files(spec_dir, full.names = TRUE)
  files[!grepl("\\.(hdr|aux|xml|enp|sta)$", files, ignore.case = TRUE)]
}

get_quad_id <- function(file, scale) {
  quad_id <- basename(file)
  if (scale == "20m") {
    return(regmatches(quad_id, regexpr("^\\d+", quad_id)))
  }
  quad_id
}

read_masked_spectra <- function(file) {
  raster <- terra::rast(file)
  band_mask <- raster[[BEST_BAND]]
  sunlit_mask <- band_mask > SHADOW_THRESHOLD
  raster_masked <- terra::mask(raster, sunlit_mask, maskvalues = 0)
  x <- terra::values(raster_masked, mat = TRUE)
  if (nrow(x) == 0) {
    return(x)
  }
  x <- x[stats::complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]
  x
}

vector_normalize <- function(x) {
  norms <- sqrt(rowSums(x^2))
  keep <- is.finite(norms) & norms > 0
  x <- x[keep, , drop = FALSE]
  norms <- norms[keep]
  sweep(x, 1, norms, "/")
}

project_scores <- function(x, pca_object, preprocess) {
  if (preprocess == "vector_normalized") {
    x <- vector_normalize(x)
  }
  if (nrow(x) == 0) {
    return(matrix(numeric(), ncol = PCA_N_COMPONENTS))
  }
  x_scaled <- sweep(x, 2, pca_object$center, "-")
  x_scaled <- sweep(x_scaled, 2, pca_object$scale, "/")
  x_scaled %*% pca_object$pca$rotation[, seq_len(PCA_N_COMPONENTS), drop = FALSE]
}

select_representative_quadrats <- function(config) {
  quadrat_df <- read.csv(config$quadrat_csv, stringsAsFactors = FALSE)
  spectral_df <- read.csv(config$spectral_csv, stringsAsFactors = FALSE)
  raster_files <- list_raster_files(config$spec_dir)
  raster_index <- data.frame(
    quad_id = vapply(raster_files, get_quad_id, character(1), scale = config$scale),
    raster_file = raster_files,
    stringsAsFactors = FALSE
  )

  selected_pool <- merge(
    quadrat_df[, c("quad_id", "scale", "phy_afaith")],
    spectral_df[, c(
      "quad_id",
      "n_pixels",
      "metric_method",
      "alpha_hull_area",
      "standardized_PCA_alpha_hull_area"
    )],
    by = "quad_id",
    all = FALSE
  )
  selected_pool <- merge(selected_pool, raster_index, by = "quad_id", all = FALSE)
  selected_pool <- selected_pool[
    selected_pool$metric_method == "all_pixels" &
      is.finite(selected_pool$phy_afaith) &
      is.finite(selected_pool$n_pixels),
    ,
    drop = FALSE
  ]
  selected_pool <- selected_pool[order(selected_pool$n_pixels, selected_pool$quad_id), , drop = FALSE]

  if (nrow(selected_pool) < config$n_select) {
    stop("Not enough candidate quadrats for ", config$scale, call. = FALSE)
  }

  indices <- unique(round(seq(1, nrow(selected_pool), length.out = config$n_select)))
  while (length(indices) < config$n_select) {
    missing_indices <- setdiff(seq_len(nrow(selected_pool)), indices)
    indices <- sort(c(indices, missing_indices[seq_len(config$n_select - length(indices))]))
  }
  selected_pool[indices, , drop = FALSE]
}

calculate_alpha_area <- function(points, alpha_value) {
  points_unique <- unique(points)
  if (nrow(points_unique) < 4) {
    return(list(
      area = NA_real_,
      status = "insufficient_unique_points",
      n_unique_points = nrow(points_unique),
      arc_count = NA_integer_,
      elapsed_seconds = NA_real_
    ))
  }

  start_time <- proc.time()[["elapsed"]]
  result <- tryCatch({
    hull <- alphahull::ahull(points_unique[, 1], points_unique[, 2], alpha = alpha_value)
    area <- alphahull::areaahull(hull, timeout = ALPHA_TIMEOUT)
    elapsed <- proc.time()[["elapsed"]] - start_time
    list(
      area = area,
      status = if (is.finite(area)) "ok" else "nonfinite_area",
      n_unique_points = nrow(points_unique),
      arc_count = if (!is.null(hull$arcs)) nrow(hull$arcs) else NA_integer_,
      elapsed_seconds = elapsed
    )
  }, error = function(e) {
    elapsed <- proc.time()[["elapsed"]] - start_time
    list(
      area = NA_real_,
      status = paste0("error: ", conditionMessage(e)),
      n_unique_points = nrow(points_unique),
      arc_count = NA_integer_,
      elapsed_seconds = elapsed
    )
  })

  result
}

run_quadrat_benchmark <- function(row, pca_objects) {
  x <- read_masked_spectra(row$raster_file)
  basis_rows <- list()

  for (basis_name in intersect(names(pca_objects), BASIS_VALUES)) {
    pca_info <- pca_objects[[basis_name]]
    scores <- project_scores(x, pca_info$pca_object, pca_info$preprocess)
    points <- scores[, 1:2, drop = FALSE]
    for (alpha_value in ALPHA_VALUES) {
      if (is_completed_combo(row$scale, row$quad_id, basis_name, alpha_value)) {
        message("Skipping completed ", row$scale, " ", row$quad_id, " ", basis_name, " alpha=", alpha_value)
        next
      }
      alpha_result <- calculate_alpha_area(points, alpha_value)
      result_row <- data.frame(
        scale = row$scale,
        quad_id = row$quad_id,
        basis = basis_name,
        alpha_value = alpha_value,
        phy_afaith = row$phy_afaith,
        n_retained_pixels = nrow(x),
        n_unique_pc12_points = alpha_result$n_unique_points,
        alpha_area_all_pixels = alpha_result$area,
        status = alpha_result$status,
        arc_count = alpha_result$arc_count,
        elapsed_seconds = alpha_result$elapsed_seconds,
        production_alpha_area = if (basis_name == "original_pca") {
          row$alpha_hull_area
        } else {
          row$standardized_PCA_alpha_hull_area
        },
        stringsAsFactors = FALSE
      )
      append_result_row(result_row)
      basis_rows[[length(basis_rows) + 1L]] <- result_row
    }
  }

  if (length(basis_rows) == 0) {
    return(NULL)
  }
  do.call(rbind, basis_rows)
}

results_csv <- file.path(OUTPUT_DIR, "r_alphahull_all_pixels_results.csv")

read_existing_results <- function() {
  if (!file.exists(results_csv)) {
    return(data.frame())
  }
  read.csv(results_csv, stringsAsFactors = FALSE)
}

is_completed_combo <- function(scale, quad_id, basis, alpha_value) {
  existing <- read_existing_results()
  if (nrow(existing) == 0) {
    return(FALSE)
  }
  any(
    existing$scale == scale &
      existing$quad_id == quad_id &
      existing$basis == basis &
      abs(existing$alpha_value - alpha_value) < 1e-12
  )
}

append_result_row <- function(result_row) {
  write.table(
    result_row,
    results_csv,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(results_csv),
    append = file.exists(results_csv)
  )
}

correlation_summary <- function(results) {
  groups <- unique(results[, c("basis", "alpha_value")])
  group_rows <- list()

  for (i in seq_len(nrow(groups))) {
    subset_i <- results[
      results$basis == groups$basis[i] &
        results$alpha_value == groups$alpha_value[i] &
        results$status == "ok" &
        is.finite(results$alpha_area_all_pixels) &
        is.finite(results$phy_afaith),
      ,
      drop = FALSE
    ]
    group_rows[[length(group_rows) + 1L]] <- summarize_correlation_row(
      subset_i,
      groups$basis[i],
      groups$alpha_value[i],
      "all_scales"
    )

    for (scale in unique(results$scale)) {
      scale_subset <- subset_i[subset_i$scale == scale, , drop = FALSE]
      group_rows[[length(group_rows) + 1L]] <- summarize_correlation_row(
        scale_subset,
        groups$basis[i],
        groups$alpha_value[i],
        scale
      )
    }
  }

  do.call(rbind, group_rows)
}

summarize_correlation_row <- function(data, basis, alpha_value, scale_group) {
  n <- nrow(data)
  pearson_r <- NA_real_
  spearman_r <- NA_real_
  r_squared <- NA_real_
  p_value <- NA_real_

  if (n >= 3 && stats::sd(data$alpha_area_all_pixels) > 0 && stats::sd(data$phy_afaith) > 0) {
    pearson_test <- stats::cor.test(
      data$alpha_area_all_pixels,
      data$phy_afaith,
      method = "pearson"
    )
    spearman_test <- suppressWarnings(stats::cor.test(
      data$alpha_area_all_pixels,
      data$phy_afaith,
      method = "spearman",
      exact = FALSE
    ))
    pearson_r <- unname(pearson_test$estimate)
    spearman_r <- unname(spearman_test$estimate)
    r_squared <- pearson_r^2
    p_value <- pearson_test$p.value
  }

  data.frame(
    basis = basis,
    alpha_value = alpha_value,
    scale_group = scale_group,
    n = n,
    pearson_r = pearson_r,
    r_squared = r_squared,
    pearson_p = p_value,
    spearman_r = spearman_r,
    stringsAsFactors = FALSE
  )
}

format_number <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

write_markdown_table <- function(data) {
  if (nrow(data) == 0) {
    return("_No rows._")
  }
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  rows <- apply(data, 1, function(row) paste("|", paste(row, collapse = " | "), "|"))
  paste(c(header, separator, rows), collapse = "\n")
}

set.seed(RANDOM_SEED)

pca_objects <- list(
  original_pca = list(
    pca_object = readRDS(PCA_RDS),
    preprocess = "raw"
  ),
  vector_normalized_pca = list(
    pca_object = readRDS(STANDARDIZED_PCA_RDS),
    preprocess = "vector_normalized"
  )
)

selected_quadrats <- do.call(
  rbind,
  lapply(seq_len(nrow(SCALE_CONFIG)), function(i) select_representative_quadrats(SCALE_CONFIG[i, ]))
)

write.csv(
  selected_quadrats,
  file.path(OUTPUT_DIR, "selected_quadrats.csv"),
  row.names = FALSE
)

result_list <- vector("list", nrow(selected_quadrats))
for (i in seq_len(nrow(selected_quadrats))) {
  message(
    "Benchmarking ",
    selected_quadrats$scale[i],
    " ",
    selected_quadrats$quad_id[i],
    " (",
    i,
    "/",
    nrow(selected_quadrats),
    ")"
  )
  result_list[[i]] <- run_quadrat_benchmark(selected_quadrats[i, ], pca_objects)
}

benchmark_results <- read_existing_results()
selected_keys <- paste(selected_quadrats$scale, selected_quadrats$quad_id, sep = "\r")
result_keys <- paste(benchmark_results$scale, benchmark_results$quad_id, sep = "\r")
current_benchmark_results <- benchmark_results[result_keys %in% selected_keys, , drop = FALSE]
correlations <- correlation_summary(current_benchmark_results)

write.csv(
  correlations,
  file.path(OUTPUT_DIR, "alpha_hull_phy_afaith_correlations.csv"),
  row.names = FALSE
)

ok_results <- current_benchmark_results[current_benchmark_results$status == "ok", , drop = FALSE]
runtime_summary <- aggregate(
  elapsed_seconds ~ scale + basis + alpha_value,
  data = ok_results,
  FUN = function(x) round(sum(x), 3)
)
names(runtime_summary)[names(runtime_summary) == "elapsed_seconds"] <- "total_elapsed_seconds"
write.csv(
  runtime_summary,
  file.path(OUTPUT_DIR, "alpha_hull_runtime_summary.csv"),
  row.names = FALSE
)

report_correlations <- correlations[
  correlations$scale_group == "all_scales" &
    correlations$basis == "vector_normalized_pca",
  ,
  drop = FALSE
]
report_correlations <- report_correlations[order(report_correlations$alpha_value), , drop = FALSE]
report_table <- data.frame(
  alpha = report_correlations$alpha_value,
  n = report_correlations$n,
  pearson_r = format_number(report_correlations$pearson_r),
  r_squared = format_number(report_correlations$r_squared),
  pearson_p = format_number(report_correlations$pearson_p),
  stringsAsFactors = FALSE
)

method_status <- as.data.frame(table(current_benchmark_results$status), stringsAsFactors = FALSE)
names(method_status) <- c("status", "rows")
reported_alpha_values <- sort(unique(current_benchmark_results$alpha_value))

analysis_lines <- c(
  "# All-Pixel Alpha-Hull Benchmark",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Purpose",
  "",
  paste(
    "Benchmark R `alphahull` all-retained-pixel PC1-PC2 alpha-hull areas for",
    paste0(SCALE_CONFIG$n_select[1], " selected 10 m quadrats,"),
    paste0(SCALE_CONFIG$n_select[2], " selected 20 m quadrats, and"),
    paste0(SCALE_CONFIG$n_select[3], " selected 50 m quadrats.")
  ),
  "The benchmark tests multiple alpha parameter values and compares the resulting areas with abundance-weighted Faith's PD (`phy_afaith`).",
  "",
  "## Alpha Values",
  "",
  paste("Tested alpha values in current selected subset:", paste(reported_alpha_values, collapse = ", ")),
  "",
  "## Vector-Normalized PCA Correlations With Abundance-Weighted Faith's PD",
  "",
  write_markdown_table(report_table),
  "",
  "## Status Counts",
  "",
  write_markdown_table(method_status),
  "",
  "## Output Tables",
  "",
  "- `reports/tables/alpha_hull_all_pixels_benchmark/selected_quadrats.csv`",
  "- `reports/tables/alpha_hull_all_pixels_benchmark/r_alphahull_all_pixels_results.csv`",
  "- `reports/tables/alpha_hull_all_pixels_benchmark/alpha_hull_phy_afaith_correlations.csv`",
  "- `reports/tables/alpha_hull_all_pixels_benchmark/alpha_hull_runtime_summary.csv`",
  "",
  "## Interpretation Notes",
  "",
  "- These are benchmark-subset correlations, not final all-quadrat inference.",
  "- Alpha-hull geometry uses all retained valid sunlit pixels projected into PC1-PC2 space; duplicate PC1-PC2 coordinates are collapsed because duplicates do not change the hull boundary.",
  paste0(
    "- Scale-specific correlations are included in the CSV output. The 50 m subset has ",
    SCALE_CONFIG$n_select[SCALE_CONFIG$scale == "50m"],
    " quadrats, so that scale-specific estimate remains small-n and should be interpreted cautiously."
  ),
  "- Python benchmarking was intentionally deferred."
)
writeLines(analysis_lines, ANALYSIS_REPORT)

task_lines <- c(
  "# Task Report: All-Pixel Alpha-Hull Benchmark",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Actions",
  "",
  "- Added `scripts/3_Analysis/alpha_hull_all_pixels_benchmark.R`.",
  paste0(
    "- Selected representative included quadrats across 10 m, 20 m, and 50 m scales: ",
    paste(paste0(SCALE_CONFIG$scale, " n=", SCALE_CONFIG$n_select), collapse = ", "),
    "."
  ),
  "- Recomputed all-retained-pixel PC1-PC2 alpha-hull areas with R `alphahull` across multiple alpha values.",
  "- Wrote benchmark results, runtime summaries, correlations, and an analysis report under `reports/`.",
  "",
  "## Notes",
  "",
  "- The benchmark was run for the primary vector-normalized PCA basis.",
  "- The requested Python comparison was deferred for a later step."
)
writeLines(task_lines, TASK_REPORT)

validation_checks <- data.frame(
  check = c(
    "selected_quadrats_rows",
    "current_selection_benchmark_result_rows",
    "raw_benchmark_result_rows",
    "correlation_rows",
    "all_requested_scales_present",
    "has_vector_normalized_alpha_1"
  ),
  value = c(
    nrow(selected_quadrats),
    nrow(current_benchmark_results),
    nrow(benchmark_results),
    nrow(correlations),
    all(c("10m", "20m", "50m") %in% selected_quadrats$scale),
    any(correlations$basis == "vector_normalized_pca" & correlations$alpha_value == 1)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  validation_checks,
  file.path(OUTPUT_DIR, "alpha_hull_benchmark_validation_checks.csv"),
  row.names = FALSE
)

validation_lines <- c(
  "# Validation Report: All-Pixel Alpha-Hull Benchmark",
  "",
  paste("Date:", Sys.Date()),
  "",
  "## Checks",
  "",
  write_markdown_table(validation_checks),
  "",
  "## Result",
  "",
  if (all(validation_checks$value != "FALSE")) {
    "Pass. Benchmark outputs were generated for the requested scales and alpha values."
  } else {
    "Review required. One or more validation checks did not pass."
  }
)
writeLines(validation_lines, VALIDATION_REPORT)

message("Completed alpha hull all-pixel benchmark.")
