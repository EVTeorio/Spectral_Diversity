# Recover area-evaluation failures using the exact same alphahull boundary.
# Run only after alpha_half_small_scales.R has finished all authorized scales.
.libPaths(c("C:/Users/PaintRock/AppData/Local/R/win-library/4.2", .libPaths()))
source("scripts/3_Analysis/alpha_hull_arc_integral.R")
OUT <- "reports/tables/20261015_alpha_half_small_scales"
CANONICAL <- "Quad_Values/Spectral_diversitySHPs"
SCALES <- c("10m", "20m", "25m")
stopifnot(file.exists(file.path(OUT, "validation_checks.csv")))
basic <- read.csv(file.path(OUT, "validation_checks.csv"))
stopifnot(all(basic$pass))
synthetic <- read.csv(file.path(OUT, "arc_integral_extended_validation.csv"))
actual10 <- read.csv(file.path(OUT, "arc_integral_actual_10m_validation.csv"))
actual20_25 <- read.csv(file.path(OUT, "arc_integral_actual_20m_25m_validation.csv"))
endpoint_synthetic <- read.csv(file.path(OUT, "endpoint_integral_validation.csv"))
endpoint_actual <- read.csv(file.path(OUT, "endpoint_integral_actual_validation.csv"))
stopifnot(nrow(synthetic) == 100, all(synthetic$pass),
          all(actual10$agreement[is.finite(actual10$reference_area)]),
          all(actual20_25$agreement), all(is.finite(actual20_25$reference_area)),
          nrow(endpoint_synthetic) == 100, all(endpoint_synthetic$pass),
          all(is.finite(endpoint_actual$endpoint_area)),
          all(endpoint_actual$absolute_difference[is.finite(endpoint_actual$reference_area)] < 1e-8))
reference_hashes <- read.csv(file.path(OUT, "protected_hashes.csv"))
stopifnot(identical(unname(tools::md5sum(reference_hashes$path)), reference_hashes$md5))
helper <- new.env()
expressions <- parse("scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R")
eval(expressions[-length(expressions)], helper)
models <- list(original_pca = readRDS(helper$PCA_RDS), vector_normalized_pca = readRDS(helper$STANDARDIZED_PCA_RDS))
recovery <- list()
all_results <- list()
backup <- file.path(OUT, "before_area_recovery")
dir.create(backup, recursive = TRUE, showWarnings = FALSE)
for (scale in SCALES) {
  path <- file.path(CANONICAL, paste0("alpha_hull_0_5_", scale, "_all_pixels.csv"))
  d <- read.csv(path, stringsAsFactors = FALSE)
  candidates <- which(!d$status %in% c("ok", "missing_raster", "skipped_manual_excluded"))
  if (length(candidates) && !file.exists(file.path(backup, basename(path)))) {
    stopifnot(file.copy(path, file.path(backup, basename(path))))
  }
  for (i in candidates) {
    row <- d[i, ]
    record <- data.frame(scale = scale, quad_id = row$quad_id, basis = row$basis,
      previous_status = row$status, previous_error = row$error, previous_warning = row$warning,
      recovered_area = NA_real_, convex_area_bound = NA_real_, recovery_status = "not_recovered", recovery_error = "",
      recovery_seconds = NA_real_, recovery_method = "arc_integral", stringsAsFactors = FALSE)
    start <- proc.time()[["elapsed"]]
    tryCatch({
      model <- models[[row$basis]]
      raster <- terra::rast(row$raster_file)
      stopifnot(identical(names(raster), names(model$center)))
      x <- helper$read_masked_spectra(row$raster_file)
      points <- unique(helper$project_scores(x, model)[, 1:2, drop = FALSE])
      stopifnot(nrow(x) == row$n_retained_pixels, nrow(points) == row$n_unique_pc12_points)
      hull <- alphahull::ahull(points[, 1], points[, 2], alpha = 0.5)
      area <- tryCatch(alpha_hull_arc_integral(hull), error = function(e) {
        record$recovery_method <<- "endpoint_integral"
        alpha_hull_endpoint_integral(hull)
      })
      convex_bound <- helper$convex_hull_area_2d(points)
      stopifnot(is.finite(convex_bound), area <= convex_bound + 1e-8 + 1e-10 * abs(convex_bound))
      row$alpha_area_all_pixels <- area
      row$status <- "ok"
      row$arc_count <- nrow(hull$arcs)
      row$source <- paste0("computed_all_pixels_alpha_0_5_", record$recovery_method, "_recovery")
      row$warning <- paste(na.omit(c(row$warning, "Recovered by exact integral of unchanged alphahull circular arcs")), collapse = " | ")
      row$error <- ""
      row$run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      record$recovered_area <- area
      record$convex_area_bound <- convex_bound
      record$recovery_status <- "recovered"
    }, error = function(e) record$recovery_error <<- conditionMessage(e))
    record$recovery_seconds <- proc.time()[["elapsed"]] - start
    recovery[[length(recovery) + 1]] <- record
    if (record$recovery_status == "recovered") {
      row$elapsed_seconds <- row$elapsed_seconds + record$recovery_seconds
      d[i, ] <- row
      checkpoint <- file.path(OUT, scale, paste0(row$quad_id, ".csv"))
      cp <- read.csv(checkpoint, stringsAsFactors = FALSE)
      old_checkpoint <- file.path(backup, paste0(scale, "_", row$quad_id, ".csv"))
      if (!file.exists(old_checkpoint)) stopifnot(file.copy(checkpoint, old_checkpoint))
      cp[cp$basis == row$basis, ] <- row
      write.csv(cp, checkpoint, row.names = FALSE)
    }
    print(record[, c("scale", "quad_id", "basis", "previous_status", "recovery_status", "recovery_error")]); flush.console()
  }
  write.csv(d, path, row.names = FALSE)
  summary <- as.data.frame(table(d$scale, d$basis, d$status))
  names(summary) <- c("scale", "basis", "status", "rows")
  write.csv(summary[summary$rows > 0, ], file.path(OUT, paste0(scale, "_status_summary.csv")), row.names = FALSE)
  all_results[[scale]] <- d
}
if (length(recovery)) {
  audit_path <- file.path(OUT, "area_recovery_audit.csv")
  audit <- do.call(rbind, recovery)
  if (file.exists(audit_path)) {
    previous <- read.csv(audit_path, stringsAsFactors = FALSE)
    if (!"recovery_method" %in% names(previous)) previous$recovery_method <- "arc_integral"
    audit <- rbind(previous[names(audit)], audit)
  }
  write.csv(audit, audit_path, row.names = FALSE)
}
combined <- do.call(rbind, all_results)
unresolved <- combined[!combined$status %in% c("ok", "missing_raster", "skipped_manual_excluded"), ]
write.csv(unresolved, file.path(OUT, "unresolved_geometry.csv"), row.names = FALSE)
coverage <- list()
finite_unchanged <- TRUE
timings <- list()
for (scale in SCALES) {
  d <- all_results[[scale]]
  grid <- foreign::read.dbf(paste0("Quad_Scale_SHPs/PR_", scale, ".dbf"))
  id_field <- if (scale == "10m") "sub_id" else "Name"
  for (basis in names(models)) {
    b <- d[d$basis == basis, ]
    stopifnot(setequal(as.character(b$quad_id), as.character(grid[[id_field]])))
    coverage[[length(coverage) + 1]] <- data.frame(scale = scale, basis = basis,
      grid_rows = nrow(b), finite_area = sum(b$status == "ok"),
      missing_raster = sum(b$status == "missing_raster"),
      existing_manual_exclusions = sum(b$status == "skipped_manual_excluded"),
      unresolved_geometry = sum(!b$status %in% c("ok", "missing_raster", "skipped_manual_excluded")),
      recovered_areas = sum(grepl("integral_recovery", b$source, fixed = TRUE)))
  }
  initial_path <- file.path(backup, paste0("alpha_hull_0_5_", scale, "_all_pixels.csv"))
  if (file.exists(initial_path)) {
    initial <- read.csv(initial_path, stringsAsFactors = FALSE)
    initial <- initial[initial$status == "ok", ]
    key <- paste(d$quad_id, d$basis)
    j <- match(paste(initial$quad_id, initial$basis), key)
    finite_unchanged <- finite_unchanged && all(initial$alpha_area_all_pixels == d$alpha_area_all_pixels[j])
  }
  timing <- read.csv(file.path(OUT, paste0(scale, "_timing.csv")))
  first_start <- file.info(file.path(OUT, paste0(scale, "_input_manifest.csv")))$mtime
  timings[[scale]] <- data.frame(scale = scale,
    elapsed_pass_minutes = as.numeric(difftime(as.POSIXct(timing$finished), first_start, units = "mins")))
}
write.csv(do.call(rbind, coverage), file.path(OUT, "final_metric_coverage.csv"), row.names = FALSE)
write.csv(do.call(rbind, timings), file.path(OUT, "overall_scale_timings.csv"), row.names = FALSE)
stopifnot(finite_unchanged)
stopifnot(all(combined$scale %in% SCALES), all(combined$alpha_value == 0.5),
          !anyDuplicated(combined[c("scale", "quad_id", "basis", "alpha_value")]),
          all(is.finite(combined$alpha_area_all_pixels[combined$status == "ok"])),
          all(combined$alpha_area_all_pixels[combined$status == "ok"] >= 0),
          identical(unname(tools::md5sum(reference_hashes$path)), reference_hashes$md5))
write.csv(data.frame(check = c("arc_integral_validation_passed", "full_boundary_ids_accounted_for",
  "unique_result_keys", "all_successful_areas_finite_nonnegative", "existing_finite_areas_unchanged",
  "protected_reference_inputs_and_50m_unchanged", "only_authorized_scales"), pass = TRUE),
  file.path(OUT, "post_recovery_validation.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("Unresolved after area recovery:", nrow(unresolved), "\n50 m was not processed.\n")
if (requireNamespace("beepr", quietly = TRUE)) beepr::beep(3)
