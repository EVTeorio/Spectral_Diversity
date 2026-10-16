# All-pixel alpha = 0.5 hulls, restricted to 10 m -> 20 m -> 25 m.
.libPaths(c("C:/Users/PaintRock/AppData/Local/R/win-library/4.2", .libPaths()))
stopifnot(all(vapply(c("terra", "sf", "alphahull"), requireNamespace, logical(1), quietly = TRUE)))
SCALES <- c("10m", "20m", "25m")
stopifnot(identical(SCALES, c("10m", "20m", "25m")))
BASES <- c("original_pca", "vector_normalized_pca")
ALPHA <- 0.5
AREA_TIMEOUT <- 300
WORKERS <- as.integer(Sys.getenv("SMALL_HULL_WORKERS", "8"))
OUT <- "reports/tables/20261015_alpha_half_small_scales"
CANONICAL <- "Quad_Values/Spectral_diversitySHPs"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
helper <- new.env()
expressions <- parse("scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R")
stopifnot(grepl("RUN_SPECTRAL_HET_WORKFLOW", paste(deparse(expressions[[length(expressions)]]), collapse = "")))
eval(expressions[-length(expressions)], helper)
models <- list(original_pca = readRDS(helper$PCA_RDS), vector_normalized_pca = readRDS(helper$STANDARDIZED_PCA_RDS))
historical_path <- "reports/tables/alpha_hull_all_pixels_values/alpha_hull_all_pixels_values.csv"
historical <- read.csv(historical_path, stringsAsFactors = FALSE)
historical <- historical[historical$scale %in% SCALES & historical$alpha_value == ALPHA, ]
stopifnot(!anyDuplicated(historical[c("scale", "quad_id", "basis", "alpha_value")]))
protected <- c(helper$PCA_RDS, helper$STANDARDIZED_PCA_RDS, historical_path,
  paste0(CANONICAL, "/spectral_heterogeneity_", c("10m", "20m", "50m"), "_smooth_masked_5nm_summary.csv"),
  paste0(CANONICAL, "/spectral_mean_rao_25m_smooth_masked_5nm_summary.csv"),
  "Quad_Values/25m_PCA_metrics_smooth_masked_5nm_summary.csv",
  "Quad_Values/25m_standardized_PCA_metrics_smooth_masked_5nm_summary.csv")
hashes <- tools::md5sum(protected)
hash_file <- file.path(OUT, "protected_hashes.csv")
if (file.exists(hash_file)) {
  old_hashes <- read.csv(hash_file)
  stopifnot(identical(old_hashes$path, names(hashes)), identical(old_hashes$md5, unname(hashes)))
} else write.csv(data.frame(path = names(hashes), md5 = unname(hashes)), hash_file, row.names = FALSE)

geometry_area <- function(points) {
  points <- unique(points)
  result <- list(n_unique_pc12_points = nrow(points), alpha_area_all_pixels = NA_real_,
                 status = "insufficient_unique_points", arc_count = NA_integer_,
                 elapsed_seconds = 0, warning = "", error = "")
  if (nrow(points) < 4) return(result)
  start <- proc.time()[["elapsed"]]
  messages <- character()
  tryCatch(withCallingHandlers({
    hull <- alphahull::ahull(points[, 1], points[, 2], alpha = ALPHA)
    area <- alphahull::areaahull(hull, timeout = AREA_TIMEOUT)
    result$alpha_area_all_pixels <- area
    result$arc_count <- nrow(hull$arcs)
    result$status <- if (is.finite(area) && area >= 0) "ok" else "nonfinite_area"
  }, warning = function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {
    result$status <<- "geometry_error"; result$error <<- conditionMessage(e)
  })
  result$elapsed_seconds <- proc.time()[["elapsed"]] - start
  result$warning <- paste(unique(messages), collapse = " | ")
  result
}

process_one <- function(i) {
  job <- jobs[i, ]
  checkpoint <- file.path(checkpoint_dir, paste0(job$quad_id, ".csv"))
  if (file.exists(checkpoint)) {
    d <- read.csv(checkpoint, stringsAsFactors = FALSE)
    stopifnot(nrow(d) == 2, all(d$scale == job$scale), all(d$quad_id == job$quad_id),
              identical(d$basis, BASES), all(d$alpha_value == ALPHA))
    return(d)
  }
  d <- data.frame(scale = job$scale, quad_id = job$quad_id, basis = BASES,
    alpha_value = ALPHA, n_retained_pixels = NA_integer_, n_unique_pc12_points = NA_integer_,
    alpha_area_all_pixels = NA_real_, status = "pending", arc_count = NA_integer_,
    elapsed_seconds = 0, raster_file = job$raster_file, source = "", run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    warning = "", error = "", stringsAsFactors = FALSE)
  if (is.na(job$raster_file)) {
    d$status <- "missing_raster"; d$source <- "input_inventory"
  } else if (job$reference_excluded) {
    d$status <- "skipped_manual_excluded"; d$source <- "existing_reference_policy"
  } else {
    for (b in seq_along(BASES)) {
      previous <- historical[historical$scale == job$scale & historical$quad_id == job$quad_id & historical$basis == BASES[b], ]
      expected <- if (b == 1) job$expected_raw_pixels else job$expected_normalized_pixels
      if (nrow(previous) == 1 && previous$status == "ok" &&
          is.finite(previous$alpha_area_all_pixels) && previous$alpha_area_all_pixels >= 0 &&
          is.finite(expected) && previous$n_retained_pixels == expected) {
        fields <- c("n_retained_pixels", "n_unique_pc12_points", "alpha_area_all_pixels", "status", "arc_count")
        d[b, fields] <- previous[fields]
        d$source[b] <- "reused_existing_alpha_0_5"
      }
    }
    needed <- which(d$status == "pending")
    if (length(needed)) {
      tryCatch({
        r <- terra::rast(job$raster_file)
        stopifnot(identical(names(r), names(models$original_pca$center)),
                  identical(names(r), names(models$vector_normalized_pca$center)))
        x <- helper$read_masked_spectra(job$raster_file)
        for (b in needed) {
          scores <- helper$project_scores(x, models[[BASES[b]]])
          expected <- if (b == 1) job$expected_raw_pixels else job$expected_normalized_pixels
          if (is.finite(expected)) stopifnot(nrow(scores) == expected)
          d$n_retained_pixels[b] <- nrow(scores)
          result <- geometry_area(scores[, 1:2, drop = FALSE])
          for (field in names(result)) d[b, field] <- result[[field]]
          d$source[b] <- "computed_all_pixels_alpha_0_5"
        }
      }, error = function(e) {
        pending <- which(d$status == "pending")
        d$status[pending] <- "processing_error"; d$error[pending] <- conditionMessage(e)
        d$source[pending] <- "computed_all_pixels_alpha_0_5"
      })
    }
  }
  temporary <- paste0(checkpoint, ".tmp")
  write.csv(d, temporary, row.names = FALSE)
  stopifnot(file.rename(temporary, checkpoint))
  d
}

cluster <- parallel::makePSOCKcluster(WORKERS)
parallel::clusterExport(cluster, c("helper", "models", "historical", "BASES", "ALPHA", "AREA_TIMEOUT", "geometry_area", "process_one"), envir = environment())
invisible(parallel::clusterEvalQ(cluster, {
  .libPaths(c("C:/Users/PaintRock/AppData/Local/R/win-library/4.2", .libPaths()))
  terra::terraOptions(progress = 0)
  NULL
}))
all_results <- list()
tryCatch({
  for (scale in SCALES) {
    started <- Sys.time()
    grid <- sf::st_read(paste0("Quad_Scale_SHPs/PR_", scale, ".shp"), quiet = TRUE)
    id_field <- if (scale == "10m") "sub_id" else "Name"
    ids <- as.character(grid[[id_field]])
    stopifnot(!anyDuplicated(ids))
    p <- if (scale == "25m") paste0(CANONICAL, "/spectral_mean_rao_25m_smooth_masked_5nm_summary.csv") else
      paste0(CANONICAL, "/spectral_heterogeneity_", scale, "_smooth_masked_5nm_summary.csv")
    reference <- read.csv(p, stringsAsFactors = FALSE)
    files <- helper$list_raster_files(paste0("Quad_Spectra/", scale, "_smooth_5nm"))
    file_ids <- vapply(files, helper$get_quad_id, character(1), scale = scale)
    stopifnot(!anyDuplicated(file_ids), all(file_ids %in% ids))
    idx <- match(ids, reference$quad_id)
    jobs <- data.frame(scale = scale, quad_id = ids, raster_file = files[match(ids, file_ids)],
      reference_excluded = if (scale == "25m") FALSE else !is.na(idx) & reference$manual_excluded[idx],
      expected_raw_pixels = reference$n_pixels[idx], expected_normalized_pixels = reference$standardized_PCA_n_pixels[idx])
    checkpoint_dir <- file.path(OUT, scale)
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    manifest <- data.frame(path = files, bytes = file.info(files)$size, modified = as.numeric(file.info(files)$mtime))
    manifest_path <- file.path(OUT, paste0(scale, "_input_manifest.csv"))
    if (file.exists(manifest_path)) {
      stopifnot(isTRUE(all.equal(read.csv(manifest_path), manifest, check.attributes = FALSE)))
    } else write.csv(manifest, manifest_path, row.names = FALSE)
    parallel::clusterExport(cluster, c("jobs", "checkpoint_dir"), envir = environment())
    cat(format(Sys.time()), "START", scale, "grid", nrow(jobs), "workers", WORKERS, "\n"); flush.console()
    # Stream single-quadrat jobs across the whole scale; a slow geometry must
    # not leave other workers idle at a small batch boundary. Checkpoints remain
    # visible to the monitoring process while this call runs.
    batches <- list(seq_len(nrow(jobs)))
    results <- list()
    for (k in seq_along(batches)) {
      results[[k]] <- do.call(rbind, parallel::parLapplyLB(cluster, batches[[k]], process_one, chunk.size = 1L))
      cat(format(Sys.time()), scale, "completed grid rasters", max(batches[[k]]), "/", nrow(jobs), "\n"); flush.console()
    }
    result <- do.call(rbind, results)
    result <- result[order(match(result$quad_id, ids), match(result$basis, BASES)), ]
    stopifnot(nrow(result) == 2 * nrow(grid), !anyDuplicated(result[c("scale", "quad_id", "basis", "alpha_value")]),
              all(result$scale == scale), all(result$alpha_value == 0.5),
              all(is.finite(result$alpha_area_all_pixels[result$status == "ok"])),
              all(result$alpha_area_all_pixels[result$status == "ok"] >= 0))
    stopifnot(identical(manifest$bytes, file.info(files)$size), identical(manifest$modified, as.numeric(file.info(files)$mtime)))
    canonical_path <- file.path(CANONICAL, paste0("alpha_hull_0_5_", scale, "_all_pixels.csv"))
    write.csv(result, canonical_path, row.names = FALSE)
    all_results[[scale]] <- result
    summary <- as.data.frame(table(result$scale, result$basis, result$status))
    names(summary) <- c("scale", "basis", "status", "rows")
    summary <- summary[summary$rows > 0, ]
    write.csv(summary, file.path(OUT, paste0(scale, "_status_summary.csv")), row.names = FALSE)
    timing <- data.frame(scale = scale, started = as.character(started), finished = as.character(Sys.time()),
                         wall_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    write.csv(timing, file.path(OUT, paste0(scale, "_timing.csv")), row.names = FALSE)
    cat(format(Sys.time()), "FINISHED", scale, "wall seconds", timing$wall_seconds, "\n"); print(summary); flush.console()
    stopifnot(identical(hashes, tools::md5sum(protected)))
  }
}, finally = parallel::stopCluster(cluster))
combined <- do.call(rbind, all_results)
failures <- combined[!combined$status %in% c("ok", "skipped_manual_excluded", "missing_raster"), ]
write.csv(failures, file.path(OUT, "unresolved_geometry.csv"), row.names = FALSE)
write.csv(data.frame(check = c("only_10m_20m_25m", "alpha_is_half", "unique_keys", "protected_inputs_and_50m_unchanged"),
  pass = c(all(combined$scale %in% SCALES), all(combined$alpha_value == 0.5),
           !anyDuplicated(combined[c("scale", "quad_id", "basis", "alpha_value")]), identical(hashes, tools::md5sum(protected)))),
  file.path(OUT, "validation_checks.csv"), row.names = FALSE)
cat("Unresolved geometry rows:", nrow(failures), "\nSTOPPED BEFORE 50 m.\n")
if (requireNamespace("beepr", quietly = TRUE)) beepr::beep(3)
