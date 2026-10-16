# Focused all-pixel PC1-PC3 workflow. Never runs hulls or fits a new PCA basis.
user_library <- "C:/Users/PaintRock/AppData/Local/R/win-library/4.2"
.libPaths(c(user_library, .libPaths()))
stopifnot(requireNamespace("terra", quietly = TRUE), requireNamespace("sf", quietly = TRUE))
REFERENCE_SCRIPT <- "scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R"
helper <- new.env()
expressions <- parse(REFERENCE_SCRIPT)
stopifnot(grepl("RUN_SPECTRAL_HET_WORKFLOW", paste(deparse(expressions[[length(expressions)]]), collapse = "")))
eval(expressions[-length(expressions)], envir = helper)

OUT <- "reports/tables/20261015_pca_mean_rao_completion"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
original <- readRDS(helper$PCA_RDS)
normalized <- readRDS(helper$STANDARDIZED_PCA_RDS)
stopifnot(original$preprocess == "raw", normalized$preprocess == "vector_normalized",
          original$n_bands == normalized$n_bands, original$n_components >= 3,
          normalized$n_components >= 3)
protected <- c(helper$PCA_RDS, helper$STANDARDIZED_PCA_RDS,
               unlist(lapply(c("10m", "20m", "50m"), function(s) c(
                 paste0("Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_", s, "_smooth_masked_5nm_summary.csv"),
                 paste0("Quad_Values/", s, "_PCA_metrics_smooth_masked_5nm_summary.csv"),
                 paste0("Quad_Values/", s, "_standardized_PCA_metrics_smooth_masked_5nm_summary.csv")))))
before <- tools::md5sum(protected)
hash_path <- file.path(OUT, "reference_hashes_before.csv")
if (file.exists(hash_path)) {
  previous_hashes <- read.csv(hash_path)
  stopifnot(identical(previous_hashes$path, names(before)), identical(previous_hashes$md5, unname(before)))
} else write.csv(data.frame(path = names(before), md5 = unname(before)), hash_path, row.names = FALSE)

files10 <- helper$list_raster_files("Quad_Spectra/10m_smooth_5nm")
files25 <- helper$list_raster_files("Quad_Spectra/25m_smooth_5nm")
stopifnot(length(files10) > 0, length(files25) > 0,
          !anyDuplicated(basename(files10)), !anyDuplicated(basename(files25)))
old10 <- read.csv("Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_10m_smooth_masked_5nm_summary.csv")
stopifnot(setequal(basename(files10), old10$quad_id))
old_scales <- setNames(lapply(c("10m", "20m", "50m"), function(s) {
  d <- read.csv(paste0("Quad_Values/Spectral_diversitySHPs/spectral_heterogeneity_", s, "_smooth_masked_5nm_summary.csv"))
  for (pre in c("", "standardized_PCA_")) {
    component <- read.csv(paste0("Quad_Values/", s, "_", ifelse(pre == "", "", "standardized_"), "PCA_metrics_smooth_masked_5nm_summary.csv"))
    stopifnot(identical(as.character(component$quad_id), as.character(d$quad_id)))
    for (field in c("pca_euclidean_mean", "rao_q_pca")) stopifnot(isTRUE(all.equal(component[[paste0(pre, field)]], d[[paste0(pre, field)]])))
    eligible <- !d[[paste0(pre, "manual_excluded")]]
    stopifnot(all(is.finite(d[[paste0(pre, "pca_euclidean_mean")]][eligible])),
              all(is.finite(d[[paste0(pre, "rao_q_pca")]][eligible])))
  }
  d
}), c("10m", "20m", "50m"))

jobs <- rbind(data.frame(scale = "10m", quad_id = basename(files10), file = files10),
              data.frame(scale = "25m", quad_id = basename(files25), file = files25))
jobs$manual_excluded <- jobs$scale == "10m" & jobs$quad_id %in% old10$quad_id[old10$manual_excluded]
manifest <- data.frame(path = jobs$file, bytes = file.info(jobs$file)$size,
                       modified = as.numeric(file.info(jobs$file)$mtime))
manifest_path <- file.path(OUT, "input_manifest.csv")
if (file.exists(manifest_path)) {
  stopifnot(isTRUE(all.equal(read.csv(manifest_path), manifest, check.attributes = FALSE)))
} else write.csv(manifest, manifest_path, row.names = FALSE)
worker <- function(i) {
  job <- jobs[i, ]
  bases <- if (job$scale == "10m") "vector_normalized_pca" else c("original_pca", "vector_normalized_pca")
  result <- data.frame(scale = job$scale, quad_id = job$quad_id, basis = bases,
                       manual_excluded = if (job$scale == "10m") job$manual_excluded else NA,
                       n_pixels = NA_integer_, metric_method = "pending", pca_euclidean_mean = NA_real_,
                       rao_q_pca = NA_real_, rao_q_distance = "squared_euclidean_pc1_pc3",
                       status = "pending", error = "", raster_file = job$file)
  if (job$manual_excluded) {
    result$metric_method <- "manual_excluded"; result$status <- "skipped_manual_excluded"
    return(result)
  }
  tryCatch({
    r <- terra::rast(job$file)
    stopifnot(identical(names(r), names(original$center)), identical(names(r), names(normalized$center)))
    x <- helper$read_masked_spectra(job$file)
    for (j in seq_along(bases)) {
      basis <- if (bases[j] == "original_pca") original else normalized
      scores <- helper$project_scores(x, basis)
      result$n_pixels[j] <- nrow(scores)
      if (nrow(scores) < helper$MIN_VALID_PIXELS) {
        result$metric_method[j] <- "insufficient_pixels"; result$status[j] <- "insufficient_pixels"
      } else {
        sq <- rowSums(sweep(scores, 2, colMeans(scores), "-")^2)
        result$pca_euclidean_mean[j] <- mean(sqrt(sq))
        result$rao_q_pca[j] <- 2 * mean(sq)
        result$metric_method[j] <- "all_pixels"; result$status[j] <- "ok"
      }
    }
    result
  }, error = function(e) {
    result$status <- "processing_error"; result$metric_method <- "processing_error"
    result$error <- conditionMessage(e); result
  })
}

# Persist each completed chunk so a failed run does not lose previous calculations.
workers <- as.integer(Sys.getenv("PCA_DISTANCE_WORKERS", "4"))
cluster <- parallel::makePSOCKcluster(workers)
parallel::clusterExport(cluster, c("jobs", "helper", "original", "normalized", "worker"), envir = environment())
parallel::clusterEvalQ(cluster, {
  .libPaths(c("C:/Users/PaintRock/AppData/Local/R/win-library/4.2", .libPaths()))
  terra::terraOptions(progress = 0)
  NULL
})
parts <- list()
chunks <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / 80))
tryCatch({
  for (k in seq_along(chunks)) {
    chunk_path <- file.path(OUT, sprintf("computed_chunk_%03d.csv", k))
    if (file.exists(chunk_path)) {
      part <- read.csv(chunk_path, stringsAsFactors = FALSE)
      chunk_jobs <- jobs[chunks[[k]], ]
      expected <- do.call(rbind, lapply(seq_len(nrow(chunk_jobs)), function(i) {
        j <- chunk_jobs[i, ]; b <- if(j$scale == "10m") "vector_normalized_pca" else c("original_pca", "vector_normalized_pca")
        data.frame(scale = j$scale, quad_id = j$quad_id, basis = b, raster_file = j$file)
      }))
      stopifnot(isTRUE(all.equal(part[names(expected)], expected, check.attributes = FALSE)),
                all(part$status %in% c("ok", "skipped_manual_excluded")))
    } else {
      part <- do.call(rbind, parallel::parLapplyLB(cluster, chunks[[k]], worker))
      write.csv(part, chunk_path, row.names = FALSE)
    }
    parts[[k]] <- part
    cat(format(Sys.time()), "completed", max(chunks[[k]]), "of", nrow(jobs), "rasters\n"); flush.console()
  }
}, finally = parallel::stopCluster(cluster))
computed <- do.call(rbind, parts)
write.csv(computed, file.path(OUT, "computed_metrics_long.csv"), row.names = FALSE)
stopifnot(!anyDuplicated(computed[c("scale", "quad_id", "basis")]),
          all(computed$status %in% c("ok", "skipped_manual_excluded")))

v <- computed[computed$scale == "10m", ]
idx <- match(v$quad_id, old10$quad_id)
v$reference_n_pixels <- old10$standardized_PCA_n_pixels[idx]
v$reference_mean <- old10$standardized_PCA_pca_euclidean_mean[idx]
v$reference_rao <- old10$standardized_PCA_rao_q_pca[idx]
v$mean_difference <- v$pca_euclidean_mean - v$reference_mean
v$rao_difference <- v$rao_q_pca - v$reference_rao
ABS_TOL <- 1e-10
REL_TOL <- 1e-10
v$mean_matches <- abs(v$mean_difference) <= ABS_TOL + REL_TOL * abs(v$reference_mean)
v$rao_matches <- abs(v$rao_difference) <= ABS_TOL + REL_TOL * abs(v$reference_rao)
v$pixel_count_matches <- v$n_pixels == v$reference_n_pixels
write.csv(v, file.path(OUT, "restored_10m_vector_normalized_comparison.csv"), row.names = FALSE)
eligible <- v$status == "ok"
comparison_pass <- all(v$mean_matches[eligible]) && all(v$rao_matches[eligible]) && all(v$pixel_count_matches[eligible])
cat("10m comparison pass:", comparison_pass, "maximum mean difference", max(abs(v$mean_difference[eligible])), "\n")

# Retain every 25 m boundary ID. The user handles exclusions manually.
grid25 <- sf::st_read("Quad_Scale_SHPs/PR_25m.shp", quiet = TRUE)
stopifnot(!anyDuplicated(grid25$Name), all(basename(files25) %in% grid25$Name))
qc <- data.frame(quad_id = as.character(grid25$Name), scale = "25m", raster_present = grid25$Name %in% basename(files25))
write.csv(qc, file.path(OUT, "25m_raster_coverage.csv"), row.names = FALSE)
d25 <- computed[computed$scale == "25m", ]
components <- list()
for (basis in c("original_pca", "vector_normalized_pca")) {
  d <- d25[d25$basis == basis, ]
  d <- d[match(qc$quad_id, d$quad_id), c("quad_id", "n_pixels", "metric_method", "pca_euclidean_mean", "rao_q_pca", "rao_q_distance", "status")]
  d$quad_id <- qc$quad_id
  missing <- !qc$raster_present
  d$metric_method[missing] <- "missing_raster"; d$status[missing] <- "missing_raster"
  if (basis == "vector_normalized_pca") d <- helper$prefix_standardized_pca_metrics(d)
  components[[basis]] <- d
}
summary25 <- cbind(qc, components$original_pca[setdiff(names(components$original_pca), "quad_id")],
                   components$vector_normalized_pca[setdiff(names(components$vector_normalized_pca), "quad_id")])
canonical25 <- "Quad_Values/Spectral_diversitySHPs/spectral_mean_rao_25m_smooth_masked_5nm_summary.csv"
targets <- c("Quad_Values/25m_PCA_metrics_smooth_masked_5nm_summary.csv",
             "Quad_Values/25m_standardized_PCA_metrics_smooth_masked_5nm_summary.csv", canonical25)
if (any(file.exists(targets))) stop("New 25 m output already exists; review before overwriting.")
write.csv(components$original_pca, targets[1], row.names = FALSE)
write.csv(components$vector_normalized_pca, targets[2], row.names = FALSE)
write.csv(summary25, canonical25, row.names = FALSE)

coverage <- list()
for (s in c("10m", "20m", "25m", "50m")) {
  grid <- sf::st_read(paste0("Quad_Scale_SHPs/PR_", s, ".shp"), quiet = TRUE)
  d <- if (s == "25m") summary25 else old_scales[[s]]
  for (pre in c("", "standardized_PCA_")) {
    excluded <- d[[paste0(pre, "manual_excluded")]]
    coverage[[length(coverage) + 1]] <- data.frame(scale = s, basis = ifelse(pre == "", "original_pca", "vector_normalized_pca"),
      grid_quadrats = nrow(grid), summary_rows = nrow(d),
      quadrats_without_rasters = if (s == "25m") sum(!qc$raster_present) else nrow(grid) - nrow(d),
      manual_excluded = if (s == "25m") NA_integer_ else sum(excluded),
      finite_mean = sum(is.finite(d[[paste0(pre, "pca_euclidean_mean")]])),
      finite_rao = sum(is.finite(d[[paste0(pre, "rao_q_pca")]])),
      exclusion_policy = ifelse(s == "25m", "user_managed_no_exclusions_applied", "existing_policy"))
  }
}
coverage <- do.call(rbind, coverage)
write.csv(coverage, file.path(OUT, "metric_coverage.csv"), row.names = FALSE)
after <- tools::md5sum(protected)
stopifnot(identical(before, after),
          identical(manifest$bytes, file.info(jobs$file)$size),
          identical(manifest$modified, as.numeric(file.info(jobs$file)$mtime)),
          nrow(summary25) == nrow(grid25),
          all(summary25$pca_euclidean_mean[qc$raster_present] >= 0),
          all(summary25$rao_q_pca[qc$raster_present] >= 0))
checks <- data.frame(check = c("restored_10m_mean_matches", "restored_10m_rao_matches", "restored_10m_pixel_counts_match",
                             "reference_files_unchanged", "25m_all_available_rasters_computed", "25m_full_grid_accounted_for"),
                     pass = c(all(v$mean_matches[eligible]), all(v$rao_matches[eligible]), all(v$pixel_count_matches[eligible]),
                              identical(before, after), sum(d25$status == "ok") == 2 * length(files25), nrow(summary25) == nrow(grid25)))
write.csv(checks, file.path(OUT, "validation_checks.csv"), row.names = FALSE)
print(coverage); print(checks)
cat("10m max abs mean difference:", format(max(abs(v$mean_difference[eligible])), digits = 17), "\n")
cat("10m max abs Rao difference:", format(max(abs(v$rao_difference[eligible])), digits = 17), "\n")
cat("25m missing raster IDs:", paste(qc$quad_id[!qc$raster_present], collapse = ", "), "\n")
if (requireNamespace("beepr", quietly = TRUE)) beepr::beep(3)
if (!comparison_pass) stop("Restored 10 m values differ; review the comparison CSV before accepting restored data.")
