# Accept only the validated, unchanged-point numerical circumcenter correction.
.libPaths(c('C:/Users/PaintRock/AppData/Local/R/win-library/4.2', .libPaths()))
source('scripts/3_Analysis/alpha_hull_arc_integral.R')
same_rows <- function(x, y) {
  identical(names(x), names(y)) && nrow(x) == nrow(y) &&
    all(vapply(names(x), function(n) identical(x[[n]], y[[n]]) ||
      (all(is.na(x[[n]])) && all(is.na(y[[n]]))), logical(1)))
}
root <- 'reports/tables/20261015_alpha_half_small_scales'
out <- file.path(root, '1000_c_numeric_geometry')
canonical <- 'Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_10m_all_pixels.csv'
checkpoint <- file.path(root, '10m/1000_c.csv')
saved <- readRDS(file.path(out, 'triangulation.rds'))
trials <- read.csv(file.path(out, 'circumcenter_recovery_trials.csv'))
stopifnot(nrow(trials) == 2, all(trials$replacements == 1), all(trials$mesh_nonfinite == 0),
  all(is.finite(as.matrix(trials[c('package_area', 'arc_area', 'endpoint_area')]))),
  diff(range(as.matrix(trials[c('package_area', 'arc_area', 'endpoint_area')]))) < 1e-8)
h <- readRDS(file.path(out, 'corrected_hull_translated_determinant.rds'))
h2 <- readRDS(file.path(out, 'corrected_hull_translated_linear_solve.rds'))
stopifnot(identical(h$arcs, h2$arcs), identical(h$xahull, h2$xahull),
  isTRUE(all.equal(h, h2)), h$alpha == .5)
dd <- readRDS(file.path(out, 'corrected_delvor_translated_determinant.rds'))
# Compare the corrected mesh against standard evaluation on the exact saved
# triangulation, bypassing any new triangulation or internal jitter retry.
private <- new.env(parent = asNamespace('alphahull'))
private$tri.mesh <- function(x, y = NULL, ...) saved$triangulation
standard_delvor <- get('delvor', asNamespace('alphahull'))
environment(standard_delvor) <- private
baseline <- standard_delvor(saved$points[, 1], saved$points[, 2])
bad <- !is.finite(baseline$mesh)
stopifnot(sum(bad) == 8, all(is.finite(dd$mesh)),
  identical(baseline$mesh[!bad], dd$mesh[!bad]),
  identical(unname(dd$x), unname(saved$points)), identical(dd$tri.obj, saved$triangulation))
area <- alpha_hull_endpoint_integral(h)
stopifnot(abs(area - alpha_hull_arc_integral(h)) < 1e-8,
  abs(area - trials$package_area[1]) < 1e-8)
p <- saved$points
ids <- chull(p)
polygon <- p[c(ids, ids[1]), , drop = FALSE]
bound <- abs(sum(polygon[-nrow(polygon), 1] * polygon[-1, 2] -
  polygon[-1, 1] * polygon[-nrow(polygon), 2])) / 2
stopifnot(is.finite(area), area > 0, area <= bound)
protected <- read.csv(file.path(root, 'protected_hashes.csv'), stringsAsFactors = FALSE)
stopifnot(identical(unname(tools::md5sum(protected$path)), protected$md5))
unchanged_paths <- file.path('Quad_Values/Spectral_diversitySHPs',
  c('alpha_hull_0_5_20m_all_pixels.csv', 'alpha_hull_0_5_25m_all_pixels.csv'))
hashes <- unname(tools::md5sum(unchanged_paths))
write.csv(data.frame(path = unchanged_paths, md5 = hashes), file.path(out, 'other_scale_hashes.csv'), row.names = FALSE)
paths <- c(canonical, checkpoint, file.path(root,
  c('final_metric_coverage.csv', 'unresolved_geometry.csv', '10m_status_summary.csv')))
for (path in paths) {
  backup <- file.path(out, paste0('before_', basename(path)))
  if (file.exists(backup)) {
    stopifnot(unname(tools::md5sum(path)) == unname(tools::md5sum(backup)))
  } else stopifnot(file.copy(path, backup))
}
d <- read.csv(canonical, stringsAsFactors = FALSE)
old <- d
i <- which(d$quad_id == '1000_c' & d$basis == 'original_pca')
stopifnot(length(i) == 1, is.na(d$alpha_area_all_pixels[i]),
  d$n_retained_pixels[i] == 7126, d$n_unique_pc12_points[i] == nrow(p))
d$alpha_area_all_pixels[i] <- area
d$status[i] <- 'ok'
d$arc_count[i] <- nrow(h$arcs)
d$source[i] <- 'computed_all_pixels_alpha_0_5_stable_circumcenter_recovery'
d$warning[i] <- paste('Recovered one nonfinite circumcenter using translated-coordinate',
  'formula on unchanged points and triangulation; alpha 0.5 unchanged; no manual jitter',
  'or orientation repair. See 1000_c_numeric_geometry audit.')
d$error[i] <- NA
d$run_timestamp[i] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')
stopifnot(identical(d[-i, ], old[-i, ]))
cp <- read.csv(checkpoint, stringsAsFactors = FALSE,
  colClasses = c(warning = 'character', error = 'character'))
old_cp <- cp; j <- which(cp$basis == 'original_pca')
cp[j, ] <- d[i, ]
stopifnot(same_rows(cp[-j, ], old_cp[-j, ]))
write.csv(d, canonical, row.names = FALSE)
write.csv(cp, checkpoint, row.names = FALSE)
coverage_path <- file.path(root, 'final_metric_coverage.csv')
coverage <- read.csv(coverage_path, stringsAsFactors = FALSE)
j <- which(coverage$scale == '10m' & coverage$basis == 'original_pca')
coverage$finite_area[j] <- sum(d$basis == 'original_pca' & is.finite(d$alpha_area_all_pixels))
coverage$unresolved_geometry[j] <- sum(d$basis == 'original_pca' &
  !d$status %in% c('ok', 'missing_raster', 'skipped_manual_excluded'))
write.csv(coverage, coverage_path, row.names = FALSE)
# Historical recovered_areas remains the count from the initial integration pass;
# later geometry corrections have dedicated provenance and audit records.
unresolved_path <- file.path(root, 'unresolved_geometry.csv')
unresolved <- read.csv(unresolved_path, stringsAsFactors = FALSE)
stopifnot(nrow(unresolved) == 1, unresolved$quad_id == '1000_c')
write.csv(unresolved[FALSE, ], unresolved_path, row.names = FALSE)
status <- as.data.frame(table(d$scale, d$basis, d$status))
names(status) <- c('scale', 'basis', 'status', 'rows')
write.csv(status[status$rows > 0, ], file.path(root, '10m_status_summary.csv'), row.names = FALSE)
diagnosis <- read.csv(file.path(out, 'nonfinite_standard_circumcenters.csv'))
audit <- data.frame(scale = '10m', quad_id = '1000_c', basis = 'original_pca', alpha_value = .5,
  method = 'stable_circumcenter_same_triangulation', affected_triangle = diagnosis$triangle,
  determinant = diagnosis$determinant, triangle_height = diagnosis$height_longest_edge,
  n_retained_pixels = 7126L, n_unique_points = nrow(p), replaced_mesh_entries = sum(bad),
  recovered_area = area, convex_bound = bound, timestamp = d$run_timestamp[i])
write.csv(audit, file.path(out, 'accepted_recovery_audit.csv'), row.names = FALSE)
readback <- read.csv(canonical, stringsAsFactors = FALSE, colClasses = c(error = 'character'))
stopifnot(same_rows(readback[-i, ], old[-i, ]),
  abs(readback$alpha_area_all_pixels[i] - area) < 1e-10,
  coverage$finite_area[j] == 1744, coverage$unresolved_geometry[j] == 0,
  sum(coverage$finite_area) == 4984,
  identical(unname(tools::md5sum(unchanged_paths)), hashes),
  identical(unname(tools::md5sum(protected$path)), protected$md5))
write.csv(data.frame(check = c('global_cloud_not_collinear', 'no_zero_area_triangles',
  'one_nonfinite_center_replaced', 'original_points_and_triangulation_unchanged',
  'only_eight_originally_nonfinite_mesh_entries_changed', 'two_center_formulas_same_hull',
  'package_and_two_independent_area_integrals_agree', 'closed_boundary_and_circle_consistency',
  'finite_nonnegative_within_convex_bound', 'other_canonical_rows_and_VN_checkpoint_unchanged',
  'other_scales_and_protected_references_unchanged', 'saved_result_round_trip',
  'all_small_scale_available_hulls_complete'), pass = TRUE),
  file.path(out, 'accepted_validation.csv'), row.names = FALSE)
print(audit, digits = 16)
cat('Accepted stable circumcenter recovery: PASS\n')
