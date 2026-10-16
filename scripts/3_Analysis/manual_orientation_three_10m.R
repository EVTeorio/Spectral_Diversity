# User-reviewed exceptions only; do not generalize to unreviewed hulls.
root <- 'reports/tables/20261015_alpha_half_small_scales'
diagnostic <- file.path(root, 'diagnose_four_10m')
out <- file.path(root, 'manual_orientation_three_10m')
canonical <- 'Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_10m_all_pixels.csv'
reviewed <- c('1402_b' = 353L, '1019_a' = 487L, '1816_c' = 545L)
dir.create(out, showWarnings = FALSE)
protected <- read.csv(file.path(root, 'protected_hashes.csv'), stringsAsFactors = FALSE)
stopifnot(identical(unname(tools::md5sum(protected$path)), protected$md5))
unchanged_files <- c(file.path('Quad_Values/Spectral_diversitySHPs',
  c('alpha_hull_0_5_20m_all_pixels.csv', 'alpha_hull_0_5_25m_all_pixels.csv')),
  file.path(root, '10m/1000_c.csv'),
  file.path(diagnostic, paste0(names(reviewed), '_hull.rds')))
unchanged_hash <- unname(tools::md5sum(unchanged_files))
write.csv(data.frame(path = unchanged_files, md5 = unchanged_hash),
  file.path(out, 'unchanged_file_hashes.csv'), row.names = FALSE)
d <- read.csv(canonical, stringsAsFactors = FALSE)
before <- d
affected <- which(d$basis == 'original_pca' & d$quad_id %in% names(reviewed))
stopifnot(length(affected) == 3, all(is.na(d$alpha_area_all_pixels[affected])))
paths <- c(canonical, file.path(root, c('final_metric_coverage.csv',
  'unresolved_geometry.csv', '10m_status_summary.csv')),
  file.path(root, '10m', paste0(names(reviewed), '.csv')))
for (path in paths) {
  backup <- file.path(out, paste0('before_', basename(path)))
  if (file.exists(backup)) {
    stopifnot(unname(tools::md5sum(path)) == unname(tools::md5sum(backup)))
  } else stopifnot(file.copy(path, backup))
}
audits <- list()
for (id in names(reviewed)) {
  h <- readRDS(file.path(diagnostic, paste0(id, '_hull.rds')))
  active <- which(h$arcs[, 3] > 0)
  a <- h$arcs[active, , drop = FALSE]
  k <- match(reviewed[[id]], active)
  stopifnot(!is.na(k), all(is.finite(a)))
  r <- a[, 3]; theta <- a[, 6]
  centers <- a[, 1:2, drop = FALSE]
  p <- h$xahull[as.integer(a[, 7]), , drop = FALSE]
  q <- h$xahull[as.integer(a[, 8]), , drop = FALSE]
  tol <- 1e-7 * pmax(1, r)
  stopifnot(all(theta >= 0), all(theta <= pi + 1e-10),
    all(abs(rowSums(a[, 4:5, drop = FALSE]^2) - 1) < 1e-8),
    all(abs(sqrt(rowSums((p - centers)^2)) - r) <= tol),
    all(abs(sqrt(rowSums((q - centers)^2)) - r) <= tol))
  flip <- a[, 4] * (p[, 2] - centers[, 2]) - a[, 5] * (p[, 1] - centers[, 1]) < 0
  s <- as.integer(a[, 7]); t <- as.integer(a[, 8])
  s[flip] <- a[flip, 8]; t[flip] <- a[flip, 7]
  old_balance <- tabulate(s, nrow(h$xahull)) - tabulate(t, nrow(h$xahull))
  bad <- which(old_balance != 0)
  stopifnot(length(bad) == 2, setequal(c(s[k], t[k]), bad))
  direction <- rep(1, nrow(a)); direction[k] <- -1
  tmp <- s[k]; s[k] <- t[k]; t[k] <- tmp
  stopifnot(all(tabulate(s, nrow(h$xahull)) == tabulate(t, nrow(h$xahull))))
  origin <- colMeans(h$xahull)
  sp <- sweep(h$xahull[s, , drop = FALSE], 2, origin)
  tp <- sweep(h$xahull[t, , drop = FALSE], 2, origin)
  terms <- (sp[, 1] * tp[, 2] - sp[, 2] * tp[, 1] -
    direction * r^2 * (2 * theta - sin(2 * theta))) / 2
  area <- sum(terms)
  vx <- a[, 4]; vy <- a[, 5]
  start <- centers + cbind(r * (vx * cos(theta) - vy * sin(theta)),
    r * (vx * sin(theta) + vy * cos(theta)))
  end <- centers + cbind(r * (vx * cos(theta) + vy * sin(theta)),
    r * (-vx * sin(theta) + vy * cos(theta)))
  delta <- (end - start) * direction
  closure <- sqrt(sum(colSums(delta)^2))
  cc <- sweep(centers, 2, origin)
  center_area <- sum((cc[, 1] * delta[, 2] - cc[, 2] * delta[, 1] -
    direction * 2 * r^2 * theta) / 2)
  hull_ids <- chull(h$xahull)
  polygon <- sweep(h$xahull[c(hull_ids, hull_ids[1]), , drop = FALSE], 2, origin)
  bound <- abs(sum(polygon[-nrow(polygon), 1] * polygon[-1, 2] -
    polygon[-1, 1] * polygon[-nrow(polygon), 2])) / 2
  stopifnot(is.finite(area), area >= 0, area <= bound + 1e-8,
    closure < 1e-8 * max(1, sum(2 * r * theta)), abs(area - center_area) < 1e-8)
  record <- data.frame(scale = '10m', quad_id = id, basis = 'original_pca',
    alpha_value = .5, reversed_arc = reviewed[[id]], repaired_area = area,
    independent_center_area = center_area, absolute_integral_difference = abs(area - center_area),
    closure_error = closure, convex_bound = bound,
    timestamp = format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'))
  audits[[id]] <- record
  write.csv(data.frame(arc = active, direction = direction, start_id = s,
    end_id = t, area_term = terms), file.path(out, paste0(id, '_repaired_traversal.csv')),
    row.names = FALSE)
  i <- which(d$quad_id == id & d$basis == 'original_pca')
  stopifnot(length(i) == 1, nrow(h$arcs) == d$arc_count[i])
  d$alpha_area_all_pixels[i] <- area
  d$status[i] <- 'ok'
  d$source[i] <- 'computed_all_pixels_alpha_0_5_manual_arc_orientation_repair'
  d$warning[i] <- paste('User-reviewed manual repair: reversed traversal of saved arc',
    reviewed[[id]], '; points, curves and alpha unchanged. Original package area was NA;',
    'see manual_orientation_three_10m audit.')
  d$error[i] <- NA
  d$run_timestamp[i] <- record$timestamp
}
# All three numerical validations finish before any canonical/checkpoint write.
audit <- do.call(rbind, audits)
write.csv(audit, file.path(out, 'manual_repair_audit.csv'), row.names = FALSE)
stopifnot(identical(d[-affected, ], before[-affected, ]))
for (id in names(reviewed)) {
  path <- file.path(root, '10m', paste0(id, '.csv'))
  cp <- read.csv(path, stringsAsFactors = FALSE)
  old_cp <- cp
  j <- which(cp$basis == 'original_pca')
  replacement <- d[d$quad_id == id & d$basis == 'original_pca', ]
  # read.csv infers all-NA checkpoint columns as logical; preserve that schema.
  for (name in names(cp)) {
    if (is.logical(cp[[name]]) && all(is.na(cp[[name]]))) {
      stopifnot(is.na(replacement[[name]]))
      replacement[[name]] <- as.logical(replacement[[name]])
    }
  }
  cp[j, ] <- replacement
  stopifnot(identical(cp[-j, ], old_cp[-j, ]))
  write.csv(cp, path, row.names = FALSE)
}
write.csv(d, canonical, row.names = FALSE)
coverage_path <- file.path(root, 'final_metric_coverage.csv')
coverage <- read.csv(coverage_path, stringsAsFactors = FALSE)
j <- which(coverage$scale == '10m' & coverage$basis == 'original_pca')
coverage$finite_area[j] <- sum(d$basis == 'original_pca' & is.finite(d$alpha_area_all_pixels))
coverage$unresolved_geometry[j] <- sum(d$basis == 'original_pca' &
  !d$status %in% c('ok', 'missing_raster', 'skipped_manual_excluded'))
write.csv(coverage, coverage_path, row.names = FALSE)
# recovered_areas retains its earlier integral-only recovery count; manual
# corrections are distinguished in canonical source/warning and dedicated audits.
unresolved_path <- file.path(root, 'unresolved_geometry.csv')
unresolved <- read.csv(unresolved_path, stringsAsFactors = FALSE)
unresolved <- unresolved[!(unresolved$scale == '10m' & unresolved$basis == 'original_pca' &
  unresolved$quad_id %in% names(reviewed)), ]
stopifnot(nrow(unresolved) == 1, unresolved$quad_id == '1000_c')
write.csv(unresolved, unresolved_path, row.names = FALSE)
status <- as.data.frame(table(d$scale, d$basis, d$status))
names(status) <- c('scale', 'basis', 'status', 'rows')
write.csv(status[status$rows > 0, ], file.path(root, '10m_status_summary.csv'), row.names = FALSE)
saved <- read.csv(canonical, stringsAsFactors = FALSE)
stopifnot(identical(saved[-affected, ], before[-affected, ]),
  max(abs(saved$alpha_area_all_pixels[affected] - d$alpha_area_all_pixels[affected])) < 1e-10,
  identical(unname(tools::md5sum(unchanged_files)), unchanged_hash),
  identical(unname(tools::md5sum(protected$path)), protected$md5),
  coverage$finite_area[j] == 1743, sum(coverage$finite_area) == 4983)
for (id in names(reviewed)) {
  cp <- read.csv(file.path(root, '10m', paste0(id, '.csv')), stringsAsFactors = FALSE)
  backup <- read.csv(file.path(out, paste0('before_', id, '.csv')), stringsAsFactors = FALSE)
  stopifnot(identical(cp[cp$basis != 'original_pca', ], backup[backup$basis != 'original_pca', ]),
    identical(lapply(cp[cp$basis == 'original_pca', ], as.character),
      lapply(saved[saved$basis == 'original_pca' & saved$quad_id == id, ], as.character)))
}
write.csv(data.frame(check = c('all_three_exact_vertex_balance', 'endpoint_circle_consistency',
  'independent_integral_agreement', 'reconstructed_boundary_closure', 'convex_area_bounds',
  'only_three_authorized_canonical_rows_changed', 'vector_normalized_checkpoints_unchanged',
  '1000_c_20m_25m_and_saved_hulls_unchanged', 'protected_references_unchanged',
  'saved_results_round_trip', 'coverage_and_unresolved_consistent'), pass = TRUE),
  file.path(out, 'validation.csv'), row.names = FALSE)
print(audit, digits = 16)
cat('Three reviewed manual repairs: PASS\n')
