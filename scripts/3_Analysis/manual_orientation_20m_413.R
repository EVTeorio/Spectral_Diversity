# User-authorized exception for one saved hull; never an automatic fallback.
source('scripts/3_Analysis/alpha_hull_arc_integral.R')
root <- 'reports/tables/20261015_alpha_half_small_scales'
retry <- file.path(root, 'retry_20m_413_original_20261015_161612')
out <- file.path(root, 'manual_orientation_20m_413')
dir.create(out, showWarnings = FALSE)
canonical <- 'Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_20m_all_pixels.csv'
checkpoint <- file.path(root, '20m/413.csv')
coverage_path <- file.path(root, 'final_metric_coverage.csv')
unresolved_path <- file.path(root, 'unresolved_geometry.csv')
status_path <- file.path(root, '20m_status_summary.csv')
for (f in c(canonical, checkpoint, coverage_path, unresolved_path, status_path)) {
  backup <- file.path(out, paste0('before_', basename(f)))
  if (file.exists(backup)) {
    stopifnot(unname(tools::md5sum(f)) == unname(tools::md5sum(backup)))
  } else stopifnot(file.copy(f, backup))
}
protected <- read.csv(file.path(root, 'protected_hashes.csv'), stringsAsFactors = FALSE)
stopifnot(identical(unname(tools::md5sum(protected$path)), protected$md5))
h <- readRDS(file.path(retry, 'retry_hull.rds'))
stopifnot(nrow(h$arcs) == 561)
# Zero-radius entries are isolated points, not boundary curves.
active <- which(h$arcs[, 3] > 0)
a <- h$arcs[active, , drop = FALSE]
loop <- which(a[, 7] %in% c(4324, 4325, 4535) |
              a[, 8] %in% c(4324, 4325, 4535))
stopifnot(identical(loop, c(293L, 422L, 423L)),
          setequal(a[423, 7:8], c(4324, 4535)))
# Verify that removing the entire isolated component leaves a closed boundary.
rest <- h; rest$arcs <- a[-loop, , drop = FALSE]
rest_area <- alpha_hull_endpoint_integral(rest)
stopifnot(abs(rest_area - alpha_hull_arc_integral(rest)) < 1e-8)
p <- h$xahull[as.integer(a[, 7]), , drop = FALSE]
q <- h$xahull[as.integer(a[, 8]), , drop = FALSE]
centers <- a[, 1:2, drop = FALSE]
tol <- 1e-7 * pmax(1, a[, 3])
stopifnot(all(abs(sqrt(rowSums((p - centers)^2)) - a[, 3]) <= tol),
          all(abs(sqrt(rowSums((q - centers)^2)) - a[, 3]) <= tol))
flip <- a[, 4] * (p[, 2] - centers[, 2]) - a[, 5] * (p[, 1] - centers[, 1]) < 0
s <- as.integer(a[, 7]); t <- as.integer(a[, 8])
s[flip] <- a[flip, 8]; t[flip] <- a[flip, 7]
# Preserve the ahull object itself. Record the orientation override separately.
direction <- rep(1, nrow(a)); direction[423] <- -1
tmp <- s[423]; s[423] <- t[423]; t[423] <- tmp
stopifnot(all(tabulate(s, nrow(h$xahull)) == tabulate(t, nrow(h$xahull))))
origin <- colMeans(h$xahull)
sp <- sweep(h$xahull[s, , drop = FALSE], 2, origin)
tp <- sweep(h$xahull[t, , drop = FALSE], 2, origin)
endpoint_terms <- (sp[, 1] * tp[, 2] - sp[, 2] * tp[, 1] -
  direction * a[, 3]^2 * (2 * a[, 6] - sin(2 * a[, 6]))) / 2
area <- sum(endpoint_terms)
theta <- a[, 6]; vx <- a[, 4]; vy <- a[, 5]; r <- a[, 3]
start <- centers + cbind(r * (vx * cos(theta) - vy * sin(theta)),
                        r * (vx * sin(theta) + vy * cos(theta)))
end <- centers + cbind(r * (vx * cos(theta) + vy * sin(theta)),
                      r * (-vx * sin(theta) + vy * cos(theta)))
delta <- (end - start) * direction
closure_error <- sqrt(sum(colSums(delta)^2))
cc <- sweep(centers, 2, origin)
center_area <- sum((cc[, 1] * delta[, 2] - cc[, 2] * delta[, 1] -
                   direction * 2 * r^2 * theta) / 2)
# Calculate the tiny loop around a local origin to reduce cancellation.
loop_origin <- colMeans(h$xahull[c(4324, 4325, 4535), , drop = FALSE])
ls <- sweep(h$xahull[s[loop], , drop = FALSE], 2, loop_origin)
lt <- sweep(h$xahull[t[loop], , drop = FALSE], 2, loop_origin)
loop_area <- sum((ls[, 1] * lt[, 2] - ls[, 2] * lt[, 1] -
  direction[loop] * r[loop]^2 * (2 * theta[loop] - sin(2 * theta[loop]))) / 2)
bound <- read.csv(file.path(retry, 'retry_result.csv'))$convex_area
stopifnot(is.finite(area), area >= 0, area <= bound,
          closure_error < 1e-8 * max(1, sum(2 * r * theta)),
          abs(area - center_area) < 1e-8, abs(area - rest_area - loop_area) < 1e-8,
          abs(loop_area) < 1e-8)
audit <- data.frame(scale = '20m', quad_id = '413', basis = 'original_pca',
  alpha_value = .5, method = 'manual_arc_423_traversal_reversal',
  repaired_area = area, independent_center_area = center_area,
  isolated_loop_signed_area = loop_area, remaining_boundary_area = rest_area,
  closure_error = closure_error, convex_bound = bound,
  timestamp = format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'))
write.csv(audit, file.path(out, 'manual_repair_audit.csv'), row.names = FALSE)
write.csv(data.frame(arc = active, direction = direction,
  start_id = s, end_id = t, area_term = endpoint_terms),
  file.path(out, 'repaired_traversal.csv'), row.names = FALSE)
d <- read.csv(canonical, stringsAsFactors = FALSE)
old <- d
i <- which(d$quad_id == '413' & d$basis == 'original_pca')
stopifnot(length(i) == 1, is.na(d$alpha_area_all_pixels[i]))
d$alpha_area_all_pixels[i] <- area
d$status[i] <- 'ok'
d$source[i] <- 'computed_all_pixels_alpha_0_5_manual_arc_orientation_repair'
d$warning[i] <- paste('User-authorized manual repair: reversed traversal of arc 423',
  'in isolated near-degenerate three-arc loop; coordinates, curves and alpha unchanged.',
  'Original package area was NA; see manual_orientation_20m_413 audit.')
d$error[i] <- NA
d$run_timestamp[i] <- audit$timestamp
stopifnot(identical(d[-i, ], old[-i, ]))
cp <- read.csv(checkpoint, stringsAsFactors = FALSE)
old_cp <- cp; j <- which(cp$basis == 'original_pca')
cp[j, ] <- d[i, ]
stopifnot(identical(cp[-j, ], old_cp[-j, ]))
write.csv(d, canonical, row.names = FALSE)
write.csv(cp, checkpoint, row.names = FALSE)
coverage <- read.csv(coverage_path, stringsAsFactors = FALSE)
k <- which(coverage$scale == '20m' & coverage$basis == 'original_pca')
coverage$finite_area[k] <- sum(d$basis == 'original_pca' & d$status == 'ok')
coverage$unresolved_geometry[k] <- sum(d$basis == 'original_pca' &
  !d$status %in% c('ok', 'missing_raster', 'skipped_manual_excluded'))
# recovered_areas counts the earlier integral recoveries only; manual repair is
# identified separately in the canonical source field and this dedicated audit.
write.csv(coverage, coverage_path, row.names = FALSE)
unresolved <- read.csv(unresolved_path, stringsAsFactors = FALSE)
unresolved <- unresolved[!(unresolved$scale == '20m' & unresolved$quad_id == '413' &
                         unresolved$basis == 'original_pca'), ]
write.csv(unresolved, unresolved_path, row.names = FALSE)
status <- as.data.frame(table(d$scale, d$basis, d$status))
names(status) <- c('scale', 'basis', 'status', 'rows')
write.csv(status[status$rows > 0, ], status_path, row.names = FALSE)
saved <- read.csv(canonical, stringsAsFactors = FALSE)
stopifnot(identical(saved[-i, ], old[-i, ]),
          abs(saved$alpha_area_all_pixels[i] - area) < 1e-10,
          identical(unname(tools::md5sum(protected$path)), protected$md5),
          nrow(unresolved) == 4)
write.csv(data.frame(check = c('exact_vertex_balance', 'endpoint_circle_consistency',
  'independent_integral_agreement', 'closed_remaining_boundary', 'negligible_loop_area',
  'finite_nonnegative_within_convex_bound', 'other_canonical_rows_unchanged',
  'vector_normalized_checkpoint_unchanged', 'protected_references_unchanged',
  'saved_result_round_trip'), pass = TRUE), file.path(out, 'validation.csv'), row.names = FALSE)
print(audit, digits = 16)
cat('Manual repair validation: PASS\n')
