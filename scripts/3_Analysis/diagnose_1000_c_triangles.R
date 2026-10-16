.libPaths(c('C:/Users/PaintRock/AppData/Local/R/win-library/4.2', .libPaths()))
root <- 'reports/tables/20261015_alpha_half_small_scales/1000_c_numeric_geometry'
dir.create(root, showWarnings = FALSE)
helper <- new.env()
e <- parse('scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R')
eval(e[-length(e)], helper)
x <- helper$read_masked_spectra('Quad_Spectra/10m_smooth_5nm/1000_c')
p <- unique(helper$project_scores(x, readRDS(helper$PCA_RDS))[, 1:2, drop = FALSE])
ns <- asNamespace('alphahull')
tm <- get('tri.mesh', ns)(p[, 1], p[, 2])
tr <- get('triangles', ns)(tm)
saveRDS(list(points = p, triangulation = tm, triangles = tr), file.path(root, 'triangulation.rds'))
singular <- svd(sweep(p, 2, colMeans(p)), nu = 0, nv = 0)$d
global <- data.frame(retained_pixels = nrow(x), unique_points = nrow(p),
  singular_value_1 = singular[1], singular_value_2 = singular[2],
  singular_ratio = singular[2] / singular[1], triangulation_points_unchanged =
    identical(as.numeric(tm$x), as.numeric(p[, 1])) && identical(as.numeric(tm$y), as.numeric(p[, 2])))
write.csv(global, file.path(root, 'global_collinearity.csv'), row.names = FALSE)
records <- vector('list', nrow(tr))
centers <- matrix(NA_real_, nrow(tr), 2)
for (i in seq_len(nrow(tr))) {
  xy <- p[tr[i, 1:3], , drop = FALSE]
  u <- xy[2, ] - xy[1, ]; v <- xy[3, ] - xy[1, ]
  determinant <- u[1] * v[2] - u[2] * v[1]
  lens <- c(sum(u^2), sum(v^2), sum((v - u)^2))
  std <- get('circum', ns)(xy[, 1], xy[, 2])
  if (determinant != 0) {
    centers[i, ] <- xy[1, ] + c(sum(u^2) * v[2] - sum(v^2) * u[2],
      u[1] * sum(v^2) - v[1] * sum(u^2)) / (2 * determinant)
  }
  rad <- sqrt(rowSums(sweep(xy, 2, centers[i, ])^2))
  records[[i]] <- data.frame(triangle = i, node1 = tr[i, 1], node2 = tr[i, 2],
    node3 = tr[i, 3], determinant = determinant, area = abs(determinant) / 2,
    height_longest_edge = abs(determinant) / sqrt(max(lens)),
    scaled_determinant = abs(determinant) / max(lens),
    standard_x = std$x, standard_y = std$y,
    stable_x = centers[i, 1], stable_y = centers[i, 2],
    stable_radius = mean(rad), stable_radius_spread = diff(range(rad)),
    center_discrepancy = sqrt(sum((c(std$x, std$y) - centers[i, ])^2)))
}
d <- do.call(rbind, records)
write.csv(d, file.path(root, 'triangle_diagnostics.csv'), row.names = FALSE)
bad <- d[!is.finite(d$standard_x) | !is.finite(d$standard_y), ]
write.csv(bad, file.path(root, 'nonfinite_standard_circumcenters.csv'), row.names = FALSE)
saveRDS(centers, file.path(root, 'stable_circumcenters.rds'))
print(global, digits = 16)
cat('Triangles:', nrow(d), '; exact zero-area:', sum(d$determinant == 0),
    '; nonfinite standard centers:', nrow(bad), '\n')
print(bad, digits = 16)
