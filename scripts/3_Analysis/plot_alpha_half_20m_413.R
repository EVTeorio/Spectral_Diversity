# Plot the saved failed hull; do not recompute or modify its geometry.
input <- paste0('reports/tables/20261015_alpha_half_small_scales/',
                'retry_20m_413_original_20261015_161612/retry_hull.rds')
output <- 'reports/figures/20261015_alpha_half_20m_413'
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
h <- readRDS(input)
a <- h$arcs[h$arcs[, 3] > 0, , drop = FALSE]
ids <- c(4324L, 4535L)
problem <- which(a[, 7] %in% ids | a[, 8] %in% ids)
stopifnot(length(problem) == 3L)
colors <- c('#D55E00', '#0072B2', '#CC79A7')
arc_xy <- function(k) {
  z <- a[k, ]
  angle <- atan2(z[5], z[4]) + seq(z[6], -z[6], length.out = 161)
  cbind(z[1] + z[3] * cos(angle), z[2] + z[3] * sin(angle))
}
draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), oma = c(4, 0, 3, 0),
      family = 'sans', las = 1)
  plot(h$xahull, pch = 16, cex = .24, col = '#88939A60', asp = 1,
       xlab = 'PC1 score', ylab = 'PC2 score', main = 'A  |  Full saved hull')
  for (k in seq_len(nrow(a))) lines(arc_xy(k), col = '#263B49', lwd = .9)
  for (j in seq_along(problem)) lines(arc_xy(problem[j]), col = colors[j], lwd = 2)
  xy <- h$xahull[ids, , drop = FALSE]
  mid <- colMeans(xy)
  rect(mid[1] - .7, mid[2] - .7, mid[1] + .7, mid[2] + .7,
       border = '#D55E00', lwd = 1.5)
  text(mid[1] + 1.2, mid[2] + 2, 'Detail in B', adj = 0, col = '#A54400', cex = .85)
  legend('bottomleft', c('PCA points', 'Returned circular arcs', 'Problem location'),
         pch = c(16, NA, 0), lty = c(NA, 1, NA),
         col = c('#88939A', '#263B49', '#D55E00'), bty = 'n', cex = .8)
  plot(h$xahull, xlim = mid[1] + c(-.075, .075), ylim = mid[2] + c(-.065, .065),
       pch = 16, cex = .5, col = '#88939A80', asp = 1,
       xlab = 'PC1 score', ylab = 'PC2 score', main = 'B  |  Magnified problem boundary',
       xaxs = 'i', yaxs = 'i')
  for (k in seq_len(nrow(a))) lines(arc_xy(k), col = '#B4BEC4', lwd = 1)
  for (j in seq_along(problem)) {
    z <- arc_xy(problem[j])
    lines(z, col = colors[j], lwd = 2, lty = j)
    # Arrows follow the clockwise traversal used by the area checks.
    if (diff(range(z[, 1])) + diff(range(z[, 2])) > .001)
      arrows(z[65, 1], z[65, 2], z[95, 1], z[95, 2], length = .10,
             col = colors[j], lwd = 2)
  }
  points(xy, pch = 21, bg = 'white', col = '#111111', cex = 1.25, lwd = 1.5)
  text(xy[1, 1], xy[1, 2] - .009, '4324', cex = .9, adj = c(.5, 1))
  text(xy[2, 1], xy[2, 2] + .009, '4535 / 4325*', cex = .85, adj = c(.5, 0))
  legend('bottomleft', paste('Arc', problem), col = colors, lty = 1:3, lwd = 2,
         bty = 'n', cex = .8)
  mtext('20 m quadrat 413 | Original PCA | Alpha = 0.5', outer = TRUE,
        side = 3, line = .8, cex = 1.2, font = 2)
  mtext('Saved retry geometry; no repair applied. Arrows show boundary traversal directions.',
        outer = TRUE, side = 1, line = 1.1, cex = .85)
  mtext('*Vertices 4535 and 4325 are nearly coincident; two long arcs also overlap at this scale. Area remains NA.',
        outer = TRUE, side = 1, line = 2.5, cex = .8)
}
png(paste0(output, '.png'), width = 2400, height = 1350, res = 180)
draw()
dev.off()
pdf(paste0(output, '.pdf'), width = 13.333, height = 7.5)
draw()
dev.off()
cat('Created ', output, '.png and .pdf\n', sep = '')
