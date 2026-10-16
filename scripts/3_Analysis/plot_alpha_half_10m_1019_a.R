# Visualize the saved, unrepaired 1019_a hull, without rerunning the geometry.
root <- 'reports/tables/20261015_alpha_half_small_scales/diagnose_four_10m'
h <- readRDS(file.path(root, '1019_a_hull.rds'))
bad <- read.csv(file.path(root, '1019_a_unbalanced_vertices.csv'))
active <- which(h$arcs[, 3] > 0)
a <- h$arcs[active, , drop = FALSE]
candidate <- which(active == 487L)
stopifnot(length(candidate) == 1, setequal(a[candidate, 7:8], bad$vertex))
incident <- which(a[, 7] %in% bad$vertex | a[, 8] %in% bad$vertex)
neighbors <- setdiff(incident, candidate)
arc_xy <- function(k) {
  z <- a[k, ]
  phi <- atan2(z[5], z[4]) + seq(z[6], -z[6], length.out = 241)
  cbind(z[1] + z[3] * cos(phi), z[2] + z[3] * sin(phi))
}
xy <- h$xahull[bad$vertex, , drop = FALSE]
local <- do.call(rbind, lapply(incident, arc_xy))
center <- colMeans(apply(local, 2, range))
half <- max(diff(range(local[, 1])), diff(range(local[, 2]))) / 2 + .35
zoom_x <- center[1] + c(-half, half)
zoom_y <- center[2] + c(-half, half)
orange <- '#D55E00'; teal <- '#007F8B'
out <- 'reports/figures/20261015_alpha_half_10m_1019_a.png'
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
png(out, width = 2400, height = 1350, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), oma = c(4, 0, 3, 0),
    family = 'sans', las = 1)
plot(h$xahull, pch = 16, cex = .33, col = '#88939A75', asp = 1,
     xlab = 'PC1 score', ylab = 'PC2 score', main = 'A  |  Full saved hull')
for (k in seq_len(nrow(a))) lines(arc_xy(k), col = '#263B49', lwd = 1)
lines(arc_xy(candidate), col = orange, lwd = 2.5)
rect(zoom_x[1], zoom_y[1], zoom_x[2], zoom_y[2], border = orange, lwd = 1.5)
text(zoom_x[1] - 1, zoom_y[2] + 1.4, 'Detail in B', adj = 1, col = '#A54400', cex = .9)
legend('bottomleft', c('PCA points', 'Returned hull boundary', 'Arc 487 / detail location'),
       pch = c(16, NA, NA), lty = c(NA, 1, 1),
       col = c('#88939A', '#263B49', orange), bty = 'n', cex = .8)
plot(h$xahull, xlim = zoom_x, ylim = zoom_y, asp = 1, pch = 16, cex = .8,
     col = '#88939A85', xaxs = 'i', yaxs = 'i',
     xlab = 'PC1 score', ylab = 'PC2 score', main = 'B  |  Local direction inconsistency')
for (k in seq_len(nrow(a))) lines(arc_xy(k), col = '#B4BEC4', lwd = 1)
for (k in incident) {
  z <- arc_xy(k)
  color <- if (k == candidate) orange else teal
  lines(z, col = color, lwd = 2.5)
  arrows(z[90, 1], z[90, 2], z[150, 1], z[150, 2],
         length = .12, col = color, lwd = 2)
}
points(xy, pch = 21, bg = 'white', col = '#111111', cex = 1.3, lwd = 1.5)
text(xy[1, 1] + .10, xy[1, 2] + .09, paste('Vertex', bad$vertex[1]), adj = c(0, 0), cex = .85)
text(xy[2, 1] + .12, xy[2, 2] - .09, paste('Vertex', bad$vertex[2]), adj = c(0, 1), cex = .85)
legend('bottomleft', c('Arc 487 (candidate for reversal)', 'Adjacent arcs', 'Imbalanced vertices'),
       col = c(orange, teal, '#111111'), lty = c(1, 1, NA), pch = c(NA, NA, 21),
       pt.bg = 'white', lwd = 2, bty = 'n', cex = .75)
mtext('10 m quadrat 1019_a | Original PCA | Alpha = 0.5', outer = TRUE,
      side = 3, line = .8, cex = 1.2, font = 2)
mtext('Arrows show the current traversal directions. The flagged arc belongs to a large, 292-arc boundary component.',
      outer = TRUE, side = 1, line = 1.1, cex = .83)
mtext('Saved geometry is shown unchanged. No repair has been applied; the recorded area remains NA.',
      outer = TRUE, side = 1, line = 2.5, cex = .83)
dev.off()
cat('Created ', out, '\n', sep = '')
