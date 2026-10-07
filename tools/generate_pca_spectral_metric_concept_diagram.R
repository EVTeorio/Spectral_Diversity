PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

OUT_PATH <- file.path(PROJECT_DIR, "Documents/Tables and Figures/21_pca_spectral_metric_concept_diagram.png")
dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)

cols <- list(
  text = rgb(0.10, 0.13, 0.16),
  subtext = rgb(0.37, 0.42, 0.46),
  axis = rgb(0.45, 0.50, 0.54),
  point_fill = rgb(0.11, 0.32, 0.46, 0.82),
  point_border = "white",
  pairwise = rgb(0.36, 0.39, 0.42, 0.52),
  centroid = rgb(0.78, 0.10, 0.10),
  centroid_line = rgb(0.78, 0.10, 0.10, 0.55),
  hull_fill = rgb(0.10, 0.49, 0.56, 0.62),
  hull_border = rgb(0.03, 0.33, 0.38)
)

points_xy <- rbind(
  c(-1.45, -0.25),
  c(-1.15,  0.42),
  c(-0.82, -0.88),
  c(-0.48,  0.12),
  c(-0.05, -0.48),
  c( 0.18,  0.64),
  c( 0.58, -0.08),
  c( 0.92,  0.56),
  c( 1.15, -0.54),
  c( 1.55,  0.02)
)
colnames(points_xy) <- c("PC1", "PC2")

centroid <- colMeans(points_xy)
alpha_outline <- rbind(
  c(-1.58, -0.32),
  c(-1.18,  0.56),
  c(-0.60,  0.32),
  c( 0.16,  0.80),
  c( 0.42,  0.28),
  c( 0.95,  0.66),
  c( 1.70,  0.05),
  c( 1.22, -0.64),
  c( 0.55, -0.22),
  c(-0.02, -0.62),
  c(-0.82, -1.02),
  c(-0.62, -0.18),
  c(-1.16, -0.78)
)

draw_axes <- function() {
  plot(
    points_xy[, 1], points_xy[, 2],
    xlim = c(-1.85, 1.90),
    ylim = c(-1.22, 1.05),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  usr <- par("usr")
  lines(c(usr[1] + 0.10, usr[2] - 0.10), c(usr[3] + 0.16, usr[3] + 0.16), col = cols$axis, lwd = 2)
  lines(c(usr[1] + 0.22, usr[1] + 0.22), c(usr[3] + 0.08, usr[4] - 0.08), col = cols$axis, lwd = 2)
  mtext("PC1", side = 1, line = 1.65, col = cols$axis, cex = 0.90)
  mtext("PC2", side = 2, line = 1.85, col = cols$axis, cex = 0.90)
}

draw_points <- function() {
  points(points_xy[, 1], points_xy[, 2], pch = 21, bg = cols$point_fill, col = cols$point_border, cex = 1.75, lwd = 0.9)
}

panel_title <- function(main, sub) {
  title(main = main, col.main = cols$text, cex.main = 1.15, font.main = 2, line = 1.25)
  mtext(sub, side = 3, line = 0.00, col = cols$subtext, cex = 0.68)
}

draw_rao_panel <- function() {
  draw_axes()
  pairs <- combn(seq_len(nrow(points_xy)), 2)
  for (i in seq_len(ncol(pairs))) {
    p <- pairs[, i]
    lines(points_xy[p, 1], points_xy[p, 2], col = cols$pairwise, lwd = 1.15)
  }
  draw_points()
  panel_title("A. Spectral Rao's Q", "Mean pairwise dissimilarity among pixels")
}

draw_mean_distance_panel <- function() {
  draw_axes()
  for (i in seq_len(nrow(points_xy))) {
    lines(c(centroid[1], points_xy[i, 1]), c(centroid[2], points_xy[i, 2]), col = cols$centroid_line, lwd = 1.65)
  }
  draw_points()
  points(centroid[1], centroid[2], pch = 21, bg = cols$centroid, col = "white", cex = 2.25, lwd = 1.05)
  text(centroid[1] + 0.19, centroid[2] + 0.12, "centroid", col = cols$centroid, cex = 0.76, font = 2, adj = 0)
  panel_title("B. Mean Euclidean distance", "Average distance from pixels to the centroid")
}

draw_alpha_panel <- function() {
  draw_axes()
  polygon(alpha_outline[, 1], alpha_outline[, 2], col = cols$hull_fill, border = cols$hull_border, lwd = 2.8)
  draw_points()
  panel_title("C. Alpha-hull area", "Occupied area in PC1-PC2 spectral space")
}

png(OUT_PATH, width = 3000, height = 1200, res = 220, bg = "transparent")
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)

par(
  family = "sans",
  fg = cols$text,
  col.axis = cols$subtext,
  mfrow = c(1, 3),
  mar = c(4.5, 4.7, 4.9, 1.1),
  oma = c(0.4, 0.5, 3.7, 0.5)
)

draw_rao_panel()
draw_mean_distance_panel()
draw_alpha_panel()

mtext(
  "Conceptual representation of PCA-based spectral heterogeneity metrics",
  outer = TRUE,
  side = 3,
  line = 1.4,
  adj = 0.03,
  cex = 1.28,
  font = 2,
  col = cols$text
)
mtext(
  "Each point represents a retained illuminated pixel projected into spectral PCA space.",
  outer = TRUE,
  side = 3,
  line = 0.20,
  adj = 0.03,
  cex = 0.78,
  col = cols$subtext
)

cat("Created", OUT_PATH, "\n")
