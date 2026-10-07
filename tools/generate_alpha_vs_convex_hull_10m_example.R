PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

OUT_PATH <- file.path(PROJECT_DIR, "Documents/Tables and Figures/20_alpha_hull_vs_convex_hull_10m_example.png")
SCORE_CSV <- file.path(PROJECT_DIR, "reports/tables/figure_sources/20_alpha_hull_vs_convex_hull_10m_example_scores.csv")
PCA_RDS <- file.path(PROJECT_DIR, "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds")
SPEC_DIR <- file.path(PROJECT_DIR, "Quad_Spectra/10m_smooth_5nm")

dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(SCORE_CSV), recursive = TRUE, showWarnings = FALSE)

samples <- 136L
lines <- 136L
bands <- 121L
n_pixels <- samples * lines
shadow_threshold <- 0.0305476
mask_band_index <- 34L

plot_cols <- list(
  axis = rgb(0.43, 0.48, 0.52),
  text = rgb(0.10, 0.13, 0.16),
  subtext = rgb(0.37, 0.42, 0.46),
  points = rgb(0.12, 0.25, 0.36, 0.62),
  convex_fill = rgb(0.70, 0.74, 0.76, 0.58),
  convex_border = rgb(0.31, 0.36, 0.39),
  alpha_fill = rgb(0.10, 0.49, 0.56, 0.68),
  alpha_border = rgb(0.03, 0.33, 0.38),
  legend_bg = rgb(1, 1, 1, 0.72)
)

read_envi_cube <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  vals <- readBin(con, what = "numeric", n = n_pixels * bands, size = 4, endian = "little")
  matrix(vals, nrow = n_pixels, ncol = bands)
}

clean_project_scores <- function(tile_name, pca_object) {
  x <- read_envi_cube(file.path(SPEC_DIR, tile_name))
  keep <- is.finite(x[, mask_band_index]) & x[, mask_band_index] > shadow_threshold
  x <- x[keep, , drop = FALSE]
  x <- x[complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]

  norms <- sqrt(rowSums(x^2))
  x <- x[is.finite(norms) & norms > 0, , drop = FALSE]
  norms <- norms[is.finite(norms) & norms > 0]
  x <- sweep(x, 1, norms, "/")

  x <- sweep(x, 2, pca_object$center, "-")
  x <- sweep(x, 2, pca_object$scale, "/")
  scores <- x %*% pca_object$pca$rotation[, 1:2, drop = FALSE]
  colnames(scores) <- c("PC1", "PC2")
  scores
}

polygon_area <- function(x, y) {
  if (length(x) < 3) return(NA_real_)
  0.5 * abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])))
}

tight_boundary <- function(scores, bins = 72L) {
  center <- colMeans(scores)
  dx <- scores[, 1] - center[1]
  dy <- scores[, 2] - center[2]
  theta <- atan2(dy, dx)
  theta[theta < 0] <- theta[theta < 0] + 2 * pi
  radius <- sqrt(dx^2 + dy^2)
  bin <- pmin(bins, floor(theta / (2 * pi) * bins) + 1L)

  idx <- tapply(seq_along(radius), bin, function(ii) ii[which.max(radius[ii])])
  idx <- as.integer(idx[!is.na(idx)])
  pts <- scores[idx, , drop = FALSE]
  ord <- order(atan2(pts[, 2] - center[2], pts[, 1] - center[1]))
  pts <- pts[ord, , drop = FALSE]
  pts[!duplicated(paste(round(pts[, 1], 8), round(pts[, 2], 8))), , drop = FALSE]
}

make_example <- function() {
  pca_object <- readRDS(PCA_RDS)
  candidates <- c("1000_a", "1001_a", "1002_a", "604_a", "500_a", "700_a", "0_a", "300_a")
  candidates <- candidates[file.exists(file.path(SPEC_DIR, candidates))]

  best <- NULL
  for (tile in candidates) {
    scores <- clean_project_scores(tile, pca_object)
    if (nrow(scores) < 250) next
    set.seed(1000 + sum(utf8ToInt(tile)))
    plot_scores <- scores[sample.int(nrow(scores), min(1200L, nrow(scores))), , drop = FALSE]
    convex <- plot_scores[chull(plot_scores), , drop = FALSE]
    tight <- tight_boundary(plot_scores)
    ratio <- polygon_area(tight[, 1], tight[, 2]) / polygon_area(convex[, 1], convex[, 2])
    if (is.null(best) || (!is.na(ratio) && ratio < best$ratio)) {
      best <- list(tile = tile, scores = plot_scores, tight = tight, convex = convex, ratio = ratio, n_retained = nrow(scores))
    }
  }
  best
}

draw_legend <- function(xlim, ylim, pad_x, pad_y) {
  usr <- par("usr")
  lx <- usr[1] + 0.69 * diff(usr[1:2])
  ly <- usr[4] - 0.09 * diff(usr[3:4])
  box_w <- 0.32 * diff(usr[1:2])
  box_h <- 0.23 * diff(usr[3:4])
  rect(lx, ly - box_h, lx + box_w, ly, col = plot_cols$legend_bg, border = NA)

  y1 <- ly - 0.06 * diff(usr[3:4])
  y_step <- 0.055 * diff(usr[3:4])
  sw <- 0.045 * diff(usr[1:2])
  sh <- 0.025 * diff(usr[3:4])
  tx <- lx + 0.075 * diff(usr[1:2])

  points(lx + sw / 2, y1, pch = 21, bg = plot_cols$points, col = "white", cex = 0.85, lwd = 0.5)
  text(tx, y1, "Pixel scores", adj = c(0, 0.5), cex = 0.78, col = plot_cols$text)

  rect(lx, y1 - y_step - sh / 2, lx + sw, y1 - y_step + sh / 2, col = plot_cols$convex_fill, border = plot_cols$convex_border, lwd = 1.3)
  text(tx, y1 - y_step, "Convex-hull area", adj = c(0, 0.5), cex = 0.78, col = plot_cols$text)

  rect(lx, y1 - 2 * y_step - sh / 2, lx + sw, y1 - 2 * y_step + sh / 2, col = plot_cols$alpha_fill, border = plot_cols$alpha_border, lwd = 1.3)
  text(tx, y1 - 2 * y_step, "Alpha-hull area", adj = c(0, 0.5), cex = 0.78, col = plot_cols$text)
}

draw_panel <- function(scores, convex, tight, main, sub, mode, show_legend = FALSE) {
  xlim <- range(scores[, 1])
  ylim <- range(scores[, 2])
  pad_x <- diff(xlim) * 0.18
  pad_y <- diff(ylim) * 0.18
  plot(
    scores[, 1], scores[, 2],
    xlim = xlim + c(-pad_x, if (show_legend) pad_x * 2.25 else pad_x),
    ylim = ylim + c(-pad_y, pad_y),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  axis_col <- plot_cols$axis
  lines(c(xlim[1] - pad_x * 0.55, xlim[2] + pad_x * 0.45), rep(ylim[1] - pad_y * 0.25, 2), col = axis_col, lwd = 2)
  lines(rep(xlim[1] - pad_x * 0.45, 2), c(ylim[1] - pad_y * 0.25, ylim[2] + pad_y * 0.45), col = axis_col, lwd = 2)
  mtext("PC1", side = 1, line = 1.8, col = axis_col, cex = 0.95)
  mtext("PC2", side = 2, line = 1.8, col = axis_col, cex = 0.95)

  title(main = main, col.main = plot_cols$text, cex.main = 1.30, font.main = 2, line = 1.4)
  mtext(sub, side = 3, line = 0.1, col = plot_cols$subtext, cex = 0.80)

  if (mode == "convex") {
    polygon(convex[, 1], convex[, 2], col = plot_cols$convex_fill, border = plot_cols$convex_border, lwd = 2.5)
    polygon(tight[, 1], tight[, 2], col = plot_cols$alpha_fill, border = plot_cols$alpha_border, lwd = 2.2)
  } else {
    polygon(tight[, 1], tight[, 2], col = plot_cols$alpha_fill, border = plot_cols$alpha_border, lwd = 2.5)
  }

  points(scores[, 1], scores[, 2], pch = 21, bg = plot_cols$points, col = "white", cex = 0.48, lwd = 0.45)
  if (show_legend) draw_legend(xlim, ylim, pad_x, pad_y)
}

example <- make_example()
if (is.null(example)) {
  stop("No suitable 10 m quadrat found for plotting.", call. = FALSE)
}

write.csv(
  data.frame(PC1 = example$scores[, 1], PC2 = example$scores[, 2]),
  SCORE_CSV,
  row.names = FALSE
)

png(OUT_PATH, width = 2400, height = 1400, res = 200, bg = "transparent")
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)

par(family = "sans", fg = plot_cols$text, col.axis = plot_cols$subtext, mar = c(4.3, 4.6, 5.2, 1.2), oma = c(1.1, 0.5, 4.4, 0.5), mfrow = c(1, 2), xaxs = "i", yaxs = "i")

draw_panel(
  example$scores, example$convex, example$tight,
  "A. Convex hull", "Outer envelope around 10 m pixel scores", "convex"
)
draw_panel(
  example$scores, example$convex, example$tight,
  "B. Alpha-hull area", "Boundary following occupied PC1-PC2 space", "alpha", show_legend = TRUE
)

mtext("Observed 10 m quadrat example in vector-normalized PCA spectral space", outer = TRUE, side = 3, line = 2.0, adj = 0.03, cex = 1.40, font = 2, col = plot_cols$text)
mtext(
  sprintf("Example tile: %s; retained illuminated pixels: %s; plotted sample: %s pixels", example$tile, format(example$n_retained, big.mark = ","), format(nrow(example$scores), big.mark = ",")),
  outer = TRUE, side = 3, line = 0.7, adj = 0.03, cex = 0.90, col = plot_cols$subtext
)

cat("Created", OUT_PATH, "\n")
cat("Selected tile:", example$tile, "\n")
cat("Retained illuminated pixels:", example$n_retained, "\n")
cat("Tight/convex area ratio:", round(example$ratio, 3), "\n")
