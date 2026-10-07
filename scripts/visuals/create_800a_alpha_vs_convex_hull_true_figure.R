PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

.libPaths(c(
  "C:/Users/PaintRock/AppData/Local/R/win-library/4.2",
  "C:/Program Files/R/R-4.2.3/library",
  .libPaths()
))

if (!requireNamespace("alphahull", quietly = TRUE)) {
  stop("Package 'alphahull' is required in the R 4.2 library.", call. = FALSE)
}

OUT_PATH <- file.path(PROJECT_DIR, "Documents/Tables and Figures/25_alpha_hull_vs_convex_hull_10m_tile_800_a_true_calculated.png")
META_CSV <- file.path(PROJECT_DIR, "reports/tables/figure_sources/25_alpha_hull_vs_convex_hull_10m_tile_800_a_true_calculated_metadata.csv")
PCA_RDS <- file.path(PROJECT_DIR, "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds")
SUMMARY_CSV <- file.path(PROJECT_DIR, "Quad_Values/10m_spectral_heterogeneity_smooth_masked_5nm_summary.csv")
SPEC_PATH <- file.path(PROJECT_DIR, "Quad_Spectra/10m_smooth_5nm/800_a")

dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(META_CSV), recursive = TRUE, showWarnings = FALSE)

SAMPLES <- 136L
LINES <- 136L
BANDS <- 121L
N_PIXELS <- SAMPLES * LINES
SHADOW_THRESHOLD <- 0.0305476
MASK_BAND_INDEX <- 34L
ALPHA_VALUE <- 1

plot_cols <- list(
  axis = rgb(0.35, 0.39, 0.42),
  text = rgb(0.09, 0.11, 0.13),
  subtext = rgb(0.36, 0.39, 0.42),
  points = rgb(0.10, 0.20, 0.28, 0.40),
  convex_fill = rgb(0.88, 0.48, 0.18, 0.26),
  convex_border = rgb(0.73, 0.27, 0.08),
  alpha_border = rgb(0.02, 0.39, 0.45),
  alpha_arc = rgb(0.02, 0.54, 0.61),
  label_bg = rgb(1, 1, 1, 0.82)
)

read_envi_cube <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  vals <- readBin(con, what = "numeric", n = N_PIXELS * BANDS, size = 4, endian = "little")
  if (length(vals) != N_PIXELS * BANDS) {
    stop("Unexpected raster size for 800_a.", call. = FALSE)
  }
  matrix(vals, nrow = N_PIXELS, ncol = BANDS)
}

project_standardized_pca_scores <- function(tile_path, pca_object) {
  x <- read_envi_cube(tile_path)
  keep <- is.finite(x[, MASK_BAND_INDEX]) & x[, MASK_BAND_INDEX] > SHADOW_THRESHOLD
  x <- x[keep, , drop = FALSE]
  x <- x[complete.cases(x), , drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[apply(x, 1, function(row) all(is.finite(row))), , drop = FALSE]

  norms <- sqrt(rowSums(x^2))
  keep_norm <- is.finite(norms) & norms > 0
  x <- x[keep_norm, , drop = FALSE]
  norms <- norms[keep_norm]
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

draw_base_panel <- function(points_unique, main_title, subtitle) {
  xlim <- range(points_unique[, 1])
  ylim <- range(points_unique[, 2])
  pad_x <- diff(xlim) * 0.12
  pad_y <- diff(ylim) * 0.12

  plot(
    points_unique[, 1], points_unique[, 2],
    xlim = xlim + c(-pad_x, pad_x),
    ylim = ylim + c(-pad_y, pad_y),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  axis(1, col = plot_cols$axis, col.axis = plot_cols$subtext, lwd = 1.4, cex.axis = 0.85)
  axis(2, col = plot_cols$axis, col.axis = plot_cols$subtext, lwd = 1.4, cex.axis = 0.85)
  box(col = rgb(0.70, 0.72, 0.74), lwd = 1.2)
  mtext("PC1", side = 1, line = 2.3, col = plot_cols$axis, cex = 0.95)
  mtext("PC2", side = 2, line = 2.3, col = plot_cols$axis, cex = 0.95)
  title(main = main_title, col.main = plot_cols$text, cex.main = 1.18, font.main = 2, line = 1.4)
  mtext(subtitle, side = 3, line = 0.10, col = plot_cols$subtext, cex = 0.78)
}

draw_note <- function(label) {
  usr <- par("usr")
  x_left <- usr[1] + 0.035 * diff(usr[1:2])
  x_right <- usr[1] + 0.60 * diff(usr[1:2])
  y_top <- usr[4] - 0.035 * diff(usr[3:4])
  y_bottom <- usr[4] - 0.155 * diff(usr[3:4])
  rect(x_left, y_bottom, x_right, y_top, col = plot_cols$label_bg, border = NA)
  text(
    x_left + 0.018 * diff(usr[1:2]),
    y_top - 0.030 * diff(usr[3:4]),
    label,
    adj = c(0, 1),
    cex = 0.78,
    col = plot_cols$text
  )
}

draw_alpha_hull <- function(alpha_hull) {
  arcs <- alpha_hull$arcs
  curved_arcs <- which(arcs[, 3] > 0)
  point_arcs <- which(arcs[, 3] == 0)

  if (length(curved_arcs) > 0) {
    for (arc_index in curved_arcs) {
      alphahull::arc(
        arcs[arc_index, 1:2],
        arcs[arc_index, 3],
        arcs[arc_index, 4:5],
        arcs[arc_index, 6],
        col = plot_cols$alpha_arc,
        lwd = 2.1
      )
    }
  }

  if (length(point_arcs) > 0) {
    points(
      arcs[point_arcs, 1],
      arcs[point_arcs, 2],
      pch = 16,
      col = plot_cols$alpha_border,
      cex = 0.55
    )
  }
}

pca_object <- readRDS(PCA_RDS)
scores <- project_standardized_pca_scores(SPEC_PATH, pca_object)
points_unique <- unique(round(scores, 10))
points_unique <- as.matrix(points_unique)
colnames(points_unique) <- c("PC1", "PC2")

convex_idx <- chull(points_unique[, 1], points_unique[, 2])
convex_polygon <- points_unique[convex_idx, , drop = FALSE]
convex_area <- polygon_area(convex_polygon[, 1], convex_polygon[, 2])

alpha_hull <- alphahull::ahull(points_unique[, 1], points_unique[, 2], alpha = ALPHA_VALUE)
alpha_area <- alphahull::areaahull(alpha_hull, timeout = 10)

summary_data <- read.csv(SUMMARY_CSV, check.names = FALSE)
summary_row <- summary_data[summary_data$quad_id == "800_a", , drop = FALSE]
if (nrow(summary_row) != 1) {
  stop("Expected exactly one row for 800_a in the 10 m spectral heterogeneity summary.", call. = FALSE)
}

metadata <- data.frame(
  quad_id = "800_a",
  pca_basis = "vector-normalized standardized PCA",
  alpha = ALPHA_VALUE,
  retained_pixels = nrow(scores),
  unique_pc1_pc2_points = nrow(points_unique),
  recalculated_convex_hull_area = convex_area,
  table_convex_hull_area = summary_row$standardized_PCA_pca_convex_hull_area,
  recalculated_alpha_hull_area = alpha_area,
  table_alpha_hull_area = summary_row$standardized_PCA_alpha_hull_area,
  table_alpha_hull_method = summary_row$standardized_PCA_alpha_hull_method,
  table_alpha_hull_n_points = summary_row$standardized_PCA_alpha_hull_n_points
)
write.csv(metadata, META_CSV, row.names = FALSE)

png(OUT_PATH, width = 2600, height = 1450, res = 220, bg = "white")
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)

par(
  family = "sans",
  fg = plot_cols$text,
  col.axis = plot_cols$subtext,
  mar = c(4.7, 4.9, 5.0, 1.1),
  oma = c(1.1, 0.2, 4.4, 0.2),
  mfrow = c(1, 2),
  xaxs = "i",
  yaxs = "i"
)

draw_base_panel(points_unique, "A. Convex hull", "All occupied PC1-PC2 scores enclosed by one convex boundary")
polygon(
  convex_polygon[, 1],
  convex_polygon[, 2],
  col = plot_cols$convex_fill,
  border = plot_cols$convex_border,
  lwd = 2.4
)
points(points_unique[, 1], points_unique[, 2], pch = 16, col = plot_cols$points, cex = 0.30)
draw_note(sprintf("Convex area: %.1f\nSingle continuous outer envelope", convex_area))

draw_base_panel(points_unique, "B. Alpha hull", "Actual alphahull boundary used for the calculated area")
points(points_unique[, 1], points_unique[, 2], pch = 16, col = plot_cols$points, cex = 0.30)
draw_alpha_hull(alpha_hull)
draw_note(sprintf("Alpha area: %.1f\nalpha = %s; %s unique points", alpha_area, ALPHA_VALUE, format(nrow(points_unique), big.mark = ",")))

mtext(
  "10 m quadrat 800_a in vector-normalized PCA spectral space",
  outer = TRUE,
  side = 3,
  line = 2.2,
  adj = 0.03,
  cex = 1.28,
  font = 2,
  col = plot_cols$text
)
mtext(
  sprintf(
    "Retained illuminated pixels: %s; alpha-hull method in metric table: %s; table alpha area: %.1f",
    format(nrow(scores), big.mark = ","),
    summary_row$standardized_PCA_alpha_hull_method,
    summary_row$standardized_PCA_alpha_hull_area
  ),
  outer = TRUE,
  side = 3,
  line = 0.8,
  adj = 0.03,
  cex = 0.85,
  col = plot_cols$subtext
)

cat("Created", OUT_PATH, "\n")
cat("Created", META_CSV, "\n")
cat("Retained pixels:", nrow(scores), "\n")
cat("Unique PC1-PC2 points:", nrow(points_unique), "\n")
cat("Convex area:", convex_area, "\n")
cat("Alpha area:", alpha_area, "\n")
