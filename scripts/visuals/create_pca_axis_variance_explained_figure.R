PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

ORIGINAL_VARIANCE_CSV <- file.path(PROJECT_DIR, "Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv")
VECTOR_VARIANCE_CSV <- file.path(PROJECT_DIR, "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.csv")
OUT_PATH <- file.path(PROJECT_DIR, "Documents/Tables and Figures/26_pca_axis_variance_explained_original_vs_vector_normalized.png")
SOURCE_CSV <- file.path(PROJECT_DIR, "reports/tables/figure_sources/26_pca_axis_variance_explained_original_vs_vector_normalized.csv")

dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(SOURCE_CSV), recursive = TRUE, showWarnings = FALSE)

MAX_AXIS <- 10L

read_variance_table <- function(path, basis_id, basis_label) {
  x <- read.csv(path, check.names = FALSE)
  x <- x[x$pc_axis <= MAX_AXIS, , drop = FALSE]
  x$basis_id <- basis_id
  x$basis_label <- basis_label
  x$pct_variance_percent <- x$pct_variance * 100
  x$cumulative_pct_variance_percent <- x$cumulative_pct_variance * 100
  x
}

plot_cols <- list(
  text = rgb(0.08, 0.10, 0.12),
  subtext = rgb(0.34, 0.38, 0.42),
  axis = rgb(0.45, 0.49, 0.52),
  original = rgb(0.78, 0.31, 0.11),
  original_light = rgb(0.78, 0.31, 0.11, 0.34),
  vector = rgb(0.02, 0.43, 0.50),
  vector_light = rgb(0.02, 0.43, 0.50, 0.34),
  cumulative = rgb(0.10, 0.13, 0.16)
)

variance_data <- rbind(
  read_variance_table(ORIGINAL_VARIANCE_CSV, "original_pca", "Original PCA"),
  read_variance_table(VECTOR_VARIANCE_CSV, "vector_normalized_pca", "Vector-normalized PCA")
)

write.csv(variance_data, SOURCE_CSV, row.names = FALSE)

draw_panel <- function(panel_data, bar_col, light_col, title_text) {
  ymax <- 105
  bar_positions <- barplot(
    panel_data$pct_variance_percent,
    names.arg = rep("", nrow(panel_data)),
    ylim = c(0, ymax),
    col = ifelse(panel_data$pc_axis <= 2, bar_col, light_col),
    border = NA,
    axes = FALSE,
    main = "",
    xlab = "",
    ylab = ""
  )
  axis(2, las = 1, col = plot_cols$axis, col.axis = plot_cols$subtext, lwd = 1.2)
  axis(1, at = bar_positions, labels = paste0("PC", panel_data$pc_axis), las = 2, cex.axis = 0.72, col = plot_cols$axis, col.axis = plot_cols$subtext, lwd = 0)
  box(col = rgb(0.78, 0.80, 0.82), lwd = 1)
  grid(nx = NA, ny = NULL, col = rgb(0.90, 0.91, 0.92), lty = 1)

  lines(
    bar_positions,
    panel_data$cumulative_pct_variance_percent,
    col = plot_cols$cumulative,
    lwd = 2.1,
    type = "b",
    pch = 16,
    cex = 0.85
  )

  title(main = title_text, col.main = plot_cols$text, cex.main = 1.15, font.main = 2, line = 1.1)
  mtext("PCA axis", side = 1, line = 3.6, col = plot_cols$axis, cex = 0.92)
  mtext("Variance explained (%)", side = 2, line = 3.1, col = plot_cols$axis, cex = 0.92)

  label_axes <- panel_data$pc_axis <= 3
  text(
    bar_positions[label_axes],
    panel_data$pct_variance_percent[label_axes] + ymax * 0.025,
    sprintf("%.1f%%", panel_data$pct_variance_percent[label_axes]),
    col = plot_cols$text,
    cex = 0.78
  )

  usr <- par("usr")
  legend_x <- usr[1] + 0.54 * diff(usr[1:2])
  legend_y <- usr[3] + 0.26 * diff(usr[3:4])
  legend(
    legend_x,
    legend_y,
    legend = c("Axis variance", "Cumulative variance"),
    fill = c(bar_col, NA),
    border = c(NA, NA),
    lty = c(NA, 1),
    lwd = c(NA, 2.1),
    pch = c(NA, 16),
    col = c(bar_col, plot_cols$cumulative),
    bty = "n",
    cex = 0.82,
    text.col = plot_cols$text
  )
}

png(OUT_PATH, width = 2600, height = 1450, res = 220, bg = "white")
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)

par(
  family = "sans",
  fg = plot_cols$text,
  mfrow = c(1, 2),
  mar = c(6.3, 5.1, 4.5, 1.2),
  oma = c(0.4, 0.2, 4.5, 0.2),
  xaxs = "i",
  yaxs = "i"
)

draw_panel(
  variance_data[variance_data$basis_id == "original_pca", ],
  plot_cols$original,
  plot_cols$original_light,
  "A. Original PCA"
)

draw_panel(
  variance_data[variance_data$basis_id == "vector_normalized_pca", ],
  plot_cols$vector,
  plot_cols$vector_light,
  "B. Vector-normalized PCA"
)

mtext(
  "PCA Axis Variance Explained",
  outer = TRUE,
  side = 3,
  line = 2.4,
  adj = 0.03,
  cex = 1.35,
  font = 2,
  col = plot_cols$text
)
mtext(
  "Bars show individual axis variance for PC1-PC10; black line shows cumulative variance.",
  outer = TRUE,
  side = 3,
  line = 1.0,
  adj = 0.03,
  cex = 0.88,
  col = plot_cols$subtext
)

cat("Created", OUT_PATH, "\n")
cat("Created", SOURCE_CSV, "\n")
cat("Original PCA PC1:", variance_data$pct_variance_percent[variance_data$basis_id == "original_pca" & variance_data$pc_axis == 1], "\n")
cat("Original PCA PC2:", variance_data$pct_variance_percent[variance_data$basis_id == "original_pca" & variance_data$pc_axis == 2], "\n")
cat("Vector-normalized PCA PC1:", variance_data$pct_variance_percent[variance_data$basis_id == "vector_normalized_pca" & variance_data$pc_axis == 1], "\n")
cat("Vector-normalized PCA PC2:", variance_data$pct_variance_percent[variance_data$basis_id == "vector_normalized_pca" & variance_data$pc_axis == 2], "\n")
