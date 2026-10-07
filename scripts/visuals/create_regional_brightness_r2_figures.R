input_csv <- "reports/tables/spectral_heterogeneity_relationships/spectral_metric_regional_illumination_correlations.csv"
out_dir <- "Documents/Tables and Figures"

d <- read.csv(input_csv, stringsAsFactors = FALSE)

region_levels <- c("blue", "green", "red", "nir")
region_labels <- c("Blue\n450-495 nm", "Green\n500-570 nm", "Red\n620-750 nm", "Near infrared\n750-998 nm")
scale_levels <- c("10m", "20m", "50m")
metric_colors <- c("#1b9e77", "#d95f02", "#7570b3")

make_plot <- function(metrics, labels, outfile, title, subtitle) {
  dd <- d[d$metric %in% metrics, ]
  dd$metric_label2 <- labels[match(dd$metric, metrics)]
  dd$region_i <- match(dd$region, region_levels)
  dd$scale <- factor(dd$scale, levels = scale_levels)

  png(file.path(out_dir, outfile), width = 2600, height = 1100, res = 220)
  par(mfrow = c(1, 3), mar = c(5.2, 4.8, 4.7, 1.2), oma = c(0, 0, 3.8, 0), family = "sans")

  y_max <- max(dd$r_squared, na.rm = TRUE) * 1.08
  for (sc in scale_levels) {
    plot(
      NA,
      xlim = c(0.75, 4.25),
      ylim = c(0, y_max),
      xaxt = "n",
      xlab = "Brightness region ordered across the electromagnetic spectrum",
      ylab = expression(R^2),
      main = paste0(gsub("m", " m", sc), " quadrats"),
      cex.main = 1.05
    )
    axis(1, at = seq_along(region_labels), labels = region_labels)
    grid(nx = NA, ny = NULL, col = "grey88", lty = 1)
    abline(v = seq_along(region_labels), col = "grey92", lty = 1)

    for (i in seq_along(metrics)) {
      sub <- dd[dd$scale == sc & dd$metric == metrics[i], ]
      sub <- sub[order(sub$region_i), ]
      lines(
        sub$region_i,
        sub$r_squared,
        type = "b",
        lwd = 2.3,
        pch = 16 + i,
        col = metric_colors[i],
        cex = 1.05
      )
    }

    if (sc == scale_levels[1]) {
      legend(
        "topright",
        legend = labels,
        col = metric_colors,
        lwd = 2.3,
        pch = 17:(16 + length(metrics)),
        bty = "n",
        cex = 0.78
      )
    }
  }

  mtext(title, outer = TRUE, side = 3, line = 2.1, font = 2, cex = 1.2)
  mtext(subtitle, outer = TRUE, side = 3, line = 0.7, cex = 0.82)
  dev.off()
}

make_plot(
  metrics = c("spec_alpha", "spec_pca_mean", "spec_rao_q"),
  labels = c("Alpha-hull area", "Mean distance", "Spectral Rao's Q"),
  outfile = "23_raw_pca_regional_brightness_r2_by_spectrum.png",
  title = "Raw PCA Metrics: Regional Brightness Driver Strength",
  subtitle = "R-squared values from Pearson relationships between raw PCA spectral heterogeneity metrics and retained-pixel regional brightness"
)

make_plot(
  metrics = c("spec_spca_alpha", "spec_spca_mean", "spec_spca_rao"),
  labels = c("Alpha-hull area", "Mean distance", "Spectral Rao's Q"),
  outfile = "24_vector_normalized_pca_regional_brightness_r2_by_spectrum.png",
  title = "Vector-Normalized PCA Metrics: Regional Brightness Driver Strength",
  subtitle = "R-squared values from Pearson relationships between vector-normalized PCA spectral heterogeneity metrics and retained-pixel regional brightness"
)
