# Show the failed case's PCA points; do not invent a hull boundary.
.libPaths(c('C:/Users/PaintRock/AppData/Local/R/win-library/4.2', .libPaths()))
expressions <- parse('scripts/2_Indices Creation/Spectral_diversity/spectral_heterogeneity_all_metrics.R')
stopifnot(grepl('RUN_SPECTRAL_HET_WORKFLOW',
  paste(deparse(expressions[[length(expressions)]]), collapse = '')))
helper <- new.env()
eval(expressions[-length(expressions)], helper)
raster_file <- 'Quad_Spectra/10m_smooth_5nm/1000_c'
model <- readRDS(helper$PCA_RDS)
stopifnot(identical(names(terra::rast(raster_file)), names(model$center)))
x <- helper$read_masked_spectra(raster_file)
scores <- helper$project_scores(x, model)[, 1:2, drop = FALSE]
stopifnot(nrow(scores) == 7126, nrow(unique(scores)) == 6890, all(is.finite(scores)))
out <- 'reports/figures/20261015_original_pca_10m_1000_c.png'
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
png(out, width = 1800, height = 1350, res = 180)
par(mar = c(5, 5, 4, 1.5), oma = c(4, 0, 1, 0), las = 1, family = 'sans')
plot(scores, asp = 1, pch = 16, cex = .65, col = '#344F6555',
  xlab = 'PC1 score', ylab = 'PC2 score',
  main = '10 m quadrat 1000_c | Original PCA', cex.main = 1.2)
legend('topright', c('7,126 retained pixels', '6,890 unique PC1-PC2 points'),
  bty = 'n', cex = .9)
mtext('Point cloud only: alpha = 0.5 hull construction failed; no boundary was returned.',
  side = 1, outer = TRUE, line = 1, cex = .9)
mtext('The location of the construction failure has not been identified in PCA space.',
  side = 1, outer = TRUE, line = 2.2, cex = .85)
dev.off()
cat('Created ', out, '\n', sep = '')
