PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
setwd(PROJECT_DIR)

PCA_RDS <- file.path(PROJECT_DIR, "Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm.rds")
SPEC_DIR <- file.path(PROJECT_DIR, "Quad_Spectra/10m_smooth_5nm")
OUT_DIR <- file.path(PROJECT_DIR, "reports/tables/figure_sources/10m_alpha_examples")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

samples <- 136L
lines <- 136L
bands <- 121L
n_pixels <- samples * lines
shadow_threshold <- 0.0305476
mask_band_index <- 34L

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

pca_object <- readRDS(PCA_RDS)
candidates <- list.files(SPEC_DIR, pattern = "^[0-9]+_a$", full.names = FALSE)
candidates <- unique(c("700_a", candidates[seq_len(min(length(candidates), 80L))]))

summary_rows <- list()
for (tile in candidates) {
  scores <- tryCatch(clean_project_scores(tile, pca_object), error = function(e) NULL)
  if (is.null(scores) || nrow(scores) < 500) next
  set.seed(2500 + sum(utf8ToInt(tile)))
  plot_scores <- scores[sample.int(nrow(scores), min(900L, nrow(scores))), , drop = FALSE]
  out_path <- file.path(OUT_DIR, paste0(tile, ".csv"))
  write.csv(data.frame(PC1 = plot_scores[, 1], PC2 = plot_scores[, 2]), out_path, row.names = FALSE)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    tile = tile,
    retained_pixels = nrow(scores),
    plotted_pixels = nrow(plot_scores),
    csv = out_path
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR, "candidate_summary.csv"), row.names = FALSE)
cat("Exported", nrow(summary_df), "candidate 10 m score tables to", OUT_DIR, "\n")
