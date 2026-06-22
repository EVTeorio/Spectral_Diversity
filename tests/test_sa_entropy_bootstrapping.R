Sys.setenv(RUN_SA_ENTROPY_WORKFLOW = "false")
source("scripts/2_Indices Creation/Spectral_diversity/SA_entropy_bootstrapping.R")

assert_equal <- function(actual, expected, tolerance = sqrt(.Machine$double.eps)) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop("Expected ", paste(expected, collapse = ", "), " but got ", paste(actual, collapse = ", "))
  }
}

assert_true <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message)
  }
}

assert_equal(get_quad_id("Quad_Spectra/10m_smooth_5nm/0_a", "10m"), "0_a")
assert_equal(get_quad_id("Quad_Spectra/20m_smooth_5nm/1000", "20m"), "1000")
assert_equal(get_quad_id("Quad_Spectra/50m_smooth_5nm/sub50_10", "50m"), "sub50_10")

spectra <- rbind(
  c(3, 4, 0),
  c(0, 0, 0),
  c(NA, 1, 1),
  c(1, 2, 2)
)
normalized <- normalize_spectra(spectra)

assert_equal(nrow(normalized), 2)
assert_equal(as.numeric(sqrt(rowSums(normalized^2))), c(1, 1), tolerance = 1e-12)

spectra_norm <- rbind(
  c(1, 0),
  c(0, 1),
  c(1, 0)
)
angles <- sort(pairwise_angles(spectra_norm, block_size = 2))

assert_equal(length(angles), 3)
assert_equal(angles, sort(c(pi / 2, 0, pi / 2)), tolerance = 1e-12)

spectra_norm <- normalize_spectra(rbind(
  c(1, 0, 0),
  c(0, 1, 0),
  c(0, 0, 1),
  c(1, 1, 0)
))
entropy <- spectral_entropy_from_norm(spectra_norm, n_bins = 4)

assert_true(is.finite(entropy), "Expected finite entropy")
assert_true(entropy > 0, "Expected entropy greater than zero")

set.seed(42)
sampled_angles <- sample_pairwise_angles(spectra_norm, n_pairs = 10)
assert_equal(length(sampled_angles), 10)
assert_true(all(is.finite(sampled_angles)), "Expected finite sampled pair angles")

sampled_entropy <- spectral_entropy_from_sampled_pairs(spectra_norm, n_pairs = 10, n_bins = 4)
assert_true(is.finite(sampled_entropy), "Expected finite sampled-pair entropy")

cat("SA entropy helper tests passed\n")
