PROJECT_DIR <- "C:/Users/PaintRock/OneDrive - Alabama A&M University/PaintRock RemoteSens/Spectral_Diversity"
TABLE_DIR <- file.path(PROJECT_DIR, "reports/tables/pca_sa_phylogenetic_only_multiple_regression")
PAPER_TABLE <- file.path(PROJECT_DIR, "Documents/Paper/20261009_phylogenetic_only_pca_sa_multiple_regression_table.md")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

required_files <- c(
  "pca_sa_phylogenetic_only_full_model_summary.csv",
  "pca_sa_phylogenetic_only_full_model_coefficients.csv",
  "pca_sa_phylogenetic_only_partial_r2.csv",
  "pca_sa_phylogenetic_only_predictor_vif.csv",
  "pca_sa_phylogenetic_only_top_full_models.csv",
  "pca_sa_phylogenetic_only_validation_checks.csv"
)

missing_files <- required_files[!file.exists(file.path(TABLE_DIR, required_files))]
assert_true(length(missing_files) == 0, paste("Missing output files:", paste(missing_files, collapse = ", ")))
assert_true(file.exists(PAPER_TABLE), "Missing paper-facing Markdown table.")

model_summary <- read.csv(file.path(TABLE_DIR, "pca_sa_phylogenetic_only_full_model_summary.csv"), stringsAsFactors = FALSE)
coefs <- read.csv(file.path(TABLE_DIR, "pca_sa_phylogenetic_only_full_model_coefficients.csv"), stringsAsFactors = FALSE)
partial_r2 <- read.csv(file.path(TABLE_DIR, "pca_sa_phylogenetic_only_partial_r2.csv"), stringsAsFactors = FALSE)
vifs <- read.csv(file.path(TABLE_DIR, "pca_sa_phylogenetic_only_predictor_vif.csv"), stringsAsFactors = FALSE)
checks <- read.csv(file.path(TABLE_DIR, "pca_sa_phylogenetic_only_validation_checks.csv"), stringsAsFactors = FALSE)

assert_true(nrow(model_summary) == 21, "Expected 21 full model rows.")
assert_true("spec_sa" %in% model_summary$spectral_metric, "Expected spectral angle entropy model rows.")
assert_true(nrow(coefs) == 84, "Expected 84 coefficient rows.")
assert_true(nrow(partial_r2) == 63, "Expected 63 partial R2 rows.")
assert_true(nrow(vifs) == 63, "Expected 63 VIF rows.")
assert_true(all(checks$passed), "Expected all workflow validation checks to pass.")
assert_true(all(is.finite(model_summary$r_squared)), "Expected finite model R2 values.")
assert_true(all(model_summary$r_squared >= 0 & model_summary$r_squared <= 1), "Expected model R2 values between 0 and 1.")
assert_true(all(vifs$vif >= 1), "Expected all VIF values to be at least 1.")

message("PCA and spectral angle phylogenetic-only multiple-regression output tests passed.")
