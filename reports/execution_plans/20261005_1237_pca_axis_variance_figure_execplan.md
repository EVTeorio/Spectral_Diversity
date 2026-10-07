# ExecPlan: PCA Axis Variance Figure and Naming Note

## Objective

Create a figure comparing the variance explained by PCA axes for the original PCA and vector-normalized PCA bases, and document that former references to "raw PCA" should be called "original PCA" in forward-facing work.

## Requested Task

- Create a figure showing axis variance explained for both PCA bases.
- Add a project-wide note that moving forward, the raw PCA will be called the original PCA.

## Files To Review

- `CODEX_AGENT_GUIDELINES.md`
- `reports/project_state.md`
- `Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm_variance_explained.csv`
- `Quad_Values/Spectral_diversitySHPs/standardized_PCA_global_pca_smooth_masked_5nm_variance_explained.csv`

## Relevant Files

- `scripts/visuals/create_pca_axis_variance_explained_figure.R`
- `Documents/Tables and Figures/26_pca_axis_variance_explained_original_vs_vector_normalized.png`
- `reports/tables/figure_sources/26_pca_axis_variance_explained_original_vs_vector_normalized.csv`

## Proposed Changes

- Add an R plotting script using the existing PCA variance-explained CSVs.
- Save a publication-support figure under `Documents/Tables and Figures/`.
- Save the plotted data under `reports/tables/figure_sources/`.
- Update governance/project-state notes to establish "original PCA" as the forward-facing term for the former raw PCA.

## Validation Plan

- Confirm the plotted PC1 and PC2 percentages match the source CSVs:
  - original PCA: PC1 = 66.23%, PC2 = 20.72%
  - vector-normalized PCA: PC1 = 45.46%, PC2 = 22.54%
- Confirm output PNG and figure-source CSV exist.

## Risks

- Existing script and column names still use raw/standardized naming for reproducibility. The note should clarify that the naming change is interpretive/manuscript-facing and does not rename existing files or columns.
