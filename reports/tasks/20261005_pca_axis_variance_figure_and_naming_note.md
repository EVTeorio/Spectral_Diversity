# Task Report: PCA Axis Variance Figure and Naming Note

## Objective

Create a figure comparing PCA axis variance explained for the original PCA and vector-normalized PCA bases, and establish "original PCA" as the forward-facing term for the former raw PCA.

## Actions Performed

- Created `scripts/visuals/create_pca_axis_variance_explained_figure.R`.
- Read the existing PCA variance-explained CSVs from `Quad_Values/Spectral_diversitySHPs/`.
- Generated a two-panel figure showing PC1-PC10 individual axis variance and cumulative variance.
- Saved a figure-source CSV with the plotted variance values.
- Updated project guidance and project state to note that future writing should call the former raw PCA the original PCA.

## Outputs

- `Documents/Tables and Figures/26_pca_axis_variance_explained_original_vs_vector_normalized.png`
- `reports/tables/figure_sources/26_pca_axis_variance_explained_original_vs_vector_normalized.csv`
- `scripts/visuals/create_pca_axis_variance_explained_figure.R`

## Naming Decision

Moving forward, manuscript text, figure titles, table labels, and interpretation notes should use **original PCA** instead of **raw PCA**. Existing filenames, column names, and code identifiers may remain unchanged where renaming would reduce reproducibility or break links to established outputs.
