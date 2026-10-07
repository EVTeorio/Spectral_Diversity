# Validation: PCA Axis Variance Figure

## Validation Steps

- Ran `scripts/visuals/create_pca_axis_variance_explained_figure.R` with R 4.2.3.
- Confirmed that the figure-source CSV was created from the two existing PCA variance-explained CSVs.
- Visually inspected the exported PNG.
- Checked the plotted PC1 and PC2 percentages against source values.

## Results

| PCA basis | PC1 variance explained | PC2 variance explained | Cumulative PC1-PC2 |
|---|---:|---:|---:|
| Original PCA | 66.23% | 20.72% | 86.95% |
| Vector-normalized PCA | 45.46% | 22.54% | 68.00% |

## Conclusion

The figure correctly displays PC1-PC10 individual and cumulative variance explained for both PCA bases. The naming convention has been recorded so future forward-facing work uses "original PCA" for the former raw PCA.
