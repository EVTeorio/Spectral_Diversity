# Bootstrap Variation Analysis for Spectral Heterogeneity

Date: 2026-06-18

Updated: 2026-06-22

## Purpose

This report evaluates whether the bootstrap variation within quadrats is small enough relative to variation among quadrats to treat the spectral heterogeneity values as valid for downstream biodiversity analyses.

The key distinction is that `boot_sd` measures how much spectral entropy changes across repeated pixel subsamples within the same quadrat, while `boot_sd / sqrt(70)` estimates the uncertainty of the reported `boot_mean` value.

## Inputs

- `Indices_SHPs/10m_SA_entropy_smooth_masked_5nm_summary.csv` and matching bootstrap files
- `Indices_SHPs/20m_SA_entropy_smooth_masked_5nm_summary.csv` and matching bootstrap files
- `Indices_SHPs/50m_SA_entropy_smooth_masked_5nm_summary.csv` and matching bootstrap files

All inputs were generated from current smoothed 5 nm spectra using the shadow mask threshold of `0.0305476` at `563 nm`.

## Main Conclusion

The bootstrap means are usable with quality-control flags, but raw within-quadrat bootstrap variation is not uniformly small enough to ignore. A subset of quadrats, especially at 10 m and 20 m, should be flagged or included in sensitivity checks before downstream modeling.

Decision: use the reported `spectral_entropy` / `boot_mean` values as the main heterogeneity estimates, but do not treat all quadrats as equally stable. Carry `boot_sd`, `boot_cv`, and `method` forward as quality-control fields, and run sensitivity checks that flag or exclude high-CV quadrats.

Across all bootstrap-estimated quadrats, the 95th percentile bootstrap CV was 10.4%, and the maximum bootstrap CV was 13.5% in quadrat `420_a` at 10 m.

The most important validation metric is the ratio of mean within-quadrat bootstrap variance to between-quadrat variance. Lower values mean the bootstrap noise is small relative to the spatial signal among quadrats.

Using a 5% bootstrap-CV flag, 372 of 1866 bootstrap-estimated 10 m quadrats, 60 of 477 20 m quadrats, and 11 of 79 50 m quadrats would be flagged.

## Scale Summary



| Scale| Valid estimates| Bootstrap quads| Exact quads| Mean entropy| Median entropy| Skewness| Between SD| Median within SD| Mean within SD| Mean boot mean SE| Within/Between variance| Between fraction| Median CV| P90 CV|  Class|
|-----:|---------------:|---------------:|-----------:|------------:|--------------:|--------:|----------:|----------------:|--------------:|-----------------:|-----------------------:|----------------:|---------:|------:|------:|
|  10 m|            1903|            1866|          37|        2.843|          2.829|    0.275|      0.108|            0.010|          0.068|             0.008|                  121.6%|            45.1%|      0.4%|   8.5%| Review|
|  20 m|             485|             477|           8|        2.835|          2.821|   -0.002|      0.098|            0.013|          0.053|             0.006|                   92.3%|            52.0%|      0.5%|   5.7%| Review|
|  50 m|              80|              79|           1|        2.847|          2.833|    1.187|      0.076|            0.014|          0.053|             0.006|                  122.4%|            45.0%|      0.5%|   5.8%| Review|

## Figures

### Within Versus Between SD

![Within versus between SD](../figures/bootstrap_variation/within_vs_between_sd.png)

### Bootstrap CV Distribution

Dashed and dotted reference lines mark 5% and 10% bootstrap CV.

![Bootstrap CV distribution](../figures/bootstrap_variation/bootstrap_cv_distribution.png)

### Bootstrap SD Versus Entropy

![Bootstrap SD versus entropy](../figures/bootstrap_variation/bootstrap_sd_vs_entropy.png)

### Spectral Entropy Histograms

These histograms show the distribution and skewness of the final spectral heterogeneity values at each scale. The dashed red line is the mean and the solid blue line is the median.

![Spectral entropy histograms by scale](../figures/bootstrap_variation/spectral_entropy_histograms_by_scale.png)

### Pixel Count Versus Bootstrap CV

![Pixel count versus bootstrap CV](../figures/bootstrap_variation/pixel_count_vs_bootstrap_cv.png)

### Variance Partition Diagnostic

![Variance partition diagnostic](../figures/bootstrap_variation/variance_partition_diagnostic.png)

### Most Variable Quadrat Replicate Distributions

![Top unstable quadrat bootstrap distributions](../figures/bootstrap_variation/top_unstable_quadrat_bootstrap_distributions.png)

## Most Variable Quadrats



| Scale|  Quadrat| Pixels| Entropy| Boot SD| Boot CV| Outlier reps|
|-----:|--------:|------:|-------:|-------:|-------:|------------:|
|  10 m|    420_a|   7179|   2.360|   0.318|   13.5%|            0|
|  10 m|    316_c|   7837|   2.559|   0.343|   13.4%|            0|
|  10 m|    511_c|   7120|   2.647|   0.343|   12.9%|            0|
|  10 m|   1414_a|   6227|   2.707|   0.345|   12.7%|            0|
|  10 m|   1711_d|   8185|   2.628|   0.334|   12.7%|            0|
|  10 m|   1208_b|   5509|   2.641|   0.333|   12.6%|            0|
|  10 m|   1024_c|   6205|   2.635|   0.331|   12.6%|            0|
|  10 m|   1316_d|   6691|   2.735|   0.342|   12.5%|            0|
|  20 m|      918|  19833|   2.761|   0.346|   12.5%|            0|
|  20 m|      804|  32325|   2.792|   0.346|   12.4%|            0|
|  20 m|      412|  27508|   2.667|   0.319|   12.0%|            0|
|  20 m|     1217|   9606|   2.831|   0.338|   12.0%|            0|
|  20 m|       24|  15470|   2.699|   0.321|   11.9%|            0|
|  20 m|      419|  25439|   2.551|   0.299|   11.7%|            0|
|  20 m|      719|  29838|   2.877|   0.328|   11.4%|            0|
|  20 m|      713|  28554|   2.780|   0.309|   11.1%|            0|
|  50 m| sub50_12| 169495|   2.811|   0.281|   10.0%|           15|
|  50 m| sub50_54| 156070|   2.854|   0.218|    7.6%|            8|
|  50 m|  sub50_4| 144375|   2.950|   0.221|    7.5%|            8|
|  50 m| sub50_58| 189724|   2.774|   0.177|    6.4%|            5|
|  50 m| sub50_44| 152090|   3.011|   0.191|    6.3%|            6|
|  50 m| sub50_60| 164492|   2.887|   0.178|    6.2%|            5|
|  50 m| sub50_62|  48957|   3.148|   0.192|    6.1%|            6|
|  50 m| sub50_53|  73073|   3.001|   0.179|    6.0%|            6|

## Interpretation

- Median bootstrap CV is the best quick read on typical within-quadrat stability.
- The mean bootstrap SD can be pulled upward by a small number of quadrats with occasional low or high bootstrap replicates.
- The SE of the bootstrap mean is much smaller than the raw bootstrap SD because each reported value averages 70 bootstrap replicates.
- The raw within-quadrat bootstrap variance is large enough in some quadrats that it should be represented in downstream model diagnostics.
- Exact all-pixel estimates were used where computationally reasonable; most quadrats use `bootstrap_mean` because all-pixel pair counts were too large.
- Some bootstrap variation may reflect the sampled-pair approximation used for computationally large bootstrap replicates; high-CV quadrats are the best candidates for reruns with larger pair samples if more precision is needed.

## Recommendation

Use `spectral_entropy` from the per-scale summary CSVs as the primary spectral heterogeneity value, but carry `boot_sd`, `boot_cv`, and `method` into downstream analysis tables as quality-control fields. In model sensitivity checks, consider excluding or flagging quadrats with bootstrap CV above 5% or with many outlying bootstrap replicates.

## Output Tables

- `reports/tables/bootstrap_variation/bootstrap_variation_scale_summary.csv`
- `reports/tables/bootstrap_variation/bootstrap_variation_quadrat_diagnostics.csv`
- `reports/tables/bootstrap_variation/bootstrap_variation_top_unstable_quadrats.csv`
