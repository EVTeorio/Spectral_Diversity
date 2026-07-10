# Validation: 50 m PC1-PC3 Mean Reflectance Correlation

Date: 2026-07-07

## Checks

| check | expected | observed |
| --- | --- | --- |
| 50m raster files found | 80.0000 | 80 |
| pixels per quadrat | 250.0000 | 250 |
| total sampled pixels per PCA basis | 20000.0000 | 20000, 20000 |
| quadrat mean rows per PCA basis | 80.0000 | 80, 80 |
| manual excluded quadrats included | 6.0000 | 6 |

## Result

Pass. The script sampled 250 valid sunlit pixels from each of 80 current 50 m quadrat rasters and projected the same sampled spectra into both regular and standardized PCA space.
