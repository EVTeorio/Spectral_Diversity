# Illumination Training Polygon Area Summary

Date: 2026-09-14

## Purpose

Summarize the shadow and sunlit training polygons used by the illumination-screening workflow. The active illumination workflow reads `Shadow_vs_sunlit_SHP/Shadow_vs_sunlit_Rev.shp`.

## Active revised training polygons

Source layer: `Shadow_vs_sunlit_SHP/Shadow_vs_sunlit_Rev.shp`

CRS: WGS 84 geographic coordinates

| Category | Polygon count | Total area (m2) | Total area (ha) | Mean polygon area (m2) |
|---|---:|---:|---:|---:|
| shadow | 30 | 229.104 | 0.022910 | 7.637 |
| sunlit | 50 | 198.879 | 0.019888 | 3.978 |
| total | 80 | 427.982 | 0.042798 | 5.350 |

## Previous training polygons

The earlier, non-revised layer is retained in the same folder. It is not the layer used by the active `Sunlit_vs_Shadow_AUC.R` workflow, but it is summarized here to document the difference.

Source layer: `Shadow_vs_sunlit_SHP/Shadow_vs_sunlit.shp`

| Category | Polygon count | Total area (m2) | Total area (ha) | Mean polygon area (m2) |
|---|---:|---:|---:|---:|
| shadow | 27 | 200.596 | 0.020060 | 7.429 |
| sunlit | 33 | 285.037 | 0.028504 | 8.637 |
| total | 60 | 485.633 | 0.048563 | 8.094 |

## Calculation note

The shapefiles are stored in WGS 84 latitude/longitude, so raw coordinate area would be in degree-squared units and is not meaningful for reporting. Areas above were calculated by reading the shapefile polygons and projecting coordinates to a local meter-based tangent plane centered on the polygon set latitude before applying the polygon shoelace area formula. Because these are small training polygons, this local-area approximation is appropriate for a protocol summary.

The active revised layer increased the number of training polygons from 60 to 80. Shadow polygons increased from 27 to 30 and sunlit polygons increased from 33 to 50. Total sampled polygon area decreased from 485.633 m2 to 427.982 m2, mostly because the revised sunlit polygons are smaller on average.
