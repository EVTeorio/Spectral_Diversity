# Alpha-Hull Average Processing Speed per Quadrat

Date: 2026-10-11

## Purpose

Summarize average processing time per quadrat for the all-retained-pixel 2D R `alphahull` benchmark by alpha parameter and spatial scale.
This report uses the same runtime source as the total-time report but divides each scale total by the number of selected quadrats at that scale: 16 selected 10 m quadrats, 8 selected 20 m quadrats, and 4 selected 50 m quadrats.

## Software, Operating System, and Hardware

| Item | Value |
|---|---|
| Operating system | Microsoft Windows 10.0.26200, x64 |
| R version | R version 4.2.3 (2023-03-15 ucrt) |
| R platform | x86_64-w64-mingw32 |
| `alphahull` version | 2.5 |
| `terra` version | 1.7.71 |
| Computer model | HP ZBook Firefly 16 inch G10 Mobile Workstation PC |
| CPU | 13th Gen Intel(R) Core(TM) i7-1365U |
| CPU cores / logical processors | 10 cores / 12 logical processors |
| Installed physical memory | 31.64 GiB |
| Benchmark script | `scripts/3_Analysis/alpha_hull_all_pixels_benchmark.R` |
| Runtime table source | `reports/tables/alpha_hull_all_pixels_benchmark/alpha_hull_runtime_summary.csv` |

## Average Processing Time per Quadrat

Times are average elapsed seconds per selected quadrat within each spatial scale.

| alpha | 10 m seconds / quadrat | 20 m seconds / quadrat | 50 m seconds / quadrat |
|---:|---:|---:|---:|
| 0.10 | 140.87 | 848.93 | 2,911.79 |
| 0.25 | 26.25 | 68.94 | 302.92 |
| 0.50 | 3.53 | 11.58 | 149.67 |
| 1.00 | 0.97 | 5.33 | 121.28 |
| 1.50 | 0.77 | 4.99 | 117.46 |
| 2.00 | 0.69 | 4.74 | 117.73 |

## Average Processing Time in Minutes

| alpha | 10 m | 20 m | 50 m |
|---:|---:|---:|---:|
| 0.10 | 2.35 min/quadrat | 14.15 min/quadrat | 48.53 min/quadrat |
| 0.25 | 26.25 sec/quadrat | 1.15 min/quadrat | 5.05 min/quadrat |
| 0.50 | 3.53 sec/quadrat | 11.58 sec/quadrat | 2.49 min/quadrat |
| 1.00 | 0.97 sec/quadrat | 5.33 sec/quadrat | 2.02 min/quadrat |
| 1.50 | 0.77 sec/quadrat | 4.99 sec/quadrat | 1.96 min/quadrat |
| 2.00 | 0.69 sec/quadrat | 4.74 sec/quadrat | 1.96 min/quadrat |

## Notes

- These are scale-level averages, not per-quadrat medians. Individual quadrats can differ substantially depending on retained-pixel count and unique PC1-PC2 point count.
- Values were calculated from total elapsed time by scale divided by the selected quadrat count for that scale.
- Alpha = 0.10 was the slowest setting at every spatial scale, averaging about 48.5 minutes per 50 m quadrat.
- The original total-runtime report remains unchanged at `reports/analysis/20261011_alpha_hull_processing_speed_by_alpha.md`.
