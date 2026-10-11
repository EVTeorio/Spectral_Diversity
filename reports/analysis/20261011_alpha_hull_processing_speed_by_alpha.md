# Alpha-Hull Processing Speed by Alpha

Date: 2026-10-11

## Purpose

Summarize processing time for the all-retained-pixel 2D R `alphahull` benchmark by alpha parameter and spatial scale.
Each timing uses the current matched benchmark subset: 16 selected 10 m quadrats, 8 selected 20 m quadrats, and 4 selected 50 m quadrats.
All alpha-hull values were calculated in vector-normalized PCA PC1-PC2 space after collapsing duplicate PC1-PC2 coordinates.

## Software and Operating System

| Item | Value |
|---|---|
| Operating system | Microsoft Windows 10.0.26200, x64 |
| R version | R version 4.2.3 (2023-03-15 ucrt) |
| R platform | x86_64-w64-mingw32 |
| `alphahull` version | 2.5 |
| `terra` version | 1.7.71 |
| Benchmark script | `scripts/3_Analysis/alpha_hull_all_pixels_benchmark.R` |
| Runtime table source | `reports/tables/alpha_hull_all_pixels_benchmark/alpha_hull_runtime_summary.csv` |

## Hardware Context

These hardware values were queried on 2026-10-11 from the same workstation used for the benchmark. Available memory and free disk space are point-in-time values.

| Item | Value |
|---|---|
| Computer model | HP ZBook Firefly 16 inch G10 Mobile Workstation PC |
| CPU | 13th Gen Intel(R) Core(TM) i7-1365U |
| CPU cores / logical processors | 10 cores / 12 logical processors |
| Reported maximum clock speed | 1,800 MHz |
| Installed physical memory | 31.64 GiB, reported as 33,974,120,448 bytes |
| RAM modules | 2 x 16 GiB Samsung modules |
| RAM speed | 5,600 MT/s module speed; 5,200 MT/s configured clock |
| Visible physical memory at query time | 31.64 GiB |
| Free physical memory at query time | 6.79 GiB |
| Total virtual memory at query time | 48.75 GiB |
| Free virtual memory at query time | 7.00 GiB |
| System drive | C: fixed NTFS drive, 920.37 GiB total, 25.01 GiB free at query time |
| Disk model | Samsung MZVL21T0HCLR-00BH1 |
| Workspace location | OneDrive-backed local path under `C:/Users/PaintRock/OneDrive - Alabama A&M University/` |

## Processing Time by Alpha and Spatial Scale

Times are total elapsed seconds for all selected quadrats within each spatial scale.

| alpha | 10 m seconds | 20 m seconds | 50 m seconds | total seconds |
|---:|---:|---:|---:|---:|
| 0.10 | 2,253.90 | 6,791.43 | 11,647.16 | 20,692.49 |
| 0.25 | 420.03 | 551.53 | 1,211.68 | 2,183.24 |
| 0.50 | 56.40 | 92.62 | 598.67 | 747.69 |
| 1.00 | 15.47 | 42.62 | 485.13 | 543.22 |
| 1.50 | 12.33 | 39.89 | 469.82 | 522.04 |
| 2.00 | 11.00 | 37.94 | 470.93 | 519.87 |

## Processing Time in Minutes and Hours

| alpha | 10 m | 20 m | 50 m | total |
|---:|---:|---:|---:|---:|
| 0.10 | 37.6 min | 1.89 hr | 3.24 hr | 5.75 hr |
| 0.25 | 7.0 min | 9.2 min | 20.2 min | 36.4 min |
| 0.50 | 56.4 sec | 1.5 min | 10.0 min | 12.5 min |
| 1.00 | 15.5 sec | 42.6 sec | 8.1 min | 9.1 min |
| 1.50 | 12.3 sec | 39.9 sec | 7.8 min | 8.7 min |
| 2.00 | 11.0 sec | 37.9 sec | 7.8 min | 8.7 min |

## Notes

- The runtimes are benchmark timings for this workstation and this selected quadrat subset, not a universal performance estimate.
- Smaller alpha values were much slower in R `alphahull`. Alpha = 0.10 took about 5.75 hours total for the 28-quadrat subset, compared with about 36.4 minutes for alpha = 0.25 and about 8.7-12.5 minutes for alpha values 0.5 through 2.
- The 50 m rows dominate runtime because they include many more retained and unique PC1-PC2 points than the smaller spatial scales.
- These timings were produced from all retained valid sunlit pixels for each included quadrat, with duplicate PC1-PC2 coordinates removed before alpha-hull calculation.
- Hardware and background system state can affect these timings. CPU frequency scaling, available RAM, disk free space, OneDrive synchronization, and other running processes may all change elapsed time for repeated runs.
