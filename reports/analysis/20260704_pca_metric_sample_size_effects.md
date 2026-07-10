# Regular PCA-Derived Spectral Metric Sample-Size Effects

Date: 2026-07-04

## Purpose

Extend the sample-size sensitivity analysis from spectral angle entropy to regular PCA mean distance, spectral Rao's Q, and alpha-hull area.

## Design

- Bootstrap iterations per quadrat x sample rule: 50
- PCA basis: Regular PCA (`Quad_Values/Spectral_diversitySHPs/global_pca_smooth_masked_5nm.rds`).
- Quadrat and sample-size rules are reused from `reports/tables/sample_size_effects/sa_entropy/sa_entropy_sample_size_design.csv` so the non-SA metrics are evaluated on the same quadrats and retained-pixel sample sizes.
- Each replicate samples retained pixels without replacement. If a rule resolves to 100% of retained pixels, the metric is calculated once from the full retained-pixel set and repeated across the 50 output rows, so full-pixel conditions have zero artificial bootstrap variation.
- PCA mean distance and spectral Rao's Q use all sampled pixels in PCA space. Alpha-hull area uses all sampled PC1-PC2 points unless the existing alpha-hull helper has to remove duplicate points or returns an internal failure method.

## Selected Quadrats



|Scale |Quadrat  |Retained pixels |Selection note                                                                              |
|:-----|:--------|:---------------|:-------------------------------------------------------------------------------------------|
|10m   |1104_b   |9,423           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1110_a   |7,650           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1124_a   |8,301           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |112_b    |7,340           |Original random bootstrap-mean quadrat selected with seed 20260703                          |
|10m   |114_b    |7,239           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |11_c     |6,199           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |120_d    |7,068           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1305_c   |8,027           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1315_b   |9,077           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1400_b   |6,806           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1414_a   |6,227           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1417_c   |7,513           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1514_a   |7,027           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1516_c   |8,358           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1604_c   |7,229           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1606_c   |7,594           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1701_b   |9,299           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1803_c   |5,295           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1816_c   |5,193           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1912_c   |5,467           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |1917_c   |5,452           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |204_d    |7,147           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |217_d    |7,850           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |22_c     |7,136           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |319_c    |8,983           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |409_d    |5,989           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |503_c    |7,403           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |614_a    |7,102           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |700_c    |7,623           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |722_a    |8,781           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|10m   |800_a    |3,976           |Original 10 m replacement: closest bootstrap-mean quadrat at or below 4,000 retained pixels |
|10m   |914_a    |7,220           |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |100      |30,808          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |1201     |24,803          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |1402     |27,061          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |1512     |30,984          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |1518     |35,713          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |1800     |29,990          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |1805     |25,420          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |1910     |25,613          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |206      |23,348          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |219      |29,380          |Original random bootstrap-mean quadrat selected with seed 20260703                          |
|20m   |302      |26,326          |Original random bootstrap-mean quadrat selected with seed 20260703                          |
|20m   |317      |31,152          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |405      |36,367          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |821      |29,194          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |905      |27,021          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|20m   |912      |29,618          |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|50m   |sub50_15 |172,782         |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|50m   |sub50_24 |166,816         |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|50m   |sub50_36 |174,203         |Original random bootstrap-mean quadrat selected with seed 20260703                          |
|50m   |sub50_49 |206,958         |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|50m   |sub50_51 |200,631         |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|50m   |sub50_56 |145,500         |Additional random bootstrap-mean quadrat selected with seed 20260703                        |
|50m   |sub50_60 |164,492         |Original random bootstrap-mean quadrat selected with seed 20260703                          |
|50m   |sub50_8  |187,893         |Additional random bootstrap-mean quadrat selected with seed 20260703                        |

## Sample Rules



|Scale |Quadrat  |Retained pixels |Sample rule                        |Sample pixels |Sample percent |
|:-----|:--------|:---------------|:----------------------------------|:-------------|:--------------|
|10m   |1104_b   |9,423           |1% of pixels = 94 pixels (1.0%)    |94            |1.0%           |
|10m   |1104_b   |9,423           |2% of pixels = 188 pixels (2.0%)   |188           |2.0%           |
|10m   |1104_b   |9,423           |3% of pixels = 283 pixels (3.0%)   |283           |3.0%           |
|10m   |1104_b   |9,423           |Fixed 1,250 (13.3%)                |1,250         |13.3%          |
|10m   |1104_b   |9,423           |Fixed 2,000 (21.2%)                |2,000         |21.2%          |
|10m   |1104_b   |9,423           |Fixed 3,000 (31.8%)                |3,000         |31.8%          |
|10m   |1104_b   |9,423           |Fixed 4,000 (42.4%)                |4,000         |42.4%          |
|10m   |1110_a   |7,650           |1% of pixels = 76 pixels (1.0%)    |76            |1.0%           |
|10m   |1110_a   |7,650           |2% of pixels = 153 pixels (2.0%)   |153           |2.0%           |
|10m   |1110_a   |7,650           |3% of pixels = 230 pixels (3.0%)   |230           |3.0%           |
|10m   |1110_a   |7,650           |Fixed 1,250 (16.3%)                |1,250         |16.3%          |
|10m   |1110_a   |7,650           |Fixed 2,000 (26.1%)                |2,000         |26.1%          |
|10m   |1110_a   |7,650           |Fixed 3,000 (39.2%)                |3,000         |39.2%          |
|10m   |1110_a   |7,650           |Fixed 4,000 (52.3%)                |4,000         |52.3%          |
|10m   |1124_a   |8,301           |1% of pixels = 83 pixels (1.0%)    |83            |1.0%           |
|10m   |1124_a   |8,301           |2% of pixels = 166 pixels (2.0%)   |166           |2.0%           |
|10m   |1124_a   |8,301           |3% of pixels = 249 pixels (3.0%)   |249           |3.0%           |
|10m   |1124_a   |8,301           |Fixed 1,250 (15.1%)                |1,250         |15.1%          |
|10m   |1124_a   |8,301           |Fixed 2,000 (24.1%)                |2,000         |24.1%          |
|10m   |1124_a   |8,301           |Fixed 3,000 (36.1%)                |3,000         |36.1%          |
|10m   |1124_a   |8,301           |Fixed 4,000 (48.2%)                |4,000         |48.2%          |
|10m   |112_b    |7,340           |1% of pixels = 73 pixels (1.0%)    |73            |1.0%           |
|10m   |112_b    |7,340           |2% of pixels = 147 pixels (2.0%)   |147           |2.0%           |
|10m   |112_b    |7,340           |3% of pixels = 220 pixels (3.0%)   |220           |3.0%           |
|10m   |112_b    |7,340           |Fixed 1,250 (17.0%)                |1,250         |17.0%          |
|10m   |112_b    |7,340           |Fixed 2,000 (27.2%)                |2,000         |27.2%          |
|10m   |112_b    |7,340           |Fixed 3,000 (40.9%)                |3,000         |40.9%          |
|10m   |112_b    |7,340           |Fixed 4,000 (54.5%)                |4,000         |54.5%          |
|10m   |114_b    |7,239           |1% of pixels = 72 pixels (1.0%)    |72            |1.0%           |
|10m   |114_b    |7,239           |2% of pixels = 145 pixels (2.0%)   |145           |2.0%           |
|10m   |114_b    |7,239           |3% of pixels = 217 pixels (3.0%)   |217           |3.0%           |
|10m   |114_b    |7,239           |Fixed 1,250 (17.3%)                |1,250         |17.3%          |
|10m   |114_b    |7,239           |Fixed 2,000 (27.6%)                |2,000         |27.6%          |
|10m   |114_b    |7,239           |Fixed 3,000 (41.4%)                |3,000         |41.4%          |
|10m   |114_b    |7,239           |Fixed 4,000 (55.3%)                |4,000         |55.3%          |
|10m   |11_c     |6,199           |1% of pixels = 62 pixels (1.0%)    |62            |1.0%           |
|10m   |11_c     |6,199           |2% of pixels = 124 pixels (2.0%)   |124           |2.0%           |
|10m   |11_c     |6,199           |3% of pixels = 186 pixels (3.0%)   |186           |3.0%           |
|10m   |11_c     |6,199           |Fixed 1,250 (20.2%)                |1,250         |20.2%          |
|10m   |11_c     |6,199           |Fixed 2,000 (32.3%)                |2,000         |32.3%          |
|10m   |11_c     |6,199           |Fixed 3,000 (48.4%)                |3,000         |48.4%          |
|10m   |11_c     |6,199           |Fixed 4,000 (64.5%)                |4,000         |64.5%          |
|10m   |120_d    |7,068           |1% of pixels = 71 pixels (1.0%)    |71            |1.0%           |
|10m   |120_d    |7,068           |2% of pixels = 141 pixels (2.0%)   |141           |2.0%           |
|10m   |120_d    |7,068           |3% of pixels = 212 pixels (3.0%)   |212           |3.0%           |
|10m   |120_d    |7,068           |Fixed 1,250 (17.7%)                |1,250         |17.7%          |
|10m   |120_d    |7,068           |Fixed 2,000 (28.3%)                |2,000         |28.3%          |
|10m   |120_d    |7,068           |Fixed 3,000 (42.4%)                |3,000         |42.4%          |
|10m   |120_d    |7,068           |Fixed 4,000 (56.6%)                |4,000         |56.6%          |
|10m   |1305_c   |8,027           |1% of pixels = 80 pixels (1.0%)    |80            |1.0%           |
|10m   |1305_c   |8,027           |2% of pixels = 161 pixels (2.0%)   |161           |2.0%           |
|10m   |1305_c   |8,027           |3% of pixels = 241 pixels (3.0%)   |241           |3.0%           |
|10m   |1305_c   |8,027           |Fixed 1,250 (15.6%)                |1,250         |15.6%          |
|10m   |1305_c   |8,027           |Fixed 2,000 (24.9%)                |2,000         |24.9%          |
|10m   |1305_c   |8,027           |Fixed 3,000 (37.4%)                |3,000         |37.4%          |
|10m   |1305_c   |8,027           |Fixed 4,000 (49.8%)                |4,000         |49.8%          |
|10m   |1315_b   |9,077           |1% of pixels = 91 pixels (1.0%)    |91            |1.0%           |
|10m   |1315_b   |9,077           |2% of pixels = 182 pixels (2.0%)   |182           |2.0%           |
|10m   |1315_b   |9,077           |3% of pixels = 272 pixels (3.0%)   |272           |3.0%           |
|10m   |1315_b   |9,077           |Fixed 1,250 (13.8%)                |1,250         |13.8%          |
|10m   |1315_b   |9,077           |Fixed 2,000 (22.0%)                |2,000         |22.0%          |
|10m   |1315_b   |9,077           |Fixed 3,000 (33.1%)                |3,000         |33.1%          |
|10m   |1315_b   |9,077           |Fixed 4,000 (44.1%)                |4,000         |44.1%          |
|10m   |1400_b   |6,806           |1% of pixels = 68 pixels (1.0%)    |68            |1.0%           |
|10m   |1400_b   |6,806           |2% of pixels = 136 pixels (2.0%)   |136           |2.0%           |
|10m   |1400_b   |6,806           |3% of pixels = 204 pixels (3.0%)   |204           |3.0%           |
|10m   |1400_b   |6,806           |Fixed 1,250 (18.4%)                |1,250         |18.4%          |
|10m   |1400_b   |6,806           |Fixed 2,000 (29.4%)                |2,000         |29.4%          |
|10m   |1400_b   |6,806           |Fixed 3,000 (44.1%)                |3,000         |44.1%          |
|10m   |1400_b   |6,806           |Fixed 4,000 (58.8%)                |4,000         |58.8%          |
|10m   |1414_a   |6,227           |1% of pixels = 62 pixels (1.0%)    |62            |1.0%           |
|10m   |1414_a   |6,227           |2% of pixels = 125 pixels (2.0%)   |125           |2.0%           |
|10m   |1414_a   |6,227           |3% of pixels = 187 pixels (3.0%)   |187           |3.0%           |
|10m   |1414_a   |6,227           |Fixed 1,250 (20.1%)                |1,250         |20.1%          |
|10m   |1414_a   |6,227           |Fixed 2,000 (32.1%)                |2,000         |32.1%          |
|10m   |1414_a   |6,227           |Fixed 3,000 (48.2%)                |3,000         |48.2%          |
|10m   |1414_a   |6,227           |Fixed 4,000 (64.2%)                |4,000         |64.2%          |
|10m   |1417_c   |7,513           |1% of pixels = 75 pixels (1.0%)    |75            |1.0%           |
|10m   |1417_c   |7,513           |2% of pixels = 150 pixels (2.0%)   |150           |2.0%           |
|10m   |1417_c   |7,513           |3% of pixels = 225 pixels (3.0%)   |225           |3.0%           |
|10m   |1417_c   |7,513           |Fixed 1,250 (16.6%)                |1,250         |16.6%          |
|10m   |1417_c   |7,513           |Fixed 2,000 (26.6%)                |2,000         |26.6%          |
|10m   |1417_c   |7,513           |Fixed 3,000 (39.9%)                |3,000         |39.9%          |
|10m   |1417_c   |7,513           |Fixed 4,000 (53.2%)                |4,000         |53.2%          |
|10m   |1514_a   |7,027           |1% of pixels = 70 pixels (1.0%)    |70            |1.0%           |
|10m   |1514_a   |7,027           |2% of pixels = 141 pixels (2.0%)   |141           |2.0%           |
|10m   |1514_a   |7,027           |3% of pixels = 211 pixels (3.0%)   |211           |3.0%           |
|10m   |1514_a   |7,027           |Fixed 1,250 (17.8%)                |1,250         |17.8%          |
|10m   |1514_a   |7,027           |Fixed 2,000 (28.5%)                |2,000         |28.5%          |
|10m   |1514_a   |7,027           |Fixed 3,000 (42.7%)                |3,000         |42.7%          |
|10m   |1514_a   |7,027           |Fixed 4,000 (56.9%)                |4,000         |56.9%          |
|10m   |1516_c   |8,358           |1% of pixels = 84 pixels (1.0%)    |84            |1.0%           |
|10m   |1516_c   |8,358           |2% of pixels = 167 pixels (2.0%)   |167           |2.0%           |
|10m   |1516_c   |8,358           |3% of pixels = 251 pixels (3.0%)   |251           |3.0%           |
|10m   |1516_c   |8,358           |Fixed 1,250 (15.0%)                |1,250         |15.0%          |
|10m   |1516_c   |8,358           |Fixed 2,000 (23.9%)                |2,000         |23.9%          |
|10m   |1516_c   |8,358           |Fixed 3,000 (35.9%)                |3,000         |35.9%          |
|10m   |1516_c   |8,358           |Fixed 4,000 (47.9%)                |4,000         |47.9%          |
|10m   |1604_c   |7,229           |1% of pixels = 72 pixels (1.0%)    |72            |1.0%           |
|10m   |1604_c   |7,229           |2% of pixels = 145 pixels (2.0%)   |145           |2.0%           |
|10m   |1604_c   |7,229           |3% of pixels = 217 pixels (3.0%)   |217           |3.0%           |
|10m   |1604_c   |7,229           |Fixed 1,250 (17.3%)                |1,250         |17.3%          |
|10m   |1604_c   |7,229           |Fixed 2,000 (27.7%)                |2,000         |27.7%          |
|10m   |1604_c   |7,229           |Fixed 3,000 (41.5%)                |3,000         |41.5%          |
|10m   |1604_c   |7,229           |Fixed 4,000 (55.3%)                |4,000         |55.3%          |
|10m   |1606_c   |7,594           |1% of pixels = 76 pixels (1.0%)    |76            |1.0%           |
|10m   |1606_c   |7,594           |2% of pixels = 152 pixels (2.0%)   |152           |2.0%           |
|10m   |1606_c   |7,594           |3% of pixels = 228 pixels (3.0%)   |228           |3.0%           |
|10m   |1606_c   |7,594           |Fixed 1,250 (16.5%)                |1,250         |16.5%          |
|10m   |1606_c   |7,594           |Fixed 2,000 (26.3%)                |2,000         |26.3%          |
|10m   |1606_c   |7,594           |Fixed 3,000 (39.5%)                |3,000         |39.5%          |
|10m   |1606_c   |7,594           |Fixed 4,000 (52.7%)                |4,000         |52.7%          |
|10m   |1701_b   |9,299           |1% of pixels = 93 pixels (1.0%)    |93            |1.0%           |
|10m   |1701_b   |9,299           |2% of pixels = 186 pixels (2.0%)   |186           |2.0%           |
|10m   |1701_b   |9,299           |3% of pixels = 279 pixels (3.0%)   |279           |3.0%           |
|10m   |1701_b   |9,299           |Fixed 1,250 (13.4%)                |1,250         |13.4%          |
|10m   |1701_b   |9,299           |Fixed 2,000 (21.5%)                |2,000         |21.5%          |
|10m   |1701_b   |9,299           |Fixed 3,000 (32.3%)                |3,000         |32.3%          |
|10m   |1701_b   |9,299           |Fixed 4,000 (43.0%)                |4,000         |43.0%          |
|10m   |1803_c   |5,295           |1% of pixels = 53 pixels (1.0%)    |53            |1.0%           |
|10m   |1803_c   |5,295           |2% of pixels = 106 pixels (2.0%)   |106           |2.0%           |
|10m   |1803_c   |5,295           |3% of pixels = 159 pixels (3.0%)   |159           |3.0%           |
|10m   |1803_c   |5,295           |Fixed 1,250 (23.6%)                |1,250         |23.6%          |
|10m   |1803_c   |5,295           |Fixed 2,000 (37.8%)                |2,000         |37.8%          |
|10m   |1803_c   |5,295           |Fixed 3,000 (56.7%)                |3,000         |56.7%          |
|10m   |1803_c   |5,295           |Fixed 4,000 (75.5%)                |4,000         |75.5%          |
|10m   |1816_c   |5,193           |1% of pixels = 52 pixels (1.0%)    |52            |1.0%           |
|10m   |1816_c   |5,193           |2% of pixels = 104 pixels (2.0%)   |104           |2.0%           |
|10m   |1816_c   |5,193           |3% of pixels = 156 pixels (3.0%)   |156           |3.0%           |
|10m   |1816_c   |5,193           |Fixed 1,250 (24.1%)                |1,250         |24.1%          |
|10m   |1816_c   |5,193           |Fixed 2,000 (38.5%)                |2,000         |38.5%          |
|10m   |1816_c   |5,193           |Fixed 3,000 (57.8%)                |3,000         |57.8%          |
|10m   |1816_c   |5,193           |Fixed 4,000 (77.0%)                |4,000         |77.0%          |
|10m   |1912_c   |5,467           |1% of pixels = 55 pixels (1.0%)    |55            |1.0%           |
|10m   |1912_c   |5,467           |2% of pixels = 109 pixels (2.0%)   |109           |2.0%           |
|10m   |1912_c   |5,467           |3% of pixels = 164 pixels (3.0%)   |164           |3.0%           |
|10m   |1912_c   |5,467           |Fixed 1,250 (22.9%)                |1,250         |22.9%          |
|10m   |1912_c   |5,467           |Fixed 2,000 (36.6%)                |2,000         |36.6%          |
|10m   |1912_c   |5,467           |Fixed 3,000 (54.9%)                |3,000         |54.9%          |
|10m   |1912_c   |5,467           |Fixed 4,000 (73.2%)                |4,000         |73.2%          |
|10m   |1917_c   |5,452           |1% of pixels = 55 pixels (1.0%)    |55            |1.0%           |
|10m   |1917_c   |5,452           |2% of pixels = 109 pixels (2.0%)   |109           |2.0%           |
|10m   |1917_c   |5,452           |3% of pixels = 164 pixels (3.0%)   |164           |3.0%           |
|10m   |1917_c   |5,452           |Fixed 1,250 (22.9%)                |1,250         |22.9%          |
|10m   |1917_c   |5,452           |Fixed 2,000 (36.7%)                |2,000         |36.7%          |
|10m   |1917_c   |5,452           |Fixed 3,000 (55.0%)                |3,000         |55.0%          |
|10m   |1917_c   |5,452           |Fixed 4,000 (73.4%)                |4,000         |73.4%          |
|10m   |204_d    |7,147           |1% of pixels = 71 pixels (1.0%)    |71            |1.0%           |
|10m   |204_d    |7,147           |2% of pixels = 143 pixels (2.0%)   |143           |2.0%           |
|10m   |204_d    |7,147           |3% of pixels = 214 pixels (3.0%)   |214           |3.0%           |
|10m   |204_d    |7,147           |Fixed 1,250 (17.5%)                |1,250         |17.5%          |
|10m   |204_d    |7,147           |Fixed 2,000 (28.0%)                |2,000         |28.0%          |
|10m   |204_d    |7,147           |Fixed 3,000 (42.0%)                |3,000         |42.0%          |
|10m   |204_d    |7,147           |Fixed 4,000 (56.0%)                |4,000         |56.0%          |
|10m   |217_d    |7,850           |1% of pixels = 78 pixels (1.0%)    |78            |1.0%           |
|10m   |217_d    |7,850           |2% of pixels = 157 pixels (2.0%)   |157           |2.0%           |
|10m   |217_d    |7,850           |3% of pixels = 236 pixels (3.0%)   |236           |3.0%           |
|10m   |217_d    |7,850           |Fixed 1,250 (15.9%)                |1,250         |15.9%          |
|10m   |217_d    |7,850           |Fixed 2,000 (25.5%)                |2,000         |25.5%          |
|10m   |217_d    |7,850           |Fixed 3,000 (38.2%)                |3,000         |38.2%          |
|10m   |217_d    |7,850           |Fixed 4,000 (51.0%)                |4,000         |51.0%          |
|10m   |22_c     |7,136           |1% of pixels = 71 pixels (1.0%)    |71            |1.0%           |
|10m   |22_c     |7,136           |2% of pixels = 143 pixels (2.0%)   |143           |2.0%           |
|10m   |22_c     |7,136           |3% of pixels = 214 pixels (3.0%)   |214           |3.0%           |
|10m   |22_c     |7,136           |Fixed 1,250 (17.5%)                |1,250         |17.5%          |
|10m   |22_c     |7,136           |Fixed 2,000 (28.0%)                |2,000         |28.0%          |
|10m   |22_c     |7,136           |Fixed 3,000 (42.0%)                |3,000         |42.0%          |
|10m   |22_c     |7,136           |Fixed 4,000 (56.1%)                |4,000         |56.1%          |
|10m   |319_c    |8,983           |1% of pixels = 90 pixels (1.0%)    |90            |1.0%           |
|10m   |319_c    |8,983           |2% of pixels = 180 pixels (2.0%)   |180           |2.0%           |
|10m   |319_c    |8,983           |3% of pixels = 269 pixels (3.0%)   |269           |3.0%           |
|10m   |319_c    |8,983           |Fixed 1,250 (13.9%)                |1,250         |13.9%          |
|10m   |319_c    |8,983           |Fixed 2,000 (22.3%)                |2,000         |22.3%          |
|10m   |319_c    |8,983           |Fixed 3,000 (33.4%)                |3,000         |33.4%          |
|10m   |319_c    |8,983           |Fixed 4,000 (44.5%)                |4,000         |44.5%          |
|10m   |409_d    |5,989           |1% of pixels = 60 pixels (1.0%)    |60            |1.0%           |
|10m   |409_d    |5,989           |2% of pixels = 120 pixels (2.0%)   |120           |2.0%           |
|10m   |409_d    |5,989           |3% of pixels = 180 pixels (3.0%)   |180           |3.0%           |
|10m   |409_d    |5,989           |Fixed 1,250 (20.9%)                |1,250         |20.9%          |
|10m   |409_d    |5,989           |Fixed 2,000 (33.4%)                |2,000         |33.4%          |
|10m   |409_d    |5,989           |Fixed 3,000 (50.1%)                |3,000         |50.1%          |
|10m   |409_d    |5,989           |Fixed 4,000 (66.8%)                |4,000         |66.8%          |
|10m   |503_c    |7,403           |1% of pixels = 74 pixels (1.0%)    |74            |1.0%           |
|10m   |503_c    |7,403           |2% of pixels = 148 pixels (2.0%)   |148           |2.0%           |
|10m   |503_c    |7,403           |3% of pixels = 222 pixels (3.0%)   |222           |3.0%           |
|10m   |503_c    |7,403           |Fixed 1,250 (16.9%)                |1,250         |16.9%          |
|10m   |503_c    |7,403           |Fixed 2,000 (27.0%)                |2,000         |27.0%          |
|10m   |503_c    |7,403           |Fixed 3,000 (40.5%)                |3,000         |40.5%          |
|10m   |503_c    |7,403           |Fixed 4,000 (54.0%)                |4,000         |54.0%          |
|10m   |614_a    |7,102           |1% of pixels = 71 pixels (1.0%)    |71            |1.0%           |
|10m   |614_a    |7,102           |2% of pixels = 142 pixels (2.0%)   |142           |2.0%           |
|10m   |614_a    |7,102           |3% of pixels = 213 pixels (3.0%)   |213           |3.0%           |
|10m   |614_a    |7,102           |Fixed 1,250 (17.6%)                |1,250         |17.6%          |
|10m   |614_a    |7,102           |Fixed 2,000 (28.2%)                |2,000         |28.2%          |
|10m   |614_a    |7,102           |Fixed 3,000 (42.2%)                |3,000         |42.2%          |
|10m   |614_a    |7,102           |Fixed 4,000 (56.3%)                |4,000         |56.3%          |
|10m   |700_c    |7,623           |1% of pixels = 76 pixels (1.0%)    |76            |1.0%           |
|10m   |700_c    |7,623           |2% of pixels = 152 pixels (2.0%)   |152           |2.0%           |
|10m   |700_c    |7,623           |3% of pixels = 229 pixels (3.0%)   |229           |3.0%           |
|10m   |700_c    |7,623           |Fixed 1,250 (16.4%)                |1,250         |16.4%          |
|10m   |700_c    |7,623           |Fixed 2,000 (26.2%)                |2,000         |26.2%          |
|10m   |700_c    |7,623           |Fixed 3,000 (39.4%)                |3,000         |39.4%          |
|10m   |700_c    |7,623           |Fixed 4,000 (52.5%)                |4,000         |52.5%          |
|10m   |722_a    |8,781           |1% of pixels = 88 pixels (1.0%)    |88            |1.0%           |
|10m   |722_a    |8,781           |2% of pixels = 176 pixels (2.0%)   |176           |2.0%           |
|10m   |722_a    |8,781           |3% of pixels = 263 pixels (3.0%)   |263           |3.0%           |
|10m   |722_a    |8,781           |Fixed 1,250 (14.2%)                |1,250         |14.2%          |
|10m   |722_a    |8,781           |Fixed 2,000 (22.8%)                |2,000         |22.8%          |
|10m   |722_a    |8,781           |Fixed 3,000 (34.2%)                |3,000         |34.2%          |
|10m   |722_a    |8,781           |Fixed 4,000 (45.6%)                |4,000         |45.6%          |
|10m   |800_a    |3,976           |1% of pixels = 40 pixels (1.0%)    |40            |1.0%           |
|10m   |800_a    |3,976           |2% of pixels = 80 pixels (2.0%)    |80            |2.0%           |
|10m   |800_a    |3,976           |3% of pixels = 119 pixels (3.0%)   |119           |3.0%           |
|10m   |800_a    |3,976           |Fixed 1,250 (31.4%)                |1,250         |31.4%          |
|10m   |800_a    |3,976           |Fixed 2,000 (50.3%)                |2,000         |50.3%          |
|10m   |800_a    |3,976           |Fixed 3,000 (75.5%)                |3,000         |75.5%          |
|10m   |800_a    |3,976           |Fixed 4,000 (100.0%)               |3,976         |100.0%         |
|10m   |914_a    |7,220           |1% of pixels = 72 pixels (1.0%)    |72            |1.0%           |
|10m   |914_a    |7,220           |2% of pixels = 144 pixels (2.0%)   |144           |2.0%           |
|10m   |914_a    |7,220           |3% of pixels = 217 pixels (3.0%)   |217           |3.0%           |
|10m   |914_a    |7,220           |Fixed 1,250 (17.3%)                |1,250         |17.3%          |
|10m   |914_a    |7,220           |Fixed 2,000 (27.7%)                |2,000         |27.7%          |
|10m   |914_a    |7,220           |Fixed 3,000 (41.6%)                |3,000         |41.6%          |
|10m   |914_a    |7,220           |Fixed 4,000 (55.4%)                |4,000         |55.4%          |
|20m   |100      |30,808          |1% of pixels = 308 pixels (1.0%)   |308           |1.0%           |
|20m   |100      |30,808          |2% of pixels = 616 pixels (2.0%)   |616           |2.0%           |
|20m   |100      |30,808          |3% of pixels = 924 pixels (3.0%)   |924           |3.0%           |
|20m   |100      |30,808          |Fixed 1,250 (4.1%)                 |1,250         |4.1%           |
|20m   |100      |30,808          |Fixed 2,000 (6.5%)                 |2,000         |6.5%           |
|20m   |100      |30,808          |Fixed 3,000 (9.7%)                 |3,000         |9.7%           |
|20m   |100      |30,808          |Fixed 4,000 (13.0%)                |4,000         |13.0%          |
|20m   |1201     |24,803          |1% of pixels = 248 pixels (1.0%)   |248           |1.0%           |
|20m   |1201     |24,803          |2% of pixels = 496 pixels (2.0%)   |496           |2.0%           |
|20m   |1201     |24,803          |3% of pixels = 744 pixels (3.0%)   |744           |3.0%           |
|20m   |1201     |24,803          |Fixed 1,250 (5.0%)                 |1,250         |5.0%           |
|20m   |1201     |24,803          |Fixed 2,000 (8.1%)                 |2,000         |8.1%           |
|20m   |1201     |24,803          |Fixed 3,000 (12.1%)                |3,000         |12.1%          |
|20m   |1201     |24,803          |Fixed 4,000 (16.1%)                |4,000         |16.1%          |
|20m   |1402     |27,061          |1% of pixels = 271 pixels (1.0%)   |271           |1.0%           |
|20m   |1402     |27,061          |2% of pixels = 541 pixels (2.0%)   |541           |2.0%           |
|20m   |1402     |27,061          |3% of pixels = 812 pixels (3.0%)   |812           |3.0%           |
|20m   |1402     |27,061          |Fixed 1,250 (4.6%)                 |1,250         |4.6%           |
|20m   |1402     |27,061          |Fixed 2,000 (7.4%)                 |2,000         |7.4%           |
|20m   |1402     |27,061          |Fixed 3,000 (11.1%)                |3,000         |11.1%          |
|20m   |1402     |27,061          |Fixed 4,000 (14.8%)                |4,000         |14.8%          |
|20m   |1512     |30,984          |1% of pixels = 310 pixels (1.0%)   |310           |1.0%           |
|20m   |1512     |30,984          |2% of pixels = 620 pixels (2.0%)   |620           |2.0%           |
|20m   |1512     |30,984          |3% of pixels = 930 pixels (3.0%)   |930           |3.0%           |
|20m   |1512     |30,984          |Fixed 1,250 (4.0%)                 |1,250         |4.0%           |
|20m   |1512     |30,984          |Fixed 2,000 (6.5%)                 |2,000         |6.5%           |
|20m   |1512     |30,984          |Fixed 3,000 (9.7%)                 |3,000         |9.7%           |
|20m   |1512     |30,984          |Fixed 4,000 (12.9%)                |4,000         |12.9%          |
|20m   |1518     |35,713          |1% of pixels = 357 pixels (1.0%)   |357           |1.0%           |
|20m   |1518     |35,713          |2% of pixels = 714 pixels (2.0%)   |714           |2.0%           |
|20m   |1518     |35,713          |3% of pixels = 1,071 pixels (3.0%) |1,071         |3.0%           |
|20m   |1518     |35,713          |Fixed 1,250 (3.5%)                 |1,250         |3.5%           |
|20m   |1518     |35,713          |Fixed 2,000 (5.6%)                 |2,000         |5.6%           |
|20m   |1518     |35,713          |Fixed 3,000 (8.4%)                 |3,000         |8.4%           |
|20m   |1518     |35,713          |Fixed 4,000 (11.2%)                |4,000         |11.2%          |
|20m   |1800     |29,990          |1% of pixels = 300 pixels (1.0%)   |300           |1.0%           |
|20m   |1800     |29,990          |2% of pixels = 600 pixels (2.0%)   |600           |2.0%           |
|20m   |1800     |29,990          |3% of pixels = 900 pixels (3.0%)   |900           |3.0%           |
|20m   |1800     |29,990          |Fixed 1,250 (4.2%)                 |1,250         |4.2%           |
|20m   |1800     |29,990          |Fixed 2,000 (6.7%)                 |2,000         |6.7%           |
|20m   |1800     |29,990          |Fixed 3,000 (10.0%)                |3,000         |10.0%          |
|20m   |1800     |29,990          |Fixed 4,000 (13.3%)                |4,000         |13.3%          |
|20m   |1805     |25,420          |1% of pixels = 254 pixels (1.0%)   |254           |1.0%           |
|20m   |1805     |25,420          |2% of pixels = 508 pixels (2.0%)   |508           |2.0%           |
|20m   |1805     |25,420          |3% of pixels = 763 pixels (3.0%)   |763           |3.0%           |
|20m   |1805     |25,420          |Fixed 1,250 (4.9%)                 |1,250         |4.9%           |
|20m   |1805     |25,420          |Fixed 2,000 (7.9%)                 |2,000         |7.9%           |
|20m   |1805     |25,420          |Fixed 3,000 (11.8%)                |3,000         |11.8%          |
|20m   |1805     |25,420          |Fixed 4,000 (15.7%)                |4,000         |15.7%          |
|20m   |1910     |25,613          |1% of pixels = 256 pixels (1.0%)   |256           |1.0%           |
|20m   |1910     |25,613          |2% of pixels = 512 pixels (2.0%)   |512           |2.0%           |
|20m   |1910     |25,613          |3% of pixels = 768 pixels (3.0%)   |768           |3.0%           |
|20m   |1910     |25,613          |Fixed 1,250 (4.9%)                 |1,250         |4.9%           |
|20m   |1910     |25,613          |Fixed 2,000 (7.8%)                 |2,000         |7.8%           |
|20m   |1910     |25,613          |Fixed 3,000 (11.7%)                |3,000         |11.7%          |
|20m   |1910     |25,613          |Fixed 4,000 (15.6%)                |4,000         |15.6%          |
|20m   |206      |23,348          |1% of pixels = 233 pixels (1.0%)   |233           |1.0%           |
|20m   |206      |23,348          |2% of pixels = 467 pixels (2.0%)   |467           |2.0%           |
|20m   |206      |23,348          |3% of pixels = 700 pixels (3.0%)   |700           |3.0%           |
|20m   |206      |23,348          |Fixed 1,250 (5.4%)                 |1,250         |5.4%           |
|20m   |206      |23,348          |Fixed 2,000 (8.6%)                 |2,000         |8.6%           |
|20m   |206      |23,348          |Fixed 3,000 (12.8%)                |3,000         |12.8%          |
|20m   |206      |23,348          |Fixed 4,000 (17.1%)                |4,000         |17.1%          |
|20m   |219      |29,380          |1% of pixels = 294 pixels (1.0%)   |294           |1.0%           |
|20m   |219      |29,380          |2% of pixels = 588 pixels (2.0%)   |588           |2.0%           |
|20m   |219      |29,380          |3% of pixels = 881 pixels (3.0%)   |881           |3.0%           |
|20m   |219      |29,380          |Fixed 1,250 (4.3%)                 |1,250         |4.3%           |
|20m   |219      |29,380          |Fixed 2,000 (6.8%)                 |2,000         |6.8%           |
|20m   |219      |29,380          |Fixed 3,000 (10.2%)                |3,000         |10.2%          |
|20m   |219      |29,380          |Fixed 4,000 (13.6%)                |4,000         |13.6%          |
|20m   |302      |26,326          |1% of pixels = 263 pixels (1.0%)   |263           |1.0%           |
|20m   |302      |26,326          |2% of pixels = 527 pixels (2.0%)   |527           |2.0%           |
|20m   |302      |26,326          |3% of pixels = 790 pixels (3.0%)   |790           |3.0%           |
|20m   |302      |26,326          |Fixed 1,250 (4.7%)                 |1,250         |4.7%           |
|20m   |302      |26,326          |Fixed 2,000 (7.6%)                 |2,000         |7.6%           |
|20m   |302      |26,326          |Fixed 3,000 (11.4%)                |3,000         |11.4%          |
|20m   |302      |26,326          |Fixed 4,000 (15.2%)                |4,000         |15.2%          |
|20m   |317      |31,152          |1% of pixels = 312 pixels (1.0%)   |312           |1.0%           |
|20m   |317      |31,152          |2% of pixels = 623 pixels (2.0%)   |623           |2.0%           |
|20m   |317      |31,152          |3% of pixels = 935 pixels (3.0%)   |935           |3.0%           |
|20m   |317      |31,152          |Fixed 1,250 (4.0%)                 |1,250         |4.0%           |
|20m   |317      |31,152          |Fixed 2,000 (6.4%)                 |2,000         |6.4%           |
|20m   |317      |31,152          |Fixed 3,000 (9.6%)                 |3,000         |9.6%           |
|20m   |317      |31,152          |Fixed 4,000 (12.8%)                |4,000         |12.8%          |
|20m   |405      |36,367          |1% of pixels = 364 pixels (1.0%)   |364           |1.0%           |
|20m   |405      |36,367          |2% of pixels = 727 pixels (2.0%)   |727           |2.0%           |
|20m   |405      |36,367          |3% of pixels = 1,091 pixels (3.0%) |1,091         |3.0%           |
|20m   |405      |36,367          |Fixed 1,250 (3.4%)                 |1,250         |3.4%           |
|20m   |405      |36,367          |Fixed 2,000 (5.5%)                 |2,000         |5.5%           |
|20m   |405      |36,367          |Fixed 3,000 (8.2%)                 |3,000         |8.2%           |
|20m   |405      |36,367          |Fixed 4,000 (11.0%)                |4,000         |11.0%          |
|20m   |821      |29,194          |1% of pixels = 292 pixels (1.0%)   |292           |1.0%           |
|20m   |821      |29,194          |2% of pixels = 584 pixels (2.0%)   |584           |2.0%           |
|20m   |821      |29,194          |3% of pixels = 876 pixels (3.0%)   |876           |3.0%           |
|20m   |821      |29,194          |Fixed 1,250 (4.3%)                 |1,250         |4.3%           |
|20m   |821      |29,194          |Fixed 2,000 (6.9%)                 |2,000         |6.9%           |
|20m   |821      |29,194          |Fixed 3,000 (10.3%)                |3,000         |10.3%          |
|20m   |821      |29,194          |Fixed 4,000 (13.7%)                |4,000         |13.7%          |
|20m   |905      |27,021          |1% of pixels = 270 pixels (1.0%)   |270           |1.0%           |
|20m   |905      |27,021          |2% of pixels = 540 pixels (2.0%)   |540           |2.0%           |
|20m   |905      |27,021          |3% of pixels = 811 pixels (3.0%)   |811           |3.0%           |
|20m   |905      |27,021          |Fixed 1,250 (4.6%)                 |1,250         |4.6%           |
|20m   |905      |27,021          |Fixed 2,000 (7.4%)                 |2,000         |7.4%           |
|20m   |905      |27,021          |Fixed 3,000 (11.1%)                |3,000         |11.1%          |
|20m   |905      |27,021          |Fixed 4,000 (14.8%)                |4,000         |14.8%          |
|20m   |912      |29,618          |1% of pixels = 296 pixels (1.0%)   |296           |1.0%           |
|20m   |912      |29,618          |2% of pixels = 592 pixels (2.0%)   |592           |2.0%           |
|20m   |912      |29,618          |3% of pixels = 889 pixels (3.0%)   |889           |3.0%           |
|20m   |912      |29,618          |Fixed 1,250 (4.2%)                 |1,250         |4.2%           |
|20m   |912      |29,618          |Fixed 2,000 (6.8%)                 |2,000         |6.8%           |
|20m   |912      |29,618          |Fixed 3,000 (10.1%)                |3,000         |10.1%          |
|20m   |912      |29,618          |Fixed 4,000 (13.5%)                |4,000         |13.5%          |
|50m   |sub50_15 |172,782         |Fixed 1,250 (0.7%)                 |1,250         |0.7%           |
|50m   |sub50_15 |172,782         |1% of pixels = 1,728 pixels (1.0%) |1,728         |1.0%           |
|50m   |sub50_15 |172,782         |2% of pixels = 3,456 pixels (2.0%) |3,456         |2.0%           |
|50m   |sub50_15 |172,782         |Fixed 4,000 (2.3%)                 |4,000         |2.3%           |
|50m   |sub50_15 |172,782         |3% of pixels = 5,000 pixels (2.9%) |5,000         |2.9%           |
|50m   |sub50_15 |172,782         |Fixed 6,000 (3.5%)                 |6,000         |3.5%           |
|50m   |sub50_15 |172,782         |Fixed 8,000 (4.6%)                 |8,000         |4.6%           |
|50m   |sub50_24 |166,816         |Fixed 1,250 (0.7%)                 |1,250         |0.7%           |
|50m   |sub50_24 |166,816         |1% of pixels = 1,668 pixels (1.0%) |1,668         |1.0%           |
|50m   |sub50_24 |166,816         |2% of pixels = 3,336 pixels (2.0%) |3,336         |2.0%           |
|50m   |sub50_24 |166,816         |Fixed 4,000 (2.4%)                 |4,000         |2.4%           |
|50m   |sub50_24 |166,816         |3% of pixels = 5,000 pixels (3.0%) |5,000         |3.0%           |
|50m   |sub50_24 |166,816         |Fixed 6,000 (3.6%)                 |6,000         |3.6%           |
|50m   |sub50_24 |166,816         |Fixed 8,000 (4.8%)                 |8,000         |4.8%           |
|50m   |sub50_36 |174,203         |Fixed 1,250 (0.7%)                 |1,250         |0.7%           |
|50m   |sub50_36 |174,203         |1% of pixels = 1,742 pixels (1.0%) |1,742         |1.0%           |
|50m   |sub50_36 |174,203         |2% of pixels = 3,484 pixels (2.0%) |3,484         |2.0%           |
|50m   |sub50_36 |174,203         |Fixed 4,000 (2.3%)                 |4,000         |2.3%           |
|50m   |sub50_36 |174,203         |3% of pixels = 5,000 pixels (2.9%) |5,000         |2.9%           |
|50m   |sub50_36 |174,203         |Fixed 6,000 (3.4%)                 |6,000         |3.4%           |
|50m   |sub50_36 |174,203         |Fixed 8,000 (4.6%)                 |8,000         |4.6%           |
|50m   |sub50_49 |206,958         |Fixed 1,250 (0.6%)                 |1,250         |0.6%           |
|50m   |sub50_49 |206,958         |1% of pixels = 2,070 pixels (1.0%) |2,070         |1.0%           |
|50m   |sub50_49 |206,958         |Fixed 4,000 (1.9%)                 |4,000         |1.9%           |
|50m   |sub50_49 |206,958         |2% of pixels = 4,139 pixels (2.0%) |4,139         |2.0%           |
|50m   |sub50_49 |206,958         |3% of pixels = 5,000 pixels (2.4%) |5,000         |2.4%           |
|50m   |sub50_49 |206,958         |Fixed 6,000 (2.9%)                 |6,000         |2.9%           |
|50m   |sub50_49 |206,958         |Fixed 8,000 (3.9%)                 |8,000         |3.9%           |
|50m   |sub50_51 |200,631         |Fixed 1,250 (0.6%)                 |1,250         |0.6%           |
|50m   |sub50_51 |200,631         |1% of pixels = 2,006 pixels (1.0%) |2,006         |1.0%           |
|50m   |sub50_51 |200,631         |Fixed 4,000 (2.0%)                 |4,000         |2.0%           |
|50m   |sub50_51 |200,631         |2% of pixels = 4,013 pixels (2.0%) |4,013         |2.0%           |
|50m   |sub50_51 |200,631         |3% of pixels = 5,000 pixels (2.5%) |5,000         |2.5%           |
|50m   |sub50_51 |200,631         |Fixed 6,000 (3.0%)                 |6,000         |3.0%           |
|50m   |sub50_51 |200,631         |Fixed 8,000 (4.0%)                 |8,000         |4.0%           |
|50m   |sub50_56 |145,500         |Fixed 1,250 (0.9%)                 |1,250         |0.9%           |
|50m   |sub50_56 |145,500         |1% of pixels = 1,455 pixels (1.0%) |1,455         |1.0%           |
|50m   |sub50_56 |145,500         |2% of pixels = 2,910 pixels (2.0%) |2,910         |2.0%           |
|50m   |sub50_56 |145,500         |Fixed 4,000 (2.7%)                 |4,000         |2.7%           |
|50m   |sub50_56 |145,500         |3% of pixels = 4,365 pixels (3.0%) |4,365         |3.0%           |
|50m   |sub50_56 |145,500         |Fixed 6,000 (4.1%)                 |6,000         |4.1%           |
|50m   |sub50_56 |145,500         |Fixed 8,000 (5.5%)                 |8,000         |5.5%           |
|50m   |sub50_60 |164,492         |Fixed 1,250 (0.8%)                 |1,250         |0.8%           |
|50m   |sub50_60 |164,492         |1% of pixels = 1,645 pixels (1.0%) |1,645         |1.0%           |
|50m   |sub50_60 |164,492         |2% of pixels = 3,290 pixels (2.0%) |3,290         |2.0%           |
|50m   |sub50_60 |164,492         |Fixed 4,000 (2.4%)                 |4,000         |2.4%           |
|50m   |sub50_60 |164,492         |3% of pixels = 4,935 pixels (3.0%) |4,935         |3.0%           |
|50m   |sub50_60 |164,492         |Fixed 6,000 (3.6%)                 |6,000         |3.6%           |
|50m   |sub50_60 |164,492         |Fixed 8,000 (4.9%)                 |8,000         |4.9%           |
|50m   |sub50_8  |187,893         |Fixed 1,250 (0.7%)                 |1,250         |0.7%           |
|50m   |sub50_8  |187,893         |1% of pixels = 1,879 pixels (1.0%) |1,879         |1.0%           |
|50m   |sub50_8  |187,893         |2% of pixels = 3,758 pixels (2.0%) |3,758         |2.0%           |
|50m   |sub50_8  |187,893         |Fixed 4,000 (2.1%)                 |4,000         |2.1%           |
|50m   |sub50_8  |187,893         |3% of pixels = 5,000 pixels (2.7%) |5,000         |2.7%           |
|50m   |sub50_8  |187,893         |Fixed 6,000 (3.2%)                 |6,000         |3.2%           |
|50m   |sub50_8  |187,893         |Fixed 8,000 (4.3%)                 |8,000         |4.3%           |

## Summary Results



|Metric            |Scale |Quadrat  |Sample rule                        |Mean value |Boot SD |Boot CV |95% CI low |95% CI high |Delta from fixed 4,000 |Calculation method                                                        |
|:-----------------|:-----|:--------|:----------------------------------|:----------|:-------|:-------|:----------|:-----------|:----------------------|:-------------------------------------------------------------------------|
|PCA Mean Distance |10m   |1104_b   |1% of pixels = 94 pixels (1.0%)    |6.293      |0.512   |8.1%    |6.148      |6.439       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1104_b   |2% of pixels = 188 pixels (2.0%)   |6.273      |0.320   |5.1%    |6.182      |6.364       |-0.022                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1104_b   |3% of pixels = 283 pixels (3.0%)   |6.308      |0.290   |4.6%    |6.225      |6.390       |0.013                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1104_b   |Fixed 1,250 (13.3%)                |6.252      |0.107   |1.7%    |6.221      |6.282       |-0.043                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1104_b   |Fixed 2,000 (21.2%)                |6.285      |0.080   |1.3%    |6.262      |6.308       |-0.010                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1104_b   |Fixed 3,000 (31.8%)                |6.276      |0.065   |1.0%    |6.258      |6.294       |-0.019                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1104_b   |Fixed 4,000 (42.4%)                |6.295      |0.061   |1.0%    |6.277      |6.312       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1110_a   |1% of pixels = 76 pixels (1.0%)    |6.476      |0.380   |5.9%    |6.368      |6.584       |-0.055                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1110_a   |2% of pixels = 153 pixels (2.0%)   |6.501      |0.312   |4.8%    |6.412      |6.589       |-0.031                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1110_a   |3% of pixels = 230 pixels (3.0%)   |6.485      |0.166   |2.6%    |6.438      |6.533       |-0.046                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1110_a   |Fixed 1,250 (16.3%)                |6.543      |0.103   |1.6%    |6.514      |6.573       |0.012                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1110_a   |Fixed 2,000 (26.1%)                |6.529      |0.068   |1.0%    |6.510      |6.548       |-0.002                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1110_a   |Fixed 3,000 (39.2%)                |6.527      |0.051   |0.8%    |6.513      |6.542       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1110_a   |Fixed 4,000 (52.3%)                |6.531      |0.035   |0.5%    |6.521      |6.541       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1124_a   |1% of pixels = 83 pixels (1.0%)    |7.007      |0.501   |7.2%    |6.865      |7.150       |-0.082                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1124_a   |2% of pixels = 166 pixels (2.0%)   |7.085      |0.369   |5.2%    |6.980      |7.190       |-0.005                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1124_a   |3% of pixels = 249 pixels (3.0%)   |7.073      |0.302   |4.3%    |6.987      |7.159       |-0.017                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1124_a   |Fixed 1,250 (15.1%)                |7.105      |0.104   |1.5%    |7.075      |7.134       |0.015                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1124_a   |Fixed 2,000 (24.1%)                |7.071      |0.093   |1.3%    |7.044      |7.097       |-0.019                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1124_a   |Fixed 3,000 (36.1%)                |7.124      |0.064   |0.9%    |7.105      |7.142       |0.034                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1124_a   |Fixed 4,000 (48.2%)                |7.090      |0.056   |0.8%    |7.074      |7.106       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |112_b    |1% of pixels = 73 pixels (1.0%)    |8.564      |1.113   |13.0%   |8.247      |8.880       |0.011                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |112_b    |2% of pixels = 147 pixels (2.0%)   |8.584      |0.860   |10.0%   |8.340      |8.829       |0.032                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |112_b    |3% of pixels = 220 pixels (3.0%)   |8.631      |0.815   |9.4%    |8.399      |8.863       |0.079                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |112_b    |Fixed 1,250 (17.0%)                |8.552      |0.274   |3.2%    |8.474      |8.630       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |112_b    |Fixed 2,000 (27.2%)                |8.551      |0.215   |2.5%    |8.490      |8.612       |-0.002                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |112_b    |Fixed 3,000 (40.9%)                |8.508      |0.149   |1.8%    |8.465      |8.550       |-0.045                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |112_b    |Fixed 4,000 (54.5%)                |8.552      |0.119   |1.4%    |8.519      |8.586       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |114_b    |1% of pixels = 72 pixels (1.0%)    |9.334      |0.738   |7.9%    |9.124      |9.544       |-0.053                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |114_b    |2% of pixels = 145 pixels (2.0%)   |9.377      |0.506   |5.4%    |9.233      |9.521       |-0.010                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |114_b    |3% of pixels = 217 pixels (3.0%)   |9.272      |0.339   |3.7%    |9.175      |9.368       |-0.115                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |114_b    |Fixed 1,250 (17.3%)                |9.389      |0.180   |1.9%    |9.337      |9.440       |0.002                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |114_b    |Fixed 2,000 (27.6%)                |9.346      |0.126   |1.3%    |9.311      |9.382       |-0.040                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |114_b    |Fixed 3,000 (41.4%)                |9.410      |0.092   |1.0%    |9.384      |9.436       |0.023                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |114_b    |Fixed 4,000 (55.3%)                |9.387      |0.072   |0.8%    |9.366      |9.407       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |11_c     |1% of pixels = 62 pixels (1.0%)    |6.955      |0.497   |7.1%    |6.814      |7.097       |-0.047                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |11_c     |2% of pixels = 124 pixels (2.0%)   |6.925      |0.323   |4.7%    |6.834      |7.017       |-0.077                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |11_c     |3% of pixels = 186 pixels (3.0%)   |7.060      |0.301   |4.3%    |6.975      |7.146       |0.058                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |11_c     |Fixed 1,250 (20.2%)                |6.969      |0.100   |1.4%    |6.940      |6.997       |-0.034                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |11_c     |Fixed 2,000 (32.3%)                |6.991      |0.068   |1.0%    |6.972      |7.010       |-0.012                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |11_c     |Fixed 3,000 (48.4%)                |6.988      |0.050   |0.7%    |6.974      |7.002       |-0.015                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |11_c     |Fixed 4,000 (64.5%)                |7.003      |0.037   |0.5%    |6.992      |7.013       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |120_d    |1% of pixels = 71 pixels (1.0%)    |8.804      |0.569   |6.5%    |8.642      |8.966       |-0.148                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |120_d    |2% of pixels = 141 pixels (2.0%)   |8.977      |0.403   |4.5%    |8.863      |9.091       |0.025                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |120_d    |3% of pixels = 212 pixels (3.0%)   |8.907      |0.381   |4.3%    |8.798      |9.015       |-0.045                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |120_d    |Fixed 1,250 (17.7%)                |8.959      |0.138   |1.5%    |8.920      |8.998       |0.007                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |120_d    |Fixed 2,000 (28.3%)                |8.947      |0.088   |1.0%    |8.922      |8.972       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |120_d    |Fixed 3,000 (42.4%)                |8.957      |0.068   |0.8%    |8.938      |8.976       |0.005                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |120_d    |Fixed 4,000 (56.6%)                |8.952      |0.054   |0.6%    |8.936      |8.967       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1305_c   |1% of pixels = 80 pixels (1.0%)    |7.200      |0.394   |5.5%    |7.088      |7.312       |-0.072                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1305_c   |2% of pixels = 161 pixels (2.0%)   |7.271      |0.328   |4.5%    |7.178      |7.364       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1305_c   |3% of pixels = 241 pixels (3.0%)   |7.276      |0.251   |3.5%    |7.205      |7.348       |0.005                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1305_c   |Fixed 1,250 (15.6%)                |7.269      |0.116   |1.6%    |7.236      |7.302       |-0.003                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1305_c   |Fixed 2,000 (24.9%)                |7.271      |0.072   |1.0%    |7.250      |7.291       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1305_c   |Fixed 3,000 (37.4%)                |7.266      |0.056   |0.8%    |7.250      |7.282       |-0.006                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1305_c   |Fixed 4,000 (49.8%)                |7.272      |0.038   |0.5%    |7.261      |7.283       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1315_b   |1% of pixels = 91 pixels (1.0%)    |8.018      |0.579   |7.2%    |7.853      |8.182       |0.020                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1315_b   |2% of pixels = 182 pixels (2.0%)   |7.995      |0.375   |4.7%    |7.888      |8.102       |-0.003                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1315_b   |3% of pixels = 272 pixels (3.0%)   |7.913      |0.328   |4.1%    |7.820      |8.007       |-0.084                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1315_b   |Fixed 1,250 (13.8%)                |7.955      |0.138   |1.7%    |7.916      |7.994       |-0.043                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1315_b   |Fixed 2,000 (22.0%)                |7.979      |0.097   |1.2%    |7.951      |8.007       |-0.019                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1315_b   |Fixed 3,000 (33.1%)                |7.967      |0.074   |0.9%    |7.946      |7.988       |-0.030                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1315_b   |Fixed 4,000 (44.1%)                |7.998      |0.075   |0.9%    |7.977      |8.019       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1400_b   |1% of pixels = 68 pixels (1.0%)    |7.799      |0.491   |6.3%    |7.660      |7.939       |-0.128                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1400_b   |2% of pixels = 136 pixels (2.0%)   |7.883      |0.350   |4.4%    |7.784      |7.983       |-0.044                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1400_b   |3% of pixels = 204 pixels (3.0%)   |7.928      |0.274   |3.5%    |7.850      |8.006       |0.001                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1400_b   |Fixed 1,250 (18.4%)                |7.937      |0.111   |1.4%    |7.905      |7.969       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1400_b   |Fixed 2,000 (29.4%)                |7.937      |0.068   |0.9%    |7.918      |7.956       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1400_b   |Fixed 3,000 (44.1%)                |7.916      |0.057   |0.7%    |7.900      |7.932       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1400_b   |Fixed 4,000 (58.8%)                |7.927      |0.039   |0.5%    |7.916      |7.938       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1414_a   |1% of pixels = 62 pixels (1.0%)    |8.896      |0.581   |6.5%    |8.731      |9.061       |0.084                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1414_a   |2% of pixels = 125 pixels (2.0%)   |8.710      |0.426   |4.9%    |8.589      |8.831       |-0.102                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1414_a   |3% of pixels = 187 pixels (3.0%)   |8.796      |0.372   |4.2%    |8.690      |8.902       |-0.016                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1414_a   |Fixed 1,250 (20.1%)                |8.801      |0.150   |1.7%    |8.758      |8.844       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1414_a   |Fixed 2,000 (32.1%)                |8.796      |0.080   |0.9%    |8.773      |8.819       |-0.016                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1414_a   |Fixed 3,000 (48.2%)                |8.823      |0.061   |0.7%    |8.806      |8.841       |0.011                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1414_a   |Fixed 4,000 (64.2%)                |8.812      |0.046   |0.5%    |8.799      |8.825       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1417_c   |1% of pixels = 75 pixels (1.0%)    |9.084      |0.616   |6.8%    |8.908      |9.259       |-0.134                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1417_c   |2% of pixels = 150 pixels (2.0%)   |9.224      |0.448   |4.9%    |9.096      |9.351       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1417_c   |3% of pixels = 225 pixels (3.0%)   |9.189      |0.375   |4.1%    |9.082      |9.296       |-0.029                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1417_c   |Fixed 1,250 (16.6%)                |9.201      |0.147   |1.6%    |9.159      |9.243       |-0.017                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1417_c   |Fixed 2,000 (26.6%)                |9.207      |0.114   |1.2%    |9.174      |9.239       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1417_c   |Fixed 3,000 (39.9%)                |9.198      |0.076   |0.8%    |9.177      |9.220       |-0.019                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1417_c   |Fixed 4,000 (53.2%)                |9.217      |0.058   |0.6%    |9.201      |9.234       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1514_a   |1% of pixels = 70 pixels (1.0%)    |5.882      |0.471   |8.0%    |5.748      |6.016       |-0.041                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1514_a   |2% of pixels = 141 pixels (2.0%)   |5.924      |0.273   |4.6%    |5.846      |6.001       |0.001                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1514_a   |3% of pixels = 211 pixels (3.0%)   |5.908      |0.271   |4.6%    |5.831      |5.985       |-0.015                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1514_a   |Fixed 1,250 (17.8%)                |5.945      |0.100   |1.7%    |5.917      |5.974       |0.022                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1514_a   |Fixed 2,000 (28.5%)                |5.928      |0.072   |1.2%    |5.907      |5.948       |0.005                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1514_a   |Fixed 3,000 (42.7%)                |5.912      |0.049   |0.8%    |5.898      |5.926       |-0.010                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1514_a   |Fixed 4,000 (56.9%)                |5.923      |0.038   |0.6%    |5.912      |5.933       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1516_c   |1% of pixels = 84 pixels (1.0%)    |8.558      |0.567   |6.6%    |8.397      |8.719       |-0.030                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1516_c   |2% of pixels = 167 pixels (2.0%)   |8.524      |0.342   |4.0%    |8.426      |8.621       |-0.065                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1516_c   |3% of pixels = 251 pixels (3.0%)   |8.634      |0.267   |3.1%    |8.559      |8.710       |0.046                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1516_c   |Fixed 1,250 (15.0%)                |8.598      |0.112   |1.3%    |8.566      |8.629       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1516_c   |Fixed 2,000 (23.9%)                |8.597      |0.105   |1.2%    |8.567      |8.626       |0.009                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1516_c   |Fixed 3,000 (35.9%)                |8.599      |0.063   |0.7%    |8.581      |8.617       |0.011                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1516_c   |Fixed 4,000 (47.9%)                |8.588      |0.057   |0.7%    |8.572      |8.604       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1604_c   |1% of pixels = 72 pixels (1.0%)    |7.268      |0.560   |7.7%    |7.109      |7.427       |0.020                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1604_c   |2% of pixels = 145 pixels (2.0%)   |7.267      |0.371   |5.1%    |7.161      |7.372       |0.019                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1604_c   |3% of pixels = 217 pixels (3.0%)   |7.283      |0.310   |4.3%    |7.195      |7.371       |0.035                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1604_c   |Fixed 1,250 (17.3%)                |7.261      |0.116   |1.6%    |7.228      |7.294       |0.013                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1604_c   |Fixed 2,000 (27.7%)                |7.245      |0.078   |1.1%    |7.223      |7.267       |-0.003                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1604_c   |Fixed 3,000 (41.5%)                |7.224      |0.051   |0.7%    |7.209      |7.238       |-0.024                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1604_c   |Fixed 4,000 (55.3%)                |7.248      |0.042   |0.6%    |7.236      |7.260       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1606_c   |1% of pixels = 76 pixels (1.0%)    |7.514      |0.377   |5.0%    |7.407      |7.621       |-0.007                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1606_c   |2% of pixels = 152 pixels (2.0%)   |7.579      |0.258   |3.4%    |7.506      |7.653       |0.058                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1606_c   |3% of pixels = 228 pixels (3.0%)   |7.572      |0.278   |3.7%    |7.493      |7.651       |0.051                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1606_c   |Fixed 1,250 (16.5%)                |7.539      |0.096   |1.3%    |7.512      |7.566       |0.018                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1606_c   |Fixed 2,000 (26.3%)                |7.521      |0.066   |0.9%    |7.502      |7.539       |-0.000                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1606_c   |Fixed 3,000 (39.5%)                |7.523      |0.043   |0.6%    |7.510      |7.535       |0.002                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1606_c   |Fixed 4,000 (52.7%)                |7.521      |0.035   |0.5%    |7.511      |7.531       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1701_b   |1% of pixels = 93 pixels (1.0%)    |7.107      |0.484   |6.8%    |6.970      |7.245       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1701_b   |2% of pixels = 186 pixels (2.0%)   |7.167      |0.363   |5.1%    |7.064      |7.270       |0.065                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1701_b   |3% of pixels = 279 pixels (3.0%)   |7.100      |0.290   |4.1%    |7.018      |7.183       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1701_b   |Fixed 1,250 (13.4%)                |7.117      |0.125   |1.8%    |7.081      |7.152       |0.015                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1701_b   |Fixed 2,000 (21.5%)                |7.112      |0.101   |1.4%    |7.083      |7.140       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1701_b   |Fixed 3,000 (32.3%)                |7.098      |0.073   |1.0%    |7.077      |7.119       |-0.003                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1701_b   |Fixed 4,000 (43.0%)                |7.101      |0.049   |0.7%    |7.087      |7.115       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1803_c   |1% of pixels = 53 pixels (1.0%)    |6.625      |0.488   |7.4%    |6.487      |6.764       |-0.164                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1803_c   |2% of pixels = 106 pixels (2.0%)   |6.755      |0.320   |4.7%    |6.664      |6.846       |-0.034                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1803_c   |3% of pixels = 159 pixels (3.0%)   |6.757      |0.239   |3.5%    |6.689      |6.825       |-0.032                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1803_c   |Fixed 1,250 (23.6%)                |6.781      |0.084   |1.2%    |6.757      |6.805       |-0.008                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1803_c   |Fixed 2,000 (37.8%)                |6.795      |0.059   |0.9%    |6.778      |6.812       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1803_c   |Fixed 3,000 (56.7%)                |6.795      |0.041   |0.6%    |6.783      |6.806       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1803_c   |Fixed 4,000 (75.5%)                |6.789      |0.025   |0.4%    |6.782      |6.796       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1816_c   |1% of pixels = 52 pixels (1.0%)    |6.982      |0.536   |7.7%    |6.830      |7.134       |-0.107                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1816_c   |2% of pixels = 104 pixels (2.0%)   |7.081      |0.316   |4.5%    |6.991      |7.171       |-0.008                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1816_c   |3% of pixels = 156 pixels (3.0%)   |7.041      |0.327   |4.6%    |6.948      |7.134       |-0.048                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1816_c   |Fixed 1,250 (24.1%)                |7.078      |0.093   |1.3%    |7.051      |7.104       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1816_c   |Fixed 2,000 (38.5%)                |7.085      |0.063   |0.9%    |7.067      |7.103       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1816_c   |Fixed 3,000 (57.8%)                |7.089      |0.044   |0.6%    |7.077      |7.102       |-0.000                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1816_c   |Fixed 4,000 (77.0%)                |7.089      |0.028   |0.4%    |7.081      |7.097       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1912_c   |1% of pixels = 55 pixels (1.0%)    |7.094      |0.515   |7.3%    |6.948      |7.241       |-0.110                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1912_c   |2% of pixels = 109 pixels (2.0%)   |7.188      |0.302   |4.2%    |7.102      |7.274       |-0.016                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1912_c   |3% of pixels = 164 pixels (3.0%)   |7.234      |0.252   |3.5%    |7.162      |7.305       |0.029                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1912_c   |Fixed 1,250 (22.9%)                |7.204      |0.106   |1.5%    |7.174      |7.234       |-0.000                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1912_c   |Fixed 2,000 (36.6%)                |7.204      |0.062   |0.9%    |7.187      |7.222       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1912_c   |Fixed 3,000 (54.9%)                |7.222      |0.046   |0.6%    |7.209      |7.235       |0.017                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1912_c   |Fixed 4,000 (73.2%)                |7.204      |0.029   |0.4%    |7.196      |7.212       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1917_c   |1% of pixels = 55 pixels (1.0%)    |7.001      |0.510   |7.3%    |6.857      |7.146       |-0.089                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1917_c   |2% of pixels = 109 pixels (2.0%)   |7.037      |0.387   |5.5%    |6.927      |7.147       |-0.053                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1917_c   |3% of pixels = 164 pixels (3.0%)   |7.072      |0.305   |4.3%    |6.985      |7.159       |-0.018                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1917_c   |Fixed 1,250 (22.9%)                |7.081      |0.098   |1.4%    |7.053      |7.109       |-0.009                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1917_c   |Fixed 2,000 (36.7%)                |7.097      |0.065   |0.9%    |7.079      |7.116       |0.007                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1917_c   |Fixed 3,000 (55.0%)                |7.088      |0.049   |0.7%    |7.074      |7.102       |-0.002                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |1917_c   |Fixed 4,000 (73.4%)                |7.090      |0.036   |0.5%    |7.080      |7.100       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |204_d    |1% of pixels = 71 pixels (1.0%)    |8.083      |0.595   |7.4%    |7.914      |8.252       |-0.120                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |204_d    |2% of pixels = 143 pixels (2.0%)   |8.192      |0.388   |4.7%    |8.082      |8.302       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |204_d    |3% of pixels = 214 pixels (3.0%)   |8.213      |0.344   |4.2%    |8.116      |8.311       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |204_d    |Fixed 1,250 (17.5%)                |8.199      |0.111   |1.4%    |8.167      |8.231       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |204_d    |Fixed 2,000 (28.0%)                |8.200      |0.089   |1.1%    |8.175      |8.225       |-0.003                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |204_d    |Fixed 3,000 (42.0%)                |8.191      |0.067   |0.8%    |8.172      |8.210       |-0.012                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |204_d    |Fixed 4,000 (56.0%)                |8.203      |0.043   |0.5%    |8.191      |8.216       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |217_d    |1% of pixels = 78 pixels (1.0%)    |7.570      |0.509   |6.7%    |7.425      |7.714       |-0.073                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |217_d    |2% of pixels = 157 pixels (2.0%)   |7.652      |0.278   |3.6%    |7.573      |7.731       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |217_d    |3% of pixels = 236 pixels (3.0%)   |7.582      |0.255   |3.4%    |7.509      |7.654       |-0.061                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |217_d    |Fixed 1,250 (15.9%)                |7.646      |0.101   |1.3%    |7.618      |7.675       |0.004                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |217_d    |Fixed 2,000 (25.5%)                |7.626      |0.082   |1.1%    |7.603      |7.649       |-0.017                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |217_d    |Fixed 3,000 (38.2%)                |7.638      |0.072   |0.9%    |7.618      |7.659       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |217_d    |Fixed 4,000 (51.0%)                |7.643      |0.051   |0.7%    |7.628      |7.657       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |22_c     |1% of pixels = 71 pixels (1.0%)    |8.808      |0.726   |8.2%    |8.601      |9.014       |0.064                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |22_c     |2% of pixels = 143 pixels (2.0%)   |8.658      |0.392   |4.5%    |8.546      |8.769       |-0.086                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |22_c     |3% of pixels = 214 pixels (3.0%)   |8.661      |0.412   |4.8%    |8.544      |8.778       |-0.083                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |22_c     |Fixed 1,250 (17.5%)                |8.733      |0.134   |1.5%    |8.695      |8.771       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |22_c     |Fixed 2,000 (28.0%)                |8.743      |0.099   |1.1%    |8.715      |8.771       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |22_c     |Fixed 3,000 (42.0%)                |8.725      |0.082   |0.9%    |8.702      |8.748       |-0.019                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |22_c     |Fixed 4,000 (56.1%)                |8.744      |0.058   |0.7%    |8.727      |8.760       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |319_c    |1% of pixels = 90 pixels (1.0%)    |8.182      |0.568   |6.9%    |8.020      |8.343       |-0.044                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |319_c    |2% of pixels = 180 pixels (2.0%)   |8.188      |0.416   |5.1%    |8.069      |8.306       |-0.038                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |319_c    |3% of pixels = 269 pixels (3.0%)   |8.247      |0.351   |4.3%    |8.147      |8.346       |0.021                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |319_c    |Fixed 1,250 (13.9%)                |8.231      |0.133   |1.6%    |8.194      |8.269       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |319_c    |Fixed 2,000 (22.3%)                |8.201      |0.111   |1.3%    |8.170      |8.232       |-0.025                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |319_c    |Fixed 3,000 (33.4%)                |8.211      |0.082   |1.0%    |8.188      |8.235       |-0.014                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |319_c    |Fixed 4,000 (44.5%)                |8.226      |0.057   |0.7%    |8.209      |8.242       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |409_d    |1% of pixels = 60 pixels (1.0%)    |8.933      |0.674   |7.5%    |8.742      |9.125       |-0.138                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |409_d    |2% of pixels = 120 pixels (2.0%)   |9.150      |0.485   |5.3%    |9.012      |9.287       |0.079                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |409_d    |3% of pixels = 180 pixels (3.0%)   |9.085      |0.367   |4.0%    |8.981      |9.189       |0.014                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |409_d    |Fixed 1,250 (20.9%)                |9.066      |0.141   |1.6%    |9.026      |9.106       |-0.005                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |409_d    |Fixed 2,000 (33.4%)                |9.095      |0.102   |1.1%    |9.066      |9.124       |0.024                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |409_d    |Fixed 3,000 (50.1%)                |9.085      |0.071   |0.8%    |9.065      |9.105       |0.015                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |409_d    |Fixed 4,000 (66.8%)                |9.071      |0.046   |0.5%    |9.058      |9.084       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |503_c    |1% of pixels = 74 pixels (1.0%)    |7.698      |0.484   |6.3%    |7.561      |7.836       |-0.100                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |503_c    |2% of pixels = 148 pixels (2.0%)   |7.693      |0.373   |4.8%    |7.587      |7.799       |-0.104                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |503_c    |3% of pixels = 222 pixels (3.0%)   |7.858      |0.262   |3.3%    |7.784      |7.933       |0.061                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |503_c    |Fixed 1,250 (16.9%)                |7.793      |0.117   |1.5%    |7.760      |7.827       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |503_c    |Fixed 2,000 (27.0%)                |7.786      |0.085   |1.1%    |7.762      |7.811       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |503_c    |Fixed 3,000 (40.5%)                |7.801      |0.058   |0.7%    |7.785      |7.818       |0.004                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |503_c    |Fixed 4,000 (54.0%)                |7.798      |0.045   |0.6%    |7.785      |7.810       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |614_a    |1% of pixels = 71 pixels (1.0%)    |6.124      |0.440   |7.2%    |5.999      |6.249       |-0.039                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |614_a    |2% of pixels = 142 pixels (2.0%)   |6.176      |0.321   |5.2%    |6.085      |6.267       |0.013                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |614_a    |3% of pixels = 213 pixels (3.0%)   |6.129      |0.270   |4.4%    |6.052      |6.206       |-0.034                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |614_a    |Fixed 1,250 (17.6%)                |6.181      |0.127   |2.1%    |6.145      |6.217       |0.018                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |614_a    |Fixed 2,000 (28.2%)                |6.168      |0.063   |1.0%    |6.150      |6.186       |0.005                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |614_a    |Fixed 3,000 (42.2%)                |6.162      |0.050   |0.8%    |6.148      |6.177       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |614_a    |Fixed 4,000 (56.3%)                |6.163      |0.036   |0.6%    |6.153      |6.173       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |700_c    |1% of pixels = 76 pixels (1.0%)    |6.937      |0.451   |6.5%    |6.808      |7.065       |-0.089                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |700_c    |2% of pixels = 152 pixels (2.0%)   |7.002      |0.263   |3.8%    |6.927      |7.077       |-0.024                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |700_c    |3% of pixels = 229 pixels (3.0%)   |7.056      |0.253   |3.6%    |6.984      |7.128       |0.030                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |700_c    |Fixed 1,250 (16.4%)                |6.986      |0.114   |1.6%    |6.954      |7.018       |-0.040                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |700_c    |Fixed 2,000 (26.2%)                |7.039      |0.079   |1.1%    |7.017      |7.062       |0.013                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |700_c    |Fixed 3,000 (39.4%)                |7.013      |0.062   |0.9%    |6.996      |7.031       |-0.013                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |700_c    |Fixed 4,000 (52.5%)                |7.026      |0.048   |0.7%    |7.012      |7.040       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |722_a    |1% of pixels = 88 pixels (1.0%)    |8.185      |0.521   |6.4%    |8.037      |8.333       |-0.019                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |722_a    |2% of pixels = 176 pixels (2.0%)   |8.143      |0.396   |4.9%    |8.030      |8.256       |-0.061                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |722_a    |3% of pixels = 263 pixels (3.0%)   |8.172      |0.309   |3.8%    |8.084      |8.259       |-0.032                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |722_a    |Fixed 1,250 (14.2%)                |8.208      |0.138   |1.7%    |8.169      |8.247       |0.004                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |722_a    |Fixed 2,000 (22.8%)                |8.194      |0.111   |1.4%    |8.162      |8.226       |-0.010                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |722_a    |Fixed 3,000 (34.2%)                |8.214      |0.082   |1.0%    |8.191      |8.237       |0.011                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |722_a    |Fixed 4,000 (45.6%)                |8.204      |0.052   |0.6%    |8.189      |8.218       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |800_a    |1% of pixels = 40 pixels (1.0%)    |9.366      |1.188   |12.7%   |9.028      |9.703       |-0.383                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |800_a    |2% of pixels = 80 pixels (2.0%)    |9.797      |1.050   |10.7%   |9.499      |10.095      |0.048                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |800_a    |3% of pixels = 119 pixels (3.0%)   |9.898      |0.740   |7.5%    |9.688      |10.109      |0.150                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |800_a    |Fixed 1,250 (31.4%)                |9.723      |0.189   |1.9%    |9.670      |9.777       |-0.026                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |800_a    |Fixed 2,000 (50.3%)                |9.772      |0.128   |1.3%    |9.736      |9.809       |0.023                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |800_a    |Fixed 3,000 (75.5%)                |9.763      |0.094   |1.0%    |9.737      |9.790       |0.014                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |800_a    |Fixed 4,000 (100.0%)               |9.749      |0.000   |0.0%    |9.749      |9.749       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |914_a    |1% of pixels = 72 pixels (1.0%)    |6.631      |0.504   |7.6%    |6.488      |6.775       |-0.073                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |914_a    |2% of pixels = 144 pixels (2.0%)   |6.678      |0.349   |5.2%    |6.579      |6.778       |-0.026                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |914_a    |3% of pixels = 217 pixels (3.0%)   |6.662      |0.274   |4.1%    |6.585      |6.740       |-0.041                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |914_a    |Fixed 1,250 (17.3%)                |6.709      |0.124   |1.8%    |6.674      |6.745       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |914_a    |Fixed 2,000 (27.7%)                |6.694      |0.082   |1.2%    |6.671      |6.717       |-0.010                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |914_a    |Fixed 3,000 (41.6%)                |6.699      |0.048   |0.7%    |6.686      |6.713       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |10m   |914_a    |Fixed 4,000 (55.4%)                |6.704      |0.040   |0.6%    |6.692      |6.715       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |100      |1% of pixels = 308 pixels (1.0%)   |8.805      |0.273   |3.1%    |8.727      |8.883       |0.035                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |100      |2% of pixels = 616 pixels (2.0%)   |8.769      |0.185   |2.1%    |8.717      |8.822       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |100      |3% of pixels = 924 pixels (3.0%)   |8.749      |0.134   |1.5%    |8.711      |8.787       |-0.021                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |100      |Fixed 1,250 (4.1%)                 |8.749      |0.124   |1.4%    |8.714      |8.784       |-0.021                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |100      |Fixed 2,000 (6.5%)                 |8.755      |0.100   |1.1%    |8.726      |8.783       |-0.015                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |100      |Fixed 3,000 (9.7%)                 |8.759      |0.080   |0.9%    |8.736      |8.781       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |100      |Fixed 4,000 (13.0%)                |8.770      |0.066   |0.7%    |8.751      |8.789       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1201     |1% of pixels = 248 pixels (1.0%)   |8.062      |0.254   |3.1%    |7.990      |8.134       |0.011                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1201     |2% of pixels = 496 pixels (2.0%)   |8.021      |0.187   |2.3%    |7.968      |8.074       |-0.030                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1201     |3% of pixels = 744 pixels (3.0%)   |8.066      |0.167   |2.1%    |8.018      |8.113       |0.015                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1201     |Fixed 1,250 (5.0%)                 |8.062      |0.112   |1.4%    |8.030      |8.094       |0.012                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1201     |Fixed 2,000 (8.1%)                 |8.068      |0.105   |1.3%    |8.038      |8.098       |0.017                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1201     |Fixed 3,000 (12.1%)                |8.076      |0.077   |1.0%    |8.054      |8.098       |0.025                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1201     |Fixed 4,000 (16.1%)                |8.051      |0.069   |0.9%    |8.031      |8.070       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1402     |1% of pixels = 271 pixels (1.0%)   |7.073      |0.255   |3.6%    |7.000      |7.146       |0.044                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1402     |2% of pixels = 541 pixels (2.0%)   |7.005      |0.154   |2.2%    |6.961      |7.049       |-0.024                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1402     |3% of pixels = 812 pixels (3.0%)   |7.045      |0.118   |1.7%    |7.012      |7.079       |0.016                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1402     |Fixed 1,250 (4.6%)                 |7.041      |0.104   |1.5%    |7.012      |7.071       |0.012                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1402     |Fixed 2,000 (7.4%)                 |7.032      |0.089   |1.3%    |7.007      |7.058       |0.003                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1402     |Fixed 3,000 (11.1%)                |7.050      |0.062   |0.9%    |7.032      |7.067       |0.021                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1402     |Fixed 4,000 (14.8%)                |7.029      |0.051   |0.7%    |7.014      |7.043       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1512     |1% of pixels = 310 pixels (1.0%)   |9.768      |0.346   |3.5%    |9.670      |9.867       |0.026                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1512     |2% of pixels = 620 pixels (2.0%)   |9.695      |0.185   |1.9%    |9.642      |9.747       |-0.048                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1512     |3% of pixels = 930 pixels (3.0%)   |9.772      |0.152   |1.6%    |9.729      |9.816       |0.030                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1512     |Fixed 1,250 (4.0%)                 |9.727      |0.107   |1.1%    |9.697      |9.758       |-0.015                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1512     |Fixed 2,000 (6.5%)                 |9.725      |0.095   |1.0%    |9.698      |9.752       |-0.018                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1512     |Fixed 3,000 (9.7%)                 |9.747      |0.101   |1.0%    |9.719      |9.776       |0.005                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1512     |Fixed 4,000 (12.9%)                |9.742      |0.088   |0.9%    |9.717      |9.767       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1518     |1% of pixels = 357 pixels (1.0%)   |9.819      |0.447   |4.6%    |9.692      |9.946       |0.052                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1518     |2% of pixels = 714 pixels (2.0%)   |9.828      |0.227   |2.3%    |9.763      |9.892       |0.061                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1518     |3% of pixels = 1,071 pixels (3.0%) |9.866      |0.265   |2.7%    |9.790      |9.941       |0.099                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1518     |Fixed 1,250 (3.5%)                 |9.771      |0.207   |2.1%    |9.712      |9.830       |0.004                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1518     |Fixed 2,000 (5.6%)                 |9.771      |0.178   |1.8%    |9.720      |9.821       |0.004                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1518     |Fixed 3,000 (8.4%)                 |9.786      |0.120   |1.2%    |9.752      |9.820       |0.019                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1518     |Fixed 4,000 (11.2%)                |9.767      |0.099   |1.0%    |9.739      |9.795       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1800     |1% of pixels = 300 pixels (1.0%)   |8.282      |0.278   |3.4%    |8.203      |8.361       |-0.049                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1800     |2% of pixels = 600 pixels (2.0%)   |8.339      |0.185   |2.2%    |8.287      |8.392       |0.008                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1800     |3% of pixels = 900 pixels (3.0%)   |8.345      |0.153   |1.8%    |8.302      |8.389       |0.014                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1800     |Fixed 1,250 (4.2%)                 |8.320      |0.132   |1.6%    |8.283      |8.358       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1800     |Fixed 2,000 (6.7%)                 |8.317      |0.091   |1.1%    |8.291      |8.343       |-0.014                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1800     |Fixed 3,000 (10.0%)                |8.329      |0.082   |1.0%    |8.306      |8.352       |-0.002                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1800     |Fixed 4,000 (13.3%)                |8.331      |0.072   |0.9%    |8.310      |8.352       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1805     |1% of pixels = 254 pixels (1.0%)   |8.223      |0.404   |4.9%    |8.108      |8.337       |-0.014                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1805     |2% of pixels = 508 pixels (2.0%)   |8.240      |0.210   |2.5%    |8.180      |8.299       |0.003                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1805     |3% of pixels = 763 pixels (3.0%)   |8.199      |0.176   |2.1%    |8.149      |8.249       |-0.038                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1805     |Fixed 1,250 (4.9%)                 |8.196      |0.136   |1.7%    |8.157      |8.234       |-0.041                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1805     |Fixed 2,000 (7.9%)                 |8.214      |0.101   |1.2%    |8.185      |8.242       |-0.023                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1805     |Fixed 3,000 (11.8%)                |8.219      |0.096   |1.2%    |8.192      |8.246       |-0.018                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1805     |Fixed 4,000 (15.7%)                |8.237      |0.066   |0.8%    |8.218      |8.256       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1910     |1% of pixels = 256 pixels (1.0%)   |6.477      |0.223   |3.4%    |6.414      |6.541       |-0.046                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1910     |2% of pixels = 512 pixels (2.0%)   |6.561      |0.146   |2.2%    |6.519      |6.602       |0.037                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1910     |3% of pixels = 768 pixels (3.0%)   |6.544      |0.116   |1.8%    |6.511      |6.577       |0.020                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1910     |Fixed 1,250 (4.9%)                 |6.522      |0.110   |1.7%    |6.491      |6.553       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1910     |Fixed 2,000 (7.8%)                 |6.516      |0.055   |0.8%    |6.501      |6.532       |-0.007                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1910     |Fixed 3,000 (11.7%)                |6.516      |0.072   |1.1%    |6.495      |6.536       |-0.008                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |1910     |Fixed 4,000 (15.6%)                |6.523      |0.055   |0.8%    |6.508      |6.539       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |206      |1% of pixels = 233 pixels (1.0%)   |8.195      |0.262   |3.2%    |8.120      |8.269       |-0.000                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |206      |2% of pixels = 467 pixels (2.0%)   |8.127      |0.258   |3.2%    |8.054      |8.200       |-0.068                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |206      |3% of pixels = 700 pixels (3.0%)   |8.187      |0.205   |2.5%    |8.129      |8.245       |-0.008                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |206      |Fixed 1,250 (5.4%)                 |8.170      |0.106   |1.3%    |8.140      |8.201       |-0.025                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |206      |Fixed 2,000 (8.6%)                 |8.192      |0.102   |1.2%    |8.163      |8.221       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |206      |Fixed 3,000 (12.8%)                |8.191      |0.091   |1.1%    |8.166      |8.217       |-0.004                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |206      |Fixed 4,000 (17.1%)                |8.195      |0.063   |0.8%    |8.178      |8.213       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |219      |1% of pixels = 294 pixels (1.0%)   |8.819      |0.351   |4.0%    |8.719      |8.918       |-0.067                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |219      |2% of pixels = 588 pixels (2.0%)   |8.852      |0.238   |2.7%    |8.785      |8.920       |-0.033                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |219      |3% of pixels = 881 pixels (3.0%)   |8.915      |0.162   |1.8%    |8.868      |8.961       |0.029                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |219      |Fixed 1,250 (4.3%)                 |8.903      |0.146   |1.6%    |8.861      |8.944       |0.017                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |219      |Fixed 2,000 (6.8%)                 |8.923      |0.130   |1.5%    |8.887      |8.960       |0.038                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |219      |Fixed 3,000 (10.2%)                |8.922      |0.104   |1.2%    |8.893      |8.952       |0.037                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |219      |Fixed 4,000 (13.6%)                |8.886      |0.063   |0.7%    |8.867      |8.904       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |302      |1% of pixels = 263 pixels (1.0%)   |6.995      |0.210   |3.0%    |6.936      |7.055       |-0.044                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |302      |2% of pixels = 527 pixels (2.0%)   |7.016      |0.186   |2.7%    |6.964      |7.069       |-0.023                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |302      |3% of pixels = 790 pixels (3.0%)   |7.002      |0.149   |2.1%    |6.960      |7.045       |-0.037                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |302      |Fixed 1,250 (4.7%)                 |7.013      |0.142   |2.0%    |6.973      |7.054       |-0.026                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |302      |Fixed 2,000 (7.6%)                 |7.032      |0.084   |1.2%    |7.008      |7.056       |-0.007                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |302      |Fixed 3,000 (11.4%)                |7.044      |0.065   |0.9%    |7.026      |7.063       |0.005                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |302      |Fixed 4,000 (15.2%)                |7.040      |0.053   |0.8%    |7.025      |7.055       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |317      |1% of pixels = 312 pixels (1.0%)   |7.705      |0.184   |2.4%    |7.653      |7.757       |-0.048                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |317      |2% of pixels = 623 pixels (2.0%)   |7.743      |0.175   |2.3%    |7.694      |7.793       |-0.009                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |317      |3% of pixels = 935 pixels (3.0%)   |7.727      |0.135   |1.8%    |7.689      |7.765       |-0.026                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |317      |Fixed 1,250 (4.0%)                 |7.760      |0.116   |1.5%    |7.727      |7.793       |0.007                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |317      |Fixed 2,000 (6.4%)                 |7.779      |0.092   |1.2%    |7.752      |7.805       |0.026                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |317      |Fixed 3,000 (9.6%)                 |7.750      |0.071   |0.9%    |7.730      |7.771       |-0.002                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |317      |Fixed 4,000 (12.8%)                |7.753      |0.056   |0.7%    |7.737      |7.769       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |405      |1% of pixels = 364 pixels (1.0%)   |8.366      |0.296   |3.5%    |8.282      |8.450       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |405      |2% of pixels = 727 pixels (2.0%)   |8.398      |0.208   |2.5%    |8.339      |8.457       |0.032                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |405      |3% of pixels = 1,091 pixels (3.0%) |8.342      |0.171   |2.1%    |8.294      |8.391       |-0.023                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |405      |Fixed 1,250 (3.4%)                 |8.371      |0.135   |1.6%    |8.333      |8.410       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |405      |Fixed 2,000 (5.5%)                 |8.374      |0.116   |1.4%    |8.341      |8.407       |0.008                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |405      |Fixed 3,000 (8.2%)                 |8.356      |0.077   |0.9%    |8.334      |8.378       |-0.010                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |405      |Fixed 4,000 (11.0%)                |8.366      |0.094   |1.1%    |8.339      |8.392       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |821      |1% of pixels = 292 pixels (1.0%)   |8.253      |0.236   |2.9%    |8.186      |8.320       |0.003                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |821      |2% of pixels = 584 pixels (2.0%)   |8.260      |0.184   |2.2%    |8.207      |8.312       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |821      |3% of pixels = 876 pixels (3.0%)   |8.270      |0.147   |1.8%    |8.228      |8.312       |0.020                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |821      |Fixed 1,250 (4.3%)                 |8.248      |0.124   |1.5%    |8.212      |8.283       |-0.002                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |821      |Fixed 2,000 (6.9%)                 |8.252      |0.097   |1.2%    |8.224      |8.280       |0.002                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |821      |Fixed 3,000 (10.3%)                |8.281      |0.089   |1.1%    |8.256      |8.306       |0.031                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |821      |Fixed 4,000 (13.7%)                |8.250      |0.069   |0.8%    |8.231      |8.270       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |905      |1% of pixels = 270 pixels (1.0%)   |10.880     |0.546   |5.0%    |10.724     |11.035      |-0.124                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |905      |2% of pixels = 540 pixels (2.0%)   |11.007     |0.561   |5.1%    |10.848     |11.167      |0.004                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |905      |3% of pixels = 811 pixels (3.0%)   |10.969     |0.426   |3.9%    |10.848     |11.090      |-0.035                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |905      |Fixed 1,250 (4.6%)                 |10.997     |0.343   |3.1%    |10.899     |11.094      |-0.007                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |905      |Fixed 2,000 (7.4%)                 |10.959     |0.272   |2.5%    |10.881     |11.036      |-0.045                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |905      |Fixed 3,000 (11.1%)                |10.949     |0.203   |1.9%    |10.891     |11.006      |-0.055                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |905      |Fixed 4,000 (14.8%)                |11.003     |0.182   |1.7%    |10.952     |11.055      |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |912      |1% of pixels = 296 pixels (1.0%)   |9.240      |0.257   |2.8%    |9.167      |9.313       |0.003                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |912      |2% of pixels = 592 pixels (2.0%)   |9.274      |0.181   |1.9%    |9.223      |9.326       |0.037                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |912      |3% of pixels = 889 pixels (3.0%)   |9.258      |0.115   |1.2%    |9.225      |9.290       |0.021                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |912      |Fixed 1,250 (4.2%)                 |9.223      |0.125   |1.4%    |9.187      |9.258       |-0.014                 |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |912      |Fixed 2,000 (6.8%)                 |9.238      |0.118   |1.3%    |9.204      |9.271       |0.001                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |912      |Fixed 3,000 (10.1%)                |9.255      |0.080   |0.9%    |9.233      |9.278       |0.018                  |all_sampled_pixels                                                        |
|PCA Mean Distance |20m   |912      |Fixed 4,000 (13.5%)                |9.237      |0.065   |0.7%    |9.218      |9.255       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_15 |Fixed 1,250 (0.7%)                 |7.369      |0.099   |1.3%    |7.341      |7.397       |-0.023                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_15 |1% of pixels = 1,728 pixels (1.0%) |7.384      |0.103   |1.4%    |7.355      |7.413       |-0.008                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_15 |2% of pixels = 3,456 pixels (2.0%) |7.396      |0.064   |0.9%    |7.378      |7.414       |0.004                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_15 |Fixed 4,000 (2.3%)                 |7.392      |0.062   |0.8%    |7.375      |7.410       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_15 |3% of pixels = 5,000 pixels (2.9%) |7.393      |0.057   |0.8%    |7.377      |7.410       |0.001                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_15 |Fixed 6,000 (3.5%)                 |7.401      |0.050   |0.7%    |7.386      |7.415       |0.009                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_15 |Fixed 8,000 (4.6%)                 |7.380      |0.043   |0.6%    |7.368      |7.392       |-0.012                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_24 |Fixed 1,250 (0.7%)                 |7.848      |0.131   |1.7%    |7.811      |7.885       |-0.021                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_24 |1% of pixels = 1,668 pixels (1.0%) |7.892      |0.118   |1.5%    |7.858      |7.925       |0.023                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_24 |2% of pixels = 3,336 pixels (2.0%) |7.864      |0.068   |0.9%    |7.845      |7.883       |-0.005                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_24 |Fixed 4,000 (2.4%)                 |7.868      |0.066   |0.8%    |7.850      |7.887       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_24 |3% of pixels = 5,000 pixels (3.0%) |7.872      |0.056   |0.7%    |7.855      |7.888       |0.003                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_24 |Fixed 6,000 (3.6%)                 |7.865      |0.055   |0.7%    |7.850      |7.881       |-0.003                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_24 |Fixed 8,000 (4.8%)                 |7.881      |0.048   |0.6%    |7.868      |7.895       |0.013                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_36 |Fixed 1,250 (0.7%)                 |9.173      |0.125   |1.4%    |9.137      |9.209       |-0.013                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_36 |1% of pixels = 1,742 pixels (1.0%) |9.224      |0.126   |1.4%    |9.188      |9.260       |0.038                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_36 |2% of pixels = 3,484 pixels (2.0%) |9.193      |0.064   |0.7%    |9.175      |9.212       |0.007                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_36 |Fixed 4,000 (2.3%)                 |9.186      |0.071   |0.8%    |9.166      |9.206       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_36 |3% of pixels = 5,000 pixels (2.9%) |9.204      |0.061   |0.7%    |9.187      |9.221       |0.018                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_36 |Fixed 6,000 (3.4%)                 |9.196      |0.065   |0.7%    |9.178      |9.215       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_36 |Fixed 8,000 (4.6%)                 |9.201      |0.058   |0.6%    |9.184      |9.217       |0.014                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_49 |Fixed 1,250 (0.6%)                 |8.619      |0.175   |2.0%    |8.570      |8.669       |-0.024                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_49 |1% of pixels = 2,070 pixels (1.0%) |8.642      |0.144   |1.7%    |8.601      |8.683       |-0.001                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_49 |Fixed 4,000 (1.9%)                 |8.644      |0.081   |0.9%    |8.621      |8.667       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_49 |2% of pixels = 4,139 pixels (2.0%) |8.625      |0.087   |1.0%    |8.600      |8.649       |-0.019                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_49 |3% of pixels = 5,000 pixels (2.4%) |8.618      |0.084   |1.0%    |8.594      |8.642       |-0.026                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_49 |Fixed 6,000 (2.9%)                 |8.637      |0.068   |0.8%    |8.618      |8.656       |-0.007                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_49 |Fixed 8,000 (3.9%)                 |8.633      |0.072   |0.8%    |8.612      |8.653       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_51 |Fixed 1,250 (0.6%)                 |8.583      |0.141   |1.6%    |8.543      |8.623       |-0.006                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_51 |1% of pixels = 2,006 pixels (1.0%) |8.587      |0.097   |1.1%    |8.559      |8.614       |-0.002                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_51 |Fixed 4,000 (2.0%)                 |8.589      |0.074   |0.9%    |8.568      |8.610       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_51 |2% of pixels = 4,013 pixels (2.0%) |8.567      |0.078   |0.9%    |8.545      |8.589       |-0.022                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_51 |3% of pixels = 5,000 pixels (2.5%) |8.556      |0.062   |0.7%    |8.538      |8.574       |-0.033                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_51 |Fixed 6,000 (3.0%)                 |8.584      |0.063   |0.7%    |8.566      |8.602       |-0.005                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_51 |Fixed 8,000 (4.0%)                 |8.577      |0.044   |0.5%    |8.564      |8.589       |-0.013                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_56 |Fixed 1,250 (0.9%)                 |8.165      |0.101   |1.2%    |8.136      |8.194       |-0.038                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_56 |1% of pixels = 1,455 pixels (1.0%) |8.234      |0.108   |1.3%    |8.203      |8.265       |0.031                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_56 |2% of pixels = 2,910 pixels (2.0%) |8.213      |0.078   |1.0%    |8.191      |8.235       |0.010                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_56 |Fixed 4,000 (2.7%)                 |8.203      |0.067   |0.8%    |8.184      |8.222       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_56 |3% of pixels = 4,365 pixels (3.0%) |8.212      |0.062   |0.8%    |8.194      |8.229       |0.009                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_56 |Fixed 6,000 (4.1%)                 |8.200      |0.060   |0.7%    |8.183      |8.217       |-0.003                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_56 |Fixed 8,000 (5.5%)                 |8.203      |0.045   |0.6%    |8.191      |8.216       |0.001                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_60 |Fixed 1,250 (0.8%)                 |9.490      |0.154   |1.6%    |9.446      |9.533       |0.012                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_60 |1% of pixels = 1,645 pixels (1.0%) |9.483      |0.130   |1.4%    |9.446      |9.520       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_60 |2% of pixels = 3,290 pixels (2.0%) |9.496      |0.088   |0.9%    |9.472      |9.521       |0.019                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_60 |Fixed 4,000 (2.4%)                 |9.477      |0.082   |0.9%    |9.454      |9.501       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_60 |3% of pixels = 4,935 pixels (3.0%) |9.483      |0.080   |0.8%    |9.460      |9.506       |0.006                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_60 |Fixed 6,000 (3.6%)                 |9.497      |0.065   |0.7%    |9.479      |9.516       |0.020                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_60 |Fixed 8,000 (4.9%)                 |9.490      |0.051   |0.5%    |9.475      |9.504       |0.013                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_8  |Fixed 1,250 (0.7%)                 |8.544      |0.145   |1.7%    |8.503      |8.585       |-0.021                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_8  |1% of pixels = 1,879 pixels (1.0%) |8.554      |0.112   |1.3%    |8.522      |8.586       |-0.011                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_8  |2% of pixels = 3,758 pixels (2.0%) |8.573      |0.088   |1.0%    |8.548      |8.598       |0.008                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_8  |Fixed 4,000 (2.1%)                 |8.565      |0.070   |0.8%    |8.545      |8.585       |0.000                  |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_8  |3% of pixels = 5,000 pixels (2.7%) |8.558      |0.049   |0.6%    |8.544      |8.572       |-0.007                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_8  |Fixed 6,000 (3.2%)                 |8.557      |0.057   |0.7%    |8.541      |8.573       |-0.008                 |all_sampled_pixels                                                        |
|PCA Mean Distance |50m   |sub50_8  |Fixed 8,000 (4.3%)                 |8.557      |0.045   |0.5%    |8.544      |8.569       |-0.008                 |all_sampled_pixels                                                        |
|Spectral Rao's Q  |10m   |1104_b   |1% of pixels = 94 pixels (1.0%)    |113.504    |25.665  |22.6%   |106.210    |120.798     |0.693                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1104_b   |2% of pixels = 188 pixels (2.0%)   |111.913    |15.939  |14.2%   |107.383    |116.442     |-0.899                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1104_b   |3% of pixels = 283 pixels (3.0%)   |111.702    |14.084  |12.6%   |107.699    |115.704     |-1.110                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1104_b   |Fixed 1,250 (13.3%)                |111.216    |5.868   |5.3%    |109.548    |112.884     |-1.595                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1104_b   |Fixed 2,000 (21.2%)                |112.732    |4.055   |3.6%    |111.579    |113.884     |-0.080                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1104_b   |Fixed 3,000 (31.8%)                |112.335    |3.118   |2.8%    |111.449    |113.221     |-0.476                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1104_b   |Fixed 4,000 (42.4%)                |112.811    |2.927   |2.6%    |111.979    |113.643     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1110_a   |1% of pixels = 76 pixels (1.0%)    |108.307    |12.048  |11.1%   |104.883    |111.731     |-1.144                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1110_a   |2% of pixels = 153 pixels (2.0%)   |108.239    |9.929   |9.2%    |105.417    |111.061     |-1.212                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1110_a   |3% of pixels = 230 pixels (3.0%)   |108.028    |5.276   |4.9%    |106.528    |109.527     |-1.423                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1110_a   |Fixed 1,250 (16.3%)                |109.653    |3.498   |3.2%    |108.659    |110.647     |0.202                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1110_a   |Fixed 2,000 (26.1%)                |109.492    |2.119   |1.9%    |108.890    |110.094     |0.041                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1110_a   |Fixed 3,000 (39.2%)                |109.290    |1.526   |1.4%    |108.856    |109.724     |-0.161                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1110_a   |Fixed 4,000 (52.3%)                |109.451    |1.140   |1.0%    |109.127    |109.775     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1124_a   |1% of pixels = 83 pixels (1.0%)    |134.905    |17.386  |12.9%   |129.964    |139.846     |-3.006                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1124_a   |2% of pixels = 166 pixels (2.0%)   |137.693    |12.671  |9.2%    |134.092    |141.294     |-0.217                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1124_a   |3% of pixels = 249 pixels (3.0%)   |137.244    |9.954   |7.3%    |134.415    |140.073     |-0.666                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1124_a   |Fixed 1,250 (15.1%)                |138.598    |3.782   |2.7%    |137.523    |139.673     |0.688                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1124_a   |Fixed 2,000 (24.1%)                |137.332    |3.232   |2.4%    |136.414    |138.251     |-0.578                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1124_a   |Fixed 3,000 (36.1%)                |139.170    |2.336   |1.7%    |138.506    |139.834     |1.260                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1124_a   |Fixed 4,000 (48.2%)                |137.910    |2.012   |1.5%    |137.338    |138.482     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |112_b    |1% of pixels = 73 pixels (1.0%)    |276.977    |104.418 |37.7%   |247.302    |306.652     |3.174                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |112_b    |2% of pixels = 147 pixels (2.0%)   |271.131    |70.663  |26.1%   |251.048    |291.213     |-2.672                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |112_b    |3% of pixels = 220 pixels (3.0%)   |279.960    |68.386  |24.4%   |260.525    |299.395     |6.157                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |112_b    |Fixed 1,250 (17.0%)                |272.933    |21.477  |7.9%    |266.829    |279.036     |-0.870                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |112_b    |Fixed 2,000 (27.2%)                |273.224    |19.529  |7.1%    |267.674    |278.774     |-0.579                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |112_b    |Fixed 3,000 (40.9%)                |269.266    |11.955  |4.4%    |265.868    |272.663     |-4.537                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |112_b    |Fixed 4,000 (54.5%)                |273.803    |10.830  |4.0%    |270.725    |276.881     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |114_b    |1% of pixels = 72 pixels (1.0%)    |233.599    |41.994  |18.0%   |221.665    |245.533     |-1.698                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |114_b    |2% of pixels = 145 pixels (2.0%)   |235.496    |27.327  |11.6%   |227.730    |243.262     |0.199                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |114_b    |3% of pixels = 217 pixels (3.0%)   |229.053    |18.058  |7.9%    |223.921    |234.185     |-6.244                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |114_b    |Fixed 1,250 (17.3%)                |235.375    |9.564   |4.1%    |232.657    |238.094     |0.078                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |114_b    |Fixed 2,000 (27.6%)                |233.211    |6.577   |2.8%    |231.341    |235.080     |-2.087                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |114_b    |Fixed 3,000 (41.4%)                |236.313    |4.870   |2.1%    |234.929    |237.698     |1.016                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |114_b    |Fixed 4,000 (55.3%)                |235.297    |4.001   |1.7%    |234.160    |236.434     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |11_c     |1% of pixels = 62 pixels (1.0%)    |124.951    |16.433  |13.2%   |120.280    |129.621     |-1.414                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |11_c     |2% of pixels = 124 pixels (2.0%)   |124.277    |11.751  |9.5%    |120.937    |127.616     |-2.087                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |11_c     |3% of pixels = 186 pixels (3.0%)   |128.397    |10.904  |8.5%    |125.298    |131.496     |2.033                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |11_c     |Fixed 1,250 (20.2%)                |125.219    |3.719   |3.0%    |124.162    |126.276     |-1.145                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |11_c     |Fixed 2,000 (32.3%)                |126.008    |2.494   |2.0%    |125.299    |126.717     |-0.356                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |11_c     |Fixed 3,000 (48.4%)                |125.823    |1.746   |1.4%    |125.327    |126.319     |-0.541                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |11_c     |Fixed 4,000 (64.5%)                |126.364    |1.391   |1.1%    |125.969    |126.760     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |120_d    |1% of pixels = 71 pixels (1.0%)    |202.670    |23.556  |11.6%   |195.976    |209.365     |-5.706                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |120_d    |2% of pixels = 141 pixels (2.0%)   |209.966    |17.587  |8.4%    |204.968    |214.964     |1.589                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |120_d    |3% of pixels = 212 pixels (3.0%)   |205.818    |17.265  |8.4%    |200.911    |210.724     |-2.559                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |120_d    |Fixed 1,250 (17.7%)                |208.917    |6.290   |3.0%    |207.130    |210.705     |0.541                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |120_d    |Fixed 2,000 (28.3%)                |208.160    |3.867   |1.9%    |207.061    |209.259     |-0.217                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |120_d    |Fixed 3,000 (42.4%)                |208.926    |3.103   |1.5%    |208.045    |209.808     |0.550                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |120_d    |Fixed 4,000 (56.6%)                |208.377    |2.364   |1.1%    |207.705    |209.049     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1305_c   |1% of pixels = 80 pixels (1.0%)    |132.578    |14.105  |10.6%   |128.569    |136.587     |-3.160                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1305_c   |2% of pixels = 161 pixels (2.0%)   |134.958    |11.412  |8.5%    |131.714    |138.201     |-0.780                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1305_c   |3% of pixels = 241 pixels (3.0%)   |135.636    |8.842   |6.5%    |133.123    |138.149     |-0.102                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1305_c   |Fixed 1,250 (15.6%)                |135.712    |3.991   |2.9%    |134.578    |136.847     |-0.026                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1305_c   |Fixed 2,000 (24.9%)                |135.859    |2.422   |1.8%    |135.170    |136.547     |0.121                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1305_c   |Fixed 3,000 (37.4%)                |135.580    |1.907   |1.4%    |135.038    |136.122     |-0.158                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1305_c   |Fixed 4,000 (49.8%)                |135.738    |1.400   |1.0%    |135.340    |136.136     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1315_b   |1% of pixels = 91 pixels (1.0%)    |170.327    |23.958  |14.1%   |163.518    |177.136     |-0.095                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1315_b   |2% of pixels = 182 pixels (2.0%)   |170.196    |16.422  |9.6%    |165.528    |174.863     |-0.227                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1315_b   |3% of pixels = 272 pixels (3.0%)   |167.167    |14.204  |8.5%    |163.131    |171.204     |-3.255                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1315_b   |Fixed 1,250 (13.8%)                |168.735    |5.949   |3.5%    |167.045    |170.426     |-1.687                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1315_b   |Fixed 2,000 (22.0%)                |169.711    |4.064   |2.4%    |168.556    |170.867     |-0.711                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1315_b   |Fixed 3,000 (33.1%)                |169.213    |3.437   |2.0%    |168.236    |170.189     |-1.210                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1315_b   |Fixed 4,000 (44.1%)                |170.422    |3.072   |1.8%    |169.549    |171.295     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1400_b   |1% of pixels = 68 pixels (1.0%)    |154.902    |18.987  |12.3%   |149.506    |160.298     |-3.735                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1400_b   |2% of pixels = 136 pixels (2.0%)   |157.350    |13.150  |8.4%    |153.612    |161.087     |-1.287                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1400_b   |3% of pixels = 204 pixels (3.0%)   |159.019    |9.864   |6.2%    |156.216    |161.823     |0.382                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1400_b   |Fixed 1,250 (18.4%)                |158.786    |4.053   |2.6%    |157.634    |159.938     |0.149                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1400_b   |Fixed 2,000 (29.4%)                |158.943    |2.598   |1.6%    |158.205    |159.681     |0.306                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1400_b   |Fixed 3,000 (44.1%)                |158.187    |2.038   |1.3%    |157.608    |158.767     |-0.450                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1400_b   |Fixed 4,000 (58.8%)                |158.637    |1.316   |0.8%    |158.263    |159.011     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1414_a   |1% of pixels = 62 pixels (1.0%)    |205.081    |28.988  |14.1%   |196.843    |213.320     |4.730                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1414_a   |2% of pixels = 125 pixels (2.0%)   |194.894    |20.139  |10.3%   |189.171    |200.617     |-5.458                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1414_a   |3% of pixels = 187 pixels (3.0%)   |199.455    |19.153  |9.6%    |194.012    |204.899     |-0.896                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1414_a   |Fixed 1,250 (20.1%)                |199.500    |7.608   |3.8%    |197.338    |201.663     |-0.851                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1414_a   |Fixed 2,000 (32.1%)                |199.369    |3.945   |2.0%    |198.248    |200.490     |-0.982                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1414_a   |Fixed 3,000 (48.2%)                |200.910    |3.011   |1.5%    |200.055    |201.766     |0.558                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1414_a   |Fixed 4,000 (64.2%)                |200.352    |2.287   |1.1%    |199.702    |201.002     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1417_c   |1% of pixels = 75 pixels (1.0%)    |218.388    |30.482  |14.0%   |209.725    |227.051     |-5.757                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1417_c   |2% of pixels = 150 pixels (2.0%)   |226.158    |22.780  |10.1%   |219.684    |232.632     |2.013                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1417_c   |3% of pixels = 225 pixels (3.0%)   |222.819    |17.493  |7.9%    |217.847    |227.790     |-1.326                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1417_c   |Fixed 1,250 (16.6%)                |223.231    |7.285   |3.3%    |221.160    |225.301     |-0.914                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1417_c   |Fixed 2,000 (26.6%)                |223.667    |5.599   |2.5%    |222.076    |225.259     |-0.478                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1417_c   |Fixed 3,000 (39.9%)                |223.078    |3.607   |1.6%    |222.053    |224.103     |-1.067                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1417_c   |Fixed 4,000 (53.2%)                |224.145    |2.830   |1.3%    |223.341    |224.949     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1514_a   |1% of pixels = 70 pixels (1.0%)    |92.449     |14.455  |15.6%   |88.341     |96.557      |-1.567                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1514_a   |2% of pixels = 141 pixels (2.0%)   |93.922     |9.050   |9.6%    |91.350     |96.494      |-0.094                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1514_a   |3% of pixels = 211 pixels (3.0%)   |93.558     |8.777   |9.4%    |91.063     |96.052      |-0.457                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1514_a   |Fixed 1,250 (17.8%)                |94.572     |3.217   |3.4%    |93.657     |95.486      |0.557                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1514_a   |Fixed 2,000 (28.5%)                |94.181     |2.253   |2.4%    |93.540     |94.821      |0.165                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1514_a   |Fixed 3,000 (42.7%)                |93.740     |1.485   |1.6%    |93.318     |94.162      |-0.275                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1514_a   |Fixed 4,000 (56.9%)                |94.015     |1.186   |1.3%    |93.678     |94.352      |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1516_c   |1% of pixels = 84 pixels (1.0%)    |183.818    |25.251  |13.7%   |176.642    |190.994     |-0.162                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1516_c   |2% of pixels = 167 pixels (2.0%)   |181.148    |14.352  |7.9%    |177.069    |185.226     |-2.833                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1516_c   |3% of pixels = 251 pixels (3.0%)   |185.970    |11.408  |6.1%    |182.728    |189.212     |1.990                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1516_c   |Fixed 1,250 (15.0%)                |184.348    |5.051   |2.7%    |182.912    |185.784     |0.368                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1516_c   |Fixed 2,000 (23.9%)                |184.222    |4.319   |2.3%    |182.995    |185.450     |0.242                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1516_c   |Fixed 3,000 (35.9%)                |184.511    |2.831   |1.5%    |183.706    |185.315     |0.530                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1516_c   |Fixed 4,000 (47.9%)                |183.980    |2.408   |1.3%    |183.296    |184.665     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1604_c   |1% of pixels = 72 pixels (1.0%)    |139.081    |20.660  |14.9%   |133.209    |144.952     |2.492                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1604_c   |2% of pixels = 145 pixels (2.0%)   |137.763    |12.179  |8.8%    |134.302    |141.225     |1.175                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1604_c   |3% of pixels = 217 pixels (3.0%)   |137.847    |10.751  |7.8%    |134.791    |140.902     |1.258                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1604_c   |Fixed 1,250 (17.3%)                |136.869    |4.118   |3.0%    |135.699    |138.039     |0.281                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1604_c   |Fixed 2,000 (27.7%)                |136.380    |2.694   |2.0%    |135.614    |137.145     |-0.208                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1604_c   |Fixed 3,000 (41.5%)                |135.807    |1.864   |1.4%    |135.277    |136.337     |-0.781                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1604_c   |Fixed 4,000 (55.3%)                |136.588    |1.588   |1.2%    |136.137    |137.040     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1606_c   |1% of pixels = 76 pixels (1.0%)    |135.352    |13.671  |10.1%   |131.467    |139.237     |-0.731                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1606_c   |2% of pixels = 152 pixels (2.0%)   |138.542    |9.549   |6.9%    |135.828    |141.256     |2.459                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1606_c   |3% of pixels = 228 pixels (3.0%)   |137.601    |9.689   |7.0%    |134.848    |140.355     |1.518                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1606_c   |Fixed 1,250 (16.5%)                |136.854    |3.430   |2.5%    |135.879    |137.829     |0.771                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1606_c   |Fixed 2,000 (26.3%)                |136.079    |2.229   |1.6%    |135.446    |136.713     |-0.004                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1606_c   |Fixed 3,000 (39.5%)                |136.066    |1.468   |1.1%    |135.649    |136.484     |-0.017                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1606_c   |Fixed 4,000 (52.7%)                |136.083    |1.190   |0.9%    |135.745    |136.421     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1701_b   |1% of pixels = 93 pixels (1.0%)    |137.957    |18.857  |13.7%   |132.598    |143.316     |1.900                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1701_b   |2% of pixels = 186 pixels (2.0%)   |138.201    |13.776  |10.0%   |134.286    |142.116     |2.144                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1701_b   |3% of pixels = 279 pixels (3.0%)   |135.676    |10.891  |8.0%    |132.581    |138.772     |-0.380                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1701_b   |Fixed 1,250 (13.4%)                |136.704    |5.121   |3.7%    |135.248    |138.159     |0.647                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1701_b   |Fixed 2,000 (21.5%)                |136.426    |3.867   |2.8%    |135.327    |137.525     |0.370                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1701_b   |Fixed 3,000 (32.3%)                |136.180    |2.922   |2.1%    |135.350    |137.011     |0.124                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1701_b   |Fixed 4,000 (43.0%)                |136.056    |1.978   |1.5%    |135.494    |136.618     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1803_c   |1% of pixels = 53 pixels (1.0%)    |107.686    |13.719  |12.7%   |103.787    |111.585     |-5.711                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1803_c   |2% of pixels = 106 pixels (2.0%)   |111.727    |10.526  |9.4%    |108.735    |114.718     |-1.671                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1803_c   |3% of pixels = 159 pixels (3.0%)   |112.767    |7.410   |6.6%    |110.661    |114.872     |-0.631                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1803_c   |Fixed 1,250 (23.6%)                |113.030    |2.663   |2.4%    |112.273    |113.786     |-0.368                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1803_c   |Fixed 2,000 (37.8%)                |113.569    |1.886   |1.7%    |113.033    |114.105     |0.172                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1803_c   |Fixed 3,000 (56.7%)                |113.535    |1.198   |1.1%    |113.194    |113.875     |0.138                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1803_c   |Fixed 4,000 (75.5%)                |113.397    |0.765   |0.7%    |113.180    |113.615     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1816_c   |1% of pixels = 52 pixels (1.0%)    |124.296    |19.203  |15.4%   |118.838    |129.753     |-3.871                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1816_c   |2% of pixels = 104 pixels (2.0%)   |128.657    |11.282  |8.8%    |125.451    |131.864     |0.491                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1816_c   |3% of pixels = 156 pixels (3.0%)   |126.184    |11.600  |9.2%    |122.887    |129.481     |-1.982                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1816_c   |Fixed 1,250 (24.1%)                |127.859    |3.372   |2.6%    |126.901    |128.817     |-0.307                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1816_c   |Fixed 2,000 (38.5%)                |127.795    |2.241   |1.8%    |127.158    |128.432     |-0.372                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1816_c   |Fixed 3,000 (57.8%)                |128.040    |1.547   |1.2%    |127.600    |128.480     |-0.127                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1816_c   |Fixed 4,000 (77.0%)                |128.167    |0.952   |0.7%    |127.896    |128.437     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1912_c   |1% of pixels = 55 pixels (1.0%)    |128.584    |18.637  |14.5%   |123.287    |133.880     |-2.443                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1912_c   |2% of pixels = 109 pixels (2.0%)   |130.320    |10.739  |8.2%    |127.268    |133.372     |-0.707                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1912_c   |3% of pixels = 164 pixels (3.0%)   |131.868    |8.394   |6.4%    |129.483    |134.254     |0.841                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1912_c   |Fixed 1,250 (22.9%)                |130.866    |3.460   |2.6%    |129.883    |131.849     |-0.161                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1912_c   |Fixed 2,000 (36.6%)                |131.144    |2.282   |1.7%    |130.496    |131.792     |0.117                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1912_c   |Fixed 3,000 (54.9%)                |131.627    |1.664   |1.3%    |131.154    |132.099     |0.600                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1912_c   |Fixed 4,000 (73.2%)                |131.027    |0.984   |0.8%    |130.747    |131.307     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1917_c   |1% of pixels = 55 pixels (1.0%)    |125.694    |19.991  |15.9%   |120.013    |131.376     |-2.257                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1917_c   |2% of pixels = 109 pixels (2.0%)   |125.720    |13.997  |11.1%   |121.742    |129.698     |-2.232                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1917_c   |3% of pixels = 164 pixels (3.0%)   |127.885    |10.988  |8.6%    |124.762    |131.008     |-0.067                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1917_c   |Fixed 1,250 (22.9%)                |127.600    |3.656   |2.9%    |126.561    |128.639     |-0.352                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1917_c   |Fixed 2,000 (36.7%)                |128.175    |2.286   |1.8%    |127.525    |128.824     |0.223                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1917_c   |Fixed 3,000 (55.0%)                |127.809    |1.789   |1.4%    |127.300    |128.317     |-0.143                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |1917_c   |Fixed 4,000 (73.4%)                |127.952    |1.246   |1.0%    |127.598    |128.306     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |204_d    |1% of pixels = 71 pixels (1.0%)    |167.665    |23.492  |14.0%   |160.989    |174.341     |-4.692                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |204_d    |2% of pixels = 143 pixels (2.0%)   |172.574    |14.731  |8.5%    |168.388    |176.761     |0.217                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |204_d    |3% of pixels = 214 pixels (3.0%)   |172.195    |13.328  |7.7%    |168.407    |175.983     |-0.162                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |204_d    |Fixed 1,250 (17.5%)                |172.511    |4.593   |2.7%    |171.206    |173.817     |0.154                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |204_d    |Fixed 2,000 (28.0%)                |172.075    |3.424   |2.0%    |171.102    |173.048     |-0.282                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |204_d    |Fixed 3,000 (42.0%)                |171.913    |2.485   |1.4%    |171.207    |172.619     |-0.444                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |204_d    |Fixed 4,000 (56.0%)                |172.357    |1.754   |1.0%    |171.859    |172.856     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |217_d    |1% of pixels = 78 pixels (1.0%)    |150.000    |20.154  |13.4%   |144.273    |155.728     |-3.536                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |217_d    |2% of pixels = 157 pixels (2.0%)   |153.261    |11.506  |7.5%    |149.991    |156.531     |-0.276                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |217_d    |3% of pixels = 236 pixels (3.0%)   |150.629    |10.354  |6.9%    |147.686    |153.571     |-2.908                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |217_d    |Fixed 1,250 (15.9%)                |153.256    |4.122   |2.7%    |152.085    |154.428     |-0.280                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |217_d    |Fixed 2,000 (25.5%)                |153.006    |3.169   |2.1%    |152.105    |153.906     |-0.531                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |217_d    |Fixed 3,000 (38.2%)                |153.205    |2.634   |1.7%    |152.457    |153.953     |-0.332                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |217_d    |Fixed 4,000 (51.0%)                |153.537    |1.904   |1.2%    |152.996    |154.078     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |22_c     |1% of pixels = 71 pixels (1.0%)    |208.563    |31.923  |15.3%   |199.490    |217.635     |3.331                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |22_c     |2% of pixels = 143 pixels (2.0%)   |202.302    |17.645  |8.7%    |197.287    |207.317     |-2.930                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |22_c     |3% of pixels = 214 pixels (3.0%)   |201.507    |17.997  |8.9%    |196.392    |206.622     |-3.725                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |22_c     |Fixed 1,250 (17.5%)                |204.529    |5.634   |2.8%    |202.927    |206.130     |-0.703                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |22_c     |Fixed 2,000 (28.0%)                |204.785    |4.478   |2.2%    |203.513    |206.058     |-0.447                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |22_c     |Fixed 3,000 (42.0%)                |203.956    |3.687   |1.8%    |202.909    |205.004     |-1.275                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |22_c     |Fixed 4,000 (56.1%)                |205.232    |2.507   |1.2%    |204.519    |205.944     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |319_c    |1% of pixels = 90 pixels (1.0%)    |177.581    |23.395  |13.2%   |170.933    |184.230     |-1.212                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |319_c    |2% of pixels = 180 pixels (2.0%)   |179.486    |17.371  |9.7%    |174.549    |184.423     |0.693                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |319_c    |3% of pixels = 269 pixels (3.0%)   |179.195    |13.870  |7.7%    |175.253    |183.136     |0.402                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |319_c    |Fixed 1,250 (13.9%)                |179.080    |5.405   |3.0%    |177.544    |180.616     |0.287                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |319_c    |Fixed 2,000 (22.3%)                |177.978    |4.494   |2.5%    |176.701    |179.255     |-0.815                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |319_c    |Fixed 3,000 (33.4%)                |178.018    |3.578   |2.0%    |177.001    |179.035     |-0.775                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |319_c    |Fixed 4,000 (44.5%)                |178.793    |2.290   |1.3%    |178.142    |179.444     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |409_d    |1% of pixels = 60 pixels (1.0%)    |213.237    |32.702  |15.3%   |203.943    |222.531     |-2.930                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |409_d    |2% of pixels = 120 pixels (2.0%)   |219.350    |23.248  |10.6%   |212.743    |225.957     |3.184                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |409_d    |3% of pixels = 180 pixels (3.0%)   |216.541    |19.881  |9.2%    |210.891    |222.191     |0.374                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |409_d    |Fixed 1,250 (20.9%)                |215.018    |6.591   |3.1%    |213.145    |216.891     |-1.148                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |409_d    |Fixed 2,000 (33.4%)                |217.095    |5.034   |2.3%    |215.664    |218.525     |0.928                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |409_d    |Fixed 3,000 (50.1%)                |216.832    |3.643   |1.7%    |215.797    |217.867     |0.665                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |409_d    |Fixed 4,000 (66.8%)                |216.167    |2.357   |1.1%    |215.497    |216.836     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |503_c    |1% of pixels = 74 pixels (1.0%)    |153.525    |17.371  |11.3%   |148.588    |158.461     |-2.680                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |503_c    |2% of pixels = 148 pixels (2.0%)   |152.245    |12.501  |8.2%    |148.692    |155.798     |-3.959                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |503_c    |3% of pixels = 222 pixels (3.0%)   |157.879    |9.358   |5.9%    |155.220    |160.538     |1.675                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |503_c    |Fixed 1,250 (16.9%)                |155.856    |4.141   |2.7%    |154.679    |157.033     |-0.349                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |503_c    |Fixed 2,000 (27.0%)                |155.616    |3.233   |2.1%    |154.697    |156.535     |-0.589                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |503_c    |Fixed 3,000 (40.5%)                |156.134    |2.052   |1.3%    |155.551    |156.717     |-0.070                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |503_c    |Fixed 4,000 (54.0%)                |156.204    |1.617   |1.0%    |155.745    |156.664     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |614_a    |1% of pixels = 71 pixels (1.0%)    |101.786    |14.731  |14.5%   |97.599     |105.972     |-0.949                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |614_a    |2% of pixels = 142 pixels (2.0%)   |102.977    |10.522  |10.2%   |99.986     |105.967     |0.242                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |614_a    |3% of pixels = 213 pixels (3.0%)   |101.038    |9.197   |9.1%    |98.425     |103.652     |-1.697                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |614_a    |Fixed 1,250 (17.6%)                |103.344    |4.249   |4.1%    |102.137    |104.552     |0.609                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |614_a    |Fixed 2,000 (28.2%)                |102.865    |2.187   |2.1%    |102.243    |103.486     |0.130                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |614_a    |Fixed 3,000 (42.2%)                |102.784    |1.405   |1.4%    |102.384    |103.183     |0.049                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |614_a    |Fixed 4,000 (56.3%)                |102.735    |1.212   |1.2%    |102.391    |103.079     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |700_c    |1% of pixels = 76 pixels (1.0%)    |127.717    |15.817  |12.4%   |123.222    |132.212     |-2.255                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |700_c    |2% of pixels = 152 pixels (2.0%)   |129.318    |9.115   |7.0%    |126.727    |131.908     |-0.655                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |700_c    |3% of pixels = 229 pixels (3.0%)   |130.865    |9.040   |6.9%    |128.296    |133.434     |0.892                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |700_c    |Fixed 1,250 (16.4%)                |128.774    |4.052   |3.1%    |127.623    |129.926     |-1.198                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |700_c    |Fixed 2,000 (26.2%)                |130.518    |2.888   |2.2%    |129.697    |131.339     |0.546                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |700_c    |Fixed 3,000 (39.4%)                |129.575    |2.228   |1.7%    |128.942    |130.209     |-0.397                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |700_c    |Fixed 4,000 (52.5%)                |129.972    |1.764   |1.4%    |129.471    |130.474     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |722_a    |1% of pixels = 88 pixels (1.0%)    |178.066    |19.659  |11.0%   |172.479    |183.653     |-1.396                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |722_a    |2% of pixels = 176 pixels (2.0%)   |177.474    |16.046  |9.0%    |172.913    |182.034     |-1.989                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |722_a    |3% of pixels = 263 pixels (3.0%)   |177.803    |12.210  |6.9%    |174.333    |181.273     |-1.659                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |722_a    |Fixed 1,250 (14.2%)                |180.002    |5.421   |3.0%    |178.461    |181.543     |0.540                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |722_a    |Fixed 2,000 (22.8%)                |179.164    |4.305   |2.4%    |177.940    |180.387     |-0.299                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |722_a    |Fixed 3,000 (34.2%)                |179.906    |3.245   |1.8%    |178.984    |180.828     |0.443                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |722_a    |Fixed 4,000 (45.6%)                |179.462    |2.101   |1.2%    |178.865    |180.059     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |800_a    |1% of pixels = 40 pixels (1.0%)    |270.011    |130.438 |48.3%   |232.941    |307.081     |-28.353                |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |800_a    |2% of pixels = 80 pixels (2.0%)    |294.542    |107.537 |36.5%   |263.981    |325.104     |-3.822                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |800_a    |3% of pixels = 119 pixels (3.0%)   |319.472    |76.980  |24.1%   |297.594    |341.349     |21.107                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |800_a    |Fixed 1,250 (31.4%)                |294.273    |21.320  |7.2%    |288.214    |300.333     |-4.091                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |800_a    |Fixed 2,000 (50.3%)                |300.877    |12.688  |4.2%    |297.271    |304.483     |2.513                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |800_a    |Fixed 3,000 (75.5%)                |299.033    |9.995   |3.3%    |296.193    |301.874     |0.669                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |800_a    |Fixed 4,000 (100.0%)               |298.364    |0.000   |0.0%    |298.364    |298.364     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |914_a    |1% of pixels = 72 pixels (1.0%)    |122.314    |18.237  |14.9%   |117.131    |127.496     |-1.789                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |914_a    |2% of pixels = 144 pixels (2.0%)   |124.833    |13.623  |10.9%   |120.961    |128.704     |0.730                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |914_a    |3% of pixels = 217 pixels (3.0%)   |122.139    |11.732  |9.6%    |118.805    |125.474     |-1.963                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |914_a    |Fixed 1,250 (17.3%)                |124.778    |4.783   |3.8%    |123.418    |126.137     |0.675                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |914_a    |Fixed 2,000 (27.7%)                |123.699    |3.384   |2.7%    |122.737    |124.661     |-0.404                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |914_a    |Fixed 3,000 (41.6%)                |124.020    |2.175   |1.8%    |123.402    |124.639     |-0.083                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |10m   |914_a    |Fixed 4,000 (55.4%)                |124.103    |1.639   |1.3%    |123.637    |124.569     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |100      |1% of pixels = 308 pixels (1.0%)   |194.353    |11.386  |5.9%    |191.117    |197.588     |1.429                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |100      |2% of pixels = 616 pixels (2.0%)   |193.310    |8.193   |4.2%    |190.982    |195.639     |0.387                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |100      |3% of pixels = 924 pixels (3.0%)   |192.411    |6.176   |3.2%    |190.656    |194.167     |-0.512                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |100      |Fixed 1,250 (4.1%)                 |192.246    |5.320   |2.8%    |190.734    |193.758     |-0.678                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |100      |Fixed 2,000 (6.5%)                 |192.399    |4.417   |2.3%    |191.144    |193.654     |-0.525                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |100      |Fixed 3,000 (9.7%)                 |192.708    |3.297   |1.7%    |191.771    |193.645     |-0.216                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |100      |Fixed 4,000 (13.0%)                |192.924    |2.731   |1.4%    |192.148    |193.700     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1201     |1% of pixels = 248 pixels (1.0%)   |165.880    |10.635  |6.4%    |162.858    |168.903     |0.909                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1201     |2% of pixels = 496 pixels (2.0%)   |163.969    |7.579   |4.6%    |161.815    |166.122     |-1.003                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1201     |3% of pixels = 744 pixels (3.0%)   |165.077    |7.308   |4.4%    |163.000    |167.154     |0.106                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1201     |Fixed 1,250 (5.0%)                 |165.656    |4.602   |2.8%    |164.348    |166.964     |0.685                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1201     |Fixed 2,000 (8.1%)                 |165.648    |4.195   |2.5%    |164.456    |166.840     |0.677                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1201     |Fixed 3,000 (12.1%)                |165.690    |3.241   |2.0%    |164.769    |166.612     |0.719                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1201     |Fixed 4,000 (16.1%)                |164.971    |2.599   |1.6%    |164.233    |165.710     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1402     |1% of pixels = 271 pixels (1.0%)   |127.164    |10.172  |8.0%    |124.273    |130.055     |1.352                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1402     |2% of pixels = 541 pixels (2.0%)   |125.274    |5.820   |4.6%    |123.620    |126.929     |-0.538                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1402     |3% of pixels = 812 pixels (3.0%)   |126.127    |4.629   |3.7%    |124.812    |127.443     |0.315                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1402     |Fixed 1,250 (4.6%)                 |126.122    |3.894   |3.1%    |125.015    |127.228     |0.309                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1402     |Fixed 2,000 (7.4%)                 |125.889    |3.291   |2.6%    |124.953    |126.824     |0.076                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1402     |Fixed 3,000 (11.1%)                |126.688    |2.428   |1.9%    |125.998    |127.378     |0.876                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1402     |Fixed 4,000 (14.8%)                |125.812    |2.044   |1.6%    |125.231    |126.393     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1512     |1% of pixels = 310 pixels (1.0%)   |249.149    |17.066  |6.8%    |244.299    |253.999     |0.706                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1512     |2% of pixels = 620 pixels (2.0%)   |245.803    |8.986   |3.7%    |243.249    |248.356     |-2.640                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1512     |3% of pixels = 930 pixels (3.0%)   |249.798    |7.138   |2.9%    |247.770    |251.827     |1.356                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1512     |Fixed 1,250 (4.0%)                 |247.277    |5.347   |2.2%    |245.757    |248.796     |-1.166                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1512     |Fixed 2,000 (6.5%)                 |247.229    |4.385   |1.8%    |245.983    |248.475     |-1.214                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1512     |Fixed 3,000 (9.7%)                 |248.362    |4.881   |2.0%    |246.975    |249.749     |-0.081                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1512     |Fixed 4,000 (12.9%)                |248.443    |4.063   |1.6%    |247.288    |249.597     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1518     |1% of pixels = 357 pixels (1.0%)   |280.889    |34.272  |12.2%   |271.149    |290.629     |3.048                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1518     |2% of pixels = 714 pixels (2.0%)   |280.010    |15.853  |5.7%    |275.504    |284.515     |2.168                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1518     |3% of pixels = 1,071 pixels (3.0%) |283.826    |20.249  |7.1%    |278.072    |289.581     |5.985                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1518     |Fixed 1,250 (3.5%)                 |276.962    |15.201  |5.5%    |272.642    |281.282     |-0.879                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1518     |Fixed 2,000 (5.6%)                 |277.656    |12.171  |4.4%    |274.197    |281.115     |-0.186                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1518     |Fixed 3,000 (8.4%)                 |278.416    |7.580   |2.7%    |276.262    |280.570     |0.575                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1518     |Fixed 4,000 (11.2%)                |277.841    |7.194   |2.6%    |275.797    |279.886     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1800     |1% of pixels = 300 pixels (1.0%)   |174.942    |11.896  |6.8%    |171.561    |178.323     |-1.708                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1800     |2% of pixels = 600 pixels (2.0%)   |176.455    |8.338   |4.7%    |174.085    |178.824     |-0.195                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1800     |3% of pixels = 900 pixels (3.0%)   |177.469    |6.220   |3.5%    |175.701    |179.236     |0.819                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1800     |Fixed 1,250 (4.2%)                 |176.406    |5.581   |3.2%    |174.820    |177.992     |-0.244                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1800     |Fixed 2,000 (6.7%)                 |176.464    |4.366   |2.5%    |175.224    |177.705     |-0.185                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1800     |Fixed 3,000 (10.0%)                |176.616    |3.652   |2.1%    |175.578    |177.654     |-0.034                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1800     |Fixed 4,000 (13.3%)                |176.650    |2.861   |1.6%    |175.837    |177.463     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1805     |1% of pixels = 254 pixels (1.0%)   |181.827    |22.940  |12.6%   |175.307    |188.346     |1.668                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1805     |2% of pixels = 508 pixels (2.0%)   |181.741    |11.872  |6.5%    |178.368    |185.115     |1.583                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1805     |3% of pixels = 763 pixels (3.0%)   |178.670    |10.345  |5.8%    |175.730    |181.610     |-1.488                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1805     |Fixed 1,250 (4.9%)                 |178.210    |6.402   |3.6%    |176.391    |180.029     |-1.948                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1805     |Fixed 2,000 (7.9%)                 |178.828    |5.322   |3.0%    |177.315    |180.340     |-1.331                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1805     |Fixed 3,000 (11.8%)                |179.288    |5.172   |2.9%    |177.818    |180.758     |-0.870                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1805     |Fixed 4,000 (15.7%)                |180.158    |3.467   |1.9%    |179.173    |181.144     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1910     |1% of pixels = 256 pixels (1.0%)   |110.618    |7.724   |7.0%    |108.423    |112.813     |-1.241                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1910     |2% of pixels = 512 pixels (2.0%)   |113.122    |4.904   |4.3%    |111.729    |114.516     |1.263                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1910     |3% of pixels = 768 pixels (3.0%)   |112.524    |4.077   |3.6%    |111.365    |113.683     |0.664                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1910     |Fixed 1,250 (4.9%)                 |111.762    |3.543   |3.2%    |110.755    |112.769     |-0.098                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1910     |Fixed 2,000 (7.8%)                 |112.041    |1.800   |1.6%    |111.529    |112.553     |0.182                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1910     |Fixed 3,000 (11.7%)                |111.745    |2.353   |2.1%    |111.076    |112.413     |-0.115                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |1910     |Fixed 4,000 (15.6%)                |111.859    |1.988   |1.8%    |111.294    |112.424     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |206      |1% of pixels = 233 pixels (1.0%)   |180.228    |12.282  |6.8%    |176.738    |183.719     |0.364                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |206      |2% of pixels = 467 pixels (2.0%)   |177.060    |12.938  |7.3%    |173.383    |180.737     |-2.804                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |206      |3% of pixels = 700 pixels (3.0%)   |180.289    |10.235  |5.7%    |177.381    |183.198     |0.425                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |206      |Fixed 1,250 (5.4%)                 |179.537    |5.434   |3.0%    |177.993    |181.081     |-0.327                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |206      |Fixed 2,000 (8.6%)                 |179.747    |4.547   |2.5%    |178.455    |181.040     |-0.117                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |206      |Fixed 3,000 (12.8%)                |180.042    |4.267   |2.4%    |178.830    |181.255     |0.178                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |206      |Fixed 4,000 (17.1%)                |179.864    |3.202   |1.8%    |178.954    |180.774     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |219      |1% of pixels = 294 pixels (1.0%)   |208.180    |17.227  |8.3%    |203.285    |213.076     |-3.559                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |219      |2% of pixels = 588 pixels (2.0%)   |208.810    |12.385  |5.9%    |205.291    |212.330     |-2.929                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |219      |3% of pixels = 881 pixels (3.0%)   |212.853    |9.637   |4.5%    |210.115    |215.592     |1.114                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |219      |Fixed 1,250 (4.3%)                 |211.905    |7.617   |3.6%    |209.740    |214.070     |0.165                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |219      |Fixed 2,000 (6.8%)                 |212.921    |7.063   |3.3%    |210.914    |214.928     |1.181                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |219      |Fixed 3,000 (10.2%)                |212.779    |5.411   |2.5%    |211.242    |214.317     |1.040                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |219      |Fixed 4,000 (13.6%)                |211.740    |3.572   |1.7%    |210.725    |212.755     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |302      |1% of pixels = 263 pixels (1.0%)   |129.753    |7.202   |5.6%    |127.707    |131.800     |-1.131                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |302      |2% of pixels = 527 pixels (2.0%)   |129.641    |6.680   |5.2%    |127.743    |131.540     |-1.243                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |302      |3% of pixels = 790 pixels (3.0%)   |129.790    |5.230   |4.0%    |128.304    |131.277     |-1.094                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |302      |Fixed 1,250 (4.7%)                 |130.132    |4.914   |3.8%    |128.736    |131.529     |-0.752                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |302      |Fixed 2,000 (7.6%)                 |130.791    |3.011   |2.3%    |129.935    |131.646     |-0.094                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |302      |Fixed 3,000 (11.4%)                |130.988    |2.320   |1.8%    |130.329    |131.647     |0.104                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |302      |Fixed 4,000 (15.2%)                |130.884    |1.776   |1.4%    |130.379    |131.389     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |317      |1% of pixels = 312 pixels (1.0%)   |151.179    |8.531   |5.6%    |148.755    |153.603     |-1.823                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |317      |2% of pixels = 623 pixels (2.0%)   |152.931    |6.893   |4.5%    |150.973    |154.890     |-0.070                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |317      |3% of pixels = 935 pixels (3.0%)   |152.217    |5.692   |3.7%    |150.600    |153.835     |-0.784                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |317      |Fixed 1,250 (4.0%)                 |153.416    |4.590   |3.0%    |152.112    |154.721     |0.415                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |317      |Fixed 2,000 (6.4%)                 |154.242    |3.991   |2.6%    |153.108    |155.376     |1.241                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |317      |Fixed 3,000 (9.6%)                 |152.757    |3.003   |2.0%    |151.904    |153.610     |-0.245                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |317      |Fixed 4,000 (12.8%)                |153.002    |2.529   |1.7%    |152.283    |153.720     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |405      |1% of pixels = 364 pixels (1.0%)   |198.524    |16.968  |8.5%    |193.701    |203.346     |-0.732                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |405      |2% of pixels = 727 pixels (2.0%)   |200.200    |11.439  |5.7%    |196.949    |203.451     |0.944                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |405      |3% of pixels = 1,091 pixels (3.0%) |198.442    |9.626   |4.9%    |195.707    |201.178     |-0.813                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |405      |Fixed 1,250 (3.4%)                 |198.863    |7.188   |3.6%    |196.821    |200.906     |-0.392                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |405      |Fixed 2,000 (5.5%)                 |199.463    |6.401   |3.2%    |197.644    |201.282     |0.208                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |405      |Fixed 3,000 (8.2%)                 |199.025    |4.425   |2.2%    |197.767    |200.282     |-0.231                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |405      |Fixed 4,000 (11.0%)                |199.255    |5.252   |2.6%    |197.763    |200.748     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |821      |1% of pixels = 292 pixels (1.0%)   |172.449    |9.114   |5.3%    |169.859    |175.039     |-0.168                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |821      |2% of pixels = 584 pixels (2.0%)   |172.619    |6.866   |4.0%    |170.668    |174.570     |0.002                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |821      |3% of pixels = 876 pixels (3.0%)   |172.968    |5.621   |3.2%    |171.371    |174.565     |0.351                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |821      |Fixed 1,250 (4.3%)                 |172.417    |5.147   |3.0%    |170.955    |173.880     |-0.200                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |821      |Fixed 2,000 (6.9%)                 |172.750    |3.961   |2.3%    |171.624    |173.875     |0.133                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |821      |Fixed 3,000 (10.3%)                |173.888    |3.437   |2.0%    |172.912    |174.865     |1.271                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |821      |Fixed 4,000 (13.7%)                |172.617    |2.823   |1.6%    |171.815    |173.420     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |905      |1% of pixels = 270 pixels (1.0%)   |464.701    |86.033  |18.5%   |440.250    |489.151     |-18.703                |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |905      |2% of pixels = 540 pixels (2.0%)   |485.172    |86.208  |17.8%   |460.672    |509.672     |1.768                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |905      |3% of pixels = 811 pixels (3.0%)   |484.233    |60.857  |12.6%   |466.938    |501.529     |0.829                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |905      |Fixed 1,250 (4.6%)                 |479.967    |46.196  |9.6%    |466.838    |493.095     |-3.437                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |905      |Fixed 2,000 (7.4%)                 |477.378    |44.293  |9.3%    |464.791    |489.966     |-6.026                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |905      |Fixed 3,000 (11.1%)                |475.429    |30.728  |6.5%    |466.696    |484.162     |-7.975                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |905      |Fixed 4,000 (14.8%)                |483.404    |26.468  |5.5%    |475.882    |490.926     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |912      |1% of pixels = 296 pixels (1.0%)   |210.572    |11.226  |5.3%    |207.381    |213.762     |-0.683                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |912      |2% of pixels = 592 pixels (2.0%)   |213.252    |8.584   |4.0%    |210.813    |215.692     |1.998                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |912      |3% of pixels = 889 pixels (3.0%)   |212.238    |5.204   |2.5%    |210.759    |213.716     |0.983                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |912      |Fixed 1,250 (4.2%)                 |210.799    |5.869   |2.8%    |209.131    |212.467     |-0.455                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |912      |Fixed 2,000 (6.8%)                 |211.044    |5.456   |2.6%    |209.493    |212.595     |-0.210                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |912      |Fixed 3,000 (10.1%)                |211.985    |3.419   |1.6%    |211.014    |212.957     |0.731                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |20m   |912      |Fixed 4,000 (13.5%)                |211.254    |2.916   |1.4%    |210.425    |212.083     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_15 |Fixed 1,250 (0.7%)                 |140.532    |3.788   |2.7%    |139.455    |141.608     |-0.975                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_15 |1% of pixels = 1,728 pixels (1.0%) |141.329    |4.092   |2.9%    |140.166    |142.492     |-0.178                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_15 |2% of pixels = 3,456 pixels (2.0%) |141.760    |2.628   |1.9%    |141.013    |142.507     |0.253                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_15 |Fixed 4,000 (2.3%)                 |141.507    |2.479   |1.8%    |140.802    |142.211     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_15 |3% of pixels = 5,000 pixels (2.9%) |141.493    |2.533   |1.8%    |140.773    |142.213     |-0.013                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_15 |Fixed 6,000 (3.5%)                 |141.911    |2.097   |1.5%    |141.315    |142.507     |0.404                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_15 |Fixed 8,000 (4.6%)                 |141.064    |1.673   |1.2%    |140.588    |141.539     |-0.443                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_24 |Fixed 1,250 (0.7%)                 |158.187    |6.560   |4.1%    |156.323    |160.052     |-1.113                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_24 |1% of pixels = 1,668 pixels (1.0%) |160.147    |5.030   |3.1%    |158.718    |161.577     |0.847                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_24 |2% of pixels = 3,336 pixels (2.0%) |159.303    |3.507   |2.2%    |158.306    |160.300     |0.002                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_24 |Fixed 4,000 (2.4%)                 |159.301    |2.780   |1.7%    |158.510    |160.091     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_24 |3% of pixels = 5,000 pixels (3.0%) |159.236    |2.691   |1.7%    |158.471    |160.001     |-0.064                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_24 |Fixed 6,000 (3.6%)                 |159.296    |2.673   |1.7%    |158.537    |160.056     |-0.004                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_24 |Fixed 8,000 (4.8%)                 |159.677    |2.187   |1.4%    |159.055    |160.299     |0.376                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_36 |Fixed 1,250 (0.7%)                 |213.112    |6.874   |3.2%    |211.158    |215.065     |-1.064                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_36 |1% of pixels = 1,742 pixels (1.0%) |214.904    |6.230   |2.9%    |213.134    |216.675     |0.729                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_36 |2% of pixels = 3,484 pixels (2.0%) |214.307    |3.519   |1.6%    |213.307    |215.307     |0.132                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_36 |Fixed 4,000 (2.3%)                 |214.175    |3.829   |1.8%    |213.087    |215.264     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_36 |3% of pixels = 5,000 pixels (2.9%) |215.133    |3.666   |1.7%    |214.091    |216.175     |0.958                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_36 |Fixed 6,000 (3.4%)                 |214.255    |3.339   |1.6%    |213.306    |215.204     |0.080                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_36 |Fixed 8,000 (4.6%)                 |214.568    |3.244   |1.5%    |213.646    |215.490     |0.392                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_49 |Fixed 1,250 (0.6%)                 |212.531    |13.226  |6.2%    |208.772    |216.289     |-1.782                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_49 |1% of pixels = 2,070 pixels (1.0%) |212.622    |10.452  |4.9%    |209.652    |215.593     |-1.690                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_49 |Fixed 4,000 (1.9%)                 |214.312    |8.195   |3.8%    |211.984    |216.641     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_49 |2% of pixels = 4,139 pixels (2.0%) |211.047    |7.205   |3.4%    |208.999    |213.095     |-3.265                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_49 |3% of pixels = 5,000 pixels (2.4%) |212.562    |6.733   |3.2%    |210.648    |214.476     |-1.751                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_49 |Fixed 6,000 (2.9%)                 |212.778    |5.343   |2.5%    |211.260    |214.297     |-1.534                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_49 |Fixed 8,000 (3.9%)                 |212.593    |5.443   |2.6%    |211.046    |214.140     |-1.719                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_51 |Fixed 1,250 (0.6%)                 |192.845    |6.625   |3.4%    |190.963    |194.728     |-0.082                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_51 |1% of pixels = 2,006 pixels (1.0%) |193.069    |4.476   |2.3%    |191.797    |194.341     |0.141                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_51 |Fixed 4,000 (2.0%)                 |192.927    |3.731   |1.9%    |191.867    |193.988     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_51 |2% of pixels = 4,013 pixels (2.0%) |192.501    |4.103   |2.1%    |191.335    |193.667     |-0.427                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_51 |3% of pixels = 5,000 pixels (2.5%) |191.646    |2.947   |1.5%    |190.808    |192.484     |-1.282                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_51 |Fixed 6,000 (3.0%)                 |192.907    |3.363   |1.7%    |191.951    |193.863     |-0.020                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_51 |Fixed 8,000 (4.0%)                 |192.710    |2.423   |1.3%    |192.021    |193.398     |-0.218                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_56 |Fixed 1,250 (0.9%)                 |169.270    |4.307   |2.5%    |168.046    |170.495     |-1.302                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_56 |1% of pixels = 1,455 pixels (1.0%) |171.837    |4.653   |2.7%    |170.515    |173.159     |1.264                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_56 |2% of pixels = 2,910 pixels (2.0%) |171.241    |3.198   |1.9%    |170.332    |172.150     |0.668                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_56 |Fixed 4,000 (2.7%)                 |170.573    |2.724   |1.6%    |169.799    |171.347     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_56 |3% of pixels = 4,365 pixels (3.0%) |170.923    |2.615   |1.5%    |170.179    |171.666     |0.350                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_56 |Fixed 6,000 (4.1%)                 |170.588    |2.621   |1.5%    |169.844    |171.333     |0.016                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_56 |Fixed 8,000 (5.5%)                 |170.840    |1.946   |1.1%    |170.287    |171.393     |0.267                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_60 |Fixed 1,250 (0.8%)                 |231.662    |8.722   |3.8%    |229.183    |234.141     |1.967                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_60 |1% of pixels = 1,645 pixels (1.0%) |230.331    |7.637   |3.3%    |228.161    |232.502     |0.636                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_60 |2% of pixels = 3,290 pixels (2.0%) |230.783    |5.045   |2.2%    |229.349    |232.216     |1.088                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_60 |Fixed 4,000 (2.4%)                 |229.695    |4.466   |1.9%    |228.426    |230.964     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_60 |3% of pixels = 4,935 pixels (3.0%) |230.648    |4.569   |2.0%    |229.349    |231.946     |0.953                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_60 |Fixed 6,000 (3.6%)                 |231.150    |3.763   |1.6%    |230.081    |232.219     |1.455                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_60 |Fixed 8,000 (4.9%)                 |230.610    |3.084   |1.3%    |229.733    |231.486     |0.915                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_8  |Fixed 1,250 (0.7%)                 |192.631    |6.592   |3.4%    |190.757    |194.504     |-0.541                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_8  |1% of pixels = 1,879 pixels (1.0%) |192.819    |4.950   |2.6%    |191.412    |194.226     |-0.352                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_8  |2% of pixels = 3,758 pixels (2.0%) |193.562    |4.009   |2.1%    |192.423    |194.702     |0.391                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_8  |Fixed 4,000 (2.1%)                 |193.171    |3.294   |1.7%    |192.235    |194.107     |0.000                  |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_8  |3% of pixels = 5,000 pixels (2.7%) |192.669    |2.432   |1.3%    |191.978    |193.360     |-0.502                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_8  |Fixed 6,000 (3.2%)                 |192.653    |2.526   |1.3%    |191.935    |193.370     |-0.519                 |squared_euclidean_pc1_pc3                                                 |
|Spectral Rao's Q  |50m   |sub50_8  |Fixed 8,000 (4.3%)                 |192.740    |2.357   |1.2%    |192.070    |193.410     |-0.432                 |squared_euclidean_pc1_pc3                                                 |
|Alpha-Hull Area   |10m   |1104_b   |1% of pixels = 94 pixels (1.0%)    |14.652     |4.503   |30.7%   |13.372     |15.932      |-411.137               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1104_b   |2% of pixels = 188 pixels (2.0%)   |53.714     |7.175   |13.4%   |51.675     |55.753      |-372.075               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1104_b   |3% of pixels = 283 pixels (3.0%)   |87.877     |7.471   |8.5%    |85.754     |90.001      |-337.912               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1104_b   |Fixed 1,250 (13.3%)                |266.305    |10.188  |3.8%    |263.410    |269.200     |-159.484               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1104_b   |Fixed 2,000 (21.2%)                |329.662    |10.044  |3.0%    |326.808    |332.517     |-96.127                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1104_b   |Fixed 3,000 (31.8%)                |380.997    |10.844  |2.8%    |377.916    |384.079     |-44.792                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1104_b   |Fixed 4,000 (42.4%)                |425.789    |9.870   |2.3%    |422.984    |428.594     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1110_a   |1% of pixels = 76 pixels (1.0%)    |7.475      |2.962   |39.6%   |6.633      |8.317       |-326.677               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1110_a   |2% of pixels = 153 pixels (2.0%)   |39.265     |6.595   |16.8%   |37.391     |41.140      |-294.886               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1110_a   |3% of pixels = 230 pixels (3.0%)   |77.279     |7.623   |9.9%    |75.113     |79.446      |-256.872               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1110_a   |Fixed 1,250 (16.3%)                |252.040    |8.613   |3.4%    |249.592    |254.487     |-82.112                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1110_a   |Fixed 2,000 (26.1%)                |289.941    |6.627   |2.3%    |288.058    |291.824     |-44.211                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1110_a   |Fixed 3,000 (39.2%)                |316.493    |5.859   |1.9%    |314.828    |318.158     |-17.658                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1110_a   |Fixed 4,000 (52.3%)                |334.151    |5.590   |1.7%    |332.563    |335.740     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1124_a   |1% of pixels = 83 pixels (1.0%)    |11.859     |4.062   |34.3%   |10.705     |13.014      |-334.899               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1124_a   |2% of pixels = 166 pixels (2.0%)   |45.034     |5.227   |11.6%   |43.549     |46.520      |-301.724               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1124_a   |3% of pixels = 249 pixels (3.0%)   |77.764     |7.704   |9.9%    |75.575     |79.954      |-268.994               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1124_a   |Fixed 1,250 (15.1%)                |250.613    |7.022   |2.8%    |248.618    |252.609     |-96.145                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1124_a   |Fixed 2,000 (24.1%)                |293.596    |6.941   |2.4%    |291.624    |295.569     |-53.162                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1124_a   |Fixed 3,000 (36.1%)                |327.214    |6.288   |1.9%    |325.427    |329.001     |-19.544                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1124_a   |Fixed 4,000 (48.2%)                |346.758    |6.065   |1.7%    |345.035    |348.482     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |112_b    |1% of pixels = 73 pixels (1.0%)    |8.849      |3.526   |39.8%   |7.847      |9.851       |-456.698               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |112_b    |2% of pixels = 147 pixels (2.0%)   |34.537     |6.219   |18.0%   |32.770     |36.305      |-431.010               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |112_b    |3% of pixels = 220 pixels (3.0%)   |56.881     |5.547   |9.8%    |55.305     |58.458      |-408.666               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |112_b    |Fixed 1,250 (17.0%)                |202.990    |9.974   |4.9%    |200.156    |205.825     |-262.557               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |112_b    |Fixed 2,000 (27.2%)                |279.078    |12.177  |4.4%    |275.618    |282.539     |-186.469               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |112_b    |Fixed 3,000 (40.9%)                |377.672    |13.323  |3.5%    |373.886    |381.459     |-87.875                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |112_b    |Fixed 4,000 (54.5%)                |465.547    |12.188  |2.6%    |462.083    |469.011     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |114_b    |1% of pixels = 72 pixels (1.0%)    |2.413      |1.425   |59.0%   |2.008      |2.818       |-576.135               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |114_b    |2% of pixels = 145 pixels (2.0%)   |17.818     |4.527   |25.4%   |16.531     |19.104      |-560.730               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |114_b    |3% of pixels = 217 pixels (3.0%)   |45.513     |6.263   |13.8%   |43.733     |47.293      |-533.035               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |114_b    |Fixed 1,250 (17.3%)                |339.582    |11.867  |3.5%    |336.210    |342.955     |-238.966               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |114_b    |Fixed 2,000 (27.6%)                |430.019    |13.201  |3.1%    |426.267    |433.771     |-148.529               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |114_b    |Fixed 3,000 (41.4%)                |511.160    |11.179  |2.2%    |507.983    |514.337     |-67.388                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |114_b    |Fixed 4,000 (55.3%)                |578.548    |13.372  |2.3%    |574.748    |582.348     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |11_c     |1% of pixels = 62 pixels (1.0%)    |2.849      |1.466   |51.5%   |2.432      |3.266       |-386.694               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |11_c     |2% of pixels = 124 pixels (2.0%)   |22.437     |4.210   |18.8%   |21.240     |23.633      |-367.106               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |11_c     |3% of pixels = 186 pixels (3.0%)   |48.681     |7.278   |14.9%   |46.612     |50.749      |-340.862               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |11_c     |Fixed 1,250 (20.2%)                |276.050    |9.958   |3.6%    |273.220    |278.880     |-113.493               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |11_c     |Fixed 2,000 (32.3%)                |325.517    |8.467   |2.6%    |323.111    |327.924     |-64.026                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |11_c     |Fixed 3,000 (48.4%)                |363.522    |7.956   |2.2%    |361.261    |365.783     |-26.021                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |11_c     |Fixed 4,000 (64.5%)                |389.543    |7.103   |1.8%    |387.524    |391.562     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |120_d    |1% of pixels = 71 pixels (1.0%)    |2.112      |1.459   |69.1%   |1.697      |2.526       |-577.026               |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |10m   |120_d    |2% of pixels = 141 pixels (2.0%)   |15.494     |3.931   |25.4%   |14.377     |16.611      |-563.644               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |120_d    |3% of pixels = 212 pixels (3.0%)   |41.886     |7.942   |19.0%   |39.629     |44.143      |-537.252               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |120_d    |Fixed 1,250 (17.7%)                |367.584    |11.576  |3.1%    |364.294    |370.874     |-211.554               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |120_d    |Fixed 2,000 (28.3%)                |454.629    |12.981  |2.9%    |450.940    |458.318     |-124.509               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |120_d    |Fixed 3,000 (42.4%)                |530.346    |11.545  |2.2%    |527.065    |533.627     |-48.792                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |120_d    |Fixed 4,000 (56.6%)                |579.138    |9.874   |1.7%    |576.332    |581.944     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1305_c   |1% of pixels = 80 pixels (1.0%)    |8.791      |3.118   |35.5%   |7.904      |9.677       |-370.775               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1305_c   |2% of pixels = 161 pixels (2.0%)   |43.184     |6.801   |15.7%   |41.251     |45.116      |-336.382               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1305_c   |3% of pixels = 241 pixels (3.0%)   |80.795     |7.638   |9.5%    |78.624     |82.965      |-298.771               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1305_c   |Fixed 1,250 (15.6%)                |254.389    |9.059   |3.6%    |251.814    |256.964     |-125.176               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1305_c   |Fixed 2,000 (24.9%)                |305.753    |7.435   |2.4%    |303.640    |307.867     |-73.812                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1305_c   |Fixed 3,000 (37.4%)                |351.192    |9.179   |2.6%    |348.583    |353.801     |-28.373                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1305_c   |Fixed 4,000 (49.8%)                |379.565    |9.055   |2.4%    |376.992    |382.139     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1315_b   |1% of pixels = 91 pixels (1.0%)    |5.732      |2.598   |45.3%   |4.993      |6.470       |-576.279               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1315_b   |2% of pixels = 182 pixels (2.0%)   |35.030     |7.506   |21.4%   |32.896     |37.163      |-546.982               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1315_b   |3% of pixels = 272 pixels (3.0%)   |75.846     |8.104   |10.7%   |73.543     |78.149      |-506.165               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1315_b   |Fixed 1,250 (13.8%)                |339.677    |14.172  |4.2%    |335.649    |343.705     |-242.334               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1315_b   |Fixed 2,000 (22.0%)                |437.796    |12.687  |2.9%    |434.191    |441.402     |-144.215               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1315_b   |Fixed 3,000 (33.1%)                |521.237    |12.938  |2.5%    |517.560    |524.914     |-60.775                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1315_b   |Fixed 4,000 (44.1%)                |582.011    |13.443  |2.3%    |578.191    |585.832     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1400_b   |1% of pixels = 68 pixels (1.0%)    |2.956      |1.285   |43.5%   |2.591      |3.321       |-457.530               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1400_b   |2% of pixels = 136 pixels (2.0%)   |20.441     |4.002   |19.6%   |19.303     |21.578      |-440.046               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1400_b   |3% of pixels = 204 pixels (3.0%)   |50.504     |7.811   |15.5%   |48.284     |52.724      |-409.982               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1400_b   |Fixed 1,250 (18.4%)                |308.596    |10.189  |3.3%    |305.700    |311.491     |-151.890               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1400_b   |Fixed 2,000 (29.4%)                |367.187    |8.905   |2.4%    |364.656    |369.717     |-93.300                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1400_b   |Fixed 3,000 (44.1%)                |419.252    |8.341   |2.0%    |416.882    |421.623     |-41.234                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1400_b   |Fixed 4,000 (58.8%)                |460.486    |10.855  |2.4%    |457.401    |463.571     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1414_a   |1% of pixels = 62 pixels (1.0%)    |1.137      |1.107   |97.3%   |0.823      |1.452       |-620.028               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1414_a   |2% of pixels = 125 pixels (2.0%)   |9.588      |3.403   |35.5%   |8.621      |10.555      |-611.577               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1414_a   |3% of pixels = 187 pixels (3.0%)   |26.591     |6.077   |22.9%   |24.864     |28.318      |-594.574               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1414_a   |Fixed 1,250 (20.1%)                |392.265    |13.865  |3.5%    |388.324    |396.205     |-228.901               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1414_a   |Fixed 2,000 (32.1%)                |484.607    |12.850  |2.7%    |480.955    |488.259     |-136.559               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1414_a   |Fixed 3,000 (48.2%)                |566.533    |10.549  |1.9%    |563.535    |569.531     |-54.633                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1414_a   |Fixed 4,000 (64.2%)                |621.166    |9.600   |1.5%    |618.437    |623.894     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1417_c   |1% of pixels = 75 pixels (1.0%)    |1.296      |0.974   |75.1%   |1.019      |1.573       |-689.556               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1417_c   |2% of pixels = 150 pixels (2.0%)   |11.783     |3.624   |30.8%   |10.753     |12.813      |-679.069               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1417_c   |3% of pixels = 225 pixels (3.0%)   |37.073     |8.478   |22.9%   |34.664     |39.483      |-653.779               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1417_c   |Fixed 1,250 (16.6%)                |413.010    |14.275  |3.5%    |408.910    |417.111     |-277.842               |all_pixels; fallback_sampled_pixels_alpha_failed                          |
|Alpha-Hull Area   |10m   |1417_c   |Fixed 2,000 (26.6%)                |539.943    |12.507  |2.3%    |536.388    |543.497     |-150.909               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1417_c   |Fixed 3,000 (39.9%)                |635.042    |15.469  |2.4%    |630.646    |639.439     |-55.810                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1417_c   |Fixed 4,000 (53.2%)                |690.852    |12.594  |1.8%    |687.235    |694.469     |0.000                  |all_pixels; fallback_sampled_pixels_alpha_failed                          |
|Alpha-Hull Area   |10m   |1514_a   |1% of pixels = 70 pixels (1.0%)    |7.253      |3.340   |46.0%   |6.304      |8.202       |-397.624               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1514_a   |2% of pixels = 141 pixels (2.0%)   |35.838     |7.348   |20.5%   |33.749     |37.926      |-369.040               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1514_a   |3% of pixels = 211 pixels (3.0%)   |67.057     |5.646   |8.4%    |65.453     |68.662      |-337.820               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1514_a   |Fixed 1,250 (17.8%)                |244.773    |10.183  |4.2%    |241.879    |247.667     |-160.105               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1514_a   |Fixed 2,000 (28.5%)                |307.335    |13.684  |4.5%    |303.446    |311.224     |-97.543                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1514_a   |Fixed 3,000 (42.7%)                |367.913    |8.600   |2.3%    |365.469    |370.357     |-36.964                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1514_a   |Fixed 4,000 (56.9%)                |404.878    |7.666   |1.9%    |402.699    |407.056     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1516_c   |1% of pixels = 84 pixels (1.0%)    |3.300      |1.937   |58.7%   |2.750      |3.850       |-583.425               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1516_c   |2% of pixels = 167 pixels (2.0%)   |21.024     |5.173   |24.6%   |19.554     |22.494      |-565.702               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1516_c   |3% of pixels = 251 pixels (3.0%)   |59.540     |8.577   |14.4%   |57.103     |61.978      |-527.185               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1516_c   |Fixed 1,250 (15.0%)                |372.095    |13.582  |3.7%    |368.235    |375.955     |-214.630               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1516_c   |Fixed 2,000 (23.9%)                |459.888    |13.040  |2.8%    |456.182    |463.594     |-126.837               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1516_c   |Fixed 3,000 (35.9%)                |536.607    |10.675  |2.0%    |533.573    |539.641     |-50.118                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1516_c   |Fixed 4,000 (47.9%)                |586.725    |9.318   |1.6%    |584.077    |589.373     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1604_c   |1% of pixels = 72 pixels (1.0%)    |4.396      |2.339   |53.2%   |3.731      |5.060       |-449.559               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1604_c   |2% of pixels = 145 pixels (2.0%)   |27.387     |6.256   |22.8%   |25.609     |29.165      |-426.568               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1604_c   |3% of pixels = 217 pixels (3.0%)   |61.981     |8.766   |14.1%   |59.489     |64.472      |-391.974               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1604_c   |Fixed 1,250 (17.3%)                |305.891    |10.218  |3.3%    |302.987    |308.795     |-148.064               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1604_c   |Fixed 2,000 (27.7%)                |370.384    |8.025   |2.2%    |368.103    |372.665     |-83.571                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1604_c   |Fixed 3,000 (41.5%)                |420.531    |7.136   |1.7%    |418.503    |422.560     |-33.423                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1604_c   |Fixed 4,000 (55.3%)                |453.955    |8.064   |1.8%    |451.663    |456.247     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1606_c   |1% of pixels = 76 pixels (1.0%)    |3.143      |1.563   |49.7%   |2.699      |3.587       |-475.284               |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |10m   |1606_c   |2% of pixels = 152 pixels (2.0%)   |23.310     |5.302   |22.7%   |21.803     |24.816      |-455.118               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1606_c   |3% of pixels = 228 pixels (3.0%)   |58.371     |9.165   |15.7%   |55.767     |60.976      |-420.056               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1606_c   |Fixed 1,250 (16.5%)                |336.839    |10.435  |3.1%    |333.874    |339.805     |-141.588               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1606_c   |Fixed 2,000 (26.3%)                |395.442    |10.476  |2.6%    |392.464    |398.419     |-82.986                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1606_c   |Fixed 3,000 (39.5%)                |441.313    |10.489  |2.4%    |438.332    |444.294     |-37.114                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1606_c   |Fixed 4,000 (52.7%)                |478.427    |8.754   |1.8%    |475.940    |480.915     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1701_b   |1% of pixels = 93 pixels (1.0%)    |10.659     |3.633   |34.1%   |9.627      |11.692      |-477.330               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1701_b   |2% of pixels = 186 pixels (2.0%)   |45.988     |6.298   |13.7%   |44.198     |47.778      |-442.002               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1701_b   |3% of pixels = 279 pixels (3.0%)   |87.352     |8.112   |9.3%    |85.047     |89.658      |-400.637               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1701_b   |Fixed 1,250 (13.4%)                |295.465    |11.590  |3.9%    |292.171    |298.759     |-192.525               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1701_b   |Fixed 2,000 (21.5%)                |374.048    |13.880  |3.7%    |370.103    |377.992     |-113.942               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1701_b   |Fixed 3,000 (32.3%)                |441.126    |10.894  |2.5%    |438.030    |444.222     |-46.863                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1701_b   |Fixed 4,000 (43.0%)                |487.989    |8.068   |1.7%    |485.696    |490.282     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1803_c   |1% of pixels = 53 pixels (1.0%)    |1.942      |1.458   |75.1%   |1.527      |2.356       |-382.192               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1803_c   |2% of pixels = 106 pixels (2.0%)   |14.647     |3.826   |26.1%   |13.559     |15.734      |-369.486               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1803_c   |3% of pixels = 159 pixels (3.0%)   |37.314     |6.137   |16.4%   |35.570     |39.058      |-346.820               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1803_c   |Fixed 1,250 (23.6%)                |274.096    |7.088   |2.6%    |272.081    |276.110     |-110.037               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1803_c   |Fixed 2,000 (37.8%)                |317.999    |9.019   |2.8%    |315.436    |320.562     |-66.134                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1803_c   |Fixed 3,000 (56.7%)                |353.936    |6.073   |1.7%    |352.210    |355.662     |-30.197                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1803_c   |Fixed 4,000 (75.5%)                |384.133    |6.019   |1.6%    |382.423    |385.844     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1816_c   |1% of pixels = 52 pixels (1.0%)    |1.745      |1.218   |69.8%   |1.399      |2.091       |-451.269               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1816_c   |2% of pixels = 104 pixels (2.0%)   |13.658     |4.361   |31.9%   |12.419     |14.898      |-439.355               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1816_c   |3% of pixels = 156 pixels (3.0%)   |34.536     |6.701   |19.4%   |32.632     |36.441      |-418.477               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1816_c   |Fixed 1,250 (24.1%)                |279.968    |11.126  |4.0%    |276.772    |283.163     |-173.046               |all_pixels; fallback_sampled_pixels_alpha_failed                          |
|Alpha-Hull Area   |10m   |1816_c   |Fixed 2,000 (38.5%)                |342.729    |6.564   |1.9%    |340.780    |344.678     |-110.285               |all_pixels; fallback_sampled_pixels; fallback_sampled_pixels_alpha_failed |
|Alpha-Hull Area   |10m   |1816_c   |Fixed 3,000 (57.8%)                |402.651    |4.163   |1.0%    |400.970    |404.332     |-50.363                |all_pixels; fallback_sampled_pixels; fallback_sampled_pixels_alpha_failed |
|Alpha-Hull Area   |10m   |1816_c   |Fixed 4,000 (77.0%)                |453.014    |9.830   |2.2%    |450.190    |455.837     |0.000                  |all_pixels; fallback_sampled_pixels; fallback_sampled_pixels_alpha_failed |
|Alpha-Hull Area   |10m   |1912_c   |1% of pixels = 55 pixels (1.0%)    |1.570      |1.148   |73.1%   |1.244      |1.896       |-401.988               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1912_c   |2% of pixels = 109 pixels (2.0%)   |14.036     |4.584   |32.7%   |12.733     |15.339      |-389.522               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1912_c   |3% of pixels = 164 pixels (3.0%)   |36.464     |6.361   |17.4%   |34.656     |38.272      |-367.094               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1912_c   |Fixed 1,250 (22.9%)                |292.647    |8.355   |2.9%    |290.272    |295.021     |-110.911               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1912_c   |Fixed 2,000 (36.6%)                |337.687    |8.796   |2.6%    |335.187    |340.187     |-65.870                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1912_c   |Fixed 3,000 (54.9%)                |376.528    |6.596   |1.8%    |374.653    |378.402     |-27.030                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1912_c   |Fixed 4,000 (73.2%)                |403.557    |5.374   |1.3%    |402.030    |405.085     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1917_c   |1% of pixels = 55 pixels (1.0%)    |1.423      |1.025   |72.1%   |1.132      |1.715       |-476.250               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1917_c   |2% of pixels = 109 pixels (2.0%)   |12.642     |4.251   |33.6%   |11.434     |13.850      |-465.032               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1917_c   |3% of pixels = 164 pixels (3.0%)   |33.355     |6.048   |18.1%   |31.636     |35.073      |-444.319               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1917_c   |Fixed 1,250 (22.9%)                |314.302    |10.030  |3.2%    |311.451    |317.152     |-163.372               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1917_c   |Fixed 2,000 (36.7%)                |386.986    |8.964   |2.3%    |384.439    |389.534     |-90.687                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1917_c   |Fixed 3,000 (55.0%)                |443.074    |9.175   |2.1%    |440.466    |445.681     |-34.600                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |1917_c   |Fixed 4,000 (73.4%)                |477.674    |6.181   |1.3%    |475.917    |479.430     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |204_d    |1% of pixels = 71 pixels (1.0%)    |2.847      |1.490   |52.3%   |2.424      |3.271       |-483.291               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |204_d    |2% of pixels = 143 pixels (2.0%)   |21.258     |4.889   |23.0%   |19.869     |22.648      |-464.880               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |204_d    |3% of pixels = 214 pixels (3.0%)   |54.031     |8.474   |15.7%   |51.623     |56.439      |-432.108               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |204_d    |Fixed 1,250 (17.5%)                |325.899    |10.234  |3.1%    |322.990    |328.807     |-160.240               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |204_d    |Fixed 2,000 (28.0%)                |394.720    |10.617  |2.7%    |391.703    |397.737     |-91.419                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |204_d    |Fixed 3,000 (42.0%)                |446.683    |9.400   |2.1%    |444.012    |449.354     |-39.456                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |204_d    |Fixed 4,000 (56.0%)                |486.138    |7.906   |1.6%    |483.891    |488.385     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |217_d    |1% of pixels = 78 pixels (1.0%)    |5.019      |2.291   |45.6%   |4.368      |5.670       |-453.278               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |217_d    |2% of pixels = 157 pixels (2.0%)   |29.565     |5.939   |20.1%   |27.877     |31.253      |-428.732               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |217_d    |3% of pixels = 236 pixels (3.0%)   |67.881     |8.245   |12.1%   |65.538     |70.225      |-390.415               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |217_d    |Fixed 1,250 (15.9%)                |310.215    |10.431  |3.4%    |307.251    |313.180     |-148.082               |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |10m   |217_d    |Fixed 2,000 (25.5%)                |370.767    |8.942   |2.4%    |368.226    |373.308     |-87.530                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |217_d    |Fixed 3,000 (38.2%)                |422.598    |9.220   |2.2%    |419.977    |425.218     |-35.699                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |217_d    |Fixed 4,000 (51.0%)                |458.297    |8.699   |1.9%    |455.825    |460.769     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |22_c     |1% of pixels = 71 pixels (1.0%)    |2.688      |1.654   |61.5%   |2.218      |3.158       |-624.810               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |22_c     |2% of pixels = 143 pixels (2.0%)   |21.326     |6.091   |28.6%   |19.595     |23.057      |-606.172               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |22_c     |3% of pixels = 214 pixels (3.0%)   |47.065     |5.885   |12.5%   |45.393     |48.738      |-580.433               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |22_c     |Fixed 1,250 (17.5%)                |333.120    |14.240  |4.3%    |329.073    |337.167     |-294.379               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |22_c     |Fixed 2,000 (28.0%)                |447.856    |14.023  |3.1%    |443.871    |451.841     |-179.642               |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |10m   |22_c     |Fixed 3,000 (42.0%)                |556.793    |14.260  |2.6%    |552.740    |560.845     |-70.705                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |22_c     |Fixed 4,000 (56.1%)                |627.498    |13.479  |2.1%    |623.668    |631.329     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |319_c    |1% of pixels = 90 pixels (1.0%)    |11.447     |4.054   |35.4%   |10.295     |12.599      |-397.934               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |319_c    |2% of pixels = 180 pixels (2.0%)   |44.772     |6.139   |13.7%   |43.028     |46.517      |-364.608               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |319_c    |3% of pixels = 269 pixels (3.0%)   |78.447     |7.895   |10.1%   |76.203     |80.691      |-330.933               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |319_c    |Fixed 1,250 (13.9%)                |281.561    |9.868   |3.5%    |278.757    |284.366     |-127.819               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |319_c    |Fixed 2,000 (22.3%)                |338.958    |11.257  |3.3%    |335.758    |342.157     |-70.423                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |319_c    |Fixed 3,000 (33.4%)                |381.078    |7.000   |1.8%    |379.089    |383.067     |-28.302                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |319_c    |Fixed 4,000 (44.5%)                |409.380    |8.387   |2.0%    |406.997    |411.764     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |409_d    |1% of pixels = 60 pixels (1.0%)    |1.030      |0.832   |80.8%   |0.794      |1.267       |-596.723               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |409_d    |2% of pixels = 120 pixels (2.0%)   |9.472      |3.292   |34.7%   |8.537      |10.408      |-588.281               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |409_d    |3% of pixels = 180 pixels (3.0%)   |26.059     |5.823   |22.3%   |24.404     |27.714      |-571.694               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |409_d    |Fixed 1,250 (20.9%)                |372.785    |12.711  |3.4%    |369.173    |376.398     |-224.968               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |409_d    |Fixed 2,000 (33.4%)                |466.459    |11.133  |2.4%    |463.295    |469.623     |-131.294               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |409_d    |Fixed 3,000 (50.1%)                |541.288    |11.477  |2.1%    |538.026    |544.549     |-56.466                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |409_d    |Fixed 4,000 (66.8%)                |597.753    |7.849   |1.3%    |595.523    |599.984     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |503_c    |1% of pixels = 74 pixels (1.0%)    |4.166      |1.884   |45.2%   |3.631      |4.702       |-383.331               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |503_c    |2% of pixels = 148 pixels (2.0%)   |29.855     |5.786   |19.4%   |28.211     |31.500      |-357.642               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |503_c    |3% of pixels = 222 pixels (3.0%)   |67.217     |8.764   |13.0%   |64.727     |69.708      |-320.280               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |503_c    |Fixed 1,250 (16.9%)                |289.615    |7.953   |2.7%    |287.354    |291.875     |-97.883                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |503_c    |Fixed 2,000 (27.0%)                |333.171    |7.541   |2.3%    |331.028    |335.315     |-54.326                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |503_c    |Fixed 3,000 (40.5%)                |365.202    |7.109   |1.9%    |363.182    |367.222     |-22.295                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |503_c    |Fixed 4,000 (54.0%)                |387.498    |5.808   |1.5%    |385.847    |389.148     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |614_a    |1% of pixels = 71 pixels (1.0%)    |9.312      |3.486   |37.4%   |8.321      |10.303      |-329.944               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |614_a    |2% of pixels = 142 pixels (2.0%)   |38.843     |5.659   |14.6%   |37.235     |40.452      |-300.413               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |614_a    |3% of pixels = 213 pixels (3.0%)   |67.287     |7.681   |11.4%   |65.104     |69.470      |-271.969               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |614_a    |Fixed 1,250 (17.6%)                |234.983    |8.675   |3.7%    |232.517    |237.448     |-104.273               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |614_a    |Fixed 2,000 (28.2%)                |276.390    |8.055   |2.9%    |274.100    |278.679     |-62.866                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |614_a    |Fixed 3,000 (42.2%)                |312.208    |7.219   |2.3%    |310.156    |314.260     |-27.048                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |614_a    |Fixed 4,000 (56.3%)                |339.256    |5.598   |1.7%    |337.665    |340.847     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |700_c    |1% of pixels = 76 pixels (1.0%)    |6.154      |2.133   |34.7%   |5.548      |6.760       |-394.071               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |700_c    |2% of pixels = 152 pixels (2.0%)   |32.526     |5.936   |18.3%   |30.839     |34.213      |-367.699               |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |10m   |700_c    |3% of pixels = 229 pixels (3.0%)   |68.874     |7.070   |10.3%   |66.865     |70.883      |-331.351               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |700_c    |Fixed 1,250 (16.4%)                |281.109    |9.372   |3.3%    |278.446    |283.773     |-119.116               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |700_c    |Fixed 2,000 (26.2%)                |333.074    |8.176   |2.5%    |330.750    |335.398     |-67.152                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |700_c    |Fixed 3,000 (39.4%)                |372.634    |7.591   |2.0%    |370.477    |374.792     |-27.591                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |700_c    |Fixed 4,000 (52.5%)                |400.226    |5.670   |1.4%    |398.614    |401.837     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |722_a    |1% of pixels = 88 pixels (1.0%)    |4.676      |2.179   |46.6%   |4.057      |5.296       |-529.713               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |722_a    |2% of pixels = 176 pixels (2.0%)   |32.706     |6.939   |21.2%   |30.734     |34.678      |-501.683               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |722_a    |3% of pixels = 263 pixels (3.0%)   |74.033     |8.496   |11.5%   |71.619     |76.448      |-460.356               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |722_a    |Fixed 1,250 (14.2%)                |354.136    |13.230  |3.7%    |350.376    |357.896     |-180.254               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |722_a    |Fixed 2,000 (22.8%)                |434.253    |12.963  |3.0%    |430.569    |437.937     |-100.136               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |722_a    |Fixed 3,000 (34.2%)                |499.688    |9.208   |1.8%    |497.071    |502.305     |-34.701                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |722_a    |Fixed 4,000 (45.6%)                |534.389    |7.219   |1.4%    |532.338    |536.441     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |800_a    |1% of pixels = 40 pixels (1.0%)    |0.415      |0.538   |129.4%  |0.263      |0.568       |-557.007               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |800_a    |2% of pixels = 80 pixels (2.0%)    |3.200      |1.628   |50.9%   |2.737      |3.663       |-554.222               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |800_a    |3% of pixels = 119 pixels (3.0%)   |10.639     |3.461   |32.5%   |9.656      |11.623      |-546.783               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |800_a    |Fixed 1,250 (31.4%)                |320.756    |10.869  |3.4%    |317.667    |323.845     |-236.666               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |800_a    |Fixed 2,000 (50.3%)                |415.190    |11.460  |2.8%    |411.933    |418.447     |-142.232               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |800_a    |Fixed 3,000 (75.5%)                |500.026    |8.397   |1.7%    |497.639    |502.412     |-57.397                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |800_a    |Fixed 4,000 (100.0%)               |557.422    |0.000   |0.0%    |557.422    |557.422     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |10m   |914_a    |1% of pixels = 72 pixels (1.0%)    |6.177      |2.799   |45.3%   |5.382      |6.973       |-441.640               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |914_a    |2% of pixels = 144 pixels (2.0%)   |31.670     |5.072   |16.0%   |30.228     |33.111      |-416.148               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |914_a    |3% of pixels = 217 pixels (3.0%)   |69.749     |8.848   |12.7%   |67.235     |72.264      |-378.068               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |914_a    |Fixed 1,250 (17.3%)                |253.963    |9.463   |3.7%    |251.274    |256.652     |-193.854               |all_pixels                                                                |
|Alpha-Hull Area   |10m   |914_a    |Fixed 2,000 (27.7%)                |324.052    |11.322  |3.5%    |320.835    |327.270     |-123.765               |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |10m   |914_a    |Fixed 3,000 (41.6%)                |391.735    |11.147  |2.8%    |388.568    |394.903     |-56.082                |all_pixels                                                                |
|Alpha-Hull Area   |10m   |914_a    |Fixed 4,000 (55.4%)                |447.817    |10.308  |2.3%    |444.888    |450.746     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |100      |1% of pixels = 308 pixels (1.0%)   |90.493     |9.823   |10.9%   |87.701     |93.284      |-482.148               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |100      |2% of pixels = 616 pixels (2.0%)   |227.738    |10.765  |4.7%    |224.679    |230.798     |-344.903               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |100      |3% of pixels = 924 pixels (3.0%)   |312.573    |11.729  |3.8%    |309.239    |315.906     |-260.069               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |100      |Fixed 1,250 (4.1%)                 |371.716    |11.503  |3.1%    |368.447    |374.985     |-200.925               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |100      |Fixed 2,000 (6.5%)                 |456.344    |12.655  |2.8%    |452.748    |459.941     |-116.297               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |100      |Fixed 3,000 (9.7%)                 |523.848    |13.134  |2.5%    |520.116    |527.581     |-48.793                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |100      |Fixed 4,000 (13.0%)                |572.641    |11.780  |2.1%    |569.293    |575.989     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1201     |1% of pixels = 248 pixels (1.0%)   |67.266     |9.233   |13.7%   |64.642     |69.890      |-419.093               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1201     |2% of pixels = 496 pixels (2.0%)   |189.732    |11.472  |6.0%    |186.472    |192.993     |-296.626               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1201     |3% of pixels = 744 pixels (3.0%)   |271.522    |13.225  |4.9%    |267.764    |275.281     |-214.836               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1201     |Fixed 1,250 (5.0%)                 |345.026    |11.172  |3.2%    |341.851    |348.201     |-141.332               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1201     |Fixed 2,000 (8.1%)                 |402.843    |11.715  |2.9%    |399.514    |406.172     |-83.515                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1201     |Fixed 3,000 (12.1%)                |452.170    |11.545  |2.6%    |448.889    |455.451     |-34.189                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1201     |Fixed 4,000 (16.1%)                |486.358    |9.005   |1.9%    |483.799    |488.917     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1402     |1% of pixels = 271 pixels (1.0%)   |91.918     |11.397  |12.4%   |88.679     |95.157      |-341.986               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1402     |2% of pixels = 541 pixels (2.0%)   |188.851    |10.038  |5.3%    |185.998    |191.704     |-245.053               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1402     |3% of pixels = 812 pixels (3.0%)   |242.991    |10.818  |4.5%    |239.917    |246.066     |-190.913               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1402     |Fixed 1,250 (4.6%)                 |298.997    |10.304  |3.4%    |296.068    |301.925     |-134.907               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1402     |Fixed 2,000 (7.4%)                 |352.906    |9.505   |2.7%    |350.205    |355.608     |-80.998                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1402     |Fixed 3,000 (11.1%)                |402.732    |10.316  |2.6%    |399.800    |405.663     |-31.172                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1402     |Fixed 4,000 (14.8%)                |433.904    |9.075   |2.1%    |431.325    |436.483     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1512     |1% of pixels = 310 pixels (1.0%)   |83.007     |9.656   |11.6%   |80.263     |85.751      |-514.448               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1512     |2% of pixels = 620 pixels (2.0%)   |232.065    |13.297  |5.7%    |228.286    |235.844     |-365.391               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1512     |3% of pixels = 930 pixels (3.0%)   |327.108    |14.090  |4.3%    |323.104    |331.113     |-270.347               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1512     |Fixed 1,250 (4.0%)                 |392.511    |12.934  |3.3%    |388.835    |396.187     |-204.945               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1512     |Fixed 2,000 (6.5%)                 |481.205    |11.332  |2.4%    |477.984    |484.425     |-116.251               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1512     |Fixed 3,000 (9.7%)                 |547.162    |11.970  |2.2%    |543.760    |550.564     |-50.293                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1512     |Fixed 4,000 (12.9%)                |597.455    |11.687  |2.0%    |594.134    |600.777     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1518     |1% of pixels = 357 pixels (1.0%)   |100.670    |12.875  |12.8%   |97.011     |104.329     |-562.727               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1518     |2% of pixels = 714 pixels (2.0%)   |234.129    |14.026  |6.0%    |230.143    |238.115     |-429.268               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1518     |3% of pixels = 1,071 pixels (3.0%) |319.537    |15.462  |4.8%    |315.143    |323.931     |-343.860               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1518     |Fixed 1,250 (3.5%)                 |355.716    |13.693  |3.8%    |351.824    |359.607     |-307.681               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1518     |Fixed 2,000 (5.6%)                 |463.687    |17.304  |3.7%    |458.769    |468.605     |-199.710               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1518     |Fixed 3,000 (8.4%)                 |572.715    |16.254  |2.8%    |568.095    |577.334     |-90.682                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1518     |Fixed 4,000 (11.2%)                |663.397    |18.938  |2.9%    |658.015    |668.779     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1800     |1% of pixels = 300 pixels (1.0%)   |87.079     |10.166  |11.7%   |84.190     |89.968      |-459.471               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1800     |2% of pixels = 600 pixels (2.0%)   |227.397    |12.885  |5.7%    |223.735    |231.059     |-319.152               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1800     |3% of pixels = 900 pixels (3.0%)   |303.791    |10.166  |3.3%    |300.902    |306.680     |-242.759               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1800     |Fixed 1,250 (4.2%)                 |364.049    |12.818  |3.5%    |360.406    |367.692     |-182.500               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1800     |Fixed 2,000 (6.7%)                 |442.557    |10.363  |2.3%    |439.612    |445.502     |-103.992               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1800     |Fixed 3,000 (10.0%)                |502.383    |9.306   |1.9%    |499.738    |505.027     |-44.167                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1800     |Fixed 4,000 (13.3%)                |546.550    |12.203  |2.2%    |543.082    |550.018     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1805     |1% of pixels = 254 pixels (1.0%)   |72.718     |7.264   |10.0%   |70.653     |74.782      |-415.938               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1805     |2% of pixels = 508 pixels (2.0%)   |175.180    |11.285  |6.4%    |171.973    |178.387     |-313.476               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1805     |3% of pixels = 763 pixels (3.0%)   |242.834    |9.454   |3.9%    |240.147    |245.520     |-245.822               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1805     |Fixed 1,250 (4.9%)                 |314.523    |11.999  |3.8%    |311.113    |317.933     |-174.133               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1805     |Fixed 2,000 (7.9%)                 |384.409    |8.467   |2.2%    |382.003    |386.815     |-104.247               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1805     |Fixed 3,000 (11.8%)                |442.925    |11.179  |2.5%    |439.748    |446.102     |-45.731                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1805     |Fixed 4,000 (15.7%)                |488.656    |12.394  |2.5%    |485.134    |492.178     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1910     |1% of pixels = 256 pixels (1.0%)   |88.462     |7.922   |9.0%    |86.211     |90.713      |-267.536               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1910     |2% of pixels = 512 pixels (2.0%)   |158.458    |8.397   |5.3%    |156.071    |160.844     |-197.540               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1910     |3% of pixels = 768 pixels (3.0%)   |194.072    |7.111   |3.7%    |192.051    |196.092     |-161.926               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1910     |Fixed 1,250 (4.9%)                 |239.722    |9.187   |3.8%    |237.111    |242.333     |-116.276               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1910     |Fixed 2,000 (7.8%)                 |283.892    |9.042   |3.2%    |281.322    |286.462     |-72.106                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1910     |Fixed 3,000 (11.7%)                |325.223    |9.101   |2.8%    |322.637    |327.810     |-30.774                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |1910     |Fixed 4,000 (15.6%)                |355.998    |8.803   |2.5%    |353.496    |358.500     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |206      |1% of pixels = 233 pixels (1.0%)   |61.637     |8.628   |14.0%   |59.184     |64.089      |-471.582               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |206      |2% of pixels = 467 pixels (2.0%)   |156.927    |10.822  |6.9%    |153.851    |160.002     |-376.292               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |206      |3% of pixels = 700 pixels (3.0%)   |223.997    |12.615  |5.6%    |220.412    |227.582     |-309.222               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |206      |Fixed 1,250 (5.4%)                 |322.426    |13.909  |4.3%    |318.473    |326.379     |-210.793               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |206      |Fixed 2,000 (8.6%)                 |407.912    |11.148  |2.7%    |404.744    |411.080     |-125.307               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |206      |Fixed 3,000 (12.8%)                |477.650    |13.903  |2.9%    |473.699    |481.602     |-55.569                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |206      |Fixed 4,000 (17.1%)                |533.219    |12.162  |2.3%    |529.762    |536.676     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |219      |1% of pixels = 294 pixels (1.0%)   |82.588     |11.104  |13.4%   |79.433     |85.744      |-486.029               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |219      |2% of pixels = 588 pixels (2.0%)   |208.177    |12.411  |6.0%    |204.650    |211.704     |-360.441               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |219      |3% of pixels = 881 pixels (3.0%)   |281.053    |12.798  |4.6%    |277.416    |284.690     |-287.565               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |219      |Fixed 1,250 (4.3%)                 |350.277    |13.029  |3.7%    |346.574    |353.980     |-218.341               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |219      |Fixed 2,000 (6.8%)                 |443.463    |13.736  |3.1%    |439.559    |447.367     |-125.155               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |219      |Fixed 3,000 (10.2%)                |518.410    |12.800  |2.5%    |514.772    |522.048     |-50.208                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |219      |Fixed 4,000 (13.6%)                |568.618    |13.894  |2.4%    |564.669    |572.566     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |302      |1% of pixels = 263 pixels (1.0%)   |84.222     |7.805   |9.3%    |82.004     |86.440      |-335.907               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |302      |2% of pixels = 527 pixels (2.0%)   |165.075    |9.632   |5.8%    |162.338    |167.812     |-255.054               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |302      |3% of pixels = 790 pixels (3.0%)   |212.510    |9.547   |4.5%    |209.797    |215.224     |-207.618               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |302      |Fixed 1,250 (4.7%)                 |263.142    |10.672  |4.1%    |260.109    |266.175     |-156.987               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |302      |Fixed 2,000 (7.6%)                 |324.349    |10.641  |3.3%    |321.325    |327.373     |-95.780                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |302      |Fixed 3,000 (11.4%)                |382.179    |11.133  |2.9%    |379.015    |385.343     |-37.950                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |302      |Fixed 4,000 (15.2%)                |420.129    |11.321  |2.7%    |416.912    |423.346     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |317      |1% of pixels = 312 pixels (1.0%)   |102.175    |11.372  |11.1%   |98.944     |105.407     |-390.970               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |317      |2% of pixels = 623 pixels (2.0%)   |226.084    |15.360  |6.8%    |221.719    |230.449     |-267.061               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |317      |3% of pixels = 935 pixels (3.0%)   |295.178    |10.573  |3.6%    |292.173    |298.183     |-197.967               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |317      |Fixed 1,250 (4.0%)                 |342.912    |9.456   |2.8%    |340.225    |345.599     |-150.233               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |317      |Fixed 2,000 (6.4%)                 |409.264    |10.462  |2.6%    |406.290    |412.237     |-83.881                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |317      |Fixed 3,000 (9.6%)                 |458.222    |9.094   |2.0%    |455.638    |460.807     |-34.923                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |317      |Fixed 4,000 (12.8%)                |493.145    |8.461   |1.7%    |490.740    |495.550     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |405      |1% of pixels = 364 pixels (1.0%)   |104.579    |9.727   |9.3%    |101.815    |107.343     |-547.567               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |405      |2% of pixels = 727 pixels (2.0%)   |215.882    |14.046  |6.5%    |211.890    |219.874     |-436.264               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |405      |3% of pixels = 1,091 pixels (3.0%) |296.290    |13.383  |4.5%    |292.487    |300.094     |-355.856               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |405      |Fixed 1,250 (3.4%)                 |330.167    |12.037  |3.6%    |326.747    |333.588     |-321.979               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |405      |Fixed 2,000 (5.5%)                 |446.205    |14.777  |3.3%    |442.006    |450.405     |-205.941               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |405      |Fixed 3,000 (8.2%)                 |559.765    |17.918  |3.2%    |554.673    |564.857     |-92.381                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |405      |Fixed 4,000 (11.0%)                |652.146    |21.122  |3.2%    |646.144    |658.149     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |821      |1% of pixels = 292 pixels (1.0%)   |91.188     |10.733  |11.8%   |88.138     |94.238      |-414.094               |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |20m   |821      |2% of pixels = 584 pixels (2.0%)   |210.015    |9.299   |4.4%    |207.373    |212.658     |-295.266               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |821      |3% of pixels = 876 pixels (3.0%)   |278.097    |11.106  |4.0%    |274.940    |281.253     |-227.185               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |821      |Fixed 1,250 (4.3%)                 |338.072    |10.917  |3.2%    |334.969    |341.175     |-167.210               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |821      |Fixed 2,000 (6.9%)                 |413.415    |12.853  |3.1%    |409.762    |417.068     |-91.867                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |821      |Fixed 3,000 (10.3%)                |472.813    |10.233  |2.2%    |469.905    |475.721     |-32.469                |all_pixels; fallback_sampled_pixels                                       |
|Alpha-Hull Area   |20m   |821      |Fixed 4,000 (13.7%)                |505.282    |9.942   |2.0%    |502.456    |508.107     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |905      |1% of pixels = 270 pixels (1.0%)   |65.236     |9.517   |14.6%   |62.532     |67.941      |-509.240               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |905      |2% of pixels = 540 pixels (2.0%)   |187.446    |14.427  |7.7%    |183.346    |191.546     |-387.030               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |905      |3% of pixels = 811 pixels (3.0%)   |266.511    |14.083  |5.3%    |262.509    |270.513     |-307.965               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |905      |Fixed 1,250 (4.6%)                 |351.976    |11.343  |3.2%    |348.753    |355.200     |-222.500               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |905      |Fixed 2,000 (7.4%)                 |437.436    |11.867  |2.7%    |434.063    |440.809     |-137.040               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |905      |Fixed 3,000 (11.1%)                |515.928    |17.562  |3.4%    |510.937    |520.920     |-58.548                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |905      |Fixed 4,000 (14.8%)                |574.476    |12.287  |2.1%    |570.984    |577.968     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |20m   |912      |1% of pixels = 296 pixels (1.0%)   |74.050     |9.119   |12.3%   |71.458     |76.642      |-534.834               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |912      |2% of pixels = 592 pixels (2.0%)   |209.368    |12.208  |5.8%    |205.898    |212.837     |-399.517               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |912      |3% of pixels = 889 pixels (3.0%)   |305.005    |11.447  |3.8%    |301.751    |308.258     |-303.879               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |912      |Fixed 1,250 (4.2%)                 |379.228    |12.975  |3.4%    |375.540    |382.915     |-229.656               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |912      |Fixed 2,000 (6.8%)                 |477.409    |14.081  |2.9%    |473.408    |481.411     |-131.475               |all_pixels                                                                |
|Alpha-Hull Area   |20m   |912      |Fixed 3,000 (10.1%)                |560.362    |12.213  |2.2%    |556.891    |563.833     |-48.522                |all_pixels                                                                |
|Alpha-Hull Area   |20m   |912      |Fixed 4,000 (13.5%)                |608.884    |13.301  |2.2%    |605.104    |612.664     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_15 |Fixed 1,250 (0.7%)                 |313.928    |12.425  |4.0%    |310.397    |317.459     |-176.076               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_15 |1% of pixels = 1,728 pixels (1.0%) |357.812    |11.524  |3.2%    |354.537    |361.087     |-132.192               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_15 |2% of pixels = 3,456 pixels (2.0%) |467.794    |13.439  |2.9%    |463.974    |471.613     |-22.210                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_15 |Fixed 4,000 (2.3%)                 |490.004    |10.855  |2.2%    |486.919    |493.089     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_15 |3% of pixels = 5,000 pixels (2.9%) |523.523    |10.670  |2.0%    |520.491    |526.556     |33.519                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_15 |Fixed 6,000 (3.5%)                 |553.926    |11.283  |2.0%    |550.720    |557.133     |63.922                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_15 |Fixed 8,000 (4.6%)                 |600.689    |11.350  |1.9%    |597.464    |603.915     |110.685                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_24 |Fixed 1,250 (0.7%)                 |328.010    |13.103  |4.0%    |324.286    |331.734     |-170.952               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_24 |1% of pixels = 1,668 pixels (1.0%) |371.023    |11.257  |3.0%    |367.824    |374.222     |-127.939               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_24 |2% of pixels = 3,336 pixels (2.0%) |473.805    |12.884  |2.7%    |470.144    |477.467     |-25.157                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_24 |Fixed 4,000 (2.4%)                 |498.962    |10.787  |2.2%    |495.896    |502.028     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_24 |3% of pixels = 5,000 pixels (3.0%) |536.865    |11.527  |2.1%    |533.589    |540.141     |37.903                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_24 |Fixed 6,000 (3.6%)                 |563.545    |12.535  |2.2%    |559.982    |567.108     |64.583                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_24 |Fixed 8,000 (4.8%)                 |611.805    |12.646  |2.1%    |608.211    |615.399     |112.843                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_36 |Fixed 1,250 (0.7%)                 |387.419    |13.508  |3.5%    |383.580    |391.258     |-229.196               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_36 |1% of pixels = 1,742 pixels (1.0%) |454.098    |11.499  |2.5%    |450.830    |457.366     |-162.517               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_36 |2% of pixels = 3,484 pixels (2.0%) |593.861    |11.844  |2.0%    |590.495    |597.226     |-22.754                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_36 |Fixed 4,000 (2.3%)                 |616.615    |11.241  |1.8%    |613.420    |619.809     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_36 |3% of pixels = 5,000 pixels (2.9%) |658.706    |13.355  |2.0%    |654.911    |662.502     |42.091                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_36 |Fixed 6,000 (3.4%)                 |690.741    |12.935  |1.9%    |687.064    |694.417     |74.126                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_36 |Fixed 8,000 (4.6%)                 |744.801    |14.839  |2.0%    |740.584    |749.019     |128.187                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_49 |Fixed 1,250 (0.6%)                 |347.955    |13.542  |3.9%    |344.106    |351.803     |-247.191               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_49 |1% of pixels = 2,070 pixels (1.0%) |457.189    |14.171  |3.1%    |453.162    |461.217     |-137.957               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_49 |Fixed 4,000 (1.9%)                 |595.146    |14.678  |2.5%    |590.975    |599.317     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_49 |2% of pixels = 4,139 pixels (2.0%) |602.564    |13.477  |2.2%    |598.734    |606.394     |7.418                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_49 |3% of pixels = 5,000 pixels (2.4%) |644.718    |11.480  |1.8%    |641.455    |647.980     |49.572                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_49 |Fixed 6,000 (2.9%)                 |682.938    |13.918  |2.0%    |678.983    |686.894     |87.793                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_49 |Fixed 8,000 (3.9%)                 |742.940    |13.029  |1.8%    |739.237    |746.642     |147.794                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_51 |Fixed 1,250 (0.6%)                 |348.930    |14.125  |4.0%    |344.915    |352.944     |-233.002               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_51 |1% of pixels = 2,006 pixels (1.0%) |445.863    |11.310  |2.5%    |442.648    |449.077     |-136.069               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_51 |Fixed 4,000 (2.0%)                 |581.931    |12.554  |2.2%    |578.364    |585.499     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_51 |2% of pixels = 4,013 pixels (2.0%) |583.295    |13.148  |2.3%    |579.558    |587.031     |1.363                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_51 |3% of pixels = 5,000 pixels (2.5%) |629.186    |13.847  |2.2%    |625.250    |633.121     |47.254                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_51 |Fixed 6,000 (3.0%)                 |668.303    |14.317  |2.1%    |664.234    |672.372     |86.372                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_51 |Fixed 8,000 (4.0%)                 |725.997    |13.623  |1.9%    |722.125    |729.868     |144.065                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_56 |Fixed 1,250 (0.9%)                 |344.798    |10.980  |3.2%    |341.678    |347.919     |-189.853               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_56 |1% of pixels = 1,455 pixels (1.0%) |370.377    |10.982  |3.0%    |367.256    |373.499     |-164.273               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_56 |2% of pixels = 2,910 pixels (2.0%) |484.616    |13.078  |2.7%    |480.899    |488.333     |-50.035                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_56 |Fixed 4,000 (2.7%)                 |534.651    |12.436  |2.3%    |531.117    |538.185     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_56 |3% of pixels = 4,365 pixels (3.0%) |549.934    |11.182  |2.0%    |546.756    |553.112     |15.283                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_56 |Fixed 6,000 (4.1%)                 |597.378    |13.480  |2.3%    |593.547    |601.209     |62.727                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_56 |Fixed 8,000 (5.5%)                 |641.158    |14.444  |2.3%    |637.053    |645.263     |106.507                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_60 |Fixed 1,250 (0.8%)                 |378.921    |11.288  |3.0%    |375.713    |382.129     |-225.018               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_60 |1% of pixels = 1,645 pixels (1.0%) |433.571    |14.336  |3.3%    |429.497    |437.646     |-170.368               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_60 |2% of pixels = 3,290 pixels (2.0%) |569.345    |11.830  |2.1%    |565.983    |572.707     |-34.594                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_60 |Fixed 4,000 (2.4%)                 |603.939    |15.748  |2.6%    |599.464    |608.414     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_60 |3% of pixels = 4,935 pixels (3.0%) |642.888    |12.418  |1.9%    |639.358    |646.417     |38.949                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_60 |Fixed 6,000 (3.6%)                 |678.660    |13.068  |1.9%    |674.946    |682.374     |74.721                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_60 |Fixed 8,000 (4.9%)                 |730.859    |12.554  |1.7%    |727.291    |734.426     |126.920                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_8  |Fixed 1,250 (0.7%)                 |348.411    |12.027  |3.5%    |344.993    |351.829     |-224.738               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_8  |1% of pixels = 1,879 pixels (1.0%) |429.608    |12.752  |3.0%    |425.983    |433.232     |-143.542               |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_8  |2% of pixels = 3,758 pixels (2.0%) |560.364    |15.768  |2.8%    |555.883    |564.846     |-12.785                |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_8  |Fixed 4,000 (2.1%)                 |573.149    |12.491  |2.2%    |569.600    |576.699     |0.000                  |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_8  |3% of pixels = 5,000 pixels (2.7%) |616.684    |13.007  |2.1%    |612.987    |620.381     |43.535                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_8  |Fixed 6,000 (3.2%)                 |655.308    |13.235  |2.0%    |651.546    |659.069     |82.159                 |all_pixels                                                                |
|Alpha-Hull Area   |50m   |sub50_8  |Fixed 8,000 (4.3%)                 |707.509    |13.447  |1.9%    |703.687    |711.330     |134.359                |all_pixels                                                                |

## Figures

Mean-by-sample-size figures include 95 percent confidence interval error bars around the 50-iteration mean. Compact overview distribution plots are retained for each metric. Longer distribution plots are split by sample-size rule, with quadrats as boxes and retained pixel counts in the x-axis labels.

## PCA Mean Distance

![PCA Mean Distance mean](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_mean_by_sample_size.png)

![PCA Mean Distance CV](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_cv_by_sample_size.png)

![PCA Mean Distance delta](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_delta_from_fixed_4000.png)

![PCA Mean Distance replicate distributions](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_replicate_distributions.png)

### PCA Mean Distance: 10m Scale Figures

![PCA Mean Distance 10m mean](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_mean_by_sample_size_10m.png)

![PCA Mean Distance 10m CV](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_cv_by_sample_size_10m.png)

![PCA Mean Distance 10m delta](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_delta_from_fixed_4000_10m.png)

### PCA Mean Distance: 10m Distribution Charts Split By Sample Size

- [10m 1%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_10m_pct_1.png)
- [10m 2%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_10m_pct_2.png)
- [10m 3%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_10m_pct_3.png)
- [10m 1,250](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_10m_fixed_1250.png)
- [10m 2,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_10m_fixed_2000.png)
- [10m 3,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_10m_fixed_3000.png)
- [10m 4,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_10m_fixed_4000.png)

### PCA Mean Distance: 20m Scale Figures

![PCA Mean Distance 20m mean](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_mean_by_sample_size_20m.png)

![PCA Mean Distance 20m CV](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_cv_by_sample_size_20m.png)

![PCA Mean Distance 20m delta](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_delta_from_fixed_4000_20m.png)

### PCA Mean Distance: 20m Distribution Charts Split By Sample Size

- [20m 1%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_20m_pct_1.png)
- [20m 2%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_20m_pct_2.png)
- [20m 3%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_20m_pct_3.png)
- [20m 1,250](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_20m_fixed_1250.png)
- [20m 2,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_20m_fixed_2000.png)
- [20m 3,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_20m_fixed_3000.png)
- [20m 4,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_20m_fixed_4000.png)

### PCA Mean Distance: 50m Scale Figures

![PCA Mean Distance 50m mean](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_mean_by_sample_size_50m.png)

![PCA Mean Distance 50m CV](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_cv_by_sample_size_50m.png)

![PCA Mean Distance 50m delta](../figures/sample_size_effects/pca_mean_distance/pca_mean_distance_delta_from_fixed_4000_50m.png)

### PCA Mean Distance: 50m Distribution Charts Split By Sample Size

- [50m 1%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_50m_pct_1.png)
- [50m 2%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_50m_pct_2.png)
- [50m 3%](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_50m_pct_3.png)
- [50m 1,250](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_50m_fixed_1250.png)
- [50m 4,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_50m_fixed_4000.png)
- [50m 6,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_50m_fixed_6000.png)
- [50m 8,000](../figures/sample_size_effects/pca_mean_distance/distributions_by_sample_size/pca_mean_distance_replicate_distribution_50m_fixed_8000.png)


Output tables:

- `reports/tables/sample_size_effects/pca_mean_distance/pca_mean_distance_sample_size_design.csv`
- `reports/tables/sample_size_effects/pca_mean_distance/pca_mean_distance_sample_size_boot_long.csv`
- `reports/tables/sample_size_effects/pca_mean_distance/pca_mean_distance_sample_size_summary.csv`

## Spectral Rao's Q

![Spectral Rao's Q mean](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_mean_by_sample_size.png)

![Spectral Rao's Q CV](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_cv_by_sample_size.png)

![Spectral Rao's Q delta](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_delta_from_fixed_4000.png)

![Spectral Rao's Q replicate distributions](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_replicate_distributions.png)

### Spectral Rao's Q: 10m Scale Figures

![Spectral Rao's Q 10m mean](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_mean_by_sample_size_10m.png)

![Spectral Rao's Q 10m CV](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_cv_by_sample_size_10m.png)

![Spectral Rao's Q 10m delta](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_delta_from_fixed_4000_10m.png)

### Spectral Rao's Q: 10m Distribution Charts Split By Sample Size

- [10m 1%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_10m_pct_1.png)
- [10m 2%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_10m_pct_2.png)
- [10m 3%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_10m_pct_3.png)
- [10m 1,250](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_10m_fixed_1250.png)
- [10m 2,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_10m_fixed_2000.png)
- [10m 3,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_10m_fixed_3000.png)
- [10m 4,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_10m_fixed_4000.png)

### Spectral Rao's Q: 20m Scale Figures

![Spectral Rao's Q 20m mean](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_mean_by_sample_size_20m.png)

![Spectral Rao's Q 20m CV](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_cv_by_sample_size_20m.png)

![Spectral Rao's Q 20m delta](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_delta_from_fixed_4000_20m.png)

### Spectral Rao's Q: 20m Distribution Charts Split By Sample Size

- [20m 1%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_20m_pct_1.png)
- [20m 2%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_20m_pct_2.png)
- [20m 3%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_20m_pct_3.png)
- [20m 1,250](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_20m_fixed_1250.png)
- [20m 2,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_20m_fixed_2000.png)
- [20m 3,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_20m_fixed_3000.png)
- [20m 4,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_20m_fixed_4000.png)

### Spectral Rao's Q: 50m Scale Figures

![Spectral Rao's Q 50m mean](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_mean_by_sample_size_50m.png)

![Spectral Rao's Q 50m CV](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_cv_by_sample_size_50m.png)

![Spectral Rao's Q 50m delta](../figures/sample_size_effects/spectral_rao_q/spectral_rao_q_delta_from_fixed_4000_50m.png)

### Spectral Rao's Q: 50m Distribution Charts Split By Sample Size

- [50m 1%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_50m_pct_1.png)
- [50m 2%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_50m_pct_2.png)
- [50m 3%](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_50m_pct_3.png)
- [50m 1,250](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_50m_fixed_1250.png)
- [50m 4,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_50m_fixed_4000.png)
- [50m 6,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_50m_fixed_6000.png)
- [50m 8,000](../figures/sample_size_effects/spectral_rao_q/distributions_by_sample_size/spectral_rao_q_replicate_distribution_50m_fixed_8000.png)


Output tables:

- `reports/tables/sample_size_effects/spectral_rao_q/spectral_rao_q_sample_size_design.csv`
- `reports/tables/sample_size_effects/spectral_rao_q/spectral_rao_q_sample_size_boot_long.csv`
- `reports/tables/sample_size_effects/spectral_rao_q/spectral_rao_q_sample_size_summary.csv`

## Alpha-Hull Area

![Alpha-Hull Area mean](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_mean_by_sample_size.png)

![Alpha-Hull Area CV](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_cv_by_sample_size.png)

![Alpha-Hull Area delta](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_delta_from_fixed_4000.png)

![Alpha-Hull Area replicate distributions](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_replicate_distributions.png)

### Alpha-Hull Area: 10m Scale Figures

![Alpha-Hull Area 10m mean](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_mean_by_sample_size_10m.png)

![Alpha-Hull Area 10m CV](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_cv_by_sample_size_10m.png)

![Alpha-Hull Area 10m delta](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_delta_from_fixed_4000_10m.png)

### Alpha-Hull Area: 10m Distribution Charts Split By Sample Size

- [10m 1%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_10m_pct_1.png)
- [10m 2%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_10m_pct_2.png)
- [10m 3%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_10m_pct_3.png)
- [10m 1,250](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_10m_fixed_1250.png)
- [10m 2,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_10m_fixed_2000.png)
- [10m 3,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_10m_fixed_3000.png)
- [10m 4,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_10m_fixed_4000.png)

### Alpha-Hull Area: 20m Scale Figures

![Alpha-Hull Area 20m mean](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_mean_by_sample_size_20m.png)

![Alpha-Hull Area 20m CV](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_cv_by_sample_size_20m.png)

![Alpha-Hull Area 20m delta](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_delta_from_fixed_4000_20m.png)

### Alpha-Hull Area: 20m Distribution Charts Split By Sample Size

- [20m 1%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_20m_pct_1.png)
- [20m 2%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_20m_pct_2.png)
- [20m 3%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_20m_pct_3.png)
- [20m 1,250](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_20m_fixed_1250.png)
- [20m 2,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_20m_fixed_2000.png)
- [20m 3,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_20m_fixed_3000.png)
- [20m 4,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_20m_fixed_4000.png)

### Alpha-Hull Area: 50m Scale Figures

![Alpha-Hull Area 50m mean](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_mean_by_sample_size_50m.png)

![Alpha-Hull Area 50m CV](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_cv_by_sample_size_50m.png)

![Alpha-Hull Area 50m delta](../figures/sample_size_effects/alpha_hull_area/alpha_hull_area_delta_from_fixed_4000_50m.png)

### Alpha-Hull Area: 50m Distribution Charts Split By Sample Size

- [50m 1%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_50m_pct_1.png)
- [50m 2%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_50m_pct_2.png)
- [50m 3%](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_50m_pct_3.png)
- [50m 1,250](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_50m_fixed_1250.png)
- [50m 4,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_50m_fixed_4000.png)
- [50m 6,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_50m_fixed_6000.png)
- [50m 8,000](../figures/sample_size_effects/alpha_hull_area/distributions_by_sample_size/alpha_hull_area_replicate_distribution_50m_fixed_8000.png)


Output tables:

- `reports/tables/sample_size_effects/alpha_hull_area/alpha_hull_area_sample_size_design.csv`
- `reports/tables/sample_size_effects/alpha_hull_area/alpha_hull_area_sample_size_boot_long.csv`
- `reports/tables/sample_size_effects/alpha_hull_area/alpha_hull_area_sample_size_summary.csv`

