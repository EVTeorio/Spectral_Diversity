# All-pixel alpha = 0.5 hulls at 10 m, 20 m, and 25 m

Date: 2026-10-15

## Outcome and coverage

Completed the authorized production passes at 10 m, then 20 m, then 25 m in
original and vector-normalized global PCA PC1-PC2 space. Stopped before 50 m.
No 33 m or 100 m processing was performed. No new exclusion flags were added;
25 m exclusion decisions remain with the user.

Counts below use available/non-excluded spectral quadrats. Existing reference
exclusions remain at 10 m and 20 m; no exclusions are applied at 25 m.

| Scale | Expected quadrats per basis | Original finite areas | Vector-normalized finite areas | Production-pass elapsed minutes |
|---|---:|---:|---:|---:|
| 10 m | 1,744 | 1,743* | 1,744 | 30.5 |
| 20 m | 436 | 436* | 436 | 19.7 |
| 25 m | 312 | 312 | 312 | 48.5 |

*Includes subsequent user-authorized manual orientation repairs: 10 m 1402_b,
1019_a, and 1816_c; 20 m 413.

Pass timings exclude subsequent recovery and validation. The 10 m elapsed time
includes the checkpointed scheduling restart. Eight local R workers were used.

Canonical tables retain every boundary ID and both bases: 4,000 rows at 10 m,
1,000 at 20 m, and 640 at 25 m. Missing-raster cells number 91, 15, and 8,
respectively, per basis. Existing manually excluded reference cells number 165
at 10 m and 49 at 20 m, per basis. These statuses are separate from geometry
failures and do not receive invented values.

## Remaining geometry failures

One original-PCA area remains NA with an explicit failure status:

| Scale | Quadrat | Reason |
|---|---|---|
| 10 m | `1000_c` | Hull construction fails inside the package's alpha-shape routine. |

`1000_c` also failed in the historical original-PCA alpha = 1 product. Three exact
point reorderings did not resolve its construction error. Reordering the points
for `1402_b` did not resolve its closure failure before manual review. The user
deferred further inspection of `1000_c`; its raster, checkpoint, and NA value
are unchanged. Manual orientation repairs for the other three 10 m cases and
20 m 413 are documented below. No alpha changes or manual jitter were applied.

## Methods and recovery

Used the saved original and vector-normalized PCA bases, the established 563 nm
shadow threshold of `0.0305476`, current spectrum cleaning and vector normalization,
and all retained valid pixels. Hull geometry uses all unique PC1-PC2 coordinates
with R `alphahull::ahull(alpha = 0.5)`. Successful existing alpha = 0.5 values were
reused only when keys and retained-pixel counts matched. Older alpha = 1 fields
were preserved.

The default area routine initially left 18 geometry/area failures. Thirteen were
recovered by exact line integration over the same returned circular-arc boundary:
10 through the center/arc representation and three through an equivalent endpoint/
chord-and-circular-segment representation. The latter enforces exact oriented
vertex-ID balance and verifies endpoint-to-circle consistency. Both use Green's
theorem; they do not sample the boundary or change alpha. Recovered values must
be finite, nonnegative, and no larger than the convex hull of the same points.
Existing finite values were preserved, first attempts were backed up, and each
recovery has explicit source provenance and an audit record.

Each integration form was checked against the package on 100 synthetic hulls
covering clouds, disconnected components, holes, and translations. Maximum
absolute discrepancies were approximately `2.12e-12` for the center form and
`7.38e-14` for the endpoint form. Actual successful hulls in both PCA bases were
also compared; detailed validation tables are in the dated run directory.

The installed triangulator can automatically retry with jitter during numerical
failure handling. Its warnings are retained in the output `warning` field.
Successful rows carrying such warnings are original-PCA 10 m `413_d`, and both
bases at 25 m `sub25_69`. These are available for the user's inclusion review;
no manual exclusion flags were added. No manual coordinate perturbation was
introduced by these scripts.

## Validation and outputs

Full-grid ID coverage, unique scale/quadrat/basis/alpha keys, positive-area bounds,
CSV round trips, unchanged prior finite areas, and protected reference hashes
were checked. Both PCA basis objects, established spectral summaries, the 25 m
distance/Rao products, and the historical hull table containing 50 m results
remain unchanged. No manuscript assets were modified. The R completion beep
was invoked.

Canonical alpha = 0.5 long-format tables:

- [10 m](../../Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_10m_all_pixels.csv)
- [20 m](../../Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_20m_all_pixels.csv)
- [25 m](../../Quad_Values/Spectral_diversitySHPs/alpha_hull_0_5_25m_all_pixels.csv)

Audit and reproducibility:

- [Final coverage](../tables/20261015_alpha_half_small_scales/final_metric_coverage.csv)
- [Unresolved cases](../tables/20261015_alpha_half_small_scales/unresolved_geometry.csv)
- [Recovery audit, including repeated attempts](../tables/20261015_alpha_half_small_scales/area_recovery_audit.csv)
- [Post-recovery validation](../tables/20261015_alpha_half_small_scales/post_recovery_validation.csv)
- [Independent saved-file validation](../tables/20261015_alpha_half_small_scales/final_independent_validation.csv)
- [Observed pass timings](../tables/20261015_alpha_half_small_scales/overall_scale_timings.csv)
- [R session](../tables/20261015_alpha_half_small_scales/sessionInfo.txt)
- [Production script](../../scripts/3_Analysis/alpha_half_small_scales.R)
- [Exact area functions](../../scripts/3_Analysis/alpha_hull_arc_integral.R)
- [Recovery script](../../scripts/3_Analysis/recover_alpha_half_area_evaluations.R)

The existing combined-table alpha fields continue to represent the earlier
alpha = 1 workflow. Current alpha = 0.5 analyses should join the new explicit
alpha tables using scale, quadrat ID, basis, and alpha value.

### Targeted retry of 20 m quadrat 413

On 2026-10-15 at 16:16 CDT, original PCA quadrat 413 was rerun alone at
alpha = 0.5. The fresh run retained 28,511 pixels (28,180 unique PC1-PC2
points) and returned 561 arcs in 13.98 seconds. The package area remained NA;
both exact integration checks failed closure, with the endpoint check again
finding two unbalanced boundary vertices. The canonical result remains NA
and its file hash was verified unchanged. No other quadrats were rerun.

The [retry diagnostics](../tables/20261015_alpha_half_small_scales/retry_20m_413_original_20261015_161612/retry_result.csv),
previous result, and returned hull are saved together in the retry directory.

Further inspection found two incident arcs at each flagged vertex: the problem
is inconsistent orientation in a nearly degenerate three-arc loop, rather than
absent endpoint IDs. The flagged vertices are about 0.08630616 PCA units apart,
so merging them would alter the geometry. Three unchanged-point ordering tests
were then attempted with the same algorithm and alpha. Lexicographic ordering
returned a package area of approximately 604.1341, but failed both independent
closure checks; it was not accepted. Reverse input and reverse lexicographic
orders failed hull construction, with package-internal jitter retry warnings.
The [ordering diagnostics](../tables/20261015_alpha_half_small_scales/retry_20m_413_original_20261015_161612/point_order_diagnostics.csv)
and constructed hull are retained. The canonical NA remains unchanged.

A [diagnostic figure](../figures/20261015_alpha_half_20m_413.png)
shows the saved full hull and a magnified view of the problematic arcs;
a [vector PDF](../figures/20261015_alpha_half_20m_413.pdf) is also available.
It plots the original retry geometry without repairs or recomputation, using
[this script](../../scripts/3_Analysis/plot_alpha_half_20m_413.R).

### User-authorized manual orientation repair of 413

After inspecting the figure, the user authorized manual closure of this case.
All three problematic arcs form a separate near-degenerate component with
vertex IDs 4324, 4325, and 4535. Every vertex has two incident arcs, but arc 423
has inconsistent traversal relative to the other two. Reversing its traversal
balances all directed vertex counts without moving any coordinates, changing
any circular curves, changing alpha, or removing points. The original `ahull`
object remains intact; the corrected traversal is stored separately. This is a
manual interpretation, not an unmodified package result.

The corrected total area is **460.7109822782754 PCA-score squared units**.
Endpoint/chord and independent center/arc integration agree within 4e-13;
the reconstructed boundary closure error is 8.37e-14. The isolated loop's signed
area is only 3.928773328883682e-11. For comparison, the closed remainder alone
has area 460.7109822782361. The repaired area is finite, nonnegative, and below
the saved convex bound of 848.118450530432.

The canonical 20 m original-PCA 413 row and its checkpoint now contain the
repaired value, with source `computed_all_pixels_alpha_0_5_manual_arc_orientation_repair`
and an explicit manual-repair warning. All other canonical rows and protected
references were verified unchanged. Finite original-PCA 20 m coverage is now
436/436; four unresolved original-PCA cases remain at 10 m. Earlier retry notes
and validation tables above describe their historical pre-repair state.

- [Manual repair audit](../tables/20261015_alpha_half_small_scales/manual_orientation_20m_413/manual_repair_audit.csv)
- [Corrected signed traversal](../tables/20261015_alpha_half_small_scales/manual_orientation_20m_413/repaired_traversal.csv)
- [Focused validation](../tables/20261015_alpha_half_small_scales/manual_orientation_20m_413/validation.csv)
- [Repair script](../../scripts/3_Analysis/manual_orientation_20m_413.R)

Original failed rows and affected-table snapshots are retained in the manual
repair directory. No other quadrats were repaired under this exception.

### Focused diagnosis of the four unresolved 10 m cases

Fresh reconstruction on October 15 confirmed that all four rasters can be read,
cleaned, and projected into the saved original PCA basis. Failure occurs during
geometry construction or area evaluation, not because their spectra are absent.
No canonical metric values were changed by this diagnostic pass.

| Quadrat | Diagnosis of freshly returned geometry |
|---|---|
| `1000_c` | Construction again throws `missing value where TRUE/FALSE needed`; no hull is returned. 7,126 retained pixels. |
| `1402_b` | Two directed vertex imbalances in an isolated three-arc/three-vertex component, spanning approximately 0.326 by 0.217 PCA units. Reversing arc 353 would balance vertex counts. |
| `1019_a` | Two directed vertex imbalances in a 292-arc/292-vertex component spanning approximately 36.52 by 16.44 PCA units. Reversing arc 487 would balance counts, but the affected component is a large boundary, not an isolated tiny loop. |
| `1816_c` | Two directed vertex imbalances in an isolated four-arc/four-vertex component, spanning approximately 0.120 by 0.048 PCA units. Reversing arc 545 would balance vertex counts. |

Single-arc reversal candidates are diagnostic only: count balance alone does
not validate the intended hull or authorize applying the 413 exception elsewhere.
Candidate signed component areas are approximately 0.01405, 330.0377, and
0.001389 PCA-score squared units for 1402_b, 1019_a, and 1816_c, respectively;
these are not accepted full-quadrat metric values. Saved hulls, affected-component
arc tables, and [diagnosis summary](../tables/20261015_alpha_half_small_scales/diagnose_four_10m/diagnosis.csv)
are retained for further inspection.

The [1019_a diagnostic figure](../figures/20261015_alpha_half_10m_1019_a.png)
shows the full saved hull and a close-up of arc 487 and its neighboring arcs.
It depicts the current inferred traversal directions without changing geometry
or stored values; [plotting script](../../scripts/3_Analysis/plot_alpha_half_10m_1019_a.R).

### Reviewed manual repairs of three 10 m quadrats

The user subsequently reviewed the arc changes and explicitly authorized
repairs of 1402_b, 1019_a, and 1816_c, deferring 1000_c. Using the saved
diagnostic hulls, reversed only the reviewed arc traversal in each case. All
input points, alpha = 0.5, circular centers, radii, and curves remain unchanged.
Original hulls and failed-table snapshots are preserved. Corrected signed
traversals are separate from the original `ahull` objects.

| Quadrat | Reversed saved arc | Repaired total area (PCA-score squared units) |
|---|---:|---:|
| `1402_b` | 353 | 339.4753276068346 |
| `1019_a` | 487 | 302.3550874568043 |
| `1816_c` | 545 | 290.6404426297782 |

All three now have exact directed vertex-ID balance. Reconstructed geometric
closure residuals are at most 6.63e-10; independent center/arc and endpoint/chord
areas differ by at most 2.78e-9 PCA-score squared units. Each total is finite,
nonnegative, and below the convex bound of its points. The 1019_a figure above
continues to depict the historical unrepaired traversal.

Canonical rows and per-quadrat checkpoints contain these values with source
`computed_all_pixels_alpha_0_5_manual_arc_orientation_repair` and explicit
manual-repair warnings. All other canonical rows, vector-normalized checkpoint
rows, 1000_c, 20 m/25 m tables, saved hulls, and protected references were
verified unchanged. Current original-PCA 10 m coverage is 1,743/1,744; the only
unresolved case is 1000_c. Earlier diagnosis and validation artifacts describe
their historical snapshots and are retained.

- [Three-repair audit](../tables/20261015_alpha_half_small_scales/manual_orientation_three_10m/manual_repair_audit.csv)
- [Focused validation](../tables/20261015_alpha_half_small_scales/manual_orientation_three_10m/validation.csv)
- [Repair script](../../scripts/3_Analysis/manual_orientation_three_10m.R)

Before-table backups and each corrected traversal CSV are in the same audit
directory. No construction or repair of 1000_c was attempted during this step.

### Subsequent targeted retry of 1000_c

At the user's request, retried only original-PCA 10 m 1000_c at alpha = 0.5
on October 15 at 16:53 CDT. All 7,126 retained spectra projected successfully,
giving the same 6,890 unique PC1-PC2 points. The run failed during hull
construction again with `missing value where TRUE/FALSE needed`; no hull was
returned for area integration or orientation repair. Captured the construction
call stack and failing conditional, which compares intermediate `mx1`/`mx2`
values. This identifies the failing operation, not its ultimate numerical cause.
Elapsed time was 7.2 seconds including input preparation. The saved canonical
table, 1000_c checkpoint, and original PCA object were hash-verified unchanged.
The original-PCA value remains NA; its vector-normalized value is unchanged.

- [Retry result](../tables/20261015_alpha_half_small_scales/retry_10m_1000_c_original_20261015_165350/retry_result.csv)
- [Construction call stack](../tables/20261015_alpha_half_small_scales/retry_10m_1000_c_original_20261015_165350/construction_error_calls.txt)

The [1000_c original-PCA point-cloud figure](../figures/20261015_original_pca_10m_1000_c.png)
shows all 7,126 retained pixels in PC1-PC2 space. No hull boundary or inferred
failure location is drawn because construction did not return a hull.
The [plotting script](../../scripts/3_Analysis/plot_original_pca_10m_1000_c.R)
uses the saved PCA basis and established spectrum filtering; metric values are unchanged.

Further tracing captured the actual failing intermediate row inside `ashape`:
candidate edge 1490 of 1542, linking unique point IDs 2402 and 2492. Both
intermediate coordinate pairs (`mx1,my1` and `mx2,my2`) are NaN, although the
input endpoints are finite: approximately (-1.9469131, 23.8768040) and
(-14.0583711, 15.7067545). The triangulation mesh contains eight nonfinite
entries. Thus the immediate failure is invalid derived geometry followed by an
unguarded comparison, not missing input PCA scores. The deeper reason those
derived coordinates become NaN has not yet been established. See the
[captured edge](../tables/20261015_alpha_half_small_scales/retry_10m_1000_c_original_20261015_165350/failing_construction_edge.csv).
