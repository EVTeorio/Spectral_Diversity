# Diagnose and repair 1000_c numerical geometry

User requests collinearity testing first, correction if implicated, otherwise
degenerate-triangle testing and correction. Scope is 10 m original-PCA 1000_c,
alpha = 0.5 only. Reviewed project governance, retry trace, failing edge,
installed alphahull delvor/ashape and interp circum/tri.mesh implementations.

Test global rank and local triangle determinants, area, heights, and edge
lengths. Identify which triangles produce nonfinite circumcenters and whether
they have collinear vertices or zero area. Attempt stable, translated-coordinate
circumcenter evaluation on the unchanged triangulation before any point removal,
jitter, or alpha changes. Preserve installed package code; any replacement
functions must live in a private evaluation environment or reconstructed mesh.

Validate improved circumcenters by equal radii, finite geometry and agreement
with standard evaluation on ordinary triangles. If a hull can be obtained,
require closure, independent area agreement, and convex-area bounds before
updating only the authorized canonical/checkpoint row and summaries. Preserve
original failures, protected references, prior repairs, and other scales.
If no mathematically justified unchanged-point fix succeeds, report the exact
remaining obstacle rather than fabricating a metric. Update documentation if
the numerical workflow or accepted dataset changes. No 50 m processing.
