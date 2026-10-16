# Reviewed manual repairs of three 10 m original-PCA hulls

The user reviewed the arc changes and authorized repairs of 1402_b, 1019_a,
and 1816_c using the approach used for 20 m 413. Leave 1000_c unchanged.
Reviewed the saved diagnostic hulls and CSVs, manual_orientation_20m_413.R,
the canonical 10 m table, coverage summary, and project governance.

Reverse only the reviewed traversal for original saved arc IDs 353, 487, and
545, respectively. Keep input points, alpha = 0.5, and all circular curves
unchanged; preserve each original hull. Store signed traversal and explicit
manual source/warning metadata. Validate all three results before writing any
canonical result: directed vertex balance, circle consistency, independent
center/arc versus endpoint integrals, closure residual, and convex area bounds.
Do not require tiny component area for these reviewed cases: 1019_a involves
the large outer component, and 1402_b/1816_c have finite small loops.

Back up affected tables/checkpoints, update only the three authorized rows,
refresh current coverage/status/unresolved summaries, and verify all other rows,
1000_c, vector-normalized checkpoint values, 20/25 m tables, and protected
reference hashes remain unchanged. Update current documentation and report.
The main risk is manual interpretation of arc orientation; distinguish the
results from unmodified package outputs and retain original failed results.
