# Manual orientation repair: 20 m original-PCA quadrat 413

User authorized a manual closure after reviewing the saved diagnostic figure.
Scope is only this quadrat/basis at alpha = 0.5; retain the original NA result.

Reviewed the saved retry hull, exact circular-arc integration helpers, canonical
20 m table, checkpoint, coverage/status summaries, and project governance.
The isolated three-arc component (293, 422, 423) has three vertices of degree
two but inconsistent traversal. Reversing arc 423 balances the component.
This changes orientation metadata, not points, alpha, centers, radii, or curves.

Save a separate signed traversal table and manual-repair audit. Validate exact
oriented vertex balance, radial endpoint consistency, agreement of center and
endpoint integrals, nonnegative total area below the saved convex bound, and
the negligible area of the small loop. Back up affected tables; update only the
413 original-PCA canonical/checkpoint row, coverage and failure summaries.
Verify all other canonical rows and protected reference hashes remain unchanged.
Update the current-state documentation and report with explicit manual provenance.

Risk: this is a manual interpretation of a numerical degeneracy, not an
unmodified package result. Preserve that distinction in source/warning fields.
Do not apply this exception automatically to other failures or recompute 50 m.
