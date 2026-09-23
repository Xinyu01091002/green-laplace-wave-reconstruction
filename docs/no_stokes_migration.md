# Pure Green--Laplace migration, 2026-09-23

This repository now maintains the selected Green--Laplace implementations
from the retired research repository. The migration excludes the historical
Stokes-corrected executors and the order-four angular/compact repairs.

## Runtime changes

The earlier public third-order executor applied pair-level eta/Phi repairs
and third-order elevation/surface-potential diagonal repairs. Its replacement,
`src/internal/gl_no_stokes_eta33.m`, contains none of those branches or helper
functions. Both the existing public API and the paper rank-study entry point
call this implementation. There is no option to enable the old repairs.

The third-order extraction preserves the selected source's no-repair path.
The fourth-order source is copied from the author's already uncorrected
implementation, with only repository paths and environment-variable names
adapted. `gl_pure_sum_order4` exposes its analytic dimensionless elevation and
true surface-potential outputs as an experimental interface. Required lower
states are generated internally from eta11.

The separate Two-Scale/Shared-Scale implementation and its ordered-triple
references are preserved under `research/two_scale`. Existing eta20
difference-frequency diagnostics retain their separate diagnostic status.

Surface Taylor terms are retained: evaluating phi at z=eta is necessary to
obtain surface potential and is not a Stokes-diagonal repair. Exact Stokes
traces and historical symbolic certification sources already present in the
repository remain mathematical references, not runtime correction functions.
The runtime third-order JSON is a metadata projection that omits the old
correction expressions; no governing equation was derived or retuned during
this migration. The historical generator can regenerate its full mathematical
record, but the current MATLAB runtime never consumes its correction fields.

## Verification and limits

MATLAB R2022b, double precision, raw same-unit comparisons:

- Third-order migration parity on a 64 by 64 grid at ranks 4, 6 and 8:
  elevation and surface-potential relative L2 differences were zero against
  the source evaluator with every repair disabled.
- Fourth-order migration parity on the same input at rank 4: both field
  differences were zero against the author's no-Stokes source.
- Two-Scale/Shared-Scale ordered-triple parity uses the independent MATLAB
  ordered reconstruction. These checks test both elevation and surface
  potential, including the mixed inner/outer surface route.
- The release suite includes a runtime source guard and a fourth-order smoke
  check. The fourth-order forcing graph was regenerated and exactly certified
  by Wolfram against the original Euler residual.

These are migration, algebraic and implementation checks. They do not establish
new broad-spectrum accuracy or measured speedups. Earlier corrected-version
accuracy, exact-diagonal and operation-count claims must not be reused for this
version. Report costs from the returned audit for the requested output.

The public total-field API keeps its previous strict-forward, pure-sum and
alias-safe domain restrictions. Mixed-sign, resonant and arbitrary-order
production support are not introduced by this migration.

## Provenance and retained material

`no_stokes_source_manifest.json` records hashes of the source files, including
uncommitted author changes. It supplements and, for replaced files, supersedes
the earlier extraction manifest. The source worktrees and their research
artifacts are preserved in a separate local archive; corrected variants,
machine-specific launchers, incomplete experiments and large field archives
are not copied into this public runtime.

Manuscript working drafts remain in the preserved manuscript worktree. This
migration does not replace the paper with a stale research snapshot. Committed
figure references remain version-specific historical assets; regenerate
figures with the current entry points before using them as current results.
