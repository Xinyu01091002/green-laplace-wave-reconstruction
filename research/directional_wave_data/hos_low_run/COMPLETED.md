# Completed directional focused-wavegroup comparison

Verified 2026-09-26. Both four-phase families reached 220 s, with 1101 finite
probe records per phase and successful completion checks. The low family used
8 MPI ranks per phase (32 total); the earlier high family used 4 per phase.

Low HOS: 13:10:34--14:32:43 UTC, batch wall time 4929 s (1:22:09).
Peak sampled aggregate RSS: 4457388 KiB (4.251 GiB).
High HOS: 04:12:14--06:45:33 UTC, batch wall time 9199 s (2:33:19).
Different process counts and background loads preclude a steepness-only timing
comparison. Full GL processing of both families finished at 14:34:28 UTC.

The reference is the four-phase second-harmonic sector, not an exact isolated
perturbation order. GL uses the separated first-harmonic record jointly with
the initial complex directional spectrum; directional allocation is assumed,
not uniquely inferred from one probe. No gain, phase or shift fitting is used.

J8, 7.5-degree directional allocation, identical sum-band filtering of reference
and prediction, main-group window around the first-harmonic envelope peak:

| y offset (m) | Akp .02 relative L2 (%) | Akp .12 relative L2 (%) |
| --- | ---: | ---: |
| 0 | 0.267241 | 4.214164 |
| -52.7823 | 0.338830 | 3.088965 |
| +52.7823 | 0.350890 | 3.131221 |
| -70.3762 | 0.454777 | 2.572421 |
| +70.3762 | 0.464956 | 2.614949 |

Full raw/common-band metrics are saved as focused_akp002_metrics.csv and
focused_akp012_metrics.csv alongside this note. Plots and JSON reports are in
artifacts/hos_ocean/directional-full/{akp002,akp012}/ locally. Remote outputs
are full-gl-comparison-20260926-v1 under each original HOS run root named in
README.md. This completion claim concerns second-harmonic reconstruction;
no new directional third-harmonic result is claimed.

## Next random-wave trial: assessment, not launched

Use one fixed realization (existing seed 20260925), four global phase shifts
of that same realization, and identical random phases between amplitude levels.
Existing prepare_compact_wavegroup.m supports phase-only randomization and
fresh MF12 order-2 eta/true-surface-psi construction. Do not reuse the historical
tapered unit-RMS prototype or independently randomize the four phase runs.

The prior user-selected amplitude rule preserves wavegroup modal magnitudes.
Its high focusing label Akp=.12 has linear Hs=0.4076151725 m and kp*Hs/2=
0.005686231656; the .02 counterpart gives Hs=0.0679358621 m and kp*Hs/2=
0.000947705276. Both are weak random seas. Whether the new request instead
means sea-state steepness .12/.02 has been asked explicitly; no rescaling or
new production run is performed before that distinction is resolved.

The focused-wavegroup postprocessor needs adaptation: use a declared interior
record window rather than selecting a strongest packet; report finite-record
edge/window sensitivity and the conditioning of complex directional allocation.
Apply any temporal window or frequency filter identically to both compared
signals. A 220 s single-seed trial is a practical first comparison, not converged
random-sea statistics. Keep initial-transient effects distinct from GL errors.

Remote read at 15:08 UTC: 48 visible CPUs, 8 OW3D processes, zero HOS processes,
load 8.22 and 605 GiB available RAM. Reuse four concurrent phases x 8 MPI ranks
for one amplitude at a time; measure wall time and aggregate memory anew.
No random-run timing claim is available yet.
