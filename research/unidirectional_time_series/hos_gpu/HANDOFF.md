# GPU HOS and GL time-series handoff — 2026-09-27

## User's current request and the main correction

The user requested a genuinely accelerated implementation using the ORIGINAL
Green--Laplace separable operators and FFTs. They explicitly rejected the
assistant's added interpolation/low-rank/TT-SVD/cache route and requested its
deletion. That route has been deleted locally and on 93, including source
bundles and the approximately 606 MiB saved plan. Historical reports remain
only as evidence of the abandoned experiment. Do not restore or extend it.

The user then requested commit/push and a handoff. This checkpoint does NOT
claim completion of the requested fast GL time-series implementation.

Critical distinction: using GL quadrature inside pair/triple enumeration is
NOT the original fast GL executor. Earlier OW3D and HOS time-series comparisons
used precisely those reference executors. Their matching small fixtures do
not establish acceleration. Read SECOND_ORDER_IMPLEMENTATION.md first; it
contains the audited call chains and historical quadrature choices.

## Authoritative workspace and original GL entry points

- Workspace: `C:/Users/spet5947/Documents/green-laplace-unidirectional-time-series`.
- Branch: `codex/unidirectional-time-series`; preserve original public APIs.
- `src/gl_spectral_surface.m` deposits the native spatial analytic spectrum.
- `src/internal/green_laplace_eta22.m` selects the native second-order graph.
- `src/internal/finite_depth_directional_pure_gl8_eta22.m` loops over GL nodes,
  applies filtered IFFTs, pointwise forcing products and FFT output operators.
- `src/internal/gl_no_stokes_eta33.m` is the retained pure third-order graph.
- `research/unidirectional_time_series/test_time_pair_reference.m` compares
  temporal pair references to native spatial GL on small fixtures.

Current input requirement remains the observed/separated first harmonic at
the probe plus the initial COMPLEX directional spectrum, not total eta(t),
not an independently supplied higher-order reference. Current frequency-test
request is every nonzero temporal bin through 4omega_p; do not silently restore
the old initial-spectrum min/max cutoff or a kh>.5 deletion. Any applicability
or convergence limitation must be stated separately.

The unresolved interface is explicit: temporal-bin dispersion wavevectors
need not lie on the native spatial FFT lattice. Do not silently round them,
replace the observed first harmonic with initial linear propagation, change
the directional assumption, or introduce another physical closure. Audit the
native interface and derive/validate the required adaptation first. If it
cannot yet be used faithfully, say that; do not label a new approximation as
the original implementation. Pair/triple references are for small checks only.

## Completed GPU HOS evidence

92 is Tesla P40 24 GiB. Connect to `root@192.168.2.92:22` through
`root@60.188.112.99:60093`, using the existing local
`~/.ssh/id_ed25519_cursor` key. Public port 60092 is not the working route.
Never change clocks, power limits, ECC, thermal protection or drivers.

Completed low family: kpHs/2=.02 (NOT monochromatic akp), Hs=1.43369 m,
h=35.84229 m, kph=1, kp=.0279/m, Tp=13.762 s. JONSWAP gamma3.3,
11128 initial parents, four global phase shifts. HOS base grid1024x512,
M5/q3, FP64 and absolute time tolerance1e-12. Every phase completed1100.8 s,
5504 accepted steps, zero rejections,5505 finite probe records.

- Original GPU run on92: `/root/hos-low-80tp-20260927T052004Z`.
- MPS replay: `/root/hos-mps-20260927T090455Z`.
- Strict ordinary simultaneous-start baseline:
  `/root/hos-ordinary-full-20260927T123931Z`.
- Same full workload: ordinary12257.2766 s; MPS9527.7277 s;
  speedup1.28648x, waiting-time reduction22.2688%. Full probe/energy/stored
  checkpoint comparisons have zero measured difference. MPS service shut down.
- Read MPS_REPLAY.md. Earlier adoption run had a solo phase-zero head start;
  do not use its parallel-stage timer as a clean batch comparison.

The GPU q3 convention matches an isolated correction of official MPI global
zero-mode gravity. It is not byte-identical to uncorrected MPI8 evolution.
The original CPU production solver/jobs were not changed. Earlier CPU8/GPU
two-second FP64 timings gave about9.1--9.2x; these are bounded application-loop
measurements. GPL-derived source notices are in NOTICE.md.

## Completed second-order postprocessing and third-order limitations

93 postprocessing root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-low-postprocess-20260927T085331Z`.
Run MATLAB as **lxy**, not root. Root failed at sign-in; failed logs were kept.
R2026a executable: `/home/lxy/Desktop/matlabr2026a/bin/matlab`.
92 has R2025b installed too, but do not imply GPU execution merely from host.

Full second-harmonic GL/HOS comparison completed. It uses pairwise GL16 with
3.75/1.875-degree directional allocations. At1.875deg, raw fixed10--70Tp L2
errors for five probes are4.21%,3.30%,3.51%,3.19%,4.56%; common-band values
are similar. The original input retained149 temporal bins,99.8065% positive
energy, but changed the center first-harmonic waveform by4.40% in full-record
L2. A10s record truncation changed GL prediction by3.76%; not fully converged.
Read POSTPROCESS_GPU_LOW.md; plots/compact metrics are in ignored local
`artifacts/hos_gpu/postprocess-20260927T085331Z/`.

Higher-harmonic pilot root on93:
`/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-higher-harmonics-20260927T145332Z`.
- User explicitly requires **R for2-minus**. R4 time-pair adapter gave12.43%
  low-frequency-band L2 at center,1.875deg; R6 gave12.56%. Raw R4 L2 was11.88%.
  Preserve same-frequency/nonzero-spatial-K contributions before explicitly
  labelled temporal-DC/oscillatory-band projection. Strict spatial zero omitted.
- Third-harmonic direct GL30deg/rank8 pilot had92.14% raw error. It is not a
  successful validation. Direct15deg and full-band jobs were intentionally
  stopped on the user's request to use fast GL.
- Pure-linear controls demonstrate material finite-record Hilbert leakage.
  The user explicitly requested mean removal: this alone did not improve the
  third-harmonic error. A fixed end taper reduced leakage but left about70%
  discrepancy. Preserve raw, demeaned and tapered results separately.
- Initial absence of eta33 does NOT prevent HOS generating third-order content;
  do not use that fact to dismiss comparison. Free/transient components and
  extraction contamination need evidence, not an assumed explanation.
- The removed numerical-separation route reproduced its direct reference and
  completed a through-4omega_p sensitivity test, but its timings are NOT the
  original GL's performance. Do not present it as the requested fast solution.

Private R adapter snapshots remain only in ignored artifacts/remote storage.
Do not vendor private SPARK source into this public repository. No raw large
fields were fetched locally. Do not delete original GL, R, HOS or data while
cleaning an abandoned route.

## Separate concurrent project work

Another workstream has committed a new kpHs/2=.06 CPU HOS campaign. Current
AGENTS.md and `hos_ocean/random/JONSWAP_KPHS006.md` describe it. Its base is
`hos-jonswap-kphs006-80tp-20260927-v2` on93. Do not restart, stop, overwrite or
claim its live status from this handoff. The .12 failure and .02 completion
are distinct. Existing OW3D jobs are also outside this task's process scope.

## Next work

1. Start from the original native GL node/FFT graph above, not time-pair/triple
   evaluation or the deleted interpolation route.
2. Make the time-record-to-native-operator input contract explicit and verify
   preservation of supplied first-harmonic information and directional prior.
3. Derive new identities in Wolfram and check MATLAB implementation against
   small direct/MF12 SECOND-order fixtures. Do not add MF12 third-order work
   contrary to the existing scope.
4. Require actual same-input/same-output/same-host timings against MF12 before
   claiming acceleration; report preparation and repeated evaluation separately.
5. Preserve the fixed scientific inputs, no-Stokes scope, native output band,
   and user-requested R2-minus path. Do not tune against HOS reference curves.

No new fast GL time-series solver has been delivered at this checkpoint.
