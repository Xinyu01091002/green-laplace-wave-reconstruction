# Compact directional OW3D redesign, 2026-09-25

This replaces the previous depth/steepness/time-step campaign proposal.
The user selected only **kpd=1, Akp=.12, a focused wavegroup and a random-wave
case**, with no dt study or solver parallelization project. New fields and
processing remain remote. Published GL code and all historical runs are unchanged.

## Current outcome

- Four focused-wavegroup initial eta/psi pairs and native input files exist
  remotely. The independent second-order initialization checks passed.
- A reduced-domain first-order geometry check and a unit-amplitude finite
  random-field boundary prototype ran remotely in MATLAB R2026a.
- Remote GL release checks passed 8/8, from the independent clone without
  SWORD paths. All four new MATLAB drivers pass Code Analyzer.
- The actual OW3D binary parsed the corrected input, built its multigrid
  hierarchy, wrote initial EP/kinematics and entered time stepping. It did
  **not complete the first advance within the 300 s smoke limit**. This is
  not a successful two-step propagation test or proof of instability.
- Native t=0 kinematics surface phi does not match the prescribed surface
  potential. Initial eta/psi must come from EP; kinematics begins at t=.2 s.
- Random-wave physical normalization remains pending the user question.
  Absorbing-zone evolution and usable record duration are not yet validated.
  No production OW3D campaign or random nonlinear initial state was launched.

## Wavegroup input

| Parameter | Adopted value |
| --- | --- |
| Gravity, peak wavenumber | g=9.81 m/s2, kp=.0279 1/m |
| Depth | h=1/kp=35.8422939068 m |
| Peak wavelength, period | 225.203774451 m; 13.7619975 s from finite-depth dispersion |
| Domain | 50 lambda by 20 lambda = 11260.1887225 by 4504.0754890 m |
| Native OW3D nodes | 1025 by 257 by 17, excluding ghost nodes |
| Unique horizontal FFT grid | 1024 by 256 |
| Horizontal intervals | dx=10.996278 m, dy=17.594045 m |
| Vertical grid | Native sine-clustered grid, with one bottom ghost |
| Linear focus amplitude | A=.12/kp=4.301075269 m |
| Focus position/time | (Lx/2,Ly/2), t=110 s |
| Proposed run / integration / output | 0--220 s; dt=.2 s; saved interval .2 s |
| Phase runs | 0, 90, 180, 270 degrees |
| Kinematics range | all x, physical y indices 125--133, time indices 2--1101 |
| Initial sample | eta and true surface psi from EP_00000.bin / initial file |
| Sparse EP cadence | 55 s, including initial and final fields |

The box uses exact multiples of 2*pi/kp. This differs by about 0.09% from
halving the historical dimensions, which used lambda=225 m. Direction is
explicitly +x, with `Re(C exp(i k dot x - i omega t))` and linear surface
potential coefficient `-i*g*C/omega`. Do not import the old positive-time
phase sign or its nominal Tp=12 s without conversion.

The prescribed focused modal-amplitude density has a Gaussian radial lower
width .004606 1/m, upper width kp/sqrt(2 log(10^8)), and Gaussian **amplitude**
angular width 25 degrees. That is not a claim of 25-degree energy RMS spread.
Cartesian quadrature includes the polar-coordinate 1/k Jacobian and strict
kx>0 support. This is a newly declared spectrum; the edited historical VWA
generator does not establish exact identity with the old input fields.

A prescribed relative modal-amplitude threshold 1e-4 leaves 3777 parents.
The removed L1 amplitude fraction is 5.7202e-5 and removed squared-amplitude
fraction 5.7741e-9. Remaining amplitudes are normalized to the declared
linear focus A before any OW3D or GL comparison. This is input construction,
not fitting to a numerical reference. The exact complex parent vector is
saved for later GL use. Second-order sums stay within the native grid.

An untruncated first-order geometry audit covered focus leads 60/90/110/130 s.
At 110 s, over 0--220 s, maximum wall analytic-eta amplitude / focus A was
3.1405e-5; the analogous potential ratio was 3.3201e-5. The full second-order
initial surface potential is broader: its boundary value reaches about
1.31% of gA/wp. Its initial surface-potential normal-coordinate gradients
are at most 1.374e-4 (x edges) and 1.13e-5 (y edges), normalized by g*kp*A/wp.
These are surface diagnostics, not exact volume boundary residuals or a
reduced-domain propagation certificate. No boundary correction was fitted.

## Independent initialization and output checks

Initialization uses the existing external MF12 **order 2** coefficients and
surface evaluator: eta11+eta20+eta22 and their true surface potentials.
It includes nonzero difference modes, including stationary nonzero spatial
K; strict spatial means are zero, with no imposed return current. Third-order
MF12 was not computed. No GL high-order result was supplied to the reference
initializer, and the old VWA/Stokes-style third-to-fifth-order additions were
not copied. Third-order startup/free-wave content is not certified by this
second-order initialization.

The preallocated MF12 order-2 path was checked against its independent direct
coefficient/evaluation path on three declared parents, covering both sum and
difference terms and surface potential. Four-phase first-sector extraction
of the prepared fields recovers the declared linear eta/psi to roundoff.
Nonfinite and output-support checks passed. The external MF12 files retain
their authorship headers and are hashed under an ignored remote dependency
directory; they are not vendored into this Git repository.

The current binary's list-directed input reader continues to the following
line when optional fields are omitted without the old `<-` comments. The
new writer explicitly supplies the acceleration limit, density and disabled
breaking-model fields. Earlier malformed input attempts are preserved.
`wavegroup-v4` is the adopted input directory; v1--v3 are not the current
run configuration. Physical field arrays are unchanged from `wavegroup-v2`.

The t=0 native-output audit found exact eta equality between EP and
kinematics, but maximum phi discrepancy 32.6479 m2/s: the native volume
potential has not been solved at that point in setup. Production selection
therefore skips the t=0 kinematics record. Later top-level phi/EP parity
still requires a completed advance; the timed smoke attempt did not supply
one. The default two-step audit retains that requirement.

## Random-wave design still pending a physical definition

The prototype uses a fixed seed 20260925 and the same modal-amplitude shape,
with independent random phases. It is a **finite random field with an interior
sea region**, not a stationary infinite or periodically evolved sea. Native
OW3D periodic walls are not assumed.

A prescribed C2 quintic edge taper of width 5 lambda along x and 3 lambda
along y is applied to first-order eta only. The resulting forward spectrum
is saved, and linear psi is recomputed from that spectrum. No independent
window is applied to nonlinear eta/phi. In the central region
x/lambda=10--40, y/lambda=5--15, the relative eta change from the untapered
prototype is 1.113e-4. Wall eta / core RMS is 3.690e-4; wall psi / (g*sigma/wp)
is about .0062. These checks concern initial fields only.

The pending user choice is whether Akp=.12 means **kp*Hs/2=.12**
(Hs=8.602150538 m, sigma_eta=2.150537634 m), or unchanged wavegroup modal
amplitudes with randomized phases and a separately reported actual Hs.
No physical normalization was chosen silently. Nonlinear random initial
conditions, absorbing-zone input and an accepted duration remain unfinished.

## Resources and provenance

The reduced extended volume grid has 4,787,874 nodes, 25.18% of the earlier
2049x513x17 proposal. Compared with the historical 2049x513x9 run it is
about 45.3%, because the new proposal doubles vertical resolution.
Four-phase 0.2--220 s native kinematics plus six EP files per phase are
estimated at **71.87 GiB**, before inputs and derived products.

The bounded OW3D smoke recorded one thread and peak RSS 24,518,468 KiB
(23.38 GiB), with 302.35 s elapsed before timeout. This includes startup
and an unfinished first advance, not steady per-step timing. Eight jobs
at that observed peak plus 25% memory allowance would use about 234 GiB,
but batch throughput and later peak memory have not been verified.

Actual remote project: `/home/lxy/green-laplace-unidirectional-time-series`.
Actual run root: `results/ow3d_redesign/20260925-compact/` beneath it.
The GitHub HTTPS clone failed with a TLS error; an independent clone was
created from a bundle of the local research branch at
`441f8bd161ab1562e945a5aa28c32c70b3695f54`, with origin set to the public GitHub
URL. No historical .git/worktree or run directory was copied or modified.

Executions used hashed source snapshots under `source-v1` through `source-v6`.
The final manifest records the committed source, snapshot hashes, solver
hash, external MF12 hashes and generated input hashes. Relevant artifacts:

- `geometry-v1/`: declared first-order spectrum and geometry screening.
- `wavegroup-v4/`: current four input pairs and effective configuration.
- `wavegroup-v2/initial_fields.mat`: frozen independent physical arrays.
- `random-boundary-v1/`: unscaled random prototype.
- `smoke-v4/`: correctly parsed but timed-out OW3D attempt and initial-output audit.
- `release-tests.log`: 8/8 remote GL checks.

All new numerical execution and field processing occurred remotely. Only
small logs/status metrics were read through SSH; no new raw fields were
downloaded. No active production OW3D process is left by this design task.

Next steps are to settle the random normalization, complete its boundary
configuration, and determine whether the actual serial OW3D first advance
finishes within an acceptable time. The two-step smoke has not passed, so
the present input package is not labelled propagation-validated or ready
for an unattended large campaign.
