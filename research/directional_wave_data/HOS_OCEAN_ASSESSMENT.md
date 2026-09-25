# HOS-Ocean assessment and phase-only random input

2026-09-25. This is a solver-selection assessment, not an HOS propagation
result or a measured speedup. The user chose to preserve the focused
wavegroup's first-order modal amplitudes and randomize phases only.

## Random input actually prepared remotely

`prepare_compact_wavegroup` now accepts an optional fourth argument pointing
to the frozen focused-wavegroup `initial_fields.mat`. It applies
`C_random=C_group*exp(i*theta)` with independent uniform theta and fixed
seed 20260925. The four phase runs share that same random realization;
their additional phase shifts are 0/90/180/270 degrees.

No amplitude renormalization or spatial taper is applied. The earlier
unit-RMS edge-taper prototype is retained as history, not selected input.
Observed maximum first-order modal-amplitude change is 1.800e-16 relative.
The common linear spatial RMS is 0.1019037931 m, giving the variance-based
Hs=4*sigma=0.4076151725 m and kp*Hs/2=0.005686231656. The label Akp=.12
is the **potential focusing steepness** kp*sum(abs(C)), not the random
sea-state steepness. This follows the user choice; no target Hs is imposed.

Four independent MF12 order-2 random eta/true-surface-psi fields were
generated, with the same finite/support, direct-MF12 parity, four-phase
linear recovery and Parseval checks as the focused family. Third-order
MF12 was not computed. Output is on the unique 1024x256 periodic grid.
The fields are generic solver inputs, not already formatted/validated HOS
input files. Closed-wall OW3D configurations were deliberately not generated
for the untapered random realization because those walls would change the
physical experiment.

Remote directory under `/home/lxy/green-laplace-unidirectional-time-series/`:
`results/ow3d_redesign/20260925-compact/random-phase-only-v1/`.
It contains `initial_fields.mat`, `initialization.json`, and per-phase
`initial_surface.mat` files. `source-phase-randomization-v1/` records the
executed source hash and base commit. New raw fields remain remote.

## What the official HOS-Ocean documentation supports

- Constant finite depth and nonlinear open-ocean directional wave fields:
  [overview](https://lheea.gitlab.io/HOS-Ocean/index.html).
- External eta and free-surface potential initialization, including a
  dimensional `Initial surface quantities` input route; physical and modal
  free-surface eta/psi output:
  [running and input/output documentation](https://lheea.gitlab.io/HOS-Ocean/running-HOS-Ocean.html).
- Surface spectral discretization, selectable HOS nonlinearity order,
  dealiasing and adaptive time-integration tolerance:
  [numerical parameter guidance](https://lheea.gitlab.io/HOS-Ocean/choice-numerical-parameters.html).
- Official Linux/Windows prebuilt binaries, or CMake/Fortran with FFTW3 and
  LAPACK; MPI/HDF5 are optional:
  [installation](https://lheea.gitlab.io/HOS-Ocean/getting-started.html).
- Current official source is linked on
  [GitLab](https://gitlab.com/lheea/HOS-Ocean). Pin a release/commit and binary
  hash; do not mix its current YAML schema with older GitHub input examples.

HOS-Ocean is a better structural match for these constant-depth, periodic
surface-wave experiments than a full-volume finite-difference solve. It
does not require an OW3D-style Nz grid. Lower memory and shorter runtimes
are an engineering expectation from the formulation, not an observed
performance result on this host. No HOS executable or Fortran compiler was
found on PATH or HOS installation in the bounded common-directory search;
that is not proof that no installation exists elsewhere. HOS has not been
installed, compiled, downloaded or run during this assessment.

## Proposed first HOS comparison

Use the same declared g, h, domain, physical units, origins, parent spectra
and independent MF12 order-2 eta/psi fields. Target the same unique physical
grid 1024x256; verify the selected version's real/Fourier storage conventions
before assigning its input grid numbers. Do not import OW3D ghost nodes or
its repeated endpoint. No Nz parameter is carried over.

Start from HOS order M=5, the documented standard choice. M controls the
evolution approximation; it is not an eta22/eta33 harmonic label. The
documented initial dealiasing choice is qx=qy=3 for directional waves;
this is partial rather than full dealiasing at M=5. Accuracy of the small
high-order components must be checked before treating them as reference
truth. This is not a request to reinstate the cancelled OW3D dt sweep.

Keep saved physical times 0,.2,.4,... s. HOS uses adaptive integration
tolerance; the saved interval is not its internal step size. Verify output
frequency units for the selected initialization type. Set and record the
initial nonlinear activation/ramp explicitly, so it does not silently
replace the prescribed already-bound initial state. No new physical
breaking/dissipation model, fitted gain, phase shift or GL correction is
introduced. Retain the same strict-zero spatial mode and potential-gauge
conventions when comparing eta20/psi components.

Prefer free-surface physical/modal output, potentially HDF5, and extract
probe eta/psi with MATLAB remotely. Do not assume the default probe file
contains surface potential without checking its actual columns. No volumic
kinematics is needed. A native t=0 input/output equality check is essential.

For the focused packet, check that periodic wraparound is absent from the
main-group window. For the random case, periodic evolution is the intended
setting; closed-wall OW3D is not automatically the same boundary-value
problem. OW3D remains an independent check where those boundary differences
can be made negligible. GL outputs still come solely from the declared
first-order input; HOS is an independent evolution reference.

## Concurrent runs

At the 20:45 UTC refresh the host had 48 visible logical CPUs and 738.6 GiB
available RAM, with no OW3D process. The earlier actual OW3D startup peak
was 23.382633 GiB with one thread, but no completed advance was measured.
With a 25% memory allowance:

| Independent OW3D jobs | Estimated aggregate RAM |
| ---: | ---: |
| 8 | 233.83 GiB |
| 16 | 467.65 GiB |
| 24 | 701.48 GiB |

Recommend eight jobs for the two four-phase families if OW3D is used; sixteen
is a conditional later target, not a tested scaling result. Twenty-four
exceeds the proposed 600 GiB calculation budget. CPU/memory bandwidth and
per-step runtime still need measurement. These numbers must not be reused
as an HOS memory estimate.

For HOS, first measure one fixed-version case, then aim for the same eight
independent jobs, initially one compute worker per job. More concurrency
depends on its measured memory/thread use and throughput. Existing MPI is
available if separately chosen, but no solver parallel rewrite is required.
