# Unidirectional deep-water four-phase HOS / MF12 experiment

Completed, including absolute-tolerance checks. See [RESULTS.md](RESULTS.md)
for final errors, measured runtime/memory and interpretation limits.

User-authorized on 2026-09-26: first run low steepness, then high; initialize
with spectral MF12 excluding 31; compare spatial second/third harmonics
reconstructed from the current first harmonic at 3Tp and 20Tp; record time
and peak memory. This supersedes earlier MF12-order-two-only limits for this
experiment only. Frozen OW3D jobs are not modified.

## Frozen setup

- Local settings reference: `ChatGPT/ESC-Time Series/configs/ow3d_esc_input_plan.json`
  and `src/matlab/+esc_ts/make_vwa_general_source_modes.m`.
- g=9.81 m/s2; kp=0.0279 1/m; asymmetric semi-Gaussian amplitude shape,
  Alpha=1, left width 0.004606 1/m, right width kp/sqrt(2 log(10)).
- Akp=0.02 and 0.12 denote potential linear focusing steepness, not current
  local steepness. Nominal focus: x=L/2 at t=40Tp after initialization.
- Domain 68 peak wavelengths, 4096 unique periodic points. No duplicate endpoint.
- Finite-depth approximation to deep water: kph=20 (h=716.8458781362 m),
  minimum retained parent kh=5.5882352941; Tp=12.01000706916 s.
  This is not a literal infinite-depth solution, especially for long difference modes.
- Retain 190 source modes using the local source rule: 99.999% of shape L1,
  not an energy fraction. Actual retained fraction 0.999990480432779. Same
  fixed parent support is used at all evaluation times; omitted evolved
  first-harmonic norm is reported, never hidden or tuned to prediction errors.
- Each phase has independently evaluated coefficients/fields, 0/90/180/270 deg.
  MF12 supplies 11+20+22+33 in eta and true surface psi. Mixed-sign order-three
  branches and primary-harmonic muStar are excluded; frequencies are linear.
- HOS M=5, qx=5 (full product dealiasing at this M), Ny=1, Ta=0,
  no breaking/filter/dissipation. Duration 20Tp, snapshots every Tp/4.

## Source and execution

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-deepwater-fourphase-20260926-v1`.

Local MF12 source:
`C:/Research/OW3D_benchmark/unidirectional_wave_generation_package/deps/mf12`.
Frozen archive SHA-256:
`65011bff7ae59dfba822430f49d157066c60accbd0968e599aca48ce781aa8bb`.
The wrapper explicitly sets muStar=0 in addition to
`disable_third_order_correction=true`, and sets linear modal frequencies.
No external MF12 source file is edited.

HOS binary SHA-256:
`21a171f03c32a49ae05b62e4896d5b081c26a895372552d6ba563ba8fd44f9b4`.
The tested IO-only source patch is documented in the parent PRECISION.md.
MATLAB is used for numerical fields and comparisons; Python only orchestrates.

## Measurement and interpretation

The four phase records are combined and separated using positive spatial
wavenumber projection for the first and third harmonics. These are harmonic
phase sectors, not exact perturbation orders. Second harmonic uses the
alternating four-phase sum. No phase, gain, offset, sign or spatial shift is fitted.

At t=0,3Tp,20Tp the current first-harmonic spectrum is projected onto the fixed
input support and supplied to spectral MF12 at that snapshot. Predicted eta22
and eta33 are compared with the simultaneous HOS sectors. Full-domain relative
L2/Linf and a +/-2 lambda_p main-group window are reported; the window center
is selected from the first-harmonic envelope alone.

All output records are checked for finite values and all three requested
snapshots are required. Native t=0 eta/psi must agree with the prescribed
initial fields to relative L2 <1e-12. The synthetic initial four-phase sector
recovery is checked separately before any HOS process is launched.

GNU time records process wall/CPU time and peak RSS for each independent
HOS case and each MATLAB preparation/analysis stage. MATLAB stage measurements
include process startup, while internal coefficient/reconstruction timers
are reported separately. Four HOS phases run concurrently, with one compute
thread per process. Initial low/high batches run sequentially. Per-process
RSS is not an observed aggregate peak, and summing case wall times is not
batch elapsed time.

## Numerical precision audit

The initial 1e-10 tolerance is HOS's default absolute-error mode, applied to
nondimensional modal variables. The low-amplitude higher harmonics are tiny;
a misleading low-steepness result triggered a numerical audit at 1e-12 and
1e-14 with byte-identical initial fields. Original results are retained under
akp002/akp012, and tighter results under tol12/ and tol14/. No model formula
is adjusted to reduce observed error. See the final results record for which
precision is supported by convergence evidence.