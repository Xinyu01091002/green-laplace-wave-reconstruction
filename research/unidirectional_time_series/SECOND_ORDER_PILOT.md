# First eta1(t) to eta22(t) pilot

## Scope and fixed choices

The user selected separated first-order time series and requested GL, spectral
MF12, VWA and Walker against the existing unidirectional OW3D data. This pilot
reconstructs elevation in the second-order positive sum sector only. It does
not include difference frequencies, mean flow, inversion or nonlinear evolution.

Read-only data source:
`C:/Research/VWA/VWA time series/unidirectional/timeseriesdata`.
The chosen family is `T_init-40_Tp_Alpha_1.0_Akp_002_kd1.0_phi_{0,90,180,270}`.
Probe index 3800 and steps 2800:4:4200 come from the existing
`paperplot_VWA_time_series.m`, not a search for better agreement. The four
records must have identical spatial grids and physical/time metadata. The
binary reader follows the existing `ReadBinFile.m` Fortran-record layout,
adds record-length/finite-value checks and reads no historical runtime code.
The file headers include ghost points; retain the historical index convention.

Read dt, depth and g from each OceanWave3D.inp, and kp/Tp from OW_readme.txt.
The old plotting script's dt=0.15 and Tp=12 are not authoritative for this case.
Plots use simulation elapsed seconds, with no fitted or inferred focus shift.

Four-phase separation uses the temporal Hilbert sign of the existing VWA
script. The first phase sector is projected onto the common parent support
and supplied as the declared eta1 input. The second phase sector is held out
as the reference. Neither sector is a pure perturbation order: the input can
contain higher-order primary-harmonic content and the reference can contain
higher-order second-harmonic content/phase aliasing. Akp=0.02 limits, but does
not certify the absence of, that contamination.

Shared support is fixed before comparison: all positive temporal FFT bins
with kh>=0.3 and 2*omega below the native temporal Nyquist. No energy-ranking,
maximum-component cap, padding, gain, offset, sign, phase or time adjustment
is applied. DC is excluded and reported; this is an explicit input projection,
not a claim that the original phase-sector record has zero mean. Retained
energy and input projection error are reported. No OW3D reference values
choose this support. The reference record is not retrospectively filtered
to improve agreement. Finite-window leakage and input projection remain
limitations. Full-window and fixed middle-half metrics are both retained.

## GL and comparison conventions

Use physical amplitudes A in `eta1(t)=Re sum A exp(-i omega (t-t0))`;
for a real temporal record, `A=2*conj(fft(eta1)/N)` on positive FFT bins.
At the probe, the spatial phase is already included in A. Solve the linear
dispersion relation for each positive k; use the true pair sums k_i+k_j and
omega_i+omega_j. Do not infer a bound-wave k from its summed frequency.

`gl_eta22_time_pairs` transfers the existing released, symmetric ordered-pair
GL kernel from `tests/finite_depth_directional_pure_gl8_ordered_pair.m` to
single-point temporal synthesis. Dimensionless kernel times A_i*A_j/h gives
physical elevation. There is no new physical closure and no Stokes correction.
The stable exponential form is certified in Wolfram against the original
sinh/cosh form. The primary rank is 6; 8 and 12 are prescribed diagnostics,
not ranks selected against OW3D or MF12. This O(number_of_parents^2) reference
executor is not a fast temporal FFT implementation or a speedup claim.

Spectral MF12 uses the external `mf12_spectral_coefficients` at
`C:/Research/spectral domain implementation of wave interaction theory`.
The adapter follows the repository's `gl_mf12_order2_comparison`, retaining
only self and positive-sum coefficients. It does not use the old time adapter's
99%-energy / 64-component defaults. All methods receive exactly the same A/k.

VWA uses `Re[z * sum B22(kh)*k*A*exp(-i*omega*t)]`, where
`B22=(3-tanh(kh)^2)/(4*tanh(kh)^3)` and z is the analytic first-order input.
This is the second-order temporal expression in the existing VWA script.
Walker uses `kp*B22(kp*h)*Re[z^2]`, equivalent to the Appendix A second-order
expression implemented in `compare_deepwater_virtual_superharmonics.m`.
Walker, Taylor and Eatock Taylor (2004), *The shape of large surface waves on
the open sea and the Draupner New Year wave*, DOI 10.1016/j.apor.2005.02.001,
is available locally as extracted text under `VWA time series/UWA_data_20260401`.
Its coefficients are an independent comparison, never a repair to GL.

## Commands

From the new project root:

```powershell
wolframscript -file research/unidirectional_time_series/verify_time_pair_algebra.wl
matlab -batch "restoredefaultpath; addpath(pwd); addpath('research/unidirectional_time_series'); test_time_pair_reference('C:/Research/spectral domain implementation of wave interaction theory');"
matlab -batch "restoredefaultpath; addpath(pwd); addpath('research/unidirectional_time_series'); run_ow3d_eta22_pilot('C:/Research/VWA/VWA time series/unidirectional/timeseriesdata','C:/Research/spectral domain implementation of wave interaction theory');"
```

Numerical results and PNG/PDF figures are in
`results/unidirectional_time_series/ow3d_kh1_alpha1_akp002/`.
`source_files.csv` hashes every consumed EP file; `metadata_files.csv` hashes
the physical input records. `pilot.mat`, `report.json`, `time_series.csv` and
`metrics.csv` retain raw numerical evidence. Original data remain read-only.
Symbolic and implementation checks are under `artifacts/unidirectional_time_series/`.

## Implementation validation

Wolfram verifies four algebraic identities: stable exponential evaluation,
dynamic-source symmetrization, kinematic-source symmetrization and pair symmetry.
MATLAB checks ranks 4/6/8/12 against the released spatial GL API at three times
and three spatial points, complex amplitudes included. The initial 1e-10
point-relative gate exposed a maximum 5.389e-10 discrepancy; the maximum
absolute discrepancy is 6.370e-14 m. The independent frozen ordered GL8
reference agrees to 3.581e-16 relative. The test therefore retains both the
spatial relative bound (1e-8) and absolute bound (1e-12 m), plus the stricter
frozen-kernel gate (1e-12 relative); no kernel or quadrature parameter was
changed to pass. MF12 temporal/spatial relative difference is 1.318e-16 and
temporal FFT amplitude recovery is 1.144e-15. These are implementation and
algebra checks, not broad physical-accuracy certification.

## Executed first-case results

MATLAB R2022b read all 1404 EP files (351 samples times four phases), checked
finite fields and matching spatial grids, and wrote per-file SHA-256 records.
Actual h=35.84229 m, dt_output=0.6005 s, kp=0.0279 rad/m, Tp=13.761997 s;
probe x=28514.208251953125 m and y=0. The elapsed-time interval is
420.35 to 630.525 s. The shared support contains 82 parents, kh=0.348685 to
24.574502, retaining 99.9999635994% of positive-frequency input energy.
The input projection relative L2 is 0.0603516%; excluded DC is 2.20207e-6 m.
No data-dependent amplitude rescaling or reference alignment was performed.

Full-window raw relative L2 against the OW3D second phase sector:

| Method | Relative L2 (%) |
|---|---:|
| GL6 (declared primary) | 0.236271 |
| GL8 | 0.107809 |
| GL12 | 0.094233 |
| Spectral MF12 | 0.094204 |
| VWA | 0.895420 |
| Walker | 11.500584 |

GL6/8/12 relative L2 against spectral MF12 on the identical first-order input
is 0.161379%, 0.024061% and 0.001181%, respectively. This is a prescribed-rank
comparison, not numerical fitting or a proof of exactness. MF12's recovered
linear record agrees with the common input to 2.142e-14 relative L2.
The complete metrics also include Linf, Q, norm ratios, signed extrema and
the fixed middle-half window. The figure displays the full window so the
whole packet remains visible. The middle-half metric is not a focus window.

All new MATLAB files passed Code Analyzer with no messages. Source execution
and input/output hashes are retained in the ignored result directory. No
remote calculation, new OW3D run or Git push was performed.

Interpretation: the first small-amplitude case supports the feasibility of
the GL pairwise time-series route and its close agreement with spectral MF12.
The approximately 0.094% OW3D discrepancy shared by GL12/MF12 is not assigned
to a single cause: phase-sector contamination, finite record effects, input
projection and OW3D numerical error have not been separately measured here.
The next validation should use the existing Alpha=8 case at the same kh and
Akp, with the same predeclared rules, followed by a depth comparison; do not
adjust the method based on the current OW3D waveform.

## Shallower-water follow-up: kph=0.5

At the user's request, the runner now accepts an optional third argument
`kph` (default 1). This selects a data family and a separate output directory;
it does not change the GL kernel or quadrature. Alpha=1, Akp=0.02, probe index
3800, steps 2800:4:4200, ranks 6/8/12, support guards and metric definitions
remain unchanged. Each depth uses its own source metadata, so the same steps
do not imply the same physical time interval or number of peak periods.

```matlab
run_ow3d_eta22_pilot( ...
    'C:/Research/VWA/VWA time series/unidirectional/timeseriesdata', ...
    'C:/Research/spectral domain implementation of wave interaction theory',0.5);
```

The kh=1 results remain intact. The new results are written under
`results/unidirectional_time_series/ow3d_kh0p5_alpha1_akp002/`.

All 1404 files were read successfully. Source metadata give h=17.92115 m,
dt_output=0.800668 s and Tp=17.667179 s; the elapsed-time interval is 560.4676
to 840.7014 s. The identical support rule retains 78 parents and 99.9947792%
of positive-frequency input energy. The input projection relative L2 is
0.723736%, and excluded DC is 5.24465e-5 m. Minimum parent kh is 0.306859.

| Method | Full-window relative L2 against OW3D (%) |
|---|---:|
| GL6 | 23.5069 |
| GL8 | 23.4659 |
| GL12 | 23.4560 |
| Spectral MF12 | 23.4545 |
| VWA | 23.7313 |
| Walker | 36.5712 |

GL6/8/12 against spectral MF12: 1.31954%, 0.407012%, 0.0525246% relative L2.
The MF12 first-order reconstruction agrees with the common input to
4.775e-14 relative L2. Both modified MATLAB files pass Code Analyzer.

The complete waveform shows an OW3D second-phase-sector oscillation after
approximately 730 s, after the first-order packet has passed. All four
bound-wave methods lack that tail. Its physical/numerical origin has not
been identified, and the raw phase sector cannot yet be treated as a pure
bound eta22 reference. Do not attribute the full-window discrepancy solely
to the GL approximation or suppress the tail to improve agreement. This
time window covers a different part of the passage than kh=1; its fixed
middle half is particularly unsuitable as a main-packet metric.

An exploratory `window_diagnostic.csv` partitions the *saved* results at
the time-window midpoint. It reports first-half errors and the second-half
share of squared residual, without changing any inputs, predictions or
primary full-window metrics. This diagnostic was added after inspecting
the complete plot and is not a predeclared validation gate.
