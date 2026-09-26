# Directional four-phase HOS production launch

Active run: `hos-directional-fourphase-20260926-v2`.
Started: 2026-09-26 04:12:14 UTC (05:12:14 Europe/London).
All four phases were observed advancing; all four initial-probe checks passed.
This is launch evidence, not completed evolution or GL accuracy validation.

## Numerical specification

| Item | Setting |
| --- | --- |
| Gravity | 9.81 m/s2 |
| Peak wavenumber | .0279 1/m |
| Peak wavelength | 225.2037744509 m |
| Peak period | 13.7619970876 s |
| Constant depth | 35.8422939068 m; kph=1 |
| Potential linear focusing steepness | Akp=.12; nominal linear focus A=4.301075269 m |
| Domain | 11260.188722544 x 4504.075489018 m = 50 x 20 lambda_p |
| Boundary | Periodic in x and y |
| Unique physical horizontal grid | 1024 x 256 |
| dx, dy | 10.9962780494 m, 17.5940448790 m |
| HOS nonlinearity order | M=5 |
| Dealiasing | qx=qy=3, partial at M=5; extended FFT grid 2048 x 512 |
| Time integration | Six-stage embedded Cash-Karp Runge-Kutta 5(4), adaptive; analytical linear evolution |
| Error control | Absolute tolerance 1e-12 on HOS nondimensional modal variables |
| Nonlinear activation | Ta=0; no ramp |
| Breaking or added dissipation | None |
| Physical duration | 0--220 s = 15.986051923 Tp |
| Nominal linear focus | t=110 s = 7.993025961 Tp, x=Lx/2, y=Ly/2 |
| Output interval | .2 s = approximately Tp/68.81, not a fixed integration step |
| Samples | 1101 expected per probe including t=0 |
| Global phases | 0,90,180,270 degrees |
| MPI | 4 ranks per phase, four phases concurrently, 16 compute ranks total |

The initial condition is the exact frozen directional OW3D wavegroup's
independent MF12 order-2 eta/true-surface-psi pair: 11+20+22, no33 or31.
It is not the earlier unidirectional third-order initialization. No GL
higher-order field is used. The input has 3777 retained parents, +x mean
direction, 25-degree Gaussian AMPLITUDE angular sigma, and the declared
semi-Gaussian radial widths .004606 1/m below kp and kp/sqrt(2*log(10^8))
above kp. The original support and normalization are preserved.

Input source, read only and copied into this run:
`/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z/inputs/wavegroup/initial_fields.mat`.
MPI inputs are contiguous 64-row y slabs, with high-precision decimal IO.
Full initialization fields and first-order spectrum are retained in inputs/.

## Probe output

All five probes share x=5630.094361272 m. Their y coordinates are:
2252.037744509 (center), 2199.255609872, 2304.819879146,
2181.661564993 and 2322.413924025 m. Actual transverse offsets are
0, +/-52.782134637, +/-70.376179516 m, or 0, +/-0.234375, +/-0.3125 lambda_p.
They are native grid points, with 1-based (x,y) indices (513,129), (513,126),
(513,132), (513,125), (513,133).

Saved probe columns are time and five eta values; they do not contain psi.
No dense spatial/volumic time series is written. First-harmonic directional
information remains available in the frozen full initialization. Native HOS
energy/volume diagnostics and runtime logs are also retained.

## MPI probe bug caught before accepted production

The first attempt, v1, was stopped automatically on an initial-probe mismatch.
Its inputs, logs and failure status remain intact. The pinned upstream MPI
probe initialization broadcast x/y coordinates but omitted nprobes; non-root
ranks consequently did not have the correct probe-loop count. This does not
invalidate the earlier MPI full-spatial-field equivalence test.

The v2 source broadcasts nprobes in both probe initialization/restart paths,
retains the existing high-precision IO formats, and flushes probe output for
prompt monitoring. Evolution equations are unchanged. A four-rank phase-90
preflight checked all five probes at t=0 and .2 s against the assembled native
spatial output: maximum absolute difference 1.59e-15 m. Full production then
passed all twenty initial probe checks, with maximum difference 5.41e-15 m.

Active binary SHA-256:
`93d64252f732449fad9d2171979c91e13b0b9de1d17a73d7587fa00721efcbb6`.
Patch SHA-256:
`747b0473aa4af380e3fa6d06db4d8294aafcf821363277ef9bfa1091243c504a`.
Do not use the older MPI binary for multi-rank probe production. Its full
spatial output remains the validated historical benchmark.

## Resource control and status

Rank CPU affinity is enforced with taskset: phase 0 uses CPUs 16--19,
phase 90 uses 20--23, phase 180 uses 24--27 and phase 270 uses 28--31.
Actual solver process affinities were checked in /proc. OpenMP, MKL and
OpenBLAS worker counts are one. Existing 16 OW3D processes were not modified.
At the launch check the four HOS jobs used approximately 3.9 GiB combined
sampled RSS; available host memory remained about 469 GiB. An 8 GiB aggregate
HOS-job RSS guard terminates only these jobs if exceeded.

Absolute tolerance is tighter than the 1e-8 MPI timing pilot, so the earlier
approximately two-hour extrapolation is not a measured completion estimate
for this batch. Use actual progress as evolution approaches focus.

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-fourphase-20260926-v2`.
Screen session: `hos-directional-4x4-v2-20260926`.
Read status.json for per-case progress, PIDs, affinity, RSS and initial checks;
settings.json and manifest.json record physical configuration and file hashes.
On completion each case must exit zero and provide all 1101 finite records
at expected times. Only then is its completion.json written. The root reaches
completed and writes finished.utc after all four validated case outputs.
This controller does not claim completion of subsequent GL reconstruction.
Per-rank GNU time files record wall/CPU time and peak RSS on exit.

Local metadata copies are under artifacts/hos_ocean/directional-production/.
No raw simulation data are fetched locally and no commit/push was performed.