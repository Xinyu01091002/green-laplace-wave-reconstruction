# Short all-device RK4 validation

The fixed-step RK4 integrator in rk4.cuh evolves the degree-five HOS operator
on an explicitly retained rectangular input band. Every stage RHS and completed
state are projected to that band; the independent CPU reference uses exactly
the same Galerkin projection. This is a small, defined semidiscrete system,
not a claim that production HOS-Ocean uses these narrow bands or RK4.

The time stepper retains state, stage state and accumulated slopes on the
GPU. Step calls have no host transfers or allocations. Six additional real
buffers are allocated once. The existing HOS RHS handles its own work arrays.
No hardware settings or existing remote simulations are changed.

## Independent reference

make_time_reference.m uses MATLAB direct Fourier convolutions for every
degree-five RHS evaluation, then classical RK4 in coefficient space. No FFT
or GPU code is reused. State modes are |mx|<=3, |my|<=2, with the same initial
two directional mode pairs and mean as the instantaneous tests. Parameters:
Lx=9 m, Ly=7 m, h=1.3 m, g=9.81 m/s^2, M=5, physical grid 64x32.

Duration is only 0.2 s: 40 steps of 0.005 s. The reference is independently
repeated with 80 steps of 0.0025 s to check local time-discretisation sensitivity.
The maximum coefficient differences were about 4.7994e-11 for eta and
7.6403e-11 for psi. Initial/final reference mean eta agrees within 1e-12;
no mean correction is applied. This is not a long-time convergence study.

hos_time_check compares the GPU final eta and psi at every physical grid
point with the matched-step MATLAB reference, reports raw maximum absolute
and relative L2 errors and observed mean-elevation change. Predeclared
absolute thresholds are 2e-5 FP32 and 2e-11 FP64. The reference includes
nonlinear evolution and the mean potential; no fitting, shifting or scaling
is used.

Run with a newly generated reference:

    python run_p40_checks.py --time-reference <directory>/hos_time.bin

This does not validate the production q=3 masks, adaptive Cash-Karp stages,
an integrating-factor scheme, long-time energy/phase behavior or performance.

## Verified result, 27 September 2026

P40 snapshot `/root/hos-rhs-check-20260927T005948Z`; hashes and compact results
are in time_results.json, with the matching raw local artifacts/hos_gpu log.
The two time-integration checks and all eight existing CTest tests passed.

| Precision | eta max abs | psi max abs | eta relative L2 | psi relative L2 |
|---|---:|---:|---:|---:|
| FP32 | 7.04438e-7 | 3.58131e-6 | 1.22677e-5 | 1.41691e-5 |
| FP64 | 3.57353e-16 | 1.60982e-15 | 5.09021e-15 | 6.25526e-15 |

Mean eta changed by -7.29808e-10 (FP32) and -1.12757e-17 (FP64), without
mean adjustment. Those observations cover only 0.2 s. All 2048 grid points
were compared against the same-step independent MATLAB reference. The
64x32 evaluator+RK allocations were 321,536 bytes FP32 and 643,072 bytes
FP64, excluding context/library overhead. No GPU/CPU timing ratio is claimed.
