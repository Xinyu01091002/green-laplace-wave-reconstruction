# Existing directional OW3D data inventory

Discovery only: no directional GL algorithm, new OW3D run, project rename,
large remote transfer or historical-source modification was performed.

## Local VWA directional fields

Root: `C:/Research/VWA/VWA Unidirectinal/Directional`.

| Folder | kh | Directory spreading parameter | Akp | Four-phase groups | Ordinary snapshots per phase | Output interval |
|---|---|---|---|---:|---:|---:|
| test1 | 1, 2, 5 | 5, 15, 25 degrees | 0.02 | 9 | 91, steps 0:10:900 | 4 s |
| test6 | 1 / 5 | 25 / 15 degrees, respectively | 0.12 | 2 | 121, steps 0:10:1200 | 4 s |

All 36 test1 and all 8 test6 phase directories have every expected ordinary EP
snapshot and an OceanWave3D.end file. Four phases are 0/90/180/270 degrees.
The integration step is 0.4 s; elapsed-time coverage is 0--360 s and 0--480 s.
These are actual two-dimensional output fields, not just initialization files.

Representative shallow-water case paths:

```text
test1/kd1.0_spread_25_Akp_0.02_phi_shift_{0,90,180,270}/EP_00300.bin
test6/kd1.0_spread_25_Akp_0.12_phi_shift_{0,90,180,270}/EP_00400.bin
```

The test1 input grid is 1025 x 257; the sample EP header is 1027 x 259.
The test6 input grid is 2049 x 513; the sample EP header is 2051 x 515.
This includes additional output points; do not silently reinterpret either
header as a periodic FFT grid. EP stores x, y, eta and potential arrays. The
same binary layout was checked in the unidirectional work, and representative
directional field checks are saved in `artifacts/directional_data_inventory/`.
Only representative snapshots, not every file, are numerically read/checked.

The directory `spread` values are reported as labels. A current generator
uses a Gaussian angular weight with the spreading parameter as its width,
but that edited generator alone does not prove the exact historical generation
settings of every dataset. Read source metadata before assigning that semantic.
Likewise the legacy README's Tp=12 should not override finite-depth dispersion
computed from the stated kp=0.0279 and depth.

The 4 s interval gives a temporal Nyquist of 0.125 Hz, insufficient for the
carrier's second harmonic in these cases. It is nonetheless usable for fixed-
time two-dimensional field comparisons. `processed_eta22/` and `processed_eta33/`
contain earlier MATLAB field products, but raw four-phase fields remain the
preferred source for a fresh comparison.

`inversetest1` advertises Akp=0.16 and `inversetest2` Akp=0.40 at kh=1/5.
None of the sixteen inspected phase directories contains EP output; these
must not be counted as completed directional simulation datasets.

## Remote crossing-wave fields, verified live

Use the existing SSH alias `60.188.112.99:60093`. No authentication files
were read or copied.

Completed kh=1 root:
`/home/lxy/gl-r2-ow3d-runs/gl-r2-ow3d-kpd1-hotstart-20260827T135129Z`.

Cases: `cases/glr2_kpd1_cross{0,60,120}_sigma30_phi{0,90,180,270}`.
OW_readme reports Akp=0.18. The run-level status is
`completed_12_of_12_kpd1_hotstarts`, exit code 0, completed_count=12,
failed_count=0. Every case is recorded at 300.00000000000563 s with exit 0;
all twelve end files exist. A representative end-file header was checked.
Checkpoint dimensions are 1025 x 513 with repeated periodic endpoints; the
documented physical FFT grid is 1024 x 512. These paired eta/psi ASCII end
fields are confirmed fixed-time data. EP time-series files were not found at
the inspected case depth; a regular output interval is not claimed.

The initial fields use historical GL/R2 construction, according to their
metadata. Do not treat this as a new no-Stokes release validation or import
its historical executors. The four-phase fields are usable as data subject
to a separate grid, finiteness, initialization and domain audit.

Review-required kh=5 root:
`/home/lxy/gl-r2-ow3d-runs/directional-ow3d-kpd5-glscale-repaired-20260904T201519Z`.
The live run-level status is `failed_kpd5_repaired`, exit 1, failed_count=12.
However all twelve summary rows have solver exit 0 and final time
300.30000000000564 s; a representative log ends `JOB IS COMPLETE` and the
checkpoint agrees with that time. A strict final-time gate may explain the
run-level failure, but this has not been established. Do not call these
validated completed cases or infer numerical blow-up from status alone.

## Recommended starting point

Use the local test1 kh=1, spread=25, Akp=0.02 four-phase family for a first
fixed-time directional field test, then test6 kh=1, spread=25, Akp=0.12.
Direction information must come from the declared directional input/spatial
spectrum or adequate observations; the single-point eta1(t) workflow cannot
simply assume every temporal frequency has a unique propagation direction.
The user subsequently accepted knowledge of the initial linear spatial
wave-number spectrum together with the observed probe eta1(t). See
[INPUT_ASSUMPTION.md](INPUT_ASSUMPTION.md) for their distinct roles and the
directional allocation assumption still to be established. Agreement of a
linearly propagated initial spectrum with the observed record is not a gate.
