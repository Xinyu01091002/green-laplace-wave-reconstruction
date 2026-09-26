# HOS-Ocean deployment and first pilot

**Update 2026-09-25 22:35 UTC:** A separately identified source build now
resolves the external-surface IO precision limitation. Three matched 2 s
pilots passed; raw t=0 eta/psi errors are below 9e-16 relative L2.
See [PRECISION.md](PRECISION.md) for the patch, binary and evidence.
The original deployment and six-digit pilot below are retained as history.

Verified on 2026-09-25 (22:18 UTC). Remote Linux deployment and a short
full-grid wavegroup pilot completed. This establishes runnable deployment,
input/output conventions and finite short-time output, not physical accuracy,
long-time stability, high-order reference accuracy or parallel scaling.

## Deployment

Official source: https://gitlab.com/lheea/HOS-Ocean, tag v2.1.0,
commit `4deb3b4913d993c4e6ea16f736e5fc5792e14f12`.
Official package: `HOS-Ocean-linux-x64-hdf5-v2.1.0.tar.gz`.
Archive SHA-256: `2bd3ee01ad8fe463ea26bb50634ac0eee9772616d82525ca1730107bea8fadb6`.
Solver SHA-256: `117f0a485538cb723f70ae8dd63fadcdf7c0f58ddf632d3d9ff9471a26eafdcf`.

Remote root:
`/home/lxy/green-laplace-unidirectional-time-series/artifacts/hos-ocean-v2.1.0`.
The `source/` checkout matches the local pinned source. The executable is
`install/bin/HOS-Ocean`. Ubuntu packages liblapack3, libblas3, libsz2, libaec0
and time were downloaded and unpacked into `deps/`, without system package
installation. `deployment-sha256.json` records package hashes.

For each solver process set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and
MKL_NUM_THREADS to 1. Its private LD_LIBRARY_PATH is the colon-separated
`deps/usr/lib/x86_64-linux-gnu`, and its `lapack` and `blas` subdirectories.
Do not apply this library path globally or to MATLAB.

The host is Ubuntu 24.04 x86_64 with 48 visible logical CPUs. Before starting,
available memory was 472 GiB and filesystem free space 843 GiB. Sixteen
OW3D processes each used approximately one CPU. They remained running after
the HOS pilot; no frozen OW3D case or controller was modified.

## Executed cases

1. `official-smoke-v2`: official Linear_Regular_RZ example, with only duration
   reduced from 200 to 0.2 periods. Exit 0, 0.06 s wall time, 9216 KiB peak RSS.
   This is a startup smoke test, not the upstream full CTest suite. Earlier
   `official-smoke/` contains the failed launcher attempt before installing
   private GNU time; its solver was never started.
2. `wavegroup-phase000-pilot`: frozen MF12 order-2 wavegroup phase 0,
   1024 x 256 unique periodic points, domain 11260.188722544061 x
   4504.075489017624 m, depth 35.842293906810035 m, g=9.81 m/s2.
   M=5, qx=qy=3 (partial dealiasing), RK tolerance 1e-8, no ramp,
   no breaking or added dissipation. Physical output interval 0.2 s,
   duration 2 s. No volumic output or modal output requested.

Pilot input:
`/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z/inputs/wavegroup/initial_fields.mat`.
Input SHA-256: `14362441397d27a4fe815c5da8325f74765d041f79e1cffc8c5ecca271025c5e`.

Measured solver wall time: 126.77 s, including initialization and surface
output. Peak RSS: 911848 KiB (0.870 GiB). One thread was observed. These are
one short run's measurements alongside 16 OW3D jobs, not a production runtime
forecast or matched-accuracy speed comparison.

`status.json` reports `short_pilot_io_validated`, solver exit 0.
`Results/3d.dat` contains all 11 expected frames at 0:0.2:2 s; MATLAB checked
all eta and true-surface-psi values for finiteness and checked coordinate order.
`surface_strip.mat` retains x and y indices 125:133 and physical times remotely.
`validation.json`, `input-audit.json`, `run.log`, and `resources.txt` hold evidence.
No raw fields were downloaded locally.

## Precision and import behavior

The unmodified v2.1.0 external input reader uses 67 header lines, a fixed
zone record, and two ES12.5 columns. X varies fastest; local MATLAB arrays
are transposed from the supplied ny-by-nx storage. Input units are dimensional.
The ordinary free-surface output also uses ES12.5. The HDF5 switch does not
convert this ordinary eta/psi surface output into HDF5.

Input quantization relative L2: eta 1.5155931427e-6, psi 1.5404442727e-6.
Maximum absolute input quantization: eta 4.99449e-6 m, psi 4.99955e-5 m2/s.
Quantization introduced means 1.79065e-10 m and -6.99470e-9 m2/s.
Pinned `sources/HOS/HOS-ocean.f90:386-390` resets both mean modes when the
initial eta mean is nonzero. This behavior explains the raw t=0 differences
against the quantized file: relative L2 2.68803e-10 and 3.01812e-10.
The initial exact-import assertion failed and is preserved in
`validate-pilot.log` and `validate_hos_pilot_v1.m`. The final validator checks
the documented zero-mode operation against per-value ES12.5 rounding bounds.
It reports raw errors and applies no fitted correction to saved results.
Raw t=0 error against original input remains approximately 1.5e-6 relative L2.

Before production high-order comparisons, implement and verify sufficient
input/output precision in a separately identified source build. Preserve the
unmodified binary and this run as evidence. The 220 s HOS production campaign,
random realizations and four-phase HOS families have not yet been launched.

## Scripts

- `prepare_hos_pilot.m`: MATLAB exporter, fixed 2 s phase-0 pilot, fails if
  output exists; records quantization rather than changing input amplitudes.
- `run_pilot.py`: orchestration and hashes for the exact remote deployment;
  runs MATLAB and then HOS, keeping its private libraries isolated.
- `validate_hos_pilot.m`: MATLAB full-field finite-output and import audit;
  keeps data remote. This is not a physics implementation or accuracy test.

The executed scripts were copied as hashed snapshots. Deployment/failed
attempt artifacts and official source remain under ignored `artifacts/`.