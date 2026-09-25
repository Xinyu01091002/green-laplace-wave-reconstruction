# Remote OW3D time-resolution study

Status on 2026-09-25: resource and output-format audit completed; no new OW3D
run launched. Previous Akp=0.12 results are committed and pushed at
`1038842c284a2d684f03766dfec4c6e774d380f0` on
`codex/unidirectional-time-series`. This plan does not change published GL
interfaces or OW3D equations. All new simulation data and MATLAB processing
will remain remote; there is no bulk-download or automatic deletion step.

## Native eta/phi output: no kinematics or source patch required

The documented input line is `StoreDataOnOff formattype`. Use a positive
integer stride and `formattype=1`, for example:

```text
1 1       <- save EP every integration step, unformatted binary
```

The existing directional test6 inputs already use `10 1`. With `dt=0.4 s`,
this means one snapshot every 4 s. The reader sets `iKinematics=0` when
the format is not one of the kinematics formats 20/21/22. No kinematics
range lines are added with format 1. If adapting an input that used those
formats, remove their range lines as well to preserve input alignment.

`StoreData` writes `EP_XXXXX.bin` with three Fortran sequential records:
`Nx,Ny`; the two coordinate arrays; and `E,P`. Here `E` is eta and `P` is
the free-surface potential psi=phi(x,y,z=eta,t), not a vertical volume.
No velocity, velocity gradient or volume potential is saved in EP.
The stock path also updates the surface restart file `OceanWave3D.end`.
EP covers the full horizontal surface, including ghost nodes in this
configuration; it is not a native point-probe file. MATLAB can extract
arbitrary on-grid probe records afterwards, remotely.

The earlier suggestion that eta/phi-only output requires a new Fortran
writer was incorrect: it conflated selecting variables with selecting
spatial probes. Neither a compiler nor a solver modification is required
for this plan. A probe-only storage optimization is not part of this task.

Evidence read on 2026-09-25:

- [Official annotated input](https://github.com/apengsigkarup/OceanWave3D-Fortran90/blob/master/examples/inputfiles/OceanWave3D.inp), data-storage line.
- [Official input reader](https://github.com/apengsigkarup/OceanWave3D-Fortran90/blob/master/src/IO/ReadInputFileParameters.f90), DATA STORAGE branch.
- [Official EP writer](https://github.com/apengsigkarup/OceanWave3D-Fortran90/blob/master/src/IO/StoreData.f90).
- Read-only local source: `C:/Research/OW3D_benchmark/OceanWave3D-Fortran90-master/`.
  Local `StoreData.f90` SHA-256 is
  `2bca0a3ee2cd6f44fe6e65f6af4b10e79fda41fba33771c6bcf3b06b1e953514`;
  `ReadInputFileParameters.f90` is
  `e4462b1107a2f0e6c2dc3b9e1c60dfddf90fb18666ae676b56bea02e6280ae9d`.

The remote executable has been fingerprinted, but its build is not yet
mapped to an exact source revision. A short native-output smoke run must
confirm this executable's format, record count and timestamps before a
production campaign; documentation/source inspection is not that runtime test.

## Resources checked live

SSH alias: `60.188.112.99:60093`, user `lxy`, host `jfm93`.
Snapshot at 2026-09-25 18:45:56 UTC:

| Resource | Observation |
| --- | --- |
| CPU | 48 visible logical CPUs in a VMware guest; load 0.22/0.19/0.18 |
| RAM | 747.6 GiB total, 738.6 GiB available |
| Storage | 907689693184 bytes available, approximately 845.4 GiB |
| OW3D | `/usr/local/bin/ow3d`; no running `ow3d`/`OceanWave3D` process found |
| MATLAB | `/home/lxy/Desktop/matlabr2026a/bin/matlab`, R2026a Update 2; batch launch and license test passed |
| Other tools | WolframScript, Python3, g++, CMake and Git present |
| GPU | `nvidia-smi` could not communicate with the driver; this CPU plan does not depend on it |

The OW3D binary SHA-256 is
`36bd612d2655c26daa3a3146d7d3d291f0864d8cd3f23c3d2dd6992f7a26f614`.
Its linked libraries show no OpenMP runtime; no parallel speedup or run-time
estimate has been measured. Start with one pilot, then at most four
independent phase runs after measuring elapsed time, RSS and disk writes.
The old MATLAB process observed in the audit is a zero-RSS zombie, not an
active calculation; it and all historical runs were left untouched.

The intended project directory
`/home/lxy/green-laplace-unidirectional-time-series` is absent and has not
been created. On deployment, use an independent clone of
`https://github.com/Xinyu01091002/green-laplace-wave-reconstruction.git`,
pin the research commit, and use a new `results/ow3d/<run-id>/` directory.
Recheck the destination before creation; never overwrite an existing run.

## Sampling and storage

`design_remote_sampling.m` evaluated the already-local initial positive-kx
spectra for the kh=1/spread25/Akp=.12 and kh=5/spread15/Akp=.12 cases.
These are initial-spectrum planning indicators, not a cutoff to impose on
the nonlinear result or evidence of full time convergence.

| Case | Peak period (s) | Frequency at 99.99% cumulative initial spectral energy (Hz) |
| --- | ---: | ---: |
| kh=1, spread label 25 | 13.7620 | 0.09379 |
| kh=5, spread label 15 | 12.0106 | 0.09929 |

Old 4 s sampling has a 0.125 Hz Nyquist frequency, below 2fp and 3fp in
both cases. Spatially separated values at the saved instants remain usable,
but a densely interpolated curve does not recover the missing temporal
information. This particularly limits temporal separation of higher harmonics.

Recommend `dt=0.1 s`, EP stride 1, for the first fine run. This gives about
20--21 samples per period at five times the tabulated frequency. At
`dt_out=0.2 s`, there are still about 17--18 samples per period at three
times that frequency, making 0.2 s a useful coarser comparison. Neither is
certified sufficient before testing. Integration step and saved sampling
must be distinguished: interpolation of a 0.2 s run is not a 0.1 s result.

For the existing strong-case EP dimensions 2051 by 515 (including ghosts),
one uncompressed EP file is `32*Nx*Ny+32 = 33800512` bytes: coordinates
are repeated in every file. For 0--480 s inclusive:

| Saved interval (s) | Files per phase | Four-phase EP storage (GiB) |
| ---: | ---: | ---: |
| 4.0 | 121 | 15.24 |
| 0.4 | 1201 | 151.23 |
| 0.2 | 2401 | 302.33 |
| 0.1 | 4801 | 604.53 |
| 0.05 | 9601 | 1208.93 |

These estimates exclude restart files, logs and processed products and
assume no compression. Four 0.1 s phases plus one full 0.2 s control phase
use 680.11 GiB, leaving about 165 GiB at the audited free-space level.
Two complete four-phase batches at 0.1 and 0.2 s would need 906.85 GiB
before overhead and do not fit. Do not start a full depth/bandwidth sweep
or duplicate the full EP archive into a second format. Recheck free space
and actual pilot I/O cost before launch; reduced output cadence is an
explicit future choice, never silently substituted.

## Bounded first campaign and validation

1. Keep the existing kh=1/spread25/Akp=.12 initial eta/psi pairs, domain
   22500 by 9000 m, 2049 by 513 horizontal nodes and Nz=9 unchanged.
   Locate the authoritative initial files, verify their time origin,
   boundaries and hashes, and freeze one input manifest. A previous output
   snapshot is not a substitute for an initial first-order spectrum.
2. Run only a short native-output smoke/pilot first, with a new run ID:
   check eta/psi record layout, no `Kinematics*` output, initial/final
   time and step indexing, finite fields, CPU/RSS and output cost. The
   stock path rewrites an ASCII restart every EP dump, so disk capacity
   alone does not predict its throughput.
3. Compare one fixed phase at integration/output 0.2 and 0.1 s over the
   original time window. Then, if stable and within the storage budget,
   complete phases 0/90/180/270 at 0.1 s. For a zero start time,
   `Nsteps=4801`, `dt=0.1`, output line `1 1` gives 480 s;
   the 0.2 s control uses `Nsteps=2401`. Check the actual initial-file
   start time instead of assuming it. A single-phase integration comparison
   is a pilot, not full four-phase harmonic convergence certification.
4. Extract and phase-separate in MATLAB remotely, retaining the existing
   spatial separation convention and audited temporal sign. Process frames
   in a stream and save compact probe/spectrum products rather than a
   duplicate full space-time array. Retain raw EP files. Reuse center and
   off-centerline probes, including dy=+/-52.734375 and +/-70.3125 m at
   x=11250 m, with same-x center y=4500 m. Apply the user's sampled eta22
   peak >= one-third centerline criterion independently of GL errors.
5. On the fine four-phase records, repeat the same analysis using 0.2,
   0.4 and 4 s subsets. This tests output sampling at fixed integration.
   Compare main-group waveform L2 and maximum absolute error, sampled
   extrema and spectrum, in physical units without alignment or fitted
   gain. Proposed sampling targets: below 0.5% relative L2 change in eta22
   and psi22, below 1% in eta33 and oscillatory eta20; include absolute
   errors for small signals. These are proposed convergence checks, not
   already achieved accuracies against the physical reference. If they
   fail, reserve a finer run before making a sufficiency claim.
6. Four phase sectors are not exact perturbation orders. In particular,
   the zero phase sector can contain fourth harmonics; denser time sampling
   helps frequency separation but does not isolate eta20 by itself.
   Preserve stationary nonzero-spatial-K contributions in raw eta20;
   report any common temporal-DC projection separately. Consider eight
   consistently generated phases only after this output/time audit,
   not by interpolating nonlinear initial conditions. No MF12 third-order
   comparison is introduced.

Remote MATLAB release tests and the relevant existing directional
reconstructions must pass before using its R2026a results for new claims.
This task checked MATLAB startup, not those remote regressions. Higher
temporal resolution also does not establish spatial/Nz convergence.
After this first group, choose one new depth or bandwidth control at a
time with a fresh storage budget. No automatic old-data cleanup is authorized.

## Provenance and commands actually run

The local estimate reads `results/directional_sweep/*/extracted.mat`, which
already existed locally; no remote fields were fetched. The two input hashes
(kh1 then kh5) are:

```text
aebb4734405e25cd5a2a10d1961fedbe764480ea083401704362794020ffc6d4
b1e2291e82e32552dcae4a88156d491950490b50d231e8591939d3588aab5624
```

Commands used include `git status --short`, `git push origin
codex/unidirectional-time-series`, `git ls-remote origin
refs/heads/codex/unidirectional-time-series`; local `Get-FileHash` and `rg`;
and MATLAB `design_remote_sampling`. Estimate JSON, CSV, log and the live
resource transcript are under ignored `artifacts/remote_campaign_design/`.
The updated estimate function ran successfully in local MATLAB R2022b;
its Code Analyzer result was empty. `git diff --check` also passed.

SSH used `-o BatchMode=yes -o ConnectTimeout=10 60.188.112.99:60093` for
`nproc`, `lscpu`, `uptime`, `free`, `df`, exact-name `pgrep`, `command -v`,
`sha256sum`, `file`, `ldd`, bounded directory checks and:

```sh
timeout 60s /home/lxy/Desktop/matlabr2026a/bin/matlab -batch "disp(version); disp(computer); disp(license('test','MATLAB'));"
```

No OW3D run, compiler installation, source patch, remote clone, raw-data
transfer, old-process termination or historical-directory modification was
performed. Future manifests must record repository commit, input and binary
hashes, actual coordinates/phases/units/time origin, integration/output steps
and actual ending time. Require successful exit, finite output and a final
time check with rounding tolerance; a detached shell's existence is not
evidence that a simulation is running or complete.
