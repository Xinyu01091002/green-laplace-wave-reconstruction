# Remote OW3D time-resolution study

**Superseded planning record:** the user cancelled the dt/depth/steepness
queue and approved a compact kpd=1, Akp=.12 wavegroup/random-wave redesign.
Use [COMPACT_OW3D_DESIGN.md](COMPACT_OW3D_DESIGN.md) for the current scope,
actual remote paths, prepared inputs, timeout and pending random normalization.
The older resource audit below is historical; it is not a launch instruction.

Status on 2026-09-25: resource and output-format audit completed; no new OW3D
run launched. Previous Akp=0.12 results are committed and pushed at
`1038842c284a2d684f03766dfec4c6e774d380f0` on
`codex/unidirectional-time-series`. This plan does not change published GL
interfaces or OW3D equations. All new simulation data and MATLAB processing
will remain remote; there is no bulk-download or automatic deletion step.

The current proposal uses native **format 20**, a limited time window and
a narrow horizontal strip. The user accepts 20 samples per third-harmonic
period, interpreted here as the period at **3fp**; output every **0.2 s**
meets that criterion for the two audited cases. Full-domain EP remains a
sparse diagnostic and a documented alternative, not the production output.

## Native EP alternative and local kinematics

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

Native format 20 instead permits `xbeg xend xstride ybeg yend ystride
tbeg tend tstride`. It writes selected horizontal nodes and time levels
to one appended file per selection. Its stock Cartesian layout contains
three surface arrays (eta and its two gradients), volume phi, three
velocities and nine velocity gradients at all Nz levels. There is no
vertical-subset or eta/phi-only switch in this format's audited reader.
The user now explicitly proposes using this format for local output;
retain the native file and extract only eta and top-level phi for analysis.
Confirm top-level phi agrees with EP surface P at matching times in the
smoke run, including any effects of end-of-step filtering.

Neither output choice requires a source patch or compiler. Format 22 is
not the selected option: its off-grid interpolation is unimplemented for
this two-horizontal-dimensional setup in the audited source.

Evidence read on 2026-09-25:

- [Official annotated input](https://github.com/apengsigkarup/OceanWave3D-Fortran90/blob/master/examples/inputfiles/OceanWave3D.inp), data-storage line.
- [Official input reader](https://github.com/apengsigkarup/OceanWave3D-Fortran90/blob/master/src/IO/ReadInputFileParameters.f90), DATA STORAGE branch.
- [Official EP writer](https://github.com/apengsigkarup/OceanWave3D-Fortran90/blob/master/src/IO/StoreData.f90).
- [Official kinematics writer](https://github.com/apengsigkarup/OceanWave3D-Fortran90/blob/master/src/IO/StoreKinematicData.f90).
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

## Proposed queue and concurrency, 19:19 UTC update

The user now prefers independent serial OW3D jobs without modifying the
solver for OpenMP/MPI/GPU. A **job** below means one phase; a physical
parameter group contains four jobs at 0/90/180/270 degrees. This is a
proposed queue, not permission inferred to launch a large campaign or
evidence that any of these jobs has started.

| Priority | Group | Jobs | Purpose and preparation |
| --- | --- | ---: | --- |
| P0 | Short full-grid resource/output pilot | 1 temporary job | Measure initialization and first kinematics-output peak RSS; verify native record layout, top-phi parity and timing. Move the output window to early steps for this pilot only. |
| P1 | kh=1, spread label 25, Akp=.12, integration dt=.1 | 4 | Reproduce current high-steepness main-group eta22/psi22 and eta20 with .2 s saved sampling. |
| P1 | kh=5, spread label 15, Akp=.12, integration dt=.1 | 4 | Revisit the deeper case, especially eta20. This changes both depth and spread relative to P1 kh=1; it is not an isolated depth comparison. |
| P2 | Repeat kh=1/spread25/Akp=.12 with integration dt=.2 | 4 | Full four-phase time-step comparison on identical .2 s output times. |
| P2 | kh=1/spread25/Akp=.02 on the same grid/domain as P1 | 4 | Controlled steepness comparison. Requires consistently regenerated initial eta/psi from the same normalized first-order spectrum. Do not scale total nonlinear initial fields or silently reuse the old smaller-grid weak case. |

The proposed production queue has **16 phase jobs / 4 groups**, excluding
the short pilot. Prepare and fingerprint existing P1 input pairs first.
P2 weak inputs are not yet generated or validated. Existing kh=5 main-group
display limits are 211.979--260.021 s, so the proposed 120--360 s output
window covers both existing strong cases. The scientific scope remains
unchanged: sampled eta22 eligibility >= one-third of same-x centerline;
raw versus oscillatory eta20 kept separate; no third-order MF12 comparison.

After that first queue, useful separate controls are kh=2 at fixed spread,
a changed directional spread at fixed depth, and a changed *frequency*
bandwidth at fixed depth/spread/steepness. Each changes one input property
with the same domain/grid checks. Do not equate the existing spread labels
with frequency bandwidth or claim the mixed kh=1/kh=5 pair isolates depth.

Live audit at 19:19:27 UTC found 48 visible logical CPUs, load
0.24/0.18/0.12, 793099804672 bytes (738.6 GiB) available RAM and
907689644032 bytes (845.4 GiB) available disk. A short vmstat sample showed
idle CPUs and no swap-in/out. The current SSH session, user slice and
parent slice have `cpu.max=max`, `memory.max=max`, `memory.high=max`;
the effective CPU set is 0--47. These are guest-visible resources, not a
guarantee of dedicated physical cores. No OW3D process appeared in the
previous exact-name process check; no substantial compute load appeared
in this refresh.

Historical controllers scheduled 12--16 jobs on other grids, but no
comparable recorded peak-RSS measurement was found in the inspected
artifacts. Their schedules do not certify the capacity of this grid.
The audited 2051*515*10 extended grid has 10562650 volume nodes. Using
8-byte reals and 4-byte integers, the GMRES workspace at cap 55 is
4.486 GiB and the fine-grid cross-derivative index/weight tables are
4.250 GiB. Their **8.735 GiB subtotal is not total or peak memory**:
multigrid matrices, preconditioner assembly, RK/state/work arrays and
kinematics temporaries must also be measured. Narrow output does not
shrink the computational domain or these core allocations.

Propose a **600 GiB aggregate OW3D RAM budget**, leaving about 139 GiB
of the currently available memory for MATLAB processing, the OS and
headroom. Let M be the largest observed per-job peak RSS across setup,
time stepping and kinematics. Use `floor(600/(1.25*M))` as the memory
ceiling, subject to fresh availability and the separate CPU/I/O checks:

| Concurrent phase jobs | Meaning | Peak RSS per job required for this budget |
| ---: | --- | ---: |
| 8 | Two four-phase groups | <=60 GiB |
| 12 | Three groups | <=40 GiB |
| 16 | Four groups | <=30 GiB |

Start with the single resource pilot, then four jobs, and target eight
concurrent jobs if the measurements support it. Twelve or sixteen are
conditional expansion points, not verified capacities. Compare elapsed
time per step/total throughput and I/O wait as concurrency grows; more
logical CPUs used does not by itself imply faster completion. Measure
actual process thread counts and keep MATLAB/BLAS thread use bounded.
Do not count swap as additional calculation RAM or oversubscribe to 48
simultaneous jobs merely because 48 CPUs are visible.

At the current full-x/nine-row output design, each four-phase group costs
about 88.67 GiB including sparse EP. Four production groups therefore
reserve about **354.69 GiB** raw output, plus inputs, checkpoints and
processed products. With a proposed 100 GiB free-disk reserve, at most
eight such raw groups fit in today's free space before other overhead.
That is a storage limit on retained results, not a concurrency limit;
finishing a job does not free its files. All data remain remote, with no
automatic old-run deletion. Rebudget larger Nz, wider output strips or
longer windows separately.

No resource pilot or production simulation has been executed for this
assessment, so run duration and the optimal concurrency remain unmeasured.
The resource transcript is `artifacts/remote_campaign_design/capacity_resource_audit.txt`.

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

The user's 20-point target at 3fp gives
`dt_out <= Tp/(3*20)`: 0.22937 s at kh=1 and 0.20018 s at kh=5.
Output every 0.2 s therefore gives 22.94 and 20.02 samples respectively.
This is a peak-frequency criterion. At three times the tabulated 99.99%
frequency it gives only 17.77 and 16.79 samples; a stricter requirement
covering that band would instead need about 0.15 s. No universal 20-point
claim is made for broad spectra or later nonlinear spectral broadening.
Use 0.2 s as the first output choice and check the actual spectrum.
Integration and output are separate: a 0.1 s integrator can save every two
steps, and a 0.2 s integrator every step. Compare integration accuracy on
their common saved times; interpolating a coarse record adds no resolution.

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

These full-domain figures are alternatives for comparison, not the revised
storage request. The current center and all four qualified lateral probes
have the same old sampled envelope-peak time, 240 s, with the established
main-group display window 212.476--267.524 s. Propose **120--360 s** output
to preserve buffers for the Fourier operations and eta20. Stop integration
at 360 s instead of 480 s. Integration must still begin from the original
initial state; starting output at 120 s does not authorize restarting the
physics there without a valid checkpoint.

For native kinematics, Nz=9 plus one bottom ghost gives 10 saved vertical
levels. With 8-byte reals and 4-byte record markers, each time sample uses
`8*(3+13*10)*Nx_out*Ny_out + 16*8` bytes. All native extra fields are
included in these estimates. At 0.2 s, 120--360 s gives 1201 samples:

| Horizontal selection | Nodes | Four-phase kinematics (GiB) |
| --- | ---: | ---: |
| Complete x, 9 near-center y rows | 2049 x 9 | 87.790 |
| Local patch, approximately +/-2 peak wavelengths along x | 81 x 9 | 3.471 |
| Nine probes on the center x column | 1 x 9 | 0.0434 |

Prefer the complete-x strip for the first reference: existing first-sector
extraction applies a positive-kx spatial projection. A full x line at each
selected y preserves that operation; a short x patch or a few probes do
not. The same-x off-centerline nodes are y indices 253--261, covering
y=4429.6875--4570.3125 m. No y Fourier decomposition is required for this
particular projection. The full two-dimensional *initial* spectrum must
still be retained for the joint directional input, independent of the
reduced later output. The smaller patch can follow after a separately
validated local temporal/phase separator; it is not a drop-in replacement.

`check_strip_projection.m` checked the existing four-phase fields at 120,
240 and 360 s. Strip versus full-field positive-kx projection differed by
at most 3.395e-16 relative L2. At the five saved probes, maximum absolute
difference was 8.90e-16 m; relative errors at the extremely small tail
signals reached 2.24e-12. These are floating-point extraction checks, not
new physical-accuracy or kinematics-binary validation results.

Sparse full EP every 60 s from 0--360 s adds about 0.882 GiB for four
phases, so the preferred total is about **88.7 GiB**, before small logs,
checkpoints and extracted products. Kinematics still computes velocity
derivatives over the full numerical domain before writing the subset;
this storage saving is not a measured CPU saving. Time-window restriction
avoids that output work outside the window. Measure its overhead in the pilot.

## Bounded first campaign and validation

1. Keep the existing kh=1/spread25/Akp=.12 initial eta/psi pairs, domain
   22500 by 9000 m, 2049 by 513 horizontal nodes and Nz=9 unchanged.
   Locate the authoritative initial files, verify their time origin,
   boundaries and hashes, and freeze one input manifest. A previous output
   snapshot is not a substitute for an initial first-order spectrum.
2. Run only a short native-output smoke/pilot first, with a new run ID:
   request a delayed-start kinematics interval and overlapping EP output.
   Check header and record lengths, indexing and actual sample count,
   finite fields, eta/top-phi parity with EP, CPU/RSS and output cost.
   The stock path rewrites an ASCII restart every EP dump; the revised
   sparse EP cadence avoids doing that at every sample.
3. Compare one fixed phase at integration steps 0.2 and 0.1 s, both saved
   every 0.2 s. Complete phases 0/90/180/270 using the verified step. A
   single-phase check is a pilot, not four-phase harmonic convergence.
   Example parameter lines for `dt=0.1 s`, start time zero and end 360 s:

   ```text
   3601 0.1 1 0.0 1       <- time integration line, not adjacent to storage
   600 20 1 1            <- storage line: sparse EP, format, kinematics on, one file
   1 2049 1 253 261 1 1201 3601 2
   ```

   Physical node indices exclude ghost nodes; the writer adds their offset.
   Time levels are one-based: `(index-1)*dt + initial_time`. Here 1201 and
   3601 give 120 and 360 s and stride 2 gives 0.2 s output. For the 0.2 s
   integration control use Nsteps=1801, EP stride 300, time indices
   601--1801 and output stride 1. Validate the real binary before deployment;
   these are proposed fragments, not a newly executed input file.
4. Extract and phase-separate in MATLAB remotely, retaining the existing
   spatial separation convention and audited temporal sign. Process frames
   in a stream and save compact probe/spectrum products rather than a
   duplicate full space-time array. Retain native raw files. Reuse center and
   off-centerline probes, including dy=+/-52.734375 and +/-70.3125 m at
   x=11250 m, with same-x center y=4500 m. Apply the user's sampled eta22
   peak >= one-third centerline criterion independently of GL errors.
5. On the four-phase records, repeat the same analysis using longer and
   shorter *retained windows*, and 0.4 and 4 s subsets. Window truncation
   changes the frequency resolution and can affect eta20 and Hilbert/FFT
   edges. Compare the common core, not a cropped display alone. The 120 s
   start requires preserving the original absolute time in the directional
   prior: the old driver assumes records begin at its initial-spectrum time
   and cannot simply be passed a cropped record with the same unshifted
   coefficients. Use the declared phase evolution, never fitted alignment.
   If eta20 is not window-stable, extending the same strip to 0--360 s
   would cost about 131.65 GiB, still far below dense full EP output.
   Compare main-group waveform L2 and maximum absolute error, sampled
   extrema and spectrum, in physical units without alignment or fitted
   gain. Proposed convergence targets: below 0.5% relative L2 change in eta22
   and psi22, below 1% in eta33 and oscillatory eta20; include absolute
   errors for small signals. These are proposed convergence checks, not
   already achieved accuracies against the physical reference. If they
   fail, extend the record or reserve a finer short sampling check before
   making a sufficiency claim. Twenty points at 3fp is the user's sampling
   target; it does not replace an integration or window-convergence check.
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
The updated estimate function and strip check ran successfully in local
MATLAB R2022b; both Code Analyzer results were empty. The projection result
is in `artifacts/remote_campaign_design/strip_projection.json`, and local
kinematics sizes in `local_kinematics_storage.csv`. `git diff --check` passed.

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
