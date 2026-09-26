# Two-amplitude JONSWAP R4-GL HOS campaign

User approved starting low then high on2026-09-26 and requested a longer record.
Interpretation: increase record duration to80 peak periods, keep Tp unchanged.

Active remote root:
`/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-r4gl-80tp-20260926-v2`.
Screen:`hos-jonswap-low-high-80tp-v2`. Whole queue:`pipeline-status.json`;
per-family:`low/status.json`, `high/status.json`; preflight:`short-high/status.json`.

## Frozen physical and numerical inputs

- kp Hs/2=.02 then.12; Hs1.4336917563 and8.6021505376m, defined from initial
  linear spatial variance. No total-field or single-probe renormalization.
- kp=.0279/m, h35.84229391m, kph1, g9.81, Tp13.76199709s.
- JONSWAP gamma3.3, initial design support.5--2.5fp, energy Gaussian direction
  sigma17.6777deg around+x, original directions+/-60deg. Preserve the11128
  parents selected by cumulative99% modal variance, actual bandwidth25.274%.
- Same seed20260925 and complex phases in both amplitudes. Four global shifts
  0/90/180/270deg; no independent randomization between those four runs.
- Domain50x20lambda =11260.1887x4504.0755m; periodic unique grid1024x512.
  M5/qx=qy3 partial dealiasing; adaptive Cash-Karp5(4), absolute tol1e-12;
  Ta0, no extra damping or breaking model.
- Duration1100.8s, rounded to a multiple of.4s close to80Tp; output.2s,
  expected5505 probe samples. Score10--70Tp and show the full record.
- Same physical five probes as before: x=Lx/2; y=Ly/2 plus0,+/-52.7821,
  +/-70.3762m. Native indices now[513,257/251/263/249/265].
- Per family4 phases x8 MPI ranks on CPUs8--39. Only one HOS family at a time.
  Existing8 OW3D jobs untouched. CPU40 for sequential MATLAB preparation/GL.
- Initial full eta/true-surface-psi retained; production writes probes only,
  not dense evolving spatial fields. This differs from the optional sparse
  snapshots discussed in the design; no nonexistent spatial output is claimed.

## Hybrid initialization

11 from exact declared analytic parents;20 from current pure SPARK R4 eta/psi;
22 from current pure SPARK GL fixed-FFT second-order graph. No initial31/33.
The source remains a private runtime snapshot in ignored artifacts/remote
storage, not vendored into public code. Source manifests contain hashes.

The already examined R4 numerical residuals3.332e-8/1.786e-8 are handled in
the RUN-LOCAL copy:1e-8..1e-7 warns and is recorded; above1e-7 still stops.
The original SPARK code, formulas, low-output repair, zero-sector convention,
and physical fields are unchanged. This policy is explicit, not a claim that
the original1e-8 gate passed. The prior MF12 differences1.879/1.774% remain.

GL12->GL16 relative changes: eta7.5063e-5, psi3.4299e-5. GL16 selected by the
declared0.5% refinement condition. This is quadrature consistency, not an
independent MF12 accuracy claim for22. High and low initialize from the same
fields by linear1/6 and quadratic1/36 scaling, preserving the exact order2
homogeneity and global-phase symmetry. Initial first-sector recovery checked.

## Execution order and resource handling

1. Generate both initial families and export8 y-slabs each.
2. High family short run ONLY2 physical seconds (four concurrent phases).
   Require all initial probe comparisons, finite11 records, zero exits.
3. Full low HOS, then full high HOS. No second full benchmark configuration.
4. Low GL time-series comparison, then high GL comparison, using the same
   reference/prediction band rules and direction resolutions3.75/1.875deg.
   Ignore EXACT-zero amplitude direction bins before the unchanged pair kernel;
   no tolerance-based pruning or spectrum tuning is introduced.
5. One completion/failure email via existing remote recipient configuration.
   Email after the full workflow; mail acceptance is not inbox confirmation.

HOS aggregate RSS guard16GiB. MATLAB initialization32GiB guard and GL128GiB
guard, minimum host available64GiB. No other jobs are terminated. Timings and
memory are recorded separately for initialization, preflight, each phase and
GL. The initial live high preflight consumed about7.7GiB aggregate and all
initial probe errors were below1e-13m. Broadband high adaptive stepping is much
slower than the previous weak narrowband case; old timing estimates must not
be presented as measured completion forecasts.

## Deployment notes

v1 stopped in preparation because CRCRLF in a copied MATLAB continuation line
was invalid; no HOS ran in v1. v2 normalizes the source newline format and
preserves original/deployed hashes. Its completion shell was normalized to LF;
`sh -n` passed. The additional exact-zero postprocessing change and final mail
script hash are recorded in postprocess-runtime-patch.json. No existing run
or public SPARK source was overwritten. Read current status files to distinguish
preflight, low production, high production and completed results.

## Verified production launch

Short high preflight completed21:25:59 UTC: all four phases passed initial
probe checks, exited0 and produced11 finite samples through2s. Slowest wall
532.805s; aggregate peak RSS8097080KiB (7.722GiB). This is only a2s check,
not a second full high campaign.

Low full production started2026-09-26 21:26:00 UTC. All four initial-probe
checks passed below1e-14m; nonzero advancement to.2s was observed. High is
queued automatically after low. Early low phase0 output cost11.64s per.2s
suggests about18h if speed persists; high short-run extrapolation is about
81h. Both are early estimates, not bounds or promises. Postprocessing uses
GL16 and full-window plots; rank/display runtime patches retain SHA256.
