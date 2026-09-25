# Larger-steepness OW3D discovery

Read-only inventory under `C:/Research/VWA/VWA time series/unidirectional`.
Directory names alone are not evidence of completed output.

- `testboundtail/T_init-40_Tp_Alpha_8.0_Akp_006_kd0.5_phi_{0,90,180,270}`:
  all four phases have every EP file at steps 0:4:2400 (601 ordinary snapshots
  per phase, plus EP_99999). Source dt=0.400334 s gives output dt=1.601336 s;
  h=17.92115 m, kp=0.0279 rad/m, Tp=17.667179 s, nominal Akp=0.06.
  The phase-zero log reports 2401 completed steps. File continuity and metadata
  were checked; full field finiteness and reconstruction have not yet been run.
  This is the next candidate for a shallow-water, higher-steepness comparison.
  Because Alpha=8 differs from the existing Alpha=1 pilot, also use the
  available Alpha=8, Akp=0.02 case at kh=0.5 as the matched lower-amplitude case.
- `timeseriesdata5`: directories advertise Akp=0.12 at kh=1/5 and Alpha=1/8.
  The four-phase kh=1 families for both Alpha values contain input/init/readme
  files but no EP output. Do not claim that this directory supplies completed
  Akp=0.12 time series. No rerun or remote search was initiated.
- `data/kd1.0_Alpha_1.0_Akp_012_phi_{0,90,180,270}`:
  all four phases contain steps 0:10:1800 (181 snapshots each). The source
  integration dt=0.4 s gives 4 s output sampling; its Nyquist is 0.125 Hz.
  The carrier implied by kp=0.027947 rad/m and h=35.78259 m has a second
  harmonic above that Nyquist. These spatial snapshots are therefore not
  directly suitable for this temporal eta22 comparison at the original rate.
  They may support spatial comparisons, which are outside this immediate task.

The user now prefers main-wave-group plots. `plot_ow3d_main_group` creates
additional PNG/PDF views from saved results, centred on the common first-order
analytic-envelope maximum, with +/-2Tp limits. It keeps elapsed time on the
axis, makes no shift of the data or prediction, and retains full-window files.
This is a display choice; earlier full-window metrics remain unchanged.
