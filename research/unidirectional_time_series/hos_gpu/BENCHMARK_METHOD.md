# Matched CPU/GPU benchmark contract

The user requested additional tests, measured acceleration and memory use.
This benchmark compares the current degree-five prototype on one jfm92 host.
It is NOT a measured speed ratio against the production HOS-Ocean MPI solver.

## Same operator and execution sequence

generate_cpu_backend.py mechanically translates the four current kernel/RHS/RK
headers: CUDA launches become ordinary calls, outer per-cell kernel loops
become OpenMP parallel loops, and CUDA allocation/copy/FFT calls are mapped by
cpu_cuda.hpp to host allocations/copies and FFTW. No physics formulas, HOS
degree loops, support projections or RK stages are replaced or simplified.
CPU uses FFTW_MEASURE plans with 8 FFTW threads and 8 OpenMP threads; phases
are sequential, not nested FFTW/OpenMP regions. Both processes are restricted
to CPU IDs 0..7. GPU uses the same original CUDA headers and cuFFT on P40.
Both binaries use Release builds without fast-math. This is a competent
threaded backend for the same prototype, not an assertion of best possible
CPU tuning, MPI scaling or optimal GPU kernels. Spectral multipliers are
currently recomputed by both implementations rather than cached.

## Fixed inputs and matrix

- Grids 512x256 and 2048x1024; M=5.
- Retained input |mx|<=20, |my|<=10; complete degree-five support fits strictly
  below working-grid Nyquist on both grids. No physical product modes alias.
- 273 positive-x modes (mx=8..20, my=-10..10) plus their conjugates.
  Fixed smooth weights and deterministic phases, initial eta RMS 0.01 m.
- Lx=90 m, Ly=60 m, h=1.3 m, g=9.81 m/s^2. Initial psi follows the linear
  progressive-wave coefficient relation for the same eta spectrum.
- Physical inputs are generated once using double-precision FFTW before the
  measurement boundary, then cast to the requested precision. CPU/GPU inputs
  use the same generator and mode/phase constants.
- FP32 versus FP32 and FP64 versus FP64 are the primary speed comparisons.
  GPU FP32 is additionally compared numerically against CPU FP64.

## Timing boundaries

One-shot time starts after prepared input and device discovery. It includes
operator/RK allocation, plan creation, input upload/copy, first RHS evaluation
and output download/copy; it excludes executable startup and input generation.
CPU MEASURE and cuFFT plan creation are both charged here.

Steady RHS: one additional warm-up, then three timed evaluations; every GPU
sample ends with cudaDeviceSynchronize. No upload/download, allocation or
reference calculation occurs inside this timed interval. Median and all
individual samples are saved; these measure a reused grid/support/plan.

RK4: one warm-up step, then three measured single steps from the same reset
initial state, dt=0.005 s. State reset and transfer are outside timing, but
all four RHS stages, stage updates and retained-band projections are inside.
This is per-step timing, not a long simulation or adaptive-integrator result.

## Correctness gates

Save full-grid eta_t, psi_t, eta after one RK4 step and psi after one step.
An independent streaming C++ comparator checks unaligned, unscaled fields.
Predeclared relative L2 tolerances: same-precision FP64 1e-9; same-precision
FP32 and GPU FP32 against CPU FP64 5e-4. All cells must be finite. Results
are separate from the earlier independent MATLAB polynomial references.
Timing ratios are reportable only after these numerical checks pass.

## Memory accounting

- Process peak host RSS is measured with Linux getrusage after computation
  and output. It includes host inputs/outputs, library/plan overhead and any
  preparation peak; it is not merely counted solver-array bytes.
- GPU memory is sampled read-only through NVML every 0.05 s. Report the
  device-wide peak and increment above the immediately preceding idle value.
  This captures context/library usage in addition to arrays, but is a sampled
  peak and includes any other device activity. Hardware snapshot and no-job
  preflight constrain that interpretation; it is not an exclusive reservation.
- Also record exact declared array/workspace allocation for transparency.
  CPU FFTW internal plan memory is not included in that declared-byte count.
- CPU RAM and GPU VRAM are different resources; no total-system memory-saving
  ratio is inferred by comparing the two columns alone.

Temperature is observed, not controlled. No power, clock, ECC, thermal,
firmware or driver setting is changed. All remote files use a new hash-bound
snapshot; existing calculations and historical results are untouched.
