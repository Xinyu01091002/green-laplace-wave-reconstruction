# Production-size JONSWAP right-hand-side check

Scope: the existing low/high JONSWAP campaign's phase-zero initial states,
1024x512 base grid and 2048x1024 q=3 working grid, M=5. This is a frozen
initial-state operator check, not a new simulation or a four-phase validation.

## Input provenance and reference

Source campaign: hos-jonswap-r4gl-80tp-20260926-v2 on 93. Low/high labels are
kp*Hs/2=.02/.12, h=35.842293906810035 m, kp*h=1, 11128 retained free parents,
JONSWAP gamma=3.3 and the saved R4(20)+GL16(22) initialization. The full eta
and true surface potential in the saved inputs are used; they are not treated
as purely first order, rescaled or regenerated.

The eight saved rank partitions of each phase-zero input are assembled into
a fresh one-rank initial file by copying the data lines without numerical
conversion. Each original file is hashed before and after reading. Only the
zone header changes from the local partition dimensions to 1024x512.

The reference executable reuses unchanged production numerical module objects.
For this test, the main-program probe is inserted AFTER the original initial
condition loading, nondimensionalisation, retained-band handling and existing
initial-volume/mean-potential handling, but BEFORE time integration. The
original initialization path therefore defines the actual state supplied to
both solvers. No new mean correction is added by the GPU checker.

The exporter writes the same eight fields as OCEAN_ALIGNMENT.md: expanded
eta/psi, base eta/psi, nonlinear eta/psi RHS and full eta/psi RHS. The reader
header stores int32 dimensions and FP64 parameters/data. The official CPU
reference uses the MPI build with one rank, bound to CPU 40, and exits after
one RHS evaluation. Existing CPU campaigns are not modified.

## Precision gates fixed before running

FP64: relative L2 <=1e-8 AND max absolute error <=2e-12 + 1e-8*reference max.
FP32: relative L2 <=5e-4 AND max absolute error <=2e-6 + 5e-4*reference max.
Both criteria apply separately to all six checked fields (expanded inputs,
nonlinear RHS and full RHS). Zero-reference fields use the absolute gate.
Errors use the solver's nondimensional units. All values must be finite.
These are preliminary implementation gates, not a long-time error budget.

FP32 uses the same FP64 saved reference, cast only at the GPU input boundary;
its measured discrepancy includes that precision conversion. No support
restriction beyond the official q=3 algorithm is added to improve scores.

run_ocean_production.py stores source/input hashes and compact errors locally.
Raw inputs and reference arrays remain on the remote machines. The first
public-network relay timed out and its incomplete file was retained but not
used. Transfer resumed over a temporary LAN endpoint restricted to the 92
client IP and a random per-run capability, serving only the two named files.
It closed after both transfers; complete files were SHA-256 verified. No raw
local data files were written.
No hardware settings, power/clock limits or protections are changed.

## Results — 27 September 2026

Both precisions passed for both saved phase-zero initial states: 24 full-grid
field checks. Reference snapshot on 93: gpu-production-probe-20260927T021020Z;
GPU snapshot on 92: /root/hos-production-check-20260927T021020Z. Compact
results, input-part checksums and unchanged module-object provenance are in
production_rhs_results.json. Local raw compact evidence is under the matching
production timestamp in artifacts/hos_gpu; neither full fields nor initial
campaign files were fetched locally.

| State | Precision | eta full relative L2 | psi full relative L2 | eta nonlinear relative L2 | psi nonlinear relative L2 |
|---|---|---:|---:|---:|---:|
| Low .02 | FP64 | 3.95056e-15 | 3.87031e-16 | 6.98420e-14 | 3.72685e-15 |
| Low .02 | FP32 | 1.94876e-6 | 4.49710e-8 | 2.60553e-5 | 2.12170e-6 |
| High .12 | FP64 | 3.91084e-15 | 4.92101e-16 | 1.11391e-14 | 3.55568e-15 |
| High .12 | FP32 | 1.94884e-6 | 1.80161e-7 | 4.52191e-6 | 2.05085e-6 |

The more sensitive low-amplitude nonlinear eta field still passes the
unchanged 5e-4 FP32 relative gate. This supports pursuing FP32 on the actual
spatial RHS, but does not make a 1e-12 adaptive time tolerance meaningful in
FP32, or certify long-time phase/energy behavior. Only phase zero and t=0
were checked here; the other global phases and evolved states remain untested.
