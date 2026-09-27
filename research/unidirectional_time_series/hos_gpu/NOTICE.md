# HOS-Ocean-derived research code

The `ocean_*.cuh`, `ocean_*.cu`, and `ocean_*.f90` adaptation/diagnostic files,
and `mpi_global_zero_gravity.patch`, follow or instrument HOS-Ocean v2.1.0,
upstream commit `4deb3b4913d993c4e6ea16f736e5fc5792e14f12`.

Original HOS-Ocean copyright:
Copyright (C) 2014 — LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598.

These HOS-Ocean-derived files are distributed under GPL-3.0-or-later, as
identified by their SPDX headers and this notice. The upstream license text
is retained in LICENSE-HOS-OCEAN.txt. The repository's root MIT license does
not replace the GPL terms for these files. The independent original GL
library and its existing license are unchanged.

Changes in this research directory include a CUDA M=5/q=3 adapter, the
linear-split Cash--Karp controller, diagnostics, and an isolated global-zero
MPI correction. The original production source was not overwritten. See
OCEAN_ALIGNMENT.md and LONG_RUN_BENCHMARK.md for provenance and scope.
