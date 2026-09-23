#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "$0")" && pwd)
cd "$here/../.."
export OMP_NUM_THREADS=2
build=artifacts/order4-build
out=artifacts/order4-smoke
"$build/wit_eta44" "$out/wit_input.bin" "$out/wit.bin" --sector 44
"$build/wit_eta44_materialized" "$out/wit_input.bin" "$out/wit_materialized.bin" --sector 44
for rank in 6 8 10; do
  python3 paper/order4/run_cpp.py gl "$out/gl_input.bin" "$out/gl${rank}.bin" --rank "$rank" --threads 2
done
python3 paper/order4/run_cpp.py gl "$out/gl_input.bin" "$out/gl6_split.bin" --rank 6 --split-workers 2
