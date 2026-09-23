#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
out=${1:-"$root/artifacts/order4-build"}
mkdir -p "$out"
cxx=${CXX:-g++}
flags=(-O3 -std=c++17 -fopenmp)
fftw_flags=()
if [[ -n ${FFTW_PREFIX:-} ]]; then
  fftw_flags=(-I"$FFTW_PREFIX/include" -L"$FFTW_PREFIX/lib")
fi
"$cxx" "${flags[@]}" "$here/cpp/gl_eta44_fftw.cpp" "${fftw_flags[@]}" \
  -lfftw3_threads -lfftw3 -lm -lpthread -o "$out/gl_eta44"
"$cxx" "${flags[@]}" "$here/cpp/directional_eta4_openmp.cpp" -o "$out/wit_eta44_materialized"
patch --silent -o "$out/wit_eta44_streaming.cpp" \
  "$here/cpp/directional_eta4_openmp.cpp" "$here/cpp/directional_eta4_true_streaming.patch"
"$cxx" "${flags[@]}" "$out/wit_eta44_streaming.cpp" -o "$out/wit_eta44"
"$cxx" -O3 -std=c++17 "$here/cpp/merge_eta44_outer.cpp" -o "$out/merge_eta44_outer"
{
  "$cxx" --version
  printf 'flags: -O3 -std=c++17 -fopenmp\n'
  printf 'FFTW_PREFIX: %s\n' "${FFTW_PREFIX:-system}"
  uname -a
  sha256sum "$out/gl_eta44" "$out/wit_eta44" "$out/wit_eta44_materialized"
} > "$out/build-report.txt"
printf 'Built GL and WIT executors in %s\n' "$out"
