#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 6 ]]; then
  echo "usage: $0 BINARY MERGER INPUT OUTPUT_DIR RANK MAX_WORKERS [FIRST_CPU]" >&2
  exit 2
fi
binary=$1
merger=$2
input=$3
output_dir=$4
rank=$5
max_workers=$6
first_cpu=${7:-0}
last_cpu=$((first_cpu+max_workers-1))
mkdir -p "$output_dir"
start_ns=$(date +%s%N)
partials=()
running=0
launch_task() {
  local output=$1
  shift
  partials+=("$output.bin")
  taskset -c "$first_cpu-$last_cpu" "$binary" "$input" "$output.bin" \
    "$rank" 1 "$@" >"$output.log" 2>&1 &
  running=$((running+1))
  if ((running>=max_workers)); then
    wait -n
    running=$((running-1))
  fi
}
for ((outer4=0; outer4<rank; ++outer4)); do
  launch_task "$output_dir/base_${outer4}" "$outer4" -1
done
for ((outer4=0; outer4<rank; ++outer4)); do
  for ((outer3=0; outer3<rank; ++outer3)); do
    launch_task "$output_dir/k31_${outer4}_${outer3}" "$outer4" "$outer3"
  done
done
while ((running>0)); do
  wait -n
  running=$((running-1))
done
end_ns=$(date +%s%N)
wall_seconds=$(awk -v a="$start_ns" -v b="$end_ns" \
  'BEGIN {printf "%.9f", (b-a)/1e9}')
"$merger" "$output_dir/merged.bin" "${partials[@]}"
printf 'CPP_GL_ETA44_SPLIT_WALL_SECONDS=%s\n' "$wall_seconds" | \
  tee "$output_dir/wall_seconds.txt"
printf 'CPP_GL_ETA44_SPLIT_MAX_WORKERS=%s\n' "$max_workers" | \
  tee "$output_dir/max_workers.txt"
if [[ ${KEEP_PARTIALS:-0} != 1 ]]; then
  rm -- "${partials[@]}"
fi
