#!/bin/sh
# Adapted from the archived, audited remote_completion_email_watcher.sh.
# The historical source is unchanged. Recipient stays in the remote run.sh.
set -eu
run_dir=${1:?remote run directory is required}
run_id=${2:?run id is required}
mail_status="$run_dir/mail-status.txt"
watcher_status="$run_dir/mail-watcher-status.txt"
if test -f "$mail_status" && grep -q '^accepted_by_local_mail$' "$mail_status"; then exit 0; fi
if ! mkdir "$run_dir/.mail-watcher-lock" 2>/dev/null; then exit 0; fi
printf '%s\n' "$$" > "$run_dir/mail-watcher.pid"
printf '%s\n' waiting > "$watcher_status.tmp"
mv "$watcher_status.tmp" "$watcher_status"
while ! test -f "$run_dir/finished.utc"; do sleep 60; done
sleep 3
state=$(cat "$run_dir/status.txt" 2>/dev/null || printf unknown)
exit_code=$(cat "$run_dir/exit-code.txt" 2>/dev/null || printf unknown)
recipient=$(sed -n 's/^EMAIL="${EMAIL:-\([^}]*\)}"$/\1/p' /usr/local/bin/run.sh | head -n 1)
if test -z "$recipient"; then
  printf '%s\n' failed_recipient_not_found > "$watcher_status"
  exit 1
fi
printf '%s\n' sending > "$watcher_status"
if {
  printf 'RunId: %s\nStatus: %s\nExit code: %s\n' "$run_id" "$state" "$exit_code"
  printf 'Remote run: %s\n' "$run_dir"
  printf '\nSixteen native-wall OW3D runs: one focused group and three random realizations, four phases each.\n'
  printf 'Randomization preserves modal amplitudes. No local raw-data download.\n'
  printf 'Completion means native final time and raw-output integrity checks; physical/model accuracy is not certified.\n'
  printf '\nSummary: %s/summary.json\nPer-case processed eta/phi: cases/*/processed/surface_strip.mat\n' "$run_dir"
  if test -f "$run_dir/runtime.log"; then printf '\nLast log lines:\n'; tail -40 "$run_dir/runtime.log"; fi
} | /usr/bin/mail -s "[Green-Laplace OW3D] $state $run_id" "$recipient"; then
  printf '%s\n' accepted_by_local_mail > "$mail_status.tmp"
  mv "$mail_status.tmp" "$mail_status"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$run_dir/mail-sent.utc.tmp"
  mv "$run_dir/mail-sent.utc.tmp" "$run_dir/mail-sent.utc"
  printf '%s\n' completed > "$watcher_status.tmp"
  mail_rc=0
else
  mail_rc=$?
  printf 'mail_failed_exit_%s\n' "$mail_rc" > "$mail_status.tmp"
  mv "$mail_status.tmp" "$mail_status"
  printf '%s\n' failed > "$watcher_status.tmp"
fi
mv "$watcher_status.tmp" "$watcher_status"
exit "$mail_rc"
