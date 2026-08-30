#!/usr/bin/env bash
set -euo pipefail

if (($# == 0)); then
  echo "usage: tools/run_fortran_targets.sh TARGET..." >&2
  exit 2
fi

build_sh="${BUILD_SH:-./build.sh}"
profile="${FPM_PROFILE:-debug}"
action="${FPM_ACTION:-test}"
target_kind="${FPM_TARGET_KIND:-target}"
timing_csv="${BEACH_FORTRAN_TIMING_CSV:-}"

case "$action" in
  test|run) ;;
  *)
    echo "FPM_ACTION must be test or run" >&2
    exit 2
    ;;
esac
case "$target_kind" in
  target)
    target_options=(--target)
    ;;
  example)
    [[ "$action" == run ]] || { echo "example targets require FPM_ACTION=run" >&2; exit 2; }
    target_options=(--example --target)
    ;;
  *)
    echo "FPM_TARGET_KIND must be target or example" >&2
    exit 2
    ;;
esac

if [[ -n "$timing_csv" ]]; then
  mkdir -p "$(dirname "$timing_csv")"
  if [[ ! -s "$timing_csv" ]]; then
    printf 'target,profile,status,elapsed_seconds\n' >"$timing_csv"
  fi
fi

record_timing() {
  local target="$1"
  local status="$2"
  local elapsed="$3"
  [[ -z "$timing_csv" ]] || printf '%s,%s,%s,%s\n' "$target" "$profile" "$status" "$elapsed" >>"$timing_csv"
}

elapsed_seconds() {
  awk -v start="$1" -v end="$2" 'BEGIN { printf "%.3f", (end - start) / 1000000000.0 }'
}

for target in "$@"; do
  echo "==> fpm $action ${target_options[*]} $target (profile=$profile)"
  start_ns="$(date +%s%N)"
  if FPM_ACTION="$action" FPM_PROFILE="$profile" \
    "$build_sh" "${target_options[@]}" "$target"; then
    end_ns="$(date +%s%N)"
    record_timing "$target" passed "$(elapsed_seconds "$start_ns" "$end_ns")"
  else
    status=$?
    end_ns="$(date +%s%N)"
    record_timing "$target" failed "$(elapsed_seconds "$start_ns" "$end_ns")"
    exit "$status"
  fi
done
