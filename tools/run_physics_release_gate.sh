#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: tools/run_physics_release_gate.sh [--dry-run] [--manifest PATH]

Run the L3, far-correction, and MPI physics gates sequentially and record a
reproducibility manifest. KUDPC login nodes are rejected outside Slurm.
EOF
}

dry_run=0
manifest="build/physics-release/manifest.txt"
while (($# > 0)); do
  case "$1" in
    --dry-run)
      dry_run=1
      shift
      ;;
    --manifest)
      [[ $# -ge 2 ]] || { echo "--manifest requires a path" >&2; exit 2; }
      manifest="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

host="${BEACH_RELEASE_GATE_HOSTNAME:-$(hostname)}"
case "$host" in
  camphor*|laurel*|cinnamon*|gardenia*)
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
      echo "physics release gate must run on a compute node inside a Slurm allocation" >&2
      exit 2
    fi
    ;;
esac

mkdir -p "$(dirname "$manifest")"
artifact_dir="$(dirname "$manifest")"
convergence_csv="$artifact_dir/convergence.csv"
test_l3_timings_csv="$artifact_dir/test_l3-target-timings.csv"
far_correction_timings_csv="$artifact_dir/far_correction-target-timings.csv"
max_rss_kb="${BEACH_RELEASE_MAX_RSS_KB:-8388608}"
time_command="${BEACH_RELEASE_TIME_COMMAND:-/usr/bin/time}"
if [[ ! "$max_rss_kb" =~ ^[1-9][0-9]*$ ]]; then
  echo "BEACH_RELEASE_MAX_RSS_KB must be a positive integer" >&2
  exit 2
fi
started_epoch="$(date +%s)"
started_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
commit="$(git rev-parse HEAD)"
dirty=false
[[ -n "$(git status --porcelain)" ]] && dirty=true
fc_command="${FPM_FC:-${FC:-gfortran}}"
export FPM_FC="$fc_command"
mpi_fc_command="${MPI_FC:-mpiifx}"
if ! command -v "$mpi_fc_command" >/dev/null 2>&1; then
  echo "MPI compiler not found: $mpi_fc_command (set MPI_FC to an available wrapper)" >&2
  exit 2
fi
export MPI_FC="$mpi_fc_command"
mpi_runner="${BEACH_RELEASE_MPI_RUNNER:-}"
if [[ -z "$mpi_runner" && -n "${SLURM_JOB_ID:-}" ]]; then
  mpi_runner="srun"
fi
fc_version="$($fc_command --version 2>/dev/null | head -n 1 || printf unavailable)"
mpi_fc_version="$($mpi_fc_command --version 2>/dev/null | head -n 1 || printf unavailable)"

cat >"$manifest" <<EOF
schema_version=2
status=$(if ((dry_run)); then printf planned; else printf running; fi)
started_utc=$started_utc
hostname=$host
slurm_job_id=${SLURM_JOB_ID:-none}
git_commit=$commit
git_dirty=$dirty
fortran_compiler=$fc_version
mpi_fortran_compiler=$mpi_fc_version
test_l3.command=make test-fortran-release-correctness
far_correction.command=make test-fortran-far-correction
mpi.command=make test-mpi$(if [[ -n "$mpi_runner" ]]; then printf ' MPI_RUNNER=%s' "$mpi_runner"; fi)
mpi_cache.command=make test-mpi-periodic-cache$(if [[ -n "$mpi_runner" ]]; then printf ' MPI_RUNNER=%s' "$mpi_runner"; fi)
budget.max_rss_kb=$max_rss_kb
cache.warm_max_cold_ratio=0.5
convergence_csv=$convergence_csv
test_l3.target_timings_csv=$test_l3_timings_csv
far_correction.target_timings_csv=$far_correction_timings_csv
convergence.required_categories=boris_dt,panel_fmm_order,rough_panel_mesh,rough_outer_grid,rough_accessibility,outer_orbit_dt
EOF

run_gate() {
  local name="$1"
  shift
  local gate_start gate_end gate_rss log_path time_path
  gate_start="$(date +%s)"
  log_path="$artifact_dir/$name.log"
  time_path="$artifact_dir/$name.time"
  if ((dry_run)); then
    printf '%s.status=planned\n%s.elapsed_seconds=0\n%s.max_rss_kb=0\n%s.log=%s\n' \
      "$name" "$name" "$name" "$name" "$log_path" >>"$manifest"
    return
  fi
  if [[ ! -x "$time_command" ]]; then
    echo "GNU time command is unavailable: $time_command" >&2
    exit 2
  fi
  if "$time_command" -v -o "$time_path" "$@" 2>&1 | tee "$log_path"; then
    gate_end="$(date +%s)"
    gate_rss="$(awk -F: '/Maximum resident set size/{gsub(/[[:space:]]/, "", $2); print $2}' "$time_path")"
    [[ "$gate_rss" =~ ^[0-9]+$ ]] || gate_rss=0
    if ((gate_rss > max_rss_kb)); then
      printf '%s.status=failed_memory_budget\n%s.elapsed_seconds=%s\n%s.max_rss_kb=%s\nstatus=failed\n' \
        "$name" "$name" "$((gate_end - gate_start))" "$name" "$gate_rss" >>"$manifest"
      echo "$name exceeded memory budget: ${gate_rss} KiB > ${max_rss_kb} KiB" >&2
      exit 1
    fi
    printf '%s.status=passed\n%s.elapsed_seconds=%s\n%s.max_rss_kb=%s\n%s.log=%s\n' \
      "$name" "$name" "$((gate_end - gate_start))" "$name" "$gate_rss" "$name" "$log_path" >>"$manifest"
  else
    gate_end="$(date +%s)"
    gate_rss="$(awk -F: '/Maximum resident set size/{gsub(/[[:space:]]/, "", $2); print $2}' "$time_path" 2>/dev/null || true)"
    [[ "$gate_rss" =~ ^[0-9]+$ ]] || gate_rss=0
    printf '%s.status=failed\n%s.elapsed_seconds=%s\n%s.max_rss_kb=%s\n%s.log=%s\nstatus=failed\n' \
      "$name" "$name" "$((gate_end - gate_start))" "$name" "$gate_rss" "$name" "$log_path" >>"$manifest"
    exit 1
  fi
}

if ((!dry_run)); then
  rm -f "$convergence_csv" "$test_l3_timings_csv" "$far_correction_timings_csv"
fi
run_gate test_l3 env BEACH_FORTRAN_TIMING_CSV="$test_l3_timings_csv" make test-fortran-release-correctness
run_gate far_correction env BEACH_FORTRAN_TIMING_CSV="$far_correction_timings_csv" make test-fortran-far-correction
if [[ -n "$mpi_runner" ]]; then
  run_gate mpi make test-mpi "MPI_RUNNER=$mpi_runner"
  run_gate mpi_cache make test-mpi-periodic-cache "MPI_RUNNER=$mpi_runner"
else
  run_gate mpi make test-mpi
  run_gate mpi_cache make test-mpi-periodic-cache
fi

if ((!dry_run)); then
  printf 'category,configuration,metric_1,metric_2,metric_3,acceptance\n' >"$convergence_csv"
  sed -n 's/^.*BEACH_CONVERGENCE,//p' \
    "$artifact_dir/test_l3.log" "$artifact_dir/far_correction.log" \
    "$artifact_dir/mpi.log" "$artifact_dir/mpi_cache.log" >>"$convergence_csv"
  for category in boris_dt panel_fmm_order rough_panel_mesh rough_outer_grid rough_accessibility outer_orbit_dt; do
    if ! grep -q "^${category}," "$convergence_csv"; then
      printf 'convergence.status=failed_missing_%s\nstatus=failed\n' "$category" >>"$manifest"
      echo "missing convergence category: $category" >&2
      exit 1
    fi
  done
  printf 'convergence.status=passed\nconvergence.rows=%s\n' \
    "$(( $(wc -l <"$convergence_csv") - 1 ))" >>"$manifest"
  finished_epoch="$(date +%s)"
  printf 'status=passed\nelapsed_seconds=%s\n' "$((finished_epoch - started_epoch))" >>"$manifest"
fi
printf 'physics release manifest: %s\n' "$manifest"
