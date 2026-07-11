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
started_epoch="$(date +%s)"
started_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
commit="$(git rev-parse HEAD)"
dirty=false
[[ -n "$(git status --porcelain)" ]] && dirty=true
fc_command="${FC:-gfortran}"
mpi_fc_command="${MPI_FC:-mpiifort}"
mpi_runner="${BEACH_RELEASE_MPI_RUNNER:-}"
if [[ -z "$mpi_runner" && -n "${SLURM_JOB_ID:-}" ]]; then
  mpi_runner="srun"
fi
fc_version="$($fc_command --version 2>/dev/null | head -n 1 || printf unavailable)"
mpi_fc_version="$($mpi_fc_command --version 2>/dev/null | head -n 1 || printf unavailable)"

cat >"$manifest" <<EOF
schema_version=1
status=$(if ((dry_run)); then printf planned; else printf running; fi)
started_utc=$started_utc
hostname=$host
slurm_job_id=${SLURM_JOB_ID:-none}
git_commit=$commit
git_dirty=$dirty
fortran_compiler=$fc_version
mpi_fortran_compiler=$mpi_fc_version
test_l3.command=make test-l3
far_correction.command=make test-fortran-far-correction
mpi.command=make test-mpi$(if [[ -n "$mpi_runner" ]]; then printf ' MPI_RUNNER=%s' "$mpi_runner"; fi)
EOF

run_gate() {
  local name="$1"
  shift
  local gate_start gate_end
  gate_start="$(date +%s)"
  if ((dry_run)); then
    printf '%s.status=planned\n%s.elapsed_seconds=0\n' "$name" "$name" >>"$manifest"
    return
  fi
  if "$@"; then
    gate_end="$(date +%s)"
    printf '%s.status=passed\n%s.elapsed_seconds=%s\n' \
      "$name" "$name" "$((gate_end - gate_start))" >>"$manifest"
  else
    gate_end="$(date +%s)"
    printf '%s.status=failed\n%s.elapsed_seconds=%s\nstatus=failed\n' \
      "$name" "$name" "$((gate_end - gate_start))" >>"$manifest"
    exit 1
  fi
}

run_gate test_l3 make test-l3
run_gate far_correction make test-fortran-far-correction
if [[ -n "$mpi_runner" ]]; then
  run_gate mpi make test-mpi "MPI_RUNNER=$mpi_runner"
else
  run_gate mpi make test-mpi
fi

if ((!dry_run)); then
  finished_epoch="$(date +%s)"
  printf 'status=passed\nelapsed_seconds=%s\n' "$((finished_epoch - started_epoch))" >>"$manifest"
fi
printf 'physics release manifest: %s\n' "$manifest"
