#!/bin/bash
#SBATCH -p gr20001a
#SBATCH --rsc p=1:t=112:c=112
#SBATCH -t 120:00:00
#SBATCH -o stdout.%J.log
#SBATCH -e stderr.%J.log

set -eu

MPI_RANKS=4
THREADS_PER_RANK=28
TOTAL_CORES=$((MPI_RANKS * THREADS_PER_RANK))

if [ "${TOTAL_CORES}" -gt 112 ]; then
  echo "MPI_RANKS * THREADS_PER_RANK must be <= 112 (now: ${TOTAL_CORES})" >&2
  exit 1
fi

export OMP_NUM_THREADS=$THREADS_PER_RANK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
CASE_DIR="${ROOT_DIR}/outputs/exp_sphere_details_mpi"
CONFIG_PATH="${CASE_DIR}/beach.toml"
OUTPUT_DIR="${CASE_DIR}/outputs/latest"

source /home/b/b36291/large0/.venv/bin/activate

cd "${ROOT_DIR}"

date

beach-estimate-workload "${CONFIG_PATH}"

date

make run-mpi \
  CONFIG="${CONFIG_PATH}" \
  MPI_FC=mpiifort \
  MPI_NP="${MPI_RANKS}" \
  MPI_RUNNER="srun -n ${MPI_RANKS} -c ${THREADS_PER_RANK} --cpu-bind=cores"

date

beach-inspect "${OUTPUT_DIR}" \
  --save-bar "${OUTPUT_DIR}/charges_bar.png" \
  --save-mesh "${OUTPUT_DIR}/charges_mesh.png" \
  --save-potential-mesh "${OUTPUT_DIR}/potential_mesh.png" \
  --potential-self-term area-equivalent

beach-animate-history "${OUTPUT_DIR}" \
  --quantity charge \
  --save-gif "${OUTPUT_DIR}/charge_history.gif" \
  --total-frames 300

beach-animate-history "${OUTPUT_DIR}" \
  --quantity potential \
  --save-gif "${OUTPUT_DIR}/potential_history.gif" \
  --potential-self-term area-equivalent \
  --total-frames 300

date
