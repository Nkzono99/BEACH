#!/bin/bash
# Generic camphor MPI+OpenMP hybrid batch job template for BEACH.
#SBATCH -p gr20001a
#SBATCH --rsc p=5:t=112:c=112
#SBATCH -t 120:00:00
#SBATCH -o stdout.%J.log
#SBATCH -e stderr.%J.log

set -eu

MPI_RANKS="${MPI_RANKS:-${SLURM_NTASKS:-112}}"
THREADS_PER_RANK="${THREADS_PER_RANK:-${SLURM_DPC_CPUS:-1}}"

if [ "${MPI_RANKS}" -le 0 ] || [ "${THREADS_PER_RANK}" -le 0 ]; then
  echo "MPI_RANKS and THREADS_PER_RANK must be > 0." >&2
  exit 1
fi

export OMP_NUM_THREADS="${THREADS_PER_RANK}"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export BEACH_PROFILE="${BEACH_PROFILE:-1}"
export BEACH_PROFILE_DETAIL="${BEACH_PROFILE_DETAIL:-0}"

CONFIG_PATH="${CONFIG_PATH:-beach.toml}"
OUTPUT_DIR="${OUTPUT_DIR:-outputs/latest}"

date

beachx workload "${CONFIG_PATH}" --mpi-ranks "${MPI_RANKS}"

date

srun -n "${MPI_RANKS}" -c "${THREADS_PER_RANK}" --cpu-bind=cores beach "${CONFIG_PATH}"

date

beachx inspect "${OUTPUT_DIR}" \
  --save-bar "${OUTPUT_DIR}/charges_bar.png" \
  --save-mesh "${OUTPUT_DIR}/charges_mesh.png" \
  --save-potential-mesh "${OUTPUT_DIR}/potential_mesh.png" \
  --potential-self-term area-equivalent

beachx animate "${OUTPUT_DIR}" \
  --quantity charge \
  --save-gif "${OUTPUT_DIR}/charge_history.gif" \
  --total-frames 300

beachx animate "${OUTPUT_DIR}" \
  --quantity potential \
  --save-gif "${OUTPUT_DIR}/potential_history.gif" \
  --potential-self-term area-equivalent \
  --total-frames 300

date
