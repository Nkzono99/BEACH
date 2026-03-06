#!/usr/bin/env bash
# Install BEACH with Intel MPI Fortran toolchain via fpm.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

: "${USE_MODULES:=1}"
: "${PRGENV_MOD:=PrgEnvIntel}"
: "${INTEL_COMPILER_MOD:=intel/2023.2}"
: "${INTEL_MPI_MOD:=intelmpi/2023.2}"

if [[ "${USE_MODULES}" == "1" ]]; then
  if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi

  if type module >/dev/null 2>&1; then
    module unload "${PRGENV_MOD}" >/dev/null 2>&1 || true
    module unload intel-python >/dev/null 2>&1 || true
    module load "${PRGENV_MOD}"
    module load "${INTEL_COMPILER_MOD}" "${INTEL_MPI_MOD}"
  else
    echo "[install.sh] warning: module command not found; continuing with current environment." >&2
  fi
fi

: "${CC:=mpiicc}"
: "${CXX:=mpiicpc}"
: "${FC:=mpiifort}"
: "${PROFILE:=release}"
: "${PREFIX:=${HOME}/.local}"

DEFAULT_OPT="-ipo -O3 -no-prec-div -xHost -fp-model precise"
DEFAULT_OPT+=" -inline-level=2 -inline-forceinline -inline-factor=300 -no-inline-max-total-size"
DEFAULT_OPT+=" -qopt-report=5 -qopenmp -traceback -fpp -mcmodel=medium"

: "${FFLAGS:=${DEFAULT_OPT}}"
: "${CFLAGS:=${DEFAULT_OPT}}"
: "${CXXFLAGS:=${DEFAULT_OPT}}"
: "${LDFLAGS:=}"

export FPM_CC="${CC}"
export FPM_CXX="${CXX}"
export FPM_FC="${FC}"
export FPM_CFLAGS="${CFLAGS}"
export FPM_CXXFLAGS="${CXXFLAGS}"
export FPM_FFLAGS="${FFLAGS}"
if [[ -n "${LDFLAGS}" ]]; then
  export FPM_LDFLAGS="${LDFLAGS}"
fi

echo "[install.sh] FPM_FC=${FPM_FC}"
echo "[install.sh] FPM_FFLAGS=${FPM_FFLAGS}"
echo "[install.sh] PREFIX=${PREFIX}"

if ! command -v "${FPM_FC}" >/dev/null 2>&1; then
  echo "[install.sh] error: compiler not found: ${FPM_FC}" >&2
  exit 1
fi

fpm install --profile "${PROFILE}" --compiler "${FPM_FC}" --flag "${FPM_FFLAGS}" --prefix "${PREFIX}"

echo "[install.sh] done: ${PREFIX}/bin/beach"
