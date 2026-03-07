#!/usr/bin/env bash
# Install BEACH via fpm with optional host-specific optimization profile.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${ROOT_DIR}/env"
cd "${ROOT_DIR}"

detect_host_profile() {
  local host_name
  host_name="$(uname -n 2>/dev/null || true)"
  case "${host_name}" in
    camphor*) echo "camphor" ;;
    fn*) echo "fugaku" ;;
    *) echo "generic" ;;
  esac
}

load_env_file() {
  local env_file="$1"
  if [[ ! -f "${env_file}" ]]; then
    return 1
  fi
  # shellcheck disable=SC1090
  source "${env_file}"
  return 0
}

load_modules_if_requested() {
  local module_name
  if [[ "${USE_MODULES}" != "1" ]]; then
    return 0
  fi

  if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if ! type module >/dev/null 2>&1; then
    echo "[install.sh] warning: module command not found; continuing with current environment." >&2
    return 0
  fi

  for module_name in ${MODULES_TO_UNLOAD:-}; do
    module unload "${module_name}" >/dev/null 2>&1 || true
  done
  for module_name in ${MODULES_TO_LOAD:-}; do
    module load "${module_name}"
  done
}

setup_generic_defaults() {
  : "${USE_MODULES:=0}"
  : "${FC:=gfortran}"
  : "${FFLAGS:=-O3 -fopenmp}"
  : "${LDFLAGS:=}"
}

: "${BUILD_PROFILE:=auto}"
REQUESTED_PROFILE="${BUILD_PROFILE}"
if [[ "${BUILD_PROFILE}" == "auto" ]]; then
  BUILD_PROFILE="$(detect_host_profile)"
fi

if [[ -f "${ENV_DIR}/common.env" ]]; then
  load_env_file "${ENV_DIR}/common.env"
fi

if [[ "${BUILD_PROFILE}" != "generic" ]]; then
  if ! load_env_file "${ENV_DIR}/${BUILD_PROFILE}.env"; then
    if [[ "${REQUESTED_PROFILE}" == "auto" ]]; then
      echo "[install.sh] warning: profile env not found (${BUILD_PROFILE}). Fallback to generic." >&2
      BUILD_PROFILE="generic"
    else
      echo "[install.sh] error: profile env not found (${ENV_DIR}/${BUILD_PROFILE}.env)." >&2
      exit 1
    fi
  fi
fi

setup_generic_defaults

: "${PROFILE:=release}"
: "${PREFIX:=${HOME}/.local}"

load_modules_if_requested

export FPM_FC="${FC}"
export FPM_FFLAGS="${FFLAGS}"
if [[ -n "${LDFLAGS}" ]]; then
  export FPM_LDFLAGS="${LDFLAGS}"
fi

echo "[install.sh] BUILD_PROFILE=${BUILD_PROFILE}"
echo "[install.sh] FPM_FC=${FPM_FC}"
echo "[install.sh] FPM_FFLAGS=${FPM_FFLAGS}"
echo "[install.sh] PREFIX=${PREFIX}"

if ! command -v "${FPM_FC}" >/dev/null 2>&1; then
  if [[ "${REQUESTED_PROFILE}" == "auto" && "${BUILD_PROFILE}" != "generic" ]]; then
    echo "[install.sh] warning: compiler not found for ${BUILD_PROFILE} profile (${FPM_FC})." >&2
    echo "[install.sh] hint: set BUILD_PROFILE=generic or prepare module/compilers." >&2
  fi
  echo "[install.sh] error: compiler not found: ${FPM_FC}" >&2
  exit 1
fi

fpm install --profile "${PROFILE}" --compiler "${FPM_FC}" --flag "${FPM_FFLAGS}" --prefix "${PREFIX}"

echo "[install.sh] done: ${PREFIX}/bin/beach"
