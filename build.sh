#!/usr/bin/env bash
# Run fpm with BEACH build metadata in a way that can keep development builds
# incrementally reusable.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

FPM="${FPM:-fpm}"
FPM_ACTION="${FPM_ACTION:-build}"
FPM_PROFILE="${FPM_PROFILE:-${PROFILE:-}}"
PREFIX="${PREFIX:-${HOME}/.local}"

read_fpm_version() {
  sed -n -E 's/^version[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/p' fpm.toml | head -n 1
}

BASE_VERSION="$(read_fpm_version)"
if [[ -z "${BASE_VERSION}" ]]; then
  echo "ERROR: failed to read version from fpm.toml." >&2
  exit 1
fi

VERSION_MODE="${BEACH_VERSION_MODE:-git}"
VERSION_OVERRIDE="${BEACH_VERSION_OVERRIDE:-}"

if [[ -n "${VERSION_OVERRIDE}" ]]; then
  FULL_VERSION="${VERSION_OVERRIDE}"
  RESOLVED_VERSION_MODE="override"
else
  case "${VERSION_MODE}" in
    git)
      if git rev-parse --git-dir >/dev/null 2>&1; then
        GIT_HASH="$(git describe --tags --always --dirty)"
      else
        GIT_HASH="nogit"
      fi
      FULL_VERSION="${BASE_VERSION}-${GIT_HASH}"
      RESOLVED_VERSION_MODE="git"
      ;;
    dev|static)
      # Keep fpm's compile-flag hash stable across ordinary development commits.
      FULL_VERSION="${BASE_VERSION}-dev"
      RESOLVED_VERSION_MODE="dev"
      ;;
    plain)
      FULL_VERSION="${BASE_VERSION}"
      RESOLVED_VERSION_MODE="plain"
      ;;
    *)
      echo "ERROR: BEACH_VERSION_MODE must be one of: git, dev, plain" >&2
      echo "       Use BEACH_VERSION_OVERRIDE to pass an exact version string." >&2
      exit 1
      ;;
  esac
fi

VERSION_FLAGS="-D__BEACH_VERSION__=\\'${FULL_VERSION}\\' -D__BEACH_VERSION_MODE__=\\'${RESOLVED_VERSION_MODE}\\'"
BASE_FFLAGS="${FPM_FFLAGS:-${FFLAGS:-}}"
if [[ -n "${BASE_FFLAGS}" ]]; then
  EFFECTIVE_FFLAGS="${BASE_FFLAGS} ${VERSION_FLAGS}"
else
  EFFECTIVE_FFLAGS="${VERSION_FLAGS}"
fi

echo "[build.sh] FPM_ACTION=${FPM_ACTION}"
echo "[build.sh] BEACH_VERSION=${FULL_VERSION}"
echo "[build.sh] BEACH_VERSION_MODE=${RESOLVED_VERSION_MODE}"
echo "[build.sh] FPM_FFLAGS=${EFFECTIVE_FFLAGS}"

case "${FPM_ACTION}" in
  build)
    : "${FPM_PROFILE:=release}"
    exec "${FPM}" build --profile "${FPM_PROFILE}" --flag "${EFFECTIVE_FFLAGS}" "$@"
    ;;
  run)
    : "${FPM_PROFILE:=release}"
    exec "${FPM}" run --profile "${FPM_PROFILE}" --flag "${EFFECTIVE_FFLAGS}" "$@"
    ;;
  test)
    : "${FPM_PROFILE:=debug}"
    exec "${FPM}" test --profile "${FPM_PROFILE}" --flag "${EFFECTIVE_FFLAGS}" "$@"
    ;;
  install)
    : "${FPM_PROFILE:=release}"
    exec "${FPM}" install --profile "${FPM_PROFILE}" --flag "${EFFECTIVE_FFLAGS}" --prefix "${PREFIX}" "$@"
    ;;
  *)
    echo "ERROR: FPM_ACTION must be one of: build, run, test, install" >&2
    exit 1
    ;;
esac
