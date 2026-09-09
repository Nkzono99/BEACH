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

prepare_fpm_cache() {
  local compiler="${FPM_FC:-gfortran}"
  local compiler_name compiler_path compiler_version context context_file cache_dir cache_suffix
  local source_dir
  local source_dirs=()

  case "${FPM_ACTION}" in
    build|run|test|install) ;;
    *) return ;;
  esac
  while (($#)); do
    case "$1" in
      --) break ;;
      --help|-h|--version|-V|--list|--list=*|--show-model) return ;;
      --compiler)
        (($# >= 2)) || return 0
        compiler="$2"
        shift
        ;;
      --compiler=*) compiler="${1#--compiler=}" ;;
    esac
    shift
  done

  # Let fpm report an unavailable compiler without discarding a usable cache.
  compiler_path="$(command -v "${compiler}")" || return 0
  compiler_path="$(readlink -f -- "${compiler_path}")" || return 0
  compiler_version="$(LC_ALL=C "${compiler}" --version 2>&1)" || true
  compiler_name="${compiler##*/}"
  context_file="build/.beach-cache-${compiler_name}.context"
  for source_dir in src app tests/fortran benchmarks/fortran; do
    [[ ! -d "${source_dir}" ]] || source_dirs+=("${source_dir}")
  done
  context="$(
    printf 'compiler_name=%s\ncompiler_path=%s\ncompiler_version:\n%s\nsources:\n' \
      "${compiler_name}" "${compiler_path}" "${compiler_version}"
    if ((${#source_dirs[@]})); then
      find "${source_dirs[@]}" -type f \
        \( -iname '*.f90' -o -iname '*.f' -o -iname '*.for' -o -iname '*.fpp' \
        -o -iname '*.c' -o -iname '*.cc' -o -iname '*.cpp' -o -iname '*.cxx' \) \
        -print | LC_ALL=C sort
    fi
  )"

  if [[ "${BEACH_REBUILD:-0}" == "1" ]] || \
    ! cmp -s "${context_file}" <(printf '%s\n' "${context}"); then
    # fpm keeps library and executable objects in separate flag-hash directories.
    # Invalidate both, including their modules and archives, after a source move.
    # Keep other compilers, dependency checkouts, logs and generated products.
    echo "[build.sh] Resetting ${compiler_name} caches: source layout/compiler changed or rebuild requested."
    for cache_dir in build/"${compiler_name}_"*; do
      [[ -d "${cache_dir}" && ! -L "${cache_dir}" ]] || continue
      cache_suffix="${cache_dir#"build/${compiler_name}_"}"
      [[ "${cache_suffix}" =~ ^[[:xdigit:]]{16}$ ]] || continue
      rm -rf -- "${cache_dir}"
    done
    mkdir -p build
    printf '%s\n' "${context}" >"${context_file}"
  fi
}

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

SOURCE_COMMIT="unknown"
BUILD_ID="${RESOLVED_VERSION_MODE}"
if [[ "${RESOLVED_VERSION_MODE}" == "git" || "${RESOLVED_VERSION_MODE}" == "override" ]]; then
  if git rev-parse --verify HEAD >/dev/null 2>&1; then
    SOURCE_COMMIT="$(git rev-parse --verify HEAD)"
    if [[ -n "$(git status --porcelain=v1 --untracked-files=all)" ]]; then
      SOURCE_STATE="dirty"
    else
      SOURCE_STATE="clean"
    fi
    BUILD_ID="${SOURCE_COMMIT}:${SOURCE_STATE}"
  else
    BUILD_ID="nogit"
  fi
fi

VERSION_FLAGS="-D__BEACH_VERSION__=\\'${FULL_VERSION}\\' -D__BEACH_VERSION_MODE__=\\'${RESOLVED_VERSION_MODE}\\'"
VERSION_FLAGS+=" -D__BEACH_SOURCE_COMMIT__=\\'${SOURCE_COMMIT}\\' -D__BEACH_BUILD_ID__=\\'${BUILD_ID}\\'"
BASE_FFLAGS="${FPM_FFLAGS:-${FFLAGS:-}}"
if [[ -n "${BASE_FFLAGS}" ]]; then
  EFFECTIVE_FFLAGS="${BASE_FFLAGS} ${VERSION_FLAGS}"
else
  EFFECTIVE_FFLAGS="${VERSION_FLAGS}"
fi

echo "[build.sh] FPM_ACTION=${FPM_ACTION}"
echo "[build.sh] BEACH_VERSION=${FULL_VERSION}"
echo "[build.sh] BEACH_VERSION_MODE=${RESOLVED_VERSION_MODE}"
echo "[build.sh] BEACH_SOURCE_COMMIT=${SOURCE_COMMIT}"
echo "[build.sh] BEACH_BUILD_ID=${BUILD_ID}"
echo "[build.sh] FPM_FFLAGS=${EFFECTIVE_FFLAGS}"

prepare_fpm_cache "$@"

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
