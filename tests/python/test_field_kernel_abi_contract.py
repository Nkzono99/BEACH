from __future__ import annotations

import ctypes
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

import beach.fortran_results.kernel as kernel_module


ROOT = Path(__file__).resolve().parents[2]
FORTRAN_SOURCE = ROOT / "src/physics/field_solver/bem_field_kernel_c.f90"
PUBLIC_HEADER = ROOT / "beach/include/beach_field_kernel.h"

C_VOID_P = ctypes.c_void_p
C_INT = ctypes.c_int
C_DOUBLE = ctypes.c_double
EVAL_SIGNATURE = (C_VOID_P, C_INT, C_VOID_P, C_VOID_P)
CTYPES_SIGNATURES = {
    "beach_kernel_get_abi_version": (C_VOID_P, C_VOID_P),
    "beach_kernel_get_build_info": (C_VOID_P, C_INT, C_VOID_P),
    "beach_kernel_create": (ctypes.POINTER(C_VOID_P),),
    "beach_kernel_destroy": (C_VOID_P,),
    "beach_kernel_set_periodic_cache": (C_VOID_P, C_VOID_P, C_INT, C_DOUBLE),
    "beach_kernel_get_periodic_cache_info": (
        C_VOID_P,
        C_VOID_P,
        C_VOID_P,
        C_VOID_P,
        C_INT,
        C_VOID_P,
        C_VOID_P,
        C_INT,
        C_VOID_P,
    ),
    "beach_kernel_build": (
        C_VOID_P,
        C_INT,
        C_VOID_P,
        C_VOID_P,
        C_VOID_P,
        C_DOUBLE,
        C_INT,
        C_INT,
        C_INT,
        C_VOID_P,
        C_VOID_P,
        C_INT,
        C_INT,
        C_DOUBLE,
        C_INT,
        C_VOID_P,
        C_VOID_P,
    ),
    "beach_kernel_update_charges": (C_VOID_P, C_INT, C_VOID_P),
    "beach_kernel_eval_e": EVAL_SIGNATURE,
    "beach_kernel_eval_phi": EVAL_SIGNATURE,
    "beach_kernel_eval_e_direct": EVAL_SIGNATURE,
    "beach_kernel_eval_phi_direct": EVAL_SIGNATURE,
    "beach_kernel_force_on_charges": (
        C_VOID_P,
        C_INT,
        C_VOID_P,
        C_VOID_P,
        C_VOID_P,
        C_VOID_P,
        C_VOID_P,
    ),
}
PUBLIC_SYMBOLS = frozenset(CTYPES_SIGNATURES)


class _FakeFunction:
    argtypes = None
    restype = None

    def __call__(self, *_args: object) -> int:
        return 0


def test_header_and_fortran_expose_the_frozen_public_symbol_set() -> None:
    fortran = FORTRAN_SOURCE.read_text(encoding="utf-8")
    header = PUBLIC_HEADER.read_text(encoding="utf-8")

    bound_symbols = set(
        re.findall(
            r"bind\s*\(\s*C\s*,\s*name\s*=\s*['\"]([^'\"]+)['\"]",
            fortran,
            flags=re.IGNORECASE,
        )
    )
    declared_symbols = set(re.findall(r"\bint\s+(beach_kernel_\w+)\s*\(", header))

    assert declared_symbols == PUBLIC_SYMBOLS
    assert bound_symbols == PUBLIC_SYMBOLS


def test_abi_constants_are_synchronized_across_layers() -> None:
    fortran = FORTRAN_SOURCE.read_text(encoding="utf-8")
    header = PUBLIC_HEADER.read_text(encoding="utf-8")

    assert (
        kernel_module.FIELD_KERNEL_ABI_MAJOR,
        kernel_module.FIELD_KERNEL_ABI_MINOR,
    ) == (
        2,
        1,
    )
    for name, value in {"major": 2, "minor": 1}.items():
        assert f"BEACH_FIELD_KERNEL_ABI_{name.upper()} {value}" in header
        assert f"beach_kernel_abi_{name} = {value}_c_int" in fortran

    status_codes = {"ok": 0, "invalid_handle": 1, "invalid_argument": 2, "not_ready": 3}
    assert set(kernel_module._STATUS_MESSAGES) == set(status_codes.values())
    for name, value in status_codes.items():
        assert f"BEACH_KERNEL_{name.upper()} = {value}" in header
        assert f"beach_kernel_{name} = {value}_c_int" in fortran

    assert kernel_module._FAR_CORRECTION_CODES == {
        "auto": 0,
        "none": 1,
        "cached_kneq0": 3,
    }


def test_public_header_compiles_against_the_frozen_c_abi(tmp_path: Path) -> None:
    compiler = shlex.split(os.environ.get("CC", "cc"))
    if not compiler or shutil.which(compiler[0]) is None:
        pytest.skip("a C compiler is required to check the public field-kernel header")

    probe = tmp_path / "field_kernel_abi_probe.c"
    probe.write_text(
        r'''
#include "beach_field_kernel.h"

_Static_assert(BEACH_FIELD_KERNEL_ABI_MAJOR == 2, "ABI major changed");
_Static_assert(BEACH_FIELD_KERNEL_ABI_MINOR == 1, "ABI minor changed");
_Static_assert(BEACH_KERNEL_OK == 0, "status code changed");
_Static_assert(BEACH_KERNEL_INVALID_HANDLE == 1, "status code changed");
_Static_assert(BEACH_KERNEL_INVALID_ARGUMENT == 2, "status code changed");
_Static_assert(BEACH_KERNEL_NOT_READY == 3, "status code changed");
_Static_assert(BEACH_KERNEL_FAR_AUTO == 0, "far-correction code changed");
_Static_assert(BEACH_KERNEL_FAR_NONE == 1, "far-correction code changed");
_Static_assert(BEACH_KERNEL_FAR_RESERVED_2 == 2, "reserved code changed");
_Static_assert(BEACH_KERNEL_FAR_CACHED_KNEQ0 == 3, "far-correction code changed");
_Static_assert(_Generic((beach_kernel_handle)0, void *: 1, default: 0),
               "handle representation changed");

#define CHECK_SIGNATURE(name, type) \
  _Static_assert(_Generic(&(name), type: 1, default: 0), #name " signature changed")

typedef int (*version_fn)(int *, int *);
typedef int (*build_info_fn)(char *, int, int *);
typedef int (*create_fn)(beach_kernel_handle *);
typedef int (*handle_fn)(beach_kernel_handle);
typedef int (*set_cache_fn)(beach_kernel_handle, const char *, int, double);
typedef int (*get_cache_fn)(beach_kernel_handle, int *, int *, char *, int,
                            int *, char *, int, int *);
typedef int (*build_fn)(beach_kernel_handle, int, const double *,
                        const double *, const double *, double, int, int, int,
                        const int *, const double *, int, int, double, int,
                        const double *, const double *);
typedef int (*update_fn)(beach_kernel_handle, int, const double *);
typedef int (*eval_fn)(beach_kernel_handle, int, const double *, double *);
typedef int (*force_fn)(beach_kernel_handle, int, const double *,
                        const double *, const double *, double *, double *);

CHECK_SIGNATURE(beach_kernel_get_abi_version, version_fn);
CHECK_SIGNATURE(beach_kernel_get_build_info, build_info_fn);
CHECK_SIGNATURE(beach_kernel_create, create_fn);
CHECK_SIGNATURE(beach_kernel_destroy, handle_fn);
CHECK_SIGNATURE(beach_kernel_set_periodic_cache, set_cache_fn);
CHECK_SIGNATURE(beach_kernel_get_periodic_cache_info, get_cache_fn);
CHECK_SIGNATURE(beach_kernel_build, build_fn);
CHECK_SIGNATURE(beach_kernel_update_charges, update_fn);
CHECK_SIGNATURE(beach_kernel_eval_e, eval_fn);
CHECK_SIGNATURE(beach_kernel_eval_phi, eval_fn);
CHECK_SIGNATURE(beach_kernel_eval_e_direct, eval_fn);
CHECK_SIGNATURE(beach_kernel_eval_phi_direct, eval_fn);
CHECK_SIGNATURE(beach_kernel_force_on_charges, force_fn);

int main(void) { return 0; }
''',
        encoding="utf-8",
    )

    completed = subprocess.run(
        [
            *compiler,
            "-std=c11",
            "-Werror",
            "-fsyntax-only",
            f"-I{PUBLIC_HEADER.parent}",
            str(probe),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr


def test_python_loader_configures_the_frozen_ctypes_abi(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    library = SimpleNamespace(
        **{symbol: _FakeFunction() for symbol in PUBLIC_SYMBOLS}
    )
    monkeypatch.setattr(kernel_module, "_validate_library_abi", lambda _lib: None)

    kernel_module._configure_library(library)  # type: ignore[arg-type]

    for symbol, expected_argtypes in CTYPES_SIGNATURES.items():
        function = getattr(library, symbol)
        assert tuple(function.argtypes) == expected_argtypes
        assert function.restype is C_INT
