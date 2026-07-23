from __future__ import annotations

import re
from pathlib import Path

import beach.fortran_results.kernel as kernel_module


ROOT = Path(__file__).resolve().parents[2]
FORTRAN_SOURCE = ROOT / "src/physics/field_solver/bem_field_kernel_c.f90"
PUBLIC_HEADER = ROOT / "beach/include/beach_field_kernel.h"


def test_public_header_covers_every_fortran_c_binding() -> None:
    fortran = FORTRAN_SOURCE.read_text(encoding="utf-8")
    header = PUBLIC_HEADER.read_text(encoding="utf-8")

    bound_symbols = set(re.findall(r"bind\(C,\s*name='([^']+)'\)", fortran))
    declared_symbols = set(re.findall(r"\bint\s+(beach_kernel_\w+)\s*\(", header))

    assert declared_symbols == bound_symbols


def test_abi_status_and_far_correction_codes_are_synchronized() -> None:
    fortran = FORTRAN_SOURCE.read_text(encoding="utf-8")
    header = PUBLIC_HEADER.read_text(encoding="utf-8")

    assert (
        f"BEACH_FIELD_KERNEL_ABI_MAJOR {kernel_module.FIELD_KERNEL_ABI_MAJOR}" in header
    )
    assert (
        f"BEACH_FIELD_KERNEL_ABI_MINOR {kernel_module.FIELD_KERNEL_ABI_MINOR}" in header
    )
    assert (
        f"beach_kernel_abi_major = {kernel_module.FIELD_KERNEL_ABI_MAJOR}_c_int"
        in fortran
    )
    assert (
        f"beach_kernel_abi_minor = {kernel_module.FIELD_KERNEL_ABI_MINOR}_c_int"
        in fortran
    )

    for status, name in {
        0: "OK",
        1: "INVALID_HANDLE",
        2: "INVALID_ARGUMENT",
        3: "NOT_READY",
    }.items():
        assert status in kernel_module._STATUS_MESSAGES
        assert f"BEACH_KERNEL_{name} = {status}" in header

    for name, code in kernel_module._FAR_CORRECTION_CODES.items():
        macro_name = name.upper()
        assert f"BEACH_KERNEL_FAR_{macro_name} = {code}" in header
    assert "BEACH_KERNEL_FAR_RESERVED_2 = 2" in header
    assert "M2L_ROOT_ORACLE" not in header
