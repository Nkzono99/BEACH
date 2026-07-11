from pathlib import Path

import numpy as np
import pytest

from beach import FieldKernelError, PeriodicZeroMode


EPS0 = 8.8541878128e-12


def _kernel_lib() -> Path:
    path = Path("build/libbeach_field_kernel.so")
    if not path.exists():
        pytest.skip("field kernel shared library is not built; run `make build-kernel`")
    return path


def test_periodic_zero_mode_matches_two_sheet_solution_and_traces() -> None:
    area = 2.0
    heights = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    charges = EPS0 * area * np.array([1.0, -1.0])

    with PeriodicZeroMode(
        heights,
        charges,
        area,
        library_path=_kernel_lib(),
    ) as zero:
        phi, ez = zero.eval(np.array([-0.5, 0.0, 0.5, 1.0, 1.5]))
        _, ez_minus = zero.eval(np.array([0.0, 1.0]), trace="minus")
        _, ez_plus = zero.eval(np.array([0.0, 1.0]), trace="plus")

    np.testing.assert_allclose(phi, [0.0, 0.0, -0.5, -1.0, -1.0], atol=1.0e-14)
    np.testing.assert_allclose(ez, [0.0, 0.5, 1.0, 0.5, 0.0], atol=1.0e-14)
    np.testing.assert_allclose(ez_minus, [0.0, 1.0], atol=1.0e-14)
    np.testing.assert_allclose(ez_plus, [1.0, 0.0], atol=1.0e-14)
    assert not phi.flags.writeable
    assert not ez.flags.writeable


def test_periodic_zero_mode_updates_charges_without_rebuilding() -> None:
    area = 2.0
    heights = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    charges = EPS0 * area * np.array([1.0, -1.0])

    with PeriodicZeroMode(
        heights,
        charges,
        area,
        library_path=_kernel_lib(),
    ) as zero:
        _, before = zero.eval(np.array([0.5]))
        zero.update_charges(2.0 * charges)
        _, after = zero.eval(np.array([0.5]))

    np.testing.assert_allclose(after, 2.0 * before)


def test_periodic_zero_mode_rejects_invalid_shape_trace_and_closed_use() -> None:
    with pytest.raises(ValueError, match="shape"):
        PeriodicZeroMode(np.zeros((2, 2)), np.ones(2), 1.0, library_path=_kernel_lib())

    zero = PeriodicZeroMode(
        np.zeros((1, 3)), np.ones(1), 1.0, library_path=_kernel_lib()
    )
    with pytest.raises(ValueError, match="trace"):
        zero.eval(np.array([0.0]), trace="surface")
    zero.close()
    with pytest.raises(FieldKernelError, match="closed"):
        zero.eval(np.array([0.0]))

