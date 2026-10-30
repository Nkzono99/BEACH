import ctypes
import gc
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import beach.native.periodic_zero_mode as zero_mode_module
from beach import FieldKernelError, PeriodicZeroMode


EPS0 = 8.8541878128e-12


def _fake_zero_mode_library() -> tuple[SimpleNamespace, list[object]]:
    destroy_calls: list[object] = []

    def create(handle_out: object) -> int:
        ctypes.cast(handle_out, ctypes.POINTER(ctypes.c_void_p))[0] = (
            ctypes.c_void_p(1234)
        )
        return 0

    def destroy(handle: object) -> int:
        destroy_calls.append(handle)
        return 0

    def ok(*_args: object) -> int:
        return 0

    library = SimpleNamespace(
        beach_zero_mode_create=create,
        beach_zero_mode_destroy=destroy,
        beach_zero_mode_build=ok,
        beach_zero_mode_update=ok,
        beach_zero_mode_eval=ok,
    )
    return library, destroy_calls


def test_periodic_zero_mode_matches_two_sheet_solution_and_traces(
    field_kernel_library: Path,
) -> None:
    area = 2.0
    heights = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    charges = EPS0 * area * np.array([1.0, -1.0])
    updated_charges = EPS0 * area * np.array([2.0, -0.5])

    with PeriodicZeroMode(
        heights,
        charges,
        area,
        library_path=field_kernel_library,
    ) as zero:
        phi, ez = zero.eval(np.array([-0.5, 0.0, 0.5, 1.0, 1.5]))
        _, ez_minus = zero.eval(np.array([0.0, 1.0]), trace="minus")
        _, ez_plus = zero.eval(np.array([0.0, 1.0]), trace="plus")
        zero.update_charges(updated_charges)
        updated_phi, updated_ez = zero.eval(np.array([-0.5, 0.5, 1.5]))

    np.testing.assert_allclose(
        phi, [0.0, 0.0, -0.5, -1.0, -1.0], rtol=0.0, atol=1.0e-14,
        err_msg="initial potential",
    )
    np.testing.assert_allclose(
        ez, [0.0, 0.5, 1.0, 0.5, 0.0], rtol=0.0, atol=1.0e-14,
        err_msg="initial field",
    )
    np.testing.assert_allclose(ez_minus, [0.0, 1.0], rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(ez_plus, [1.0, 0.0], rtol=0.0, atol=1.0e-14)
    assert not phi.flags.writeable
    assert not ez.flags.writeable
    # The charge update keeps the height plan and gives the new two-sheet solution.
    np.testing.assert_allclose(
        updated_phi, [0.0, -1.0, -2.75], rtol=0.0, atol=1.0e-14,
        err_msg="potential after charge update",
    )
    np.testing.assert_allclose(
        updated_ez, [0.0, 2.0, 1.5], rtol=0.0, atol=1.0e-14,
        err_msg="field after charge update",
    )


def test_periodic_zero_mode_rejects_invalid_shape_without_native_library() -> None:
    with pytest.raises(ValueError, match="shape"):
        PeriodicZeroMode(np.zeros((2, 2)), np.ones(2), 1.0)


def test_periodic_zero_mode_rejects_invalid_trace_and_closed_use(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    library, destroy_calls = _fake_zero_mode_library()

    def forbidden_eval(*_args: object) -> int:
        raise AssertionError("rejected evaluations must not reach the native library")

    library.beach_zero_mode_eval = forbidden_eval
    monkeypatch.setattr(
        zero_mode_module,
        "_load_kernel_library",
        lambda _path: library,
    )

    zero = PeriodicZeroMode(np.zeros((1, 3)), np.ones(1), 1.0)
    with pytest.raises(ValueError, match="trace"):
        zero.eval(np.array([0.0]), trace="surface")
    zero.close()
    zero.close()
    assert len(destroy_calls) == 1
    with pytest.raises(FieldKernelError, match="closed"):
        zero.eval(np.array([0.0]))


def test_periodic_zero_mode_finalizer_releases_unclosed_native_handle(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    library, destroy_calls = _fake_zero_mode_library()
    monkeypatch.setattr(
        zero_mode_module,
        "_load_kernel_library",
        lambda _path: library,
    )

    zero = PeriodicZeroMode(np.zeros((1, 3)), np.ones(1), 1.0)
    del zero
    gc.collect()

    assert len(destroy_calls) == 1


def test_periodic_zero_mode_update_override_is_local_and_preserves_gauge(
    field_kernel_library: Path,
) -> None:
    z = np.array([-1.0, 1.0])
    with PeriodicZeroMode(
        np.zeros((1, 3)),
        np.zeros(1),
        1.0,
        e_bottom_V_m=1.25,
        z_gauge_m=0.25,
        phi_gauge_V=3.0,
        library_path=field_kernel_library,
    ) as zero:
        default_phi, default_ez = zero.eval(z)
        zero.update_charges(np.zeros(1), e_bottom_V_m=-2.5)
        override_phi, override_ez = zero.eval(z)
        zero.update_charges(np.zeros(1))
        restored_phi, restored_ez = zero.eval(z)

    np.testing.assert_allclose(
        default_phi, [4.5625, 2.0625], rtol=0.0, atol=1.0e-14
    )
    np.testing.assert_allclose(default_ez, [1.25, 1.25], rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(
        override_phi, [-0.125, 4.875], rtol=0.0, atol=1.0e-14
    )
    np.testing.assert_allclose(override_ez, [-2.5, -2.5], rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(
        restored_phi, [4.5625, 2.0625], rtol=0.0, atol=1.0e-14
    )
    np.testing.assert_allclose(restored_ez, [1.25, 1.25], rtol=0.0, atol=1.0e-14)
