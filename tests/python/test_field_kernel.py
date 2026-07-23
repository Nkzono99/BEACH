import ctypes
from pathlib import Path

import numpy as np
import pytest

import beach.fortran_results.kernel as kernel_module
from beach.cli.kernel_forces import main as kernel_forces_main
from beach import (
    BeachScene,
    FieldKernel,
    FieldKernelBuildInfo,
    FieldKernelDiagnostics,
    FieldKernelError,
    FieldKernelOptions,
    FortranRunResult,
    calc_object_forces_kernel,
    field_kernel_build_info,
    field_kernel_options_from_result,
)
from beach.fortran_results.constants import K_COULOMB
from beach.fortran_results.panel_quadrature import panel_target_quadrature
from beach.fortran_results.potential import (
    _auto_periodic2_from_result,
    _coerce_periodic2,
)


def _int_value(value: object) -> int:
    return int(getattr(value, "value", value))


class _FakeFunction:
    def __init__(self, name: str, callback, events: list[str]) -> None:  # type: ignore[no-untyped-def]
        self.name = name
        self.callback = callback
        self.events = events
        self.calls: list[tuple[object, ...]] = []
        self.argtypes = None
        self.restype = None

    def __call__(self, *args: object) -> int:
        self.calls.append(args)
        self.events.append(self.name)
        return int(self.callback(*args))


class _FakeKernelLibrary:
    def __init__(
        self,
        *,
        cache_setter: bool,
        cache_getter: bool,
        diagnostics: tuple[int, int, bytes, bytes] = (0, 0, b"", b""),
        build_info: bytes | None = None,
        abi_version: tuple[int, int] | None = None,
    ) -> None:
        self.events: list[str] = []
        self.cache_settings: list[tuple[bytes, float]] = []
        self.diagnostics_value = diagnostics

        def ok(*_args: object) -> int:
            return 0

        def create(handle_out: object) -> int:
            ctypes.cast(handle_out, ctypes.POINTER(ctypes.c_void_p))[0] = (
                ctypes.c_void_p(1234)
            )
            return 0

        def set_cache(
            _handle: object, path_ptr: object, path_len: object, tolerance: object
        ) -> int:
            path = ctypes.string_at(path_ptr, _int_value(path_len))
            self.cache_settings.append(
                (path, float(getattr(tolerance, "value", tolerance)))
            )
            return 0

        def get_cache(
            _handle: object,
            hit_ptr: object,
            count_ptr: object,
            fingerprint_ptr: object,
            fingerprint_capacity: object,
            fingerprint_length_ptr: object,
            path_ptr: object,
            path_capacity: object,
            path_length_ptr: object,
        ) -> int:
            hit, count, fingerprint, path = self.diagnostics_value
            ctypes.cast(hit_ptr, ctypes.POINTER(ctypes.c_int))[0] = hit
            ctypes.cast(count_ptr, ctypes.POINTER(ctypes.c_int))[0] = count
            ctypes.cast(fingerprint_length_ptr, ctypes.POINTER(ctypes.c_int))[0] = len(
                fingerprint
            )
            ctypes.cast(path_length_ptr, ctypes.POINTER(ctypes.c_int))[0] = len(path)
            if (
                _int_value(fingerprint_capacity) <= len(fingerprint)
                or _int_value(path_capacity) <= len(path)
            ):
                return 2
            ctypes.memmove(fingerprint_ptr, fingerprint + b"\0", len(fingerprint) + 1)
            ctypes.memmove(path_ptr, path + b"\0", len(path) + 1)
            return 0

        def get_build_info(
            buffer_ptr: object,
            buffer_capacity: object,
            length_ptr: object,
        ) -> int:
            assert build_info is not None
            ctypes.cast(length_ptr, ctypes.POINTER(ctypes.c_int))[0] = len(build_info)
            if _int_value(buffer_capacity) <= len(build_info):
                return 2
            ctypes.memmove(buffer_ptr, build_info + b"\0", len(build_info) + 1)
            return 0

        def get_abi_version(major_ptr: object, minor_ptr: object) -> int:
            assert abi_version is not None
            major, minor = abi_version
            ctypes.cast(major_ptr, ctypes.POINTER(ctypes.c_int))[0] = major
            ctypes.cast(minor_ptr, ctypes.POINTER(ctypes.c_int))[0] = minor
            return 0

        self.beach_kernel_create = _FakeFunction("create", create, self.events)
        self.beach_kernel_destroy = _FakeFunction("destroy", ok, self.events)
        self.beach_kernel_build = _FakeFunction("build", ok, self.events)
        self.beach_kernel_build_panel = _FakeFunction("build_panel", ok, self.events)
        self.beach_kernel_update_charges = _FakeFunction("update", ok, self.events)
        self.beach_kernel_eval_e = _FakeFunction("eval_e", ok, self.events)
        self.beach_kernel_eval_phi = _FakeFunction("eval_phi", ok, self.events)
        self.beach_kernel_force_on_charges = _FakeFunction("force", ok, self.events)
        if cache_setter:
            self.beach_kernel_set_periodic_cache = _FakeFunction(
                "set_cache", set_cache, self.events
            )
        if cache_getter:
            self.beach_kernel_get_periodic_cache_info = _FakeFunction(
                "get_cache", get_cache, self.events
            )
        if build_info is not None:
            self.beach_kernel_get_build_info = _FakeFunction(
                "get_build_info", get_build_info, self.events
            )
        if abi_version is not None:
            self.beach_kernel_get_abi_version = _FakeFunction(
                "get_abi_version", get_abi_version, self.events
            )


def _install_fake_kernel(
    monkeypatch: pytest.MonkeyPatch, lib: _FakeKernelLibrary
) -> None:
    monkeypatch.setattr(kernel_module, "_load_kernel_library", lambda _path: lib)


def _one_source() -> tuple[np.ndarray, np.ndarray]:
    return np.array([[0.0, 0.0, 0.0]], dtype=float), np.array([0.0], dtype=float)


def _kernel_lib() -> Path:
    path = Path("build/libbeach_field_kernel.so")
    if not path.exists():
        pytest.skip("field kernel shared library is not built; run `make build-kernel`")
    return path


def test_field_kernel_diagnostics_is_a_frozen_top_level_api() -> None:
    assert FieldKernelDiagnostics.__dataclass_params__.frozen


def test_field_kernel_accepts_matching_or_newer_minor_abi(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(
        cache_setter=False,
        cache_getter=False,
        abi_version=(kernel_module.FIELD_KERNEL_ABI_MAJOR, 99),
    )
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()

    with FieldKernel(source_pos, source_q):
        pass

    assert "get_abi_version" in lib.events


def test_field_kernel_rejects_incompatible_major_abi(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(
        cache_setter=False,
        cache_getter=False,
        abi_version=(kernel_module.FIELD_KERNEL_ABI_MAJOR + 1, 0),
    )
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()

    with pytest.raises(FieldKernelError, match="ABI is incompatible"):
        FieldKernel(source_pos, source_q)


def test_field_kernel_build_info_is_a_frozen_top_level_api(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = (
        b"build_info_schema_version=1\n"
        b"build_version=1.5.0-test\n"
        b"build_version_mode=git\n"
        b"build_source_commit=0123456789abcdef0123456789abcdef01234567\n"
        b"build_id=0123456789abcdef0123456789abcdef01234567:clean"
    )
    lib = _FakeKernelLibrary(
        cache_setter=False,
        cache_getter=False,
        build_info=payload,
    )
    _install_fake_kernel(monkeypatch, lib)

    info = field_kernel_build_info("unused.so")

    assert FieldKernelBuildInfo.__dataclass_params__.frozen
    assert info == FieldKernelBuildInfo(
        schema_version=1,
        version="1.5.0-test",
        version_mode="git",
        source_commit="0123456789abcdef0123456789abcdef01234567",
        build_id="0123456789abcdef0123456789abcdef01234567:clean",
    )


def test_field_kernel_build_info_rejects_reordered_payload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = (
        b"build_version=1.5.0-test\n"
        b"build_info_schema_version=1\n"
        b"build_version_mode=git\n"
        b"build_source_commit=0123456789abcdef0123456789abcdef01234567\n"
        b"build_id=0123456789abcdef0123456789abcdef01234567:clean"
    )
    lib = _FakeKernelLibrary(
        cache_setter=False,
        cache_getter=False,
        build_info=payload,
    )
    _install_fake_kernel(monkeypatch, lib)

    with pytest.raises(FieldKernelError, match="unexpected fields"):
        field_kernel_build_info("reordered.so")


def test_field_kernel_build_info_requires_native_attestation_symbol(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=False, cache_getter=False)
    _install_fake_kernel(monkeypatch, lib)

    with pytest.raises(FieldKernelError, match="build-info ABI"):
        field_kernel_build_info("legacy.so")


def test_field_kernel_cache_options_are_set_before_geometry_build(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=True, cache_getter=True)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()
    options = FieldKernelOptions(
        periodic_cache_dir="cache/\N{LATIN SMALL LETTER E WITH ACUTE}",
        periodic_generation_tolerance=2.5e-9,
    )

    with FieldKernel(source_pos, source_q, options=options):
        pass

    assert lib.cache_settings == [
        ("cache/\N{LATIN SMALL LETTER E WITH ACUTE}".encode(), 2.5e-9)
    ]
    assert lib.events.index("set_cache") < lib.events.index("build")


def test_field_kernel_legacy_library_keeps_default_free_build(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=False, cache_getter=False)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()

    with FieldKernel(source_pos, source_q):
        pass

    assert "build" in lib.events


def test_field_kernel_legacy_library_keeps_existing_eval_but_direct_is_optional(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=False, cache_getter=False)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()

    with FieldKernel(source_pos, source_q) as kernel:
        np.testing.assert_array_equal(kernel.eval_e(np.array([[1.0, 0.0, 0.0]])), 0.0)
        with pytest.raises(FieldKernelError, match="exact-direct"):
            kernel.eval_e_direct(np.array([[1.0, 0.0, 0.0]]))
        with pytest.raises(FieldKernelError, match="exact-direct"):
            kernel.eval_phi_direct(np.array([[1.0, 0.0, 0.0]]))


def test_field_kernel_direct_methods_reject_periodic_plan_before_native_call(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=False, cache_getter=False)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()
    options = FieldKernelOptions(
        periodic2=((0, 1), (2.0, 2.0), (0.0, 0.0), 1, "none", 0.0, 4),
        box_min=(0.0, 0.0, -1.0),
        box_max=(2.0, 2.0, 1.0),
    )

    with FieldKernel(source_pos, source_q, options=options) as kernel:
        with pytest.raises(FieldKernelError, match="non-periodic"):
            kernel.eval_e_direct(np.array([[1.0, 0.0, 0.0]]))
        with pytest.raises(FieldKernelError, match="non-periodic"):
            kernel.eval_phi_direct(np.array([[1.0, 0.0, 0.0]]))


def test_field_kernel_legacy_library_keeps_default_finite_periodic_build(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=False, cache_getter=False)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()
    options = FieldKernelOptions(
        periodic2=((0, 1), (2.0, 2.0), (0.0, 0.0), 1, "none", 0.0, 4),
        box_min=(0.0, 0.0, -1.0),
        box_max=(2.0, 2.0, 1.0),
    )

    with FieldKernel(source_pos, source_q, options=options):
        pass

    assert "build" in lib.events


@pytest.mark.parametrize(
    "options, error_model",
    [
        (
            FieldKernelOptions(
                periodic2=((0, 1), (2.0, 2.0), (0.0, 0.0), 1, "cached_kneq0", 0.0, 4),
                box_min=(0.0, 0.0, -1.0),
                box_max=(2.0, 2.0, 1.0),
            ),
            "cached_kneq0",
        ),
        (
            FieldKernelOptions(periodic_cache_dir="custom-cache"),
            "periodic cache configuration",
        ),
        (
            FieldKernelOptions(periodic_generation_tolerance=2.5e-9),
            "periodic cache configuration",
        ),
    ],
)
def test_field_kernel_legacy_library_rejects_models_requiring_cache_setter(
    monkeypatch: pytest.MonkeyPatch,
    options: FieldKernelOptions,
    error_model: str,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=False, cache_getter=False)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()

    with pytest.raises(FieldKernelError, match=error_model):
        FieldKernel(source_pos, source_q, options=options)

    assert "build" not in lib.events


@pytest.mark.parametrize(
    "native, expected_hit, expected_count, expected_fingerprint, expected_path",
    [
        (
            (0, 1, b"c0ffee1234abcdef", b"/tmp/cache/operator.bin"),
            False,
            1,
            "c0ffee1234abcdef",
            Path("/tmp/cache/operator.bin"),
        ),
        (
            (1, 0, b"c0ffee1234abcdef", b"/tmp/cache/operator.bin"),
            True,
            0,
            "c0ffee1234abcdef",
            Path("/tmp/cache/operator.bin"),
        ),
        ((0, 0, b"", b""), None, 0, None, None),
    ],
)
def test_field_kernel_diagnostics_decodes_cache_metadata(
    monkeypatch: pytest.MonkeyPatch,
    native: tuple[int, int, bytes, bytes],
    expected_hit: bool | None,
    expected_count: int,
    expected_fingerprint: str | None,
    expected_path: Path | None,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=True, cache_getter=True, diagnostics=native)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()

    with FieldKernel(source_pos, source_q) as kernel:
        diagnostics = kernel.diagnostics()

    assert diagnostics.periodic_cache_hit is expected_hit
    assert diagnostics.periodic_operator_build_count == expected_count
    assert diagnostics.periodic_cache_fingerprint == expected_fingerprint
    assert diagnostics.periodic_cache_path == expected_path
    assert len(lib.beach_kernel_get_periodic_cache_info.calls) == 2


def test_field_kernel_diagnostics_requires_new_symbol(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lib = _FakeKernelLibrary(cache_setter=False, cache_getter=False)
    _install_fake_kernel(monkeypatch, lib)
    source_pos, source_q = _one_source()

    with FieldKernel(source_pos, source_q) as kernel:
        with pytest.raises(FieldKernelError, match="field-kernel diagnostics"):
            kernel.diagnostics()


def test_removed_root_oracle_is_rejected() -> None:
    with pytest.raises(ValueError, match='was removed; use "cached_kneq0"'):
        kernel_module._far_correction_code("m2l_root_oracle")


def test_coerce_periodic2_validates_legacy_tuple_and_cached_policy() -> None:
    cached = ((0, 1), (2.0, 2.0), (0.0, 0.0), 1, "cached_kneq0", 0.0, 4)

    with pytest.raises(ValueError, match="far_correction"):
        _coerce_periodic2(cached)
    assert _coerce_periodic2(cached, allow_cached_kneq0=True) == cached
    with pytest.raises(ValueError, match="distinct axis"):
        _coerce_periodic2(((0, 0), (2.0, 2.0), (0.0, 0.0), 1, "none", 0.0, 4))
    with pytest.raises(ValueError, match="positive"):
        _coerce_periodic2(((0, 1), (-2.0, 2.0), (0.0, 0.0), 1, "none", 0.0, 4))


@pytest.mark.parametrize(
    "periodic2, message",
    [
        (
            ((0, 1), (float("nan"), 2.0), (0.0, 0.0), 1, "none", 0.0, 4),
            "finite and positive",
        ),
        (
            ((0, 1), (2.0, 2.0), (float("inf"), 0.0), 1, "none", 0.0, 4),
            "origins must be finite",
        ),
        (
            ((0, 1), (2.0, 2.0), (0.0, 0.0), -1, "none", 0.0, 4),
            "image_layers",
        ),
        (
            ((0, 1), (2.0, 2.0), (0.0, 0.0), 1, "none", float("nan"), 4),
            "ewald_alpha",
        ),
        (
            ((0, 1), (2.0, 2.0), (0.0, 0.0), 1, "none", 0.0, -1),
            "ewald_layers",
        ),
    ],
)
def test_coerce_periodic2_tuple_uses_mapping_validation(
    periodic2: tuple[object, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        _coerce_periodic2(periodic2)  # type: ignore[arg-type]


def test_native_field_kernel_resolver_allows_cached_config_but_python_auto_fails_closed(
    tmp_path: Path,
) -> None:
    (tmp_path / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [2.0, 2.0, 1.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                'field_periodic_far_correction = "cached_kneq0"',
                'field_periodic_cache_dir = "custom-cache"',
                "field_periodic_generation_tolerance = 2.5e-9",
            ]
        ),
        encoding="utf-8",
    )
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0]),
        triangles=np.array(
            [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]
        ),
    )

    options = field_kernel_options_from_result(result)

    assert options.periodic2 is not None
    assert options.periodic2[4] == "cached_kneq0"
    assert options.periodic_cache_dir == "custom-cache"
    assert options.periodic_generation_tolerance == pytest.approx(2.5e-9)
    with pytest.raises(ValueError, match="far_correction"):
        _auto_periodic2_from_result(result)


def test_field_kernel_eval_e_matches_direct_sum() -> None:
    lib = _kernel_lib()
    source_pos = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    source_q = np.array([1.0e-9, -2.0e-9], dtype=float)
    target = np.array([[0.0, 1.0, 0.0]], dtype=float)

    with FieldKernel(source_pos, source_q, library_path=lib) as kernel:
        field = kernel.eval_e(target)

    delta = target[0] - source_pos
    expected = K_COULOMB * np.sum(
        source_q[:, None] * delta / (np.linalg.norm(delta, axis=1)[:, None] ** 3),
        axis=0,
    )
    np.testing.assert_allclose(field[0], expected, rtol=1.0e-14, atol=1.0e-14)


def test_field_kernel_exact_direct_matches_softened_sum_and_charge_refresh() -> None:
    lib = _kernel_lib()
    axis = (-0.6, 0.0, 0.6)
    source_pos = np.array(
        [[x, y, z] for x in axis for y in axis for z in axis], dtype=float
    )
    source_q = (
        np.where(np.arange(source_pos.shape[0]) % 2 == 0, 1.0, -1.0)
        * (1.0 + 0.07 * np.arange(source_pos.shape[0]))
        * 1.0e-9
    )
    targets = np.array([[0.13, -0.08, 0.04], source_pos[0]], dtype=float)
    softening = 0.15
    options = FieldKernelOptions(
        theta=20.0,
        leaf_max=1,
        order=0,
        softening=softening,
    )

    delta = targets[:, None, :] - source_pos[None, :, :]
    r2 = np.sum(delta * delta, axis=2) + softening**2
    expected_e = K_COULOMB * np.sum(
        source_q[None, :, None] * delta / (r2[:, :, None] ** 1.5), axis=1
    )
    expected_phi = K_COULOMB * np.sum(source_q[None, :] / np.sqrt(r2), axis=1)

    with FieldKernel(source_pos, source_q, options=options, library_path=lib) as kernel:
        ordinary_e = kernel.eval_e(targets)
        direct_e = kernel.eval_e_direct(targets)
        direct_phi = kernel.eval_phi_direct(targets)
        kernel.update_charges(-0.25 * source_q)
        refreshed_e = kernel.eval_e_direct(targets)
        refreshed_phi = kernel.eval_phi_direct(targets)

    np.testing.assert_allclose(direct_e, expected_e, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(direct_phi, expected_phi, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(refreshed_e, -0.25 * expected_e, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(refreshed_phi, -0.25 * expected_phi, rtol=2.0e-14, atol=1.0e-14)
    assert np.linalg.norm(ordinary_e[0] - expected_e[0]) > 1.0e-6 * np.linalg.norm(
        expected_e[0]
    )


def test_field_kernel_from_panel_result_dispatches_triangle_geometry(tmp_path: Path) -> None:
    lib = _kernel_lib()
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
        dtype=float,
    )
    charges = np.array([2.5e-9], dtype=float)
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        field_source_model="triangle_p0",
        field_kernel_id="triangle_p0_exact_p2m_near",
    )
    target = np.array([[0.25, 0.2, 0.4]], dtype=float)

    with FieldKernel.from_result(result, step=None, library_path=lib) as panel_kernel:
        panel_field = panel_kernel.eval_e(target)
    with FieldKernel(triangles.mean(axis=1), charges, library_path=lib) as point_kernel:
        point_field = point_kernel.eval_e(target)

    assert np.linalg.norm(panel_field - point_field) > 1.0e-6 * np.linalg.norm(panel_field)


def test_panel_direct_and_ordinary_use_same_principal_value_at_target_quadrature() -> None:
    lib = _kernel_lib()
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
        dtype=float,
    )
    charges = np.array([2.0e-9])
    points, _, _ = panel_target_quadrature(triangles, charges, order=7)

    with FieldKernel(
        triangles.mean(axis=1),
        charges,
        source_triangles=triangles,
        library_path=lib,
    ) as kernel:
        ordinary_e = kernel.eval_e(points)
        direct_e = kernel.eval_e_direct(points)
        ordinary_phi = kernel.eval_phi(points)
        direct_phi = kernel.eval_phi_direct(points)

    np.testing.assert_allclose(ordinary_e, direct_e, rtol=2.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(ordinary_phi, direct_phi, rtol=2.0e-14, atol=1.0e-14)


def test_field_kernel_adds_uniform_external_e0() -> None:
    lib = _kernel_lib()
    source_pos = np.array([[0.0, 0.0, 0.0]], dtype=float)
    source_q = np.array([0.0], dtype=float)
    target = np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0]], dtype=float)
    target_q = np.array([2.0, 3.0], dtype=float)

    with FieldKernel(
        source_pos,
        source_q,
        options=FieldKernelOptions(external_e0=(1.0, 0.0, 0.0)),
        library_path=lib,
    ) as kernel:
        field = kernel.eval_e(target)
        direct_field = kernel.eval_e_direct(target)
        force, torque = kernel.force_on_charges(target, target_q, origin=(0.0, 0.0, 0.0))

    np.testing.assert_allclose(field, np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    np.testing.assert_array_equal(direct_field, np.zeros((2, 3)))
    np.testing.assert_allclose(force, np.array([5.0, 0.0, 0.0]))
    np.testing.assert_allclose(torque, np.array([0.0, 0.0, -8.0]))


def test_calc_object_forces_kernel_excludes_self_sources(tmp_path: Path) -> None:
    lib = _kernel_lib()
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    records = calc_object_forces_kernel(result, step=None, library_path=lib)

    centers = triangles.mean(axis=1)
    delta_12 = centers[0] - centers[1]
    expected_1 = K_COULOMB * charges[0] * charges[1] * delta_12 / (
        np.linalg.norm(delta_12) ** 3
    )
    assert [record.mesh_id for record in records] == [1, 2]
    np.testing.assert_allclose(records[0].force_N, expected_1, rtol=1.0e-14)
    np.testing.assert_allclose(records[1].force_N, -expected_1, rtol=1.0e-14)
    assert records[0].total_charge_C == pytest.approx(charges[0])


def test_calc_object_forces_kernel_uses_explicit_config_softening(tmp_path: Path) -> None:
    lib = _kernel_lib()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    config_path = tmp_path / "beach.toml"
    config_path.write_text("[sim]\nsoftening = 0.5\ntree_leaf_max = 8\n", encoding="utf-8")
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
    result = FortranRunResult(
        directory=run_dir,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    records = calc_object_forces_kernel(
        result,
        step=None,
        library_path=lib,
        config_path=config_path,
    )

    centers = triangles.mean(axis=1)
    delta_12 = centers[0] - centers[1]
    soft2 = 0.5**2
    expected_1 = K_COULOMB * charges[0] * charges[1] * delta_12 / (
        (np.dot(delta_12, delta_12) + soft2) ** 1.5
    )
    np.testing.assert_allclose(records[0].force_N, expected_1, rtol=1.0e-14)


def test_scene_move_updates_kernel_force_geometry(tmp_path: Path) -> None:
    lib = _kernel_lib()
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9], dtype=float)
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=2,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=np.array([1, 2], dtype=np.int64),
    )

    scene = BeachScene.from_result(result, step=None)
    moved = scene.move(2, by=[1.0, 0.0, 0.0])
    records = moved.calc_object_forces_kernel(
        target_mesh_ids=1,
        library_path=lib,
    )

    delta_12 = moved.centers[0] - moved.centers[1]
    expected_1 = K_COULOMB * charges[0] * charges[1] * delta_12 / (
        np.linalg.norm(delta_12) ** 3
    )
    np.testing.assert_allclose(records[0].force_N, expected_1, rtol=1.0e-14)
    np.testing.assert_allclose(scene.centers[1], triangles[1].mean(axis=0))


def test_scene_rotate_keeps_charges_attached_to_mesh_elements(tmp_path: Path) -> None:
    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.3]],
            [[1.0, 0.0, 0.0], [1.0, 0.3, 0.0], [1.0, 0.0, 0.3]],
            [[3.0, 0.0, 0.0], [3.0, 0.3, 0.0], [3.0, 0.0, 0.3]],
        ],
        dtype=float,
    )
    charges = np.array([1.0e-9, 2.0e-9, 3.0e-9], dtype=float)
    mesh_ids = np.array([1, 2, 2], dtype=np.int64)
    result = FortranRunResult(
        directory=tmp_path,
        mesh_nelem=3,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=charges,
        triangles=triangles,
        mesh_ids=mesh_ids,
    )

    scene = BeachScene.from_result(result, step=None)
    rotated = scene.rotate(2, axis=[0.0, 0.0, 1.0], angle_deg=90.0)

    mask = mesh_ids == 2
    origin = scene.centers[mask].mean(axis=0)
    angle = np.deg2rad(90.0)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    expected_centers = (scene.centers[mask] - origin) @ rotation.T + origin
    np.testing.assert_allclose(rotated.centers[mask], expected_centers)
    np.testing.assert_allclose(rotated.charges, scene.charges)


def test_kernel_forces_cli_writes_csv(tmp_path: Path) -> None:
    lib = _kernel_lib()
    out = tmp_path / "run"
    out.mkdir()
    config_path = tmp_path / "beach.toml"
    config_path.write_text("[sim]\nsoftening = 0.0\n", encoding="utf-8")
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=0",
                "absorbed=0",
                "escaped=0",
                "batches=0",
                "last_rel_change=0.0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-9\n2,2.0e-9\n",
        encoding="utf-8",
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0.0,0.0,0.0,0.0,0.3,0.0,0.0,0.0,0.3,1.0e-9,1\n"
        "2,1.0,0.0,0.0,1.0,0.3,0.0,1.0,0.0,0.3,2.0e-9,2\n",
        encoding="utf-8",
    )
    save_csv = tmp_path / "forces.csv"

    kernel_forces_main(
        [
            str(out),
            "--config",
            str(config_path),
            "--library",
            str(lib),
            "--save-csv",
            str(save_csv),
        ]
    )

    text = save_csv.read_text(encoding="utf-8")
    assert "mesh_id,step,total_charge_C" in text
    assert "1,-1,1e-09" in text
