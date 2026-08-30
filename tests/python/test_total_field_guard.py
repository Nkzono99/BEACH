from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from beach.fortran_results.context import RunContext
from beach.fortran_results.kernel import (
    FieldKernelOptions,
    _require_total_field_reconstruction,
    field_kernel_options_from_result,
)
from beach.fortran_results.periodic import Periodic2Config
from beach.fortran_results.types import (
    FieldReconstructionReceipt,
    FortranRunResult,
)


def _result(
    directory: Path,
    *,
    field_reconstruction: FieldReconstructionReceipt | None = None,
) -> FortranRunResult:
    return FortranRunResult(
        directory=directory,
        mesh_nelem=1,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.array([0.0]),
        triangles=np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]),
        field_reconstruction=field_reconstruction,
    )


def _field_receipt() -> FieldReconstructionReceipt:
    return FieldReconstructionReceipt(
        schema_version=2,
        resolved_field_solver="fmm",
        fmm_expansion_order=4,
        field_bc_mode="free",
        tree_theta=0.375,
        tree_leaf_max=23,
        external_e0_v_m=(1.0, 2.0, 3.0),
        use_box=True,
        box_min_m=(-1.0, -2.0, -3.0),
        box_max_m=(4.0, 5.0, 6.0),
        boundary_low=(0, 0, 0),
        boundary_high=(0, 0, 0),
        periodic_image_layers=2,
        periodic_far_correction="none",
        periodic_nonzero_mode_backend="not_applicable",
        periodic_zero_mode_policy="not_applicable",
        periodic_lower_boundary_model="not_applicable",
        periodic_reference_mode_layers=4,
        periodic_panel_quadrature_order=12,
        periodic_ewald_alpha=0.25,
        periodic_ewald_layers=5,
        periodic_cache_dir="verified-cache",
        periodic_generation_tolerance=2.0e-9,
    )


def _periodic2(*, far_correction: str = "none") -> Periodic2Config:
    return Periodic2Config(
        axes=(0, 1),
        lengths=(1.0, 1.0),
        origins=(0.0, 0.0),
        image_layers=1,
        far_correction=far_correction,
        ewald_alpha=0.0,
        ewald_layers=4,
    )


def _write_free_field_config(path: Path) -> None:
    path.write_text(
        '[sim]\nfield_solver = "fmm"\n[field_boundary]\nmode = "free"\n',
        encoding="utf-8",
    )


def test_total_field_guard_requires_config(tmp_path: Path) -> None:
    context = RunContext.from_value(_result(tmp_path))

    with pytest.raises(ValueError, match=r"requires .*beach\.toml"):
        _require_total_field_reconstruction(
            context,
            FieldKernelOptions(),
            operation="test operation",
        )


@pytest.mark.parametrize(
    ("options", "match"),
    [
        pytest.param(
            FieldKernelOptions(periodic2=_periodic2(far_correction="cached_kneq0")),
            "only the k!=0 component",
            id="cached-component",
        ),
        pytest.param(
            FieldKernelOptions(
                resolved_field_solver="direct",
                periodic2=_periodic2(),
            ),
            "exact-direct reconstruction.*periodic field plan",
            id="direct-periodic",
        ),
    ],
)
def test_total_field_guard_rejects_incomplete_kernel(
    tmp_path: Path,
    options: FieldKernelOptions,
    match: str,
) -> None:
    _write_free_field_config(tmp_path / "beach.toml")
    context = RunContext.from_value(_result(tmp_path))

    with pytest.raises(ValueError, match=match):
        _require_total_field_reconstruction(
            context,
            options,
            operation="test operation",
        )


def test_total_field_guard_accepts_configured_complete_kernel(tmp_path: Path) -> None:
    _write_free_field_config(tmp_path / "beach.toml")
    context = RunContext.from_value(_result(tmp_path))

    _require_total_field_reconstruction(
        context,
        FieldKernelOptions(),
        operation="test operation",
    )


def test_total_field_guard_accepts_resolved_receipt_without_config(
    tmp_path: Path,
) -> None:
    context = RunContext.from_value(
        _result(tmp_path, field_reconstruction=_field_receipt())
    )

    _require_total_field_reconstruction(
        context,
        FieldKernelOptions(),
        operation="test operation",
    )


def test_total_field_guard_rejects_panel_spectral_receipt(
    tmp_path: Path,
) -> None:
    receipt = replace(
        _field_receipt(),
        field_bc_mode="periodic2",
        boundary_low=(2, 2, 0),
        boundary_high=(2, 2, 0),
        periodic_nonzero_mode_backend="panel_spectral_reference",
        periodic_zero_mode_policy="exclude_k0",
        periodic_lower_boundary_model="e_bottom_zero",
    )
    context = RunContext.from_value(
        _result(tmp_path, field_reconstruction=receipt)
    )

    with pytest.raises(ValueError, match="panel_spectral_reference"):
        _require_total_field_reconstruction(
            context,
            FieldKernelOptions(),
            operation="test operation",
        )
    with pytest.raises(ValueError, match="panel_spectral_reference"):
        field_kernel_options_from_result(context)


@pytest.mark.parametrize("explicit", [False, True])
def test_resolved_receipt_ignores_mismatched_external_config(
    tmp_path: Path,
    explicit: bool,
) -> None:
    config = tmp_path / "beach.toml"
    config.write_text(
        "[sim]\n"
        "e0 = [99.0, 98.0, 97.0]\n"
        "tree_theta = 0.9\n"
        "tree_leaf_max = 2\n"
        "[domain]\n"
        "box_min = [10.0, 10.0, 10.0]\n"
        "box_max = [20.0, 20.0, 20.0]\n",
        encoding="utf-8",
    )
    result = _result(tmp_path, field_reconstruction=_field_receipt())
    config_path = config if explicit else None
    options = field_kernel_options_from_result(
        result,
        config_path=config_path,
    )

    assert options.external_e0 == (1.0, 2.0, 3.0)
    assert options.resolved_field_solver == "fmm"
    assert options.order == 4
    assert options.theta == pytest.approx(0.375)
    assert options.leaf_max == 23
    assert options.box_min == (-1.0, -2.0, -3.0)
    assert options.box_max == (4.0, 5.0, 6.0)
    assert options.periodic_cache_dir == "verified-cache"
    assert options.periodic_generation_tolerance == pytest.approx(2.0e-9)


def test_resolved_receipt_order_is_default_and_explicit_order_overrides(
    tmp_path: Path,
) -> None:
    receipt = replace(_field_receipt(), fmm_expansion_order=7)
    result = _result(tmp_path, field_reconstruction=receipt)

    automatic = field_kernel_options_from_result(result)
    explicit = field_kernel_options_from_result(result, order=5)

    assert automatic.order == 7
    assert explicit.order == 5


def test_resolved_direct_receipt_routes_to_exact_direct_options(
    tmp_path: Path,
) -> None:
    receipt = replace(_field_receipt(), resolved_field_solver="direct")
    result = _result(tmp_path, field_reconstruction=receipt)

    options = field_kernel_options_from_result(result)

    assert options.resolved_field_solver == "direct"
    assert options.periodic2 is None


def test_resolved_treecode_receipt_fails_closed(tmp_path: Path) -> None:
    receipt = replace(_field_receipt(), resolved_field_solver="treecode")
    context = RunContext.from_value(
        _result(tmp_path, field_reconstruction=receipt)
    )

    with pytest.raises(ValueError, match="treecode"):
        field_kernel_options_from_result(context)
    with pytest.raises(ValueError, match="treecode"):
        field_kernel_options_from_result(
            context,
            periodic2={"axes": (0, 1), "lengths": (1.0, 1.0)},
        )
    with pytest.raises(ValueError, match="treecode"):
        _require_total_field_reconstruction(
            context,
            FieldKernelOptions(resolved_field_solver="treecode"),
            operation="test operation",
        )


@pytest.mark.parametrize(
    ("mesh_nelem", "expected"),
    [(255, "direct"), (256, "fmm")],
)
def test_legacy_auto_solver_resolves_from_mesh_threshold(
    tmp_path: Path,
    mesh_nelem: int,
    expected: str,
) -> None:
    (tmp_path / "beach.toml").write_text(
        "[sim]\n"
        'field_solver = "auto"\n'
        'field_bc_mode = "free"\n'
        "tree_min_nelem = 256\n",
        encoding="utf-8",
    )
    result = replace(_result(tmp_path), mesh_nelem=mesh_nelem)

    options = field_kernel_options_from_result(result)

    assert options.resolved_field_solver == expected


def test_total_field_guard_rejects_legacy_panel_spectral_config(
    tmp_path: Path,
) -> None:
    config = tmp_path / "beach.toml"
    config.write_text(
        '[field_boundary]\nmode = "periodic2"\n'
        '[domain]\nbox_min = [0.0, 0.0, 0.0]\n'
        'box_max = [1.0, 1.0, 1.0]\nperiodic_axes = ["x", "y"]\n'
        '[periodic2]\nnonzero_mode_backend = "panel_spectral_reference"\n',
        encoding="utf-8",
    )
    context = RunContext.from_value(_result(tmp_path))

    with pytest.raises(ValueError, match="panel_spectral_reference"):
        _require_total_field_reconstruction(
            context,
            FieldKernelOptions(),
            operation="test operation",
        )
    with pytest.raises(ValueError, match="panel_spectral_reference"):
        field_kernel_options_from_result(context)
