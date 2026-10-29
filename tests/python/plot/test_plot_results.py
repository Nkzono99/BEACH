"""Result plot contracts."""

from dataclasses import replace

from beach.io.types import RunField, RunHistory, RunMesh, RunSummary

from pathlib import Path
import numpy as np
import pytest
from beach import (
    Beach,
    RunResult,
    animate_history_mesh,
    calc_coulomb,
    plot_mesh_source_boxplot,
    plot_charge_mesh,
    plot_coulomb_force_matrix,
    plot_potential_slices,
    plot_potential_mesh,
    compute_potential_mesh,
)
from beach.plot.animation import _select_frame_columns
from result_fixtures import (
    _make_history,
    _write_complete_free_field_config,
    _write_coulomb_matrix_fixture,
    _write_three_mesh_fixture,
    _unit_panel_kernel as _unit_panel_kernel,
)


def test_potential_animation_rejects_cached_kneq0_component(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    (output_dir / "beach.toml").write_text(
        "\n".join(
            [
                "[sim]",
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -2.0]",
                "box_max = [4.0, 5.0, 12.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                'field_periodic_far_correction = "cached_kneq0"',
            ]
        ),
        encoding="utf-8",
    )
    triangles = np.array([[[0.2, 0.2, 0.0], [0.8, 0.2, 0.0], [0.2, 0.8, 0.0]]])
    result = RunResult(
        directory=output_dir,
        charges=np.array([2.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=1,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
        history=RunHistory(
            charge=_make_history(
                mesh_nelem=1,
                history=np.array([[1.0e-9, 2.0e-9]]),
            )
        ),
    )
    with pytest.raises(ValueError, match="only the k!=0 component"):
        animate_history_mesh(
            result,
            quantity="potential",
            reference_point=[2.0, 2.5, 10.0],
        )


def test_plot_potential_slices_returns_figure_and_axes(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    fig, axes = plot_potential_slices(
        result,
        box_min=[-1.0, -1.0, -1.0],
        box_max=[1.0, 1.0, 1.0],
        grid_n=10,
        xy_z=0.3,
        yz_x=0.3,
        xz_y=0.3,
    )

    assert fig is not None
    assert len(axes) == 3
    fig.clf()


def test_plot_potential_slices_applies_vmin_vmax(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    fig, axes = plot_potential_slices(
        result,
        box_min=[-1.0, -1.0, -1.0],
        box_max=[1.0, 1.0, 1.0],
        grid_n=8,
        vmin=-5.0,
        vmax=12.0,
    )

    for ax in axes:
        assert ax.collections
        assert ax.collections[0].get_clim() == (-5.0, 12.0)
    fig.clf()


def test_plot_potential_slices_rejects_invalid_vmin_vmax(tmp_path: Path) -> None:
    _write_complete_free_field_config(tmp_path)
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    with pytest.raises(ValueError, match="vmin"):
        plot_potential_slices(
            result,
            box_min=[-1.0, -1.0, -1.0],
            box_max=[1.0, 1.0, 1.0],
            grid_n=8,
            vmin=1.0,
            vmax=1.0,
        )


def test_plot_charge_mesh_can_apply_periodic2_mesh_wrapping() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    triangles = np.array([[[1.10, 0.10, 0.0], [1.20, 0.10, 0.0], [1.10, 0.20, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    fig, ax = plot_charge_mesh(
        result,
        periodic2={"axes": (0, 1), "lengths": (1.0, 1.0), "origins": (0.0, 0.0)},
        apply_periodic2_mesh=True,
    )

    x_center = 0.5 * sum(ax.get_xlim())
    y_center = 0.5 * sum(ax.get_ylim())
    assert x_center == pytest.approx(0.15)
    assert y_center == pytest.approx(0.15)
    fig.clf()


def test_plot_potential_mesh_returns_figure_and_axes(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9, -1.0e-9]),
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    fig, ax = plot_potential_mesh(result, reference_point=None)

    assert fig is not None
    assert ax is not None
    fig.clf()

    fig, ax = plot_potential_mesh(
        result,
        reference_point=None,
        axis_unit="um",
    )
    assert ax.get_xlabel() == "x [um]"
    fig.clf()


def test_plot_potential_mesh_can_apply_periodic2_mesh_wrapping(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array([[[1.10, 0.10, 0.0], [1.20, 0.10, 0.0], [1.10, 0.20, 0.0]]])
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    fig, ax = plot_potential_mesh(
        result,
        reference_point=None,
        periodic2={"axes": (0, 1), "lengths": (1.0, 1.0), "origins": (0.0, 0.0)},
        apply_periodic2_mesh=True,
    )

    x_center = 0.5 * sum(ax.get_xlim())
    y_center = 0.5 * sum(ax.get_ylim())
    assert x_center == pytest.approx(0.15)
    assert y_center == pytest.approx(0.15)
    fig.clf()


def test_plot_potential_mesh_accepts_custom_view_angles(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9, -1.0e-9]),
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    fig, ax = plot_potential_mesh(
        result,
        reference_point=None,
        view_elev=11.0,
        view_azim=37.0,
    )

    assert ax.elev == pytest.approx(11.0)
    assert ax.azim == pytest.approx(37.0)
    ax.view_init(elev=5.0, azim=15.0)
    assert ax.elev == pytest.approx(5.0)
    assert ax.azim == pytest.approx(15.0)
    fig.clf()


def _right_triangle(height: float, z: float) -> list[list[float]]:
    """Triangle with unit base, so its area is ``height / 2``."""

    return [[0.0, 0.0, z], [1.0, 0.0, z], [0.0, height, z]]


def test_plot_mesh_source_boxplot_charge_uses_area_weighted_statistics() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    # Areas 1, 1, 0.5, 1, 0.5 put the weighted CDF at 0.25, 0.5, 0.625, 0.875, 1 for
    # the sorted values 10, 20, 30, 40, 1000. Unweighted quartiles would be 20/30/40.
    values = np.array([40.0, 10.0, 1000.0, 30.0, 20.0])
    heights = np.array([2.0, 2.0, 1.0, 1.0, 2.0])
    triangles = np.array(
        [_right_triangle(height, float(index)) for index, height in enumerate(heights)]
    )
    result = RunResult(
        directory=Path("dummy"),
        charges=values,
        summary=RunSummary(
            mesh_nelem=5,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles, mesh_ids=np.full(5, 11, dtype=np.int64)),
    )

    fig, ax = plot_mesh_source_boxplot(result, quantity="charge")

    stats = getattr(ax, "_beach_box_stats")
    assert len(stats) == 1
    assert {key: stats[0][key] for key in ("whislo", "q1", "med", "q3", "whishi")} == (
        pytest.approx({"whislo": 10.0, "q1": 10.0, "med": 20.0, "q3": 40.0, "whishi": 40.0})
    )
    np.testing.assert_array_equal(stats[0]["fliers"], [1000.0])
    fig.clf()


def test_plot_mesh_source_boxplot_supports_potential_quantity(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _write_complete_free_field_config(tmp_path)

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 0.8, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9, -2.0e-9]),
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(
            triangles=triangles,
            mesh_ids=np.array([1, 2], dtype=np.int64),
            mesh_potential_v=np.array([5.0, -7.0]),
        ),
    )

    # step=None reads the saved mesh_potential.csv values.
    fig, ax = plot_mesh_source_boxplot(result, quantity="potential", step=None)
    saved = getattr(ax, "_beach_box_stats")
    assert [item["med"] for item in saved] == [5.0, -7.0]
    assert [(item["q1"], item["q3"]) for item in saved] == [(5.0, 5.0), (-7.0, -7.0)]
    fig.clf()

    # A history step recomputes the potential from that step's charges.
    fig, ax = plot_mesh_source_boxplot(result, quantity="potential")
    recomputed = getattr(ax, "_beach_box_stats")
    expected = compute_potential_mesh(result, source="kernel")
    assert np.all(np.abs(expected) > 1.0)
    np.testing.assert_allclose([item["med"] for item in recomputed], expected)
    assert ax.get_ylabel() == "potential [V]"
    fig.clf()


def test_plot_mesh_source_boxplot_rejects_invalid_quantity() -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
    )

    with pytest.raises(ValueError, match="quantity"):
        plot_mesh_source_boxplot(result, quantity="invalid")


def test_beach_plot_mesh_source_boxplot_uses_mesh_source_labels(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    out = tmp_path / "run_source_boxplot"
    out.mkdir()
    _write_three_mesh_fixture(out)

    beach = Beach(out)
    fig, ax = beach.plot_mesh_source_boxplot(quantity="charge")

    stats = getattr(ax, "_beach_box_stats")
    labels = [str(item["label"]) for item in stats]
    assert "id=1 (template/plane/insulator/eps=1)" in labels
    assert "id=2 (template/box/insulator/eps=1)" in labels
    assert "id=3 (template/sphere/insulator/eps=1)" in labels
    fig.clf()


def _force_x(beach: Beach, target: int, sources: list[int]) -> float:
    return float(calc_coulomb(beach, target=target, source=sources).force_on_a_N[0])


def test_plot_coulomb_force_matrix_auto_labels_targets_from_config(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    out = tmp_path / "run_coulomb_matrix"
    out.mkdir()
    _write_coulomb_matrix_fixture(out)

    beach = Beach(out)
    fig, ax = beach.plot_coulomb_force_matrix(component="x")

    matrix_info = getattr(ax, "_beach_coulomb_matrix")
    assert matrix_info["target_labels"] == ("plane", "sphere1", "sphere2")
    assert matrix_info["source_labels"] == ("plane", "sphere1", "sphere2", "net")
    matrix = matrix_info["matrix"]
    # Rows are sources and columns are targets; each entry is calc_coulomb's Fx.
    mesh_ids = [1, 2, 3]
    expected = np.zeros((4, 3))
    for col, target in enumerate(mesh_ids):
        others = [mesh_id for mesh_id in mesh_ids if mesh_id != target]
        for source in others:
            expected[mesh_ids.index(source), col] = _force_x(beach, target, [source])
        expected[3, col] = _force_x(beach, target, others)
    np.testing.assert_allclose(matrix, expected, rtol=1.0e-12, atol=0.0)
    # Charges +1, +2, -3 nC at x = 0, 1, 2: like charges repel, unlike attract.
    assert matrix[0, 1] > 0.0 and matrix[1, 0] < 0.0
    assert matrix[0, 2] < 0.0 and matrix[2, 0] > 0.0
    np.testing.assert_allclose(np.diag(matrix[:3]), np.zeros(3))
    np.testing.assert_allclose(matrix[-1], np.sum(matrix[:-1], axis=0))
    np.testing.assert_array_equal(np.asarray(ax.images[0].get_array()), matrix)
    fig.clf()

    fig, ax = beach.plot_coulomb_force_matrix(component="x", workers=2)
    worker_info = getattr(ax, "_beach_coulomb_matrix")
    assert worker_info["workers"] == 2
    np.testing.assert_allclose(worker_info["matrix"], matrix_info["matrix"])
    fig.clf()


def test_plot_coulomb_force_matrix_accepts_explicit_target_kinds(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    out = tmp_path / "run_coulomb_matrix_plane"
    out.mkdir()
    _write_coulomb_matrix_fixture(out)
    beach = Beach(out)

    fig, ax = plot_coulomb_force_matrix(
        beach,
        component="x",
        target_kinds=("plane",),
        annotate=False,
    )

    matrix_info = getattr(ax, "_beach_coulomb_matrix")
    assert matrix_info["target_labels"] == ("plane",)
    assert matrix_info["source_labels"] == ("plane", "sphere1", "sphere2", "net")
    expected = [0.0, _force_x(beach, 1, [2]), _force_x(beach, 1, [3]), _force_x(beach, 1, [2, 3])]
    np.testing.assert_allclose(matrix_info["matrix"][:, 0], expected, rtol=1.0e-12, atol=0.0)
    # The +2 nC panel pushes the +1 nC plane toward -x; the -3 nC panel pulls it to +x.
    assert expected[1] < 0.0 < expected[2]
    fig.clf()


@pytest.mark.parametrize(
    ("snapshot_count", "frame_stride", "total_frames", "expected"),
    [
        pytest.param(10, 3, None, [0, 3, 6, 9], id="stride"),
        pytest.param(10, 1, 4, [0, 3, 6, 9], id="even-total"),
        pytest.param(3, 1, 10, [0, 1, 2], id="total-exceeds-snapshots"),
        pytest.param(7, 1, 1, [0], id="one-frame"),
    ],
)
def test_select_frame_columns(
    snapshot_count: int,
    frame_stride: int,
    total_frames: int | None,
    expected: list[int],
) -> None:
    cols = _select_frame_columns(
        snapshot_count,
        frame_stride=frame_stride,
        total_frames=total_frames,
    )

    np.testing.assert_array_equal(cols, np.array(expected, dtype=np.int64))


def test_animate_history_mesh_requires_charge_history(tmp_path: Path) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
        history=RunHistory(charge=None),
    )

    with pytest.raises(ValueError, match="charge_history.csv"):
        animate_history_mesh(result, tmp_path / "history.gif")


def test_animate_history_mesh_rejects_invalid_quantity(tmp_path: Path) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
        history=RunHistory(
            charge=_make_history(
                mesh_nelem=1,
                history=np.array([[1.0e-9, 1.2e-9]]),
            )
        ),
    )

    with pytest.raises(ValueError, match="quantity"):
        animate_history_mesh(result, tmp_path / "history.gif", quantity="invalid")


@pytest.mark.parametrize("source_model", ["point", "unknown"])
def test_potential_animation_rejects_legacy_or_missing_source_receipt(
    tmp_path: Path,
    source_model: str,
) -> None:
    pytest.importorskip("matplotlib")
    triangles = np.array(
        [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]],
        dtype=float,
    )
    result = RunResult(
        directory=tmp_path,
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
        history=RunHistory(
            charge=_make_history(
                mesh_nelem=1,
                history=np.array([[1.0e-9, 1.2e-9]]),
            )
        ),
        field=RunField(source_model=source_model),
    )

    with pytest.raises(ValueError, match="field_source_model='triangle_p0'"):
        animate_history_mesh(result, quantity="potential")


def test_animate_history_mesh_rejects_invalid_total_frames(tmp_path: Path) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
        history=RunHistory(
            charge=_make_history(
                mesh_nelem=1,
                history=np.array([[1.0e-9, 1.2e-9]]),
            )
        ),
    )

    with pytest.raises(ValueError, match="total_frames"):
        animate_history_mesh(result, tmp_path / "history.gif", total_frames=0)


def test_animate_history_mesh_rejects_stride_and_total_frames_combination(
    tmp_path: Path,
) -> None:
    triangles = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    result = RunResult(
        directory=Path("dummy"),
        charges=np.array([1.0e-9]),
        summary=RunSummary(
            mesh_nelem=1,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
        history=RunHistory(
            charge=_make_history(
                mesh_nelem=1,
                history=np.array([[1.0e-9, 1.2e-9]]),
            )
        ),
    )

    with pytest.raises(ValueError, match="cannot be used together"):
        animate_history_mesh(
            result,
            tmp_path / "history.gif",
            frame_stride=2,
            total_frames=2,
        )


def _history_result(directory: Path) -> RunResult:
    """Two unit-area panels whose three history frames all differ."""

    triangles = np.array(
        [
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [2.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        ]
    )
    charge_history = np.array(
        [
            [1.0e-10, 2.0e-10, 3.0e-10],
            [-2.0e-10, -1.0e-10, -0.5e-10],
        ]
    )
    return RunResult(
        directory=directory,
        charges=charge_history[:, -1],
        summary=RunSummary(
            mesh_nelem=2,
            processed_particles=0,
            absorbed=0,
            escaped=0,
            batches=0,
            escaped_boundary=0,
            survived_max_step=0,
            last_rel_change=0.0,
        ),
        mesh=RunMesh(triangles=triangles),
        history=RunHistory(
            charge=_make_history(
                mesh_nelem=2,
                history=charge_history,
                batch_indices=np.array([1, 2, 3]),
                processed_particles_by_batch=np.array([10, 20, 30]),
                rel_change_by_batch=np.array([1.0e-1, 1.0e-2, 1.0e-3]),
            )
        ),
    )


def _sorted_rows(values: np.ndarray) -> np.ndarray:
    return values[np.lexsort(values.T[::-1])]


def test_animate_history_mesh_writes_gif(tmp_path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    image_module = pytest.importorskip("PIL.Image")
    matplotlib.use("Agg")
    out = tmp_path / "charge.gif"

    written = animate_history_mesh(
        _history_result(tmp_path),
        out,
        quantity="charge",
        fps=4,
        frame_stride=1,
    )

    assert written == out
    frames = []
    with image_module.open(out) as image:
        assert image.format == "GIF"
        assert image.n_frames == 3
        for index in range(image.n_frames):
            image.seek(index)
            assert image.info["duration"] == 250
            frames.append(np.asarray(image.convert("RGB")))
    assert all(
        not np.array_equal(frames[first], frames[second])
        for first in range(3)
        for second in range(first + 1, 3)
    )


@pytest.mark.parametrize("quantity", ["charge", "potential"])
def test_animate_history_mesh_frames_follow_history_values(
    tmp_path: Path, quantity: str
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import AbstractMovieWriter
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    _write_complete_free_field_config(tmp_path)
    result = _history_result(tmp_path)
    history = result.history.charge.as_array()
    if quantity == "charge":
        expected = history  # unit-area panels: density equals charge
    else:
        expected = np.column_stack(
            [
                compute_potential_mesh(replace(result, charges=charges), source="kernel")
                for charges in history.T
            ]
        )
    limit = float(np.max(np.abs(expected)))
    cmap = matplotlib.colormaps["viridis"]
    norm = matplotlib.colors.Normalize(vmin=-limit, vmax=limit)

    class _FaceColorRecorder(AbstractMovieWriter):
        """Record each saved frame's face colors instead of encoding images."""

        def __init__(self) -> None:
            super().__init__()
            self.frames: list[np.ndarray] = []

        def setup(self, fig, outfile, dpi=None) -> None:
            super().setup(fig, outfile, dpi)

        def grab_frame(self, **_savefig_kwargs) -> None:
            self.fig.canvas.draw()
            (mesh,) = [
                item
                for item in self.fig.axes[0].collections
                if isinstance(item, Poly3DCollection)
            ]
            self.frames.append(np.asarray(mesh.get_facecolor())[:, :3])

        def finish(self) -> None:
            pass

    _, mesh_ax = plot_charge_mesh(result)
    animation = animate_history_mesh(result, quantity=quantity, cmap="viridis")
    # An omitted view uses the same defaults as the static mesh plot.
    anim_ax = plt.gcf().axes[0]
    assert (anim_ax.elev, anim_ax.azim) == (mesh_ax.elev, mesh_ax.azim)
    recorder = _FaceColorRecorder()
    animation.save(tmp_path / "unused.gif", writer=recorder)

    assert len(recorder.frames) == history.shape[1]
    for frame, colors in enumerate(recorder.frames):
        # Depth sorting may reorder faces, so compare the face colors as sorted rows.
        np.testing.assert_allclose(
            _sorted_rows(colors), _sorted_rows(cmap(norm(expected[:, frame]))[:, :3])
        )
