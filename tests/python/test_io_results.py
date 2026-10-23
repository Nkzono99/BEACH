"""Result io contracts."""

from pathlib import Path
import numpy as np
import pytest
from beach import (
    MatchingPlaneHistoryEntry as PublicMatchingPlaneHistoryEntry,
    MatchingPlaneState as PublicMatchingPlaneState,
)
from beach import (
    Beach,
    FortranChargeHistory,
    MatchingPlaneHistoryEntry,
    MatchingPlaneState,
    list_fortran_runs,
    load_fortran_result,
)
from result_fixtures import (
    _charge_ledger_row,
    _matching_history_row,
    _matching_plane_summary_extra,
    _write_charge_history,
    _write_minimal_result_fixture,
    _write_three_mesh_fixture,
    _unit_panel_kernel as _unit_panel_kernel,
)


def test_load_fortran_result(tmp_path: Path) -> None:
    out = tmp_path / "run1"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=10",
                "absorbed=7",
                "escaped=3",
                "batches=1",
                "escaped_boundary=2",
                "survived_max_step=1",
                "multiple_box_events_retry_attempted=6",
                "multiple_box_events_retry_resolved=2",
                "multiple_box_events_soft_discarded=4",
                "multiple_box_events_soft_discarded_abs_charge_C=2.5e-15",
                "multiple_box_events_soft_discard_fraction=4.0e-1",
                "last_rel_change=1.0e-8",
                "simulated_time_s=7.5e-6",
                "adaptive_nonzero_mode_rejected_trials=3",
                "adaptive_nonzero_mode_last_batch_duration_s=2.5e-6",
                "adaptive_nonzero_mode_last_potential_step_V=8.0e-3",
                "adaptive_nonzero_mode_omp_threads=16",
                "periodic2_max_nonzero_mode_potential_step_V=1.0e-2",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-10\n2,-2.0e-10\n", encoding="utf-8"
    )
    (out / "mesh_triangles.csv").write_text(
        "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
        "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
        "2,0,0,1,1,0,1,0,1,1,-2.0e-10,2\n",
        encoding="utf-8",
    )
    (out / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,plane,insulator,1.0,1\n"
        "2,template,sphere,conductor,2.5,1\n",
        encoding="utf-8",
    )
    (out / "mesh_potential.csv").write_text(
        "elem_idx,potential_V\n1,1.5\n2,-2.5\n",
        encoding="utf-8",
    )

    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n"
        "1,5,3.0e-1,1,2.0e-11\n"
        "1,5,3.0e-1,2,-1.0e-11\n"
        "3,10,1.0e-8,1,1.0e-10\n"
        "3,10,1.0e-8,2,-2.0e-10\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.mesh_nelem == 2
    assert result.absorbed == 7
    assert result.escaped_boundary == 2
    assert result.survived_max_step == 1
    assert result.multiple_box_events_retry_attempted == 6
    assert result.multiple_box_events_retry_resolved == 2
    assert result.multiple_box_events_soft_discarded == 4
    assert result.multiple_box_events_soft_discarded_abs_charge_c == 2.5e-15
    assert result.multiple_box_events_soft_discard_fraction == 0.4
    assert result.simulated_time_s == 7.5e-6
    assert result.adaptive_nonzero_mode_rejected_trials == 3
    assert result.adaptive_nonzero_mode_last_batch_duration_s == 2.5e-6
    assert result.adaptive_nonzero_mode_last_potential_step_v == 8.0e-3
    assert result.adaptive_nonzero_mode_omp_threads == 16
    assert result.periodic2_max_nonzero_mode_potential_step_v == 1.0e-2
    assert result.triangles is not None
    assert result.triangles.shape == (2, 3, 3)
    assert result.history is not None
    assert result.mesh_ids is not None
    np.testing.assert_array_equal(result.mesh_ids, np.array([1, 2]))
    assert result.mesh_sources is not None
    assert result.mesh_sources[1].template_kind == "plane"
    assert result.mesh_sources[1].surface_model == "insulator"
    assert result.mesh_sources[1].epsilon_r == 1.0
    assert result.mesh_sources[2].surface_model == "conductor"
    assert result.mesh_sources[2].epsilon_r == 2.5
    np.testing.assert_allclose(result.charges, np.array([1.0e-10, -2.0e-10]))
    assert result.mesh_potential_v is not None
    np.testing.assert_allclose(result.mesh_potential_v, np.array([1.5, -2.5]))
    np.testing.assert_allclose(
        result.history.as_array(),
        np.array([[2.0e-11, 1.0e-10], [-1.0e-11, -2.0e-10]]),
    )
    np.testing.assert_array_equal(
        result.history.processed_particles_by_batch, np.array([5, 10])
    )
    np.testing.assert_allclose(
        result.history.rel_change_by_batch, np.array([3.0e-1, 1.0e-8])
    )
    np.testing.assert_array_equal(result.history.batch_indices, np.array([1, 3]))


@pytest.mark.parametrize("field_values", [None, "-2.7,-2.7", "NaN,NaN"])
def test_load_fortran_result_matching_plane_state_and_history(tmp_path: Path, field_values: str | None) -> None:
    out = tmp_path / "matching"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "checkpoint_schema_version=9",
            *_matching_plane_summary_extra(),
        ],
    )
    (out / "matching_plane_history.csv").write_text(
        (f"{_MATCHING_HISTORY_HEADER},secondary_share,secondary_budget_residual_a_m2,E_top_sheath_V_m,E_top_cell_mean_V_m\n"
         f"{_matching_history_row()},0,0,{field_values}\n") if field_values is not None else
        f"{_MATCHING_HISTORY_HEADER}\n{_matching_history_row()}\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.matching_plane_state is not None
    assert PublicMatchingPlaneState is MatchingPlaneState
    assert isinstance(result.matching_plane_state, PublicMatchingPlaneState)
    assert result.matching_plane_state.phi_v == 1.5
    assert result.matching_plane_state.photoelectron_escape_flux_m2_s == 1.5
    assert result.matching_plane_history is not None
    assert len(result.matching_plane_history) == 1
    assert PublicMatchingPlaneHistoryEntry is MatchingPlaneHistoryEntry
    assert isinstance(result.matching_plane_history[0], PublicMatchingPlaneHistoryEntry)
    assert result.matching_plane_history[0].batch == 1
    assert result.matching_plane_history[0].state.iterations == 3

    expected_field = -2.7 if field_values == "-2.7,-2.7" else None
    assert result.matching_plane_history[0].state.top_sheath_field_v_m == expected_field
    assert result.matching_plane_history[0].state.top_cell_mean_field_v_m == expected_field


def test_load_fortran_result_ignores_stale_matching_history(tmp_path: Path) -> None:
    out = tmp_path / "nonmatching_reused_output"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=["matching_plane_state_valid=F"],
    )
    (out / "matching_plane_history.csv").write_text(
        "stale,history,from,a,previous,run\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.matching_plane_state is None
    assert result.matching_plane_history is None


def test_load_fortran_result_rejects_inconsistent_matching_budget(
    tmp_path: Path,
) -> None:
    out = tmp_path / "matching_bad_budget"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=_matching_plane_summary_extra(
            photoelectron_escape_flux_m2_s=0.5
        ),
    )

    with pytest.raises(ValueError, match="photoelectron budget"):
        load_fortran_result(out)


@pytest.mark.parametrize(
    ("history_header", "history_rows", "match"),
    [
        pytest.param(
            "batch,simulated_time_s",
            (_matching_history_row(),),
            "header does not match",
            id="wrong-header",
        ),
        pytest.param(
            _MATCHING_HISTORY_HEADER,
            (_matching_history_row().rsplit(",", maxsplit=1)[0],),
            "exactly 17 columns",
            id="short-row",
        ),
        pytest.param(
            _MATCHING_HISTORY_HEADER,
            (_matching_history_row(electron_inward_flux="nan"),),
            "electron_inward_flux_m2_s must be finite",
            id="nonfinite-value",
        ),
        pytest.param(
            _MATCHING_HISTORY_HEADER,
            (_matching_history_row(), _matching_history_row()),
            "positive and strictly increasing",
            id="nonincreasing-batch",
        ),
        pytest.param(
            _MATCHING_HISTORY_HEADER,
            (_matching_history_row(photoelectron_escape_flux="0.5"),),
            "photoelectron budget",
            id="inconsistent-photoelectron-budget",
        ),
    ],
)
def test_load_fortran_result_rejects_corrupt_matching_plane_history(
    tmp_path: Path,
    history_header: str,
    history_rows: tuple[str, ...],
    match: str,
) -> None:
    out = tmp_path / "matching_bad_history"
    out.mkdir()
    _write_minimal_result_fixture(out, summary_extra=_matching_plane_summary_extra())
    (out / "matching_plane_history.csv").write_text(
        history_header + "\n" + "\n".join(history_rows) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


@pytest.mark.parametrize(
    ("charges_text", "match"),
    [
        pytest.param(
            "elem_idx,charge_C\n1,1.0e-10\n",
            "charges.csv row count",
            id="row-count",
        ),
        pytest.param(
            "elem_idx,charge_C\n1,1.0e-10\n3,-2.0e-10\n",
            "charges.csv elem_idx",
            id="elem-index",
        ),
        pytest.param(
            "elem_idx,charge_C\n1,1.0e-10\n2,nan\n",
            "charge_C",
            id="nonfinite-charge",
        ),
    ],
)
def test_load_fortran_result_rejects_corrupt_charges(
    tmp_path: Path,
    charges_text: str,
    match: str,
) -> None:
    out = tmp_path / "run_bad_charges"
    out.mkdir()
    _write_minimal_result_fixture(out, charges_text=charges_text)

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


def test_load_fortran_result_rejects_mesh_triangle_row_count_mismatch(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_bad_triangles_count"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_triangles_text=(
            "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
            "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
        ),
    )

    with pytest.raises(ValueError, match="mesh_triangles.csv row count"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_nonfinite_triangle_vertex(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_triangle_nan"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_triangles_text=(
            "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
            "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
            "2,0,0,1,1,nan,1,0,1,1,-2.0e-10,2\n"
        ),
    )

    with pytest.raises(ValueError, match="vertex values"):
        load_fortran_result(out)


@pytest.mark.parametrize(
    ("mesh_id", "match"),
    [
        pytest.param("1.5", "mesh_id values must be integer", id="fractional"),
        pytest.param("0", "mesh_id values must be positive", id="nonpositive"),
    ],
)
def test_load_fortran_result_rejects_invalid_triangle_mesh_id(
    tmp_path: Path,
    mesh_id: str,
    match: str,
) -> None:
    out = tmp_path / "run_bad_triangle_mesh_id"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_triangles_text=(
            "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
            "1,0,0,0,1,0,0,0,1,0,1.0e-10,1\n"
            f"2,0,0,1,1,0,1,0,1,1,-2.0e-10,{mesh_id}\n"
        ),
    )

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


@pytest.mark.parametrize(
    ("second_source", "match"),
    [
        pytest.param(
            "2,template,sphere,dielectric,nan,1",
            "epsilon_r",
            id="nonfinite-epsilon",
        ),
        pytest.param(
            "2,template,sphere,resistive,2.0,1",
            "surface_model",
            id="unknown-surface-model",
        ),
        pytest.param(
            "2,template,sphere,dielectric,2.0,-1",
            "elem_count",
            id="negative-element-count",
        ),
        pytest.param(
            "0,template,sphere,dielectric,2.0,1",
            "mesh_id",
            id="nonpositive-mesh-id",
        ),
    ],
)
def test_load_fortran_result_rejects_invalid_mesh_source_fields(
    tmp_path: Path,
    second_source: str,
    match: str,
) -> None:
    out = tmp_path / "run_bad_mesh_source"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_sources_text=(
            "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
            "1,template,plane,insulator,1.0,1\n"
            f"{second_source}\n"
        ),
    )

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


@pytest.mark.parametrize(
    ("mesh_potential_text", "match"),
    [
        pytest.param(
            "elem_idx,potential_V\n1,1.5\n",
            "mesh_potential.csv row count",
            id="row-count",
        ),
        pytest.param(
            "elem_idx,potential_V\n1,1.5\n3,-2.5\n",
            "mesh_potential.csv elem_idx",
            id="elem-index",
        ),
        pytest.param(
            "elem_idx,potential_V\n1,1.5\n2,inf\n",
            "potential_V",
            id="nonfinite-potential",
        ),
    ],
)
def test_load_fortran_result_rejects_corrupt_mesh_potential(
    tmp_path: Path,
    mesh_potential_text: str,
    match: str,
) -> None:
    out = tmp_path / "run_bad_mesh_potential"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        mesh_potential_text=mesh_potential_text,
    )

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


def test_load_fortran_result_rejects_negative_summary_stat(tmp_path: Path) -> None:
    out = tmp_path / "run_bad_summary"
    out.mkdir()
    _write_minimal_result_fixture(out, summary_extra=["escaped_boundary=-1"])

    with pytest.raises(ValueError, match="escaped_boundary"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_inconsistent_retry_counts(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_bad_retry_counts"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "multiple_box_events_retry_attempted=1",
            "multiple_box_events_retry_resolved=2",
        ],
    )

    with pytest.raises(ValueError, match="retry_resolved must not exceed"):
        load_fortran_result(out)


def test_load_fortran_result_derives_missing_soft_discard_fraction(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_legacy_soft_discard_fraction"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=["multiple_box_events_soft_discarded=4"],
    )

    result = load_fortran_result(out)

    assert result.multiple_box_events_soft_discard_fraction == 0.4


def test_load_fortran_result_accepts_zero_soft_discard_fraction(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_zero_soft_discard_fraction"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        processed_particles=0,
        absorbed=0,
        escaped=0,
        summary_extra=[
            "multiple_box_events_soft_discarded=0",
            "multiple_box_events_soft_discard_fraction=0.0",
        ],
    )

    result = load_fortran_result(out)

    assert result.multiple_box_events_soft_discard_fraction == 0.0


def test_load_fortran_result_accepts_soft_discard_fraction_roundoff(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_rounded_soft_discard_fraction"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "multiple_box_events_soft_discarded=1",
            "multiple_box_events_soft_discard_fraction=1.000000000000001e-1",
        ],
    )

    result = load_fortran_result(out)

    assert result.multiple_box_events_soft_discard_fraction == 0.1


def test_load_fortran_result_rejects_soft_discard_count_above_processed(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_bad_soft_discard_count"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=["multiple_box_events_soft_discarded=11"],
    )

    with pytest.raises(ValueError, match="must not exceed processed_particles"):
        load_fortran_result(out)


def test_load_fortran_result_rejects_inconsistent_soft_discard_fraction(
    tmp_path: Path,
) -> None:
    out = tmp_path / "run_inconsistent_soft_discard_fraction"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "multiple_box_events_soft_discarded=4",
            "multiple_box_events_soft_discard_fraction=0.3",
        ],
    )

    with pytest.raises(ValueError, match="must match"):
        load_fortran_result(out)


@pytest.mark.parametrize("value", ["-1.0e-6", "1.000001", "nan", "inf"])
def test_load_fortran_result_rejects_invalid_soft_discard_fraction(
    tmp_path: Path,
    value: str,
) -> None:
    out = tmp_path / "run_bad_soft_discard_fraction"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[f"multiple_box_events_soft_discard_fraction={value}"],
    )

    with pytest.raises(ValueError, match="multiple_box_events_soft_discard_fraction"):
        load_fortran_result(out)


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*missing elem_idx=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,2,2.0e-9",
                "1,2,1.0e-1,2,9.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*duplicate elem_idx=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,2,2.0e-9",
                "1,2,1.0e-1,4,4.0e-9",
            ],
            r"batch=1.*out-of-range elem_idx=4",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,1.0e-1,2,nan",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*elem_idx=2.*non-finite charge_C",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,3,1.0e-1,2,2.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*inconsistent processed_particles.*elem_idx=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,2.0e-1,2,2.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*inconsistent rel_change.*elem_idx=2",
        ),
        (
            [
                "2,2,1.0e-1,1,1.0e-9",
                "2,2,1.0e-1,2,2.0e-9",
                "2,2,1.0e-1,3,3.0e-9",
                "1,3,5.0e-2,1,4.0e-9",
                "1,3,5.0e-2,2,5.0e-9",
                "1,3,5.0e-2,3,6.0e-9",
            ],
            r"batch=1.*previous batch=2",
        ),
        (
            [
                "1,2,1.0e-1,1,1.0e-9",
                "1,2,nan,2,2.0e-9",
                "1,2,1.0e-1,3,3.0e-9",
            ],
            r"batch=1.*elem_idx=2.*non-finite rel_change",
        ),
    ],
    ids=[
        "missing-element",
        "duplicate-element",
        "out-of-range-element",
        "non-finite-charge",
        "inconsistent-processed-particles",
        "inconsistent-rel-change",
        "decreasing-batch",
        "non-finite-rel-change",
    ],
)
def test_charge_history_rejects_corrupt_batch_when_indexed(
    tmp_path: Path,
    rows: list[str],
    message: str,
) -> None:
    history = _write_charge_history(tmp_path / "charge_history.csv", rows)

    for _ in range(2):
        with pytest.raises(ValueError, match=message):
            _ = history.batch_indices


def test_load_fortran_result_defers_corrupt_history_validation(tmp_path: Path) -> None:
    out = tmp_path / "run_corrupt_lazy_history"
    out.mkdir()
    _write_minimal_result_fixture(out)
    (out / "charge_history.csv").write_text(
        "batch,processed_particles,rel_change,elem_idx,charge_C\n1,2,1.0e-1,1,1.0e-9\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.history is not None
    with pytest.raises(ValueError, match=r"batch=1.*missing elem_idx=2"):
        _ = result.history.batch_indices


def test_charge_history_accepts_permuted_dense_rows(tmp_path: Path) -> None:
    history = _write_charge_history(
        tmp_path / "charge_history.csv",
        [
            "1,2,1.0e-1,3,3.0e-9",
            "1,2,1.0e-1,1,1.0e-9",
            "1,2,1.0e-1,2,2.0e-9",
        ],
    )

    np.testing.assert_array_equal(history.batch_indices, np.array([1]))
    np.testing.assert_allclose(history.get_step(1), np.array([1.0e-9, 2.0e-9, 3.0e-9]))


def test_load_fortran_result_history_supports_step_access(tmp_path: Path) -> None:
    out = tmp_path / "run_lazy"
    out.mkdir()
    _write_three_mesh_fixture(out)

    result = load_fortran_result(out)

    assert result.history is not None
    np.testing.assert_allclose(result.history[1], np.array([5.0e-10, 1.0e-9, -1.5e-9]))
    np.testing.assert_allclose(
        result.history_at(10), np.array([1.0e-9, 2.0e-9, -3.0e-9])
    )
    history = result.history.as_array()
    np.testing.assert_allclose(
        history,
        np.array(
            [
                [5.0e-10, 1.0e-9],
                [1.0e-9, 2.0e-9],
                [-1.5e-9, -3.0e-9],
            ]
        ),
    )
    np.testing.assert_array_equal(result.history.batch_indices, np.array([1, 10]))


def test_charge_history_materialization_preserves_cached_edits(tmp_path: Path) -> None:
    _write_three_mesh_fixture(tmp_path)
    history = FortranChargeHistory(tmp_path / "charge_history.csv", mesh_nelem=3)
    first = history[1]
    first[0] = 7.0e-9
    with pytest.raises(ValueError, match="step=2 is not found"):
        history[2]

    values = history.as_array()
    np.testing.assert_array_equal(values[:, 0], first)
    np.testing.assert_array_equal(values[:, 1], [1.0e-9, 2.0e-9, -3.0e-9])
    np.testing.assert_array_equal(history[-1], values[:, -1])
    assert history.as_array() is values
    with pytest.raises(ValueError, match="step=2 is not found"):
        history[2]


def test_list_fortran_runs(tmp_path: Path) -> None:
    valid = tmp_path / "valid"
    valid.mkdir()
    (valid / "summary.txt").write_text(
        "mesh_nelem=1\nprocessed_particles=1\nabsorbed=1\nescaped=0\nbatches=1\nlast_rel_change=0.0\n",
        encoding="utf-8",
    )
    (valid / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    invalid = tmp_path / "invalid"
    invalid.mkdir()
    (invalid / "summary.txt").write_text("mesh_nelem=1\n", encoding="utf-8")

    runs = list_fortran_runs(tmp_path)
    assert runs == [valid]


def test_beach_uses_outputs_latest_by_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "outputs" / "latest"
    out.mkdir(parents=True)
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=1",
                "absorbed=1",
                "escaped=0",
                "batches=1",
                "last_rel_change=0.0",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    monkeypatch.chdir(tmp_path)

    beach = Beach()
    assert beach.output_dir == Path("outputs/latest")
    assert beach.result.mesh_nelem == 1


def test_beach_default_lazy_history_keeps_step_access(tmp_path: Path) -> None:
    out = tmp_path / "run_beach_lazy"
    out.mkdir()
    _write_three_mesh_fixture(out)

    beach = Beach(out)
    result = beach.result

    assert result.history is not None
    np.testing.assert_allclose(beach.get_mesh_charge(2, step=10), np.array([2.0e-9]))
    np.testing.assert_allclose(
        result.history.get_step(10), np.array([1.0e-9, 2.0e-9, -3.0e-9])
    )


def test_load_fortran_result_without_new_summary_keys(tmp_path: Path) -> None:
    out = tmp_path / "legacy"
    out.mkdir()
    (out / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=1",
                "processed_particles=4",
                "absorbed=2",
                "escaped=2",
                "batches=1",
                "last_rel_change=0.5",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    result = load_fortran_result(out)

    assert result.escaped_boundary == 0
    assert result.survived_max_step == 0
    assert result.checkpoint_schema_version is None
    assert result.charge_ledger is None
    assert result.simulated_time_s is None
    assert result.adaptive_nonzero_mode_rejected_trials == 0
    assert result.adaptive_nonzero_mode_last_batch_duration_s is None
    assert result.adaptive_nonzero_mode_last_potential_step_v is None
    assert result.adaptive_nonzero_mode_omp_threads is None
    assert result.periodic2_max_nonzero_mode_potential_step_v is None
    assert result.multiple_box_events_soft_discard_fraction == 0.0


@pytest.mark.parametrize(
    ("extra_line", "match"),
    [
        (
            None,
            "adaptive nonzero-mode output requires adaptive_nonzero_mode_omp_threads",
        ),
        (
            "adaptive_nonzero_mode_omp_threads=0",
            "adaptive_nonzero_mode_omp_threads must be > 0",
        ),
    ],
)
def test_load_fortran_result_rejects_incomplete_adaptive_summary(
    tmp_path: Path, extra_line: str | None, match: str
) -> None:
    out = tmp_path / "incomplete_adaptive"
    out.mkdir()
    lines = [
        "mesh_nelem=1",
        "processed_particles=4",
        "absorbed=2",
        "escaped=2",
        "batches=1",
        "last_rel_change=0.5",
        "field_source_model=triangle_p0",
        "periodic2_max_nonzero_mode_potential_step_V=1.0e-2",
        "simulated_time_s=1.0e-3",
        "adaptive_nonzero_mode_rejected_trials=1",
        "adaptive_nonzero_mode_last_batch_duration_s=5.0e-4",
        "adaptive_nonzero_mode_last_potential_step_V=8.0e-3",
    ]
    if extra_line is not None:
        lines.append(extra_line)
    (out / "summary.txt").write_text("\n".join(lines), encoding="utf-8")
    (out / "charges.csv").write_text("elem_idx,charge_C\n1,0.0\n", encoding="utf-8")

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


def test_load_fortran_result_model_contract_and_charge_ledger(tmp_path: Path) -> None:
    out = tmp_path / "contract"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "checkpoint_schema_version=4",
            "model_fingerprint=0123456789ABCDEF",
            "mesh_fingerprint=1111111122222222",
            "species_fingerprint=3333333344444444",
            "charge_ledger_residual_C=1.0e-18",
            "periodic2_cache_hit=T",
            "periodic2_operator_build_count=0",
            "periodic2_cache_fingerprint=ABCDEF0123456789",
            "periodic2_cache_path=.beach_cache/periodic2/operator.bin",
            "periodic2_generation_tolerance=1.0e-8",
        ],
    )
    (out / "charge_ledger.csv").write_text(
        "batch,species_idx,injected_from_remote_C,emitted_from_surface_C,"
        "absorbed_on_surface_C,escaped_to_infinity_C,discarded_unresolved_C,"
        "neutral_return_correction_C,neutral_return_weight_scale,"
        "neutral_return_unresolved_fraction,fixed_absorbed_target_charge_C,"
        "fixed_absorbed_weight_scale,fixed_emission_target_charge_C,"
        "fixed_emission_weight_scale,fixed_current_correction_C,"
        "fixed_absorbed_applied_charge_C,fixed_emission_applied_charge_C,"
        "fixed_escape_target_charge_C,fixed_escape_applied_charge_C,"
        "fixed_escape_correction_C,injected_count,"
        "emitted_count,absorbed_count,escaped_count,discarded_unresolved_count\n"
        "1,1,-3,0,-2,-1,-0.1,-0.1,1.05,0.05,-4,2,3,1.5,-2,-4,3,-1,-1,-0.5,3,0,2,1,1\n",
        encoding="utf-8",
    )

    result = load_fortran_result(out)

    assert result.checkpoint_schema_version == 4
    assert result.model_fingerprint == "0123456789ABCDEF"
    assert result.mesh_fingerprint == "1111111122222222"
    assert result.species_fingerprint == "3333333344444444"
    assert result.charge_ledger_residual_c == pytest.approx(1.0e-18)
    assert result.charge_ledger is not None
    assert result.charge_ledger[0].species_idx == 1
    assert result.charge_ledger[0].injected_from_remote_c == pytest.approx(-3.0)
    assert result.charge_ledger[0].escaped_count == 1
    assert result.charge_ledger[0].neutral_return_correction_c == pytest.approx(-0.1)
    assert result.charge_ledger[0].neutral_return_weight_scale == pytest.approx(1.05)
    assert result.charge_ledger[0].neutral_return_unresolved_fraction == pytest.approx(
        0.05
    )
    assert result.charge_ledger[0].fixed_absorbed_target_charge_c == pytest.approx(-4.0)
    assert result.charge_ledger[0].fixed_absorbed_weight_scale == pytest.approx(2.0)
    assert result.charge_ledger[0].fixed_emission_target_charge_c == pytest.approx(3.0)
    assert result.charge_ledger[0].fixed_emission_weight_scale == pytest.approx(1.5)
    assert result.charge_ledger[0].fixed_absorbed_applied_charge_c == pytest.approx(
        -4.0
    )
    assert result.charge_ledger[0].fixed_emission_applied_charge_c == pytest.approx(3.0)
    assert result.charge_ledger[0].fixed_escape_target_charge_c == pytest.approx(-1.0)
    assert result.charge_ledger[0].fixed_escape_applied_charge_c == pytest.approx(-1.0)
    assert result.charge_ledger[0].fixed_escape_correction_c == pytest.approx(-0.5)
    assert result.charge_ledger[0].fixed_current_correction_c == pytest.approx(-2.0)
    assert result.periodic2_cache_hit is True
    assert result.periodic2_operator_build_count == 0
    assert result.periodic2_cache_fingerprint == "ABCDEF0123456789"
    assert result.periodic2_cache_path == ".beach_cache/periodic2/operator.bin"
    assert result.periodic2_generation_tolerance == pytest.approx(1.0e-8)


@pytest.mark.parametrize(
    ("ledger_text", "match"),
    [
        pytest.param(
            "batch,species_idx,injected_from_remote_C\n1,1,-3.0\n",
            "missing required columns",
            id="missing-required-column",
        ),
        pytest.param(
            f"{_CHARGE_LEDGER_REQUIRED_HEADER}\n"
            f"{_charge_ledger_row(species_idx='0')}\n",
            "indices and counts are invalid",
            id="invalid-species-index",
        ),
        pytest.param(
            f"{_CHARGE_LEDGER_REQUIRED_HEADER}\n"
            f"{_charge_ledger_row(injected_charge='nan')}\n",
            "charge values must be finite",
            id="nonfinite-charge",
        ),
    ],
)
def test_load_fortran_result_rejects_corrupt_charge_ledger(
    tmp_path: Path,
    ledger_text: str,
    match: str,
) -> None:
    out = tmp_path / "bad_charge_ledger"
    out.mkdir()
    _write_minimal_result_fixture(out)
    (out / "charge_ledger.csv").write_text(ledger_text, encoding="utf-8")

    with pytest.raises(ValueError, match=match):
        load_fortran_result(out)


@pytest.mark.parametrize("switch_active", [False, True])
def test_load_fortran_result_field_reconstruction_receipt(tmp_path: Path, switch_active: bool) -> None:
    out = tmp_path / "field_receipt"
    out.mkdir()
    _write_minimal_result_fixture(
        out,
        summary_extra=[
            "field_reconstruction_schema_version=2",
            "field_reconstruction_resolved_field_solver=fmm",
            "field_reconstruction_fmm_expansion_order=4",
            "field_reconstruction_field_bc_mode=periodic2",
            "field_reconstruction_tree_theta=3.75e-1",
            "field_reconstruction_tree_leaf_max=23",
            "field_reconstruction_e0_V_m=1.0 2.0 3.0",
            "field_reconstruction_use_box=T",
            "field_reconstruction_box_min_m=-1.0 -2.0 -3.0",
            "field_reconstruction_box_max_m=4.0 5.0 6.0",
            "field_reconstruction_boundary_low=2 2 0",
            "field_reconstruction_boundary_high=2 2 0",
            "field_reconstruction_periodic_image_layers=2",
            "field_reconstruction_periodic_far_correction=" + ("cached_kneq0" if switch_active else "none"),
            "field_reconstruction_periodic_nonzero_mode_backend=" + ("cached_kneq0" if switch_active else "legacy_finite_images"),
            "field_reconstruction_periodic_zero_mode_policy=" + ("exclude_k0" if switch_active else "legacy_not_decomposed"),
            "field_reconstruction_periodic_lower_boundary_model=" + ("e_top_from_sheath" if switch_active else "legacy_implicit"),
            "field_reconstruction_periodic_reference_mode_layers=4",
            "field_reconstruction_periodic_panel_quadrature_order=12",
            "field_reconstruction_periodic_ewald_alpha=2.5e-1",
            "field_reconstruction_periodic_ewald_layers=5",
            "field_reconstruction_periodic_cache_dir=verified-cache",
            "field_reconstruction_periodic_generation_tolerance=2.0e-9",
            *(["field_reconstruction_upper_fourier_switch_active=T",
               "field_reconstruction_upper_fourier_switch_height_m=5.5",
               "field_reconstruction_upper_fourier_switch_mode_layers=7"] if switch_active else []),
        ],
    )

    result = load_fortran_result(out)

    assert result.field_reconstruction is not None
    assert result.field_reconstruction.schema_version == 2
    assert result.field_reconstruction.resolved_field_solver == "fmm"
    assert result.field_reconstruction.fmm_expansion_order == 4
    assert result.field_reconstruction.external_e0_v_m == (1.0, 2.0, 3.0)
    assert result.field_reconstruction.boundary_low == (2, 2, 0)
    assert result.field_reconstruction.periodic_cache_dir == "verified-cache"

    assert result.field_reconstruction.upper_fourier_switch_active == switch_active
    assert result.field_reconstruction.upper_fourier_switch_height_m == (5.5 if switch_active else 0.0)
    assert result.field_reconstruction.upper_fourier_switch_mode_layers == (7 if switch_active else 0)


def test_load_fortran_result_secondary_emission_and_combined_matching_source(tmp_path: Path) -> None:
    out = tmp_path / "secondary"
    out.mkdir()
    _write_minimal_result_fixture(out, summary_extra=[
        *_matching_plane_summary_extra(),
        "matching_plane_secondary_share=0.25",
        "secondary_budget_residual_a_m2=-1e-8",
    ])
    (out / "matching_plane_history.csv").write_text(
        f"{_MATCHING_HISTORY_HEADER},secondary_share,secondary_budget_residual_a_m2\n"
        f"{_matching_history_row()},0.25,-1e-8\n", encoding="utf-8",
    )
    (out / "secondary_emission_history.csv").write_text(
        "batch,incident_species,hits,emitted_count,emitted_charge_c,mean_yield,generation_capped_count\n"
        '1,"incident,electron",3,6,-9.6e-19,2,0\n'
        '1,"secondary",1,0,0,1.5,2\n', encoding="utf-8",
    )
    result = load_fortran_result(out)
    assert result.matching_plane_state.secondary_share == 0.25
    assert result.matching_plane_history[0].state.secondary_budget_residual_a_m2 == -1e-8
    assert result.secondary_emission_history[0].incident_species == "incident,electron"
    assert result.secondary_emission_history[0].emitted_count == 6
    assert result.secondary_emission_history[0].emitted_charge_c == -9.6e-19
    assert result.secondary_emission_history[1].generation_capped_count == 2


@pytest.mark.parametrize("charge", ["nan", "1e-19"])
def test_secondary_emission_reader_rejects_invalid_emitted_charge(tmp_path: Path, charge: str) -> None:
    _write_minimal_result_fixture(tmp_path)
    (tmp_path / "secondary_emission_history.csv").write_text(
        "batch,incident_species,hits,emitted_count,emitted_charge_c,mean_yield,generation_capped_count\n"
        f"1,electron,3,6,{charge},2,0\n", encoding="utf-8",
    )
    with pytest.raises(ValueError, match="emitted_charge_c must be finite and nonpositive"):
        load_fortran_result(tmp_path)


@pytest.mark.parametrize("archived", [False, True])
def test_saved_config_mesh_parts_and_archived_templates_remain_readable(tmp_path: Path, archived: bool) -> None:
    import tomli_w
    from beach.config import default_config
    from beach.config._layout import to_grouped_layout
    from beach.analysis.context import load_config_for_output

    config = to_grouped_layout(default_config())
    if archived:
        config["mesh"]["templates"] = config["mesh"].pop("parts")
        config["mesh"]["mode"] = "template"
    path = tmp_path / "beach.toml"
    path.write_text(tomli_w.dumps(config), encoding="utf-8")
    loaded = load_config_for_output(tmp_path, config_path=path)
    assert loaded["mesh"] == config["mesh"]
    assert loaded["sim"]["dt"] == default_config()["sim"]["dt"]
