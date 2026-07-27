from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from beach.cli import inspect_fortran_output as inspect_cli


def _write_inspect_fixture(
    output_dir: Path,
    *,
    include_triangles: bool = True,
    potential_values: tuple[float, float] | None = None,
) -> None:
    output_dir.mkdir()
    (output_dir / "summary.txt").write_text(
        "\n".join(
            [
                "mesh_nelem=2",
                "processed_particles=10",
                "absorbed=7",
                "escaped=3",
                "batches=1",
                "last_rel_change=1.0e-8",
                "field_source_model=triangle_p0",
            ]
        ),
        encoding="utf-8",
    )
    (output_dir / "charges.csv").write_text(
        "elem_idx,charge_C\n1,1.0e-9\n2,-2.0e-9\n",
        encoding="utf-8",
    )
    if include_triangles:
        (output_dir / "mesh_triangles.csv").write_text(
            "elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id\n"
            "1,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,1.0e-9,1\n"
            "2,1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,-2.0e-9,1\n",
            encoding="utf-8",
        )
    if potential_values is not None:
        (output_dir / "mesh_potential.csv").write_text(
            "elem_idx,potential_V\n"
            f"1,{potential_values[0]}\n"
            f"2,{potential_values[1]}\n",
            encoding="utf-8",
        )


def _fail_if_potential_is_computed(*_args: object, **_kwargs: object) -> np.ndarray:
    raise AssertionError("ordinary inspect must not call Beach.compute_potential")


def test_inspect_without_precomputed_potential_does_not_compute(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output_dir = tmp_path / "run_no_potential"
    _write_inspect_fixture(output_dir)
    monkeypatch.setattr(
        inspect_cli.Beach,
        "compute_potential",
        _fail_if_potential_is_computed,
    )

    inspect_cli.main([str(output_dir)])

    output = capsys.readouterr().out
    assert "mesh_nelem=2" in output
    assert "potential_min=" not in output
    assert "potential_max=" not in output


def test_inspect_reports_precomputed_potential_without_compute(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output_dir = tmp_path / "run_precomputed"
    _write_inspect_fixture(
        output_dir,
        include_triangles=False,
        potential_values=(1.5, -2.5),
    )
    monkeypatch.setattr(
        inspect_cli.Beach,
        "compute_potential",
        _fail_if_potential_is_computed,
    )

    inspect_cli.main([str(output_dir)])

    output = capsys.readouterr().out
    assert "potential_min=-2.500000e+00" in output
    assert "potential_max=1.500000e+00" in output


@pytest.mark.parametrize("with_precomputed", [False, True])
def test_inspect_recompute_potential_calls_compute_and_uses_its_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    with_precomputed: bool,
) -> None:
    output_dir = tmp_path / f"run_recompute_{with_precomputed}"
    _write_inspect_fixture(
        output_dir,
        potential_values=(100.0, 200.0) if with_precomputed else None,
    )
    calls: list[dict[str, object]] = []

    def fake_compute(_self: object, **kwargs: object) -> np.ndarray:
        calls.append(kwargs)
        return np.array([7.0, -8.0])

    monkeypatch.setattr(inspect_cli.Beach, "compute_potential", fake_compute)

    inspect_cli.main([str(output_dir), "--recompute-potential"])

    output = capsys.readouterr().out
    assert calls == [
        {
            "reference_point": "species1_injection_center",
            "library_path": None,
        }
    ]
    assert "potential_min=-8.000000e+00" in output
    assert "potential_max=7.000000e+00" in output
    assert "potential_min=1.000000e+02" not in output
    assert "potential_max=2.000000e+02" not in output


def test_inspect_potential_plot_is_independent_of_recompute_flag(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output_dir = tmp_path / "run_plot"
    save_path = tmp_path / "potential.png"
    _write_inspect_fixture(output_dir)
    monkeypatch.setattr(
        inspect_cli.Beach,
        "compute_potential",
        _fail_if_potential_is_computed,
    )
    plot_calls: list[dict[str, object]] = []

    class FakeFigure:
        def savefig(self, path: Path, *, dpi: int) -> None:
            assert path == save_path
            assert dpi == 150

    def fake_plot(_self: object, **kwargs: object) -> tuple[FakeFigure, object]:
        plot_calls.append(kwargs)
        return FakeFigure(), object()

    monkeypatch.setattr(inspect_cli.Beach, "plot_potential", fake_plot)

    inspect_cli.main([str(output_dir), "--save-potential-mesh", str(save_path)])

    output = capsys.readouterr().out
    assert len(plot_calls) == 1
    assert "potential_min=" not in output
    assert "potential_max=" not in output
    assert f"saved_potential_mesh={save_path}" in output


def test_inspect_parser_documents_recompute_potential() -> None:
    parser = inspect_cli.build_parser()

    assert parser.parse_args([]).recompute_potential is False
    assert parser.parse_args(["--recompute-potential"]).recompute_potential is True
    help_text = parser.format_help()
    assert "--recompute-potential" in help_text
    assert "may be expensive" in help_text
    assert "potential plots" in help_text
