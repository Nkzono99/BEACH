from pathlib import Path

import pytest

from beach.cli.main import main as beachx_main


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        pytest.param(
            ["inspect", "no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="inspect",
        ),
        pytest.param(
            ["animate", "no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="animate",
        ),
        pytest.param(
            ["slices", "no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="slices",
        ),
        pytest.param(
            ["coulomb", "no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="coulomb",
        ),
        pytest.param(
            ["profile", "no_such_dir"],
            r'Performance profile file is missing under "no_such_dir/performance_profile\.csv"',
            id="profile",
        ),
        pytest.param(
            ["workload", "no_such_file.toml"],
            r"config file not found: no_such_file\.toml",
            id="workload",
        ),
    ],
)
def test_cli_missing_input_exits_with_friendly_message(
    argv: list[str],
    message: str,
) -> None:
    with pytest.raises(SystemExit, match=message):
        beachx_main(argv)


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        pytest.param(
            ["coulomb", "--softening", "nan"],
            "unrecognized arguments: --softening",
            id="coulomb-removed-softening",
        ),
        pytest.param(
            ["kernel-forces", "--theta", "nan"],
            "argument --theta: must be finite",
            id="kernel-theta",
        ),
        pytest.param(
            ["inspect", "--view-elev", "nan"],
            "argument --view-elev: must be finite",
            id="inspect-view",
        ),
        pytest.param(
            ["animate", "--periodic2-repeat", "-1"],
            "argument --periodic2-repeat: must be >= 0",
            id="animate-repeat",
        ),
        pytest.param(
            ["slices", "--vmin", "nan"],
            "argument --vmin: must be finite",
            id="slices-vmin",
        ),
    ],
)
def test_cli_rejects_invalid_option(
    capsys: pytest.CaptureFixture[str],
    argv: list[str],
    message: str,
) -> None:
    with pytest.raises(SystemExit) as exc:
        beachx_main(argv)

    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_cli_reports_library_value_error_as_message(tmp_path: Path) -> None:
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    (output_dir / "summary.txt").write_text(
        "mesh_nelem=1\nprocessed_particles=0\nabsorbed=0\nescaped=0\n"
        "batches=1\nlast_rel_change=0.0\n",
        encoding="utf-8",
    )
    (output_dir / "charges.csv").write_text(
        "elem_idx,charge_C\n1,0.0\n", encoding="utf-8"
    )

    with pytest.raises(SystemExit) as exc:
        beachx_main(["animate", str(output_dir), "--fps", "0"])

    assert exc.value.code == "fps must be > 0."
