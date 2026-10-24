from pathlib import Path

import pytest

from beach.cli.main import main as beachx_main


MISSING_OUTPUT = (
    'Fortran output files are missing under "no_such_dir". '
    "Expected at least summary.txt and charges.csv."
)


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        *(
            pytest.param([command, "no_such_dir"], MISSING_OUTPUT, id=command)
            for command in (
                "inspect", "animate", "slices", "coulomb", "kernel-forces", "mobility",
            )
        ),
        pytest.param(
            ["profile", "no_such_dir"],
            'Performance profile file is missing under "no_such_dir/performance_profile.csv"',
            id="profile",
        ),
        pytest.param(
            ["workload", "no_such_file.toml"],
            "config file not found: no_such_file.toml",
            id="workload",
        ),
        pytest.param(
            ["model", "close-pack", "--layers", "4", "--radius", "0.2", "--cells-x", "2"],
            "base config file not found: beach.toml",
            id="close-pack",
        ),
    ],
)
def test_cli_missing_input_exits_with_friendly_message(
    argv: list[str],
    message: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        beachx_main(argv)

    assert str(excinfo.value).startswith(message)


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
