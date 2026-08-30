import pytest

from beach.cli_animate_fortran_history import main as animate_main
from beach.cli_estimate_fortran_workload import main as workload_main
from beach.cli_inspect_fortran_output import main as inspect_main
from beach.cli_plot_coulomb_force_matrix import main as plot_coulomb_main
from beach.cli_plot_performance_profile import main as plot_performance_main
from beach.cli.analyze_coulomb_mobility import main as mobility_main
from beach.cli.kernel_forces import main as kernel_forces_main
from beach.cli_plot_fortran_potential_slices import main as plot_slices_main


@pytest.mark.parametrize(
    ("cli_main", "argv", "message"),
    [
        pytest.param(
            inspect_main,
            ["no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="inspect",
        ),
        pytest.param(
            animate_main,
            ["no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="animate",
        ),
        pytest.param(
            plot_slices_main,
            ["no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="slices",
        ),
        pytest.param(
            plot_coulomb_main,
            ["no_such_dir"],
            r'Fortran output files are missing under "no_such_dir"',
            id="coulomb",
        ),
        pytest.param(
            plot_performance_main,
            ["no_such_dir"],
            r'Performance profile file is missing under "no_such_dir/performance_profile\.csv"',
            id="profile",
        ),
        pytest.param(
            workload_main,
            ["no_such_file.toml"],
            r"config file not found: no_such_file\.toml",
            id="workload",
        ),
    ],
)
def test_cli_missing_input_exits_with_friendly_message(
    cli_main,
    argv: list[str],
    message: str,
) -> None:
    with pytest.raises(SystemExit, match=message):
        cli_main(argv)


@pytest.mark.parametrize(
    ("cli_main", "argv", "message"),
    [
        pytest.param(
            plot_coulomb_main,
            ["--softening", "nan"],
            "unrecognized arguments: --softening",
            id="coulomb-removed-softening",
        ),
        pytest.param(
            kernel_forces_main,
            ["--theta", "nan"],
            "--theta must be finite",
            id="kernel-theta",
        ),
        pytest.param(
            mobility_main,
            ["--density-kg-m3", "inf"],
            "--density-kg-m3 must be finite",
            id="mobility-density",
        ),
        pytest.param(
            inspect_main,
            ["--view-elev", "nan"],
            "--view-elev must be finite",
            id="inspect-view",
        ),
        pytest.param(
            animate_main,
            ["--fps", "0"],
            "--fps must be > 0",
            id="animate-fps",
        ),
        pytest.param(
            plot_slices_main,
            ["--vmin", "nan"],
            "--vmin must be finite",
            id="slices-vmin",
        ),
    ],
)
def test_cli_rejects_invalid_option(
    capsys: pytest.CaptureFixture[str],
    cli_main,
    argv: list[str],
    message: str,
) -> None:
    with pytest.raises(SystemExit) as exc:
        cli_main(argv)

    assert exc.value.code == 2
    assert message in capsys.readouterr().err
