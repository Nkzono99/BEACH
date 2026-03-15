import pytest

from beach.cli_animate_fortran_history import main as animate_main
from beach.cli_inspect_fortran_output import main as inspect_main
from beach.cli_plot_performance_profile import main as plot_performance_main
from beach.cli_plot_fortran_potential_slices import main as plot_slices_main


def test_inspect_missing_output_dir_exits_with_friendly_message() -> None:
    with pytest.raises(
        SystemExit,
        match=r'Fortran output files are missing under "no_such_dir"\.',
    ):
        inspect_main(["no_such_dir"])


def test_animate_missing_output_dir_exits_with_friendly_message() -> None:
    with pytest.raises(
        SystemExit,
        match=r'Fortran output files are missing under "no_such_dir"\.',
    ):
        animate_main(["no_such_dir"])


def test_plot_slices_missing_output_dir_exits_with_friendly_message() -> None:
    with pytest.raises(
        SystemExit,
        match=r'Fortran output files are missing under "no_such_dir"\.',
    ):
        plot_slices_main(["no_such_dir"])


def test_plot_performance_missing_profile_exits_with_friendly_message() -> None:
    with pytest.raises(
        SystemExit,
        match=r'Performance profile file is missing under "no_such_dir/performance_profile.csv"\.',
    ):
        plot_performance_main(["no_such_dir"])
