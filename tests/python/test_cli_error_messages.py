import pytest

from beach.cli_animate_fortran_history import main as animate_main
from beach.cli_inspect_fortran_output import main as inspect_main


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
