from pathlib import Path

import pytest

from beach.cli_animate_fortran_history import main as animate_main
from beach.cli_inspect_fortran_output import main as inspect_main
from beach.cli_plot_coulomb_force_matrix import main as plot_coulomb_main
from beach.cli_plot_performance_profile import main as plot_performance_main
from beach.cli.plot_fortran_potential_slices import _load_sim_box
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


def test_plot_coulomb_missing_output_dir_exits_with_friendly_message() -> None:
    with pytest.raises(
        SystemExit,
        match=r'Fortran output files are missing under "no_such_dir"\.',
    ):
        plot_coulomb_main(["no_such_dir"])


def test_plot_performance_missing_profile_exits_with_friendly_message() -> None:
    with pytest.raises(
        SystemExit,
        match=r'Performance profile file is missing under "no_such_dir/performance_profile.csv"\.',
    ):
        plot_performance_main(["no_such_dir"])


def test_plot_slices_load_sim_box_defaults_periodic2_far_correction_to_oracle(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "beach.toml"
    config_path.write_text(
        "\n".join(
            [
                "[sim]",
                'field_solver = "fmm"',
                'field_bc_mode = "periodic2"',
                "box_min = [0.0, 0.0, -1.0]",
                "box_max = [1.0, 1.0, 1.0]",
                'bc_x_low = "periodic"',
                'bc_x_high = "periodic"',
                'bc_y_low = "periodic"',
                'bc_y_high = "periodic"',
                'bc_z_low = "open"',
                'bc_z_high = "open"',
                "field_periodic_image_layers = 2",
            ]
        ),
        encoding="utf-8",
    )

    _, _, _, periodic2 = _load_sim_box(config_path)

    assert periodic2 is not None
    assert periodic2["far_correction"] == "m2l_root_oracle"
    assert periodic2["ewald_layers"] == 4
