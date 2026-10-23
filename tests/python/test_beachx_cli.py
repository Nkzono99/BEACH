from __future__ import annotations

import pytest

from beach.cli.main import main as beachx_main


def test_beachx_help_lists_all_subcommands(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["--help"])

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "inspect" in captured.out
    assert "animate" in captured.out
    assert "coulomb" in captured.out
    assert "mobility" in captured.out
    assert "kernel-forces" in captured.out
    assert "object-detachment" in captured.out
    assert "slices" in captured.out
    assert "workload" in captured.out
    assert "profile" in captured.out
    assert "lint" in captured.out
    assert "config" in captured.out
    assert "model" in captured.out
    assert "kinetic-response" not in captured.out


def test_beachx_model_help_lists_available_models(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["model", "--help"])

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "close-pack" in captured.out


def test_beachx_config_help_lists_available_subcommands(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["config", "--help"])

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "init" in captured.out
    assert "diff" in captured.out


def test_beachx_model_close_pack_missing_base_config_exits_with_friendly_message(
    tmp_path,
    monkeypatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["model", "close-pack", "--layers", "4", "--radius", "0.2", "--cells-x", "2"])

    assert str(excinfo.value) == "base config file not found: beach.toml"


@pytest.mark.parametrize("command", ["kernel-forces", "mobility"])
def test_beachx_analysis_command_missing_output_exits_with_friendly_message(
    command: str,
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        beachx_main([command, "no_such_dir"])

    assert (
        str(excinfo.value)
        == 'Fortran output files are missing under "no_such_dir". Expected at least summary.txt and charges.csv.'
    )
