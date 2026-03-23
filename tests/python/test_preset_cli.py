from __future__ import annotations

from pathlib import Path

import pytest

from beach.cli.main import main as beachx_main


def test_preset_cli_new_from_builtin_creates_user_preset(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)

    beachx_main(
        [
            "preset",
            "new",
            "sim/lab/periodic2_fast",
            "--from",
            "sim/periodic2_fmm",
        ]
    )
    created = capsys.readouterr()

    preset_path = beachx_home / "presets" / "sim" / "lab" / "periodic2_fast.toml"
    assert preset_path.exists()
    assert "scope=user" in created.out
    assert "from=sim/periodic2_fmm" in created.out
    assert (
        preset_path.read_text(encoding="utf-8").splitlines()[0]
        == "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.preset.schema.json"
    )

    beachx_main(["preset", "path", "sim/lab/periodic2_fast"])
    path_streams = capsys.readouterr()
    assert path_streams.out.strip() == str(preset_path)

    beachx_main(["preset", "validate", "sim/lab/periodic2_fast"])
    validate_streams = capsys.readouterr()
    assert "status=ok" in validate_streams.out


def test_preset_cli_new_local_creates_project_local_skeleton(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.chdir(tmp_path)

    beachx_main(["preset", "new", "output/project/debug", "--local"])
    created = capsys.readouterr()

    preset_path = (
        tmp_path / ".beachx" / "presets" / "output" / "project" / "debug.toml"
    )
    assert preset_path.exists()
    assert "scope=project-local" in created.out
    assert (
        preset_path.read_text(encoding="utf-8").splitlines()[0]
        == "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.preset.schema.json"
    )

    beachx_main(["preset", "show", "output/project/debug"])
    show_streams = capsys.readouterr()
    assert "[output]" in show_streams.out


def test_preset_cli_list_prefers_project_local_over_user_and_builtin(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)

    user_preset = beachx_home / "presets" / "output" / "standard.toml"
    user_preset.parent.mkdir(parents=True)
    user_preset.write_text(
        """
[output]
dir = "outputs/user"
""".strip()
        + "\n",
        encoding="utf-8",
    )

    local_preset = tmp_path / ".beachx" / "presets" / "output" / "standard.toml"
    local_preset.parent.mkdir(parents=True)
    local_preset.write_text(
        """
[output]
dir = "outputs/local"
""".strip()
        + "\n",
        encoding="utf-8",
    )

    beachx_main(["preset", "list"])
    listed = capsys.readouterr()
    output_lines = [line for line in listed.out.splitlines() if line.startswith("output/standard")]

    assert len(output_lines) == 1
    assert "\tproject-local\t" in output_lines[0]
    assert str(local_preset) in output_lines[0]


def test_preset_cli_show_warns_when_shadowing_occurs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)

    user_preset = beachx_home / "presets" / "mesh" / "plane_basic.toml"
    user_preset.parent.mkdir(parents=True)
    user_preset.write_text(
        """
[mesh]
mode = "template"
""".strip()
        + "\n",
        encoding="utf-8",
    )

    local_preset = tmp_path / ".beachx" / "presets" / "mesh" / "plane_basic.toml"
    local_preset.parent.mkdir(parents=True)
    local_preset.write_text(
        """
[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
""".strip()
        + "\n",
        encoding="utf-8",
    )

    beachx_main(["preset", "show", "mesh/plane_basic"])
    streams = capsys.readouterr()

    assert "[mesh]" in streams.out
    assert 'warning: preset "mesh/plane_basic"' in streams.err


def test_builtin_presets_have_schema_directive() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    preset_root = repo_root / "beach" / "config" / "presets"

    for path in preset_root.rglob("*.toml"):
        assert (
            path.read_text(encoding="utf-8").splitlines()[0]
            == "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.preset.schema.json"
        )
