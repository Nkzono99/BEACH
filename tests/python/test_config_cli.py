from __future__ import annotations

from pathlib import Path

import pytest

from beach.cli.main import main as beachx_main
from beach.config import DEFAULT_PRESET_NAMES, render_case_file


def test_render_case_file_merges_builtin_presets_and_override(tmp_path: Path) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = [
  "sim/periodic2_fmm",
  "species/solarwind_electron",
  "species/solarwind_ion",
  "mesh/plane_basic",
  "output/standard",
]

[override.sim]
dt = 1.5e-8
b0 = [1.0, 0.0, 0.0]

[override.output]
dir = "outputs/custom"
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    assert result.case.use_presets == DEFAULT_PRESET_NAMES
    assert result.config["sim"]["dt"] == pytest.approx(1.5e-8)
    assert result.config["sim"]["b0"] == [1.0, 0.0, 0.0]
    assert result.config["output"]["dir"] == "outputs/custom"
    assert len(result.config["particles"]["species"]) == 2
    assert len(result.config["mesh"]["templates"]) == 1


def test_render_case_file_appends_override_species(tmp_path: Path) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = [
  "sim/periodic2_fmm",
  "species/solarwind_electron",
  "mesh/plane_basic",
  "output/standard",
]

[[override.particles.species]]
source_mode = "volume_seed"
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
npcls_per_step = 10
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    assert len(result.config["particles"]["species"]) == 2
    assert result.config["particles"]["species"][1]["source_mode"] == "volume_seed"


def test_config_cli_init_render_validate_save_and_from(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)

    beachx_main(["config", "init", "--title", "CLI Test Case"])
    init_streams = capsys.readouterr()
    assert "saved=case.toml" in init_streams.out
    assert (tmp_path / "case.toml").exists()

    beachx_main(["config", "validate"])
    validate_streams = capsys.readouterr()
    assert "status=ok" in validate_streams.out

    beachx_main(["config", "render"])
    render_streams = capsys.readouterr()
    assert "saved=beach.toml" in render_streams.out
    assert (tmp_path / "beach.toml").exists()

    beachx_main(["config", "save", "baseline"])
    save_streams = capsys.readouterr()
    assert "baseline.toml" in save_streams.out

    beachx_main(["config", "list-saved"])
    listed = capsys.readouterr()
    assert listed.out.strip() == "baseline"

    copied_dir = tmp_path / "copied"
    copied_dir.mkdir()
    monkeypatch.chdir(copied_dir)
    beachx_main(["config", "init", "--from", "baseline"])
    copied_streams = capsys.readouterr()
    assert "source_saved_case" in copied_streams.out
    assert (copied_dir / "case.toml").exists()


def test_config_cli_render_uses_project_local_preset_shadow(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.chdir(tmp_path)
    local_preset = tmp_path / ".beachx" / "presets" / "output" / "standard.toml"
    local_preset.parent.mkdir(parents=True)
    local_preset.write_text(
        """
[output]
write_files = true
dir = "outputs/project-local"
history_stride = 5
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "case.toml").write_text(
        """
schema_version = 1
use_presets = [
  "sim/periodic2_fmm",
  "species/solarwind_electron",
  "species/solarwind_ion",
  "mesh/plane_basic",
  "output/standard",
]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    beachx_main(["config", "render"])
    streams = capsys.readouterr()
    rendered = (tmp_path / "beach.toml").read_text(encoding="utf-8")

    assert "project-local" in rendered
    assert 'warning: preset "output/standard"' in streams.err


def test_config_cli_diff_rendered_reports_changed_value(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    left = tmp_path / "left.toml"
    right = tmp_path / "right.toml"
    left.write_text(
        """
schema_version = 1
use_presets = [
  "sim/periodic2_fmm",
  "species/solarwind_electron",
  "species/solarwind_ion",
  "mesh/plane_basic",
  "output/standard",
]
""".strip()
        + "\n",
        encoding="utf-8",
    )
    right.write_text(
        """
schema_version = 1
use_presets = [
  "sim/periodic2_fmm",
  "species/solarwind_electron",
  "species/solarwind_ion",
  "mesh/plane_basic",
  "output/standard",
]

[override.sim]
dt = 1.0e-8
""".strip()
        + "\n",
        encoding="utf-8",
    )

    beachx_main(["config", "diff", "--rendered", str(left), str(right)])
    streams = capsys.readouterr()

    assert "sim.dt" in streams.out
    assert "1e-08" in streams.out
