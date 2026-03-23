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


def test_preset_cli_save_extracts_sim_from_case_document(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)
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

[override.sim]
dt = 9.9e-9
""".strip()
        + "\n",
        encoding="utf-8",
    )

    beachx_main(["preset", "save", "sim/lab/run_baseline", "--section", "sim"])
    saved = capsys.readouterr()

    preset_path = beachx_home / "presets" / "sim" / "lab" / "run_baseline.toml"
    text = preset_path.read_text(encoding="utf-8")
    assert preset_path.exists()
    assert "saved=" in saved.out
    assert "[sim]" in text
    assert "dt = 9.9e-09" in text


def test_preset_cli_save_extracts_species_from_rendered_config(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "beach.toml").write_text(
        """
[sim]
dt = 2.0e-8
batch_duration_step = 60000.0
batch_count = 100
max_step = 10000
softening = 1.0e-6
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_far_correction = "auto"

[particles]
[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 5000
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
target_macro_particles_per_batch = -1
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]

[output]
write_files = true
write_mesh_potential = false
dir = "outputs/latest"
history_stride = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    beachx_main(
        [
            "preset",
            "save",
            "species/lab/ion",
            "--section",
            "species",
            "--index",
            "2",
            "--rendered",
        ]
    )
    saved = capsys.readouterr()

    preset_path = beachx_home / "presets" / "species" / "lab" / "ion.toml"
    assert preset_path.exists()
    beachx_main(["preset", "show", "species/lab/ion"])
    shown = capsys.readouterr()
    assert "section=species" in saved.out
    assert "[[particles.species]]" in shown.out
    assert "q_particle = 1.602176634e-19" in shown.out


def test_preset_cli_edit_opens_editor_for_user_preset(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)
    editor_script = tmp_path / "fake-editor.sh"
    editor_script.write_text(
        "#!/usr/bin/env bash\nprintf '\\n# edited\\n' >> \"$1\"\n",
        encoding="utf-8",
    )
    editor_script.chmod(0o755)
    monkeypatch.setenv("EDITOR", str(editor_script))

    beachx_main(["preset", "new", "output/lab/edit-me"])
    capsys.readouterr()

    beachx_main(["preset", "edit", "output/lab/edit-me"])
    edited = capsys.readouterr()

    preset_path = beachx_home / "presets" / "output" / "lab" / "edit-me.toml"
    assert "edited=" in edited.out
    assert "# edited" in preset_path.read_text(encoding="utf-8")


def test_preset_cli_edit_rejects_builtin_preset(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        beachx_main(["preset", "edit", "sim/periodic2_fmm"])

    assert "read-only" in str(excinfo.value)
    capsys.readouterr()


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
