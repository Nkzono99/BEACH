from __future__ import annotations

import json
from pathlib import Path

import pytest

from beach.cli.main import main as beachx_main
from beach.config import (
    CaseSpecError,
    RenderValidationError,
    render_case_file,
)


def _write_project_preset(root: Path, name: str, body: str) -> None:
    path = root / ".beachx" / "presets" / f"{name}.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body.strip() + "\n", encoding="utf-8")


def _write_renderable_test_presets(root: Path) -> None:
    _write_project_preset(
        root,
        "sim/periodic2_fmm",
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
field_periodic_far_correction = "none"
""",
    )
    _write_project_preset(
        root,
        "species/solarwind_electron",
        """
[particles]
[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 10
""",
    )
    _write_project_preset(
        root,
        "species/solarwind_ion",
        """
[particles]
[[particles.species]]
source_mode = "volume_seed"
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
npcls_per_step = 10
""",
    )
    _write_project_preset(
        root,
        "mesh/plane_basic",
        """
[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.0]
""",
    )
    _write_project_preset(
        root,
        "mesh/closepack3",
        """
[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
size_x = 0.4
size_y = 0.4
nx = 8
ny = 8
center = [0.5, 0.5, 0.2]
""",
    )
    _write_project_preset(
        root,
        "output/standard",
        """
[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
""",
    )


def test_render_case_file_merges_builtin_presets_and_override(tmp_path: Path) -> None:
    _write_renderable_test_presets(tmp_path)
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

    assert result.case.use_presets == (
        "sim/periodic2_fmm",
        "species/solarwind_electron",
        "species/solarwind_ion",
        "mesh/plane_basic",
        "output/standard",
    )
    assert result.config["sim"]["dt"] == pytest.approx(1.5e-8)
    assert result.config["sim"]["b0"] == [1.0, 0.0, 0.0]
    assert result.config["output"]["dir"] == "outputs/custom"
    assert len(result.config["particles"]["species"]) == 2
    assert len(result.config["mesh"]["templates"]) == 1


def test_render_case_file_appends_override_species(tmp_path: Path) -> None:
    _write_renderable_test_presets(tmp_path)
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


def test_render_case_file_resolves_high_level_box_face_fraction_and_anchor(
    tmp_path: Path,
) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = ["output/standard"]

[override.sim]
dt = 2.0e-8
batch_duration_step = 10.0
batch_count = 2
max_step = 100
softening = 1.0e-6
use_box = true
box_origin = [1.0, 2.0, 3.0]
box_size = [2.0, 4.0, 6.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"

[override.mesh]
mode = "template"

[[override.mesh.templates]]
kind = "plane"
placement_mode = "box_anchor"
anchor = "z_low_face_center"
offset = [0.0, 0.0, 0.05]
size_mode = "box_fraction"
size_frac = [1.0, 0.5]
nx = 4
ny = 2

[[override.particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 100
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.25, 0.25]
uv_high = [0.75, 0.5]
drift_velocity = [0.0, 0.0, -4.0e5]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)
    template = result.config["mesh"]["templates"][0]
    species = result.config["particles"]["species"][0]

    assert result.config["sim"]["box_min"] == [1.0, 2.0, 3.0]
    assert result.config["sim"]["box_max"] == [3.0, 6.0, 9.0]
    assert "box_origin" not in result.config["sim"]
    assert "box_size" not in result.config["sim"]
    assert template["center"] == pytest.approx([2.0, 4.0, 3.05])
    assert template["size_x"] == pytest.approx(2.0)
    assert template["size_y"] == pytest.approx(2.0)
    assert "placement_mode" not in template
    assert "size_mode" not in template
    assert species["pos_low"] == pytest.approx([1.5, 3.0, 9.0])
    assert species["pos_high"] == pytest.approx([2.5, 4.0, 9.0])
    assert "inject_region_mode" not in species
    assert "uv_low" not in species


def test_render_case_file_resolves_mesh_group_scaling(tmp_path: Path) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = ["output/standard"]

[override.sim]
dt = 2.0e-8
batch_duration_step = 5.0
batch_count = 1
max_step = 50
softening = 1.0e-6
use_box = true
box_origin = [0.0, 0.0, 0.0]
box_size = [2.0, 4.0, 6.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"

[override.mesh]
mode = "template"

[override.mesh.groups.cavity_unit]
placement_mode = "box_anchor"
anchor = "box_center"
scale_from = "box_x"
scale_factor = 0.25

[[override.mesh.templates]]
group = "cavity_unit"
kind = "plane"
center_local = [0.0, 0.0, -1.0]
size_x = 4.0
size_y = 2.0
nx = 8
ny = 4

[[override.particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 10
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)
    template = result.config["mesh"]["templates"][0]

    assert template["center"] == pytest.approx([1.0, 2.0, 2.5])
    assert template["size_x"] == pytest.approx(2.0)
    assert template["size_y"] == pytest.approx(1.0)
    assert "group" not in template
    assert "center_local" not in template
    assert "groups" not in result.config["mesh"]


def test_render_case_file_rejects_box_origin_with_box_min_in_same_fragment(
    tmp_path: Path,
) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = ["output/standard"]

[override.sim]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0, 1.0, 1.0]
box_min = [0.0, 0.0, 0.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(CaseSpecError, match="sim.box_origin and sim.box_min"):
        render_case_file(case_path)


def test_render_case_file_rejects_high_level_inject_region_for_volume_seed(
    tmp_path: Path,
) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = ["output/standard"]

[override.sim]
dt = 2.0e-8
batch_duration_step = 1.0
batch_count = 1
max_step = 10
softening = 1.0e-6

[override.mesh]
mode = "template"

[[override.mesh.templates]]
kind = "plane"
size_x = 1.0
size_y = 1.0
nx = 1
ny = 1
center = [0.0, 0.0, 0.0]

[[override.particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 1
inject_region_mode = "face_fraction"
uv_low = [0.0, 0.0]
uv_high = [1.0, 1.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        CaseSpecError,
        match='source_mode must be "reservoir_face" or "photo_raycast"',
    ):
        render_case_file(case_path)


def test_render_case_file_rejects_negative_rendered_mesh_dimensions(
    tmp_path: Path,
) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = ["output/standard"]

[override.sim]
dt = 2.0e-8
batch_duration_step = 1.0
batch_count = 1
max_step = 10
softening = 1.0e-6
use_box = true
box_origin = [0.0, 0.0, 0.0]
box_size = [2.0, 2.0, 2.0]

[override.mesh]
mode = "template"

[[override.mesh.templates]]
kind = "plane"
placement_mode = "box_anchor"
anchor = "box_center"
size_mode = "box_fraction"
size_frac = [-0.5, 0.5]
nx = 1
ny = 1

[[override.particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(RenderValidationError, match="size_x must be > 0"):
        render_case_file(case_path)


def test_render_case_file_rejects_invalid_annulus_dimensions(tmp_path: Path) -> None:
    case_path = tmp_path / "case.toml"
    case_path.write_text(
        """
schema_version = 1
use_presets = ["output/standard"]

[override.sim]
dt = 2.0e-8
batch_duration_step = 1.0
batch_count = 1
max_step = 10
softening = 1.0e-6

[override.mesh]
mode = "template"

[[override.mesh.templates]]
kind = "annulus"
radius = 0.2
inner_radius = 0.3
n_theta = 8
n_r = 2
center = [0.0, 0.0, 0.0]

[[override.particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(RenderValidationError, match="inner_radius must be smaller than radius"):
        render_case_file(case_path)


def test_config_cli_init_render_validate_save_and_from(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    beachx_home = tmp_path / "beachx-home"
    monkeypatch.setenv("BEACHX_HOME", str(beachx_home))
    monkeypatch.chdir(tmp_path)
    _write_renderable_test_presets(tmp_path)

    beachx_main(["config", "init", "--title", "CLI Test Case"])
    init_streams = capsys.readouterr()
    assert "saved=case.toml" in init_streams.out
    assert (tmp_path / "case.toml").exists()
    assert (
        (tmp_path / "case.toml").read_text(encoding="utf-8").splitlines()[0]
        == "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.case.schema.json"
    )

    beachx_main(["config", "validate"])
    validate_streams = capsys.readouterr()
    assert "status=ok" in validate_streams.out

    beachx_main(["config", "render"])
    render_streams = capsys.readouterr()
    assert "saved=beach.toml" in render_streams.out
    assert (tmp_path / "beach.toml").exists()
    assert (
        (tmp_path / "beach.toml").read_text(encoding="utf-8").splitlines()[0]
        == "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json"
    )

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
    _write_renderable_test_presets(tmp_path)
    local_preset = tmp_path / ".beachx" / "presets" / "output" / "standard.toml"
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


def test_render_case_file_uses_parent_project_local_preset(
    tmp_path: Path,
) -> None:
    _write_renderable_test_presets(tmp_path)
    parent_preset = tmp_path / ".beachx" / "presets" / "output" / "standard.toml"
    parent_preset.write_text(
        """
[output]
write_files = true
dir = "outputs/parent-local"
history_stride = 3
""".strip()
        + "\n",
        encoding="utf-8",
    )
    case_dir = tmp_path / "runs" / "exp1"
    case_dir.mkdir(parents=True)
    case_path = case_dir / "case.toml"
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
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    assert result.config["output"]["dir"] == "outputs/parent-local"


def test_render_case_file_prefers_closest_parent_project_local_preset(
    tmp_path: Path,
) -> None:
    _write_renderable_test_presets(tmp_path)
    top_preset = tmp_path / ".beachx" / "presets" / "output" / "standard.toml"
    top_preset.write_text(
        """
[output]
write_files = true
dir = "outputs/top-local"
history_stride = 7
""".strip()
        + "\n",
        encoding="utf-8",
    )
    middle_dir = tmp_path / "runs"
    middle_preset = middle_dir / ".beachx" / "presets" / "output" / "standard.toml"
    middle_preset.parent.mkdir(parents=True)
    middle_preset.write_text(
        """
[output]
write_files = true
dir = "outputs/middle-local"
history_stride = 9
""".strip()
        + "\n",
        encoding="utf-8",
    )
    case_dir = middle_dir / "exp1"
    case_dir.mkdir(parents=True)
    case_path = case_dir / "case.toml"
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
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    assert result.config["output"]["dir"] == "outputs/middle-local"
    assert any(str(top_preset) in warning for warning in result.warnings)


def test_config_cli_diff_rendered_reports_changed_value(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_renderable_test_presets(tmp_path)
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


def test_case_schema_files_exist_and_example_case_has_schema_directive() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    case_schema = repo_root / "schemas" / "beach.case.schema.json"
    preset_schema = repo_root / "schemas" / "beach.preset.schema.json"
    example_case = repo_root / "examples" / "periodic2_basic" / "case.toml"

    loaded_case_schema = json.loads(case_schema.read_text(encoding="utf-8"))
    loaded_preset_schema = json.loads(preset_schema.read_text(encoding="utf-8"))

    assert loaded_case_schema["title"] == "BEACH Case File"
    assert loaded_preset_schema["title"] == "BEACH Preset Fragment"
    assert "box_origin" in loaded_case_schema["$defs"]["simFragment"]["properties"]
    assert "scale_from" in loaded_preset_schema["$defs"]["meshGroup"]["properties"]
    assert "allOf" in loaded_case_schema["$defs"]["simFragment"]
    assert "allOf" in loaded_preset_schema["$defs"]["speciesFragment"]
    assert (
        example_case.read_text(encoding="utf-8").splitlines()[0]
        == "#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.case.schema.json"
    )


# --- ID-based merge tests ---


def test_override_merges_template_by_id(tmp_path: Path) -> None:
    """Override with matching id deep-merges into preset template."""
    _write_renderable_test_presets(tmp_path)
    # Re-write mesh preset with an id field
    _write_project_preset(
        tmp_path,
        "mesh/plane_basic",
        """
[mesh]
mode = "template"

[[mesh.templates]]
id = "main_plane"
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.0]
""",
    )
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

[[override.mesh.templates]]
id = "main_plane"
nx = 40
ny = 40
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    # Should still have exactly 1 template (merged, not appended)
    assert len(result.config["mesh"]["templates"]) == 1
    tmpl = result.config["mesh"]["templates"][0]
    # Overridden values
    assert tmpl["nx"] == 40
    assert tmpl["ny"] == 40
    # Preserved original values
    assert tmpl["size_x"] == pytest.approx(1.0)
    assert tmpl["kind"] == "plane"
    assert tmpl["center"] == [0.5, 0.5, 0.0]
    # id should be stripped from rendered output
    assert "id" not in tmpl


def test_override_merges_species_by_id(tmp_path: Path) -> None:
    """Override with matching id deep-merges into preset species."""
    _write_renderable_test_presets(tmp_path)
    # Re-write species preset with an id field
    _write_project_preset(
        tmp_path,
        "species/solarwind_electron",
        """
[particles]
[[particles.species]]
id = "electron"
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 10
""",
    )
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
id = "electron"
npcls_per_step = 100
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    assert len(result.config["particles"]["species"]) == 1
    sp = result.config["particles"]["species"][0]
    assert sp["npcls_per_step"] == 100
    assert sp["q_particle"] == pytest.approx(-1.602176634e-19)
    assert "id" not in sp


def test_override_without_id_still_appends(tmp_path: Path) -> None:
    """Override without id uses classic append semantics."""
    _write_renderable_test_presets(tmp_path)
    _write_project_preset(
        tmp_path,
        "mesh/plane_basic",
        """
[mesh]
mode = "template"

[[mesh.templates]]
id = "main_plane"
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.0]
""",
    )
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

[[override.mesh.templates]]
kind = "plane"
size_x = 0.5
size_y = 0.5
nx = 5
ny = 5
center = [0.25, 0.25, 1.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    # id-less override appends a new template
    assert len(result.config["mesh"]["templates"]) == 2


def test_id_based_merge_with_unmatched_id_appends(tmp_path: Path) -> None:
    """Override with an id that doesn't match any base element appends."""
    _write_renderable_test_presets(tmp_path)
    _write_project_preset(
        tmp_path,
        "mesh/plane_basic",
        """
[mesh]
mode = "template"

[[mesh.templates]]
id = "main_plane"
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.0]
""",
    )
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

[[override.mesh.templates]]
id = "extra_plane"
kind = "plane"
size_x = 0.3
size_y = 0.3
nx = 3
ny = 3
center = [0.0, 0.0, 5.0]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = render_case_file(case_path)

    # Unmatched id appends a new template
    assert len(result.config["mesh"]["templates"]) == 2
