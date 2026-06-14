from __future__ import annotations

from pathlib import Path

import pytest

from beach.cli.main import main as beachx_main
from beach.config import ConfigError, RenderValidationError, load_config_file


def _write_base_config(path: Path, *, field_bc_mode: str = "periodic2") -> None:
    path.write_text(
        f"""
[sim]
dt = 2.0e-8
batch_duration_step = 10.0
batch_count = 2
max_step = 100
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
field_bc_mode = "{field_bc_mode}"

[particles]
[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 10

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

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )


def test_load_config_file_accepts_direct_beach_toml(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)

    result = load_config_file(config_path)

    assert result["sim"]["field_bc_mode"] == "periodic2"
    assert result["particles"]["species"][0]["npcls_per_step"] == 10
    assert result["mesh"]["templates"][0]["kind"] == "plane"


def test_load_config_file_resolves_high_level_notation(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    config_path.write_text(
        """
[sim]
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

[particles]
[[particles.species]]
source_mode = "reservoir_face"
number_density_m3 = 1.0
temperature_k = 0.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
w_particle = 1.0
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.25, 0.5]
uv_high = [0.75, 1.0]
drift_velocity = [0.0, 0.0, -1.0]

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
size_mode = "box_fraction"
size_frac = [0.5, 0.25]
nx = 20
ny = 20
placement_mode = "box_anchor"
anchor = "box_center"

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )

    result = load_config_file(config_path)

    assert result["sim"]["box_min"] == [1.0, 2.0, 3.0]
    assert result["sim"]["box_max"] == [3.0, 6.0, 9.0]
    species = result["particles"]["species"][0]
    assert species["pos_low"] == [1.5, 4.0, 9.0]
    assert species["pos_high"] == [2.5, 6.0, 9.0]
    template = result["mesh"]["templates"][0]
    assert template["center"] == [2.0, 4.0, 6.0]
    assert template["size_x"] == pytest.approx(1.0)
    assert template["size_y"] == pytest.approx(1.0)


def test_load_config_file_rejects_legacy_composition_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "legacy_composition.toml"
    config_path.write_text(
        """
schema_version = 1
use_presets = ["sim/periodic2_fmm"]
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigError, match="reserved top-level key"):
        load_config_file(config_path)


def test_load_config_file_rejects_conductor_with_periodic2(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8")
    text = text.replace(
        "center = [0.5, 0.5, 0.0]",
        'center = [0.5, 0.5, 0.0]\nsurface_model = "conductor"',
    )
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(RenderValidationError, match='surface_model="conductor"'):
        load_config_file(config_path)


def test_load_config_file_rejects_nonfinite_template_scalar(tmp_path: Path) -> None:
    config_path = tmp_path / "beach.toml"
    _write_base_config(config_path)
    text = config_path.read_text(encoding="utf-8").replace("size_x = 1.0", "size_x = inf")
    config_path.write_text(text, encoding="utf-8")

    with pytest.raises(RenderValidationError, match="mesh.templates\\[1\\].size_x"):
        load_config_file(config_path)


def test_config_cli_init_render_validate_and_diff(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.chdir(tmp_path)

    beachx_main(["config", "init"])
    init_streams = capsys.readouterr()
    assert "saved=beach.toml" in init_streams.out
    assert (tmp_path / "beach.toml").exists()

    beachx_main(["config", "validate"])
    validate_streams = capsys.readouterr()
    assert "status=ok" in validate_streams.out

    beachx_main(["config", "render", "--stdout"])
    render_streams = capsys.readouterr()
    assert "[sim]" in render_streams.out
    assert "use_presets" not in render_streams.out

    modified = tmp_path / "modified.toml"
    modified.write_text(
        (tmp_path / "beach.toml").read_text(encoding="utf-8").replace(
            "batch_count = 100", "batch_count = 101"
        ),
        encoding="utf-8",
    )
    beachx_main(["config", "diff", "beach.toml", str(modified)])
    diff_streams = capsys.readouterr()
    assert "sim.batch_count: 100 -> 101" in diff_streams.out
