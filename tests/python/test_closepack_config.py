from __future__ import annotations

import subprocess
import shlex
from pathlib import Path

import pytest

from beach.closepack_config import (
    ClosePackSpec,
    build_closepack_config,
    default_base_config,
    generate_closepack_sphere_templates,
    load_base_config,
    render_closepack_toml,
    unit_cell_pitch,
)


def test_generate_closepack_sphere_templates_matches_exp_bench_tree_geometry() -> None:
    spec = ClosePackSpec(layers=3, radius=0.292893, cells_x=1, cells_y=1)

    templates = generate_closepack_sphere_templates(spec)

    centers = [template["center"] for template in templates]
    expected = [
        [0.292893, 0.292893, 0.312893],
        [0.7071069999999999, 0.7071069999999999, 0.312893],
        [0.7071069999999999, 0.292893, 0.7271066956986759],
        [0.292893, 0.7071069999999999, 0.7271066956986759],
        [0.292893, 0.292893, 1.1413203913973517],
        [0.7071069999999999, 0.7071069999999999, 1.1413203913973517],
    ]
    for actual, expected_center in zip(centers, expected, strict=True):
        assert actual == pytest.approx(expected_center, abs=1.0e-6)


def test_build_closepack_config_updates_box_mesh_and_injection_face() -> None:
    spec = ClosePackSpec(layers=4, radius=0.2, cells_x=3, cells_y=2)

    config = build_closepack_config(spec)

    pitch = unit_cell_pitch(0.2)
    assert config["sim"]["box_min"] == [0.0, 0.0, 0.0]
    assert config["sim"]["box_max"] == [3.0 * pitch, 2.0 * pitch, 10.0]
    assert config["mesh"]["mode"] == "template"
    assert len(config["mesh"]["templates"]) == 1 + 2 * 4 * 3 * 2
    assert config["mesh"]["templates"][0]["size_x"] == 3.0 * pitch
    assert config["mesh"]["templates"][0]["size_y"] == 2.0 * pitch
    for species in config["particles"]["species"]:
        assert species["pos_low"] == [0.0, 0.0, 10.0]
        assert species["pos_high"] == [3.0 * pitch, 2.0 * pitch, 10.0]


def test_build_closepack_config_preserves_base_sim_and_output_values() -> None:
    base = default_base_config()
    base["sim"]["dt"] = 3.5e-8
    base["sim"]["field_solver"] = "fmm"
    base["sim"]["field_bc_mode"] = "periodic2"
    base["output"]["history_stride"] = 77
    base["particles"]["species"][0]["temperature_ev"] = 3.0

    spec = ClosePackSpec(layers=2, radius=0.18, cells_x=2, cells_y=2)
    config = build_closepack_config(spec, base_config=base)

    assert config["sim"]["dt"] == 3.5e-8
    assert config["sim"]["field_solver"] == "fmm"
    assert config["sim"]["field_bc_mode"] == "periodic2"
    assert config["output"]["history_stride"] == 77
    assert config["particles"]["species"][0]["temperature_ev"] == 3.0


def test_render_closepack_toml_roundtrips_with_toml_loader(tmp_path: Path) -> None:
    spec = ClosePackSpec(
        layers=2,
        radius=0.15,
        cells_x=2,
        cells_y=1,
        top_clearance=2.5,
    )
    config = build_closepack_config(spec, output_dir="outputs/generated")

    path = tmp_path / "generated.toml"
    path.write_text(render_closepack_toml(config, spec=spec), encoding="utf-8")
    reloaded = load_base_config(path)

    assert reloaded["sim"]["box_max"][0] == 2.0 * unit_cell_pitch(0.15)
    assert reloaded["output"]["dir"] == "outputs/generated"
    assert len(reloaded["mesh"]["templates"]) == 1 + 2 * 2 * 2 * 1


def test_generate_closepack_config_script_writes_expected_file(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    base_path = tmp_path / "base.toml"
    output_path = tmp_path / "beach.toml"
    base_path.write_text(
        """
[sim]
dt = 9.9e-8
batch_count = 100
max_step = 100000
tol_rel = 1.0e-8
softening = 1.0e-6
b0 = [0.0, 0.0, 0.0]
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
field_solver = "treecode"
batch_duration_step = 60000

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

[output]
write_files = true
dir = "./outputs/latest"
history_stride = 13
""".strip()
        + "\n",
        encoding="utf-8",
    )

    command = " ".join(
        [
            "python",
            shlex.quote(str(repo_root / "examples" / "generate_closepack_config.py")),
            shlex.quote(str(base_path)),
            "--layers",
            "3",
            "--radius",
            "0.292893",
            "--cells-x",
            "1",
            "--output",
            shlex.quote(str(output_path)),
            "--output-dir",
            "outputs/exp_closepack_generated",
        ]
    )
    completed = subprocess.run(
        ["bash", "-lc", command],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    generated = load_base_config(output_path)
    assert generated["sim"]["dt"] == 9.9e-8
    assert generated["output"]["dir"] == "outputs/exp_closepack_generated"
    assert generated["output"]["history_stride"] == 13
    assert len(generated["mesh"]["templates"]) == 7
