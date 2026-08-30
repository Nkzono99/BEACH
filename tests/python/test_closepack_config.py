from __future__ import annotations

import copy
import math
import subprocess
import sys
from pathlib import Path

import pytest

from beach.closepack_config import (
    ClosePackSpec,
    build_closepack_config,
    default_base_config,
    generate_closepack_sphere_templates,
    load_base_config,
    render_closepack_toml,
)
from beach.config.schema import load_schema, schema_errors


def test_generate_closepack_sphere_templates_follows_analytic_lattice() -> None:
    spec = ClosePackSpec(layers=3, radius=0.292893, cells_x=1, cells_y=1)

    templates = generate_closepack_sphere_templates(spec)

    sqrt2 = math.sqrt(2.0)
    pitch = spec.radius * (2.0 + sqrt2)
    upper = pitch - spec.radius
    z0 = spec.floor_z + spec.radius
    dz = spec.radius * sqrt2
    centers = [template["center"] for template in templates]
    expected = [
        [spec.radius, spec.radius, z0],
        [upper, upper, z0],
        [upper, spec.radius, z0 + dz],
        [spec.radius, upper, z0 + dz],
        [spec.radius, spec.radius, z0 + 2.0 * dz],
        [upper, upper, z0 + 2.0 * dz],
    ]
    for actual, expected_center in zip(centers, expected, strict=True):
        assert actual == pytest.approx(expected_center, rel=1.0e-14, abs=1.0e-15)
    assert all(
        template["kind"] == "sphere"
        and template["surface_side"] == "outward_closed"
        and template["radius"] == spec.radius
        and template["n_lon"] == spec.sphere_n_lon
        and template["n_lat"] == spec.sphere_n_lat
        for template in templates
    )


def test_build_closepack_config_updates_geometry_and_injection_without_mutating_base() -> None:
    base = default_base_config()
    base["domain"]["box_min"] = [1.0, -2.0, 0.0]
    base["domain"]["box_max"] = [9.0, 9.0, 100.0]
    base["sim"]["dt"] = 3.5e-8
    base["sim"]["field_solver"] = "fmm"
    base["output"]["history_stride"] = 77
    base["mesh"]["obj_path"] = "obsolete.obj"
    reservoir_species = base["particles"]["species"][0]
    reservoir_species["source_mode"] = "reservoir_face"
    reservoir_species.pop("npcls_per_step")
    reservoir_species.pop("boundary_inflow")
    reservoir_species["inject_face"] = "z_high"
    reservoir_species["pos_low"] = [0.0, 0.0, 100.0]
    reservoir_species["pos_high"] = [1.0, 1.0, 100.0]
    base_before = copy.deepcopy(base)

    spec = ClosePackSpec(
        layers=4,
        radius=0.2,
        cells_x=3,
        cells_y=2,
        top_clearance=0.8,
    )
    config = build_closepack_config(spec, base_config=base)
    preserved_top = build_closepack_config(
        ClosePackSpec(layers=4, radius=0.2, cells_x=3, cells_y=2),
        base_config=base,
    )

    pitch = spec.radius * (2.0 + math.sqrt(2.0))
    expected_top = (
        spec.floor_z
        + 2.0 * spec.radius
        + (spec.layers - 1) * spec.radius * math.sqrt(2.0)
        + 0.8
    )
    expected_box_max = [1.0 + 3.0 * pitch, -2.0 + 2.0 * pitch, expected_top]

    assert base == base_before
    assert preserved_top["domain"]["box_max"][2] == 100.0
    assert config["domain"]["box_min"] == [1.0, -2.0, 0.0]
    assert config["domain"]["box_max"] == pytest.approx(expected_box_max)
    assert config["sim"]["dt"] == 3.5e-8
    assert config["sim"]["field_solver"] == "fmm"
    assert config["field_boundary"] == base["field_boundary"]
    assert config["output"]["history_stride"] == 77
    assert config["mesh"]["mode"] == "template"
    assert "obj_path" not in config["mesh"]
    assert len(config["mesh"]["templates"]) == 1 + 2 * 4 * 3 * 2
    floor = config["mesh"]["templates"][0]
    assert floor["kind"] == "plane"
    assert floor["surface_side"] == "normal_plus"
    assert floor["size_x"] == pytest.approx(3.0 * pitch)
    assert floor["size_y"] == pytest.approx(2.0 * pitch)
    assert floor["center"] == pytest.approx(
        [1.0 + 1.5 * pitch, -2.0 + pitch, spec.floor_z]
    )
    generated_reservoir = config["particles"]["species"][0]
    assert generated_reservoir["pos_low"] == pytest.approx(
        [1.0, -2.0, expected_top]
    )
    assert generated_reservoir["pos_high"] == pytest.approx(expected_box_max)
    boundary_species = config["particles"]["species"][1]
    assert boundary_species == base["particles"]["species"][1]
    assert boundary_species["source_mode"] == "volume_seed"
    assert boundary_species["npcls_per_step"] == 0
    assert boundary_species["boundary_inflow"] == {"z_high": "reservoir"}
    assert "pos_low" not in boundary_species and "pos_high" not in boundary_species


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
    schema, _ = load_schema()

    assert reloaded == config
    assert schema_errors(reloaded, schema) == []


def test_generate_closepack_config_script_writes_expected_file(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    base_path = tmp_path / "base.toml"
    output_path = tmp_path / "beach.toml"
    base = default_base_config()
    base["sim"]["dt"] = 9.9e-8
    base["output"]["history_stride"] = 13
    base_path.write_text(render_closepack_toml(base), encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            str(repo_root / "examples" / "generate_closepack_config.py"),
            str(base_path),
            "--layers",
            "3",
            "--radius",
            "0.292893",
            "--cells-x",
            "1",
            "--output",
            str(output_path),
            "--output-dir",
            "outputs/exp_closepack_generated",
        ],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert "sphere_count=6" in completed.stdout
    generated = load_base_config(output_path)
    assert generated["sim"]["dt"] == 9.9e-8
    assert generated["output"]["dir"] == "outputs/exp_closepack_generated"
    assert generated["output"]["history_stride"] == 13
    assert len(generated["mesh"]["templates"]) == 7
