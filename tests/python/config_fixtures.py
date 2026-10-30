"""beach.toml writers and base configs shared by the config and lint tests."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def free_space_beam_config() -> dict[str, Any]:
    """Return a small valid free-space config in the internal layout.

    One volume-seeded electron species, a conductor plane, ``fields.boundary="free"`` and the
    direct solver: the base that config and preflight tests mutate one key at a time.
    """

    return {
        "sim": {
            "dt": 2.0e-10,
            "batch_count": 150,
            "batch_duration": 1.0e-8,
            "max_step": 1500,
            "rng_seed": 12345,
            "field_solver": "direct",
        },
        "domain": {"box_min": [0.0, 0.0, 0.0], "box_max": [0.04, 0.04, 0.07], "periodic_axes": []},
        "field_boundary": {"mode": "free"},
        "particle_boundary": {"z_low": "open", "z_high": "open"},
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "q_particle": -1.602176634e-19,
                    "m_particle": 9.1093837139e-31,
                    "w_particle": 4000.0,
                    "npcls_per_step": 200,
                    "pos_low": [0.008, 0.008, 0.06],
                    "pos_high": [0.032, 0.032, 0.06],
                    "drift_velocity": [0.0, 0.0, -1.3262e6],
                    "temperature_ev": 0.01,
                },
            ]
        },
        "mesh": {
            "parts": [
                {
                    "kind": "plane",
                    "enabled": True,
                    "surface_model": "conductor",
                    "surface_side": "normal_plus",
                    "size_x": 0.02,
                    "size_y": 0.02,
                    "nx": 10,
                    "ny": 10,
                    "center": [0.02, 0.02, 0.01],
                }
            ],
        },
        "output": {
            "write_files": True,
            "dir": "outputs/latest",
            "history_stride": 1,
            "write_mesh_potential": True,
            "write_potential_history": True,
        },
    }


def _write_base_config(path: Path, *, field_bc_mode: str = "periodic2") -> None:
    path.write_text(
        f"""
[run]
batch_count = 2
rng_seed = 12345

[run.batch]
duration_steps = 10.0

[domain]
box_min = [
    0.0,
    0.0,
    0.0,
]
box_max = [
    1.0,
    1.0,
    10.0,
]
periodic_axes = [
    "x",
    "y",
]

[mesh]

[[mesh.parts]]
kind = "plane"
enabled = true
surface_side = "normal_plus"
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [
    0.5,
    0.5,
    0.0,
]

[particles.tracking]
dt_s = 2e-08
max_steps_per_particle = 100

[particles.boundary]
z_low = "open"
z_high = "open"

[[particles.species]]
charge_c = -1.602176634e-19
mass_kg = 9.10938356e-31

[particles.species.source]
mode = "volume_seed"

[particles.species.sampling]
volume_macro_particles_per_batch = 10

[fields]
boundary = "{field_bc_mode}"

[fields.solver]
method = "fmm"

[output]
enabled = true
dir = "outputs/latest"

[output.history]
stride_batches = 1
""".strip()
        + "\n",
        encoding="utf-8",
    )
