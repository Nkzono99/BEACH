"""beach.toml writers shared by the config and lint tests."""

from __future__ import annotations

from pathlib import Path


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
