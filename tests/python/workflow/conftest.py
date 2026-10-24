"""Fixtures shared by the periodic validation tool tests."""

from __future__ import annotations

from pathlib import Path

import pytest
from validation_support import (
    _write_legacy_estimator_fixture,
    _write_release_mechanics_fixture,
)


@pytest.fixture
def archive_run(tmp_path: Path) -> Path:
    run = tmp_path / "archive" / "R20260625-0002"
    (run / "input").mkdir(parents=True)
    (run / "work" / "latest").mkdir(parents=True)
    (run / "input" / "beach.toml").write_text(
        """
[sim]
dt = 1.979898987322333e-12
batch_duration = 6.060915267313266
batch_count = 280000
max_step = 100000
tol_rel = 1.0e-8
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [9.899494936611664e-05, 9.899494936611664e-05, 0.0009899494936611664]
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

[particles]
[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 5000
pos_low = [0.0, 0.0, 0.0009899494936611664]
pos_high = [9.899494936611664e-05, 9.899494936611664e-05, 0.0009899494936611664]
drift_velocity = [0.0, 0.0, -4.0e5]
[particles.species.boundary_inflow]
z_high = "reservoir"

[mesh]
[[mesh.parts]]
kind = "sphere"
enabled = true
radius = 3.5e-05
center = [2.474873734152916e-05, 2.474873734152916e-05, 0.00013597484835343895]
n_lon = 24
n_lat = 10

[output]
write_files = true
dir = "work/latest"
history_stride = 1000
write_mesh_potential = true
resume = true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (run / "manifest.toml").write_text(
        """
[simulator_source]
executable = "/nonexistent/archived/beach"
exe_hash = "sha256:b82ac96cdcb6d0d5b200d4816908714a18432f52ad5d220242193ecfed2f7c63"

[job]
processes = 6
threads = 112
cores = 112
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (run / "work" / "SIMULATOR_VERSION.txt").write_text(
        "executable: /nonexistent/archived/beach\n1.3.0-v1.3.0\n",
        encoding="utf-8",
    )
    (run / "work" / "latest" / "mesh_sources.csv").write_text(
        "mesh_id,source_kind,template_kind,surface_model,epsilon_r,elem_count\n"
        "1,template,sphere,insulator,1.0,1\n",
        encoding="utf-8",
    )
    _write_release_mechanics_fixture(run)
    _write_legacy_estimator_fixture(run)
    return run


@pytest.fixture
def binary(tmp_path: Path) -> Path:
    path = tmp_path / "build" / "beach"
    path.parent.mkdir()
    path.write_bytes(b"deterministic-beach-test-binary\n")
    path.chmod(0o755)
    return path
