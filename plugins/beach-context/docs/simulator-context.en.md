# BEACH Simulator Context

BEACH combines BEM-based surface charging with test-particle tracking. The simulation core is Fortran (`src/`, `app/`), while the Python package (`beach/`) provides configuration helpers, output loading, analysis, and visualization.

Core v1 behavior:

- Store charge on triangular boundary elements.
- Evaluate the electrostatic field from current element charges.
- Advance particles with a Boris pusher.
- Treat the first segment-triangle intersection as the collision.
- Absorb collided particles and deposit `q_particle * w_particle` onto the hit element.
- Commit charge deltas, statistics, and histories once per batch.
- The standard surface model is insulator accumulation. Conductor and resistive models are reserved extension points.

Configuration entry points:

- Usually edit `case.toml` plus presets, then generate `beach.toml` with `beachx config render`.
- The Fortran runtime reads only the rendered `beach.toml`.
- `beachx config validate` performs pre-run checks including preset resolution and high-level notation expansion.
- The final key specification is `references/fortran_parameter_file.md`.

Primary outputs:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv`
- `potential_history.csv`
- `mesh_potential.csv`
- `rng_state.txt`
- `macro_residuals.csv`
- `performance_profile.csv` when `BEACH_PROFILE=1`

Analysis entry points:

- CLI: `beachx inspect`, `beachx animate`, `beachx slices`, `beachx coulomb`, `beachx mobility`, `beachx kernel-forces`
- Python: `from beach import Beach`
- Fortran field kernel: run `make build-kernel`, then use `beachx kernel-forces`

Operational notes:

- `sim.batch_count` is the number of batches to run now.
- `sim.max_step` is the per-particle step limit.
- `sim.tol_rel` is a monitoring metric, not an early-stop condition.
- `output.resume = true` reads `summary.txt`, `charges.csv`, and `rng_state.txt` from the same `output.dir`. MPI resume requires the same world size.
- `field_bc_mode = "periodic2"` handles two periodic axes and requires care with FMM and box boundary settings.
