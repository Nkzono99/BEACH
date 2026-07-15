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

- Usually edit `beach.toml` directly, and let the Fortran parser normalize high-level notation while loading.
- The Fortran runtime reads `beach.toml`.
- `beachx lint beach.toml` performs pre-run checks including TOML parsing, JSON Schema, high-level notation, and known constraints.
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

- `sim.batch_count` is the number of batches for a normal run; with `output.resume=true`, it is the cumulative target batch count.
- `sim.max_step` is the per-particle step limit.
- `sim.tol_rel` is a monitoring metric, not an early-stop condition.
- With `output.resume = true`, the checkpoint is read from `output.dir` unless `restart_from` is set. MPI resume requires the same world size.
- `field_bc_mode = "periodic2"` handles two periodic axes. It normally uses FMM; the constrained direct `triangle_p0` split reference is the exception.
- `beachx config init` generates x/y periodic boundaries, FMM, one image layer, and a finite $3\times3$ image sum with no far correction.
