# BEACH Known Failure Modes

## Install/Build

- Missing Fortran compiler or `fpm`: `pip install beach-bem` builds the Fortran binary and needs build tools.
- `make` profile mismatch: pip builds use `INSTALL_PROFILE=auto` and may fall back to `generic` unless `BEACH_PIP_FALLBACK_GENERIC=0`.
- OpenMP flags differ by compiler/environment. Prefer repository `make` targets in a development checkout.

## Config Parsing

- Unknown sections or keys are errors in the Fortran TOML parser.
- `"$schema"` before the first section is not accepted by the Fortran parser; use `#:schema ...` comments.
- `run.batch.duration_s` and `run.batch.duration_steps` are mutually exclusive.
- `inflow` and `photo_raycast` require resolved `run.batch.duration_s > 0`.
- `distribution.temperature_k` and `distribution.temperature_ev` are mutually exclusive.
- `fields.external.electric_v_m` and the `fields.external.electric_magnitude_v_m` angle form are mutually exclusive.
- `fields.boundary = "periodic2"` requires a consistent two-axis periodic box boundary setup and is designed for FMM mode.

## Runtime/Output

- `output.diagnostics.charge_rel_change_threshold` does not cause early stop; it is only reported.
- Frequent history output, `output.history.potential`, large meshes, and many batches can create large CSV outputs.
- `[run.restart]` needs compatible `summary.txt`, `charges.csv`, and RNG state files in the same output directory.
- MPI resume requires the same world size as the previous run.
- Missing `charge_history.csv` usually means `output.history.stride_batches <= 0`.
- Missing `mesh_potential.csv` usually means `output.final.mesh_potential = false`.

## Physics/Numerics

- Growing `survived_max_step` can indicate particles do not hit or escape within `particles.tracking.max_steps_per_particle`.
- Unexpectedly low absorption can come from injection face direction, box boundary, mesh placement, or drift velocity sign.
- Strong field spikes can come from under-resolved panel geometry, sampling very near a charged edge or vertex, or excessive accumulated charge.
- Python field and potential post-processing uses the native P0 panel kernel; make sure the shared field-kernel library and run configuration match the output.
