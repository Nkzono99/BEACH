# BEACH Known Failure Modes

## Install/Build

- Missing Fortran compiler or `fpm`: `pip install beach-bem` builds the Fortran binary and needs build tools.
- `make` profile mismatch: pip builds use `INSTALL_PROFILE=auto` and may fall back to `generic` unless `BEACH_PIP_FALLBACK_GENERIC=0`.
- OpenMP flags differ by compiler/environment. Prefer repository `make` targets in a development checkout.

## Config Parsing

- Unknown sections or keys are errors in the Fortran TOML parser.
- `"$schema"` before the first section is not accepted by the Fortran parser; use `#:schema ...` comments.
- `batch_duration` and `batch_duration_step` are mutually exclusive.
- `reservoir_face` and `photo_raycast` require resolved `batch_duration > 0`.
- `temperature_k` and `temperature_ev` are mutually exclusive.
- `e0` and the `e0_abs` angle form are mutually exclusive.
- `field_bc_mode = "periodic2"` requires a consistent two-axis periodic box boundary setup and is designed for FMM mode.

## Runtime/Output

- `tol_rel` does not cause early stop; it is only reported.
- Frequent history output, `write_potential_history`, large meshes, and many batches can create large CSV outputs.
- `output.resume = true` needs compatible `summary.txt`, `charges.csv`, and RNG state files in the same output directory.
- MPI resume requires the same world size as the previous run.
- Missing `charge_history.csv` usually means `history_stride <= 0`.
- Missing `mesh_potential.csv` usually means `output.write_mesh_potential = false`.

## Physics/Numerics

- Growing `survived_max_step` can indicate particles do not hit or escape within `sim.max_step`.
- Unexpectedly low absorption can come from injection face direction, box boundary, mesh placement, or drift velocity sign.
- Strong field spikes can come from tiny softening, close approach to charged element centers, or excessive accumulated charge.
- Python post-processing uses direct point-charge reconstruction and does not reproduce every Fortran FMM or periodic far-correction detail.
