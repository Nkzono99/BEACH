# BEACH known failure modes

## install/build

- Fortran compiler or `fpm` missing: `pip install beach-bem` builds the Fortran binary and needs build tools.
- `make` profile mismatch: pip build uses `INSTALL_PROFILE=auto` and may fall back to `generic` unless `BEACH_PIP_FALLBACK_GENERIC=0`.
- OpenMP flags differ by compiler/environment. Prefer repository `make` targets in development checkout.

## config parsing

- Unknown section/key is an error in the Fortran TOML parser.
- `"$schema"` before the first section is not accepted by the Fortran parser; use `#:schema ...` comments.
- `batch_duration` and `batch_duration_step` are mutually exclusive.
- `reservoir_face` and `photo_raycast` require resolved `batch_duration > 0`.
- `temperature_k` and `temperature_ev` are mutually exclusive.
- `e0` and `e0_abs` angle form are mutually exclusive.
- `field_bc_mode = "periodic2"` requires consistent 2-axis periodic box boundary setup and is designed for FMM mode.

## runtime/output

- No early stop occurs from `tol_rel`; it is only reported.
- High `history_stride` frequency, `write_potential_history`, large mesh, and many batches can create large CSV outputs.
- `output.resume = true` needs compatible `summary.txt`, `charges.csv`, and RNG state files in the same output directory.
- MPI resume requires the same world size as the previous run.
- Missing `charge_history.csv` usually means `history_stride <= 0`.
- Missing `mesh_potential.csv` usually means `output.write_mesh_potential = false`.

## physics/numerics

- `survived_max_step` growth can indicate particles do not hit/escape within `sim.max_step`.
- Unexpectedly low absorption can come from injection face direction, box boundary, mesh placement, or drift velocity sign.
- Strong field spikes can come from tiny softening, close approach to charged element centers, or excessive accumulated charge.
- Python post-processing uses direct point-charge reconstruction and does not reproduce every Fortran FMM/periodic far correction detail.
