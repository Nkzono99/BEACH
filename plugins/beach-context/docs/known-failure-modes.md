# BEACH known failure modes

## install/build

- Fortran compiler or `fpm` missing: `pip install beach-bem` builds the Fortran binary and needs build tools.
- `make` profile mismatch: pip build uses `INSTALL_PROFILE=auto` and may fall back to `generic` unless `BEACH_PIP_FALLBACK_GENERIC=0`.
- OpenMP flags differ by compiler/environment. Prefer repository `make` targets in development checkout.

## config parsing

- Unknown section/key is an error in the Fortran TOML parser.
- `"$schema"` before the first section is not accepted by the Fortran parser; use `#:schema ...` comments.
- `run.batch.duration_s` and `run.batch.duration_steps` are mutually exclusive.
- `inflow` and `photo_raycast` require resolved `run.batch.duration_s > 0`.
- `distribution.temperature_k` and `distribution.temperature_ev` are mutually exclusive.
- `fields.external.electric_v_m` and `fields.external.electric_magnitude_v_m` angle form are mutually exclusive.
- `fields.boundary = "periodic2"` requires consistent 2-axis periodic box boundary setup and is designed for FMM mode.

## runtime/output

- No early stop occurs from `output.diagnostics.charge_rel_change_threshold`; it is only reported.
- High `output.history.stride_batches` frequency, `output.history.potential`, large mesh, and many batches can create large CSV outputs.
- `[run.restart]` needs compatible `summary.txt`, `charges.csv`, and RNG state files in the same output directory.
- MPI resume requires the same world size as the previous run.
- `cached_kneq0` の初回実行で最初の batch 行がなかなか出ない: 遠方演算子の生成中か、同じ cache の lock（`operator-<fingerprint>.bin.lock`）待ち。`summary.txt` の `periodic2_cache_hit` / `periodic2_operator_build_count` で確認し、次からは `beach --prepare-periodic-cache` を別 job で先に実行する。
- Missing `charge_history.csv` usually means `output.history.stride_batches <= 0`.
- Missing `mesh_potential.csv` usually means `output.final.mesh_potential = false`.

## physics/numerics

- `survived_max_step` growth can indicate particles do not hit/escape within `particles.tracking.max_steps_per_particle`.
- Unexpectedly low absorption can come from injection face direction, box boundary, mesh placement, or drift velocity sign.
- Strong field spikes can come from under-resolved panel geometry, sampling very near a charged edge or vertex, or excessive accumulated charge.
- Python field and potential post-processing uses the native P0 panel kernel; make sure the shared field-kernel library and run configuration match the output.
