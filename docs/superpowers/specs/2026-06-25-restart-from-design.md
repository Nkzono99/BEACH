# restart_from Design

## Goal

Add `[output].restart_from` so a resumed run can read checkpoint files from one directory while writing all new outputs to `[output].dir`.

## Requirements

- Keep existing `[output].resume = true` behavior backward compatible when `restart_from` is omitted.
- When `restart_from` is present, read `summary.txt`, `charges.csv`, `rng_state*.txt`, and optional `macro_residuals*.csv` from `restart_from`.
- Always write new `summary.txt`, `charges.csv`, history, `mesh_potential.csv`, `performance_profile.csv`, `rng_state*.txt`, and `macro_residuals*.csv` to `output.dir`.
- Treat `sim.batch_count` as the cumulative target batch count, unchanged from existing resume behavior.
- Preserve existing checkpoint validation, including required-file failures and MPI world-size checks.
- Reject `restart_from` unless `output.resume = true`.
- Update schema, lint-visible validation, Python workload estimation, and user docs.

## Design

Use `[output].restart_from` rather than a new top-level `[restart]` section. This keeps restart input selection next to the existing `resume` flag and avoids adding another top-level config namespace.

The Fortran `app_config` type gains a fixed-length `output_restart_from` string. Its default is empty. Parsing `[output].restart_from` fills that string. Validation rejects a non-empty value when `resume_output` is false.

Runtime startup computes the checkpoint input directory as:

- `output_restart_from` when `output.resume = true` and the string is non-empty.
- `output_dir` otherwise.

Only the `load_restart_checkpoint` call uses this checkpoint directory. All output writer calls continue to use `output_dir`, so a continuation run does not need copied parent output files in its new output directory.

The JSON schema adds `output.restart_from` and validates that it implies `resume = true` and `write_files = true`. `beachx workload` reads checkpoint batch counts from `restart_from` when present so estimates match the Fortran runtime.

## Tests

- Fortran config parser test verifies `output.restart_from` is accepted and stored.
- Python workload test verifies resume batch counting reads `restart_from/summary.txt`, not `output.dir/summary.txt`.
- Schema/lint tests verify `restart_from` is accepted with `resume = true` and rejected without it.
- Existing restart tests continue to cover checkpoint validation and MPI rank-specific filenames.
