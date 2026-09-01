title: Troubleshooting

Lang: [English](Troubleshooting.en.md) | [日本語](Troubleshooting.md)

# Troubleshooting

Choose the section closest to the symptom, run its checks, and compare the result with the stated expected state.
Do not force an unsupported model or loosen numerical tolerances merely to bypass an error.

Do not run `beach`, `mpirun`, or long analyses directly on a KUDPC login node.
Run the execution commands below either locally or inside a compute-node allocation.

The general examples below use `outputs/latest`. Replace it with `output.dir` from the actual `beach.toml`.
The official tutorial writes to `outputs/tutorial`.

## `beach` or `beachx` is not found

**Check:** Identify the active Python, build tool, Fortran compiler, and BEACH commands.

```bash
command -v python
command -v make
command -v gfortran
command -v beach
command -v beachx
python -m site --user-base
beach --version
beachx --help
```

`gfortran` is an example. If `FC` or `FPM_FC` selects another Fortran compiler, check that command instead.
Check `command -v fpm` only when installing a checkout with `--no-build-isolation` or invoking `make` directly.

**Expected:** Each required command prints an absolute executable path. `beach --version` and `beachx --help` exit with
code `0`.

**Safe action:**

- If only `beach` or `beachx` is missing, reinstall into the active Python environment by following
  [Installation](Installation.en.html).
- If the package was installed under the user base, add the displayed user-base `bin` directory to `PATH`.
- If `make` or the compiler is missing, install the OS package or load the HPC site's compiler module before reinstalling.
- Avoid mixing `pip` from another Python environment; use `python -m pip ...`.

## `beachx lint` fails

**Check:** Validate the same file that will be passed to the simulator.

```bash
beachx lint beach.toml
```

**Expected:** The TOML, JSON Schema, and BEACH semantic checks pass, ending with `status=ok`.

**Safe action:** Fix the first reported unknown key, type, range, or mutual-exclusion error, then run lint again.
Use [Common configuration mistakes](Configuration.en.html#4-common-mistakes) and
[Input parameters](Parameters.en.html) to check the key's table, type, and unit.

## Lint passes but `beach` stops during execution

**Check:** Save standard output and standard error, then inspect the exit code and final diagnostics.

```bash
beach beach.toml > beach-run.log 2>&1
run_status=$?
printf 'exit_code=%s\n' "$run_status"
tail -n 50 beach-run.log
```

**Expected:** The command reports `exit_code=0`, and the log ends with run statistics and `results written to ...`.

**Safe action:** Fix the first `ERROR STOP` message or the diagnostic immediately before it.
Fortran physics validation can be stricter than the Python schema, so lint success alone does not guarantee that a model
combination is runnable. Make the configuration satisfy the documented compatibility and prerequisite rules; do not
loosen a tolerance solely to suppress the failure.

## The output directory is missing or `beachx inspect` cannot read it

**Check:** Inspect the launch directory, `[output]` settings, and required final files.

```bash
pwd
grep -A 10 '^\[output\]' beach.toml
ls -l outputs/latest/summary.txt outputs/latest/charges.csv
beachx inspect outputs/latest
```

**Expected:** `write_files=true`; `summary.txt` and `charges.csv` exist under `dir`; and `beachx inspect` exits with
code `0`.

**Safe action:** Resolve a relative `output.dir` from the working directory where `beach` was launched.
Fix the simulator exit status and log before retrying. Do not create empty `summary.txt` or `charges.csv` placeholders.

## The run finishes but particle statistics differ from expectations

**Check:** Read completion and terminal outcomes separately.

```bash
beachx inspect outputs/latest
grep -E '^(processed_particles|absorbed|escaped|escaped_boundary|survived_max_step|multiple_box_events_soft_discarded|batches)=' \
  outputs/latest/summary.txt
```

**Expected:** `batches` equals `sim.batch_count`, and the particle counts can be explained from the source, absorption,
boundary escape, `survived_max_step`, and unresolved soft discards. There is no universal pass value for
`survived_max_step` across all cases.

**Safe action:** If `survived_max_step` affects the quantity of interest, vary `dt`, `max_step`, box size, and injection
velocity one at a time and check dependence. Do not relabel unresolved particles as absorbed or escaped after the run.
Use [Validating Simulation Results](ValidationGuide.en.html) for acceptance.

## History is empty or too large

**Check:** Compare the output controls with the actual file sizes and row counts.

```bash
grep -A 10 '^\[output\]' beach.toml
ls -lh outputs/latest/*history.csv
wc -l outputs/latest/*history.csv
```

**Expected:** `history_stride > 0` creates `charge_history.csv`.
`potential_history.csv` is created only with `write_potential_history=true` and `history_stride > 0`.

**Safe action:** For missing history, check `history_stride`, `batch_count`, and the actual `output.dir`.
For oversized history, increase `history_stride` and enable `write_potential_history` only for cases that require
potential history.

## Cannot resume from a checkpoint

**Check:** First inspect the directory named by `restart_from`. For a same-directory resume without `restart_from`,
inspect `output.dir` itself. The commands below use a serial `outputs/latest` example.

```bash
checkpoint_dir=outputs/latest
grep -A 12 '^\[sim\]' resume.toml
grep -A 12 '^\[output\]' resume.toml
ls -l "$checkpoint_dir/summary.txt" \
  "$checkpoint_dir/charges.csv" \
  "$checkpoint_dir/checkpoint_complete.txt"
ls -l "$checkpoint_dir"/rng_state*.txt
grep -E '^(checkpoint_schema_version|batches|model_fingerprint|mesh_fingerprint|species_fingerprint|mpi_world_size)=' \
  "$checkpoint_dir/summary.txt"
sed -n '1,80p' "$checkpoint_dir/checkpoint_complete.txt"
```

**Expected:** A current-schema checkpoint satisfies all of the following:

- `summary.txt`, `charges.csv`, and either serial `rng_state.txt` or every MPI `rng_state_rankNNNNN.txt` exist.
- For schema v8 and later, `checkpoint_complete.txt` has `state=complete`, and its `batches` and `mpi_world_size`
  match `summary.txt`.
- `macro_residuals.csv` and `charge_ledger.csv` exist when declared by the manifest.
- The ordered-mesh fingerprint matches. Model or species fingerprint changes may continue with a warning.
- For an MPI resume, the saved `mpi_world_size` matches the current rank count.
- In `resume.toml`, `output.write_files=true` and `output.resume=true`. When set, `output.restart_from` names the
  checkpoint being inspected and `output.dir` is the intended new destination; when omitted, `output.dir` itself is
  the checkpoint source.
- The new `sim.batch_count` is at least the saved `batches`; use a larger target to advance by one or more batches.

**Safe action:**

- Use `restart_from` only as the checkpoint input and write new results to a separate `output.dir` to preserve the source.
- When resuming with changed physical or numerical settings, retain the warning, use a separate output directory, and check continuity across the change.
- Do not manually combine files from different generations or rewrite `checkpoint_complete.txt`.
- Do not use `resume=false` with the same `output.dir` as a workaround for a restart error.
- When periodic slots exist, BEACH automatically selects the newest complete candidate among the final output,
  `checkpoints/slot0`, and `slot1`.

After correcting the cause, follow [Resume once from a checkpoint](Execution.en.html#resume-once-from-a-checkpoint), write to a separate output directory,
and verify the recovery:

```bash
beachx lint resume.toml
beach resume.toml > resume.log 2>&1
grep '^resuming_from_batches=' resume.log
beachx inspect outputs/resumed
```

Recovery is complete when `resuming_from_batches` equals the selected checkpoint's `batches` and the new output's
`batches` equals `sim.batch_count` in `resume.toml`.
[Files used for resume](OutputReference.en.html#files-used-for-resume) is the source of truth for required files and automatic
selection.

## A surface model or periodic configuration is rejected

**Check:** First inspect input validation for the rejected combination.

```bash
beachx lint beach.toml
```

Only when lint succeeds but runtime rejects the combination, capture a runtime log:

```bash
beach beach.toml > beach-run.log 2>&1
run_status=$?
printf 'exit_code=%s\n' "$run_status"
tail -n 50 beach-run.log
```

**Expected:** The configuration satisfies the compatibility rules in [Input parameters](Parameters.en.html) and
[Field evaluation](FieldSolvers.en.html). Dielectric polarization is not implemented in the current version.

**Safe action:** Return to a supported combination of surface model, field solver, and field boundary.
Do not substitute an unimplemented physical model with another key or a looser numerical tolerance.

## FMM far correction is slow

**Check:** Inspect workload and phase timings on a compute node.

```bash
beachx workload beach.toml --threads 8
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach beach.toml
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

**Expected:** Repeated `cached_kneq0` runs reuse a compatible warm cache instead of rebuilding the cold operator every
time.

**Safe action:** Check the cache fingerprint and cache directory, then reuse a compatible warm cache.
Do not silently change the physical zero mode or far-correction model for speed. Repeat
[Validating Simulation Results](ValidationGuide.en.html) after such a model change.

## Report an unresolved problem

For a reproducible issue, provide:

- The configuration file and a minimal mesh.
- `beach --version`, compiler/MPI, and rank/thread counts.
- The command, exit code, and complete error.
- For restart problems, `summary.txt` and `checkpoint_complete.txt`; the RNG-state contents themselves are unnecessary.
- Reproduction steps with secrets removed from inputs and paths.
