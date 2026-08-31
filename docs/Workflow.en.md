title: Development Workflow

Lang: [English](Workflow.en.md) | [日本語](Workflow.md)

# Development Workflow

This contributor guide explains how to modify BEACH and select tests for the affected area.
For installation, normal case execution, OpenMP or MPI runs, workload estimation, and resume operations,
see [Installation](Installation.en.html) and [Run a Simulation](Execution.en.html).
Before reading the source, use the [Developer Architecture Overview](Architecture.en.html) to identify the
execution flow and the owner of each state.

## Create a development environment

**Prerequisites:** Install Python, `make`, a Fortran compiler, and `fpm`.

```bash
make --version
gfortran --version
fpm --version
python --version
```

**Action:** From the checkout root, install the editable package and run the lightweight build check.

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
make check
```

**Expected result:** The Python package refers to the checkout and the Fortran source compiles.
`make check` uses `BEACH_VERSION_MODE=dev`, allowing fpm to reuse incremental compilation when only the Git hash changes.

Use Make targets through `build.sh` in normal development. Run one Fortran test target with this form:

```bash
FPM_ACTION=test ./build.sh --target test_particle_stepper
```

Do not run multiple `fpm test` commands or `build.sh` test targets concurrently. They share the same `build/`
directory and can corrupt or race on compilation artifacts.

## Select tests from the change

Run the directly related tests first, then complete at least the gate shown in the table. Combine rows when a change
crosses subsystem boundaries. Direct tests alone are not sufficient when a change affects a public contract or numerical result.

| Changed area | Direct tests | Minimum gate / additional check |
| --- | --- | --- |
| Documentation, navigation, or Japanese/English parity | `pytest -q tests/python/test_docs_sync.py tests/python/test_documentation_contracts.py` | `make test-l1`; after a site-structure change, run `python tools/sync_starlight_docs.py` and then `npm --prefix docs-site run check` |
| TOML, schema, parser, or defaults | `test_app_config_parser`, `test_physics_config_types`, `tests/python/test_config_schema.py`, `tests/python/test_config_cli.py` | `make test-l1` and `make schema-check` |
| Mesh templates, OBJ import, or panel geometry | `test_templates_importers_runtime`, `test_panel_geometry_near`, `test_panel_kernel` | `make test-l1`; use `make test-l3` when panel FMM is affected |
| Particle sources, reservoirs, or photoelectron injection | `test_injection_sampling`, `test_reservoir_injection`, `test_external_field_velocity_grid` | `make test-l1` |
| Boris update, collision, box boundary, or particle events | `test_particle_stepper`, `test_boundary`, `test_dynamics_basic` | `make test-l1`; also use `make test-mpi` when the MPI or OpenMP path changes |
| Surface model, charge closure, or charge ledger | `test_surface_models`, `test_surface_current_model`, `test_charge_ledger`, `test_simulator` | `make test-l1`; include `test_matching_plane_simulator` for matching-plane work |
| Field snapshot, Direct, or Treecode | `test_electrostatic_snapshot`, `test_dynamics_field_solver`, `test_panel_kernel` | `make test-l1` |
| FMM, periodic2, zero mode, or nonzero mode | Relevant `test_coulomb_fmm_*`, `test_periodic_*`, and `test_dynamics_fmm` targets | `make test-l3`; use `make test-fortran-far-correction` for the far correction |
| C ABI or native field kernel | `test_field_kernel_c`, `test_periodic_zero_mode_c` | `make test-l2`; use `make test-field-kernel-cache` for the cache receipt |
| Output, checkpoint, fingerprint, or restart | `test_output_writer_io`, `test_output_writer_potential`, `test_restart`, `test_model_fingerprint` | `make test-l1` and the relevant Python reader tests |
| Python reader, analysis, or CLI | The corresponding `tests/python/test_*.py` files | `make test-l1` |

`fpm.toml` is canonical for Fortran test names, and `Makefile` is canonical for tier membership. The
[Developer Architecture Overview](Architecture.en.html#move-from-a-subsystem-to-its-implementation-and-tests)
maps subsystems to source, direct tests, and related documentation.

## Select a test tier

```bash
make test-l0      # static / schema / build
make test         # L1: Python + lightweight Fortran; alias of test-l1
make test-l2      # L1 + C / kernel contracts
make test-l3      # L2 + heavy FMM / panel tests
```

- L0 checks `git diff --check`, source text, JSON schemas, and `make check`.
- L1 adds the complete Python suite and the normal Fortran test targets.
- L2 adds the C ABI and periodic zero-mode C contracts.
- L3 adds `test_dynamics_fmm`, FMM core, panel near-correction, and other heavy targets.

The following gates are not all included in the normal tiers. Run them explicitly when required by the change or release decision.

```bash
make test-heavy
make test-fortran-far-correction
make test-field-kernel-cache
make test-mpi
make test-fortran-benchmark
make test-physics-release
make test-full
```

`test-field-kernel-cache` checks the native cache and plane-oracle receipt. Keep performance comparisons separate from
correctness tests and use the release-profile `make test-fortran-benchmark`. Follow
[Physics Release Verification](PhysicsReleaseVerification.en.html) for release decisions and convergence fixtures.

Before handoff, run `make fmt-check-fortran` after a Fortran change, `ruff check .` after a Python change, and finish with
`git diff --check` for whitespace errors. A test tier does not replace these format and lint checks when they apply.

## Run development tests on KUDPC

Before running a test payload on KUDPC, inspect `hostname`, `module list`, and, when available, `spartition` and `qgroup`.
Use them to determine the host role, active `Sys*` module, and allocation.

On login nodes (`camphor*`, `laurel*`, `cinnamon*`, and `gardenia*`), limit work to editing, diff inspection,
short log reads, `make check`, job submission, and monitoring. Do not run `make test*`, `pytest`, `fpm test`,
MPI or OpenMP tests, or benchmarks directly on a login node.

- Obtain a compute-node allocation with `tssrun` for a short development test.
- Run tier tests, MPI tests, and heavy or release gates with `srun` inside an `sbatch` job.
- Do not run multiple `fpm test` commands concurrently inside a job.
- Keep the checkout, inputs, cache, and logs on `/home`, `/LARGE0`, `/LARGE1`, or an agreed `/FAST` path visible
  from compute nodes. Do not depend on a login-node `/TMP` path.

Production simulation job design, thread and rank layout, and workload estimation belong in
[Run a Simulation](Execution.en.html).

## Change a public contract

When a public contract changes, update its canonical description and consumers in the same change as the implementation.

| Changed contract | Update and verify together |
| --- | --- |
| Fortran simulation behavior | `SPEC.md`, the relevant model or numerical-method page, regression tests, and convergence tests when needed |
| TOML table, key, type, default, or constraint | `schemas/beach.schema.json` and its distributed copies, Fortran parser, Python validator, `Parameters.md` / `.en.md`, and examples |
| CLI or Python / Fortran API | Help, public signatures, examples, Japanese and English API documentation, and consumer tests |
| Output file, column, or generation condition | `schemas/beach.output-manifest.json`, Fortran writer, Python reader, `OutputGuide.md` / `.en.md`, and restart compatibility |
| Checkpoint state or fingerprint | Schema version and compatibility rule, writer, loader, restart tests, and the restart contract in `SPEC.md` |
| Numerical method or physical model | Scope, fail-closed conditions, Direct or oracle comparison, numerical convergence, and `ValidationGuide.md` / `.en.md` |

For every documentation change, keep commands, identifiers, constraints, warnings, and expected results equivalent in
Japanese and English. In the current implementation, `sim.tol_rel` is a monitoring and output value, not an early-stop
condition. Report execution success, test completion, numerical convergence, and physical validity separately.

See [Canonical ownership](Architecture.en.html#distinguish-canonical-ownership) for the responsibility of each source,
and [Validate Simulation Results](ValidationGuide.en.html) for validation of a user case.
