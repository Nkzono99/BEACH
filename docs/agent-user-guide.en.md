title: BEACH Agent User Guide

Lang: [English](agent-user-guide.en.md) | [日本語](agent-user-guide.md)

# BEACH Agent User Guide

This operational guide helps AI agents start configuring, running, and verifying
BEACH. It links to canonical references instead of repeating full parameter and API
catalogs.

## Check first

BEACH combines boundary-element electric-field evaluation with particle tracking,
primarily to compute charge accumulation on insulator surfaces. The current baseline is:

- no solved outer-plasma field or outer-region particle transport
- ambient particles injected through nonperiodic box faces with `[particles.species.boundary_inflow]`
- photoelectrons closed with `photo_raycast`, species reflection on `inject_face`, and `neutral_return`
- absorbing mesh hits with surface-charge deltas committed once per batch
- execution through `sim.batch_count`, with `sim.max_step` as the per-particle limit
- `sim.tol_rel` recorded as a monitoring metric, not used for early stopping

[`SPEC.md`](../SPEC.md) is canonical for implementation behavior, and
[Input parameters](Parameters.en.html) is canonical for configuration keys. Removed
`[outer_plasma]`, `[coupling]`, and legacy selectors are rejected as unknown input.

## Shortest run

```bash
pip install beach-bem
beachx config init
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

For development from a repository checkout:

```bash
python -m pip install -e . --no-build-isolation
make check
make run CONFIG=examples/beach.toml
```

Passing `beachx lint` confirms the configuration contract, not convergence or physical
validity. Inspect outputs and convergence of the quantity of interest after the run.

## Choose a starting case

| Goal | Starting point | Read next |
|---|---|---|
| Beginner charging update | `examples/tutorial_insulator.toml` | [Ten-minute tutorial](Tutorial.en.html) |
| One-particle developer smoke | `examples/smoke_single_electron.toml` | [Configuration recipes](ConfigurationRecipes.en.html) |
| External plasma reservoir | `[particles.species.boundary_inflow]` | [Reservoir injection](ReservoirInjection.en.html) |
| Internal rectangular source | `source_mode="plane_source"` | [Choose where particles enter](ParticleSourcesBoundaries.en.html) |
| Closed photoelectron redistribution | `examples/periodic2_closed_photoelectron.toml` | [Photoelectron emission](PhotoelectronEmission.en.html) |
| Free-space field evaluation | `field_boundary.mode="free"` | [Field evaluation](FieldSolvers.en.html) |
| Two-axis-periodic field evaluation | `field_boundary.mode="periodic2"` | [Periodic electrostatics](PeriodicElectrostatics.en.html) |
| Output analysis | `beachx inspect` / `Beach(...)` | [Post-processing tutorial](PostprocessTutorial.en.html) |
| Resume from a checkpoint | `output.resume=true` | [Execution and resume](Execution.en.html) |

The boundary-reservoir plus closed-PE example is not a self-consistent outer-sheath
model. Converge box height, periodic images, `dt`, `max_step`, `batch_duration`,
particle count, and ray count against the quantity of interest.

## Edit a configuration

Use this four-stage workflow:

1. Copy the nearest `examples/*.toml`.
2. Select the required tables using [Configuration recipes](ConfigurationRecipes.en.html).
3. Check units, defaults, and exclusions in [Input parameters](Parameters.en.html).
4. Run `beachx lint beach.toml`, then check the Fortran parser diagnostics at runtime.

Check these constraints early:

| Feature | Required conditions |
|---|---|
| `boundary_inflow` | a `[domain]` box, resolved `sim.batch_duration>0`, a nonperiodic face, and an external VDF |
| `plane_source` | a `[domain]` box, positive `batch_duration`, internal rectangular `pos_low/high`, and `source_normal` |
| `photo_raycast` | a `[domain]` box, resolved `sim.batch_duration>0`, positive current density, and `rays_per_batch>=1` |
| closed PE | negative charge, `deposit_opposite_charge_on_emit=true`, `reflect` or `redistributed_reflect` on `inject_face` in `[particles.species.boundary]`, and `surface_charge_closure="neutral_return"` |
| `periodic2` | a `[domain]` box, exactly two entries in `domain.periodic_axes`, and a `[periodic2]` zero-mode policy |
| resume | `output.write_files=true`, required checkpoints, and the same MPI size |

The normal `periodic2` path uses FMM. Direct is limited to the validation split
reference `triangle_p0 + panel_spectral_reference + exclude_k0`. Distinguish the
finite-image approximation `field_periodic_far_correction="none"` from
`"cached_kneq0"`, which supplies the infinite-periodic nonzero modes. The canonical
compatibility table is in
[Field evaluation](FieldSolvers.en.html#solver-and-field-boundary-compatibility).

## Run and verify

Common development commands:

```bash
make check       # Fortran build check
make test-l0     # static/schema/build
make test        # L1: Python + quick Fortran
pytest -q        # Python only
ruff check .     # Python lint
```

Long FMM, MPI, and release checks are separated from the normal loop:

```bash
make test-l2
make test-l3
make test-heavy
make test-fortran-far-correction
make test-field-kernel-cache
make test-mpi
```

Run one Fortran target with
`FPM_ACTION=test ./build.sh --target <name>`. On KUDPC, classify the host first;
never run test payloads directly on a login node. Use `tssrun` or `sbatch` + `srun`.
See [Development workflow](Workflow.en.html).

Report test completion separately from physical validation. For a physics case,
check the quantity of interest against time step, spatial discretization, box or
periodic-image extent, particle statistics, and unresolved-particle fraction.

## Inspect outputs

Start with:

```bash
beachx inspect outputs/latest
```

The main files are `summary.txt`, `charges.csv`, `mesh_triangles.csv`,
`mesh_sources.csv`, `charge_ledger.csv`, and `rng_state.txt`. History, potential,
reservoir-residual, and performance-profile files are conditional. Their production
conditions and restart contract are canonical in the
[Output guide](OutputGuide.en.html).

The main Python entry point is:

```python
from beach import Beach

run = Beach("outputs/latest")
fig, ax = run.plot_mesh()
fig, ax = run.plot_charges(step=-1)
```

Use `beachx <command> --help` for complete CLI options and
[Python post-processing API](PythonPostprocessAPI.en.html) for the API reference.

## Change checklist

- Fortran behavior and `SPEC.md` agree.
- Public configuration changes update the schema, parser, examples, paired docs, and tests.
- Field, collision, boundary, injection, and resume changes include regression tests.
- `tol_rel` is not described as an early-stop condition.
- A completed run is not presented as numerical or physical validation.
- Paired Japanese and English pages preserve commands, identifiers, warnings, and ownership links.

## Documentation map

| Question | Canonical page |
|---|---|
| What does the simulator compute? | [`SPEC.md`](../SPEC.md), [Algorithms](Algorithms.en.html) |
| How should a case be assembled? | [Configuration recipes](ConfigurationRecipes.en.html) |
| What are the keys, types, units, and constraints? | [Input parameters](Parameters.en.html) |
| How do execution, workload estimation, and resume work? | [Execution and resume](Execution.en.html) |
| How should sources and boundaries be selected? | [Choose where particles enter](ParticleSourcesBoundaries.en.html) |
| How are collisions and box events handled? | [Particle events](ParticleEvents.en.html) |
| How should field solvers and periodic settings be selected? | [Field evaluation](FieldSolvers.en.html), [FMM](FMM.en.html) |
| How are outputs read and visualized? | [Output guide](OutputGuide.en.html), [Post-processing tutorial](PostprocessTutorial.en.html) |
| Which development tests should run? | [Development workflow](Workflow.en.html), [Validation guide](ValidationGuide.en.html) |
