title: BEACH Documentation
ordered_subpage: OutputGuide.en.md
ordered_subpage: ConfigurationRecipes.en.md
ordered_subpage: Parameters.en.md
ordered_subpage: Configuration.en.md
ordered_subpage: PostprocessTutorial.en.md
ordered_subpage: PythonPostprocessAPI.en.md
ordered_subpage: Algorithms.en.md
ordered_subpage: FieldSolvers.en.md
ordered_subpage: ParticleChargeLoop.en.md
ordered_subpage: FMMCore.en.md
ordered_subpage: BatchDurationStability.en.md
ordered_subpage: Workflow.en.md
ordered_subpage: FortranDependencyMap.en.md
ordered_subpage: agent-user-guide.en.md

Lang: [English](index.en.md) | [日本語](index.md)

# BEACH Documentation

BEACH (BEM + Accumulated CHarge) is a surface-charging simulator that couples Coulomb fields from charges on triangular boundary elements with test-particle tracking in batches.
The current release focuses on charge accumulation on insulator surfaces. The Fortran runtime `beach` performs the simulation, while the Python CLI/API (`beachx` and the `beach` package) handles configuration checks, post-processing, and visualization.

## Three-Minute Quick Start

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem

mkdir run_periodic2
cd run_periodic2

beachx config init
beachx lint beach.toml
beachx config render beach.toml --output beach.rendered.toml

beach beach.rendered.toml
beachx inspect outputs/latest
```

`beach` reads a TOML file expanded to final runtime keys. When using high-level notation, run `beachx config render` before simulation.
Because `render` overwrites the input file when no output is specified, first-time users should prefer `--output beach.rendered.toml`.

## Success Check

A normally completed run writes `outputs/latest/summary.txt`, and `batches` reaches the configured `sim.batch_count`.
Then run `beachx inspect outputs/latest` to check absorbed and escaped particles, element charges, and mesh metadata.
See [Reading Output Files](OutputGuide.en.html) for the meaning of each output file.

## What to Read Next

| Goal | Page |
| --- | --- |
| Check outputs after the first run | [Reading Output Files](OutputGuide.en.html) |
| Adapt a configuration from examples | [Configuration Recipes](ConfigurationRecipes.en.html) |
| Look up every `beach.toml` key | [Input Parameters Reference](Parameters.en.html) |
| Understand `beachx config render` | [`beachx config` / High-Level Notation Guide](Configuration.en.html) |
| Make the first plots | [Post-processing Tutorial](PostprocessTutorial.en.html) |
| Use the full Python API | [Python Post-processing API Reference](PythonPostprocessAPI.en.html) |
| Understand the numerical model | [BEACH Algorithm Overview](Algorithms.en.html) |

## Use Cases

| Goal | Entry point |
| --- | --- |
| Run a small template mesh | [Configuration Recipes](ConfigurationRecipes.en.html), "Minimal plane-mesh run" |
| Choose a particle source model | [Configuration Recipes](ConfigurationRecipes.en.html) and `particles` in [Input Parameters Reference](Parameters.en.html) |
| Use two-periodic-axis boundaries | `periodic2` in [Configuration Recipes](ConfigurationRecipes.en.html) and [Field Solvers and Boundary Conditions](FieldSolvers.en.html) |
| Tune `batch_duration` | [`batch_duration` Stability](BatchDurationStability.en.html) |
| Inspect implementation APIs | [Fortran API](https://nkzono99.github.io/BEACH/fortran/) and [Fortran Dependency Map](FortranDependencyMap.en.html) |

## Documentation Index

- [Reading Output Files](OutputGuide.en.html)
- [Configuration Recipes](ConfigurationRecipes.en.html)
- [Input Parameters Reference](Parameters.en.html)
- [`beachx config` / High-Level Notation Guide](Configuration.en.html)
- [Post-processing Tutorial](PostprocessTutorial.en.html)
- [Python Post-processing API Reference](PythonPostprocessAPI.en.html)
- [BEACH Algorithm Overview](Algorithms.en.html)
- [Field Solvers and Boundary Conditions](FieldSolvers.en.html)
- [Particle Tracking and Charge Accumulation](ParticleChargeLoop.en.html)
- [Coulomb FMM Core Details](FMMCore.en.html)
- [`batch_duration` Stability](BatchDurationStability.en.html)
- [Execution and Development Workflow](Workflow.en.html)
- [Fortran Dependency Map](FortranDependencyMap.en.html)
- [BEACH Agent User Guide](agent-user-guide.en.html)

For API-level implementation details, use the FORD-generated [Fortran API](https://nkzono99.github.io/BEACH/fortran/).
