title: BEACH Documentation
ordered_subpage: Installation.en.md
ordered_subpage: Tutorial.en.md
ordered_subpage: OutputGuide.en.md
ordered_subpage: Troubleshooting.en.md
ordered_subpage: ConfigurationRecipes.en.md
ordered_subpage: Configuration.en.md
ordered_subpage: PostprocessTutorial.en.md
ordered_subpage: ValidationGuide.en.md
ordered_subpage: Parameters.en.md
ordered_subpage: PythonPostprocessAPI.en.md
ordered_subpage: Algorithms.en.md
ordered_subpage: FieldSolvers.en.md
ordered_subpage: ParticleChargeLoop.en.md
ordered_subpage: FMMCore.en.md
ordered_subpage: BatchDurationStability.en.md
ordered_subpage: PhysicsReleaseVerification.en.md
ordered_subpage: PhysicsRedesignCompletionAudit.en.md
ordered_subpage: Workflow.en.md
ordered_subpage: FortranDependencyMap.en.md
ordered_subpage: agent-user-guide.en.md

Lang: [English](index.en.md) | [日本語](index.md)

# BEACH

BEACH simulates test-particle trajectories and charge accumulation on insulating surfaces in batches.

## First use

1. [Check requirements and install](Installation.en.html)
2. [Run the official beginner case](Tutorial.en.html)
3. [Read the output](OutputGuide.en.html)
4. [Validate the physical and numerical result](ValidationGuide.en.html)

```bash
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

## Command and data flow

```text
beach.toml
    ├─ beachx lint ── configuration checks
    ▼
beach ────────────── Fortran simulation
    ▼
outputs/<case>/
    ├─ beachx inspect / animate ── inspection and plots
    └─ Python package beach ───── custom analysis
```

## Supported scope

| Feature | Status | Note |
| --- | --- | --- |
| Insulator charge accumulation | Supported | Primary scope |
| Floating conductors | Conditional | Check field-boundary and surface-model compatibility |
| Dielectric polarization | Not implemented | `epsilon_r` is metadata, not an independent polarization boundary condition |
| Two-axis periodic fields | Supported | Distinguish finite images, cached infinite operators, and zero modes |
| Outer plasma | Conditional | Requires model applicability and numerical error contracts |
| Automatic convergence stop | Not implemented | `tol_rel` is a monitoring metric |

Unsupported combinations fail closed instead of silently selecting another model.

## Entry points

| Goal | Page |
| --- | --- |
| Build a case | [Configuration Recipes](ConfigurationRecipes.en.html) |
| Plot and analyze output | [Post-processing Tutorial](PostprocessTutorial.en.html) |
| Look up input keys | [Input Parameters](Parameters.en.html) |
| Inspect numerical methods | [Algorithm Overview](Algorithms.en.html) |
| Resolve a problem | [Troubleshooting](Troubleshooting.en.html) |

Use the generated [Fortran API](https://nkzono99.github.io/BEACH/fortran/) for procedure-level details.
