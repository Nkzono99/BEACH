title: Inspect Output Files

Lang: [English](OutputGuide.en.md) | [日本語](OutputGuide.md)

# Inspect Output Files

Use this page to locate a value in a BEACH output directory. Use [Validating Simulation Results](ValidationGuide.en.html)
to decide whether a run is acceptable, and use the [Post-processing Tutorial](PostprocessTutorial.en.html) to create plots
or animations.

## Start with Four Outputs

```bash
beachx inspect outputs/latest
```

If `output.dir` was changed, replace `outputs/latest` with the actual output directory.

| Question | Read | Main contents |
| --- | --- | --- |
| How many batches ran, and how did particles terminate? | `summary.txt` | `batches`, `absorbed`, `escaped`, `survived_max_step` |
| Where did each species' charge enter and leave? | `charge_ledger.csv` | injection, emission, absorption, escape, unresolved discard |
| Where is the final surface charge? | `charges.csv` + `mesh_triangles.csv` | element index, charge, triangle coordinates, `mesh_id` |
| Is batch-by-batch evolution available? | `charge_history.csv`, `potential_history.csv` | element charge and centroid-potential history |

`beachx inspect` prints a run summary. For completion, charge conservation, time convergence, and model-applicability
criteria, follow [Validating Simulation Results](ValidationGuide.en.html).

## When Files Are Updated

This page assumes `output.write_files = true`. With `false`, BEACH does not create or update `summary.txt`, `charges.csv`,
or other files under `output.dir`. Any existing files in that directory belong to an earlier run.

The machine-readable [`beach.output-manifest.json`](../schemas/beach.output-manifest.json) is canonical for generation
conditions and restart roles.

## `summary.txt`: Whole-run Summary

`summary.txt` uses `key=value` lines. Start with these groups.

| Group | Main keys | Meaning |
| --- | --- | --- |
| Run identity | `build_version`, `build_source_commit`, `model_fingerprint`, `mesh_fingerprint`, `species_fingerprint` | build, configuration, mesh, and species that produced the output |
| Size | `mesh_nelem`, `mesh_count`, `mpi_world_size` | element, mesh, and MPI-rank counts |
| Particle processing | `processed_particles`, `absorbed`, `escaped`, `escaped_boundary`, `survived_max_step`, `multiple_box_events_soft_discarded`, `multiple_box_events_soft_discarded_abs_charge_C` | particle-event totals and recorded soft discards |
| Progress | `batches`, `last_rel_change` | completed batches and final-batch charge-change monitor |
| Field evaluation | `field_backend`, `field_source_model`, `field_kernel_id` | field solver and source kernel used for the output |

`absorbed` is an event count; it does not include charge sign or macro-particle weight. Read charge amounts from
`charge_ledger.csv` and the final distribution from `charges.csv`.

`last_rel_change` is a monitoring value. It is not an early-stop condition in the current implementation.

## `charge_ledger.csv`: Charge Transfer by Species

`charge_ledger.csv` has one row per species and separates signed charge [C] from event counts.

| Column | Recorded charge |
| --- | --- |
| `injected_from_remote_C` | charge entering through `volume_seed` or `reservoir_face` |
| `emitted_from_surface_C` | tracked charge emitted from a surface through `photo_raycast` |
| `absorbed_on_surface_C` | charge absorbed by the mesh |
| `escaped_to_infinity_C` | charge escaping through an open boundary or outer model |
| `discarded_unresolved_C` | charge discarded alive at `max_step` |
| `interface_outward_gross_C` | gross charge transferred from the local to outer region |
| `interface_returned_gross_C` | gross charge returned from the outer to local region |

The corresponding `*_count` columns contain event counts. `charge_ledger_residual_C` and
`charge_ledger_discarded_unresolved_abs_C` in `summary.txt` aggregate all species. See
[Check particle and charge balance](ValidationGuide.en.html#2-check-particle-and-charge-balance) for acceptance criteria.

## Match Surface Values to the Mesh

| File | One row represents | Main columns |
| --- | --- | --- |
| `charges.csv` | one triangle element | `elem_idx`, `charge_C` |
| `mesh_triangles.csv` | one triangle element | `elem_idx`, three vertices, `charge_C`, `mesh_id` |
| `mesh_sources.csv` | one mesh | `mesh_id`, `source_kind`, `surface_model`, element count |
| `mesh_potential.csv` | one triangle element | `elem_idx`, `potential_V` |

Join `charges.csv` to `mesh_triangles.csv` by `elem_idx`. For multiple meshes, use `mesh_id` in `mesh_triangles.csv`
and the corresponding row in `mesh_sources.csv`.

`mesh_potential.csv` is generated only with `output.write_mesh_potential = true` and contains final element-centroid
potential. For `field_source_model="triangle_p0"`, use this file instead of reconstructing potential as point charges in Python.

## History Files

Set a positive `output.history_stride` before the run to write history.

```toml
[output]
history_stride = 1
write_potential_history = true
```

| File | Generation condition | Contents |
| --- | --- | --- |
| `charge_history.csv` | `output.history_stride > 0` | element charge at recorded batches |
| `potential_history.csv` | above plus `output.write_potential_history = true` | element-centroid potential at recorded batches |

Batch 1 is always included, followed by every `history_stride` batches. Potential history performs another field evaluation
when written, so increase `history_stride` to reduce output frequency for large meshes.

See [History Animation](PostprocessTutorial.en.html#history-animation) for plots and animations.

## Find Other Files by Purpose

| Purpose | File | Generation condition or role |
| --- | --- | --- |
| Break down runtime | `performance_profile.csv` | `BEACH_PROFILE=1` |
| Read the outer-sheath grid profile | `outer_plasma_profile.csv` | a `kinetic_1d` / `unified_linear_response` outer state is ready |
| Read the photoelectron energy distribution | `photoelectron_histogram.csv` | `outer_plasma.photoelectron_histogram_enabled=true` and the state is ready |
| Resume random-number state | `rng_state.txt` / `rng_state_rankNNNNN.txt` | former for serial, latter per MPI rank |
| Restore fractional macro-particles | `macro_residuals.csv` | residual state is allocated; one file even with MPI |

Every generated file also requires `output.write_files = true`.

## Locate Configuration-specific Values

Detailed acceptance criteria stay with each model page; this table only points to the output.

| Configuration | Main output | Details |
| --- | --- | --- |
| finite-image `periodic2` | periodic2 configuration in `summary.txt`, `charges.csv` | [Finite-image Configuration](FinitePeriodicConfiguration.en.html) |
| `cached_kneq0` | `periodic2_cache_hit`, `periodic2_operator_build_count`, `periodic2_cache_fingerprint`, `periodic2_cache_path` | [Periodic Far Correction](PeriodicFarCorrection.en.html) |
| `kinetic_1d` | `outer_plasma_profile.csv`, `interface_potential_V`, `gauss_residual_C`, `last_outer_update_batch` | [Standard 1-D Kinetic Outer Sheath](KineticOuterPlasma.en.html) |
| `unified_linear_response` | above plus `outer_accessible_fraction_min`, `outer_accessible_fraction_max`, `outer_accessible_fraction_refinement_error` | [Advanced Rough-surface Linear Screening](UnifiedLinearResponse.en.html) |
| photoelectron histogram | `photoelectron_histogram.csv`, `photoelectron_previous_charge_ratio`, `photoelectron_linear_applicability_status` | [Photoelectron Emission and Lifecycle](PhotoelectronEmission.en.html) |
| outer particle transfer | `interface_outward_gross_C`, `interface_returned_gross_C`, `max_outer_flight_time_s`, `max_outer_frozen_field_ratio`, `max_outer_energy_relative_error` | [Particle Escape and Return](ParticleEscapeReturn.en.html) |

When a mesh contains `dielectric` elements, `summary.txt` records `surface_model_dielectric_elem_count` and
`surface_model_note=metadata_only_dielectric_present`. In the current implementation, `dielectric` is metadata; this note
does not mean that a dielectric boundary condition was solved.

## Files Used for Resume

Resume state shares the output directory with analysis files, but should not be edited individually.

| Checkpoint | Resume role |
| --- | --- |
| `summary.txt`, `charges.csv` | always required |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | former for serial; every rank file for MPI |
| `macro_residuals.csv` | restores global fractional state when present |
| `charge_ledger.csv` | required when `summary.txt` contains ledger metadata |
| `outer_plasma_profile.csv` | required when resuming a ready `kinetic_1d` / `unified_linear_response` state |
| `photoelectron_histogram.csv` | required when resuming photoelectron histogram state |

With `output.restart_from`, BEACH reads checkpoints from `restart_from` and writes new output under `output.dir`.
See [Run a Simulation](Execution.en.html#resume-a-run) for the resume procedure and fingerprint checks.

## Where to Go Next

- Decide whether the result is numerically and physically acceptable: [Validating Simulation Results](ValidationGuide.en.html)
- Create plots, animations, or Python analysis: [Post-processing Tutorial](PostprocessTutorial.en.html)
- Look up `[output]` settings: [Input Parameters Reference](Parameters.en.html#output-output-and-resume)
