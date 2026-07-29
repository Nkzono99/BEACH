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
| How many batches and seconds advanced, and how did particles terminate? | `summary.txt` | `batches`, `simulated_time_s`, `absorbed`, `escaped`, `survived_max_step` |
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
| Progress | `batches`, `simulated_time_s`, `last_rel_change` | accepted batches, sum of accepted physical widths, and final-batch charge-change monitor |
| Adaptive nonzero-mode progression | `periodic2_max_nonzero_mode_potential_step_V`, `adaptive_nonzero_mode_rejected_trials`, `adaptive_nonzero_mode_last_batch_duration_s`, `adaptive_nonzero_mode_last_potential_step_V`, `adaptive_nonzero_mode_omp_threads` | configured limit, cumulative rejected trials, last accepted width and measured $k\ne0$ potential change, and thread count required for restart |
| Field evaluation | `field_backend`, `field_source_model`, `field_kernel_id` | field solver and source kernel used for the output |
| Resolved external boundary | `coupling_update_mode`, `external_inflow_map`, `external_ordinary_open_model`, `external_interface_transport`, `outer_particle_mode_resolved` | update, inflow, ordinary-open handling, z-high transport, and timing resolved from the public configuration |

`absorbed` is an event count; it does not include charge sign or macro-particle weight. Read charge amounts from
`charge_ledger.csv` and the final distribution from `charges.csv`.

`last_rel_change` is a monitoring value. It is not an early-stop condition in the current implementation.

`periodic2_max_nonzero_mode_potential_step_V=0` means that adaptive
progression was disabled. When it is positive, `sim.batch_duration` is the
maximum trial width and `batches` counts only accepted batches.
`simulated_time_s` is the sum of accepted trial widths and therefore does not
generally equal `batches * sim.batch_duration`. A rejected trial rolls back the
RNG, macro-particle residuals, and outer/mean transaction, and does not enter
particle totals, `charge_history.csv`, `potential_history.csv`, or
`charge_ledger.csv`. `adaptive_nonzero_mode_last_potential_step_V` is the
maximum $k\ne0$ potential change measured over all panel centroids for the last
accepted trial; it is not a local-truncation-error estimate.
`adaptive_nonzero_mode_rejected_trials` is the total rejection count on the
shared ladder, including both $k\ne0$ bound failures and finite `implicit_mean`
Zhao ambient or interface field/barrier trust-region rejections. A measured
source-normalization change does not contract with trial width; it and other
$k=0$ closure failures stop without retry and are not included. Use the
execution log for the split.
To reproduce the reduction order and accepted ladder, an adaptive restart is
accepted only when every MPI rank has the same actual OpenMP team size and it
equals `adaptive_nonzero_mode_omp_threads`.
The value is `0` for a run without adaptive progression.

During the run, each `BEACH adaptive-kneq0 reject` standard-output record
reports the batch and trial width. A record containing `max_delta_phi_V`
denotes a $k\ne0$ bound rejection. A recoverable `implicit_mean` Zhao
trust-region rejection instead contains `implicit_status` and `reason`, so the
two causes can be counted separately. The
`BEACH adaptive-kneq0 accept` record reports the post-acceptance `time_end_s`,
trial width, measured potential change, and halving count.

The five external-boundary keys are not configuration inputs. They are a
receipt showing how the `[external_boundary]` facade resolved at runtime.

| Key | Output values |
| --- | --- |
| `coupling_update_mode` | `explicit` / `implicit_mean` |
| `external_inflow_map` | `source_vdf` / `infinity_barrier` / `kinetic_profile` |
| `external_ordinary_open_model` | `escape` / `potential_barrier` |
| `external_interface_transport` | `none` / `kinetic_1d` |
| `outer_particle_mode_resolved` | `local_source` / `same_batch` / `zhao_queue` |

For example, `particles.inflow_model="auto"` resolves to
`external_inflow_map=source_vdf` or `kinetic_profile`
according to the field and particle mode.
`ambient_linear_debye + same_batch` with negative `photo_raycast` automatically resolves to
`coupling_update_mode=implicit_mean`; do not write this internal value in public TOML.
Preserve both the input facade and this receipt when checking reproducibility.

During an implicit mean run, the `BEACH implicit-mean` progress record on standard output reports batch surface charge,
interface potential and field, species currents, additional measured surface current `J_other_A_m2`, zero-sum
`transaction_residual_C`, `mean_solver_iterations` for the scalar root, `sample_escape_fraction` for the ray sample, and
`return_weight_scale`. For the ambient linear-Debye closure, the former is the raw macro-charge sample fraction and the
latter is $R_{\mathrm{analytic}}/R_{\mathrm{sample}}$, not a probability, and can exceed one.
`mean_solver_iterations` is the scalar-root iteration count for this closure.

For the Zhao closure, the following `BEACH implicit-zhao` record reports the branch, virtual-cathode potential `phi_min_V`,
the complete-profile `barrier_J`, `source_scale` inferred from the measured interface current, marginal energy and its escape
fraction, nonlinear charge residual, recross charge fraction, and terminal-mismatch charge fraction. This path applies
measured-CDF weights directly to energy groups, so `return_weight_scale=1`. Under Zhao, the shared
`mean_solver_iterations` field counts connected candidate solves requested by the endpoint certificates, order-statistic
binary search, and marginal bisection; it is not the internal pseudo-arclength Newton-iteration count. Neither progress
record is itself a convergence criterion; inspect it together with `summary.txt`, `charge_ledger.csv`, and `charges.csv`.

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

`injected_count`, `emitted_count`, `absorbed_count`, `escaped_count`, and `discarded_unresolved_count` contain the
corresponding terminal-event counts. There is no independent count column for interface gross crossings.
`charge_ledger_residual_C` and
`charge_ledger_discarded_unresolved_abs_C` in `summary.txt` aggregate all species. See
[Check particle and charge balance](ValidationGuide.en.html#2-check-particle-and-charge-balance) for acceptance criteria.

For photoelectrons under `coupling_update_mode=implicit_mean`, interpret charge and count separately. Replacement applies
only to the component deferred after an outward z-high crossing. Under the ambient linear-Debye closure,
`escaped_to_infinity_C` and the post-return contribution to `absorbed_on_surface_C` follow the continuous-Maxwellian escape
and return totals from the scalar closure. Returned charge is normalized with `return_weight_scale` onto the actual-hit
distribution from one ray-sampling pass. Under the Zhao closure, the measured interface-energy CDF and nonlinear
$Q(\Phi_I)$ solution assign escape and return weights directly to each ray, with no common analytic scale.
Local reabsorption before the interface, escape through another open face such as z-low, and
`discarded_unresolved_C` retain their tracked values. Ledger terminal totals therefore combine these tracked
components by category: `escaped_to_infinity_C` combines tracked other-open-face escape with closure-derived z-high escape, while
`absorbed_on_surface_C` combines tracked local reabsorption with closure-derived post-return absorption.
`escaped_count` and the interface-return contribution to `absorbed_count` are terminal classifications of the ray sample,
not the source of truth for charge fractions. Terminal counts outside analytic replacement remain normal tracked-particle
updates.

`interface_outward_gross_C` accumulates signed macro charge for actual outward z-high crossings from the local region, while
`interface_returned_gross_C` accumulates signed macro charge for crossings back from the 1-D region. Ray crossings use weights
normalized to the terminal return and escape totals. If
$Q_{\mathrm{esc,z-high}}^{\mathrm{closure}}$ is the signed closure escape charge of the z-high deferred component, they
preserve `interface_returned_gross_C` = `interface_outward_gross_C` -
$Q_{\mathrm{esc,z-high}}^{\mathrm{closure}}$. Do not substitute the terminal total `escaped_to_infinity_C` into this
identity because that column can also include tracked escape through other open faces. A ray that recrosses z-high after
return can contribute again to these columns under the ambient linear-Debye path. The Zhao measured-CDF path treats such a
recrossing as an applicability violation and stops if it carries significant charge. In either case, gross charge is not a
terminal total.

For the transient Zhao queue, `charge_ledger_outer_flight_charge_before_C` and
`charge_ledger_outer_flight_charge_after_C` record outer-flight stock before and after each batch.

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
potential. Read this file to reuse the recorded values, or use
`compute_potential_mesh` / `compute_potential_points` to evaluate arbitrary
points with the native P0 panel kernel.

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

Accepted batch 1 is always included, followed by every `history_stride`
accepted batches. Potential history performs another field evaluation
when written, so increase `history_stride` to reduce output frequency for large meshes.

See [History Animation](PostprocessTutorial.en.html#history-animation) for plots and animations.

## Find Other Files by Purpose

| Purpose | File | Generation condition or role |
| --- | --- | --- |
| Break down runtime | `performance_profile.csv` | `BEACH_PROFILE=1` |
| Read the outer-sheath grid profile | `outer_plasma_profile.csv` | a `kinetic_1d` outer state is ready |
| Resume transient Zhao events | `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | `external_boundary.particles.mode="zhao_queue"`; former for serial, latter per MPI rank |
| Resume random-number state | `rng_state.txt` / `rng_state_rankNNNNN.txt` | former for serial, latter per MPI rank |
| Restore fractional macro-particles | `macro_residuals.csv` | residual state is allocated; one file even with MPI |

Every generated file also requires `output.write_files = true`.

## Locate Configuration-specific Values

Detailed acceptance criteria stay with each model page; this table only points to the output.

| Configuration | Main output | Details |
| --- | --- | --- |
| finite-image `periodic2` | periodic2 configuration in `summary.txt`, `charges.csv` | [Finite-image Configuration](FinitePeriodicConfiguration.en.html) |
| `cached_kneq0` | `periodic2_cache_hit`, `periodic2_operator_build_count`, `periodic2_cache_fingerprint`, `periodic2_cache_path` | [Periodic Far Correction](PeriodicFarCorrection.en.html) |
| adaptive $k\ne0$ progression | `simulated_time_s`, `adaptive_nonzero_mode_rejected_trials`, `adaptive_nonzero_mode_last_batch_duration_s`, `adaptive_nonzero_mode_last_potential_step_V` | [`batch_duration` stability and steady value](BatchDurationStability.en.html#adaptive-periodic2-nonzero-mode-progression) |
| `kinetic_1d` | `outer_plasma_profile.csv`, `interface_potential_V`, `gauss_residual_C`, `last_outer_update_batch` | [Standard 1-D Kinetic Outer Sheath](KineticOuterPlasma.en.html) |
| transient Zhao queue | `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv`, `outer_photoelectron_population_fraction`, `outer_photoelectron_column_per_area_m2`, `outer_photoelectron_column_target_per_area_m2`, `outer_photoelectron_column_residual_per_area_m2`, `outer_queue_event_count`, `outer_queue_signed_charge_C`, `outer_queue_fingerprint` | [Particle Escape and Return](ParticleEscapeReturn.en.html#queue-outer-flight-for-the-transient-zhao-closure) |
| outer particle transfer | `interface_outward_gross_C`, `interface_returned_gross_C`, `max_outer_flight_time_s`, `max_outer_frozen_field_ratio`, `max_outer_energy_relative_error` | [Particle Escape and Return](ParticleEscapeReturn.en.html) |
| `implicit_mean` return shadow | `implicit_mean_last_returned_outer_flight_time_mean_s`, `implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2` | [Standard 1-D Kinetic Outer Sheath](KineticOuterPlasma.en.html) |

`outer_infinity_potential_V` is an internal infinity-gauge diagnostic, not an input key.
It is fixed at zero for the current kinetic state.
`max_outer_energy_relative_error` is the maximum normalized conservation residual of normal kinetic plus electrostatic
energy in the kinetic 1-D return/escape mapping.
Because `implicit_mean` photoelectrons also use the individual 1-D profile map, they contribute to outer flight time,
frozen-field ratio, and energy-conservation residual. Recrossing after return contributes to the same diagnostics. These
photoelectrons are quasistationary shadows, so `max_outer_frozen_field_ratio` exceeding the configured limit does not by
itself fail the run. The ratio exposes the shadow-orbit timescale; it does not mean that delayed return current during UV
turn-on was resolved. Ordinary `same_batch` particles and ambient species continue to stop fail-closed on an over-limit ratio.

The two `implicit_mean_last_*` values are neither maxima nor cumulative values; they describe only the last batch completed
by the current invocation. A no-op resume that advances no batch omits both keys.
The first is the mean outer round-trip time of analytically weighted return excursions, weighted by positive charge
magnitude. The second is the positive photoelectron-column shadow estimate obtained from the same return excursions as
$\sum W_j\tau_j/(A\Delta t)$ by Little's law. It is not actual queue or ledger stock. The `escaped_to_infinity` outcome
itself contributes no residence time, but completed return excursions before eventual escape are included. When
`outer_integrated_charge_per_area_C_m2` is nonzero, dividing the shadow estimate by its absolute value gives
$\chi_{\mathrm{PE,shadow}}$, the relative scale of the omitted returning-photoelectron column and the integrated charge of
the 1-D outer profile. This is an interpretation ratio, not a built-in acceptance threshold.

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
| `outer_plasma_profile.csv` | required when resuming a ready `kinetic_1d` state |
| `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | required with the transient Zhao queue; former for serial and every rank file from the same MPI world size for MPI |

With `output.restart_from`, BEACH reads checkpoints from `restart_from` and writes new output under `output.dir`.
See [Run a Simulation](Execution.en.html#resume-a-run) for the resume procedure and fingerprint checks.
`summary.txt` also restores `simulated_time_s`, the cumulative rejected-trial
count, the last accepted width and potential change, and the adaptive OpenMP
team size. After resume,
`batch_count` therefore remains the cumulative target count of accepted
batches.
The queue checkpoint retains active phase-space records, terminal outcomes, due times, and `next_event_id`, and fails closed on
schema, rank, world-size, or completed-batch mismatches.

## Where to Go Next

- Decide whether the result is numerically and physically acceptable: [Validating Simulation Results](ValidationGuide.en.html)
- Create plots, animations, or Python analysis: [Post-processing Tutorial](PostprocessTutorial.en.html)
- Look up `[output]` settings: [Input Parameters Reference](Parameters.en.html#output-output-and-resume)
