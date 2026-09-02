title: Output Format Reference

# Output Format Reference

Lang: [English](OutputReference.en.md) | [日本語](OutputReference.md)

Use this reference to determine **which BEACH files are generated, which columns or receipts to read, and which files
are required for restart**. A receipt is a resolved setting or runtime state written to `summary.txt` as `key=value`,
not merely a copy of the input value.

For a first inspection, follow [Inspect Output Files](OutputGuide.en.html) before using this page. The machine-readable
source of truth for file-generation conditions is `schemas/beach.output-manifest.json`.

## File-generation conditions

Each row keeps a file's generation condition, format, purpose, and reader or restart role together. For a complex file,
the linked section repeats this information before expanding its columns and decision rules.

| File | Generation condition | Format and purpose | Python reader / restart |
| --- | --- | --- | --- |
| `summary.txt` | `output.write_files=true` | `key=value` receipts. [Completion and major keys](#locate-configuration-specific-values) | Fields of `FortranRunResult`. Required for restart |
| `charges.csv` | `output.write_files=true` | `elem_idx`, `charge_C`. Final signed charge per element | `result.charges`. Required for restart |
| `mesh_triangles.csv` | `output.write_files=true` | `elem_idx`, `v0x`, `v0y`, `v0z`, `v1x`, `v1y`, `v1z`, `v2x`, `v2y`, `v2z`, `charge_C`, `mesh_id` | `result.triangles`, `result.mesh_ids`. Not used for restart |
| `mesh_sources.csv` | OBJ or template mesh with `output.write_files=true` | `mesh_id`, `source_kind`, `template_kind`, `surface_model`, `epsilon_r`, `elem_count` | `result.mesh_sources`. Not used for restart |
| `mesh_potential.csv` | `output.write_files=true` and `output.write_mesh_potential=true` | `elem_idx`, `potential_V`. [Potential composition](#mesh-potential) | `result.mesh_potential_v`. Not used for restart |
| `charge_history.csv` | `output.write_files=true` and `output.history_stride>0` | `batch`, `processed_particles`, `rel_change`, `elem_idx`, `charge_C` | `result.history`. Not used for restart |
| `potential_history.csv` | `output.write_files=true`, `output.write_potential_history=true`, and `output.history_stride>0` | `batch`, `elem_idx`, `potential_V`. [Join to the reference](#history) | No dedicated `FortranRunResult` attribute; read the CSV directly |
| `top_reference_history.csv` | `output.write_files=true`, `output.write_potential_history=true`, `output.history_stride>0`, and a `[domain]` box | `batch`, `simulated_time_s`, `z_high_m`, `sample_n`, `potential_mean_V`, `potential_std_V`, `potential_min_V`, `potential_max_V` | No dedicated `FortranRunResult` attribute; read the CSV directly |
| `matching_plane_history.csv` | `output.write_files=true`, `output.history_stride>0`, and `surface_current_model.model=matching_plane_quasistatic` | [Seventeen-column accepted state](#matching_plane_quasistatic) | `result.matching_plane_history`. Not used for restart |
| `charge_ledger.csv` | `output.write_files=true` and charge-ledger state is present | [Twenty-five columns per species](#charge-ledger) | `result.charge_ledger`. Required for restart when summary has ledger metadata |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | `output.write_files=true`; the former for serial, the latter per MPI rank | Internal RNG state | No `FortranRunResult` attribute. Required for restart |
| `macro_residuals.csv` | `output.write_files=true` and macro-particle remainder state is allocated | `species_idx`, `face`, `residual` | No `FortranRunResult` attribute. Schema v8+: required when declared by the manifest |
| `checkpoint_complete.txt` | `output.write_files=true` for final output and each periodic slot | `key=value` completion manifest. [Meaning and selection](#files-used-for-resume) | No `FortranRunResult` attribute. Required for schema v8+ restart |
| `checkpoint_latest.txt` | `output.write_files=true` and after a matching accepted-batch commit with `output.checkpoint_stride>0` | Periodic-slot index | Candidate-selection hint, not evidence of completeness |
| `performance_profile.csv` | `BEACH_PROFILE=1` and `output.write_files=true` | After comment lines: `region`, `calls_sum`, `calls_mean`, `rank_min_s`, `rank_mean_s`, `rank_max_s`, `imbalance_ratio` | `beach-plot-performance-profile <output_dir>`. Not used for restart |

The Python reader requires `summary.txt` and `charges.csv`. Other files may be absent according to configuration or
runtime state.

## Locate configuration-specific values

`summary.txt` is a `key=value` file generated with `output.write_files=true`.

`load_fortran_result(...)` converts major values to `FortranRunResult` attributes; use
`beach.summary.load_summary_file(...)` for untyped receipts.

### Completion checks

Run completion and checkpoint completeness are separate conditions.

| Check | Condition | This condition does not establish |
| --- | --- | --- |
| Requested batches completed | `batches == sim.batch_count` using `summary.txt` and the input configuration | Numerical convergence or physical validity |
| Schema v8+ checkpoint is complete | `checkpoint_complete.txt` has `state=complete`, and the manifest agrees with every required file | Completion through `sim.batch_count` |
| Python reader succeeds | `summary.txt` and `charges.csv` satisfy the reader contract | Checkpoint completeness or physical validity |

`checkpoint_complete.txt` marks completion of a write transaction. An intermediate periodic checkpoint can also have
`state=complete`, so do not use it as evidence that the run reached its requested target. See
[Validating Simulation Results](ValidationGuide.en.html) for numerical and physical acceptance. For a current-schema
normal-completion output, require both of the first two conditions in the table.

### Common `summary.txt` receipts

The tables below list major receipts used to interpret results and decide restart compatibility. They do not claim that
every internal diagnostic key in `summary.txt` is a fixed API. The canonical producer is
[`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90); [`io.py`](../beach/fortran_results/io.py) owns the
keys converted to typed Python attributes.

| Receipt | Meaning |
| --- | --- |
| `mesh_nelem`, `mesh_count`, `mpi_world_size` | Mesh and MPI configuration |
| `batches`, `processed_particles`, `absorbed`, `escaped`, `escaped_boundary`, `survived_max_step` | Executed work and major particle outcomes |
| `simulated_time_s`, `last_rel_change` | Accepted simulated time and final relative change. `last_rel_change` is not an early-stop condition |
| `checkpoint_schema_version`, `checkpoint_stride` | Checkpoint format generation and resolved output interval |
| `model_fingerprint`, `mesh_fingerprint`, `species_fingerprint` | Configuration identities checked on restart |
| `reservoir_inflow_map`, `particle_ordinary_open_model` | Resolved reservoir and ordinary-open boundary models |

The build receipt consists of these five keys.

| Receipt | Meaning |
| --- | --- |
| `build_info_schema_version` | Build-receipt schema |
| `build_version` | BEACH version |
| `build_version_mode` | How the version was resolved |
| `build_source_commit` | Source commit |
| `build_id` | Build identifier |

Diagnose multiple box-boundary events with these five keys.

| Receipt | Meaning |
| --- | --- |
| `multiple_box_events_retry_attempted` | Number of particles for which retry was attempted |
| `multiple_box_events_retry_resolved` | Number of particles resolved by retry |
| `multiple_box_events_soft_discarded` | Number of soft-discarded particles |
| `multiple_box_events_soft_discard_fraction` | Soft-discard fraction of `processed_particles` |
| `multiple_box_events_soft_discarded_abs_charge_C` | Total absolute soft-discarded charge [C] |

### Field and boundary receipts

Field-reconstruction schema v2 provides the receipts needed to reconstruct the runtime electric field. For new output,
use `result.field_reconstruction` instead of inferring settings from a nearby `beach.toml`.

| Receipt | Meaning |
| --- | --- |
| `field_reconstruction_schema_version` | Field-reconstruction receipt schema; currently 2 |
| `field_reconstruction_resolved_field_solver` | Actual `direct`, `treecode`, or `fmm` selection |
| `field_reconstruction_fmm_expansion_order` | FMM expansion order. Schema v2 requires this and the preceding key |
| `field_reconstruction_field_bc_mode` | Resolved `free` or `periodic2` mode |
| `field_reconstruction_e0_V_m` | Resolved uniform external electric field `E0` |
| `field_reconstruction_use_box`, `field_reconstruction_box_min_m`, `field_reconstruction_box_max_m` | Box presence and bounds |
| `field_reconstruction_boundary_low`, `field_reconstruction_boundary_high` | Resolved states of the six faces |
| `field_reconstruction_tree_theta`, `field_reconstruction_tree_leaf_max` | Tree settings actually used |

The periodic-field keys are listed below. See [periodic2 Field Model](PeriodicElectrostatics.en.html) for their physical
meaning and solver selection.

| Receipt | Meaning |
| --- | --- |
| `field_reconstruction_periodic_image_layers` | Number of near-image shells |
| `field_reconstruction_periodic_far_correction` | Far-correction model |
| `field_reconstruction_periodic_nonzero_mode_backend` | Nonzero-mode backend |
| `field_reconstruction_periodic_zero_mode_policy` | Zero-mode policy |
| `field_reconstruction_periodic_lower_boundary_model` | Lower-boundary model |
| `field_reconstruction_periodic_reference_mode_layers` | Image layers for the reference nonzero mode |
| `field_reconstruction_periodic_panel_quadrature_order` | Panel quadrature order |
| `field_reconstruction_periodic_ewald_alpha` | Ewald splitting parameter |
| `field_reconstruction_periodic_ewald_layers` | Number of Ewald layers |
| `field_reconstruction_periodic_cache_dir` | Operator-cache directory |
| `field_reconstruction_periodic_generation_tolerance` | Operator-generation tolerance |

### `zhao_stationary`

See [Zhao stationary closure](ZhaoStationaryClosure.en.html) for the `[surface_current_model]` physical model,
Type A/B/C, and validity limits.
This section defines how to read its receipts.

| Receipt | Interpretation |
| --- | --- |
| `surface_current_model`, `surface_current_model_zhao_branch` | Selected model and branch |
| `surface_current_model_kinetic_contract` | Velocity-space boundary-map contract paired with the Zhao currents |
| `surface_current_model_photoelectron_active` | Whether PE channels exist. PE receipts are zero when false |
| `surface_current_model_reference_area_m2` | Reference area for current targets |
| `surface_current_model_phi0_V`, `surface_current_model_phi_m_V` | Resolved surface and minimum potentials |
| `surface_current_model_ambient_electron_density_m3` | Resolved ambient electron density |
| `surface_current_model_electron_current_density_A_m2`, `surface_current_model_ion_current_density_A_m2` | Signed electron and ion current densities |
| `surface_current_model_pe_emission_current_density_A_m2`, `surface_current_model_pe_escape_current_density_A_m2`, `surface_current_model_pe_return_current_density_A_m2` | Signed PE current densities. Emission and escape are positive; return is negative |
| `surface_current_model_pe_escape_particle_current_A` | Negative because it is the signed current carried by outward PE particles |
| `surface_current_model_net_current_density_A_m2` | Net surface current density over all channels |
| `surface_current_model_pe_budget_residual_current_density_A_m2` | PE emission-return-escape continuity residual |
| `surface_current_model_surface_budget_residual_current_density_A_m2` | Stationary surface-current residual |
| `surface_current_model_electron_inflow_reservoir_potential_V`, `surface_current_model_electron_inflow_access_potential_V`, `surface_current_model_electron_inflow_face` | Electron-inflow reservoir, access potential, and face |
| `surface_current_model_pe_outflow_barrier_potential_V`, `surface_current_model_pe_outflow_barrier_face` | PE outflow-barrier potential and face; zero when PE is inactive |

Face indices are `1..6 = x_low, x_high, y_low, y_high, z_low, z_high`.

### `matching_plane_quasistatic`

For this model, treat the accepted-batch matching-plane receipts listed below as state values instead of using static
surface-current targets.
See the [matching-plane numerical and response-table reference](MatchingPlaneReference.en.html) for the canonical
fixed-point equations, response CSV, and `implicit_zero_mode` contract.

| Item | `matching_plane_history.csv` contract |
| --- | --- |
| Generation | `output.write_files=true`, `output.history_stride>0`, and `surface_current_model.model=matching_plane_quasistatic` |
| One row | Leading `batch`, `simulated_time_s`, followed by the 15 state columns below; only accepted states are written |
| Interpretation | Read the fluxes, barriers, and iteration receipt on one row as one batch's fixed-point solution |
| Python | `result.matching_plane_history`; `result.matching_plane_state` holds the last accepted state |
| Restart | The CSV is not used. Schema v9 stores the last accepted state in `summary.txt` |

#### Accepted state

| `summary.txt` receipt | `matching_plane_history.csv` column | Meaning |
| --- | --- | --- |
| `matching_plane_displacement_C_m2` | `D_H_C_m2` | Displacement charge density $D_H$ at the matching plane |
| `matching_plane_phi_V` | `phi_H_V` | Matching-plane potential $\Phi_H$ |
| `matching_plane_electron_inward_flux_m2_s` | `electron_inward_flux_m2_s` | Inward electron flux |
| `matching_plane_ion_inward_flux_m2_s` | `ion_inward_flux_m2_s` | Inward ion flux |
| `matching_plane_electron_access_potential_V` | `electron_access_potential_V` | Electron access potential |
| `matching_plane_ion_access_potential_V` | `ion_access_potential_V` | Ion access potential |
| `matching_plane_photoelectron_barrier_potential_V` | `photoelectron_barrier_potential_V` | PE barrier potential |
| `matching_plane_photoelectron_outward_flux_m2_s` | `photoelectron_outward_flux_m2_s` | Outward PE flux fed back to the fixed point |
| `matching_plane_photoelectron_mean_normal_energy_eV` | `photoelectron_mean_normal_energy_eV` | PE mean normal energy fed back to the fixed point |
| `matching_plane_electron_outward_flux_m2_s` | `electron_outward_flux_m2_s` | Outward electron flux fed back to the fixed point |
| `matching_plane_ion_outward_flux_m2_s` | `ion_outward_flux_m2_s` | Outward ion flux fed back to the fixed point |
| `matching_plane_photoelectron_return_flux_m2_s` | `photoelectron_return_flux_m2_s` | PE return flux produced by the backend |
| `matching_plane_photoelectron_escape_flux_m2_s` | `photoelectron_escape_flux_m2_s` | PE escape flux produced by the backend |
| `matching_plane_iterations` | `iterations` | Fixed-point iteration count |
| `matching_plane_residual` | `residual` | Effective relative residual at acceptance |

The first two history columns are `batch` and `simulated_time_s`, giving 17 columns in total. When
`matching_plane_state_valid` is false, do not treat the state receipts listed above as an accepted
state.

Check the PE balance within the same batch:

$$
\mathtt{photoelectron\_outward\_flux\_m2\_s}
=\mathtt{photoelectron\_return\_flux\_m2\_s}
+\mathtt{photoelectron\_escape\_flux\_m2\_s}.
$$

The $D_H$ and $\Phi_H$ in each record correspond to the pre-commit surface-charge state used to track that batch.
In contrast, `simulated_time_s` is the time after the trial was accepted and advanced. Do not interpret the record as
the post-commit field at the start of the next batch.

#### Provenance and convergence

| Receipt | Meaning |
| --- | --- |
| `surface_current_model_response_backend` | `table` or `zhao_online` |
| `surface_current_model_zhao_root_selection` | Online-Zhao root policy: `require_unique`, `minimum_energy`, or `continuation` |
| `surface_current_model_implicit_zero_mode` | `T` advances the plane-mean $D_H$ with backward Euler; `F` holds its batch-start value fixed |
| `surface_current_model_matching_plane_z_m` | Matching-plane z coordinate |
| `surface_current_model_electron_species`, `surface_current_model_ion_species`, `surface_current_model_photoelectron_species` | Species assigned to each channel |
| `surface_current_model_coupling_rtol` | Relative tolerance for active components |
| `surface_current_model_coupling_atol` | Four componentwise absolute tolerances |
| `surface_current_model_coupling_max_iterations` | Maximum iteration count |
| `surface_current_model_coupling_relaxation` | Fixed-point relaxation factor |
| `surface_current_model_dynamic_state_source` | `surface_current_model_dynamic_state_source=accepted_batch_fixed_point` for committed matching trials. This compatibility string alone does not imply convergence |

The four values of `surface_current_model_coupling_atol` are ordered as follows:

1. Outward PE flux [m^-2 s^-1]
2. PE mean normal energy [eV]
3. Outward electron flux [m^-2 s^-1]
4. Outward ion flux [m^-2 s^-1]

All default to zero, and inactive components must also remain zero. Each active component uses
`max(coupling_rtol * backend_scale, coupling_atol)` as its threshold. An absolute-tolerance-dominated component is
converted to an effective residual, so a converged state's `matching_plane_residual` remains no greater than
`surface_current_model_coupling_rtol`. A state committed with an iteration-limit warning exceeds that value and has
`matching_plane_iterations` equal to the configured limit.

#### Backend-specific receipts

| Backend | Receipt |
| --- | --- |
| `table` | `surface_current_model_response_table_path`, `surface_current_model_response_content_fingerprint`. The fingerprint identifies the loaded response operator—$H$, canonical axes, and response values—not its path |
| `zhao_online` | `surface_current_model_response_contract=matching_plane_zhao_online_v1` |
| `zhao_online` | `surface_current_model_zhao_branch` |
| `zhao_online` | `surface_current_model_zhao_root_selection` |
| `zhao_online` | `surface_current_model_outer_solver=charge_driven_finite_h_sagdeev` |
| `zhao_online` | `surface_current_model_photoelectron_closure=moment_matched_half_maxwellian` |
| `zhao_online` | `surface_current_model_ambient_outward_feedback=transparent` |
| `zhao_online` with `require_unique` / `minimum_energy` | `surface_current_model_outer_solver_state=stateless` |
| `zhao_online` with `continuation` | `surface_current_model_outer_solver_state=accepted_endpoint_continuation_v1` |

`accepted_endpoint_continuation_v1` is a provenance receipt for the configured ownership rule: only an accepted endpoint
seeds the next batch. It is not an execution receipt proving that local Newton alone tracked the root, nor does it report
the number of full-multistart or step-subdivision evaluations or the root displacement. Use
`matching_plane_state_valid` to determine whether an accepted state exists.

## History

| File | Generation | Columns and interpretation | Python |
| --- | --- | --- | --- |
| `charge_history.csv` | `output.write_files=true` and `output.history_stride>0` | `batch`, `processed_particles`, `rel_change`, `elem_idx`, `charge_C`; one snapshot contains all element charges | `result.history` |
| `potential_history.csv` | `output.write_files=true`, `output.write_potential_history=true`, and `output.history_stride>0` | `batch`, `elem_idx`, `potential_V`; one snapshot contains every element-centroid potential | No dedicated attribute; read the CSV directly |
| `top_reference_history.csv` | `output.write_files=true`, `output.write_potential_history=true`, `output.history_stride>0`, and a `[domain]` box | `batch`, `simulated_time_s`, `z_high_m`, `sample_n`, `potential_mean_V`, `potential_std_V`, `potential_min_V`, `potential_max_V` | No dedicated attribute; read the CSV directly |
| `matching_plane_history.csv` | Matching-plane with stride enabled | [Seventeen-column accepted state](#matching_plane_quasistatic) | `result.matching_plane_history` |

The reference in `top_reference_history.csv` is the mean over the box z-high face, not infinity or plasma potential.
Join it to the same batch in `potential_history.csv` and compute `potential_V - potential_mean_V` for element-relative
potential.

For matching-plane output, typed `result.matching_plane_state` and `result.matching_plane_history` receipts avoid
manual column indexing.

## Mesh potential

| Item | `mesh_potential.csv` contract |
| --- | --- |
| Generation | `output.write_files=true` and `output.write_mesh_potential=true` |
| Columns | `elem_idx`, `potential_V` |
| Interpretation | Potential [V] at each element centroid |
| Python | `result.mesh_potential_v` |
| Restart | Not used |

| Field model | Terms included in the potential |
| --- | --- |
| Common | Analytic P0 triangle-panel self term |
| Finite periodic2 | Configured image shell |
| `cached_kneq0` | Cached nonzero mode plus the configured physical zero mode |

## Charge ledger

| Item | `charge_ledger.csv` contract |
| --- | --- |
| Generation | `output.write_files=true` and charge-ledger state is present |
| One row | Signed charges, corrections, weights, and counts for one species; all 25 columns appear below |
| Python | `ChargeLedgerEntry` objects in `result.charge_ledger` |
| Restart | Required when `summary.txt` contains ledger metadata |

The column order is:

| Column | Meaning |
| --- | --- |
| `batch` | Batch represented by the ledger state |
| `species_idx` | One-based species index |
| `injected_from_remote_C` | Signed charge injected from outside [C] |
| `emitted_from_surface_C` | Signed charge emitted from surfaces [C] |
| `absorbed_on_surface_C` | Signed charge absorbed by surfaces [C] |
| `escaped_to_infinity_C` | Raw trajectory escape charge [C] |
| `discarded_unresolved_C` | Signed charge removed as unresolved [C] |
| `neutral_return_correction_C` | Closed-PE neutral-return correction [C] |
| `neutral_return_weight_scale` | Weight scale used for neutral return |
| `neutral_return_unresolved_fraction` | Unresolved neutral-return fraction |
| `fixed_absorbed_target_charge_C` | Fixed-current absorption target [C] |
| `fixed_absorbed_weight_scale` | Weight scale applied to absorption samples |
| `fixed_emission_target_charge_C` | Fixed-current emission target [C] |
| `fixed_emission_weight_scale` | Weight scale applied to emission samples |
| `fixed_current_correction_C` | Charge added to the surface by the external closure [C] |
| `fixed_absorbed_applied_charge_C` | Compatibility alias equal to the absorption target |
| `fixed_emission_applied_charge_C` | Compatibility alias equal to the emission target |
| `fixed_escape_target_charge_C` | External escape target [C] |
| `fixed_escape_applied_charge_C` | Compatibility alias equal to the escape target |
| `fixed_escape_correction_C` | Difference between target and raw escape [C] |
| `injected_count` | Injected macro-particle count |
| `emitted_count` | Surface-emitted macro-particle count |
| `absorbed_count` | Surface-absorbed macro-particle count |
| `escaped_count` | Escaped macro-particle count |
| `discarded_unresolved_count` | Unresolved-discard macro-particle count |

### Companion receipts in `summary.txt`

| Receipt | Meaning |
| --- | --- |
| `charge_ledger_nspecies`, `charge_ledger_batch_count` | Species count and batch represented by the ledger |
| `charge_ledger_surface_charge_before_C`, `charge_ledger_surface_charge_after_C` | Surface-charge stock before and after the update |
| `charge_ledger_local_flight_charge_before_C`, `charge_ledger_local_flight_charge_after_C` | Local-flight stock before and after the update |
| `charge_ledger_unresolved_stock_before_C`, `charge_ledger_unresolved_stock_after_C` | Unresolved stock before and after the update |
| `charge_ledger_residual_C`, `charge_ledger_discarded_unresolved_abs_C` | Conservation residual and absolute unresolved charge |
| `charge_ledger_neutral_return_correction_C`, `charge_ledger_fixed_current_correction_C` | Total closure corrections |
| `charge_ledger_fixed_absorbed_applied_charge_C`, `charge_ledger_fixed_emission_applied_charge_C` | Total surface fixed-current targets |
| `charge_ledger_raw_escape_charge_C`, `charge_ledger_fixed_escape_target_charge_C`, `charge_ledger_fixed_escape_applied_charge_C`, `charge_ledger_fixed_escape_correction_C` | Raw escape, target, compatibility alias, and correction |
| `charge_ledger_fixed_applied_surface_net_charge_C` | Net charge applied to the surface by fixed current |
| `charge_ledger_fixed_pe_continuity_residual_C` | Fixed-PE continuity residual when PE is active |

### Decision rules

| Check | Read together | Caution |
| --- | --- | --- |
| External escape | `escaped_to_infinity_C`, `fixed_escape_target_charge_C`, `fixed_escape_correction_C` | The raw trajectory value is preserved. The escape correction is not deposited onto surface elements |
| Conservation residual | `charge_ledger_residual_C` and the surface / local-flight / unresolved before / after pairs above | Includes stocks, external fluxes, neutral return, and fixed-current corrections |
| `fixed_current` sampling variance | `absorbed_count`, `emitted_count`, raw charge, `fixed_absorbed_weight_scale`, `fixed_emission_weight_scale` | A count of one localizes the whole target onto the element represented by one sample |
| Closed PE | Raw absorption and unresolved values, `neutral_return_correction_C`, `neutral_return_unresolved_fraction` | Check the unresolved fraction separately even when total surface charge closes |

The output aliases for existing readers map to the equal in-memory targets as follows.

| Output alias | In-memory target |
| --- | --- |
| `fixed_absorbed_applied_charge_C` | `fixed_absorbed_target_charge_C` |
| `fixed_emission_applied_charge_C` | `fixed_emission_target_charge_C` |
| `fixed_escape_applied_charge_C` | `fixed_escape_target_charge_C` |

A small conservation residual does not establish elementwise sampling convergence. Test the charge distribution
separately across particle or ray counts, batch widths, and RNG seeds. See
[Check the charge ledger](OutputGuide.en.html#charge-ledger) for the basic procedure.

## Adaptive-batch diagnostics

With `periodic2.max_nonzero_mode_potential_step > 0`, `summary.txt` records:

| Receipt | Meaning |
| --- | --- |
| `simulated_time_s` | Simulated time through the accepted trial |
| `periodic2_max_nonzero_mode_potential_step_V` | Configured nonzero-mode potential-step limit |
| `adaptive_nonzero_mode_rejected_trials` | Cumulative rejected-trial count |
| `adaptive_nonzero_mode_last_batch_duration_s` | Duration of the last accepted batch |
| `adaptive_nonzero_mode_last_potential_step_V` | Nonzero-mode potential step of the last accepted batch |
| `adaptive_nonzero_mode_omp_threads` | OpenMP thread count fixed for adaptive restart |

Rejected trials do not appear in history files or the charge ledger.

## Files used for resume

With `output.checkpoint_stride > 0`, BEACH updates this structure after an accepted-batch commit. Even with
`output.checkpoint_stride=0`, a normal completion writes the final checkpoint.

Here, “complete” means that the **transaction containing the restart files is complete**. Separately establish
completion of the requested run with `batches == sim.batch_count` using `summary.txt` and the input configuration.

```text
<output.dir>/
├── checkpoint_complete.txt
├── checkpoint_latest.txt
└── checkpoints/
    ├── slot0/
    └── slot1/
```

### Contents of each checkpoint

| File | Required when | Format and role |
| --- | --- | --- |
| `summary.txt` | Always | `key=value` statistics, schema, fingerprints, and ledger stocks |
| `charges.csv` | Always | Element charges as `elem_idx`, `charge_C` |
| `rng_state.txt` | Serial | RNG array length and state |
| `rng_state_rankNNNNN.txt` | MPI; one for every rank | Per-rank RNG array length and state |
| `macro_residuals.csv` | Declared by a schema v8+ manifest | `species_idx`, `face`, `residual` global macro-particle remainders |
| `charge_ledger.csv` | Summary / manifest contains ledger metadata | [Twenty-five-column species ledger](#charge-ledger) |
| `checkpoint_complete.txt` | Schema v8+ | `schema_version`, `state`, `batches`, `mpi_world_size`, `macro_residuals_present`, `charge_ledger_present` |

### Update and automatic selection

| Stage | Contract |
| --- | --- |
| Update | Finish the inactive slot before atomically switching `checkpoint_latest.txt`. At most two generations are retained |
| Search | Inspect the final output under `output.restart_from` and both `slot0` / `slot1` directories |
| Selection | Choose the loadable candidate with all required files and the largest `batches` value |
| Index failure | Recover a slot with a complete manifest even if `checkpoint_latest.txt` is missing, malformed, or stale |
| Mismatch | Stop when required files, the mesh fingerprint, saved-state array shapes, or MPI world size differ. Warn and continue on model or species fingerprint changes |

`output.restart_from` changes only the read source; new output is written to `output.dir`.

### Schema generations

| Schema | Contract |
| --- | --- |
| v6 | `macro_residuals.csv` uses `species_idx,face,residual`; `face=0` is the legacy source and `1..6` are boundary faces. The older two-column `species_idx,residual` form remains readable |
| v8+ | Publish `state=in_progress` before replacing state. Atomically publish `state=complete` only after closing summary, charges, every rank's RNG, and manifest-declared residual / ledger files |
| v9 | Store accepted matching-plane feedback, potential, return / escape fluxes, and iteration receipts in `summary.txt` |

For schema v8 and later, `checkpoint_complete.txt` itself is required. If the final-output directory contains files
from different generations, BEACH rejects it and falls back to a complete periodic slot.

The model fingerprint includes canonical response-table contents for the table backend, or the Zhao contract, branch
policy, and implicit mode for the online backend. Model and species fingerprints remain provenance checks, but a
mismatch alone does not stop resume. BEACH warns and carries the saved charge and statistics into the current
configuration. This is a changed-condition continuation, so use a separate output directory and assess continuity of
the observables. A mesh fingerprint mismatch remains fatal because BEACH cannot safely map `charges.csv` rows to
different elements.

Online implicit bracket and step-subdivision nodes are neither a persistent table nor checkpoint state. With
`continuation`, BEACH tries to reconstruct the restart seed from the accepted response saved in `summary.txt` and falls
back to the minimum-energy bootstrap if reconstruction fails. Every policy evaluates the needed response directly from
the same Zhao contract, so no `matching_query.csv` or fingerprint of automatic search points is required.

See [Execution and Resume](Execution.en.html) for the actual procedure.

## Reading from Python

The [file-generation table](#file-generation-conditions) gives each file's attribute or explicitly states when no
dedicated attribute exists. Read representative typed attributes as follows.

```python
from beach import load_fortran_result

result = load_fortran_result("outputs/tutorial")
print(result.charges)
print(result.charge_ledger)
print(result.field_reconstruction)
print(result.matching_plane_state)
```

| Output | Python attribute |
| --- | --- |
| `charges.csv` | `result.charges` |
| `charge_ledger.csv` | `result.charge_ledger` |
| Field-reconstruction receipt | `result.field_reconstruction` |
| Matching-plane summary | `result.matching_plane_state` |
| `matching_plane_history.csv` | `result.matching_plane_history` |

`summary.txt` and `charges.csv` are required. Mesh, history, ledger, field-reconstruction, and matching-plane state are
loaded only when the corresponding output is present.

When `matching_plane_state_valid` is false, the reader does not
expose `matching_plane_state` or `matching_plane_history` as an accepted state.

For another case, replace `outputs/tutorial` with its actual `output.dir`. See the
[Python Post-processing API](PythonPostprocessAPI.en.html) for class and attribute details, and
[Validating Simulation Results](ValidationGuide.en.html) for physical and numerical acceptance.
