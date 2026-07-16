title: Reading Output Files

Lang: [English](OutputGuide.en.md) | [日本語](OutputGuide.md)

# Reading Output Files

After a first run, inspect completion, charge balance, element charge, and history under `outputs/latest/` in that order.
For plotting commands, see [Post-processing Tutorial](PostprocessTutorial.en.html). For all configuration keys, see [Input Parameters Reference](Parameters.en.html).

## First Checks

The checks on this page assume `output.write_files = true`. With `false`, BEACH does not create or update `summary.txt`,
`charges.csv`, or the other files under `output.dir` for that run, so `beachx inspect` cannot inspect it. Existing files,
if any, belong to an earlier run. Only the summary printed to the terminal remains available.

First confirm that the run reached the configured batch count and wrote the required files. After completion is established, evaluate
the physical model and discretization with [Validating Simulation Results](ValidationGuide.en.html).

1. `outputs/latest/summary.txt` exists.
2. `batches` in `summary.txt` reaches the configured `sim.batch_count`.
3. `charges.csv` exists and contains final element charges.
4. `beachx inspect outputs/latest` prints a run summary without errors.

```bash
beachx inspect outputs/latest
```

If `output.dir` is changed, replace `outputs/latest` with that output directory.

## Main Files

| File | When it appears | What to inspect first |
| --- | --- | --- |
| `summary.txt` | Always | batch count, absorbed and escaped counts, last relative charge change, MPI rank count |
| `outer_plasma_profile.csv` | An outer state is ready for `kinetic_1d` / `unified_linear_response` | Converged outer-grid `z, phi, E, rho`; a conditional checkpoint |
| `photoelectron_histogram.csv` | `outer_plasma.photoelectron_closure="individual_return"` and the histogram state is ready | Previous-batch and cumulative photoelectron histogram; a conditional checkpoint |
| `charges.csv` | Always | final charge per element |
| `mesh_triangles.csv` | Always | triangle vertices, element IDs, `mesh_id` |
| `mesh_sources.csv` | OBJ or template mesh | `mesh_id`, `source_kind`, `surface_model`, and element count |
| `mesh_potential.csv` | `output.write_mesh_potential=true` | final centroid potential |
| `charge_history.csv` | `output.history_stride > 0` | element charge history by batch |
| `potential_history.csv` | `output.write_potential_history=true` and `output.history_stride > 0` | centroid potential history by batch |
| `performance_profile.csv` | `BEACH_PROFILE=1` | phase timing |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | Always | Serial or MPI rank-local random-number state for resume |
| `macro_residuals.csv` | Macro-particle residual state is allocated | One MPI-global fractional macro-particle state file |
| `charge_ledger.csv` | Charge-ledger state is present | per-species signed-charge flux and particle counts |

Every file in this table also requires `output.write_files = true`. The machine-readable
[`beach.output-manifest.json`](../schemas/beach.output-manifest.json) is canonical for additional file-generation conditions and restart roles.

## Species-resolved charge balance

`charge_ledger.csv` records charge and particle-count inflow, internal transfer, and outflow for each species.

| Item | Meaning |
| --- | --- |
| `injected_from_remote_C` | Charge entering from `volume_seed` or `reservoir_face` |
| `emitted_from_surface_C` | Tracked charge leaving a surface through `photo_raycast` |
| `absorbed_on_surface_C` | Charge absorbed by the mesh |
| `escaped_to_infinity_C` | Charge escaping through open or outer models |
| `discarded_unresolved_C` | Charge discarded alive at `max_step` |
| `interface_outward_gross_C` / `interface_returned_gross_C` | Gross outward and return charge transfer across the outer ownership interface |

The charge-conservation residual combines before/after changes of surface, local-flight, outer-flight, and unresolved stocks with
remote injection, escape, and discard. Surface emission and absorption are internal transfers between surface and flight stocks,
so they are not counted twice as independent external sources.

Separately from `charge_ledger_residual_C`, where opposite species can cancel,
`charge_ledger_discarded_unresolved_abs_C` reports $\sum_s|Q_{s,\mathrm{discard}}|$. A small conservation residual does not
validate a run with large unresolved discard.

## Interpreting Success and Warnings

Start with these quantities in `summary.txt`.

| Item | Meaning |
| --- | --- |
| `batches` | In a normal run, completion means this reaches `sim.batch_count` |
| `absorbed` | Number of absorbed macro-particle events; it does not include charge sign or particle weight |
| `escaped` | Number of particles leaving through open boundaries; useful for checking injection and boundary conditions |
| `survived_max_step` | Particles that remained alive until `sim.max_step`; if large, revisit `dt`, the box, or injection conditions |
| `last_rel_change` | Monitoring value for the final batch charge change; it is not an early-stop condition in the current implementation |
| `charge_ledger_residual_C` | Charge-conservation residual; unresolved discard must still be checked separately |
| `charge_ledger_discarded_unresolved_abs_C` | Non-cancelling absolute max-step discard charge across species |

`absorbed` alone does not measure accumulated charge. Use `absorbed_on_surface_C` for signed, weighted absorbed charge,
the difference between `charge_ledger_surface_charge_before_C` and `charge_ledger_surface_charge_after_C` for net
surface-charge change, and `charges.csv` for the final distribution.

### Model-Specific Diagnostics

`field_source_model` and `field_kernel_id` identify the element kernel that produced the output. `triangle_p0_exact_tree_near` denotes Treecode with all-vertex node radii, analytic panel-near evaluation, and monopole far evaluation. `triangle_p0_exact_p2m_near` denotes the all-vertex topology, analytic panel-near, exact-panel-P2M FMM. `FieldKernel.from_result` dispatches `triangle_p0` output to the panel C ABI. The other Python potential/field/force/field-line estimators remain point-only and fail closed.

For split periodic2 runs, `summary.txt` records `interface_potential_V`, `interface_normal_field_V_m`,
`interface_eta_phi_kneq0`, `interface_eta_field_kneq0`, `interface_eta_gap`, `interface_eta_local_charge`,
`gauss_residual_C`, `outer_integrated_charge_C`, and `last_outer_update_batch`. These values cover the interface
potential and normal field, the Gauss residual, integrated outer charge, and update point. They diagnose physical-model
applicability and are part of the restart contract; a split checkpoint missing its outer state is rejected.

For `unified_linear_response`, also inspect `outer_accessible_fraction_min`, `outer_accessible_fraction_max`, and
`outer_accessible_fraction_refinement_error`. The latter is the maximum change in accessible fraction when the
surface-height sample count is doubled along each periodic axis. Initialization stops if this value exceeds
`outer_plasma.accessible_fraction_tolerance`.

For `cached_kneq0`, inspect these cache diagnostics.

| Summary key | Cold miss | Warm reuse |
| --- | --- | --- |
| `periodic2_cache_hit` | `F` | `T` |
| `periodic2_operator_build_count` | `1` on the root rank | `0` |
| `periodic2_cache_fingerprint` | generated identity | matches the reused identity |
| `periodic2_cache_path` | publication target | read source |

A cold build distributes target slices across MPI ranks and proxy right-hand sides across OpenMP threads within each rank.
Only the root rank performs cache I/O. The regenerable operator payload is not stored in the checkpoint.
The configured cache directory and generation tolerance are also recorded in the summary.

With particle transfer enabled, `interface_outward_gross_C` and `interface_returned_gross_C` in `charge_ledger.csv` record gross crossings and are not added twice to the conservation residual. `max_outer_flight_time_s`, `max_outer_frozen_field_ratio`, and `max_outer_energy_relative_error` in `summary.txt` are MPI-global run maxima.

When a mesh includes `dielectric`, `summary.txt` records `surface_model_dielectric_elem_count` and
`surface_model_note=metadata_only_dielectric_present`. `dielectric` is metadata in the current implementation; it is
not a solved dielectric boundary model. A mesh containing only `conductor` does not emit this note.

## Using History

Set `output.history_stride` to a positive value to inspect time evolution.

```toml
[output]
history_stride = 1
write_potential_history = true
```

`charge_history.csv` is written when

$$
(\texttt{stats.batches}-1)\bmod\texttt{history\_stride}=0,
$$

so batch 1 is always included. With `output.write_potential_history=true`, the current `q_elem` refreshes the field at the same stride
and writes element-centroid `potential_history.csv`.

Potential history requires additional field evaluations, so it can increase runtime for large meshes or frequent history output.
Start with `history_stride = 1` on a small case, then thin the history for larger runs.

## Useful Next Commands

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh.png

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif
```

Only results with `field_source_model="point"` can reconstruct and plot potential through the Python estimator.

```bash
beachx inspect outputs/latest \
  --save-potential-mesh outputs/latest/potential_mesh.png
```

For `triangle_p0`, `--recompute-potential`, `--save-potential-mesh`, and `--show` fail closed. To obtain
element-centroid potential, set
`output.write_mesh_potential = true` before the run and use the simulator-written `mesh_potential.csv`.

From Python:

```python
from beach import Beach

b = Beach("outputs/latest")
print(b.result.absorbed, b.result.escaped)
b.plot_mesh()

if b.result.field_source_model == "point":
    b.plot_potential()
```

`Beach.plot_potential()` is also point-source only.

## Resume Outputs

With `output.resume=true`, files are loaded according to the configuration and saved state:

| Checkpoint file | Resume behavior |
| --- | --- |
| `summary.txt`, `charges.csv` | Always required |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | The former is required for serial runs; every rank file is required for MPI |
| `macro_residuals.csv` | Restores the global residual when present. MPI still uses one file; legacy rank-local files are rejected |
| `charge_ledger.csv` | Required when the summary contains ledger checkpoint metadata |
| `outer_plasma_profile.csv` | Required when resuming a ready outer state for `kinetic_1d` / `unified_linear_response` |
| `photoelectron_histogram.csv` | Required when resuming `outer_plasma.photoelectron_closure="individual_return"` |

Schema v3 restores matching model, mesh, and species fingerprints plus the complete outer solver profile/state. A schema-v2 three-column outer profile remains readable, but it forces a new outer solve at the next refresh instead of being treated as a complete held state.
When `output.restart_from` is set, checkpoint files are read from `restart_from`, while new outputs are written under `output.dir`.
