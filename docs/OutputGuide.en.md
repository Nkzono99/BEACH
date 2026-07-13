title: Reading Output Files

Lang: [English](OutputGuide.en.md) | [日本語](OutputGuide.md)

# Reading Output Files

After a first run, inspect completion, charge balance, element charge, and history under `outputs/latest/` in that order.
For plotting commands, see [Post-processing Tutorial](PostprocessTutorial.en.html). For all configuration keys, see [Input Parameters Reference](Parameters.en.html).

## First Checks

First confirm that the run reached the configured batch and wrote the required files. After completion is established, evaluate
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
| `outer_plasma_profile.csv` | `kinetic_1d` / `unified_linear_response` | Converged outer-grid `z, phi, E, rho`; also used as the restart profile |
| `charges.csv` | Always | final charge per element |
| `mesh_triangles.csv` | Always | triangle vertices, element IDs, `mesh_id` |
| `mesh_sources.csv` | Template mesh runs | mapping from `mesh_id` to template kind / surface model / element count |
| `mesh_potential.csv` | `output.write_mesh_potential=true` | final centroid potential |
| `charge_history.csv` | `output.history_stride > 0` | element charge history by batch |
| `potential_history.csv` | `write_potential_history=true` and `history_stride > 0` | centroid potential history by batch |
| `performance_profile.csv` | `BEACH_PROFILE=1` | phase timing |
| `rng_state*.txt` | Always | random-number state for resume |
| `macro_residuals*.csv` | Reservoir-style injection | fractional macro-particle state for resume |
| `charge_ledger.csv` | Always | per-species signed-charge flux and particle counts |

## Interpreting Success and Warnings

Start with these quantities in `summary.txt`.

`field_source_model` and `field_kernel_id` identify the element kernel that produced the output. `triangle_p0_exact_p2m_near` denotes the all-vertex topology, analytic panel-near, exact-panel-P2M FMM. `FieldKernel.from_result` dispatches `triangle_p0` output to the panel C ABI. The other Python potential/field/force/field-line estimators remain point-only and fail closed.

For split periodic2 runs, `summary.txt` records interface potential and normal field, `eta_phi_kneq0`, `eta_field_kneq0`, `eta_gap`, `eta_local_charge`, the Gauss residual, integrated outer charge, and the last outer-update batch. These values are part of the applicability and restart contract; a split checkpoint missing its outer state is rejected.

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

| Item | Meaning |
| --- | --- |
| `batches` | In a normal run, completion means this reaches `sim.batch_count` |
| `absorbed` | Number of particles absorbed by surfaces; the main indicator for charge accumulation |
| `escaped` | Number of particles leaving through open boundaries; useful for checking injection and boundary conditions |
| `survived_max_step` | Particles that remained alive until `sim.max_step`; if large, revisit `dt`, the box, or injection conditions |
| `last_rel_change` | Monitoring value for the final batch charge change; it is not an early-stop condition in the current implementation |
| `charge_ledger_residual_C` | Transactional charge-conservation residual; unresolved discard must still be checked separately |
| `charge_ledger_discarded_unresolved_abs_C` | Non-cancelling absolute max-step discard charge across species |

When meshes include `conductor` or `dielectric`, `summary.txt` may include notes.
`dielectric` is metadata in the current implementation; it is not a solved dielectric boundary model.

## Using History

Set `output.history_stride` to a positive value to inspect time evolution.

```toml
[output]
history_stride = 1
write_potential_history = true
```

Potential history requires additional field evaluations, so it can increase runtime for large meshes or frequent history output.
Start with `history_stride = 1` on a small case, then thin the history for larger runs.

## Useful Next Commands

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif
```

From Python:

```python
from beach import Beach

b = Beach("outputs/latest")
print(b.result.absorbed, b.result.escaped)
b.plot_mesh()
b.plot_potential()
```

## Resume Outputs

With `output.resume=true`, `summary.txt`, `charges.csv`, `rng_state*.txt`, `macro_residuals*.csv`, and `charge_ledger.csv` when enabled are used as checkpoint files. Schema v3 restores matching model, mesh, and species fingerprints plus the complete outer solver profile/state. A schema-v2 three-column outer profile remains readable, but it forces a new outer solve at the next refresh instead of being treated as a complete held state.
When `output.restart_from` is set, checkpoint files are read from `restart_from`, while new outputs are written under `output.dir`.
