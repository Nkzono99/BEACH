title: Inspect Output Files

# Inspect Output Files

Lang: [English](OutputGuide.en.md) | [日本語](OutputGuide.md)

BEACH outputs are grouped into final state, history, and restart state. The
machine-readable source of truth for file-generation conditions is
`schemas/beach.output-manifest.json`.

## Files to inspect first

| Goal | File |
| --- | --- |
| Check completion and processed batches | `summary.txt` |
| Read final element charges | `charges.csv` |
| Read triangle geometry and mesh IDs | `mesh_triangles.csv` |
| Map mesh IDs to input meshes | `mesh_sources.csv` |
| Inspect per-species charge accounting | `charge_ledger.csv` |
| Inspect performance breakdown | `performance_profile.csv` |

`summary.txt` records statistics, resolved configuration, build information,
checkpoint schema, and model / mesh / species fingerprints. The resolved
boundary reservoir and ordinary-open state is reported as
`reservoir_inflow_map` and `particle_ordinary_open_model`.
The `field_reconstruction_*` receipt records the resolved `E0`, box, boundary,
periodic nonzero/zero-mode model, tree settings, actual field solver, and FMM
expansion order. Schema v2 requires
`field_reconstruction_resolved_field_solver` and
`field_reconstruction_fmm_expansion_order`. Python field reconstruction uses
this receipt for new outputs rather than guessing from a nearby `beach.toml`.

The `[surface_current_model]` receipt is recorded as `surface_current_model`. For `zhao_stationary`,
`surface_current_model_*` fields report the selected branch, reference area, $\phi_0$, $\phi_m$, resolved ambient-electron
density, and signed electron / ion / PE-emission / PE-escape / PE-return / net current densities. In particular,
`surface_current_model_pe_return_current_density_A_m2` is negative, while emission and escape are positive.
`surface_current_model_pe_escape_particle_current_A` is negative because it is the signed current carried outward by
PE particles. The two budget residuals independently check PE emission-return-escape continuity and the stationary
surface-current balance.
`surface_current_model_kinetic_contract` and the inflow-access / outflow-barrier potential and face receipts record the
velocity-space boundary map paired with the Zhao currents. Face indices are
`1..6 = x_low, x_high, y_low, y_high, z_low, z_high`.

## History

| File | Condition | Main columns |
| --- | --- | --- |
| `charge_history.csv` | `history_stride > 0` | batch, element, charge |
| `potential_history.csv` | `write_potential_history=true` | batch, element, potential |
| `top_reference_history.csv` | Above plus a `[domain]` box | batch, time, z-high potential statistics |

The reference in `top_reference_history.csv` is the mean over the box z-high
face, not infinity or plasma potential. Join it to the same batch in
`potential_history.csv` and compute `potential_V - potential_mean_V` for
element-relative potential.

## Mesh potential

With `write_mesh_potential=true`, BEACH writes `mesh_potential.csv`.

- Potential [V] at each element centroid
- Analytic P0 triangle-panel self term
- Configured image shell for finite periodic2
- Cached nonzero mode plus the configured physical zero mode for `cached_kneq0`

## Charge ledger

`charge_ledger.csv` records the following for each species:

- Signed charge for injection, surface emission, surface absorption, escape,
  and unresolved discard
- Counts for each terminal outcome
- Closed-PE `neutral_return_correction_C`
- `neutral_return_weight_scale`
- `neutral_return_unresolved_fraction`
- `fixed_absorbed_target_charge_C` and `fixed_absorbed_weight_scale`
- `fixed_emission_target_charge_C` and `fixed_emission_weight_scale`
- `fixed_current_correction_C` added to surface charge by the external closure
- applied absorption/emission values in `fixed_absorbed_applied_charge_C` / `fixed_emission_applied_charge_C`
- external-escape values in `fixed_escape_target_charge_C` / `fixed_escape_applied_charge_C`
- the difference from raw escape in `fixed_escape_correction_C`

`escaped_to_infinity_C` remains the raw trajectory result. Zhao fixed-current operation preserves it and stores the
external escape target alongside it. Because the escape correction is not deposited onto surface elements, analysis
can keep the surface update separate from the external-boundary budget.

`charge_ledger_residual_C` in `summary.txt` combines changes in surface,
local-flight, and unresolved stocks with external fluxes, the neutral-return
correction, and the fixed-current correction.

For `fixed_current` statistical diagnostics, read the channel's `absorbed_count` or `emitted_count` together with its raw
charge and `fixed_*_weight_scale`. A count of one localizes the entire target onto the element represented by that one
sample. A small conservation residual cannot detect this sampling variance; separately test convergence of the elementwise
charge distribution across particle or ray counts, batch widths, and RNG seeds.

For closed PE, BEACH preserves raw absorption and unresolved values. The
correction and scale are separate, so the unresolved fraction remains
independently inspectable even when total surface charge closes.

## Adaptive-batch diagnostics

With `periodic2.max_nonzero_mode_potential_step > 0`, `summary.txt` records:

- `simulated_time_s`
- `periodic2_max_nonzero_mode_potential_step_V`
- `adaptive_nonzero_mode_rejected_trials`
- `adaptive_nonzero_mode_last_batch_duration_s`
- `adaptive_nonzero_mode_last_potential_step_V`
- `adaptive_nonzero_mode_omp_threads`

Rejected trials do not appear in history files or the charge ledger.

## Files used for resume

With `output.checkpoint_stride > 0`, BEACH updates this structure after an accepted-batch commit:

```text
outputs/latest/
├── checkpoint_complete.txt
├── checkpoint_latest.txt
└── checkpoints/
    ├── slot0/
    └── slot1/
```

Each slot contains the restart files in the table below. BEACH finishes the inactive slot before atomically switching
`checkpoint_latest.txt`, so at most two generations are retained.

| File | Role |
| --- | --- |
| `summary.txt` | Statistics, schema, fingerprints, and ledger stocks |
| `charges.csv` | Element charges |
| `rng_state.txt` | Serial RNG state |
| `rng_state_rankNNNNN.txt` | Per-rank MPI RNG state |
| `macro_residuals.csv` | Global macro-particle residuals by species and face; required when declared by a schema v8+ manifest |
| `charge_ledger.csv` | Restored when summary contains ledger metadata |
| `checkpoint_complete.txt` | Completion manifest published last for schema v8 and later checkpoints |

`output.restart_from` changes only the checkpoint read source. New output is
written to `output.dir`. BEACH compares the final output and both periodic slots, then selects the loadable complete
checkpoint with the largest `batches` value. A slot with a complete manifest remains recoverable when
`checkpoint_latest.txt` is missing, malformed, or stale. BEACH stops instead of falling back to a new run when a
required file is missing or when fingerprints, mesh size, species count, or MPI
world size differ.

For schema v8 and later, BEACH changes `checkpoint_complete.txt` to `in_progress` before replacing restart state. It
atomically publishes `complete` only after closing summary, charges, every rank's RNG state, and the residual and ledger
declared by the manifest. `checkpoint_complete.txt` itself is a required restart file for these schemas. If final-output
writing is interrupted and leaves files from different generations, BEACH
rejects that directory and falls back to the complete periodic slot.

Checkpoint schema v6 writes `macro_residuals.csv` as `species_idx,face,residual`. `face=0` denotes the legacy source and
`1..6` denote boundary faces. The older two-column `species_idx,residual` form remains readable.

## Reading from Python

```python
from beach import FortranResults

result = FortranResults("outputs/latest")
print(result.summary)
print(result.charges)
```

See [Python Post-processing API](PythonPostprocessAPI.en.md) for details.
For physical and numerical acceptance, see
[Validating Simulation Results](ValidationGuide.en.html).
