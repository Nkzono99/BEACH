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

`charge_ledger_residual_C` in `summary.txt` combines changes in surface,
local-flight, and unresolved stocks with external fluxes and the neutral-return
correction.

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
| `macro_residuals.csv` | Global macro-particle residuals distinguished by species and face, restored when present |
| `charge_ledger.csv` | Restored when summary contains ledger metadata |

`output.restart_from` changes only the checkpoint read source. New output is
written to `output.dir`. When the source contains `checkpoint_latest.txt`, BEACH automatically selects the complete
checkpoint with the largest `batches` value from the final output and periodic slot. BEACH stops instead of falling back to a new run when a
required file is missing or when fingerprints, mesh size, species count, or MPI
world size differ.

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
