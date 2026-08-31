title: Inspect Output Files

# Inspect Output Files

Lang: [English](OutputGuide.en.md) | [日本語](OutputGuide.md)

> **Question:** Did the official tutorial finish successfully, and is its surface charge consistent with
> charge conservation?
>
> **One-sentence answer:** Start with `beachx inspect` for run counts, then compare `charges.csv` with
> `charge_ledger.csv` for final surface charge and particle outcomes.

After reading this page, you should be able to establish completion of the official case and then inspect final charge,
meshes, per-species charge accounting, history, and potential. Model-specific receipts, complete column definitions,
and checkpoint schemas are separated into the [Output Format Reference](OutputReference.en.html).

## Check the official tutorial result

From the directory where you completed the [10-minute tutorial](Tutorial.en.html), inspect the path configured as
`output.dir` in `beach.toml`. The official tutorial uses `outputs/tutorial`.

```bash
beachx inspect outputs/tutorial
```

When the tutorial is run with one OpenMP thread as instructed, the output includes at least these lines:

```text
directory=outputs/tutorial
processed_particles=4000
absorbed=3720 escaped=280
batches=20 last_rel_change=...
charge_sum=-1.192019e-10
```

The minimum pass conditions are:

- Both `beach` and `beachx inspect` exit with code `0`.
- `batches=20` agrees with `sim.batch_count=20` in `beach.toml`.
- `processed_particles=4000` agrees with 200 particles per batch × 20 batches.
- `outputs/tutorial/summary.txt` and `outputs/tutorial/charges.csv` exist.

`absorbed=3720`, `escaped=280`, and `charge_sum=-1.192019e-10` are reference values for the current version with
`rng_seed=12345` and one OpenMP thread. If the thread count or RNG implementation differs, do not assume that
random-sequence-dependent values must equal the reference.

This check establishes run completion only. `last_rel_change` and `tol_rel` are not automatic stop conditions.
Assess numerical convergence and physical validity separately with
[Validating Simulation Results](ValidationGuide.en.html). If the command fails or a minimum condition is not met,
continue with [Troubleshooting](Troubleshooting.en.html).

## Files to inspect first

| Question | First place to inspect |
| --- | --- |
| How many particles and batches completed, and what was the final change? | `summary.txt` and `beachx inspect` |
| What is the final charge on each element? | `charges.csv` |
| What are the triangle geometry and mesh IDs? | `mesh_triangles.csv` |
| Which input mesh produced each mesh ID? | `mesh_sources.csv` |

### `summary.txt`

`summary.txt` records run statistics and resolved settings as `key=value` lines. Start with `processed_particles`,
`absorbed`, `escaped`, `batches`, and `last_rel_change`. When reading model-specific keys, use them as the runtime
receipt rather than guessing from a nearby `beach.toml`; see
[Configuration-specific output](OutputReference.en.html#locate-configuration-specific-values).

### `charges.csv`

`charges.csv` contains final-state `elem_idx,charge_C`. `charge_C` is the total charge [C] on each triangle element,
not surface charge density. The official case has 288 rows, and their sum agrees with
`charge_sum=-1.192019e-10` from `beachx inspect`.

```bash
head -n 3 outputs/tutorial/charges.csv
```

### Mesh files

`mesh_triangles.csv` contains triangle vertices, element charge, and `mesh_id`. `mesh_sources.csv` maps each `mesh_id`
to its input mesh, template, and surface model. The official case has one plane mesh with 288 triangles.

```bash
head -n 3 outputs/tutorial/mesh_triangles.csv
head -n 2 outputs/tutorial/mesh_sources.csv
```

## Charge ledger

`charge_ledger.csv` accumulates injection, surface absorption, escape, and unresolved discard for each species as signed
charge and counts. The official case has one electron species, so this short check compares particle counts, charge,
and final surface charge together.

```bash
python - <<'PY'
import csv
from math import fsum
from pathlib import Path

out = Path("outputs/tutorial")
with (out / "charge_ledger.csv").open(newline="", encoding="utf-8") as f:
    row = next(csv.DictReader(f))
with (out / "charges.csv").open(newline="", encoding="utf-8") as f:
    surface = fsum(float(item["charge_C"]) for item in csv.DictReader(f))

terminal_count = sum(int(row[name]) for name in (
    "absorbed_count", "escaped_count", "discarded_unresolved_count"
))
terminal_charge = fsum(float(row[name]) for name in (
    "absorbed_on_surface_C", "escaped_to_infinity_C", "discarded_unresolved_C"
))
print(f"counts: injected={row['injected_count']} terminal={terminal_count}")
print(f"charge_C: injected={float(row['injected_from_remote_C']):.12e} "
      f"terminal={terminal_charge:.12e}")
print(f"surface_C: charges={surface:.12e} "
      f"absorbed={float(row['absorbed_on_surface_C']):.12e}")
PY

grep -E '^(charge_ledger_surface_charge_after_C|charge_ledger_residual_C)=' \
  outputs/tutorial/summary.txt
```

The current reference result satisfies these relationships:

```text
counts: injected=4000 terminal=4000
charge_C: injected=-1.281741307200e-10 terminal=-1.281741307200e-10
surface_C: charges=-1.192019415696e-10 absorbed=-1.192019415696e-10
charge_ledger_surface_charge_after_C= -1.1920194156960000E-10
charge_ledger_residual_C= 5.5497729723797089E-25
```

In this case, every injected particle ends as absorbed, escaped, or unresolved discard. There is no surface emission or
external correction, so absorbed charge equals the sum of `charges.csv`. Confirm that `charge_ledger_residual_C` is
near zero at roundoff scale, but do not treat a small conservation residual as evidence of statistical convergence or
physical validity. See the [charge-ledger reference](OutputReference.en.html#charge-ledger) for the general balance with
corrections and every column.

## History

The official case saves charge and potential at every batch.

| File | What it contains |
| --- | --- |
| `charge_history.csv` | Per-element charge and `rel_change` by batch |
| `potential_history.csv` | Per-element potential by batch |
| `top_reference_history.csv` | Potential statistics on the box z-high face |

```bash
head -n 3 outputs/tutorial/charge_history.csv
head -n 3 outputs/tutorial/potential_history.csv

beachx animate outputs/tutorial \
  --quantity charge \
  --save-gif outputs/tutorial/charge_history.gif \
  --total-frames 20
```

The official case has 20 snapshots from batch 1 through 20, showing negative charge accumulating over the run.
Potential-history reconstruction requires the native field kernel; use the optional steps in the
[Post-processing Tutorial](PostprocessTutorial.en.html) when needed.
See the [history reference](OutputReference.en.html#history) for generation conditions and the complete
`matching_plane_history.csv` contract.

## Mesh potential

`mesh_potential.csv` contains the final potential [V] at each element centroid. The official case has 288 rows; its
one-thread reference range is `potential_min=-4.671330e+00` to `potential_max=-2.579807e+00`.
`potential_history.csv` is the per-batch history, whereas `mesh_potential.csv` is the final state.

See the [mesh-potential reference](OutputReference.en.html#mesh-potential) for potential reference, periodic2, and field
reconstruction conditions.

## Files used for resume

Do not manually combine checkpoint files. Set the read source with `output.restart_from`. To advance the official
20-batch result to batch 21, follow
[Resume once from a checkpoint](Execution.en.html#resume-once-from-a-checkpoint).

Consult the [checkpoint output reference](OutputReference.en.html#files-used-for-resume) only when you need required
files, periodic-slot selection, or schema compatibility.

## Locate configuration-specific values

The ordinary official case does not require the following details. Select only the row for the model or diagnostic used.

| Need | Reference |
| --- | --- |
| File conditions | [Generation](OutputReference.en.html#file-generation-conditions) |
| Field, periodic2, and boundary | [Receipts](OutputReference.en.html#locate-configuration-specific-values) |
| Zhao currents and corrections | [`zhao_stationary`](OutputReference.en.html#zhao_stationary) |
| Accepted matching-plane state | [`matching_plane_quasistatic`](OutputReference.en.html#matching_plane_quasistatic) |
| Every `charge_ledger.csv` column | [Charge ledger](OutputReference.en.html#charge-ledger) |
| Adaptive periodic2 trials | [Adaptive-batch diagnostics](OutputReference.en.html#adaptive-batch-diagnostics) |
| Checkpoint schemas and required files | [Files used for resume](OutputReference.en.html#files-used-for-resume) |
| Python-reader attributes | [Reading from Python](OutputReference.en.html#reading-from-python) |

### Adaptive-batch diagnostics

For `periodic2.max_nonzero_mode_potential_step > 0`, inspect its receipt and trials excluded from histories in the
[adaptive-batch reference](OutputReference.en.html#adaptive-batch-diagnostics).

## Continue

| Goal | Next page |
| --- | --- |
| Visualize charge and potential | [Post-processing Tutorial](PostprocessTutorial.en.html) |
| Resume from a checkpoint | [Execution and Resume](Execution.en.html#resume-once-from-a-checkpoint) |
| Decide whether a research result is acceptable | [Validating Simulation Results](ValidationGuide.en.html) |
| Diagnose missing or inconsistent values | [Troubleshooting](Troubleshooting.en.html) |
| Search the output contract | [Output Format Reference](OutputReference.en.html) |

### Reading from Python

For a first Python load and visualization, follow the [Post-processing Tutorial](PostprocessTutorial.en.html).
For the complete class and function contract, see the [Python Post-processing API](PythonPostprocessAPI.en.html).
