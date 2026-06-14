# BEACH Usage Workflows

## New Case

```bash
mkdir run_periodic2
cd run_periodic2
beachx config init
beachx config validate
beachx config render
beach beach.toml
```

`beach.toml` is both the everyday editing entry point and the file read by the Fortran runtime. When high-level notation is used, `beachx config render` expands it to final runtime keys.

## Direct Rendered-Config Run

```bash
beach examples/beach.toml
```

With no argument, `beach` reads `beach.toml` from the current directory. In a development checkout, after `make`, `fpm run --profile release --flag "-fopenmp" -- examples/beach.toml` is also available.

## Output Inspection

```bash
beachx inspect outputs/latest
beachx animate outputs/latest --quantity charge --save-gif outputs/latest/charge_history.gif
beachx slices outputs/latest --grid-n 200 --save outputs/latest/potential_slices.png
beachx coulomb outputs/latest --component z --save outputs/latest/coulomb_force_z.png
beachx mobility outputs/latest --density-kg-m3 2500 --mu-static 0.4 --save-csv outputs/latest/mobility_summary.csv
```

## Python Analysis

```python
from beach import Beach

run = Beach("outputs/latest")
print(run.result.absorbed, run.result.last_rel_change)
run.plot_mesh()
run.plot_potential()
```

## Workload Estimate

```bash
beachx estimate-workload beach.toml --threads 8
```

For `reservoir_face` and `photo_raycast`, `batch_duration` / `batch_duration_step`, macro-particle count, mesh element count, and history output strongly affect compute cost and output size.
