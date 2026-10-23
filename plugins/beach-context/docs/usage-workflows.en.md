# BEACH Usage Workflows

## New Case

```bash
mkdir run_periodic2
cd run_periodic2
beachx config init
beachx lint beach.toml
beach beach.toml
```

`beach.toml` is both the everyday editing entry point and the file read by the Fortran runtime. When high-level notation is used, the Fortran parser normalizes it while loading.

## Direct beach.toml Run

```bash
beach examples/beach.toml
```

With no argument, `beach` reads `beach.toml` from the current directory. In a development checkout, after `make`, `fpm run --target beach --profile release --flag "-fopenmp" -- examples/beach.toml` is also available.

## Output Inspection

```bash
beachx inspect outputs/latest
beachx animate outputs/latest --quantity charge --save outputs/latest/charge_history.gif
beachx slices outputs/latest --grid-n 200 --save outputs/latest/potential_slices.png
beachx coulomb outputs/latest --component z --save outputs/latest/coulomb_force_z.png
beachx mobility outputs/latest --density-kg-m3 2500 --mu-static 0.4 --save outputs/latest/mobility_summary.csv
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
beachx workload beach.toml --threads 8
```

For `inflow` and `photo_raycast`, `run.batch.duration_s` / `run.batch.duration_steps`, macro-particle count, mesh element count, and history output strongly affect compute cost and output size.

## Build the periodic2 far-correction cache in a separate job

With `fields.periodic.backend="cached_kneq0"`, the first run for a new box, expansion order, or Ewald setting builds the far-correction operator during initialization.
Depending on size this takes minutes to tens of minutes and consumes the time limit and the rank count of the production job (fixed on restart). Build the cache in a separate parallel job.

```bash
export BEACH_CACHE_DIR=/path/to/shared/beach-cache   # shared FS visible from compute nodes; same for all jobs
srun -n 16 beach --prepare-periodic-cache beach.toml # cache job: more ranks than production, build only
sbatch --dependency=afterok:<cache job id> run.sh    # production job after the cache job
```

- Check `periodic2_cache_hit=T` and `periodic2_operator_build_count=0` in the production `summary.txt`.
- Runs with the same box and settings share one cache even with different object sets; only missing cells are added. For sweeps over object sets, run `--prepare-periodic-cache` for every object set sequentially in one job before production.
- Builds into the same cache are serialized by a lock. Submitting many production jobs without a cache wastes the time limit of the waiting jobs.
- Changing the box size, `fields.periodic.image_layers`, `ewald_layers`, or the expansion order selects a different cache.
