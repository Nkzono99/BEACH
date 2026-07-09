title: Post-processing Tutorial

Lang: [English](PostprocessTutorial.en.md) | [日本語](PostprocessTutorial.md)

# Post-processing Tutorial

This short tutorial shows how to inspect the first run with the CLI and Python.
For file meanings, see [Reading Output Files](OutputGuide.en.html). For the full API, see [Python Post-processing API Reference](PythonPostprocessAPI.en.html).

## Inspect a Run with the CLI

```bash
beachx inspect outputs/latest
```

Save figures:

```bash
beachx inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png
```

For `sim.field_bc_mode="periodic2"`, draw the mesh wrapped into the periodic cell:

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --apply-periodic2-mesh
```

## Make a History Animation

When `charge_history.csv` exists:

```bash
beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200
```

If `potential_history.csv` is enabled, use `--quantity potential` as well.

## Read Results from Python

```python
from beach import Beach

b = Beach("outputs/latest")
print("absorbed:", b.result.absorbed)
print("escaped:", b.result.escaped)
print("batches:", b.result.batches)

b.plot_bar()
b.plot_mesh()
b.plot_potential()
```

If the config file is not near the output directory, pass it explicitly.

```python
b = Beach("outputs/latest", config_path="beach.toml")
```

## Select a Mesh

Use `mesh_sources.csv` to find `mesh_id`, then select the target mesh.

```python
from beach import Beach

b = Beach("outputs/latest")
mesh1 = b.get_mesh(1)
charge1 = b.get_mesh_charge(1)
print(mesh1.centers.shape, charge1.shape)
```

When history exists, pass a batch index.

```python
mesh1_step10 = b.get_mesh(1, step=10)
charge1_step10 = b.get_mesh_charge(1, step=10)
```

## Slices, Forces, and Mobility

For more advanced analysis, start with these CLI commands.

```bash
beachx slices outputs/latest \
  --grid-n 200 \
  --save outputs/latest/potential_slices.png

beachx coulomb outputs/latest \
  --component z \
  --save outputs/latest/coulomb_force_z.png

beachx mobility outputs/latest \
  --density-kg-m3 2500 \
  --mu-static 0.4 \
  --save-csv outputs/latest/mobility_summary.csv
```

These commands resolve geometry and periodic2 settings from nearby `beach.toml`.
If no config is found, pass the corresponding `--config` option.
