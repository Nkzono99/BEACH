title: Post-processing Tutorial

Lang: [English](PostprocessTutorial.en.md) | [日本語](PostprocessTutorial.md)

# Post-processing Tutorial

BEACH post-processing can use the `beach` Python package directly or the `beachx` commands that package common workflows.
This page introduces the Python API first and then the ready-made visualization and analysis commands. See
[Python Post-processing API Reference](PythonPostprocessAPI.en.html) for every class and function and
[Inspect Output Files](OutputGuide.en.html) for file meanings.

## Post-process with the Python API

### Read a result and make basic plots

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

### Select a mesh

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

### Evaluate detachment force while retaining periodic images

This minimal example freezes the saved charge snapshot and moves only mesh 6 upward.

```python
import numpy as np
from beach import AdhesionProfile, Beach

run = Beach("outputs/latest", config_path="beach.toml")
with run.object_interaction_snapshot(
    periodic_model="infinite_physical",
) as snapshot:
    probe = snapshot.object_probe(6)
    wrench = probe.wrench()
    path = probe.vertical_path(np.linspace(0.0, 2.0e-4, 65))

release = path.evaluate_release(
    mass_kg=2.0e-12,
    gravity_m_s2=9.80665,
    adhesion=AdhesionProfile.none(),
)
```

See
[`examples/analyze_periodic_object_detachment.py`](https://github.com/Nkzono99/BEACH/blob/main/examples/analyze_periodic_object_detachment.py)
for a complete example.

## Use ready-made visualization and analysis with `beachx`

`beachx` exposes common Python API workflows as commands. Use it to create standard plots and CSV files without writing a script.

### Inspect a run and make basic plots

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

For `field_boundary.mode="periodic2"`, draw the mesh wrapped into the cell defined by `domain.periodic_axes`:

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --apply-periodic2-mesh
```

### Make a history animation

When `charge_history.csv` exists:

```bash
beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200
```

If `potential_history.csv` is enabled, use `--quantity potential` as well.

### Slices, forces, and mobility

For more advanced analysis, start with these commands.

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

These commands read geometry, `field_boundary.mode`, `domain.periodic_axes`, and the domain box from nearby `beach.toml` to
resolve periodic2 settings.
If no config is found, pass the corresponding `--config` option.

### Evaluate detachment force while retaining periodic images

This command runs the same flow as the Python example above and writes its CSV, JSON, and report outputs.

```bash
beachx object-detachment outputs/latest \
  --config beach.toml \
  --target-mesh-id 6 \
  --periodic-model infinite-physical \
  --z-max-m 2.0e-4 \
  --z-points 65 \
  --mass-kg 2.0e-12 \
  --gravity-m-s2 9.80665 \
  --adhesion-force-n 1.0e-10 \
  --adhesion-range-m 2.0e-6 \
  --output-dir outputs/latest/object_detachment
```

The `object-detachment` CLI defaults to lunar gravity, `1.62 m/s^2`. This
example explicitly assumes Earth gravity, `9.80665 m/s^2`; change it for the
target environment.

`configured` preserves the run's finite or cached field configuration.
`infinite-physical` combines cached `k != 0` with the physical `k = 0` mode for
an x/y-periodic run. Only the target's central-cell primary self field is
excluded, so force from the target's own periodic images remains. The command
writes `instantaneous_wrench.csv`, `path.csv`, `summary.json`, and `report.md`.

A zero exit status and generated CSV/JSON files establish execution success,
not physical qualification. Check `path.status`, work versus potential
difference, mesh/quadrature, finite shells or periodic cache, path endpoint,
charge snapshot, and stochastic-seed sensitivity using
[Validating Simulation Results](ValidationGuide.en.html). Do not interpret a
finite-height speed in a non-neutral periodic cell as escape speed at infinity.
