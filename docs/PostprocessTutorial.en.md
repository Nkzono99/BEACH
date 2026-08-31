title: Post-processing Tutorial

Lang: [English](PostprocessTutorial.en.md) | [日本語](PostprocessTutorial.md)

# Post-processing Tutorial

The ordinary BEACH post-processing path is to check the run with the CLI, animate charge history, and use Python only
for custom analysis. The standard BEACH Python package is sufficient to inspect an existing official tutorial, visualize
charge, display its saved final absolute potential, and select a stored batch. Reconstructing potential history or
potential slices is an additional path that requires the native field kernel.

**Prerequisite:** Complete the [10-Minute Tutorial](Tutorial.en.html). The current directory must contain the
`beach.toml` used for that run and `outputs/tutorial`. You do not need a repository clone, the repository root, or
`examples/tutorial_insulator.toml`.

`outputs/tutorial` must contain at least `summary.txt`, `charges.csv`, `mesh_triangles.csv`, `mesh_sources.csv`,
`charge_history.csv`, `mesh_potential.csv`, and `potential_history.csv`. If these files are missing, return to the
[10-Minute Tutorial](Tutorial.en.html) and rerun it. See [Inspect Output Files](OutputGuide.en.html) for file meanings and
the [Python Post-processing API Reference](PythonPostprocessAPI.en.html) for every class and function.

## 1. Check completion and distributions with `beachx inspect`

Read the run summary before making any plots.

```bash
beachx inspect outputs/tutorial
```

For the official tutorial, check at least the following:

- `directory=outputs/tutorial`
- `processed_particles=4000` and `batches=20`
- `absorbed`, `escaped`, and `charge_sum` are reported
- `potential_min` and `potential_max` are read from the saved `mesh_potential.csv`
- `charge_history_shape` and the stored batch indices are reported
- `mesh_ids` and `mesh_source` are reported

Reference values for `rng_seed=12345` and one OpenMP thread are listed in
[Inspect Output Files](OutputGuide.en.html#check-the-official-tutorial-result). Absorption and escape counts can change
with a different thread count because they depend on the random sequence.

Next, save the standard plots from the same output.

```bash
beachx inspect outputs/tutorial \
  --save-bar outputs/tutorial/charges_bar.png \
  --save-mesh outputs/tutorial/charge.png
```

A successful command writes a per-element charge bar chart and a mesh colored by surface charge density. Generated
images establish post-processing success; they do not establish numerical convergence or physical validity.

## 2. View changes between batches with `beachx animate`

The official tutorial uses `history_stride=1`, so it has 20 charge snapshots. Render the surface-charge-density GIF,
which works with the standard distribution alone.

```bash
beachx animate outputs/tutorial \
  --quantity charge \
  --save-gif outputs/tutorial/charge_history.gif \
  --total-frames 20
```

On success, the command reports `saved_gif=outputs/tutorial/charge_history.gif`, `snapshots=20`, and
`rendered_frames=20`.

`animate` uses the snapshots stored in `charge_history.csv`. Potential mode also requires the native field kernel to
reconstruct potential from those charges and the field-reconstruction information. For another run without history,
set `output.history_stride` to a positive value and rerun the simulation.

## 3. Select plots and stored batches with the Python API

After the CLI checks pass, use Python for custom plots and calculations. Set `config_path` to the `beach.toml` used for
the run in the current directory.

```python
from beach import Beach

run = Beach(
    "outputs/tutorial",
    config_path="beach.toml",
)

print("absorbed:", run.result.absorbed)
print("escaped:", run.result.escaped)
print("batches:", run.result.batches)
print("mesh ids:", run.mesh_ids)

bar_fig, _ = run.plot_bar()
bar_fig.savefig("outputs/tutorial/python_charges_bar.png", dpi=150)

charge_fig, _ = run.plot_mesh()
charge_fig.savefig("outputs/tutorial/python_charge.png", dpi=150)

potential_fig, _ = run.plot_potential(reference_point=None)
potential_fig.savefig("outputs/tutorial/python_potential.png", dpi=150)
```

`reference_point=None` explicitly selects the final absolute potential saved in `mesh_potential.csv`, so this case does
not need the native field kernel. Re-evaluating the field with a different reference requires the kernel in the next section.

### Select a mesh and batch

Find `mesh_id` with `beachx inspect` or `mesh_sources.csv` before selecting it. `step=None` uses final charges from
`charges.csv`; an integer `step` is a batch index stored in `charge_history.csv`.

```python
mesh_id = run.mesh_ids[0]

final_mesh = run.get_mesh(mesh_id, step=None)
final_charge = run.get_mesh_charge(mesh_id, step=None)
print(final_mesh.triangles.shape, final_charge.shape)

# The official case stores batch 10 because history_stride=1.
mesh_at_10 = run.get_mesh(mesh_id, step=10)
charge_at_10 = run.get_mesh_charge(mesh_id, step=10)
print(mesh_at_10.triangles.shape, charge_at_10.shape)
```

## 4. Reconstruct potential history and slices with the native field kernel

The standard Python package installation does not contain `libbeach_field_kernel.so`. Only when you need
`animate --quantity potential` or `slices`, build the library in a BEACH source checkout with a Fortran compiler and
`fpm`. Replace `/path/to/BEACH` with the actual checkout.

```bash
make -C /path/to/BEACH build-kernel
export BEACH_FIELD_KERNEL_LIB=/path/to/BEACH/build/libbeach_field_kernel.so
```

In the shell where that environment variable is set, first render potential history.

```bash
beachx animate outputs/tutorial \
  --quantity potential \
  --save-gif outputs/tutorial/potential_history.gif \
  --total-frames 20
```

On success, the command reports `saved_gif=outputs/tutorial/potential_history.gif`, `snapshots=20`, and
`rendered_frames=20`. Then reconstruct the centered XY, YZ, and XZ slices of the simulation box.

```bash
beachx slices outputs/tutorial \
  --config beach.toml \
  --grid-n 200 \
  --save outputs/tutorial/potential_slices.png
```

On success, the command reports `saved_potential_slices=outputs/tutorial/potential_slices.png`. These analyses
reconstruct a field from saved charges. They do not replace mesh- or solver-convergence checks. Follow
[Validating Simulation Results](ValidationGuide.en.html) before accepting the result for research.

## 5. Continue to multi-object analysis

The official tutorial contains one fixed plane, so it cannot be used for object-to-object forces or detachment. If you
have a separate multi-object run, continue with [Analyze Object Forces and Detachment](ObjectForcesDetachment.en.html).
