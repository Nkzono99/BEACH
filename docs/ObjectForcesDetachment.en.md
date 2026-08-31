title: Analyze Object Forces and Detachment

Lang: [English](ObjectForcesDetachment.en.md) | [日本語](ObjectForcesDetachment.md)

# Analyze Object Forces and Detachment

This task guide checks object-to-object forces and mobility indicators at saved positions, then evaluates a vertical
detachment path against frozen saved charges. By the end, you can save review artifacts for force, path work, and a
release decision that includes gravity and adhesion.

**This analysis cannot use the official tutorial's `outputs/tutorial`.** That case contains only one fixed plane, so it
does not define object-to-object force, mobility, or detachment. You need a separate run with multiple `mesh_id` values
and the complete `beach.toml` that was actually used for that run.

**The native field kernel is also required.** A standard Python package installation does not contain
`libbeach_field_kernel.so`, so the object-force commands and Python analysis after `inspect` do not run as installed.
Build the library in a BEACH source checkout with a Fortran compiler and `fpm`, then set its actual path.

```bash
make -C /path/to/BEACH build-kernel
export BEACH_FIELD_KERNEL_LIB=/path/to/BEACH/build/libbeach_field_kernel.so
```

The examples below assume output in `outputs/latest`, the run configuration as `beach.toml` in the current directory,
and movable `mesh_id=6`. Replace `/path/to/BEACH` and all of these values with those from your environment and run.

## 1. Check the run and target object

Read the output first and select only an object that exists.

```bash
beachx inspect outputs/latest
```

Confirm that `mesh_ids` contains multiple IDs and that `mesh_source` or `mesh_sources.csv` explains each shape. If ID
`6` is absent, replace every later `6` with an existing ID. Do not continue with a run that does not contain multiple
objects.

Only when an x/y-periodic mesh needs to be displayed inside its primary cell, use:

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --apply-periodic2-mesh
```

## 2. Check forces and mobility at the saved positions

Run these three analyses before moving an object.

```bash
beachx coulomb outputs/latest \
  --config beach.toml \
  --component z \
  --save outputs/latest/coulomb_force_z.png

beachx kernel-forces outputs/latest \
  --config beach.toml \
  --target-mesh-ids 6 \
  --save-csv outputs/latest/object_forces_kernel.csv

beachx mobility outputs/latest \
  --config beach.toml \
  --density-kg-m3 2500 \
  --mu-static 0.4 \
  --save-csv outputs/latest/mobility_summary.csv
```

`coulomb` writes a component-wise force matrix between objects. `kernel-forces` writes per-object resultants from the
field kernel. `mobility` combines force with gravity, support, and friction assumptions to produce lift, slide, and roll
indicators. Output files establish analysis completion, not the onset of motion. Validate density, friction, support,
mesh, and field-model assumptions first.

## 3. Save frozen-charge detachment artifacts

Detachment analysis holds saved source geometry and charge at their initial positions and translates only the selected
central-cell target vertically. It is not a mechanical simulation that recomputes charging or moves surrounding
objects while the target moves.

Choose the field definition first.

| `--periodic-model` | When to use it |
| --- | --- |
| `configured` | Use the run's saved free, finite, or cached setting unchanged |
| `infinite-physical` | On a compatible x/y-periodic run, compose cached `k != 0` with the physical `k = 0` mode |

The following example selects `infinite-physical`, so it requires a compatible periodic run and cache.

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

On success, the command writes `instantaneous_wrench.csv`, `path.csv`, `summary.json`, and `report.md`.
`object-detachment` defaults to lunar gravity, `1.62 m/s^2`; this example explicitly selects Earth gravity. When using
adhesion, supply `--adhesion-force-n` and `--adhesion-range-m` together.

`configured` uses the run's finite-image or cached setting. `infinite-physical` uses cached `k != 0` and the physical
`k = 0` mode. Only the target's central-cell primary self field is excluded, so force from the target's own periodic
images remains.

## 4. Extend the detachment path in Python

Use the Python API for a custom displacement array or adhesion model. This example freezes the same saved charges and
moves only mesh 6 upward.

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

print("force [N]:", wrench.force_N)
print("path status:", path.status)
print("barrier free:", release.barrier_free_from_rest)
```

See
[`examples/analyze_periodic_object_detachment.py`](https://github.com/Nkzono99/BEACH/blob/main/examples/analyze_periodic_object_detachment.py)
for a complete example. See the
[Python Post-processing API Reference](PythonPostprocessAPI.en.html#104-objectinteractionsnapshot-and-frozen-source-paths)
for every class and result field and the detailed meanings of `configured` and `infinite_physical`.

## 5. Separate completion from validity

A zero exit status and generated artifacts establish analysis completion only. Check at least:

- `path.status` and the mismatch between integrated force work and potential-difference work
- Dependence on mesh and quadrature, and on the finite image shell or periodic cache
- Dependence on the path endpoint and selected charge snapshot
- Variation in force and release decisions across stochastic seeds
- Density, friction, and support supplied to `mobility`, and mass, gravity, and adhesion supplied for detachment

Set tolerances with [Validating Simulation Results](ValidationGuide.en.html). Do not interpret a finite-height speed in
a non-neutral periodic cell as escape speed at infinity.
