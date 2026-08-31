title: Case Design Workflow

Lang: [English](ConfigurationRecipes.en.md) | [日本語](ConfigurationRecipes.md)

# Case Design Workflow

This page gives the decision order for turning the tested official tutorial into a research case. It does not reproduce
the complete configuration. Each step shows only the minimum change, while the
[input parameter reference](Parameters.en.html) remains canonical for all keys and combination constraints.

**Starting point:** Complete the [10-minute tutorial](Tutorial.en.html). In the working directory that contains
`beach.toml` and `outputs/tutorial`, copy the baseline configuration.

```bash
cp beach.toml case.toml
beachx lint case.toml
```

Change one decision at a time and rerun `beachx lint case.toml` after every change. Keep the tutorial output as the
baseline; step 7 assigns a different output directory to the new case.

## 1. Define the purpose and acceptance criteria

Before editing the configuration, write down one sentence for each of these questions:

- Which particles arrive from which environment, and which surface do they reach?
- Is the quantity of interest charge, potential, absorption, escape, or runtime?
- Which conservation check, reference solution, or convergence tolerance must the result pass?

The official tutorial is a teaching case for inter-batch charging feedback. Its particle weight and
`batch_duration=0` do not describe the time evolution of a particular physical environment. Defining acceptance
criteria first limits mesh resolution, particle count, and solver accuracy to what the study actually needs.

## 2. Choose the geometry

Use built-in shapes with `mode="template"` for simple geometry. Use a triangular OBJ for measured or CAD-derived
geometry. When separate objects, especially separate conductors, need distinct identities, use separate templates so
that they receive separate `mesh_id` values.

| `kind` | Shape | Main dimensions | Main resolution keys |
| --- | --- | --- | --- |
| `plane` | XY rectangle | `size_x`, `size_y` | `nx`, `ny` |
| `plate_hole`, `plane_hole` | Rectangle with a circular hole; the latter is an alias | `size_x`, `size_y`, `radius` | `n_theta`, `n_r` |
| `disk` | Disk | `radius` | `n_theta`, `n_r` |
| `annulus` | Annulus | `radius`, `inner_radius` | `n_theta`, `n_r` |
| `box` | Closed rectangular surface | `size` | `nx`, `ny`, `nz` |
| `cylinder` | z-axis cylinder | `radius`, `height`, `cap` | `n_theta`, `n_z` |
| `sphere` | Sphere | `radius` | `n_lon`, `n_lat` |

For example, to replace the tutorial plane with a sphere, replace only the existing `[[mesh.templates]]` block:

```toml
[[mesh.templates]]
kind = "sphere"
enabled = true
surface_side = "outward_closed"
center = [0.5, 0.5, 0.5]
radius = 0.2
n_lon = 24
n_lat = 12
```

For OBJ input, remove the existing `[[mesh.templates]]` block and replace `[mesh]` with:

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj" # replace with the actual OBJ path
surface_side = "outward_closed"
```

Use `outward_closed` only for a consistently oriented, closed two-manifold. For an open surface, select `normal_plus`
or `normal_minus` from the OBJ face normals and the vacuum side. See the
[built-in shape and OBJ reference](Parameters.en.html#mesh-mesh-input) for placement transforms, element counts, and
shape constraints. [`examples/beach.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/beach.toml) shows
multiple built-in shapes.

## 3. Choose the surface model

| Goal | `surface_model` | Scope of the current model |
| --- | --- | --- |
| Retain charge at the hit location | `insulator` | Does not solve surface conduction, bulk leakage, or dielectric polarization |
| Conserve object charge and equalize its potential | `conductor` | Floating conductor in free space; requires `field_boundary.mode="free"` |

The tutorial uses `insulator`. To change a target to a conductor, change this one line in that template or in `[mesh]`
for OBJ input:

```toml
surface_model = "conductor"
```

`dielectric`, `epsilon_r`, and resistive surfaces are not available in the current input. A conductor case needs its own
validation of object potential and total charge; the tutorial result is not its baseline. See
[How surfaces charge](SurfaceModels.en.html) for model meanings and constraints.

## 4. Choose the particle source

Choose a path from where particles are supplied. When changing `source_mode`, replace that species entry with the
minimum configuration from the dedicated guide so that keys specific to the old mode do not remain.

| `source_mode` | Use in a new case | Physical supply | Macro-sampling control |
| --- | --- | --- | --- |
| `volume_seed` | Place a specified count inside the box for initial populations or orbit tests | Does not represent a surface flux | `npcls_per_step` |
| `plane_source` | Supply particles one way from an internal rectangle | Flux × area × `batch_duration` | `w_particle` or `target_macro_particles_per_batch` |
| `photo_raycast` | Emit from the first surface hit by illumination rays | Current density × projected area × `batch_duration` | `rays_per_batch` and ray hit rate |
| `reservoir_face` | Compatibility input for deprecated existing configurations only | Do not use for a new case | Do not use for a new case |

To keep the tutorial `volume_seed` path but change its supply and region, adjust only these values in the existing
species entry:

```toml
npcls_per_step = 500
pos_low = [0.2, 0.2, 0.8]
pos_high = [0.8, 0.8, 0.8]
drift_velocity = [0.0, 0.0, -1.0e6]
```

Inflow from external plasma through a complete nonperiodic box face uses `[particles.species.boundary_inflow]`, not
another `source_mode`. The current schema also requires `source_mode="volume_seed"` and `npcls_per_step=0` for a species
that uses only boundary inflow. See [Choose where particles enter](ParticleSourcesBoundaries.en.html) for selection,
[Inflow through a simulation boundary](ReservoirInjection.en.html) for boundary inflow, and
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for photoelectrons.

## 5. Choose the box, field boundary, and particle boundaries

First make the `[domain]` box contain the mesh and every reachable particle trajectory. Then select the field closure
with `field_boundary.mode` and the outward particle action on nonperiodic faces with `[particle_boundary]`.

| Decision | Selection |
| --- | --- |
| Free-space field of finite objects | `periodic_axes=[]`, `field_boundary.mode="free"` |
| Field repeated infinitely in x/y | `periodic_axes=["x", "y"]`, `field_boundary.mode="periodic2"` |
| Particle reaching a nonperiodic face | `open`, `reflect`, or `redistributed_reflect` |

`domain.periodic_axes` is the only public setting for periodic topology. Neither `[particle_boundary]` nor a species
override can add, remove, or replace a periodic axis. Because `periodic2` currently requires periodic x/y and
nonperiodic z, start a new periodic case from the complete
[`examples/periodic2_basic/beach.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/periodic2_basic/beach.toml)
and [periodic2 electrostatics](PeriodicElectrostatics.en.html). See
[Particle escape and return](ParticleEscapeReturn.en.html) for particle-boundary behavior.

## 6. Choose the solver

| `sim.field_solver` | When to select it | How to check it |
| --- | --- | --- |
| `direct` | Small cases and reference solutions | Use as the baseline for approximate solvers |
| `treecode` | Medium free-space cases | Vary the opening condition and compare observables |
| `fmm` | Large cases, many targets, and ordinary `periodic2` production | Vary `tree_theta` and `tree_leaf_max`, then compare with Direct |
| `auto` | Select by element count for a free-space case | Read the resolved solver from output |

Keep the tutorial reference run at `field_solver="direct"`. Before adopting a faster solver, measure its difference from
Direct in a reduced case with the same mesh and particle conditions. See [Field evaluation](FieldSolvers.en.html) for
compatibility and selection, and [Use FMM](FMM.en.html) for FMM configuration and accuracy tuning.

## 7. Separate the output and run the case

Finally, change only `dir` in the existing `[output]` table so that the baseline result is not overwritten.

```toml
dir = "outputs/case"
```

On HPC systems, do not run the simulation directly on a login node; run `beach case.toml` inside a compute-node
allocation.

```bash
beachx lint case.toml
beach case.toml
beachx inspect outputs/case
```

See [Run and resume a simulation](Execution.en.html) for local execution, MPI, checkpoints, and resume.

The minimum pass condition is `status=ok` from `lint`, process exit code 0, `batches` from `inspect` equal to
`sim.batch_count`, and the presence of `summary.txt` and `charges.csv`. This proves completion, not physical correctness.
See [Inspect output files](OutputGuide.en.html) for interpreting each output.

## 8. Validate the result

For the observable defined in step 1, check at least the following:

- `processed_particles = absorbed + escaped`; its `escaped_boundary`, `survived_max_step`, and
  `multiple_box_events_soft_discarded` components are explained; and the surface/emission charge ledger closes
- Refining the mesh does not change the conclusion
- Varying `dt`, `max_step`, `batch_duration`, batch count, and particle count keeps the conclusion within tolerance
- A small Treecode or FMM result agrees with the Direct reference
- Statistical variation across random seeds or repeated runs remains within tolerance
- The conclusion stays within the scope of the selected surface model, boundaries, and particle source

Compare the same physical quantity at the same physical time before and after a change, and do not vary several
numerical controls at once. Follow [Validating simulation results](ValidationGuide.en.html) for concrete commands and
acceptance criteria.
