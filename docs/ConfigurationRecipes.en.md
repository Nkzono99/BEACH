title: Case Design Workflow

Lang: [English](ConfigurationRecipes.en.md) | [日本語](ConfigurationRecipes.md)

# Case Design Workflow

Configuration files accept only the seven-group layout (`run / domain / mesh / particles / fields / sheath / output`). See `examples/` for complete input files.

This page turns the tested official tutorial into a research case: choose geometry, sources, and boundaries, then set
time steps and particle sampling to obtain comparable results. Look up types, defaults, and combination constraints in
the [input parameter reference](Parameters.en.html).

**Starting point:** Complete the [10-minute tutorial](Tutorial.en.html). In the working directory that contains
`beach.toml` and `outputs/tutorial`, copy the baseline configuration.

```bash
cp beach.toml case.toml
beachx lint case.toml
```

Change one decision at a time and rerun `beachx lint case.toml` after every change. Keep the tutorial output as the
baseline; step 8 assigns a different output directory to the new case.

## 1. Define the purpose and acceptance criteria

Before editing the configuration, write down one sentence for each of these questions:

- Which particles arrive from which environment, and which surface do they reach?
- Is the quantity of interest charge, potential, absorption, escape, or runtime?
- Which conservation check, reference solution, or convergence tolerance must the result pass?

The official tutorial is a teaching case for inter-batch charging feedback. Its particle weight and
`run.batch.duration_s=0` do not describe the time evolution of a particular physical environment. Defining acceptance
criteria first limits mesh resolution, particle count, and solver accuracy to what the study actually needs.

## 2. Choose the geometry

Use built-in template shapes for simple geometry and triangular OBJ files for measured or CAD-derived geometry.
List both in `[[mesh.parts]]`. When separate objects, especially separate conductors, need distinct identities, use separate parts so
that they receive separate `mesh_id` values.

| `kind` | Shape | Main dimensions | Main resolution keys |
| --- | --- | --- | --- |
| `obj` | OBJ file | `path`, `scale`, `rotation`, `offset` | Faces in the file |
| `plane` | XY rectangle | `size_x`, `size_y` | `nx`, `ny` |
| `plate_hole`, `plane_hole` | Rectangle with a circular hole; the latter is an alias | `size_x`, `size_y`, `radius` | `n_theta`, `n_r` |
| `disk` | Disk | `radius` | `n_theta`, `n_r` |
| `annulus` | Annulus | `radius`, `inner_radius` | `n_theta`, `n_r` |
| `box` | Closed rectangular surface | `size` | `nx`, `ny`, `nz` |
| `cylinder` | z-axis cylinder | `radius`, `height`, `cap` | `n_theta`, `n_z` |
| `sphere` | Sphere | `radius` | `n_lon`, `n_lat` |

For example, to replace the tutorial plane with a sphere, replace only the existing `[[mesh.parts]]` block:

```toml
[[mesh.parts]]
kind = "sphere"
enabled = true
surface_side = "outward_closed"
center = [0.5, 0.5, 0.5]
radius = 0.2
n_lon = 24
n_lat = 12
```

Write OBJ input as another `[[mesh.parts]]` block. Append it to mix it with templates.

```toml
[[mesh.parts]]
kind = "obj"
path = "examples/tetrahedron.obj"
surface_side = "outward_closed"
surface_model = "insulator"
```

Use `outward_closed` only for a consistently oriented, closed two-manifold. For an open surface, select `normal_plus`
or `normal_minus` from the OBJ face normals and the vacuum side. See the
[built-in shape and OBJ reference](Parameters.en.html#mesh) for placement transforms, element counts, and
shape constraints. [`examples/beach.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/beach.toml) shows
multiple built-in shapes.

## 3. Choose the surface model

| Goal | `surface_model` | Scope of the current model |
| --- | --- | --- |
| Retain charge at the hit location | `insulator` | Does not solve surface conduction, bulk leakage, or dielectric polarization |
| Conserve object charge and equalize its potential | `conductor` | Floating conductor; available in free space and `periodic2` |

The tutorial uses `insulator`. To change a target to a conductor, change this one line in that template or in `[mesh]`
for OBJ input:

```toml
# Change this in the existing mesh.parts entry
surface_model = "conductor"
```

`dielectric`, `epsilon_r`, and resistive surfaces are not available in the current input. A conductor case needs its own
validation of object potential and total charge; the tutorial result is not its baseline. See
[How surfaces charge](SurfaceModels.en.html) for model meanings and constraints.

## 4. Choose the particle source

Choose a path from where particles are supplied. When changing `source.mode`, replace that species entry with the
minimum configuration from the dedicated guide so that keys specific to the old mode do not remain.

| `source.mode` | Use in a new case | Physical supply | Macro-sampling control |
| --- | --- | --- | --- |
| `volume_seed` | Place a specified count inside the box for initial populations or orbit tests | Does not represent a surface flux | `sampling.volume_macro_particles_per_batch` |
| `photo_raycast` | Emit from the first surface hit by illumination rays | Current density × projected area × `run.batch.duration_s` | `sampling.rays_per_batch` and ray hit rate |
| `secondary_emission` | Created by incident-species absorption and tracked in the same batch | [Secondary electron emission](SecondaryElectronEmission.en.html) |

To keep the tutorial `volume_seed` path but change its supply and region, adjust only these values in the existing
species entry:

```toml
[particles.species.sampling]
volume_macro_particles_per_batch = 500

[particles.species.source]
pos_low = [0.2, 0.2, 0.8]
pos_high = [0.8, 0.8, 0.8]

[particles.species.distribution]
drift_velocity_m_s = [0.0, 0.0, -1.0e6]
```

Inflow from external plasma through a complete nonperiodic box face uses `[particles.species.inflow]`, not
another `source.mode`. For boundary inflow alone, specify the physical parameters and inflow faces, and omit
`source.mode` and `sampling.volume_macro_particles_per_batch`. See [Choose where particles enter](ParticleSourcesBoundaries.en.html) for selection,
[Inflow through a simulation boundary](ReservoirInjection.en.html) for boundary inflow, and
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for photoelectrons.

## 5. Choose the box, field boundary, and particle boundaries

First make the `[domain]` box contain the mesh and every reachable particle trajectory. Then select the field closure
with `fields.boundary` and the outward particle action on nonperiodic faces with `[particles.boundary]`.

| Decision | Selection |
| --- | --- |
| Free-space field of finite objects | `periodic_axes=[]`, `fields.boundary="free"` |
| Field repeated infinitely in x/y | `periodic_axes=["x", "y"]`, `fields.boundary="periodic2"` |
| Particle reaching a nonperiodic face | `open`, `reflect`, or `redistributed_reflect` |

`domain.periodic_axes` is the only public setting for periodic topology. Neither `[particles.boundary]` nor a species
override can add, remove, or replace a periodic axis. Because `periodic2` currently requires periodic x/y and
nonperiodic z, start a new periodic case from the complete
[`examples/periodic2_basic.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/periodic2_basic.toml)
and [periodic2 electrostatics](PeriodicElectrostatics.en.html). See
[Particle escape and return](ParticleEscapeReturn.en.html) for particle-boundary behavior.

## 6. Choose the solver

| `fields.solver.method` | When to select it | How to check it |
| --- | --- | --- |
| `direct` | Small cases and reference solutions | Use as the baseline for approximate solvers |
| `treecode` | Medium free-space cases | Vary the opening condition and compare observables |
| `fmm` | Large cases, many targets, and ordinary `periodic2` production | Vary `fields.solver.tree.theta` and `fields.solver.tree.leaf_max`, then compare with Direct |
| `auto` | Select by element count for a free-space case | Read the resolved solver from output |

Keep the tutorial reference run at `fields.solver.method="direct"`. Before adopting a faster solver, measure its difference from
Direct in a reduced case with the same mesh and particle conditions. See [Field evaluation](FieldSolvers.en.html) for
compatibility and selection, and [Use FMM](FMM.en.html) for FMM configuration and accuracy tuning.

A `secondary_emission` species has no independent injection source. Pair it with the emission settings of an incident
species ([Secondary electron emission](SecondaryElectronEmission.en.html)).

## 7. Set time scales and particle sampling

First choose physical conditions such as density, temperature, and flux, then decide how many macro-particles will
represent them. The weight `sampling.weight` is the number of physical particles represented by one simulated particle.

| Decision | Setting | How to choose it |
| --- | --- | --- |
| Particle orbit step | `particles.tracking.dt_s` | Halve it and compare impact locations and absorption rates |
| Tracking limit per particle | `particles.tracking.max_steps_per_particle` | Increase it until unresolved particles, `survived_max_step`, have a sufficiently small effect |
| Surface-charge update interval | `run.batch.duration_s` | Use positive seconds for flux-driven sources and compare with half the interval |
| Simulation duration | `run.batch_count` | With a fixed interval, `run.batch.duration_s × batch_count` is the physical end time |

When halving the orbit step, double `particles.tracking.max_steps_per_particle` to preserve the available tracking time.
With `run.batch.duration_steps`, the relation is `run.batch.duration_s = particles.tracking.dt_s × run.batch.duration_steps`, so changing `particles.tracking.dt_s` alone also
changes the charge-update interval. Specifying `run.batch.duration_s` in seconds initially makes these two controls independent
and their effects easier to compare. The two ways to specify the batch interval are mutually exclusive.

How to increase sampling depends on the selected source.

| Source | Change to improve statistics | Physical conditions to hold fixed |
| --- | --- | --- |
| `volume_seed` | Increase `sampling.volume_macro_particles_per_batch` and decrease `sampling.weight` by the same factor | Physical particles per batch, `sampling.volume_macro_particles_per_batch × sampling.weight` |
| `inflow` | Increase `sampling.target_macro_particles_per_batch`, or decrease fixed `sampling.weight`; do not specify both | Distribution, density or flux, area, and `run.batch.duration_s` |
| `photo_raycast` | Increase `sampling.rays_per_batch` | Photoelectric current density, illumination direction, geometry, and `run.batch.duration_s` |

`sampling.volume_macro_particles_per_batch` is a count **per batch**. For example, changing the tutorial's
`sampling.volume_macro_particles_per_batch=200, sampling.weight=2.0e5` to `400, 1.0e5` doubles the sample count while preserving the number of physical
particles supplied. Increasing only the particle count also increases the supplied charge and changes the physical case.
The photoelectron macro-particle count also depends on the ray hit fraction.

For the tutorial, which does not assign physical seconds, keep `run.batch.duration_s=0` and compare results by batch.
See [Choosing `run.batch.duration_s`](BatchDurationStability.en.html) for interval comparisons and adaptive stepping in
flux-driven research cases.

## 8. Separate the output and run the case

Finally, change only `dir` in the existing `[output]` table so that the baseline result is not overwritten.

```toml
[output]
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
`run.batch_count`, and the presence of `summary.txt` and `charges.csv`. This proves completion, not physical correctness.
See [Inspect output files](OutputGuide.en.html) for interpreting each output.

## 9. Validate the result

For the observable defined in step 1, check at least the following:

- `processed_particles = absorbed + escaped`; its `escaped_boundary`, `survived_max_step`, and
  `multiple_box_events_soft_discarded` components are explained; and the surface/emission charge ledger closes
- Refining the mesh does not change the conclusion
- Varying `particles.tracking.dt_s`, `particles.tracking.max_steps_per_particle`, `run.batch.duration_s`, batch count, and particle count keeps the conclusion within tolerance
- A small Treecode or FMM result agrees with the Direct reference
- Statistical variation across random seeds or repeated runs remains within tolerance
- The conclusion stays within the scope of the selected surface model, boundaries, and particle source

Compare the same physical quantity at the same physical time before and after a change, and do not vary several
numerical controls at once. Follow [Validating simulation results](ValidationGuide.en.html) for concrete commands and
acceptance criteria.
