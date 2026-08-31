title: Inject particles through a boundary

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# Inject particles through a boundary

> **Question:** How do I bring particles corresponding to a density, temperature, or velocity distribution from plasma
> outside the simulation box across a box boundary?
>
> **One-sentence answer:** Make a nonperiodic box face `open`, then set that face to `"reservoir"` in the species'
> `[particles.species.boundary_inflow]` table.

After reading this page, you should be able to create the minimum configuration, choose Maxwell or velocity-grid input,
choose a boundary distribution or a potential correction from infinity, and confirm the actual injected amount in the
outputs. If particles should originate on a surface inside the box, first compare this path with `plane_source` in
[Choose where particles enter](ParticleSourcesBoundaries.en.html).

## 1. Create the minimum configuration

The following is the minimum change that adds electron inflow through z-high to an otherwise valid case. The values
demonstrate the configuration; they are not reference values for a calibrated plasma environment.

```toml
[sim]
batch_duration = 1.0e-6

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]
periodic_axes = []

[particle_boundary]
z_high = "open"
ordinary_open_model = "escape"

[[particles.species]]
species_key = "electron"
source_mode = "volume_seed"
npcls_per_step = 0
velocity_distribution = "maxwellian"
number_density_m3 = 5.0e6
temperature_ev = 10.0
drift_velocity = [0.0, 0.0, -4.0e5]
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
w_particle = 1.0e5

[particles.species.boundary_inflow]
z_high = "reservoir"
```

The inward normal of z-high points in the negative z direction, so the negative z drift in this example points into the
box. You can enable more than one of `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, and `z_high`.

The configuration must satisfy these conditions:

- The resolved `sim.batch_duration` is positive.
- Each inflow face is absent from `domain.periodic_axes`, and its effective outward action, including any species
  override, is `open`.
- `boundary_inflow` uses a complete box face. `pos_low`, `pos_high`, and `inject_face` do not select an aperture.
- The current schema requires `source_mode="volume_seed"`. For boundary-only inflow, set `npcls_per_step=0`; this does
  not create a volume population.

After saving the configuration, check its structure and combinations before the normal run.

```bash
beachx lint beach.toml
beach beach.toml
```

A successful lint does not establish that the flux, macro-particle weight, or potential reference is physically
suitable.

## 2. Choose Maxwell or velocity-grid input

| External-plasma information available | Setting |
| --- | --- |
| Density, temperature, and drift velocity | `velocity_distribution="maxwellian"` |
| Measured or separately calculated velocity points and distribution values | `velocity_distribution="grid"` |

### Maxwell distribution

For a Maxwell distribution, BEACH obtains the one-way flux $\Gamma_\mathrm{in}$ on each face from
`number_density_m3` or `number_density_cm3`, temperature, and `drift_velocity`. Because the probability of crossing a
surface is proportional to inward normal speed, normal velocities are sampled from a flux-weighted distribution rather
than by placing the volume Maxwell distribution directly on the boundary.

The expected macro-particle count in one batch follows the relationship needed for case design:

$$
N_\mathrm{macro,expected}
=\frac{\Gamma_\mathrm{in} A\,\Delta t_\mathrm{batch}}{w},
$$

where $A$ is the selected-face area, $\Delta t_\mathrm{batch}$ is `batch_duration`, and $w$ is `w_particle`.
Fractional counts carry into the next batch, so individual batch counts need not be equal. See
[Input Parameters](Parameters.en.html#particlesspeciesboundary_inflow) and [`SPEC.md`](../SPEC.md) for the complete
flux expression and constraints.

Set `w_particle` when the physical macro-particle weight is known. To choose the sample count per batch from
computational cost instead, set `target_macro_particles_per_batch`; BEACH then resolves the weight without changing the
physical flux. Do not specify both.

### Velocity grid

For a velocity grid, remove Maxwell density and temperature keys and supply a CSV and physical flux.

```toml
velocity_distribution = "grid"
velocity_grid_path = "inflow_vdf.csv"
velocity_grid_pdf_kind = "phase_space"
velocity_grid_sampling = "auto"
particle_flux_m2_s = 1.0e12
```

The CSV requires `vx_m_s,vy_m_s,vz_m_s,f` columns with nonnegative `f`. State which distribution `f` represents.

| `velocity_grid_pdf_kind` | Meaning of CSV `f` | BEACH weight |
| --- | --- | --- |
| `phase_space` | Phase-space distribution | $\max(v_n,0)f$, including inward normal-speed weighting |
| `flux_weighted` | Distribution already weighted for particles crossing the boundary | $f$ |

Specify exactly one of `particle_flux_m2_s` and `current_density_a_m2`. Current density is converted to particle flux as
$|J/q|$, so its sign does not select the velocity direction. Direction comes from the CSV velocities and the selected
face's inward normal. A relative CSV path is resolved from the process working directory.

## 3. Choose `source_vdf` or `infinity_barrier`

| Where the external VDF is defined | `reservoir.inflow_model` | Behavior |
| --- | --- | --- |
| At the simulation boundary | `"source_vdf"` | Use the configured VDF directly at the boundary; this is the default |
| At infinity or upstream | `"infinity_barrier"` | Correct flux and normal speed for the potential difference |

When using `infinity_barrier`, define the external reference potential explicitly.

```toml
[reservoir]
inflow_model = "infinity_barrier"
phi_infty = 0.0
face_potential_grid_n = 5
```

Let $\phi_\infty$ be the upstream potential and $\phi_f$ the inflow-face mean potential at batch start. The normal speed
at the boundary is

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}.
$$

Here $q$ is the signed particle charge. The same potential difference therefore gives opposite signs of $B$ for an
electron and a positive ion.

| $B$ | Accessibility and change at the boundary |
| ---: | --- |
| $B>0$ | Only particles with $v_{n,\infty}\ge\sqrt B$ arrive, and they decelerate toward the face |
| $B=0$ | Normal speed is unchanged |
| $B<0$ | Every inward particle is accessible and accelerates toward the face |

Tangential velocity is unchanged. $\phi_f$ is the mean of face samples set by `face_potential_grid_n`, not a local
potential evaluated separately for each incoming particle.

Inflow and outflow are separate decisions. To apply a potential barrier to outward particles as well, select
`ordinary_open_model="potential_barrier"` in `[particle_boundary]`. That model uses the local potential where each
particle actually crosses the boundary to decide return or escape.

## 4. Confirm injection in the outputs

For `output.dir="outputs/latest"`, inspect the following after the run.

```bash
beachx inspect outputs/latest
grep -E '^(reservoir_inflow_map|particle_ordinary_open_model|charge_ledger_residual_C)=' \
  outputs/latest/summary.txt
head -n 2 outputs/latest/charge_ledger.csv
```

In `summary.txt`, `reservoir_inflow_map` is `source_vdf` or `infinity_barrier` according to the selected model.
Independently, `particle_ordinary_open_model` is `escape` or `potential_barrier`.

Check at least these species-resolved columns in `charge_ledger.csv`:

- `injected_count`: for a boundary-only species such as this example, macro-particles created from outside the box
- `injected_from_remote_C`: their signed injected charge, whose sign matches `q_particle`
- `absorbed_count`, `escaped_count`, and `discarded_unresolved_count`: their main tracked outcomes

When the expected count is much smaller than one, early batches can have `injected_count=0` until a fractional remainder
accumulates. If the run produces too few samples, review `w_particle` or `target_macro_particles_per_batch` before
changing the physical flux. Changing `batch_duration` changes both particle count and the field-update interval, so
compare cases as described in [Batch Duration and Stability](BatchDurationStability.en.html).

Generated outputs establish run completion only. Evaluate flux, charge conservation, statistical error, and numerical
and physical validity separately with [Inspect Output Files](OutputGuide.en.html) and
[Validating Simulation Results](ValidationGuide.en.html).

## 5. Check the model scope

- `boundary_inflow` is a local map from conditions outside the box to particles at its boundary. It does not solve
  outside-box trajectories, intermediate $E(z)$, turning positions, flight times, space charge, or a self-consistent
  external sheath.
- A uniform external electric field has no finite potential at infinity. When using it with `infinity_barrier`, define
  `phi_infty` as an effective reservoir reference and validate that interpretation separately.
- Inflow positions span the complete selected box face. Use `source_mode="plane_source"` for a rectangle placed inside
  the box.
- `boundary_inflow` does not override outward boundary actions and does not enable `[outer_plasma]` or `[coupling]`.
  These removed tables are rejected as input.

The canonical sources for all keys, defaults, exclusions, remainder state, and checkpoint format are
[Input Parameters](Parameters.en.html) and [`SPEC.md`](../SPEC.md). MPI distribution, RNG state, remainder ownership,
and implementation responsibilities are separated into the [Developer Architecture](Architecture.en.html).

## 6. Migrate from legacy `reservoir_face`

`source_mode="reservoir_face"` is a deprecated compatibility input for existing cases. Move a new external-plasma
condition to `[particles.species.boundary_inflow]`, and move a rectangular emission surface inside the box to
`source_mode="plane_source"`.

The legacy form can define an aperture with `inject_face` and `pos_low` / `pos_high`, whereas `boundary_inflow` uses the
complete box face. If migration changes the area, compare both `injected_count` and `injected_from_remote_C`, rather
than only the configured density or flux, to confirm the intended total inflow. BEACH does not silently convert the
legacy setting to either new form. See [Input Parameters](Parameters.en.html#source_mode--reservoir_face-deprecated) and
[`SPEC.md`](../SPEC.md) for compatibility behavior and the complete migration constraints.
