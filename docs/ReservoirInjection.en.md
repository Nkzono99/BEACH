title: Inject particles through a boundary

Lang: [English](ReservoirInjection.en.md) | [日本語](ReservoirInjection.md)

# Inject particles through a boundary

This procedure injects particles through a box face from a plasma outside the box (the external plasma, or reservoir), with a
density and temperature or a velocity distribution. Make a non-periodic face open and set that face to `"reservoir"` in the
species `[particles.species.inflow]`. After reading this page you can build a minimal configuration, choose the distribution and
the inflow mapping, and check the injected amount in the output. To emit particles from a face inside the box, use
`plane_source` in [Choose a particle source](ParticleSourcesBoundaries.en.html).

## 1. Build a minimal configuration

This is the difference that injects electrons through the top face of an existing case. The values are examples and do not
represent a particular plasma environment.

```toml
[run.batch]
duration_s = 1.0e-6

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]
periodic_axes = []

[particles.boundary]
z_high = "open"

[[particles.species]]
species_key = "electron"
charge_c = -1.602176634e-19
mass_kg = 9.1093837139e-31

[particles.species.distribution]
model = "maxwellian"
number_density_m3 = 5.0e6
temperature_ev = 10.0
drift_velocity_m_s = [0.0, 0.0, -4.0e5]

[particles.species.sampling]
weight = 1.0e5

[particles.species.inflow]
z_high = "reservoir"
```

The inward normal of the top face is $-z$, so a negative z drift points into the box. Several of the six faces can be selected at once.

- Make the batch duration (`run.batch.duration_s`) positive. The number of injected particles is proportional to it.
- The inflow face must be non-periodic and open for that species.
- Particles enter through the whole selected face. A part of a face cannot be selected.
- A species that uses only boundary inflow does not have `[particles.species.source]`.

After saving the configuration, check it and run.

```bash
beachx lint beach.toml
beach beach.toml
```

Passing `lint` does not mean that the flux, weight, and potential reference are physically appropriate.

## 2. Choose the distribution

| External plasma information you have | Setting |
|---|---|
| Density, temperature, and drift velocity | `distribution.model = "maxwellian"` |
| Velocity points and distribution values from measurements or another calculation | `distribution.model = "grid"` |

### Maxwellian

From the density (`number_density_m3` or `number_density_cm3`), temperature, and drift velocity, BEACH computes the flux
$\Gamma_\mathrm{in}$ crossing each face inward. The probability of crossing a face is proportional to the inward normal
velocity, so normal velocities are sampled from the flux-weighted distribution. The expected number of macro-particles per batch is

$$
N_\mathrm{macro}
=\frac{\Gamma_\mathrm{in} A\,\Delta t_\mathrm{batch}}{w}
$$

where $A$ is the face area, $\Delta t_\mathrm{batch}$ the batch duration, and $w$ the weight. Fractions are carried to the next
batch, so the count per batch is not constant.

To set the weight physically, give `sampling.weight`; to set it from the number of samples per batch, give
`sampling.target_macro_particles_per_batch`. The latter is a target used to resolve the weight and does not change the physical
flux. The two cannot be given together.

### Velocity grid

Instead of the Maxwellian density and temperature, give a CSV file and the physical flux.

```toml
[particles.species.distribution]
model = "grid"
grid_path = "inflow_vdf.csv"
grid_pdf_kind = "phase_space"
particle_flux_m2_s = 1.0e12

[particles.species.sampling]
velocity_grid_sampling = "auto"
```

The CSV columns are `vx_m_s,vy_m_s,vz_m_s,f`. `f` is non-negative, and `grid_pdf_kind` states what it represents.

| `grid_pdf_kind` | `f` in the CSV | Weight used by BEACH |
|---|---|---|
| `phase_space` | Phase-space distribution | $\max(v_n,0)f$, multiplied by the inward normal velocity |
| `flux_weighted` | Distribution already weighted for particles crossing the face | $f$ |

Give the inflow amount with exactly one of `particle_flux_m2_s` and `current_density_a_m2`. A current density is converted to
a flux by $|J/q|$, so its sign does not set the direction; the direction comes from the CSV velocities and the inward normal of
the face. A relative CSV path is resolved against the working directory at run time.

## 3. Choose the inflow mapping

Choose `particles.reservoir.inflow_model` from where the external distribution is defined.

| Where the distribution is defined | `inflow_model` | Behavior |
|---|---|---|
| On the inflow face | `"source_vdf"` (default) | Use the configured distribution as the distribution on the inflow face |
| At infinity (upstream) | `"infinity_barrier"` | Select the particles that reach the face and change their normal velocity from the difference between the upstream potential and the mean potential of the inflow face |

```toml
[particles.reservoir]
inflow_model = "infinity_barrier"
phi_infty_v = 0.0
face_potential_grid_n = 5
```

With upstream potential $\phi_\infty$ and mean inflow-face potential $\phi_f$ at the start of the batch, the normal velocity on
the inflow face is

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}
$$

$q$ is the signed charge, so electrons and positive ions have opposite signs of $B$ for the same potential difference.

| $B$ | Particles that reach the face and the change of velocity |
|---:|---|
| $B>0$ | Only particles with $v_{n,\infty}\ge\sqrt B$ reach the face, and they decelerate on the way |
| $B=0$ | The normal velocity is unchanged |
| $B<0$ | Every particle reaches the face, and it accelerates on the way |

The tangential velocity is unchanged. $\phi_f$ is the mean over `face_potential_grid_n` × `face_potential_grid_n` points on the
inflow face, not a per-particle local potential. To test outgoing particles against the same upstream potential, set the
treatment of open faces to the potential barrier ([Particles at box boundaries](ParticleEscapeReturn.en.html#potential-barrier-potential_barrier)).
With the outer-sheath connection, the outer sheath sets the inflow mapping and `infinity_barrier` cannot be used.

## 4. Check the injected amount in the output

```bash
beachx inspect outputs/latest
grep -E '^(reservoir_inflow_map|particle_ordinary_open_model|charge_ledger_residual_C)=' \
  outputs/latest/summary.txt
head -n 2 outputs/latest/charge_ledger.csv
```

`reservoir_inflow_map` in `summary.txt` is the selected mapping (`source_vdf` or `infinity_barrier`), and
`particle_ordinary_open_model` is the treatment of open faces. In `charge_ledger.csv`, check the following for each species.

- `injected_count`: macro-particles that entered from outside the box
- `injected_from_remote_C`: their charge, with the sign of the particle charge
- `absorbed_count`, `escaped_count`, `discarded_unresolved_count`: where the injected particles went

If the expected count is much smaller than 1, `injected_count` can be 0 in the first batches until fractions accumulate. When
samples are too few, revise the weight or the target sample count before changing the physical flux. Changing the batch duration
also changes the interval between field updates, so compare according to [Choose the batch duration](BatchDurationStability.en.html).

## 5. Scope

- Boundary inflow is a local model that replaces conditions outside the box with particles on the boundary. Orbits outside the box,
  the field on the way, turning positions, flight times, space charge, and the outer sheath are not solved.
- A uniform external field has no potential at infinity. When combining it with `infinity_barrier`, define `phi_infty_v` as the
  reference of the external plasma and check its meaning separately.

All keys and constraints are in [Input parameters](Parameters.en.html#particles).
