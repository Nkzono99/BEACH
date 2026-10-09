title: Set up a periodic surface in plasma

Lang: [English](PeriodicPlasmaSurface.en.md) | [日本語](PeriodicPlasmaSurface.md)

# Set up a periodic surface in plasma

This page builds a case in which one x/y-periodic cell represents part of a surface exposed to the solar wind and
ultraviolet light, such as lunar regolith. The cell is much shorter than the Debye length, so you must decide how the
plasma and sheath beyond the top face (z-high) are treated. The page chooses one of four ways to close the top face,
shows the configuration difference for each, and lists the outputs to check after the run.

To obtain charging consistent with the outer sheath, use the outer-sheath connection (Zhao stationary sheath).
To compare only the redistribution of charge within the surface, closed photoelectrons are the simplest choice.

## Common setup

Every approach starts from the following setup. The complete configuration is
[`examples/grouped/closed_photoelectron.toml`](../examples/grouped/closed_photoelectron.toml).

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0e-4, 1.0e-4, 1.0e-3]
periodic_axes = ["x", "y"]

[fields]
boundary = "periodic2"

[fields.solver]
method = "fmm"

[particles.boundary]
z_low = "open"
z_high = "open"

[[particles.species]]          # solar-wind electrons (ions have the same form)
species_key = "solar_wind_electron"

[particles.species.inflow]
z_high = "reservoir"

[[particles.species]]          # photoelectrons
species_key = "photoelectron"

[particles.species.source]
mode = "photo_raycast"
inject_face = "z_high"
deposit_opposite_charge_on_emit = true
```

- **Solar wind:** electrons and ions enter through the top face by boundary inflow ([Inject through a boundary](ReservoirInjection.en.html)).
- **Photoelectrons:** illumination enters through the top face, and photoelectrons are emitted from the surfaces it hits ([Photoelectron emission and charge closures](PhotoelectronEmission.en.html)).
- **Field:** the x/y-periodic field is computed with FMM ([periodic2 electrostatics](PeriodicElectrostatics.en.html)).

## Choose how to close the top face

| Approach | What happens at the top face | Total surface charge | Potential reference | When to use |
|---|---|---|---|---|
| 1. Open | Every particle that reaches the top escapes | As tracked | None (read relative to the top-face mean) | Baseline for comparison |
| 2. Closed photoelectrons | Only photoelectrons are reflected at the top, and the net photoelectron current is matched to zero | Changes only by the solar wind | None (same as above) | Study redistribution within the surface |
| 3. Scalar barrier | Incoming and outgoing particles are filtered by the potential difference between the upstream plasma and the top face | As tracked | Upstream potential `phi_infty_v` | Comparison with a single barrier |
| 4. Outer sheath | The zero-current root of the outer sheath sets the currents, barriers, and potential reference | Constrained to zero by the floating condition | Upstream plasma at 0 V | Charging consistent with the outer sheath |

Choose one approach. Do not combine closed photoelectrons with the scalar barrier in one run, and
`particles.reservoir.inflow_model="infinity_barrier"` cannot be used with the outer sheath.

### 1. Open

Keep the common setup. Particles that reach the open top face escape (`particles.boundary.open_model="escape"`, the default).
Because the cell is shorter than the Debye length, many photoelectrons that reach the top would in reality be pushed back
to the surface by the outer sheath. Use the open top only as a baseline for the other approaches.

### 2. Closed photoelectrons

In the photoelectron species, make the top face reflecting and set the charge closure to `neutral_return`.

```toml
[[particles.species]]
species_key = "photoelectron"

[particles.species.charging]
closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"          # use "redistributed_reflect" to return at uniform in-plane positions
```

The charge of returned photoelectrons is matched to the emitted charge of the same batch so that the net photoelectron
current is zero. The change of total charge by the solar wind is not constrained. The mechanism and the stop conditions are in
[Photoelectron emission and charge closures](PhotoelectronEmission.en.html#closed-photoelectrons-neutral_return).

### 3. Scalar barrier

Set the upstream potential of the external plasma, and test both inflow and outflow against it.

```toml
[particles.reservoir]
inflow_model = "infinity_barrier"
phi_infty_v = 0.0
face_potential_grid_n = 5

[particles.boundary]
open_model = "potential_barrier"
```

Inflow is tested with the mean potential of the top face, and outflow with the potential at the crossing point
([inflow mapping](ReservoirInjection.en.html#3-choose-the-inflow-mapping), [potential barrier](ParticleEscapeReturn.en.html#potential-barrier-potential_barrier)).
The potential of the top face is set only by the charge inside BEACH and does not include the potential drop of the outer sheath.

### 4. Outer sheath (Zhao stationary sheath)

Solve the zero-current root of the outer sheath and fix the solar-wind and photoelectron currents to the root.
Set `charging.closure="fixed_current"` for the three role species, and use `cached_kneq0`, which treats the zero mode
of the field separately.

```toml
[sheath]
closure = "zero_current"

[sheath.species]
electron = "solar_wind_electron"
ion = "solar_wind_ion"
photoelectron = "photoelectron"

[sheath.photoelectrons]
solar_elevation_deg = 60.0
ref_density_m3 = 6.4e7

[sheath.coupling]
outflow_refresh_batches = 50   # omit to keep the initial root

[fields.periodic]
backend = "cached_kneq0"
lower_boundary_model = "e_bottom_zero"
```

Set the solar-wind electron drift to zero and the ion drift to the normal solar-wind velocity (inward, negative z).
Write the same solar-wind density for electrons and ions. The complete configuration is
[`examples/grouped/zero_current.toml`](../examples/grouped/zero_current.toml); the example that refreshes the outer root is
[`zero_current_refresh.toml`](../examples/grouped/zero_current_refresh.toml), and the example without photoelectrons (Type C) is
[`zero_current_no_photo.toml`](../examples/grouped/zero_current_no_photo.toml).
The model assumptions and known limitations are in [Connecting to the outer sheath](ZhaoStationaryClosure.en.html).

## Check the results

| Approach | Output | Decision |
|---|---|---|
| All | `charge_ledger_residual_C` in `charge_ledger.csv` and the charge of unresolved particles | The residual is at round-off level and the unresolved charge does not affect the conclusion |
| All | `potential_mean_V` and `potential_std_V` in `top_reference_history.csv` | For approaches 1–3, read each element potential relative to this mean. A large spread weakens the treatment of the top face as one plane |
| 2 | `neutral_return_weight_scale` and `neutral_return_unresolved_fraction` in `charge_ledger.csv` | The scale is near 1 and the unreturned fraction is small (the run stops above 5%) |
| 3 | `reservoir_inflow_map` and `particle_ordinary_open_model` in `summary.txt` | The selected mapping and test are in use |
| 4 | `surface_current_model_zhao_branch` and `surface_current_model_phi0_V` in `summary.txt` | The branch and wall potential are as expected |
| 4 | `target_over_tracked` in `fixed_current_history.csv` | The scale of each species stays near 1 |
| 4 (outer-root refresh) | `phi_H_V` and the photoelectron outflow in `matching_plane_history.csv` | The wall potential and the outflow have settled |

`top_reference_history.csv` is written when `output.history.potential=true` and `output.history.stride_batches>0`.
The columns are defined in [Output formats](OutputReference.en.html#history).

## Check the dependence on numerical settings

For every approach, vary at least the following and confirm that the quantity of interest does not change.

1. Periodic field: with finite images, increase `fields.periodic.image_layers` to $N, N+1, N+2$. With `cached_kneq0`, check
   [the accuracy of the far correction](PeriodicFarCorrection.en.html).
2. Cell height: move the top face up and down.
3. Batch duration: compare $T, T/2, T/4$ at the same physical time ([Choose the batch duration](BatchDurationStability.en.html)).
4. Particle tracking: halve `particles.tracking.dt_s` and double `max_steps_per_particle`.
5. Sample counts: vary the macro-particle count, the ray count (`sampling.rays_per_batch`), and the random seed.

How to set the acceptance criteria is described in [Validate results](ValidationGuide.en.html).

## Scope

- The cell represents part of a flat surface. Approaches 1–3 do not solve the sheath outside the cell; approach 4 represents it by a planar stationary solution.
- No approach solves particle orbits, flight times, or space charge outside the box.
- In approaches 1–3, potentials are meaningful only relative to the top-face mean. If you need potentials referenced to the upstream plasma, use approach 4.
