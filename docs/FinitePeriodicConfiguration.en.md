title: Finite-image periodic2 configuration

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# Finite-image periodic2 configuration

This page is the canonical `periodic2` setup for studying local regolith charging and photoelectron-driven surface
redistribution without solving a 1-D outer sheath. The solar-wind VDF is defined at z-high, only photoelectrons are closed,
and potential is read relative to the mean on the z-high plane rather than infinity.

The purpose is not to approximate a complete outer sheath. It is to separate net charging by the solar wind from local charge
transfer by photoelectrons with a small, stable set of assumptions.

## Distinguish the three configurations first

| Configuration | Solar-wind inflow | Photoelectrons | Potential reference | Use |
| --- | --- | --- | --- | --- |
| **Local reservoir + closed PE** | uncorrected VDF at z-high | z-high reflection + `neutral_return` | z-high plane mean | baseline on this page; local redistribution and batch sensitivity |
| Scalar barrier | `infinity_barrier` | common `potential_barrier` | `phi_infty` | comparison using one scalar barrier |
| Infinite periodic + outer plasma | map through a 1-D profile | profile return/escape | $\phi_\infty=0$ | self-consistent mean sheath |

Do not stack these alternatives in one run. In particular, do not add `infinity_barrier`, `potential_barrier`, or
`kinetic_1d` particle transfer to the closed-PE baseline.

## Recommended integrated configuration

The complete runnable example is
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml).
Its central settings are:

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
use_box = true
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
injection_face_phi_grid_n = 5

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "source_vdf"

[external_boundary.ordinary_open]
model = "escape"
```

Use ordinary `reservoir_face` species for solar-wind electrons and ions. Their density, temperature, and drift describe a
**local boundary VDF on z-high**. They are not an infinity VDF, and surface potential does not filter their accessibility or
map their speed.

Only the photoelectron species gets the closed surface-charge closure:

```toml
[[particles.species]]
species_key = "photoelectron"
source_mode = "photo_raycast"
deposit_opposite_charge_on_emit = true
z_high_boundary = "reflect"
surface_charge_closure = "neutral_return"
```

`z_high_boundary="reflect"` reverses only normal velocity at z-high. Solar-wind species retain the default `"inherit"` and
escape through the globally open boundary.

## What one batch does

1. Build a finite-image field snapshot from batch-start surface charge.
2. Inject solar wind from the local VDF at z-high.
3. Emit photoelectrons from ray hits and record opposite reaction charge on source elements.
4. Track all particles in the frozen snapshot, reflecting only photoelectrons at z-high.
5. Use the resolved photoelectron return distribution to close unresolved surface charge statistically.
6. Commit solar-wind absorption and photoelectron emission and return together.
7. Optionally write element potential and z-high mean from the same post-commit snapshot.

The field is not refreshed inside a batch. `batch_duration` is therefore both physical duration and explicit charge-update
width. Closed PE removes divergence from net photoelectron current; it does not remove solar-wind batch-width dependence.
[<sup>1</sup>](FieldSolvers.en.html)

## Quantity closed by `neutral_return`

For a negative photoelectron species, define signed macro charge within one batch:

$$
S=\sum_{\mathrm{emitted}}qw<0,\qquad
R=\sum_{\mathrm{resolved\ absorption}}qw<0.
$$

`neutral_return` multiplies deposits at measured return destinations by

$$
s_{\mathrm{return}}=\frac{S}{R}.
$$

The source reaction charge is $-S$, hence

$$
(-S)+s_{\mathrm{return}}R=0.
$$

This species does not change total surface charge; it only builds a distribution from emission sources to resolved return
destinations. Solar-wind net charge remains unconstrained.

The statistical assumption is that unresolved particles have the same return-destination distribution as particles resolved
in the same batch. When `s_return` differs substantially from one, the closure assumption rather than resolved orbits controls
the local distribution.
This fixed contract closes only a small long-lived population. BEACH stops without correction when the unresolved fraction
exceeds 5%. The limit is not configurable; converge below it by revisiting `max_step`, `dt`, and box height.

`neutral_return` does not mathematically remove every $k_\parallel=0$ structure. It makes the surface-charge monopole
increment zero, but transfer between surfaces at different heights can retain a plane-averaged vertical dipole.

The run fails without producing an accepted batch when:

- emission is nonzero but resolved return charge is zero;
- a photoelectron actually escapes through an open face;
- a `soft_discard`, nonfinite value, or charge-sign mismatch occurs;
- outer particle transfer or `implicit_mean` is combined with the closure.

Only particles still alive at `max_step` are statistically closed. The raw unresolved charge, correction, scale, and
unresolved fraction remain separate in `charge_ledger.csv`.

## Read potential relative to the z-high plane

With `output.write_files=true`, `output.write_potential_history=true`,
`output.history_stride>0`, and `sim.use_box=true`,
`top_reference_history.csv` is written at the same batches as `potential_history.csv`.

$$
\bar\phi_{\mathrm{top}}
=\frac{1}{N^2}\sum_{a,b}\phi(x_a,y_b,z_{\mathrm{high}}^-),\qquad
\phi_{\mathrm{rel}}(\mathbf r)=\phi(\mathbf r)-\bar\phi_{\mathrm{top}}.
$$

$N$ is `sim.injection_face_phi_grid_n`, and the samples are cell centered across the full periodic top face. Join the two
histories on `batch` and subtract `potential_mean_V` from element `potential_V` at that same batch.

```toml
[output]
history_stride = 1000
write_mesh_potential = true
write_potential_history = true
```

The plane mean is neither infinity potential nor plasma potential, and it does not feed back into solar-wind injection. It
removes a constant gauge but not dependence on box height, zero mode, or finite-image truncation. A large
`potential_std_V` or min/max span weakens the interpretation of z-high as one reservoir plane.

With either `write_mesh_potential=true` or `write_potential_history=true`, `summary.txt` also records top statistics for the
final state. Use those values to reference final `mesh_potential.csv`. If the final batch is not on the history stride, do
not reuse a top value from another batch.

## Meaning of finite images

`field_periodic_image_layers=N` adds source cells through shell $N$ around the primary cell.

| $N$ | Cells included |
| --- | --- |
| 0 | primary only, $1\times1$ |
| 1 | one surrounding shell, $3\times3$ |
| 2 | two surrounding shells, $5\times5$ |

With `field_periodic_far_correction="none"`, BEACH does not replace cells beyond that shell with an Ewald or cached operator.
Top-relative potential does not turn this into an infinite-periodic solution. Increase image layers until target quantities
stop changing. Use
[Infinite-periodic periodic2 with outer plasma](InfinitePeriodicOuterConfiguration.en.html) when physical nonzero and zero
modes require an infinite-periodic closure.

## Acceptance order

1. Check that `abs(neutral_return_weight_scale-1)` and `neutral_return_unresolved_fraction` are small.
2. When reducing `dt`, maintain or increase `max_step*dt`; also increase `rays_per_batch` to converge the return distribution.
3. Vary `batch_duration` by at least $T,T/2,T/4$ and compare charge, relative potential, and force histories.
4. Move z-high and increase `injection_face_phi_grid_n` to check top mean and variation.
5. Increase image shells through $N,N+1,N+2$.
6. Check solar-wind macro-particle count and random-seed uncertainty.

Full reflection is an artificial top mirror, not a self-consistent sheath or quasineutral solution. The result means:
“surface redistribution produced when net photoelectron current is set to zero under the specified local solar-wind flux.”

`sim.softening` and the old `[field]` table found in historical outputs have been removed from current input. Port an older
0p5x configuration from the current example rather than copying it literally.

## Scalar-barrier comparison

For a scalar-barrier comparison, remove the closed-PE settings and use:

```toml
[external_boundary.particles]
mode = "local_source"
inflow_model = "infinity_barrier"

[external_boundary.ordinary_open]
model = "potential_barrier"

[sim]
phi_infty = 0.0
```

This filters an upstream VDF with a face-average scalar and classifies reflection or escape from energy at each open crossing.
It is not forward and reverse transport through one continuous sheath profile. See
[reservoir inflow](ReservoirInjection.en.html) and
[particle escape and return](ParticleEscapeReturn.en.html) for equations and constraints.
