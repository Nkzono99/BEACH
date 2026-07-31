title: Finite-image periodic2 configuration

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# Finite-image periodic2 configuration

This page is the canonical `periodic2` setup for studying regolith charging and photoelectron surface redistribution with a
boundary reservoir + closed PE. The solar-wind VDF is defined at z-high, only photoelectrons are closed, and potential is read
relative to the z-high plane mean. BEACH does not solve a self-consistent outer-plasma/sheath model.

## Distinguish the two configurations first

| Configuration | Solar-wind inflow | Photoelectrons | Potential reference | Use |
| --- | --- | --- | --- | --- |
| **Boundary reservoir + closed PE** | uncorrected VDF at z-high | z-high reflection + `neutral_return` | z-high plane mean | baseline on this page; surface redistribution and batch sensitivity |
| Scalar barrier | `infinity_barrier` | common `potential_barrier` | `phi_infty` | comparison using one scalar barrier |

Do not stack these alternatives in one run. In particular, do not add
`infinity_barrier` or `potential_barrier` to the closed-PE baseline.

## Recommended integrated configuration

The complete runnable example is
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml).

```toml
[sim]
field_solver = "fmm"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"

[domain]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0e-4, 1.0e-4, 1.0e-3]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"

[reservoir]
inflow_model = "source_vdf"
face_potential_grid_n = 5
```

`[domain]` owns the cell and periodic topology; `[field_boundary]` owns field closure. The current `periodic2` mode accepts
only `periodic_axes=["x", "y"]` with a nonperiodic z axis. Periodicity is common to all species and cannot be overridden by a
species or `[particle_boundary]`.
The former `[external_boundary]` table has been removed. Configure external-reservoir conditions in `[reservoir]`,
per-species inflow faces in `[particles.species.boundary_inflow]`, and outward particle actions on nonperiodic faces in
`[particle_boundary]`.

Use `boundary_inflow="reservoir"` on z-high for solar-wind electrons and ions. Density, temperature, and drift describe a
**local boundary VDF on z-high**; surface potential does not filter accessibility or map speed.

Only the photoelectron species gets the closed surface-charge closure:

```toml
[[particles.species]]
species_key = "photoelectron"
source_mode = "photo_raycast"
deposit_opposite_charge_on_emit = true
surface_charge_closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"
```

The species-level `z_high="reflect"` reverses only normal velocity at nonperiodic z-high. Solar-wind species retain `"inherit"` and
follow `particle_boundary.z_high="open"` with `ordinary_open_model="escape"`.

This baseline preserves the in-plane position of the boundary event. For a sensitivity comparison that uniformly redistributes
return positions over the x-y plane, replace it with `z_high="redistributed_reflect"`. The velocity action is unchanged; only
position is resampled. This is an alternative return-destination model for closed PE, not a self-consistent sheath. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for context and simultaneous-event rules.

## What one batch does

1. Build a finite-image field snapshot from batch-start surface charge.
2. Inject solar wind from the z-high local VDF and photoelectrons from ray hits.
3. Track all particles in the frozen snapshot, reflecting only photoelectrons at z-high.
4. Close unresolved photoelectrons from the resolved return distribution and commit all charge deltas.
5. Optionally write element potential and z-high mean from the same post-commit charge.

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

The batch is rejected when:

- emission is nonzero but resolved return charge is zero;
- a photoelectron actually escapes or is `soft_discard`ed;
- a nonfinite value or charge-sign mismatch occurs.

Only particles still alive at `max_step` are statistically closed. The raw unresolved charge, correction, scale, and
unresolved fraction remain separate in `charge_ledger.csv`.

## Read potential relative to the z-high plane

With `[domain]`, `output.write_files=true`, `output.write_potential_history=true`, and `output.history_stride>0`,
`top_reference_history.csv` is written at the same batches as `potential_history.csv`.

$$
\bar\phi_{\mathrm{top}}
=\frac{1}{N^2}\sum_{a,b}\phi(x_a,y_b,z_{\mathrm{high}}^-),\qquad
\phi_{\mathrm{rel}}(\mathbf r)=\phi(\mathbf r)-\bar\phi_{\mathrm{top}}.
$$

$N$ is `reservoir.face_potential_grid_n`, and the samples are cell centered across the full periodic top face. Join the two
histories on `batch` and subtract `potential_mean_V` from element `potential_V` at that same batch.

```toml
[output]
history_stride = 1000
write_mesh_potential = true
write_potential_history = true
```

The plane mean is neither infinity potential nor plasma potential, and it does not feed back into solar-wind injection. It
removes a constant gauge but not dependence on box height, zero mode, or finite-image truncation. A large `potential_std_V` or
min/max span weakens the interpretation of z-high as one reservoir plane.

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
stop changing. Use `field_periodic_far_correction="cached_kneq0"` when the infinite-periodic nonzero mode is required, and
follow the convergence procedure in
[periodic2 electrostatics](PeriodicElectrostatics.en.html).

## Acceptance order

1. Check that `abs(neutral_return_weight_scale-1)` and `neutral_return_unresolved_fraction` are small.
2. When reducing `dt`, maintain or increase `max_step*dt`; also increase `rays_per_batch` to converge the return distribution.
3. Vary `batch_duration` by at least $T,T/2,T/4$ and compare charge, relative potential, and force histories.
4. Move z-high and increase `reservoir.face_potential_grid_n` to check top mean and variation.
5. Increase image shells through $N,N+1,N+2$.
6. Check solar-wind macro-particle count and random-seed uncertainty. With `redistributed_reflect`, also vary particle count and
   seed for the additional return-position sampling.

Full reflection is an artificial top mirror, not a self-consistent sheath or quasineutral solution. The result describes
surface redistribution with zero net photoelectron current under the specified boundary solar-wind flux.

## Scalar-barrier comparison

For a scalar-barrier comparison, remove the closed-PE settings and combine `reservoir.inflow_model="infinity_barrier"`,
`reservoir.phi_infty`, and `particle_boundary.ordinary_open_model="potential_barrier"`. This comparison filters an upstream VDF
with a face-average scalar and classifies reflection or escape from energy at each open crossing. See
[reservoir inflow](ReservoirInjection.en.html) and [particle escape and return](ParticleEscapeReturn.en.html) for configuration,
equations, and constraints.
