title: Understand the Zhao Stationary Closure

Lang: [English](ZhaoStationaryClosure.en.md) | [日本語](ZhaoStationaryClosure.md)

# Understand the Zhao Stationary Closure

`surface_current_model.model="zhao_stationary"` assumes a stationary sheath outside the BEACH domain and connects its
currents and velocity-space barriers to the top of the box. It does not evolve the outer sheath in time. BEACH solves
the zero-current root once at run startup, then applies the fixed total currents to hit and return distributions tracked
inside the box.

This page explains how the root, the 0 V reservoir, and photoelectron (PE) return and escape fit together.
See the [input-parameter reference](Parameters.en.html) for every key and constraint, and the
[output-format reference](OutputReference.en.html#zhao_stationary) for the complete receipt list.

## Problem solved at startup

The Zhao stationary closure assumes a planar, collisionless, unmagnetized outer sheath. It uses ambient electrons,
cold ions, and photoelectrons when PE is active to solve a Zhao Type A, B, or C zero-current root. The solution provides
the branch, surface potential $\phi_0$, potential minimum $\phi_m$, ambient-electron density, and species current densities.

BEACH does not resolve the root again during the run. The branch, $\phi_0$, $\phi_m$, and current targets remain fixed as
the surface charge changes between batches. This differs from
[quasistatic matching-plane coupling](MatchingPlaneCoupling.en.html), which iterates the outer response against the
outward fluxes in every batch.

### No PE means Type C

With `photoelectron_source_scale=0.0`, omit the PE species and PE-specific inputs. `zhao_branch="auto"` or `"c"` selects
Type C and solves

$$
J_e + J_i = 0.
$$

The result contains only electron and ion absorption targets and the z-high kinetic map. It creates no PE particles and
no PE emission, return, or escape targets.

### With PE, close the large channels separately

The PE source density follows from solar elevation $\alpha$, reference density $n_{pe,ref}$, and scale $s_{UV}$:

$$
n_{pe,0}=s_{UV}n_{pe,ref}\sin\alpha.
$$

The ion-species density is the ion density at infinity. The configured electron density and PE
`emit_current_density_a_m2` sample the raw particle distributions; the Zhao root determines the fixed current targets.
Current signs below denote contributions to surface charging.

| Channel | Sign | BEACH treatment |
|---|---:|---|
| Ambient-electron absorption $J_e$ | negative | Surface absorption target |
| Ion absorption $J_i$ | positive | Surface absorption target |
| PE emission $J_{emit}$ | positive | Emission-reaction target |
| PE return $J_{return}$ | nonpositive | Surface reabsorption target |
| PE escape $J_{escape}$ | positive | External-boundary target, not deposited on the surface |

The Zhao root supplies two balances:

$$
J_{return}=J_{escape}-J_{emit}\le0,
\qquad
J_e+J_i+J_{escape}=0.
$$

The PE channels and surface channels therefore close independently:

$$
J_{return}+J_{emit}-J_{escape}=0,
\qquad
J_e+J_i+J_{return}+J_{emit}=0.
$$

BEACH rescales emission and return separately instead of dividing by their small net PE current. The signed current
carried outward by escaping PE particles is $-A J_{escape}$ and is not added to a surface element.

## Map the 0 V reservoir to the current box top

The configured ambient Maxwell VDF represents a reservoir at the plasma potential at infinity, defined as 0 V. At the
start of each batch, BEACH evaluates the mean z-high potential $\phi_f$ from the current surface charge. It selects an
inflow tail that can reach both the outer access bottleneck and $\phi_f$, then maps the normal speed to the injection face by

$$
\frac12 m v_{n,f}^{2}=\frac12 m v_{n,\infty}^{2}-q(\phi_f-0).
$$

Inspect the potential variation across z-high when deciding whether this face-mean approximation is adequate.

| Branch | Ambient-electron access | Electron / PE outward barrier | Ion |
|---|---:|---:|---:|
| Type A | $\phi_m$ | $\phi_m$ | 0 V |
| Type B / C | 0 V | 0 V | 0 V |

For an outward z-high crossing, BEACH uses the local crossing potential rather than the face mean, together with the
particle's normal kinetic energy. An electron or PE that cannot reach the fixed barrier is specularly reflected at z-high;
only a particle with enough energy is classified as escape. The PE launch VDF remains the configured surface
half-Maxwellian.

## Fixed quantities and batch-dependent quantities

| Fixed at run startup | Re-evaluated or resampled each batch |
|---|---|
| Branch, $\phi_0$, $\phi_m$, ambient-electron density | Surface charge and field inside the BEACH domain |
| Signed species current densities and reference area | Mean z-high potential $\phi_f$ and selected inflow tail |
| Absorption, emission, and escape current targets | Trajectories, hit positions, and local return / escape classification |
| Branch-dependent access and barrier potentials | Target charge $I\Delta t$ and scaling of the raw distribution |

This separation enforces the stationary total currents, while the elementwise charging pattern remains a Monte Carlo
result of the trajectories in each batch.

## Validity domain

- BEACH does not solve the outer electric field, space charge, Debye shielding, turning-point distance, or flight time.
- Reflection at z-high is an adiabatic boundary contraction of the outer return trajectory.
- The uniform magnetic field must be zero because the closure is unmagnetized.
- A planar stationary solution does not respond self-consistently to curvature, collisions, outer transients, or plasma
  and illumination conditions that change during the run.
- A small zero-current residual does not establish convergence of the finite-particle hit and return maps.

See the [input-parameter reference](Parameters.en.html) for the complete species, boundary, charge, and temperature
contract. See [Photoelectron emission and return](PhotoelectronEmission.en.html) for the raycast reaction charge and VDF.

## Run an example

The complete PE case is `examples/periodic2_zhao_fixed_current.toml`.

```bash
beach examples/periodic2_zhao_fixed_current.toml
```

Use `examples/periodic2_zhao_no_photo_fixed_current.toml` for the no-PE Type C case. Both examples include the required
reservoir, open z-high face, and fixed-current species in addition to `[surface_current_model]`.

After completion, inspect `summary.txt` for the model, selected branch, and two budget residuals. Then use
`charge_ledger.csv` to compare raw, target, and applied charge, hit counts, and `fixed_*_weight_scale`. For a PE case,
vary ray count, batch width, and RNG seed until the elementwise return distribution converges. A completed run and a
closed zero-current budget do not by themselves validate the outer-sheath approximation or the spatial distribution.

The implementation-level definitions are in [SPEC section 7.7](../SPEC.md#77-自動表面電流model) and
[`bem_surface_current_model.f90`](../src/physics/sheath/bem_surface_current_model.f90).
