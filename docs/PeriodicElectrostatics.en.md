title: periodic2 electrostatics

Lang: [日本語](PeriodicElectrostatics.md) | [English](PeriodicElectrostatics.en.md)

# periodic2 electrostatics

For an x/y-periodic, z-open slab, `field_bc_mode="periodic2"` computes finite images, infinite-periodic `k\ne0`, and
plane-average `k=0` separately and adds each component exactly once.

## Decompose the field into three components

| Component | Physical meaning | Production path |
| --- | --- | --- |
| Primary and near images | Strong local field near the primary cell | Finite-image Direct/FMM sum |
| Far `k\ne0` | Infinite-periodic far field varying in x/y | `cached_kneq0` operator |
| Surface `k=0` | Plane-average field from total charge below each height | Triangle-height cumulative polynomial |

`cached_kneq0` returns only the nonzero modes; field composition adds the physical boundary-conditioned surface `k=0` exactly
once. See [Finite-image periodic2 configuration](FinitePeriodicConfiguration.en.html) for a complete run configuration.

## Select the computation path for nonzero modes

| Nonzero-mode construction | Use | Constraint |
| --- | --- | --- |
| Finite images | Small comparison and finite-image model | Contains nothing outside the image range |
| `panel_spectral_reference` | Small triangle-P0 reference | Direct and mode/quadrature convergence |
| `cached_kneq0` | Infinite-periodic nonzero modes for FMM production | x/y periodic, z nonperiodic, and `exclude_k0` |

`field_periodic_far_correction="auto"` has compatibility behavior equivalent to `none`. Infinite-periodic production must
select `cached_kneq0` explicitly. Startup validation checks high-level settings against typed `[periodic2]` configuration and
rejects contradictory zero-mode ownership. Removed `[outer_plasma]` and `[coupling]` tables are unknown input.
`m2l_root_oracle` has been removed and is rejected at startup.

The Ewald2P teacher, root-multipole-to-local operator, cache, and FMM-state connection are documented separately in
[periodic2 Far Correction](PeriodicFarCorrection.en.html). This page treats that implementation as the nonzero component of the
complete field.

## Define the range included by a finite-image sum

Sources are added explicitly for primary cell and configured image layer $N$,

$$
(i,j)\in[-N,N]^2.
$$

This evaluates near interactions with the original kernel but omits the smooth field of infinitely many images outside the range.
Therefore `field_periodic_far_correction="none"` is a finite-image model and must not be interpreted as an infinite-periodic
solution without convergence as $N$ increases.

The near-image layer in FMM must match the shell subtracted when fitting the far operator. The cache fingerprint includes image
layer for this reason.

## Separate the infinite-periodic far field with Ewald2P

`cached_kneq0` applies the difference between an Ewald2P teacher and finite-image shell as an operator and removes the
teacher's symmetric `k=0`. Field composition adds the selected physical `k=0` exactly once:

$$
K_\mathrm{surface}
=\left(K_\mathrm{shell}+R_\mathrm{Ewald}^{\mathrm{full}}-K_0^\mathrm{sym}\right)
+K_0^\mathrm{physical}.
$$

The expression in parentheses belongs to the nonzero backend. `zero_mode_policy="exclude_k0"` prevents double counting; it
does not discard the mean field. See [periodic2 Far Correction](PeriodicFarCorrection.en.html) for Ewald splitting, operator
fitting, the FMM insertion point, and cache lifecycle.

## Add the physical `k=0` component exactly once

For triangle total charge $q_i$, let $F_i(z)$ be the fraction of its area at or below height $z$. Plane-average cumulative charge is

$$
C(z)=\sum_iq_iF_i(z).
$$

$F_i$ is piecewise quadratic between the three vertex heights. A geometry plan sorts all vertex heights into breakpoints and
stores quadratic coefficients contributed by each triangle in each interval. Horizontal triangles are separate sheet charges;
evaluation on a sheet distinguishes minus trace, plus trace, and principal value.

When $q_i$ changes, stored geometry coefficients are multiplied by charge into interval differences, then a prefix sum forms

$$
C(z)=a_0+a_1z+a_2z^2
$$

and its primitive. Geometry plan is not rebuilt.

## Derive zero-mode field and potential from Gauss's law

For cell area $A=L_xL_y$ and lower far field $E_\mathrm{bottom}$, Gauss's law gives

$$
E_0(z)=E_\mathrm{bottom}+\frac{C(z)}{\epsilon_0A}.
$$

Potential from gauge point $(z_g,\phi_g)$ is

$$
\phi_0(z)=\phi_g-E_\mathrm{bottom}(z-z_g)
-\frac1{\epsilon_0A}\int_{z_g}^zC(\zeta)\,d\zeta.
$$

A binary search locates the breakpoint interval and evaluates its quadratic and cubic primitive, giving $O(\log N_z)$ per point.

A nonneutral cell can retain constant far field and linearly growing potential in z. Zero mode is a physical Gauss-law component,
not a numerical term that may be deleted.

## Close the mean field with a z-boundary condition

For total surface charge $Q=\sum_iq_i$, current choices are:

| `lower_boundary_model` | $E_\mathrm{bottom}$ | $E_\mathrm{top}$ | Meaning |
| --- | ---: | ---: | --- |
| `symmetric_vacuum` | $-Q/(2\epsilon_0A)$ | $+Q/(2\epsilon_0A)$ | Equal vacuum half-spaces above and below with no external field |
| `e_bottom_zero` | $0$ | $Q/(\epsilon_0A)$ | Legacy closure fixing lower flux to zero |

Neither solves dielectric screening or polarization. `symmetric_vacuum` is the minimal symmetric closure without an additional
interface or permittivity; `e_bottom_zero` exists for legacy reproduction and is not a universal physical default.

## Search periodic images reachable by particle trajectories

Field targets are wrapped into the primary periodic cell, while physical trajectory-event positions are retained. Mesh collision
searches geometric periodic images required by a segment. Field near-image layer and collision image bound are distinct numerical
concepts and do not share one fixed count. See [Particle collision and boundary events](ParticleEvents.en.html).

## Converge each component separately

- Increase image layers and converge the observable for a finite-image model.
- For cached models, compare cache miss/hit and thread/MPI configurations.
- Vary Ewald $\alpha$, real/reciprocal layers, and proxy/check settings to bound teacher/operator error.
- Compare against an oracle to exclude duplicate primary, near, far, symmetric `k=0` subtraction, or physical `k=0` addition.
- Inspect Gauss residual and lower/upper closure.
- Do not interpret a finite-height potential difference in a nonneutral cell directly as escape energy to infinity.

See [FMM internals](FMMCore.en.html) for internal Ewald formulas and operator APIs.

## Code reference

- Periodic FMM plan, state, and evaluation: [`bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90)
- Ewald teacher and cached root operator: [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
- Cached symmetric `k=0` subtraction: [`bem_coulomb_fmm_eval_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90)
- Surface zero-mode plan/state: [`bem_periodic_zero_mode_plan.f90`](../src/physics/periodic_zero_mode/bem_periodic_zero_mode_plan.f90)
- Zero-mode evaluation: [`bem_periodic_zero_mode_eval.f90`](../src/physics/periodic_zero_mode/bem_periodic_zero_mode_eval.f90)
- Component ownership and snapshot composition: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
