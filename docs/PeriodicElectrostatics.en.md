title: periodic2 electrostatics

Lang: [日本語](PeriodicElectrostatics.md) | [English](PeriodicElectrostatics.en.md)

# periodic2 electrostatics

The `field_bc_mode="periodic2"` field is the sum of near finite images, laterally varying infinite-periodic modes, the
plane-average `k=0` component, and outer-plasma response. In an x/y-periodic, z-open slab, one computation path owns each of
these four components and adds it exactly once to the electrostatic snapshot.

## Decompose the field into four components

| Component | Physical meaning | Production path |
| --- | --- | --- |
| Primary and near images | Strong local field near the primary cell | Finite-image Direct/FMM sum |
| Far `k\ne0` | Infinite-periodic far field varying in x/y | `cached_kneq0` operator |
| Surface `k=0` | Plane-average field from total charge below each height | Triangle-height cumulative polynomial |
| Plasma response | Outer-plasma response of zero and nonzero modes | Selected outer model |

`cached_kneq0` returns only the nonzero modes; the snapshot adds the physical boundary-conditioned surface `k=0` exactly once.
This separates lateral infinite-periodic correction from z-direction boundary and sheath choices.

See [Finite-image periodic2 configuration](FinitePeriodicConfiguration.en.html) and
[Infinite-periodic periodic2 with outer plasma](InfinitePeriodicOuterConfiguration.en.html) for complete combinations.

## Select the computation path for nonzero modes

| Nonzero-mode construction | Use | Constraint |
| --- | --- | --- |
| Legacy finite images | Small comparison and finite-image model | Contains nothing outside the image range |
| `panel_spectral_reference` | Small triangle-P0 reference | Direct, zero softening, and mode/quadrature convergence |
| `cached_kneq0` | Infinite-periodic nonzero modes for FMM production | x/y periodic, z nonperiodic, and `exclude_k0` |
| `m2l_root_oracle` | Far-correction diagnostics | Not a production particle hot path |

`field_periodic_far_correction="auto"` currently has compatibility behavior equivalent to `none`. Infinite-periodic production
must select `cached_kneq0` explicitly. Startup validation checks consistency between high-level settings and typed `[periodic2]`
configuration and rejects contradictory zero-mode ownership or unsupported outer models.

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

The direct periodic Coulomb sum converges slowly, and a nonneutral slab has no unique plane-average field without a z-boundary
closure. Ewald splitting uses numerical parameter $\alpha$:

$$
\frac1r=
\frac{\operatorname{erfc}(\alpha r)}r
+\frac{\operatorname{erf}(\alpha r)}r.
$$

The first term decays rapidly in real space; the smooth second term is expanded in reciprocal modes. Ewald2P transforms only x/y
to the reciprocal lattice and leaves z open.

| Ewald part | Evaluation |
| --- | --- |
| Real space | Finite-image sum of screened Coulomb terms |
| Reciprocal `k\ne0` | Finite layers of x/y Fourier modes |
| `k=0` | Separate plane-average term along open z |

$\alpha$ distributes numerical work between real and reciprocal space; it is neither a Debye length nor physical screening rate.
A sufficiently converged result should not depend on it. The implemented Ewald reference is a high-accuracy teacher at configured
real and reciprocal cutoffs, not an analytically exact infinite sum.

## Turn the far correction into a reusable operator

The cold build samples the residual between Ewald2P teacher and finite image shell,

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N),
$$

using proxy sources and check points. With fixed geometry, periodic length, FMM order, and target topology, the map from root source
multipole to far local expansion is linear:

$$
\mathbf L_t^\mathrm{far}=\mathbf A_t\mathbf M_\mathrm{root}.
$$

$\mathbf A_t$ is fit once by QR and cached. When charge changes, the operator is not refit; it is applied to current
$\mathbf M_\mathrm{root}$. Ewald sums occur only on the cold teacher path, never as an all-source sum during warm particle
evaluation.

## Add the physical `k=0` component exactly once

The full Ewald residual used for fitting contains a symmetric-vacuum `k=0` required by the operator construction. Physical zero
mode, however, is selected by `lower_boundary_model` or an outer model. The FMM core reconstructs symmetric `k=0` from the same
source state and subtracts it from the cached result:

$$
K_{k\ne0}=K_\mathrm{shell}+R_\mathrm{Ewald}^{\mathrm{full}}-K_0^\mathrm{sym}.
$$

The snapshot then adds selected $K_0^\mathrm{physical}$:

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}.
$$

`zero_mode_policy="exclude_k0"` does not discard the mean field; it assigns ownership so the nonzero backend does not add it
twice.

## Reuse a cold-built operator across batches

| Phase | Performed | Not performed |
| --- | --- | --- |
| Cache miss | Place proxy/check points, evaluate Ewald teacher, QR fit, publish with checksum | Particle tracking |
| Cache hit | Validate fingerprint, shape, and checksum; load operator | Ewald reevaluation and refit |
| Charge refresh | P2M/M2M, apply cached matrix, refresh zero-mode state | Rebuild operator |
| Particle evaluation | Near terms, local expansion, subtract cached symmetric `k=0`, add physical `k=0` | All-source Ewald sum |

The fingerprint includes geometry, periodic length, FMM order, image/Ewald layers, source kernel, target topology, and generator and
build versions. A mismatched cache is never reused approximately. Cold and warm results must agree to roundoff.

## Build `k=0` from surface charge

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

An outer model uses surface zero-mode field and interface conditions to form plasma response. See
[Kinetic 1-D outer plasma](KineticOuterPlasma.en.html) for the split kinetic model and
[Unified linear response](UnifiedLinearResponse.en.html) for the surface-to-far linear model.

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
