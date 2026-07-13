title: periodic2 Zero Mode and Outer Plasma

Lang: [English](PeriodicZeroModeOuterPlasma.en.md) | [日本語](PeriodicZeroModeOuterPlasma.md)

# periodic2 Zero Mode and Outer Plasma

This document explains how the two-periodic, one-open electrostatic problem is
split into laterally varying `k!=0` modes, a plane-averaged `k=0` mode, and an
outer-plasma profile. See [FMMCore](FMMCore.en.md) for FMM expansions,
[FieldSolvers](FieldSolvers.en.md) for configuration choices, and ADR 0001/0002
for the physical applicability contracts.

## 1. End-to-end view

Let x/y be periodic and z be open. A lateral Fourier decomposition gives three
numerical responsibilities:

| Component | Meaning | Production implementation |
| --- | --- | --- |
| `k!=0` | lateral charge variations | finite-image FMM plus the `cached_kneq0` far operator |
| surface `k=0` | plane-averaged field from charge below each height | triangle-height cumulative polynomials |
| outer profile | plasma response from the interface to infinity | `kinetic_1d` or `unified_linear_response` |

`cached_kneq0` is not a complete field solver by itself. It deliberately
excludes the surface `k=0` mode. `electrostatic_snapshot` adds the physical zero
mode exactly once. This separates the infinite-periodic lateral correction from
the open-axis boundary condition and sheath model.

## 2. What `cached_kneq0` stores

### 2.1 Why a finite image shell is insufficient

The ordinary periodic FMM explicitly adds the primary cell and image shifts

$$
(i,j)\in[-N,N]^2.
$$

This resolves the near field but omits the smooth field of the infinitely many
images outside the shell. Evaluating an all-source Ewald sum for every particle
would recover it, but is too expensive for the hot path.

### 2.2 What the Ewald2P reference means

A direct sum of periodic Coulomb $1/r$ images converges slowly, and a non-neutral
slab also requires an explicit mean-field boundary condition. Ewald splitting
introduces numerical parameter $\alpha$:

$$
\frac{1}{r}=\frac{\operatorname{erfc}(\alpha r)}{r}
+\frac{\operatorname{erf}(\alpha r)}{r}.
$$

The first term decays rapidly in real space. The second is smooth and can be
represented by relatively few reciprocal modes. BEACH Ewald2P transforms only
x/y into reciprocal space and leaves z as an open coordinate.

| Part | BEACH evaluation |
| --- | --- |
| real space | finite image sum of screened Coulomb terms |
| reciprocal `k!=0` | finite layer of x/y reciprocal-lattice modes |
| `k=0` | separate open-z mean-field term |

$\alpha$ is a **numerical work-splitting parameter**, not a Debye length or a
physical screening rate. A converged result is independent of $\alpha$. What the
implementation calls `exact Ewald` is the high-accuracy teacher evaluated at the
configured real/reciprocal cutoffs, not an analytic infinite-sum value.

During a cold build, `cached_kneq0` samples only the difference

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N).
$$

Ewald2P is therefore a teacher for building the operator. It is absent from warm
particle evaluation. The implemented real-space, reciprocal-space, and zero-mode
formulas are given in [FMMCore 8.2](FMMCore.en.md#82-periodic2-ewaldewald2p-correction).

### 2.3 Turn the far correction into a linear operator

For fixed geometry, periods, FMM order, and target topology, the map from source
multipoles to far local expansions is linear. During a cold build, proxy sources
and check points fit the difference between the Ewald reference and the finite
image shell. The cache stores this matrix, not charge state or particle positions:

$$
\mathbf L_t^{\mathrm{far}}=\mathbf A_t\mathbf M_{\mathrm{root}}.
$$

$\mathbf M_{\mathrm{root}}$ is rebuilt from current charges, $\mathbf A_t$ is the
cached operator for a target anchor, and $\mathbf L_t^{\mathrm{far}}$ is its local
expansion. Charge changes therefore reuse $\mathbf A_t$.

### 2.4 Why subtract `k=0`

The full-periodic Ewald teacher contains a convenient symmetric-vacuum `k=0`
component. The physical zero mode, however, must be selected by
`lower_boundary_model` and the outer-plasma closure. The FMM core reconstructs
the same symmetric zero mode from current source state and subtracts it:

$$
K_{k\ne0}=K_{\mathrm{shell}}+R_{\mathrm{Ewald}}^{\mathrm{full}}-K_0^{\mathrm{sym}}.
$$

The snapshot then adds the selected physical zero mode:

$$
K_{\mathrm{surface}}=K_{k\ne0}+K_0^{\mathrm{physical}}.
$$

Thus `exclude_k0` means that the nonzero backend must not double count the mean
field. It does not mean that the physical mean field is discarded.

### 2.5 Cold and warm paths

| Phase | Work performed | Work avoided |
| --- | --- | --- |
| cache miss | proxy/check construction, Ewald teacher, QR fit, checked atomic publish | particle tracking |
| cache hit | fingerprint, shape, and checksum validation; operator load | Ewald evaluation and refit |
| charge refresh | P2M/M2M, cached matrix application, zero-mode state update | operator generation |
| particle evaluation | local expansion, near direct, cached symmetric-`k=0` subtraction, physical-`k=0` addition | all-source Ewald sum |

The cache identity includes geometry, periods, order, image/Ewald layers, source
kernel, target topology, and generator/build versions. A near-match is never
accepted as a reusable operator.

## 3. Surface `k=0` algorithm

### 3.1 Project triangles onto cumulative height charge

For total P0 panel charge $q_i$, let $F_i(z)$ be the fraction of triangle area
at or below height z. The plane-averaged cumulative charge is

$$
C(z)=\sum_i q_iF_i(z).
$$

$F_i$ is piecewise quadratic between its three vertex heights. Plan construction
sorts all vertex heights into breakpoints and records each triangle's quadratic
contribution. Horizontal panels are stored as sheet charges so minus, plus, and
principal-value traces remain distinct.

### 3.2 Charge refresh

When $q_i$ changes, charge-weighted coefficients are accumulated through an
interval-difference array. A prefix sum forms the interval polynomial
$C(z)=a_0+a_1z+a_2z^2$, and its primitive is precomputed at every breakpoint.
The geometry plan is not rebuilt.

### 3.3 Field and potential

For cell area $A=L_xL_y$ and lower far field $E_{\mathrm{bottom}}$, Gauss' law gives

$$
E_0(z)=E_{\mathrm{bottom}}+\frac{C(z)}{\epsilon_0A}.
$$

With gauge $(z_g,\phi_g)$,

$$
\phi_0(z)=\phi_g-E_{\mathrm{bottom}}(z-z_g)
-\frac{1}{\epsilon_0A}\int_{z_g}^{z}C(\zeta)\,d\zeta.
$$

Evaluation uses binary search for the interval and polynomial evaluation inside
it, giving $O(\log N_z)$ work per point.

### 3.4 Lower-boundary closure

For total surface charge $Q=\sum_iq_i$:

| Model | $E_{\mathrm{bottom}}$ | $E_{\mathrm{top}}$ | Meaning |
| --- | ---: | ---: | --- |
| `symmetric_vacuum` | $-Q/(2\epsilon_0A)$ | $+Q/(2\epsilon_0A)$ | equal vacuum half-spaces and no imposed field |
| `e_bottom_zero` | $0$ | $Q/(\epsilon_0A)$ | legacy closure fixing lower flux to zero |

Neither option solves dielectric screening inside regolith. `symmetric_vacuum`
is the minimal symmetric closure without an interface location or permittivity;
`e_bottom_zero` is a legacy-reproduction option rather than the physical default.

## 4. Coupling to outer-plasma models

### 4.1 Split and unified models are different

In a split model the local and outer domains meet at an interface:

| Region | Field used |
| --- | --- |
| mesh to interface | `k!=0` surface field plus surface `k=0` |
| above interface | 1D outer profile; the split contract requires lateral modes to have decayed sufficiently |

The outer profile is therefore not added throughout the local domain. Potential
and normal field match at the interface, and only particles handed to the outer
domain use the same profile for return/escape decisions.

`unified_linear_response` is different: one zero-mode Poisson grid spans the
surface projection to the far boundary, while nonzero modes are continued into
plasma tails. The particle interface is only an ownership boundary.

See [Outer Sheath and Reservoir Particle Boundaries](SheathReservoirBoundary.en.md)
for inflow acceleration/deceleration, outgoing turning/return, and the distinction
from Zhao injection corrections.

### 4.2 Numerical path by model

| `outer_plasma.model` | Zero-mode treatment | Main applicability |
| --- | --- | --- |
| `none` | surface `k=0` with the selected boundary closure | no outer plasma |
| `linear_debye` | $\Delta\phi\exp[-(z-z_I)/\lambda_D]$ | small-amplitude split reference |
| `kinetic_1d` | nonlinear 1D Poisson solve with VDF closures | monotone, collisionless, unmagnetized sheath |
| `unified_linear_response` | rough-surface and plasma sources in one 1D Poisson grid | linear response without a valid split window |

The production `cached_kneq0` path currently accepts `none`, `kinetic_1d`, and
`unified_linear_response`. `linear_debye` belongs to the small-system
`panel_spectral_reference` path.

### 4.3 `linear_debye`

For interface $z_I$ and $\Delta\phi=\phi_I-\phi_\infty$,

$$
\phi(z)=\phi_\infty+\Delta\phi e^{-(z-z_I)/\lambda_D},\qquad
E(z)=\frac{\Delta\phi}{\lambda_D}e^{-(z-z_I)/\lambda_D}.
$$

The model is rejected when $|\Delta\phi|/V_T$ exceeds its linearity limit.

### 4.4 `kinetic_1d`

#### 4.4.1 Domain and unknown

`kinetic_1d` represents the plane-averaged region from the rough-surface
interface $z=z_I$ to the infinity reservoir as a one-dimensional electrostatic,
collisionless, unmagnetized plasma. Its unknown is the mean outer potential
$\phi(z)$. The surface zero mode supplies interface field $E_I$ through

$$
-\phi'(z_I)=E_I,
$$

and the solver constructs $\rho(\phi)$ from ambient electrons, cold drifting
ions, and optional photoelectron flux before solving

$$
\frac{d^2\phi}{dz^2}=-\frac{\rho(\phi)}{\epsilon_0}.
$$

The infinity gauge is $\phi_\infty=0$. Interface potential $\phi_I$ is not an
input: it is the value $\phi(z_I)$ that satisfies the interface field, VDF
closures, and far condition together.

This solve does not impose floating balance $J_\mathrm{total}=0$ as a root
equation. Batch updates of surface charge change $E_I$, which changes $\phi_I$
and the species currents. Electron, ion, photoelectron, and external-circuit
current densities are diagnostics evaluated after the profile has converged.

#### 4.4.2 Density closures from VDFs

The first negative and positive z-high `reservoir_face` species define the
infinity ambient electrons and ions. With
`photoelectron_closure="kinetic_mean"`, the first negative `photo_raycast`
species supplies temperature and emitted current density for a plane-averaged
photoelectron source.

| Population | Infinity/surface inputs | Outer density construction |
| --- | --- | --- |
| ambient electron | $n_{e,\infty},T_e,q_e,m_e$ | map a half-Maxwellian by total-energy conservation, including absorbed and potential-reflected trajectories |
| ion | $n_{i,\infty},T_i,q_i,m_i,u_{i,\infty}$ | cold beam with speed and density from energy and flux conservation |
| photoelectron | $T_{pe},q_{pe},m_{pe},\Gamma_{pe,0}$ | surface half-Maxwellian with mean outgoing and post-turning return populations |

The cold-ion closure is

$$
u_i(z)=\sqrt{u_{i,\infty}^2-\frac{2q_i\phi(z)}{m_i}},\qquad
n_i(z)=n_{i,\infty}\frac{u_{i,\infty}}{u_i(z)}.
$$

A profile for which the square root ceases to be real is rejected because ions
cannot access the interface. The photoelectron escape fraction is

$$
f_{pe,\mathrm{esc}}=
\exp\left[-\frac{\max\{0,q_{pe}(\phi_\infty-\phi_I)\}}{T_{pe}}\right],
$$

and the remainder is the return population. The Poisson source is

$$
\rho(\phi)=q_en_e(\phi)+q_in_i(\phi)+q_{pe}n_{pe}(\phi).
$$

Temperatures are represented internally in joules. Each analytic closure also
returns $\partial n_s/\partial\phi$ and, where required,
$\partial n_s/\partial\phi_I$ for the Newton Jacobian.

The outgoing/returning density in `kinetic_mean` is a stationary closure for
outer space charge. It does not replace surface deposition by tracked particles
or add a second statistical return current. Individual escape and return are
handled by `kinetic_1d_profile_return` using the same converged profile; see
[Outer Sheath and Reservoir Particle Boundaries](SheathReservoirBoundary.en.md).

#### 4.4.3 Grid and boundary conditions

The interior uses conservative finite-volume residuals on a stretchable
nonuniform 1D grid. Current runtime values are listed below; only `debye_length`
is currently exposed as an individual input parameter.

| Item | Current value |
| --- | ---: |
| grid points | 128 |
| domain length | $10\lambda_D$ |
| grid stretch | 2 |
| maximum Newton iterations | 40 |
| residual tolerance | $10^{-8}$ |

Beyond the finite upper endpoint $L$, the model assumes an exponential tail
with $\lambda_\mathrm{tail}=\lambda_D$ and imposes

$$
\phi'(L)+\frac{\phi(L)}{\lambda_{\mathrm{tail}}}=0.
$$

This Robin condition represents exponential relaxation toward the infinity
gauge without forcing potential abruptly to zero at $L$. The remaining tail is
also used in particle flight-time integration.

#### 4.4.4 Nonlinear solve and acceptance

Analytic density derivatives form a bordered-tridiagonal Jacobian, so one Newton
step is $O(N_z)$. Dependence of interior densities on
$\phi_I=\phi_1$ produces the border column in addition to the ordinary
tridiagonal stencil. Backtracking keeps trial steps on the monotone branch.
Pseudo-transient regularization is used when ordinary Newton steps stall. When
a previous-batch profile is available, interface-field continuation advances
from its old field to the current value and halves a failed increment.

These methods alter only the convergence path. Final acceptance requires all of
the following:

| Condition | Meaning |
| --- | --- |
| original Poisson residual | the unregularized discrete Poisson residual is below tolerance |
| monotone branch | remain on the supported electron-repelling profile |
| ion accessibility | $u_i^2(z)>0$ at every grid point |
| kinetic Bohm entry | $u_{i,\infty}\ge\sqrt{(T_e+\gamma_iT_i)/m_i}$ |
| infinity quasineutrality | $q_en_{e,\infty}+q_in_{i,\infty}\simeq0$ |

Nonmonotone virtual-cathode profiles, trapped populations, and sub-Bohm ion
inflow are outside this model. A numerical iterate is not accepted when a
physical condition fails. Status distinguishes `not_applicable`,
`no_physical_solution`, and `numerical_failure`; the solver does not fall back
to another sheath or an unvalidated previous-batch solution.

#### 4.4.5 Batch update, MPI, and output

On batches selected by `outer_update_stride`, the profile is refreshed using
interface field reconstructed from committed surface charge. A previous profile
is a Newton initial guess only for the same model identity and grid. The MPI root
performs the 1D solve and broadcasts status, profile, and current diagnostics.
All particles in a batch share the updated immutable snapshot, so the profile is
not solved again after each impact.

Converged $z,\phi,E,\rho$ are written to `outer_plasma_profile.csv` and can seed
restart. Inspect `interface_potential`, `interface_field`,
`outer_integrated_charge`, species and total current densities, Newton iteration
count, and the original nonlinear residual. The quasistatic scope of immediate
return is documented under
[Why return is immediate and where the approximation applies](SheathReservoirBoundary.en.md#why-return-is-immediate-and-where-the-approximation-applies).

### 4.5 `unified_linear_response`

For plasma-accessible fraction $f_{\mathrm{access}}(z)$ and
$\kappa=1/\lambda_D$, the same nonuniform 1D grid receives the plane-averaged
surface source and

$$
\rho_{\mathrm{plasma}}(z)=-\epsilon_0f_{\mathrm{access}}(z)\kappa^2\phi(z).
$$

A tridiagonal Poisson system includes the bottom field and far Robin condition.
Above the response start, each nonzero mode continues with
$\alpha=\sqrt{k^2+\kappa^2}$ while preserving potential, normal field, and
tangential field. Violations of linearity, height-field geometry, accessible-area
convergence, or mode truncation fail closed.

## 5. Per-batch order

| Order | Operation |
| ---: | --- |
| 1 | update FMM multipoles and surface zero-mode state from new `q_elem` |
| 2 | solve the outer profile on configured batches; reuse a previous profile only as an initial guess with the same identity |
| 3 | set the zero-mode potential gauge to match outer interface potential |
| 4 | update Gauss residual and interface diagnostics |
| 5 | compose nonzero, zero, and prescribed fields exactly once according to ownership |

All particle evaluations within a batch use the same immutable snapshot. The
operator and outer profile are updated after batch charge commit, not after each
particle impact.

## 6. Diagnostics and failure

Important outputs include `interface_potential`, `interface_field`,
`gauss_residual`, `outer_integrated_charge`, nonlinear residual, linearity ratio,
and cache fingerprint/hit state.

The surface flux entering the upper domain is

$$
Q_{\mathrm{upper\ flux}}=Q+\epsilon_0AE_{\mathrm{bottom}},
$$

and the Gauss diagnostic adds the integrated outer charge. Failure to satisfy
physical applicability such as the Bohm condition, monotonicity, or linearity is
fail-closed even if a numerical iterate could be produced.

## 7. Implementation map

| Operation | Fortran implementation |
| --- | --- |
| cached nonzero operator | `fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90` |
| cached symmetric-`k=0` subtraction | `fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90` |
| surface zero-mode plan/state | `periodic_zero_mode/bem_periodic_zero_mode_plan.f90` |
| zero-mode evaluation | `periodic_zero_mode/bem_periodic_zero_mode_eval.f90` |
| field ownership/composition | `bem_electrostatic_snapshot.f90` |
| linear outer model | `outer_plasma/bem_outer_plasma_linear.f90` |
| kinetic outer model | `outer_plasma/bem_outer_plasma_kinetic.f90` |
| unified linear model | `outer_plasma/bem_outer_plasma_unified.f90` |
