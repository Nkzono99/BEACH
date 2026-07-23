title: Advanced rough-surface linear screening

Lang: [日本語](UnifiedLinearResponse.md) | [English](UnifiedLinearResponse.en.md)

# Advanced rough-surface linear screening (`unified_linear_response`)

`external_boundary.field.model="unified_linear_response"` combines the rough-surface height range, plane-averaged surface-charge source,
plasma-accessible area, and linear Debye response in one field model. It is the linear model for cases without a vacuum window in
which lateral modes have decayed between the surface and ownership interface.

This is not the standard outer-sheath model. It is an advanced linear-screening model for cases where roughness and plasma
response occupy the same region. Use `kinetic_1d` for a self-consistent outer sheath that includes species VDFs, the Bohm
condition, or current balance. The unified model owns only the field, not the source VDF; select an inflow correction
independently with `external_boundary.particles.inflow_model`. Do not treat
`unified_linear_response` as a higher-accuracy replacement or automatic fallback for `kinetic_1d`.

## Replace the split window with one field solve

`linear_debye` and `kinetic_1d` join a surface field below an ownership interface to a 1-D outer profile above it. The unified
model solves one zero-mode Poisson grid from the surface projection to the far boundary and continues nonzero modes into plasma
tails.

The interface derived from `sim.box_max[2]` is the particle-ownership boundary on the z-high box face. It is neither a field-solve boundary nor the
start of plasma response. With unchanged mesh geometry and Debye length, moving only the interface must not change the local field
profile.

The plasma response uses the Debye–Hückel relation

$$
\rho_\mathrm{closure}(z)=-\epsilon_0\kappa^2\phi(z),
\qquad
\kappa=\lambda_D^{-1}.
$$

A configuration that carries species VDFs, the Bohm condition, floating-current balance, or a photoelectron mean density is
[Kinetic 1-D outer plasma](KineticOuterPlasma.en.html).

## Derive plasma-accessible area from rough geometry

Let $h(x,y)$ be the highest surface first visible from the plasma at each periodic $(x,y)$. The horizontal area available to
plasma at height $z$ is

$$
f_\mathrm{access}(z)=\frac{1}{A_{xy}}
\int_{A_{xy}}I[z>h(x,y)]\,dx\,dy.
$$

The full-cell mean plasma charge is

$$
\rho_\mathrm{plasma}(z)
=f_\mathrm{access}(z)\rho_\mathrm{closure}(z)
=-\epsilon_0f_\mathrm{access}(z)\kappa^2\phi(z).
$$

$f_\mathrm{access}=0$ below the surface, becomes one above all surface heights, and lies between zero and one across the
roughness range, preventing plasma charge from being placed in horizontal area occupied by solid.

$h(x,y)$ is sampled with `interface_sample_n` cell-centered vertical rays per axis and recomputed at twice that resolution. The
model requires

$$
\max_z|\Delta f_\mathrm{access}(z)|
\le\texttt{accessible\_fraction\_tolerance}.
$$

Overhangs, closed cavities, multiple plasma-facing intersections on one vertical ray, and pores connected to the reservoir only by
lateral paths are outside this single-valued height representation.

## Solve the zero mode from surface charge and plasma response

The zero-mode grid spans

$$
z_{\min}=z_{\mathrm{mesh,min}}-\lambda_D,
\qquad
z_{\max}=z_{\mathrm{mesh,max}}+10\lambda_D
$$

with `unified_grid_points` points, default 129 and minimum 17. Exact triangle-height projection evaluates surface zero-mode field
$E_s(z)$ on the same grid. Its discrete derivative gives surface source

$$
\rho_s(z)=\epsilon_0\frac{dE_s}{dz}
$$

for

$$
\frac{d^2\phi}{dz^2}
-f_\mathrm{access}(z)\kappa^2\phi(z)
=-\frac{\rho_s(z)}{\epsilon_0}.
$$

The lower boundary uses $E_\mathrm{bottom}$ selected by `lower_boundary_model`; the upper boundary is

$$
\phi'(z_{\max})+\frac{\phi(z_{\max})}{\lambda_D}=0.
$$

A tridiagonal discretization supporting nonuniform spacing solves in $O(N_z)$. The same Debye length extrapolates exponentially
beyond the upper grid point, so zero-mode field remains continuous to an ownership plane above the grid.

For a flat surface with $f_\mathrm{access}=1$ and no surface source, the solution returns to
$\phi\propto e^{-z/\lambda_D}$. Regression checks use this analytic limit, grid refinement, and Gauss closure of surface plus
plasma charge.

## Connect vacuum nonzero modes to a screened tail

Vacuum $k\ne0$ modes from the periodic solver are not continued as $e^{-kz}$ to infinity. They join linear plasma response just
above the highest surface at

$$
z_r=z_{\mathrm{mesh,max}}+
\max(\epsilon_\mathrm{roundoff},10^{-6}\lambda_D).
$$

Response start $z_r$ is also independent of the ownership interface. For each Fourier mode,

$$
k=\sqrt{k_x^2+k_y^2},
\qquad
\alpha=\sqrt{k^2+\kappa^2},
$$

with reflection and transmission factors

$$
R_k=\frac{k-\alpha}{k+\alpha},
\qquad
T_k=\frac{2k}{k+\alpha}.
$$

For $z\le z_r$, $R_kI_ke^{k(z-z_r)}$ is added to the base vacuum field. Above $z_r$, base continuation
$I_ke^{-k(z-z_r)}$ is removed and replaced by $T_kI_ke^{-\alpha(z-z_r)}$. Potential, normal field, and tangential field are
continuous at $z_r$. As $\kappa\to0$, $R_k\to0$, $T_k\to1$, and the correction vanishes.

Mode amplitudes integrate `triangle_p0` panels with Duffy quadrature and retain modes through
`periodic2.reference_mode_layers`. The implementation gives no automatic omitted-mode error bound, so increase
`reference_mode_layers` and `panel_quadrature_order` and demonstrate convergence of potential, field, force, or the relevant
observable.

## Check that the response remains in the linear regime

The zero-mode measure is

$$
\eta_0=\frac{\max_z|\phi_0(z)|}{V_T},
$$

and retained nonzero modes use transmitted amplitude in $\eta_k=|q\phi_{k,\mathrm{transmitted}}|/T$. If
$\max(\eta_0,\max_k\eta_k)$ exceeds `max_linearity_ratio`, the model stops because nonlinear lateral coupling is required. This
is an operational small-perturbation gate, not an error bound.

## Geometry, source, and transfer scope

| Item | Contract |
| --- | --- |
| Geometry | x/y periodic, z open, and a single-valued plasma-facing height field |
| Source kernel | `triangle_p0` |
| Prescribed field | `sim.e0=0` |
| Mean plasma | Scalar linear Debye response only |
| Species model | No species VDF, Bohm condition, or photoelectron mean closure |
| Particle transfer | `external_boundary.particles.mode="local_source"`, or `"same_batch"` for a 3-D outer orbit |
| Magnetic field | `particles.mode="same_batch"` requires `sim.b0=0` |
| Failure | Stop without fallback to a nonlinear or legacy sheath model |

`max_gap_ratio` and `max_local_charge_ratio` are diagnostics for split scalar-interface models. Unified acceptance centers on
accessible-area convergence and zero/nonzero linearity.

## Add zero and nonzero responses exactly once to the snapshot

Each field evaluation combines once each:

1. the production FMM or spectral-reference periodic $k\ne0$ base field;
2. the unified Poisson $k=0$ profile;
3. reflection/transmission corrections for retained $k\ne0$ modes.

Every snapshot refresh after surface-charge commit updates the surface zero mode, unified linear solve, and nonzero-tail amplitudes.
The current unified path does not skip solves according to `outer_update_stride`. The MPI root performs the tridiagonal solve and
broadcasts status and $z,\phi,E,\rho$.

With `external_boundary.particles.mode="local_source"`, z-high particles follow the ordinary open boundary. With
`"same_batch"`, they move through the same zero/nonzero outer field. See
[Particle escape and return](ParticleEscapeReturn.en.html).

## Converge geometry, grid, and modes independently

`outer_plasma_profile.csv` stores zero-mode $z,\phi,E,\rho_\mathrm{plasma}$. In `summary.txt`, inspect
`outer_linearity_ratio`, `outer_nonzero_tail_linearity`, `outer_accessible_fraction_min/max`,
`outer_accessible_fraction_refinement_error`, `outer_response_start_z_m`, and the Gauss residual.

For production, vary at least:

- `unified_grid_points` for zero-mode grid refinement;
- `interface_sample_n` and `accessible_fraction_tolerance` for geometry sampling;
- `reference_mode_layers` and `panel_quadrature_order` for the nonzero tail;
- `outer_orbit_dt` when explicit particle transfer is active.

## Code reference

- Accessible fraction and unified Poisson solve: [`bem_outer_plasma_unified.f90`](../src/physics/outer_plasma/bem_outer_plasma_unified.f90)
- 1-D grid operations: [`bem_outer_plasma_grid.f90`](../src/physics/outer_plasma/bem_outer_plasma_grid.f90)
- Zero/nonzero field snapshot composition: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- Explicit outer orbit: [`bem_outer_plasma_orbit.f90`](../src/physics/outer_plasma/bem_outer_plasma_orbit.f90)
- Periodic nonzero operator: [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
