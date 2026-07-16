title: Kinetic 1-D outer plasma

Lang: [日本語](KineticOuterPlasma.md) | [English](KineticOuterPlasma.en.md)

# Kinetic 1-D outer plasma

`outer_plasma.model="kinetic_1d"` solves the plasma above a periodic surface as a plane-averaged, one-dimensional,
electrostatic, collisionless, and unmagnetized profile. It connects the interface field produced by surface charge to velocity
distributions at infinity and obtains interface potential, density, and currents self-consistently.

BEACH recommends `kinetic_1d` as the standard model that connects an external reservoir to a mean sheath. With the matching
`return_model` and `particle_transfer_mode`, the same profile also controls particle inflow and escape/return. Start outer-sheath
cases with this model unless the case has a specific requirement to resolve linear lateral-field screening near a rough surface.

The converged outer profile is used in both directions. Mapping particles from the infinity VDF to the interface is covered by
[`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html); mapping particles leaving the interface to escape or return is covered by
[Particle escape and return](ParticleEscapeReturn.en.html).

## Split local and outer regions at the interface

The local particle region meets the outer region at the z-high ownership interface $z=z_I$.

| Region | Field used |
| --- | --- |
| Mesh to interface | Periodic $k\ne0$ surface field and surface $k=0$ field |
| Beyond the interface | Plane-averaged `kinetic_1d` profile |

The outer profile is not superposed at every point in the local region. The split contract assumes lateral modes have decayed
sufficiently by the interface and matches potential and normal field there. Only when roughness and plasma response must be
solved linearly in the same region, use [advanced rough-surface linear screening](UnifiedLinearResponse.en.html).

## Determine potential from interface field and infinity VDFs

The unknown is plane-averaged potential $\phi(z)$ from the interface to the infinity reservoir. Interface field $E_I$ from the
surface zero mode supplies the Neumann condition

$$
-\phi'(z_I)=E_I,
$$

and particle VDFs construct $\rho(\phi)$ for

$$
\frac{d^2\phi}{dz^2}=-\frac{\rho(\phi)}{\epsilon_0}.
$$

Infinity potential is fixed as gauge $\phi_\infty=0$. Interface potential $\phi_I=\phi(z_I)$ is not an input; it is determined by
the interface field, VDF closures, and far boundary together.

The root equations are Poisson's equation and its boundary conditions. Surface charge changes $E_I$ between batches, which
changes $\phi_I$ and species currents. Electron, ion, photoelectron, and external-circuit current densities are diagnostics
evaluated from the converged profile.

## Map VDFs to potential-dependent charge density

The first negative and positive z-high `reservoir_face` species define ambient electrons and ions at infinity. With
`photoelectron_density_model="kinetic_mean"`, the first negative `photo_raycast` species supplies temperature and emission current
density for a plane-averaged photoelectron source.

| Population | Inputs | Outer density construction |
| --- | --- | --- |
| Ambient electron | $n_{e,\infty},T_e,q_e,m_e$ | Map a half-Maxwellian by total-energy conservation, including absorbed and potential-reflected trajectories |
| Ion | $n_{i,\infty},T_i,q_i,m_i,u_{i,\infty}$ | Map a cold beam by energy and flux conservation |
| Photoelectron | $T_{pe},q_{pe},m_{pe},\Gamma_{pe,0}$ | Mean outgoing and post-turning return populations from a surface half-Maxwellian |

The cold-ion closure is

$$
u_i(z)=\sqrt{u_{i,\infty}^2-\frac{2q_i\phi(z)}{m_i}},
\qquad
n_i(z)=n_{i,\infty}\frac{u_{i,\infty}}{u_i(z)}.
$$

A profile for which the square root is not real is rejected because ions cannot access the interface.

Photoelectron escape fraction at infinity is

$$
f_{pe,\mathrm{esc}}
=\exp\left[-\frac{\max\{0,q_{pe}(\phi_\infty-\phi_I)\}}{T_{pe}}\right],
$$

with the remainder forming the return population. Temperature is represented internally in joules. The Poisson source is

$$
\rho(\phi)=q_en_e(\phi)+q_in_i(\phi)+q_{pe}n_{pe}(\phi).
$$

Each analytic closure also returns $\partial n_s/\partial\phi$ and, where needed,
$\partial n_s/\partial\phi_I$ for the Newton Jacobian.

Outgoing and returning densities in `kinetic_mean` are a stationary closure for outer space charge. They do not replace surface
deposition by tracked particles and do not add a second statistical return current. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for source charge and tracked reabsorption.

## Connect the finite grid to infinity with a Robin tail

Interior points use conservative finite-volume residuals on a stretchable nonuniform 1-D grid. Current runtime values are below;
only `debye_length` is currently exposed as a separate input.

| Item | Current value |
| --- | ---: |
| Grid points | 128 |
| Domain length | $10\lambda_D$ |
| Grid stretch | 2 |
| Maximum Newton iterations | 40 |
| Residual tolerance | $10^{-8}$ |

Beyond finite-grid endpoint $L$, an exponential tail with $\lambda_\mathrm{tail}=\lambda_D$ gives the Robin condition

$$
\phi'(L)+\frac{\phi(L)}{\lambda_\mathrm{tail}}=0.
$$

This gives exponential relaxation toward the infinity gauge. The remaining tail is also used for return-particle flight-time
integration.

## Follow the physical branch with continued Newton solves

Analytic density derivatives form a bordered-tridiagonal Jacobian, making one Newton step $O(N_z)$. Dependence of interior
densities on $\phi_I=\phi_1$ produces the border column in addition to the ordinary tridiagonal stencil.

1. Form a Newton step.
2. Use backtracking to keep the trial state on the supported monotone branch.
3. Apply pseudo-transient regularization if ordinary Newton stalls.
4. With a previous profile, continue from its interface field to the current field and halve a failed increment.

Regularization and continuation change only the convergence path. Final validation always returns to the original discrete
Poisson residual.

## Accept only the supported sheath branch

| Condition | Requirement |
| --- | --- |
| Original Poisson residual | Unregularized residual is below tolerance |
| Monotone branch | Remain on the supported electron-repelling profile |
| Ion accessibility | $u_i^2(z)>0$ at every grid point |
| Kinetic Bohm entry | $u_{i,\infty}\ge\sqrt{(T_e+\gamma_iT_i)/m_i}$ |
| Infinity quasineutrality | $q_en_{e,\infty}+q_in_{i,\infty}\simeq0$ |

Nonmonotone virtual-cathode profiles, trapped populations, and sub-Bohm ion inflow are outside this model. A numerical iterate is
not accepted if these conditions fail. The run stops with status `not_applicable`, `no_physical_solution`, or
`numerical_failure`.

## Refresh the profile at its stride and share it across ranks

On batches selected by `outer_update_stride`, the interface field is rebuilt from committed surface charge and the profile is
updated. Only a previous profile with the same model identity and grid can seed Newton. A skipped batch keeps the previous outer
state, while the surface-side field snapshot still refreshes from current committed charge.

The MPI root performs the 1-D solve and broadcasts status, profile, and current diagnostics. All particles in a batch share the
updated immutable snapshot; the outer profile is not solved again after each impact.

## Converge residuals, currents, and charging together

Converged $z,\phi,E,\rho$ are written to `outer_plasma_profile.csv` and can seed restart. At minimum inspect
`interface_potential`, `interface_field`, `outer_integrated_charge`, species and total current densities, Newton iteration count,
and original nonlinear residual.

Also vary `debye_length`, interface location, and when needed source sampling, and confirm convergence of interface potential,
currents, and surface charging. With return active, inspect flight time, frozen-field ratio, and quasisteady applicability in
[Particle escape and return](ParticleEscapeReturn.en.html).

## Code reference

- VDF closures and nonlinear Poisson solve: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- Build solver options from runtime species: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)
- Surface-field coupling and MPI collective solve: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- Profile output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- Profile restart: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
