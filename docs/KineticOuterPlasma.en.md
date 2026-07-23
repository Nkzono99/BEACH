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

Exactly one enabled negative and one enabled positive z-high `reservoir_face` species are required for the ambient electron and
ion populations at infinity; BEACH does not silently choose the first of several candidates. With
`photoelectron_density_model="kinetic_mean"`, exactly one enabled negative `photo_raycast` species supplies temperature and
emission current density for a plane-averaged photoelectron source.

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

## Connect Zhao populations to accumulated charge

`kinetic_closure="zhao_charge_driven"` uses Zhao populations of free/reflected ambient electrons, free/captured
photoelectrons, and cold ions, and solves a profile that matches $E_I$ from the currently accumulated surface charge.
Select `zhao_branch="auto"` or `a`, `b`, or `c`. The default
`kinetic_closure="absorbing_maxwellian"` remains the existing finite-grid Poisson model.
The Zhao closure requires exactly one enabled negative z-high `reservoir_face` ambient electron and one positive z-high
`reservoir_face` ion. With `photoelectron_source_scale>0`, it also requires exactly one negative `photo_raycast`
photoelectron. With `photoelectron_source_scale=0`, it rejects an enabled photoelectron source.
The current implementation accepts only `sim.sheath_electron_drift_mode="normal"` and
`sim.sheath_ion_drift_mode="normal"`. It also requires `photo_raycast.normal_drift_speed=0` and the cold-ion condition
$T_i\le0.1T_e$.

Because the charge-driven Zhao closure prescribes the current interface field, it does not impose the legacy Zhao zero-current
equation as a root. Types B and C solve infinity quasineutrality and

$$
2\int_{\psi_0}^{0}\hat\rho(\psi)\,d\psi=\hat E_I^2.
$$

Type A solves infinity quasineutrality, the upper-branch far-field condition, and

$$
-2\int_{\psi_m}^{\psi_0}\hat\rho_{\mathrm{lower}}(\psi)\,d\psi=\hat E_I^2.
$$

Here $\psi=\phi/T_{pe}$ and $\hat E_I=E_I\lambda_{D,pe}/T_{pe}$. The legacy zero-current equation is instead evaluated as
species and total current-density diagnostics, so total current may be nonzero while charge evolves. The transient population
scale $\eta$ scales photoelectron density, infinity quasineutrality, and Sagdeev terms, but it does not scale the raw
photoelectron emission-current term in the current diagnostic; that term retains the full tracked source.
If the initial $E_I=0$ state has no quasineutral root with a strong photoelectron population, an ambient-only flat state is used
once and recorded as branch `0`, representing an outer population that has not formed yet. After the first tracked current
creates charge, ordinary Zhao roots are solved; this bootstrap is never a fallback at nonzero field.

Select the UV-free path with `outer_plasma.photoelectron_source_scale=0`. Zhao photoelectron density and raw emission current
then vanish exactly, and the normalization uses ambient $n_\infty,T_e$. The configured ambient electron density is the total
quasineutral-region density and is checked against the ion density; the solver derives the incoming-electron density needed by
the free/reflected VDF. Reservoir macro-particle counts and velocity cutoffs use that derived density and the same profile map.
The flat $E_I=0$ state is the Type-B/C junction, and a negative $E_I$ selects Type C.

The quasisteady A/B/C branches with the full photoelectron population need not connect continuously to $E_I=0$. A requested
field outside their solvable range stops with `no_physical_solution` under the default `outer_queue_enabled=false` behavior.

### Start stationary studies from the Zhao zero-current root

When the target is a strong-UV stationary observable rather than the turn-on transient from an uncharged state, select
`coupling.steady_start_mode="zhao_floating"`. Before the first batch, this mode uses the configured infinity quasineutral
conditions, temperatures, drifts, and UV source to solve the legacy Zhao zero-current stationary root. It constructs the
`phi(infinity)=0` profile from that root and initializes the uniform plane charge required by its interface field $E_I$.

For horizontal area $A$, the charge follows the zero-mode lower-boundary condition:

$$
Q_{seed}=2\epsilon_0AE_I
\quad\text{for \texttt{symmetric\_vacuum}},
$$

or

$$
Q_{seed}=\epsilon_0AE_I
\quad\text{for \texttt{e\_bottom\_zero}}.
$$

The charge is distributed by triangle area over the horizontal plane selected by `steady_start_mesh_id`. That plane must be
coplanar, cover the full horizontal periodic-cell area, and lie below the outer interface. The current implementation accepts
only `mesh.mode="template"`, whose generated plane guarantees a nonoverlapping, gap-free tiling; arbitrary OBJ planes are
rejected. Other meshes retain zero initial charge. For a plane-plus-sphere model with the plane as mesh 1,
`steady_start_mesh_id=1` seeds only the plane and starts the sphere neutral.

The first outer state, density and velocity map from the infinity reservoir to the interface, and instant return or escape
outside the interface all use the same stationary-root profile. Zero current is imposed only while constructing the initial
state. Subsequent charge-driven refreshes retain the ordinary contract: current surface charge determines $E_I$, and current
remains a diagnostic. Analytic current is not added to surface charge, so it does not double count tracked emission, inflow,
or reabsorption.

This initialization is not a physical transient. It is available only on the instant path with
`outer_queue_enabled=false`, and it does not overwrite existing charge on a fresh run. With the same configuration and
`output.resume=true`, BEACH restores checkpoint mesh charge and a complete outer state without reseeding from the zero-current
root. Combining the warm start with the queue transient closure is rejected. UV turn-on, delayed return current, and transient
cloud inventory still require the queue or a dynamic outer model.
Even for a stationary publication result, a warm start alone does not demonstrate uniqueness or dynamic stability. Confirm
that an independently relaxed or perturbed seed returns to the same stationary observables.

[ADR 0006](adr/0006-zhao-stationary-warm-start.md) records this scope and its separation from the queue transient closure.

With `outer_queue_enabled=true`, tracked photoelectrons retain their macro-particle weights as queue inventory while they occupy
the outer region. The MPI-global photoelectron number divided by horizontal area is the target column. The closure solves a
population scale $\eta$ so the Zhao integral of $n_{pe,f}+n_{pe,c}$ over $0\le z\le10\lambda_{D,pe}$ matches that target.
Despite the output name `outer_photoelectron_population_fraction`, $\eta$ is an occupancy scale relative to the stationary
reference population, not a probability. The solver follows the physical path connected to $\eta=0$ over
$0\le\eta\le16$ and permits transient overshoot above one. It does not clamp to `[0,1]`, fall back to a target-independent
full-population solution or branch `0`, or jump to a disconnected branch. Queue mode requires `zhao_branch="auto"` and allows
only continuous A/B or other branch transitions satisfying the degeneracy condition. The current bisection accepts only paths
whose column increases monotonically with $\eta$; a continuous path containing a fold is unsupported and stops. If no connected,
monotone path reaches the target, the solve stops with `no_physical_solution`.

### Diagnose a Zhao continuation failure

When Zhao continuation fails closed, the MPI root writes five stderr records prefixed with
`BEACH zhao-continuation` before the generic error. `call_stage` is `pre_batch` or `post_enqueue`, and `batch` identifies the
failing batch. The remaining records contain the solver stage, reason code, underlying and return status, target field,
`attempt`, `attempted_step`, $\eta$, and column, the previous root, the rejected candidate, root residuals,
`normalized_potential_jump`, `log_density_jump`, and `normalized_root_jump`.
An unavailable branch is encoded as `from_branch=-` or `candidate_branch=-`. Reals use a full-range scientific format with a
three-digit exponent. The root flushes all five records before every MPI rank stops.

The public Fortran procedure `trace_zhao_branch_atlas` is a diagnostic pseudo-arclength tracer for one requested A, B, or C
branch. It is not connected to production continuation, root selection, or fallback. Reaching a finite density floor or point
limit is a `search_limit`; without tracing other branches and degeneracy connections, it does not establish that the target is
unreachable on every stationary Zhao root.

The public Fortran procedure `diagnose_zhao_ab_degeneracy` examines a Type-B density-zero endpoint in the coordinate
$q=\sqrt{-\phi_m/T_{pe}}$. It records the $q^3$ coefficient of the Type-A far-field residual along the quasineutral curve, the
Type-A/B interface-field integral difference, the supplied Type-B quasineutral and field residuals, and a finite-$q$ probe in
`zhao_ab_degeneracy_diagnostics_type`. `regular_connection_conditions_met` represents necessary connection conditions only; it
does not establish the absence of a disconnected Type-A component or its dynamic stability. This API also leaves production
branch selection unchanged.

The public Fortran procedure `trace_zhao_field_column_homotopy` linearly interpolates a previous state $(E_0,N_0)$ and a target
$(E_1,N_1)$, then pseudo-arclength traces the Zhao residuals together with the finite-length column residual on one fixed Type B
or C branch with a nonzero photoelectron source. When accepted points bracket $\lambda=1$, a fixed-$\lambda=1$ corrector lands exactly on the target.
`target_reached` is true only after that corrector converges, and `homotopy_fold_detected` records a sign reversal of the tangent's
$\lambda$ component. The API does not switch branches or alter production root selection or fallback.
The nonmonotone five-coordinate Type-A corrector is explicitly unsupported until it is adequately conditioned, and no-photo
Type C is rejected because its column equation is identically zero. Type A remains available through the fixed-field atlas and
A/B endpoint diagnostics.

For the strong-UV fixture, neither the forward Type-B atlas to its density floor nor the reverse atlas to its $\eta$ lower bound
brackets the target column at $E_I=0.9072962759\,\mathrm{V/m}$ and $L=10\lambda_{D,pe}$. At the density-zero limit, the
quasineutral $q^3$ coefficient is approximately $2.9\times10^{-2}$ rather than zero, so a necessary condition for a regular
local Type-A tangent continuous with that limit is not met. The target is therefore unreachable on the refined Type-B component containing the runtime neighborhood.
This result does not exclude a different $L$ or a disconnected Type-A root, and the production solver continues to fail closed.
[ADR 0005](adr/0005-zhao-continuation-and-dynamic-outer.md) records the numerical evidence, scope, and decision to proceed toward
a dynamic outer model.

The same straight homotopy was traced from a Type-B batch-15 state recovered by a successful 15-batch replay to the failed
batch-16 target. The forward curve starts at $E_0=0.8424570666\,\mathrm{V/m}$ and
$N_0=9.3202065681\times10^7\,\mathrm{m^{-2}}$, then reaches the finite ambient-density floor at
$\lambda\simeq0.33179$ and $E_I\simeq0.86397\,\mathrm{V/m}$. Changing
$\log(n_{e,\infty}/n_{pe,ref})$ from $-27$ to $-24$ moves the endpoint by less than $10^{-5}$ in $\lambda$, consistent with
convergence toward the density-zero limit. At every accepted point, the root and normalized-column residuals are at most
$5\times10^{-10}$, and the row-rank indicator is finite and above the numerical rejection threshold of $10^{-12}$. It does not reach the target
$E_1=0.9072962759\,\mathrm{V/m}$ and $N_1=9.9455765203\times10^7\,\mathrm{m^{-2}}$, and no $\lambda$ fold is detected.
The quasineutral $q^3$ coefficient at this endpoint is approximately $2.9\times10^{-2}$, so a necessary condition for a regular
local Type-A tangent is not met. The fixed-branch quasisteady path in the increasing-$\lambda$ direction from the batch-15 root therefore ends at
the finite floor before the target. This diagnostic does not cover the global curve in the reverse direction, disconnected components, or a different
$L$.

The same $0\le z\le10\lambda_{D,pe}$ interval is the finite queue-particle control volume. A particle that does not turn inside
it is absorbed by the exterior reservoir and escapes at $L=10\lambda_{D,pe}$; queue mode does not use a Robin tail outside $L$
to classify return. For each event, it applies the frozen-field bound to
`tau_outer + delta_poll + batch_duration/2`, where `delta_poll` is the quantization delay from `t_due` to the first batch-start
poll and the half-batch term bounds uncertainty from approximating the crossing time by the batch midpoint.
Configuration validation applies the same bound to `batch_duration`. See
[Particle escape and return](ParticleEscapeReturn.en.html#queue-outer-flight-for-the-transient-zhao-closure) for time
discretization, the end-of-batch corrector, and checkpoint rules.

[`periodic2_zhao_charge_driven_outer.toml`](../examples/periodic2_zhao_charge_driven_outer.toml) is a warm-start smoke that
does not solve the transient. It uses `zhao_floating` to begin ordinary small batches from the stationary Type-A root.
[`periodic2_zhao_no_photo_outer.toml`](../examples/periodic2_zhao_no_photo_outer.toml) is a no-photo smoke that charges from
the flat state toward Type C using only the same ambient inputs.
[`periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml) is a frozen-field guard fixture for the
strong-UV queue closure. With its stated physical timescale, a long outer flight is expected to stop fail-closed; it is not a
successful validation run. Interpreting return current requires a separately justified slower field timescale and convergence
checks in `batch_duration`, particle count, and control-volume length.

This first implementation treats the z-high interface as the effective Zhao emitting plane. With
`photoelectron_source_scale>0`, it reuses `sim.sheath_alpha_deg` and `sim.sheath_photoelectron_ref_density_cm3`, and obtains
mass and temperature from the first negative `photo_raycast` species. The Zhao solver uses that temperature $T_{pe}$ as its
potential scale and the photoelectron Debye length
$\lambda_{D,pe}$ derived from the temperature and reference density as its length scale. It writes the derived length as
`outer_debye_length_m` in the converged state.

`outer_plasma.debye_length` and `outer_plasma.thermal_voltage` are not physical scales for the Zhao root or profile. They remain
reference inputs for the split-interface applicability diagnostics: the former scales `interface_eta_gap` and the local-charge
estimate, while the latter scales the lateral `interface_eta_phi_kneq0` and `interface_eta_field_kneq0` diagnostics. Do not test
Zhao convergence by varying these inputs. Refine the profile grid, effective interface location, tracked-particle count, and
`dt` or batch resolution instead. The current 128-point profile grid is fixed at runtime, so a production grid-convergence study
requires exposing its point count as an input.

The tracked `photo_raycast.emit_current_density_a_m2` must agree within 1% with the analytic raw source at the effective plane,

$$
J_{pe,\mathrm{raw}}=s_{UV}
\frac{|q_{pe}|n_{\mathrm{ref}}\sin(\alpha)v_{\mathrm{th},pe}}{2\sqrt{\pi}},
\qquad
v_{\mathrm{th},pe}=\sqrt{\frac{2T_{pe}}{m_{pe}}}.
$$

Here $s_{UV}$ is `photoelectron_source_scale` and $T_{pe}$ in the speed formula is thermal energy in joules. The runtime rejects
a mismatch. The analytic raw current enters
the tracked-source consistency check and current-density diagnostics, but not the root, surface charge, or ledger; only tracked
emission and reabsorption update the latter two. The population scale $\eta$ does not scale that raw photoelectron
emission-current term in the current diagnostic.
The closure also rejects a legacy Zhao `sheath_injection_model`,
`reservoir_potential_model`, and `photoelectron_density_model="kinetic_mean"`.

This effective-plane approximation does not self-consistently connect tracked-ray directions or a VDF reaching the interface
from a rough surface to the Zhao outer population. `ray_direction` controls illumination-ray sampling of emitting surfaces, while $\alpha$
independently controls the analytic Zhao source. See [ADR 0003](adr/0003-zhao-charge-driven-outer-closure.md) for scope and the
boundary VDF needed by a future generalization, [ADR 0004](adr/0004-zhao-transient-photoelectron-column-queue.md) for the
transient-queue decision, [ADR 0005](adr/0005-zhao-continuation-and-dynamic-outer.md) for staged continuation diagnostics
and a possible dynamic outer model, and [ADR 0006](adr/0006-zhao-stationary-warm-start.md) for the stationary warm start.

## Connect the `absorbing_maxwellian` finite grid to infinity with a Robin tail

This section and the next Newton section are specific to the default `absorbing_maxwellian`. Zhao constructs its root and profile
from the Sagdeev conditions described above. `absorbing_maxwellian` interior points use conservative finite-volume residuals on
a stretchable nonuniform 1-D grid. Current runtime values are below;
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

## Follow the physical branch with continued Newton solves (`absorbing_maxwellian`)

Analytic density derivatives form a bordered-tridiagonal Jacobian, making one Newton step $O(N_z)$. Dependence of interior
densities on $\phi_I=\phi_1$ produces the border column in addition to the ordinary tridiagonal stencil.

1. Form a Newton step.
2. Use backtracking to keep the trial state on the supported monotone branch.
3. Apply pseudo-transient regularization if ordinary Newton stalls.
4. With a previous profile, continue from its interface field to the current field and halve a failed increment.

Regularization and continuation change only the convergence path. Final validation always returns to the original discrete
Poisson residual.

## Accept only branches supported by each closure

| Condition | Requirement |
| --- | --- |
| Original Poisson residual | Unregularized residual is below tolerance |
| Branch | The supported monotone branch for `absorbing_maxwellian`; the requested A/B/C sign and population conditions for Zhao |
| Ion accessibility | $u_i^2(z)>0$ at every grid point |
| Kinetic Bohm entry | $u_{i,\infty}\ge\sqrt{(T_e+\gamma_iT_i)/m_i}$ |
| Infinity quasineutrality | $q_en_{e,\infty}+q_in_{i,\infty}\simeq0$ |

Nonmonotone virtual-cathode profiles and trapped populations remain outside `absorbing_maxwellian`.
`zhao_charge_driven` supports the single potential minimum of Type A and only the reflected/captured populations represented by
the Zhao formulas. Sub-Bohm ion inflow and states outside the selected closure are rejected. The run stops with status
`not_applicable`, `no_physical_solution`, or `numerical_failure`.

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

For `absorbing_maxwellian`, vary `debye_length`, interface location, and when needed source sampling, and confirm convergence of
interface potential, currents, and surface charging. Zhao uses the derived $\lambda_{D,pe}$, so refine the profile grid,
effective interface location, tracked-particle count, and time resolution instead. With return active, inspect flight time,
frozen-field ratio, and quasisteady applicability in [Particle escape and return](ParticleEscapeReturn.en.html).

## Code reference

- VDF closures and nonlinear Poisson solve: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- Charge-driven Zhao roots and nonmonotone profiles: [`bem_outer_plasma_zhao.f90`](../src/physics/outer_plasma/bem_outer_plasma_zhao.f90)
- Build solver options from runtime species: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)
- Stationary-root and plane-charge warm start: [`bem_zhao_steady_start.f90`](../src/runtime/coupling/bem_zhao_steady_start.f90)
- Surface-field coupling and MPI collective solve: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- Profile output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- Profile restart: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
