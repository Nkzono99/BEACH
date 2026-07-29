title: Kinetic 1-D outer plasma

Lang: [日本語](KineticOuterPlasma.md) | [English](KineticOuterPlasma.en.md)

# Kinetic 1-D outer plasma

`external_boundary.field.model="kinetic_1d"` solves the plasma above a periodic surface as a plane-averaged, one-dimensional,
electrostatic, collisionless, and unmagnetized profile. It connects the interface field produced by surface charge to velocity
distributions at infinity and obtains interface potential, density, and currents self-consistently.

BEACH recommends `kinetic_1d` as the standard model that connects an external reservoir to a mean sheath. With
`external_boundary.particles.mode="same_batch"`, the same profile also controls particle inflow and ordinary escape/return.
For the `implicit_mean` photoelectron paths described below, the initial local trace stops at the interface crossing and holds
that crossing until the updated mean field is available. The two sibling paths use an analytic-Maxwellian backward-Euler mean
update for `ambient_linear_debye` and a generalized empirical-interface-energy root for `zhao_charge_driven`. Both use the
kinetic 1-D profile mapper after the mean update, but they do not share the rule that determines mean escape/return charge or
ray weights.
Start outer-sheath cases with this model unless the case has a specific requirement to resolve linear lateral-field screening
near a rough surface.
The current model set does not support that requirement. [ADR 0010](adr/0010-remove-unified-linear-response.md) records the
removed approximation and the criteria for a redesign.

The converged outer profile is used in both directions for ordinary particle transfer. Mapping particles from the infinity VDF to the interface is covered by
[`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html); mapping particles leaving the interface to escape or return is covered by
[Particle escape and return](ParticleEscapeReturn.en.html). This page explains how `implicit_mean` places that common map
in a one-pass orbit sample that separates fast mean charging from slow local redistribution.

## Split local and outer regions at the interface

The local particle region meets the outer region at the z-high ownership interface $z=z_I$.

| Region | Field used |
| --- | --- |
| Mesh to interface | Periodic $k\ne0$ surface field and surface $k=0$ field |
| Beyond the interface | Plane-averaged `kinetic_1d` profile |

The outer profile is not superposed at every point in the local region. The split contract assumes lateral modes have decayed
sufficiently by the interface and matches potential and normal field there.

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

## Separate the ambient linear Debye response from tracked photoelectrons

`kinetic_closure="ambient_linear_debye"` analytically constructs

$$
\phi_I=\lambda_D E_I,\qquad
\phi(z)=\phi_I\exp(-z/\lambda_D),\qquad
\rho_{\mathrm{amb}}(z)=-\frac{\epsilon_0}{\lambda_D^2}\phi(z)
$$

from the interface field $E_I$ supplied by the current surface zero mode. It performs no spatial Poisson root or Newton
iteration. The infinity gauge is $\phi(\infty)=0$. The same profile controls ambient-reservoir accessibility and velocity
mapping, ordinary z-high return or escape, and `implicit_mean` photoelectron classification after the mean solve.
Diagnostic ambient inflow currents use the one-sided
drifting-Maxwellian flux of each `reservoir_face` species with the same interface barrier.

Do not write `photoelectron_density_model` in public TOML for this closure. The facade resolves its internal state to `none`
and rejects the explicit key even when its value is `"none"`. It does not reject enabled `photo_raycast` sources:
photoelectrons are emitted and tracked within the 3-D region until reabsorption or an outward z-high interface crossing.
Their mean density and outer space charge are excluded from the closure.

When `external_boundary.particles.mode="same_batch"` and an enabled negative `photo_raycast` species are also present,
BEACH automatically selects internal `coupling.update_mode="implicit_mean"`. Do not write an update mode in public TOML.
This multirate update separates fast mean charging, resolved within the current batch by a continuous Maxwellian closure,
from slow local charge redistribution sampled by finitely many tracked rays and carried into the next batch.

| Component | Update |
| --- | --- |
| Local $k\ne0$ | Hold the batch-start nonzero-mode operator fixed and share it between the initial 3-D trace and post-return local traces |
| Fast mean mode | Backward-Euler mean solver determines total charge, $\phi_I$, and continuous-Maxwellian escape/return fractions |
| Energy-resolved rays | Trace each crossing once on the solved mean profile to sample its orbit and terminal destination |
| Slow local redistribution | Normalize analytic return charge over the source/destination distribution of reabsorbed samples and commit it once |

Analytic escape/return replacement applies only to the component deferred after an outward z-high crossing. Local
reabsorption before the interface, escape through another open face such as z-low, and unresolved discard retain their
tracked values from the first trace. The mean solve and return normalization below do not overwrite them.

A photoelectron's initial local trace ends when it crosses z-high. Its full macro weight, source element, crossing position,
and velocity pass to the mean update. BEACH measures $J_{pe,\mathrm{out}}^n$ from crossing charge in that batch.
It also measures $J_{\mathrm{pending}}^n$ from every surface-charge delta staged by the first trace and defines

$$
J_{\mathrm{other}}^n
=J_{\mathrm{pending}}^n-J_{e,\mathrm{tracked}}^n-J_{i,\mathrm{tracked}}^n-J_{pe,\mathrm{out}}^n.
$$

The mean interface potential then advances according to

$$
C_A\frac{\phi_I^{n+1}-\phi_I^n}{\Delta t}
=J_{e,\mathrm{tracked}}^n+J_{i,\mathrm{tracked}}^n
+J_{\mathrm{other}}^n
+J_{pe,\mathrm{out}}^n f_{\mathrm{esc}}(\phi_I^{n+1})+J_{\mathrm{ext}}.
$$

Here $C_A=\epsilon_0/\lambda_D$ for `e_bottom_zero` and
$C_A=2\epsilon_0/\lambda_D$ for `symmetric_vacuum`. Ambient current is measured from the particles actually absorbed by the
surface during that batch, restricted to the unique z-high `reservoir_face` electron and ion species. Measured surface
changes from extra species, other-face inflow, another open-face escape, or unresolved particles remain in
$J_{\mathrm{other}}^n$. `emit_current_density_a_m2` determines tracked 3-D emission weights, but its
value is not reused as an independent PE mean source. The mean-source amplitude is measured
$J_{pe,\mathrm{out}}^n$. Configured temperature supplies Maxwellian barrier transmission $f_{\mathrm{esc}}$, evaluated at
$\phi_I^{n+1}$ so the mean retarding potential acts within the same step.
The backward-Euler solution is authoritative for batch-end mean total charge and $\phi_I$. In terms of nonnegative
charge magnitudes, its analytic return and escape charges are

$$
Q_{\mathrm{ret}}^{\mathrm{ana}}
=A\Delta t\,J_{pe,\mathrm{out}}^n(1-f_{\mathrm{esc}}),
\qquad
Q_{\mathrm{esc}}^{\mathrm{ana}}
=A\Delta t\,J_{pe,\mathrm{out}}^n f_{\mathrm{esc}}.
$$

Finite-ray classification does not replace these totals.
For the element charge used during ray tracing, BEACH multiplies each crossing's source countercharge by the analytic return
fraction and temporarily neutralizes that amount at the source. This temporary distribution maps the mean-solved total into
the $k=0$ snapshot; it is not the physical reabsorption destination.

Each crossing keeps its full macro weight. The mapper evaluates its normal energy,

$$
v_n^2(z)
=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m},
$$

along the discrete profile and far tail. A zero of $v_n^2$ gives a turning return; a particle that traverses the profile
escapes. For a return, the same mapper integrates outer flight time, advances the lateral position by
$\mathbf v_t\tau_{\mathrm{outer}}$ with periodic wrapping, and places the packet just inside the interface. BEACH retraces it
through the local 3-D region with batch-start $k\ne0$ and mean-solved $k=0$. If the returned packet crosses z-high again, the
same driver and profile mapper handle the recrossing until it is absorbed by a surface or escapes to infinity.
A shadow that leaves through another local open face after return is not reclassified as analytic upward escape; it fails closed.

This deferred packet is a quasistationary shadow that samples the orbit and destination of the analytic Maxwellian return
total. The mapper still computes each shadow's outer flight time and frozen-field ratio and retains them as diagnostics, but
an over-limit ratio does not stop the run. Ordinary `same_batch` particles and ambient species are not shadows and continue
to stop fail-closed when they exceed the same limit.

For this shadow, `summary.txt` reduces only the analytically weighted return excursions from the last batch completed by the
current invocation to two values. A no-op resume that advances no batch omits these two keys because they are run-local
derived diagnostics, not solver state that must be restored. Let $W_j>0$ be the positive charge magnitude of return
excursion $j$ after analytic weighting, let $\tau_j$ be its
outer round-trip time, let $A$ be the horizontal area, and let $\Delta t$ be the batch duration. Then

$$
\bar{\tau}_{\mathrm{ret},Q}
=\frac{\sum_j W_j\tau_j}{\sum_j W_j},
\qquad
\widehat{\sigma}_{\mathrm{PE,ret}}
=\frac{\sum_j W_j\tau_j}{A\Delta t}
=J_{\mathrm{ret,gross}}^{(|Q|)}\bar{\tau}_{\mathrm{ret},Q}.
$$

`implicit_mean_last_returned_outer_flight_time_mean_s` is the charge-magnitude-weighted
$\bar{\tau}_{\mathrm{ret},Q}$, while
`implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2` is the Little's-law quasistationary
shadow estimate $\widehat{\sigma}_{\mathrm{PE,ret}}$. Both are nonnegative and are zero when there is no return excursion.
The latter is not actual cloud stock in a queue or ledger. No residence time is assigned to the `escaped_to_infinity`
outcome itself, but completed return excursions before that ray eventually escapes are included. Neither value is a
cumulative or maximum value over all batches.

When `outer_integrated_charge_per_area_C_m2` is nonzero,

$$
\chi_{\mathrm{PE,shadow}}
=\frac{\widehat{\sigma}_{\mathrm{PE,ret}}}
       {\left|\texttt{outer\_integrated\_charge\_per\_area\_C\_m2}\right|}
$$

compares the magnitude of the omitted returning-photoelectron shadow column with the magnitude of integrated charge in the
1-D outer profile. It is an interpretation aid for the model limitation, not an acceptance threshold enforced by BEACH.

Orbit tracing runs once after the mean solve; its classification is not fed back into that solve. Let $\mathcal R$ be the
samples ultimately reabsorbed by a surface and define

$$
W_{\mathcal R}=\sum_{j\in\mathcal R}|q_j|w_j,
\qquad
s_{\mathrm{ret}}=\frac{Q_{\mathrm{ret}}^{\mathrm{ana}}}{W_{\mathcal R}}.
$$

The same $s_{\mathrm{ret}}$ multiplies each sample's source leg—the temporary neutralization of its emission
countercharge—and destination leg at the actual hit. The transaction therefore remains zero-sum while its total return charge
matches the analytic value. The final state retains the emission countercharge already present in pending charge and adds the
normalized destination deposit. The sampled escape fraction is an orbit-statistics diagnostic, not an input to mean total charge.

In this one-pass composition, mean total charge and retarding potential respond within the current batch, while the sampled
reabsorption pattern changes local $k\ne0$ only after commit in the next batch. Not feeding discrete classification back into
the same mean solve prevents a few rays near the separatrix from producing a return/escape two-cycle.

This path requires
`outer_update_stride=1`, positive `batch_duration`, exactly one negative `photo_raycast` species,
`deposit_opposite_charge_on_emit=true`, and `photo_raycast.normal_drift_speed=0`. The analytic Maxwellian escape fraction
does not include normal emission drift, so a nonzero value is rejected fail-closed. Omit `photoelectron_density_model` from
public TOML; it resolves internally to
`none`. No public update-mode or return-kernel key is added, and the path cannot be combined with an outer queue. This path
adds no mesh-mode-specific requirement. The run fails closed if the mean solver does not converge, analytic return charge is
positive but no reabsorbed sample is available, transaction charge does not balance, or a deferred packet does not reach
absorption or escape within the allowed local trace.

The applicability range is weak photoemission for which the ambient linear response dominates. The model cannot represent a
photoelectron-space-charge virtual cathode, a space-charge-limited or inverse sheath, trapped populations, or a nonmonotone
profile. It also does not evolve a nonlinear photoelectron mean density or outer-cloud inventory. In particular, it does not
solve UV turn-on from an uncharged state or delayed return-current transients. That use requires a separately designed
delayed inventory or queue compatible with the time-dependent `ambient_linear_debye` closure.
Check `coupling_update_mode=implicit_mean` in `summary.txt`. The `BEACH implicit-mean` progress record reports $\phi_I$,
species currents, `J_other_A_m2`, `transaction_residual_C`, `mean_solver_iterations`, `sample_escape_fraction`, and
`return_weight_scale`. In the photoelectron charge ledger, `escaped_to_infinity_C` combines tracked other-open-face escape
with analytic z-high deferred escape, while `absorbed_on_surface_C` combines tracked local reabsorption with analytic
post-return absorption. `discarded_unresolved_C` remains tracked. The corresponding integer counts and
`sample_escape_fraction` describe terminal ray classifications.
`interface_outward_gross_C` includes initial crossings and later recrossings, while `interface_returned_gross_C` records
returned-event charge carried back into the local region with normalized weights. If
$Q_{\mathrm{esc,z-high}}^{\mathrm{analytic}}$ is the signed analytic escape charge of the z-high deferred component, then
`interface_returned_gross_C` = `interface_outward_gross_C` -
$Q_{\mathrm{esc,z-high}}^{\mathrm{analytic}}$. Do not substitute total `escaped_to_infinity_C` into this identity because it
also includes tracked escape through other open faces.

A UV-only case with no ambient electron/ion reservoirs and
`external_boundary.field.model="none"` that lets particles escape at z-high is not a substitute for this closure. Treat it as
a finite-box transient control for local emission and reabsorption, not as an infinity-closed quasineutral sheath or a
stationary solution.

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

Because the charge-driven Zhao closure prescribes the current interface field, it does not impose the stationary Zhao zero-current
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

Select the UV-free path with `external_boundary.field.photoelectron_source_scale=0`. Zhao photoelectron density and raw emission current
then vanish exactly, and the normalization uses ambient $n_\infty,T_e$. The configured ambient electron density is the total
quasineutral-region density and is checked against the ion density; the solver derives the incoming-electron density needed by
the free/reflected VDF. Reservoir macro-particle counts and velocity cutoffs use that derived density and the same profile map.
The flat $E_I=0$ state is the Type-B/C junction, and a negative $E_I$ selects Type C.

The quasisteady A/B/C branches with the full photoelectron population need not connect continuously to $E_I=0$. A requested
field outside their solvable range stops with `no_physical_solution` in
`external_boundary.particles.mode="same_batch"`.

### Close strong photoemission with the measured interface-energy CDF

Combining photoemitting `zhao_charge_driven` with `same_batch` makes BEACH select internal `implicit_mean` automatically.
This path uses a different solver from the ambient backward-Euler mean update above. It requires a committed A/B/C branch
from the start, so `external_boundary.particles.steady_start_mode="zhao_floating"` and the corresponding
`steady_start_mesh_id` that receives the uniform seed charge are mandatory.

The initial 3-D trace freezes the batch-start field and records the normal energy and macro-charge magnitude when a
photoelectron first crosses z-high outward:

$$
{\cal E}_{n,r}=\frac{1}{2}m_{pe}v_{z,r}^2,\qquad w_r=|\Delta q_r|.
$$

Samples $({\cal E}_{n,r},w_r)$ from every rank are gathered to the MPI root. The root stably sorts energies in descending
order and combines charge weights only for energies that compare exactly equal in floating-point arithmetic. It uses no
configured bins, histogram interpolation, or smoothing, so the result is an exact charge-weighted empirical CDF of the
macro-particle samples.

Let $Q_{\rm base}$ be the total staged surface charge from the first trace with all interface photoelectron source charge
removed. A candidate cell charge $Q$ maps to interface field according to the lower boundary:

$$
Q=\beta\epsilon_0AE_I,\qquad
\beta=
\begin{cases}
1,&\texttt{e\_bottom\_zero},\\
2,&\texttt{symmetric\_vacuum}.
\end{cases}
$$

BEACH solves the 1-D profile on the same Zhao branch for each candidate. The normal energy needed for photoelectron escape is
not just the endpoint difference between the interface and infinity. It is the full-profile barrier

$$
B(Q)=\max\left[
0,\ \max_z q_{pe}\{\phi(z;Q)-\phi_I(Q)\},\
q_{pe}\{\phi_\infty-\phi_I(Q)\}
\right].
$$

Because $q_{pe}<0$, the potential minimum of a Type-A profile contributes the virtual-cathode barrier to this maximum.

The empirical CDF defines $C_{\rm esc}(B)$ by assigning groups with ${\cal E}_{n,r}>B$ to escape and equality to
turning/return. BEACH solves

$$
Q=Q_{\rm base}+C_{\rm esc}[B(Q)].
$$

Only when the barrier crosses one energy group does the solver introduce an escape weight $0\le\theta\le1$ for that group
and solve $B(Q)={\cal E}_k$ as a fractional marginal root. Equal-energy macro particles are therefore not split
arbitrarily into separate groups.

For $M$ groups, let $Q_k=Q_{\rm base}+C_k$. The solver follows connected paths from one common final-source root to
$Q_0$ and $Q_M$, numerically checking that the barrier is nondecreasing over the complete candidate-charge interval.
It then locates the first true value of

$$
P_k=
\begin{cases}
[B(Q_k)\ge{\cal E}_{k+1}], & 0\le k<M,\\
[B(Q_M)\ge{\cal E}_M], & k=M
\end{cases}
$$

by binary search. Pure-root selection therefore requires $O(\log M)$ connected candidate solves.
If $P_0$ is true, the result is all-return; if $P_M$ is false, it is all-escape. At an internal first-true index $k$,
$B(Q_k)<{\cal E}_k$ selects a pure root. Otherwise, the solver resolves a fractional marginal root in
$[Q_{k-1},Q_k]$. Marginal bisection continues to the charge tolerance and retains its upper endpoint so that
turning-point equality remains on the return side. A candidate Zhao root is not reselected by an independent Newton solve;
it is followed from the common root along this connected parameter path:

$$
\mathbf G\!\left(\mathbf y;E_I(\lambda),n_{pe,0}(\lambda)\right)=\mathbf 0,
\qquad
E_I(\lambda)=(1-\lambda)E_{I,0}+\lambda E_{I,1},
\qquad
n_{pe,0}(\lambda)=(1-\lambda)n_{pe,0}^{(0)}+\lambda n_{pe,0}^{(1)}.
$$

Here $\mathbf y$ contains the log-encoded root variables for one fixed A/B/C branch. Type A/B requires $E_I>0$, Type C
requires $E_I<0$, and this chart also rejects the nonregular $E_I=0$ endpoint. An adaptive pseudo-arclength
predictor/corrector advances from the preceding root while checking local correction distance, root jump, tangent
orientation, and Jacobian rank. A vanishing or reversed tangent in the $\lambda$ direction, or loss of rank in the
fixed-parameter Jacobian, is treated as a fold before the target and fails closed. A $\lambda=1$ event corrector lands at the
target. These checks are a local numerical guard against a nearby-root jump, not a mathematical proof that no arbitrarily
close sheet exists. This internal path adds no TOML key.

For charge candidates at fixed final source, BEACH also evaluates the tangent slope of $B$ with respect to prescribed
charge at every adaptively accepted point and stops on a negative slope or oppositely directed barrier secant. Endpoint
paths, order-statistic midpoints, and marginal bisection also check barrier ordering with charge. This remains a
finite-precision guard rather than an interval-arithmetic proof over the continuum between accepted points. Final
acceptance recomputes the measured-CDF escape charge and requires $Q-Q_{\rm base}-C_{\rm esc}=0$. Loss of rank, a fold,
a barrier-slope reversal, an inconsistent order-predicate bracket, or failure to bracket a marginal energy stops fail
closed without branch fallback.

Surface emission and interface-source normalization are separate quantities.
`photo_raycast.emit_current_density_a_m2` determines the number and macro weights of rays emitted from the rough surface.
The per-batch source amplitude passed to the Zhao density closure is instead resolved from charge that actually reaches the
interface:

$$
J_{pe,I}^{\rm meas}=\frac{\sum_r w_r}{A\Delta t},\qquad
s_{\rm eff}=
\frac{J_{pe,I}^{\rm meas}}
{|q_{pe}|n_{\rm ref}\sin(\alpha)\sqrt{2T_{pe}/m_{pe}}/(2\sqrt{\pi})}.
$$

Configured `photoelectron_source_scale` constructs the initial `zhao_floating` branch anchor; it is not reused as the
interface source for the same batch. The accepted $s_{\rm eff}$ is stored as the resolved source scale in the outer state
and checkpoint.

If the measured source differs from its previous-batch value, the connected path above first advances the previous root to
a common anchor on the final $s_{\rm eff}$ slice. Every charge candidate then holds the source fixed at that final value and
is followed from the common anchor on the committed A/B/C root sheet. A branch-label change, bootstrap branch, disconnected
root, or fallback to another closure is not allowed. Each continuation step limits dimensionless changes in $\phi_0$,
$\phi_{\min}$, and the infinity electron density to 0.25.

Because one cohort measured in the batch-start field is reused with the updated profile, the final state must also satisfy

$$
\frac{|q_e\Delta\phi_I|}{T_e}\le0.25,\quad
\frac{|\Delta B_{e,\mathrm{in}}|}{T_e}\le0.25,\quad
\frac{|q_i\Delta\phi_I|}{E_{i,n}}\le0.25,\quad
\left|\log\frac{n_{e,\infty}^{\mathrm{new}}}{n_{e,\infty}^{\mathrm{old}}}\right|\le0.25,\quad
|\log(s_{\rm eff}/s_{\rm old})|\le0.25,\quad
\frac{\lambda_{D,pe}|\Delta E_I|}{T_{pe}/|q_{pe}|}\le0.25,\quad
\frac{|\Delta B|}{T_{pe}}\le0.25.
$$

An over-limit update stops fail-closed; BEACH does not retrace the cohort in a new field or restore the old profile.
Here $B_{e,\mathrm{in}}=\max(0,q_e\phi_I,q_e\phi_{\min})$ and
$E_{i,n}=\max(T_i,m_i u_{i,n}^2/2)$. The first four terms monitor the absolute potential shift, the electron cutoff
including Type A, the ion drift-energy map, and the Zhao-resolved electron inflow density frozen into the ambient
reservoir cohort. Photoelectrons instead respond to the profile-relative barrier $B$ and field deformation, so the last
two terms retain $T_{pe}$ and $\lambda_{D,pe}$. The same $\Delta\phi_I$ is not tested a second time against $T_{pe}$.

Because this CDF is the raw measured cohort from one batch, convergence must use not only the total ray count but also the
number reaching the interface and the effective sample count left in the escaping tail above the barrier. If the
ambient-potential trust guard fails statistically while source, field, and barrier changes remain small, do not relax the
guard. Increase
`rays_per_batch` and confirm that $\phi_I$, escape fraction, and total surface charge agree at no less than twice the ray count.
If the same trust violation converges instead of disappearing as the ray count increases, the physical update in one batch
is too large rather than statistically noisy. Keep the guard fixed, reduce `batch_duration` and the corresponding emitted
charge per batch, and demonstrate time-step convergence.

After acceptance, the escape/return weight assigned by each measured energy group is applied directly to the emission
element and actual reabsorption element. The common `return_weight_scale` from the ambient path is not used. Each ray must
normally have exactly one outward crossing and at most one return. A return-weighted ray that is not absorbed, an escape
ray with a mismatched terminal outcome, or a recrossing stops the run instead of transferring weight to another ray.

The measurements provide only interface-source amplitude and the normal-energy CDF. Zhao free/reflected ambient-electron,
free/captured photoelectron, and cold-ion density closures remain analytic. This is not a model for tangential velocity
distributions, arbitrary trapped VDFs, collisions, magnetization, time-dependent Vlasov--Poisson, or PIC. Even Type A is
limited to its single virtual-cathode minimum and the reflected/captured populations represented by the Zhao equations.

The current MPI implementation gathers every interface sample with `MPI_Gatherv`, groups samples and performs the Zhao
solve on the root, then broadcasts the accepted step and profile. Root temporary memory and gather traffic scale with the
number of photoelectrons crossing the interface. Account for this root-local limit when increasing ray count for a
convergence study.

### Start stationary studies from the Zhao zero-current root

When the target is a strong-UV stationary observable rather than the turn-on transient from an uncharged state, select
`external_boundary.particles.steady_start_mode="zhao_floating"`. Before the first batch, this mode uses the configured infinity quasineutral
conditions, temperatures, drifts, and UV source to solve the Zhao zero-current stationary root. It constructs the
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
`external_boundary.particles.mode="same_batch"`, and it does not overwrite existing charge on a fresh run. With the same configuration and
`output.resume=true`, BEACH restores checkpoint mesh charge and a complete outer state without reseeding from the zero-current
root. Combining the warm start with the queue transient closure is rejected. UV turn-on, delayed return current, and transient
cloud inventory still require the queue or a dynamic outer model.
Even for a stationary publication result, a warm start alone does not demonstrate uniqueness or dynamic stability. Confirm
that an independently relaxed or perturbed seed returns to the same stationary observables.

[ADR 0006](adr/0006-zhao-stationary-warm-start.md) records this scope and its separation from the queue transient closure.

With `external_boundary.particles.mode="zhao_queue"`, tracked photoelectrons retain their macro-particle weights as queue
inventory while they occupy the outer region. The MPI-global photoelectron number divided by horizontal area is the target column. The closure solves a
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

[`periodic2_zhao_charge_driven_outer.toml`](../examples/periodic2_zhao_charge_driven_outer.toml) anchors at the stationary
Type-A `zhao_floating` root and exercises the strong-photoemission empirical-CDF update over small batches. Successful
completion checks the implementation path; it does not establish physical convergence in ray count, batch width, or
interface location.
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

`external_boundary.field.debye_length` and `external_boundary.field.thermal_voltage` are not physical scales for the Zhao root or profile. They remain
reference inputs for the split-interface applicability diagnostics: the former scales `interface_eta_gap` and the local-charge
estimate, while the latter scales the lateral `interface_eta_phi_kneq0` and `interface_eta_field_kneq0` diagnostics. Do not test
Zhao convergence by varying these inputs. Refine the profile grid, effective interface location, tracked-particle count, and
`dt` or batch resolution instead. The current 128-point profile grid is fixed at runtime, so a production grid-convergence study
requires exposing its point count as an input.

The analytic raw source associated with tracked `photo_raycast` at the effective plane is

$$
J_{pe,\mathrm{raw}}=s_{UV}
\frac{|q_{pe}|n_{\mathrm{ref}}\sin(\alpha)v_{\mathrm{th},pe}}{2\sqrt{\pi}},
\qquad
v_{\mathrm{th},pe}=\sqrt{\frac{2T_{pe}}{m_{pe}}}.
$$

For fixed-source Zhao paths other than `implicit_mean`, this value must agree within 1% with
`photo_raycast.emit_current_density_a_m2`. Here $s_{UV}$ is `photoelectron_source_scale` and $T_{pe}$ in the speed formula
is thermal energy in joules. Strong-photoemission `implicit_mean` does not impose the 1% surface-current match; it resolves
$s_{\rm eff}$ from measured interface current. The analytic raw current enters
the tracked-source consistency check and current-density diagnostics, but not the root, surface charge, or ledger; only tracked
emission and reabsorption update the latter two. The population scale $\eta$ does not scale that raw photoelectron
emission-current term in the current diagnostic.
The closure also rejects
`external_boundary.particles.inflow_model="infinity_barrier"` and
`external_boundary.field.photoelectron_density_model="kinetic_mean"`.

The fixed-source Zhao effective-plane approximation does not self-consistently connect tracked-ray directions or a VDF
reaching the interface from a rough surface to the Zhao outer population. Strong-photoemission `implicit_mean` measures only
source amplitude and the normal-energy CDF; Zhao density populations and the tangential VDF remain analytic.
`ray_direction` controls illumination-ray sampling of emitting surfaces, while $\alpha$ independently controls the analytic
Zhao density closure. See [ADR 0003](adr/0003-zhao-charge-driven-outer-closure.md) for scope and the
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
frozen-field ratio, and quasisteady applicability for species using individual profile return in
[Particle escape and return](ParticleEscapeReturn.en.html). For ambient `implicit_mean`, additionally test ray-count
dependence of mean-solver iterations, transaction residual, sampled escape fraction, and the common return weight scale.
For strong-photoemission empirical Zhao, inspect branch, $\phi_{\min}$, $n_{e,\infty}$, the full-profile barrier, resolved source scale,
marginal energy/fraction, empirical charge residual, recross fraction, and terminal-mismatch fraction.
For both `implicit_mean` paths, the frozen-field ratio diagnoses the shadow-orbit timescale and, unlike ordinary same-batch
particles, does not make an over-limit value alone a run failure.

## Code reference

- VDF closures and nonlinear Poisson solve: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- Charge-driven Zhao roots and nonmonotone profiles: [`bem_outer_plasma_zhao.f90`](../src/physics/outer_plasma/bem_outer_plasma_zhao.f90)
- Build solver options from runtime species: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)
- Ambient-mean backward-Euler update: [`bem_mean_charge_integrator.f90`](../src/physics/outer_plasma/bem_mean_charge_integrator.f90), [`bem_dynamic_k0_mean.f90`](../src/runtime/coupling/bem_dynamic_k0_mean.f90)
- Strong-photoemission empirical-energy CDF and nonlinear Zhao update: [`bem_dynamic_k0_zhao.f90`](../src/runtime/coupling/bem_dynamic_k0_zhao.f90)
- One-pass rays and normalized transaction: [`bem_mean_charge_transaction.f90`](../src/runtime/coupling/bem_mean_charge_transaction.f90), [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Stationary-root and plane-charge warm start: [`bem_zhao_steady_start.f90`](../src/runtime/coupling/bem_zhao_steady_start.f90)
- Surface-field coupling and MPI collective solve: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- Profile output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- Profile restart: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
