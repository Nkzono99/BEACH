title: Matching-plane numerical and response-table reference

Lang: [日本語](MatchingPlaneReference.md) | [English](MatchingPlaneReference.en.md)

# Matching-plane numerical and response-table reference

This reference defines the response CSV, implicit mean-charge update, and fixed-point convergence contract for
`surface_current_model.model="matching_plane_quasistatic"`. For model selection, the first four-batch run, and output
diagnosis, start with [Couple an outer sheath at a matching plane](MatchingPlaneCoupling.en.html).

## Find a contract

| Need | Section |
|---|---|
| Update only the mean charge implicitly at a seconds-scale batch width | [`implicit_zero_mode`](#implicit_zero_mode) |
| Look up online-Zhao multiplicity and root-family policies | [`zhao_root_selection`](#zhao_root_selection) |
| Reproduce an implicit-endpoint failure with the same distribution | [Failure diagnostics and the PE distribution](#failure-diagnostics-and-the-pe-distribution) |
| Look up the exact header, 11 columns, units, and Cartesian grid | [Table-backend response CSV v1](#table-backend-response-csv-v1) |
| Build a response table with `beach-zhao-response` | [Build a response table for the table backend](#build-a-response-table-for-the-table-backend) |
| Look up the fixed-point acceptance and relaxation equations | [Fixed-point numerical contract](#fixed-point-numerical-contract) |
| Converge the grid, batch width, and matching-plane height | [Validate convergence and applicability](#validate-convergence-and-applicability) |

## Zhao Type B electron density

Type B ambient electrons transport the upstream distribution by energy conservation, retaining the lower velocity
bound of particles accelerated into positive potential. At zero drift, BEACH's normalization gives

$$
\hat n_{e,f}=\frac{\hat n_{e,\infty}}{2}
\operatorname{erfcx}\!\left(\sqrt{\hat\phi/\tau}\right).
$$

Here $\hat\phi=e(\Phi-\Phi_\infty)/(k_BT_{ph})$, $\tau=T_e/T_{ph}$,
and densities are normalized by the reference PE density. Nonzero drift uses the
[upstream-VDF orbit integral](#consistency-between-the-ambient-vdf-and-sheath-density), not a shifted erfc argument.
Regenerate tables produced before the density correction with `beach-zhao-response`; the reader uses existing CSV values as supplied.

## `photoelectron_closure`

For online Zhao, `moment_matched_half_maxwellian` (default) preserves the existing PE outward-flux and mean-normal-energy approximation
at H. `energy_spectrum` uses the outward normal-energy flux distribution $F_H(K)=d\Gamma_H/dK$ at H. It requires
`model="matching_plane_quasistatic"`, `response_backend="zhao_online"`, and an active `photoelectron_species` role.
It cannot be selected with the five-input table, stationary Zhao, or a PE-free case. Ambient electron/ion density
laws and A/B/C branch constraints remain unchanged.

### Distribution and currents

Outward crossings at H are recorded before the outer-barrier decision. Macro weights are divided by area and batch
duration to obtain flux. Recrossings count too. PEs absorbed below H do not enter this supply; exterior returns
reflect at H and resume the existing interior trajectory tracking. Surface emission is not added back to the H flux.

With $K$ in eV and $B=\max(0,\Phi_H-\Phi_{pe,barrier})$ the barrier energy in the same units,

$$
\Gamma_{escape}=\int_B^\infty F_H(K)\,dK,\qquad
\Gamma_{return,H}=\Gamma_H-\Gamma_{escape}.
$$

Exterior density uses this same distribution, dividing the accessible flux at each energy by its local speed.
The segment between H and a Type-A minimum, and Type B, include the returning population; beyond the minimum only
the transmitted population is retained. Binwise density and its potential integral are evaluated analytically for
piecewise-constant $F_H$. The existing finite search checks Sagdeev integrals, far neutrality, and profile E².
Failure to detect a root does not generally prove nonexistence.

For `implicit_zero_mode=true`, the PE escape target uses the integral above throughout; it does not retain a
mean-energy exponential approximation. The surface emission target remains the configured current, and the total
return target is surface emission minus this escape. Exterior return and total return to the surface are distinct.

### Grid, iteration, and restart

For `photoelectron_spectrum_bins_per_decade=N` (default 32, positive int32), bin edges are

$$
K_j=T_{pe,config}\left(10^{j/N}-1\right)\quad[j=0,1,\ldots].
$$

$T_{pe,config}$ is the configured PE temperature in eV. The grid extends to cover measured energies. Each bin stores
its integrated flux $f_j$, representing constant $F_H=f_j/(K_{j+1}-K_j)$ over that interval. The configured temperature
sets the grid scale; it does not reset the measured distribution to a Maxwell distribution. Spectrum-mode mean
energy is computed from these bins, so it differs from the former sample mean by the discretization error.

Initial startup or restart from an older checkpoint without a spectrum builds a Maxwell initial guess from
the available PE moments. Each trial relaxes the whole distribution with `coupling_relaxation`. In addition to the
existing moment convergence criteria, $\sum_j|f_{j,observed}-f_{j,input}|$ must be no larger than the PE-flux tolerance
$\max(\mathtt{coupling\_rtol}\,s_\Gamma,\mathtt{coupling\_atol}[1])$, where `[1]` denotes the first component.
The existing warning and acceptance policy for finite unconverged trials is unchanged.

`matching_plane_spectrum_history.csv` saves observed and response-input distributions separately. Checkpoint summary
records additionally retain both distributions, their grid, and the response input; existing scalar-history columns
are preserved. If the configured PE temperature and bin resolution match the saved grid, restart restores the saved
distribution, including its non-Maxwell shape. A changed grid produces a warning and a new Maxwell initial guess
with the saved flux and mean energy; subsequent trials measure the distribution again at H. The default moment
mode also continues when only its diagnostic grid changes. See the
[output reference](OutputReference.en.html#pe-spectrum-observations-and-response-inputs) for columns and units.

Check convergence and balance by changing resolution, ray count, batch duration, and H. This remains a planar,
collisionless, unmagnetized 1D exterior approximation. It adds neither correlations beyond normal energy, delayed
outer returns, nor PE volume space charge inside BEACH.

## `implicit_zero_mode`

When the explicit mean-current update becomes stiff at a seconds-scale `batch_duration`, set
`implicit_zero_mode=true` to update only the mean $D_H$ with backward Euler. Both table and online Zhao backends support it.

### Configuration contract

| Backend | $D_H$ search domain | Feedback |
|---|---|---|
| `table` | Within the CSV $D_H$ axis, which needs at least two nodes | Positive PE-flux and PE-energy singletons with PEs; both zero without PEs. Ambient-outward axes are zero singletons |
| `zhao_online` | Searched from the current state within the selected branch | PE moments update during each particle fixed-point iteration. Ambient-outward feedback is transparent |

Both require `periodic2.lower_boundary_model="e_bottom_zero"`. A table provides an audited finite domain; online Zhao
solves the built-in model directly without a CSV.

```toml
[periodic2]
lower_boundary_model = "e_bottom_zero"

[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
implicit_zero_mode = true
```

This online run needs neither `response_table_path` nor `matching_query.csv`. A query CSV is used only when
`beach-zhao-response` creates a separate fixed table snapshot.

### Backward-Euler endpoint

BEACH solves

$$
D_H^{n+1}=D_H^n+hJ(D_H^{n+1})
$$

by first bracketing a sign change. A table uses bisection between its CSV endpoints and stops rather than extrapolating
when they do not bracket a root.

Online Type B solves the implicit endpoint directly in the boundary potential $\Phi_H$ when `zhao_branch="b"`, or
when `auto` with positive inward electron drift permits only Type B physically. At each potential, upstream neutrality
determines the ambient density and the Sagdeev integral gives $D_H(\Phi_H)$. With `energy_spectrum`, the search also
samples bin edges and quarter points within each bin, avoiding a coarse displacement step that skips a narrow
$D_H$ solution interval. Finite search does not guarantee enumeration of every mathematical root. With a valid seed,
`continuation` selects a unique nearest root; without an initial seed it requires uniqueness. It does not choose
between multiple implicit endpoints by discovery order.
Type-B bisection keeps the residual tolerance and continues until no new midpoint is representable (at most 96 iterations), rather than stopping on interval width alone.

Other online conditions search in $D_H$ with a guarded secant step and midpoint fallback.
The previous outer iteration's endpoint (initially $D_H^n$) seeds the solve. For a valid seed it starts
from the explicit endpoint correction, capped by the natural sheath scale

$$
D_{ref}=\sqrt{\epsilon_0 n_i e T_e},
$$

and doubles the width at most 64 times. If the seed lies outside an explicit A, B, or C solution interval, BEACH scans
the branch-compatible sign in $D_{ref}/32$ increments out to $8D_{ref}$. It never forms a bracket across an uncertified
gap. This searches only the current batch endpoint; it does not persistently extend a table. The run stops if the branch
ends before the endpoint, the scan or numeric range is exceeded, or no sign change is found.
With `auto` + `continuation`, an unresolved initial point triggers a scan in both signs. A valid root found on an
interval seeds local searches there; unaccepted local seeds do not become the next batch's accepted state.

After finite endpoints establish a sign-changing bracket, a roundoff-scale residual miss no longer stops the run.
BEACH accepts the endpoint with the smaller residual and emits a warning. An absent bracket, non-finite response, or
missing physical solution still stops the run.

Implicit integration does not remove strong-PE A/B coexistence; branchwise solvability still has to be validated.

### `zhao_root_selection`

| Value | Scope | Rule |
|---|---|---|
| `require_unique` | Online Zhao | Require one physical root at each query. If `auto` cannot certify uniqueness, stop rather than select a branch |
| `minimum_energy` | Online Zhao | Choose the detected multistart candidate with the lowest full-sheath potential energy |
| `continuation` | Online Zhao and `implicit_zero_mode=true` | Locally track from the last accepted endpoint; fallback accepts the detected root nearest to the seed within the configured branch policy |

The default is `require_unique`. `continuation` is an opt-in, history-dependent policy for `auto` / `a` / `b` / `c`.
It is unavailable with explicit mean-charge updates, tables, and stationary Zhao.

Each branch uses at most eight initial guesses. Potential candidates use the PE flux-distribution median, PE mean
energy, $\epsilon_0 E_H^2/(en_i)$, electron temperature, and the cold-ion kinetic-energy limit. The initial ambient
density follows upstream neutrality for each candidate. Guesses do not use fixed values in V or m$^{-3}$.

With `zhao_root_selection="minimum_energy"`, BEACH evaluates every candidate detected by the multistart search over
the full profile from the surface to infinity:

$$
U=-\frac{\epsilon_0}{2}\int_0^\infty E^2\,dx.
$$

It selects the lowest $U$ within an explicit branch, or among numerically certified A, B, and C candidates for `auto`.
The solve stops when a numerical failure prevents certification of the candidate set or when the lowest values are
tied within a relative $10^{-6}$. This candidate-selection rule follows the sheath potential-energy comparison of
[Mishra et al. (2023)](https://academic.oup.com/mnras/article/520/1/233/6987684); finite multistart search does not
guarantee enumeration of every root or prove time-dependent stability.
The response can be discontinuous where the minimum-energy root switches. If the backward-Euler residual crosses that
discontinuity without an ordinary zero, BEACH reports numerical failure instead of mixing the two roots.

For a new run, `continuation` performs finite multistart with the same uniqueness requirement as `require_unique`.
It does not rank bootstrap roots by energy. Each batch trial starts from the last accepted endpoint. Once a valid
endpoint is found within the trial, that root seeds the next feedback iteration, including in the first batch.
BEACH returns to full multistart when Newton fails to converge, its result cannot be decoded, profile certification fails,
or the candidate makes a large jump. The local Newton result must have distance

$$
d=\max\left(\frac{|\Delta\Phi_H|}{T_{pe}},\frac{|\Delta\Phi_{min}|}{T_{pe}},
\left|\log\frac{n_{e,\infty}}{n_{e,\infty}^{seed}}\right|\right)
$$

at most 0.25. Potential differences use the PE temperature scale; the path minimum is $\phi_m$ for Type A, zero for
Type B, and $\Phi_H$ for Type C. After full multistart, BEACH selects the closest detected root with the same distance. Let the two
smallest distances be $d_1$ and $d_2$. If
$|d_2-d_1|\le10^{-6}\max(1,d_1)$, they are numerically indistinguishable; BEACH reports ambiguity instead of using
initial-guess order. Otherwise it accepts the closest root whether its distance is below or above 0.25. The 0.25 bound is
only the fast-path acceptance limit for local Newton; it is not a physical distance limit on the root family after full
multistart. The solve stops when full multistart finds no root or when root search or profile certification fails
numerically. If a probe immediately after a valid root reports no physical solution or a numerical failure, the implicit
solver still subdivides that interval so a coarse scan does not skip a root near a branch endpoint.
It preserves an explicit branch; `auto` considers certified A / B / C candidates.

This is not pseudo-arclength continuation. Full multistart and branch-boundary subdivision can reacquire a root lost by
local Newton, but they do not prove retention of the same physical family, locate a fold, or prove passage through one.
The bootstrap and fallback retain the finite-multistart
limitation: they do not enumerate every mathematical root.

Rejected implicit probes, unaccepted fixed-point trials, and rejected adaptive-batch trials do not commit their root
candidates to the accepted continuation state. Rejecting a trial discards its local seed; only an accepted endpoint
seeds the next batch. On restart, BEACH tries
to reconstruct the seed from the saved accepted response and falls back to the initial unique-root search if
reconstruction fails. See the [Output format reference](OutputReference.en.html#matching_plane_quasistatic) for saved
state and receipts.

With PEs and either moment closure or a table, the half-Maxwellian reduction gives

$$
\Gamma_{pe}^{escape}(D)=\Gamma_{pe}^{out}
\exp\left[-\frac{\max(0,\Phi_H(D)-\Phi_{pe,barrier}(D))}
{\langle K_{pe,n}^{out}\rangle}\right].
$$

`energy_spectrum` instead integrates the measured distribution above the barrier. In either case, with $q_{pe}<0$,
BEACH solves the endpoint of

$$
D_H^{n+1}=D_H^n+h\left[
q_e\Gamma_e^{in}(D_H^{n+1})+q_i\Gamma_i^{in}(D_H^{n+1})
-q_{pe}\Gamma_{pe}^{escape}(D_H^{n+1})\right].
$$

Without PEs, it removes the final PE term and uses

$$
J=q_e\Gamma_e^{in}+q_i\Gamma_i^{in}
$$

and creates no PE target.

Table implicit mode fixes the PE moments at the CSV singleton values. Online implicit mode solves this endpoint with the
current PE feedback $X^m$, tracks the same trial batch, relaxes the measured PE moments, and re-solves the endpoint on the
next iteration. The PE-return feedback and $D_H^{n+1}$ are therefore reconciled as a nested fixed point.

Only the mean $k=0$ charge is implicit. The elementwise $k\ne0$ distribution still comes from the batch-start field.
Therefore, using a width such as 6 s requires separate checks of local potential change, particle sampling,
root bracketing, and the physical range. See [How to choose `batch_duration`](BatchDurationStability.en.html)
for the comparison workflow.

## Failure diagnostics and the PE distribution

When an implicit endpoint cannot be found, inspect the search result together with the distribution actually supplied
to the response. A finite search detecting no root, a located algebraic root violating physical conditions, and numerical
failure to certify a root are different outcomes. Rejecting every detected candidate does not prove nonexistence in
unsearched intervals.

If the Type-B potential search cannot certify a root, it reports one of these classifications. Both return
`numerical_failure` and include `finite search, no absence proof`.

| Classification | Meaning |
|---|---|
| `search_unresolved` | No root was detected, or bisection failure, non-finite evaluation, or numerical profile-certification failure remains |
| `detected_roots_all_nonphysical` | Every detected candidate violated a physical condition, with no bisection, non-finite, or response-evaluation failures |

`grid` / `valid` count evaluated/valid main-grid points. `neutral`, `E2neg`, `nonfinite`, and `response_fail` classify
those main-grid rejections, excluding evaluations inside bisection. `brackets` / `bisect_fail`, `located`,
`physical_reject` / `profile_numeric`, and `accepted` count sign changes/bisection failures, detected candidates,
physical/numerical certification rejections, and distinct physical roots. `located` precedes deduplication and must not
be compared directly with `accepted`. The log also records the potential search range, smallest `min_abs_F` [C/m2]
and its potential, and the first rejected potential with its `reason`. `min_abs_F=-1` means no valid evaluable point existed.
`neutral` means no positive ambient density satisfies upstream neutrality; `E2neg` means the interface field squared is negative.
`bisect_fail` separates invalid midpoints (`mid_invalid`) from residual-tolerance misses (`tol_miss`). `bisect_best_F`
is the smallest absolute residual during bisection and `F_tol` is the acceptance threshold, both in C/m2; a tolerance
miss is distinct from physical rejection.
`bisect_best_F=-1` means no bisection was performed. `min_abs_F` also includes valid evaluations inside bisection.

On failure, standard error records $D_{before}$, $D_{seed}$, duration, feedback, and the continuation seed's potentials
and density when valid. If a spectrum is attached, all bins are written as CSV between
`matching-plane failed spectrum begin` and `matching-plane failed spectrum end`, with this header:

```csv
energy_low_ev,energy_high_ev,flux_m2_s
```

`flux_m2_s` is the integrated number flux within a bin, not $d\Gamma/dK$. Empty bins are retained. Extracting the header
and numeric rows reproduces the failed query's PE distribution in double precision. The last accepted history or
checkpoint can contain a different distribution and must not be substituted. Obtain the other inputs, including
upstream conditions, from that run's `beach.toml` and failure log. Moment closure with no attached spectrum writes no
CSV block. These diagnostics do not update the checkpoint or accept the failed trial.

To recheck Type B independently, run the following from the repository root. This optional analysis tool requires
Python 3.10 or later, NumPy, and SciPy (plus `tomli` on Python 3.10). On HPC systems, run it on a compute node.

```bash
python tools/diagnose_matching_plane_failure.py run.err --config beach.toml --output diagnosis/batch14
```

It reads the last failed input and writes `diagnosis/batch14.json` plus `-spectrum.csv`, `-scan.csv`, `-profiles.csv`,
and `-intervals.csv`. JSON `diagnosis=physical_BE_endpoint_found` means a candidate satisfying the physical checks was
found for the same input. `type_B_absence_certified_in_model_with_roundoff_margin` sets
`absence_certificate.complete=true` only when conservative bounds exclude the entire Type-B potential domain.
This uses float64 and explicit roundoff margins, not a rigorous directed-rounding interval proof. Check
`configuration.all_configured_branches_covered` before extending that conclusion to every configured branch.
The other diagnoses, candidate rejection or no root located by finite search, do not prove nonexistence.

## Consistency between the ambient VDF and sheath density

Ambient-electron density transports the same upstream drifting Maxwellian used for inward flux by energy conservation.
With potentials in V, $T_e$ in eV, and inward drift $u=v_d/\sqrt{2eT_e/m_e}$, the Type-A population crossing
the potential minimum has density

$$
\frac{n_e(\phi)}{n_{e,\infty}}=
\frac{1}{\sqrt\pi}\int_{\sqrt{-\phi_m/T_e}}^\infty
\frac{s\,e^{-(s-u)^2}}{\sqrt{s^2+\phi/T_e}}\,ds.
$$

Here $s$ is the upstream inward normal speed divided by $\sqrt{2eT_e/m_e}$, with $\phi_m<0$ and $\phi\ge\phi_m$.
Type B sets the lower bound to zero; reflecting intervals add the reflected population from the same upstream VDF.
Nonzero drift is not approximated by a Boltzmann factor times shifted erfc. PE density likewise preserves the source
orbits and passing/reflected velocity ranges. Type A requires $\phi_m<\min(0,\Phi_H)$ and permits $\Phi_H<0$.

Roots of neutrality and Sagdeev integral constraints are insufficient: the whole profile must satisfy $E^2\ge0$.
Type B also checks the $\sqrt{\phi}$ term of the charge density near infinity, rejecting a negative-$E^2$ interval
that a finite profile grid can miss. This coefficient uses the PE flux density at the barrier, with distinct left and
right values at a bin edge.
Positive inward drift with complete reflection of slow ambient electrons in Type A / C violates this condition near
infinity under strict neutral, zero-field upstream conditions. A finite upstream boundary defines a different boundary-value
problem and is not part of the current online closure.

## Table-backend response CSV v1

Declare the matching-plane height once before the header. It must equal the z component of `domain.box_max`.

```csv
# matching_plane_z_m=1.0e-3
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,electron_outward_number_flux_m2_s,ion_outward_number_flux_m2_s,matching_potential_v,electron_inward_number_flux_m2_s,ion_inward_number_flux_m2_s,electron_access_potential_v,ion_access_potential_v,photoelectron_barrier_potential_v
```

The first five columns are input axes and the last six are responses.

| Column | Unit | Meaning |
|---|---|---|
| `displacement_c_m2` | C/m2 | Mean $D_z$ immediately below the interface; +z is positive |
| `photoelectron_outward_number_flux_m2_s` | 1/(m2 s) | Outward PE flux reaching the interface |
| `photoelectron_outward_mean_normal_energy_ev` | eV | Mean normal kinetic energy of outward PEs |
| `electron_outward_number_flux_m2_s` | 1/(m2 s) | Outward ambient-electron flux |
| `ion_outward_number_flux_m2_s` | 1/(m2 s) | Outward ion flux |
| `matching_potential_v` | V | Interface potential $\Phi_H$ returned by the outer sheath |
| `electron_inward_number_flux_m2_s` | 1/(m2 s) | Total electron flux into BEACH, optionally including outer return |
| `ion_inward_number_flux_m2_s` | 1/(m2 s) | Total ion flux into BEACH, optionally including outer return |
| `electron_access_potential_v` | V | Access bottleneck from the electron reservoir to the interface |
| `ion_access_potential_v` | V | Access bottleneck from the ion reservoir to the interface |
| `photoelectron_barrier_potential_v` | V | Maximum outer barrier faced by outward PEs |

### Grid and value contract

- Include the complete Cartesian product of the five input axes, with no duplicate or missing points; row order is arbitrary.
- Fluxes, PE mean energy, and output fluxes are nonnegative, and every value is finite.
- A feedback axis with at least two nodes includes zero for the initial query. BEACH does not extrapolate beyond its range.
- A response-independent feedback axis is a singleton. It accepts any finite query and disables dependence on that input.
- All four potential columns use the same upstream-0-V gauge.
- Numeric tokens are decimal reals; do not use Fortran list-directed controls such as `/`, `2*0`, or null fields.

Interpolation is multilinear over at most 32 corners. Load memory is linear in row count, and every MPI rank keeps the table.

### Build a response table for the table backend

Build a production table with an independent Zhao or 1-D PIC sweep that uses the same $H$, upstream distributions, and
sign conventions. Feed it PE flux and normal energy that actually cross the matching plane, not emission at the wall.
For a nonmonotonic potential, use the maximum barrier over the complete outer profile.
Store the table-generation code, upstream conditions, solver version, and unit conversions with the production data.

To pre-evaluate the built-in online Zhao response and write a table-format snapshot, run:

```console
beach-zhao-response \
  examples/periodic2_matching_plane_zhao_online.toml \
  examples/matching_plane_zhao_query_grid.csv \
  response.csv
```

The configuration must be a complete matching case with `response_backend="zhao_online"` and no `response_table_path`.
Because a table has no accepted-endpoint history, generation rejects `zhao_root_selection="continuation"`. Use
`require_unique` or `minimum_energy` when generating a table.
The query CSV permits blank lines and `#` comments. Its first noncomment line is the exact header below. Every value is
finite, and fluxes and PE energy are nonnegative.

```csv
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,electron_outward_number_flux_m2_s,ion_outward_number_flux_m2_s
```

For the v1 generator, include zero on the PE-flux axis and make PE energy a singleton. If the PE-flux axis has a positive
node, the energy must also be positive. Use zero singletons for the two transparent ambient-outward axes.

Supply the complete five-axis product. The CLI writes the 11-column `response.csv` only after every query solves
successfully. The sample grid uses a fixed 3 eV and checks wiring; it is not a production range. Generate a production
table with PE-energy dependence directly from an independent outer solver.

If `beach-zhao-response` is unavailable, install the
[current version described by this site](Installation.en.html#install-the-version-described-by-this-site). For a run with
the generated table, create another configuration with `response_backend="table"` and `response_table_path="response.csv"`.
A run that continues to use the online backend needs no response table.

### Map Zhao solvability

`beach-zhao-atlas` evaluates Zhao A, B, and C independently at each supplied matching-plane moment. Use this offline
diagnostic before selecting a simulation branch when you need to distinguish multiple roots, absence of a physical
root, and numerical solver failure. It does not build a response table or change the BEACH runtime configuration.

```console
beach-zhao-atlas \
  examples/periodic2_matching_plane_zhao_online.toml \
  query_grid.csv \
  atlas.csv
```

Supply a complete matching case with `response_backend="zhao_online"`. Regardless of its configured `zhao_branch`
and `zhao_root_selection`, the atlas evaluates A, B, and C separately with `require_unique`. The query CSV permits blank
lines and `#` comments. Its first noncomment line is the exact header below. A full Cartesian product is not required;
you may list only the points to diagnose.

```csv
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev
```

`atlas.csv` contains three rows per query, one for each branch. Interpret `status` as follows.

| `status` | Meaning |
|---|---|
| `ok` | The solver certified one physical root in that branch |
| `no_physical_solution` | The current solver certified that branch as physically inadmissible |
| `numerical_failure` | The solver could not certify either existence or absence |
| `ambiguous_within_branch` | More than one root remains within the same branch |
| `invalid_input` | A flux or PE-energy value violates the input contract |

Two or more `ok` branches make the query `multiple`. If exactly one branch is `ok` but another is
`numerical_failure` or `ambiguous_within_branch`, uniqueness remains uncertified. Classify a query as `no_root` only
when all three branches are `no_physical_solution`; never merge numerical failures into `no_root`.
This is a solver certificate from the current finite set of initial guesses and profile checks, not a mathematical
proof that no root exists.

## Fixed-point numerical contract

The feedback vector is ordered as

$$
X=(\Gamma_{pe}^{out},\langle K_{z,pe}\rangle^{out},\Gamma_e^{out},\Gamma_i^{out}).
$$

For active component $j$, with backend scale $s_j$, relative tolerance $r$, and absolute tolerance $a_j$, BEACH accepts
a trial when every component satisfies

$$
|X_{raw,j}^{m+1}-X_j^m|\le\max(r s_j,a_j).
$$

Otherwise, with `coupling_relaxation` $\alpha$, it updates

$$
X^{m+1}=(1-\alpha)X^m+\alpha X_{raw}^{m+1}.
$$

Inactive components are excluded and their `coupling_atol` entries must be zero.

The backend defines the scale and inactive components as follows.

| Backend | $s_j$ | Inactive components |
|---|---|---|
| `table` | Maximum minus minimum of the corresponding active feedback axis | Singleton feedback axes |
| `zhao_online` | Reference flux or reference energy defined by the Zhao model | Transparent ambient-electron and ion outward axes |

With $\Delta_j=X_{raw,j}^{m+1}-X_j^m$, BEACH reports `matching_plane_residual` as

$$
\max_j \rho_j,\qquad
\rho_j=
\begin{cases}
r|\Delta_j|/a_j, & a_j>r s_j,\\
|\Delta_j|/s_j, & a_j\le r s_j.
\end{cases}
$$

This normalization keeps `matching_plane_residual <= coupling_rtol` for a converged trial even when an absolute
tolerance dominates a component.

History response columns contain values evaluated at the accepted trial's $X^m$; feedback columns contain the observed
$X_{raw}^{m+1}$ from that same trial. BEACH does not record an unexecuted relaxation update $X^{m+1}$ after convergence.

When `coupling_max_iterations` is exhausted, BEACH commits the final trial with a warning if its feedback and response
remain finite. `matching_plane_residual > coupling_rtol` together with the maximum iteration count is the
nonconvergence receipt for that batch, and the next batch starts from the observed feedback. BEACH still stops when a
table query leaves an active range, an online solve fails, or a non-finite value leaves no valid trial. See the
[Output format reference](OutputReference.en.html#matching_plane_quasistatic) for the state and residual output contract.

## Validate convergence and applicability

1. Vary `coupling_rtol`, `coupling_atol`, relaxation, and particle count and compare accepted observables.
2. Independently vary table-grid resolution and range, or the explicit online Zhao branch.
3. Converge `batch_duration`, mesh, and periodic-cell resolution.
4. Move $H$ within the overlap region and test invariance of grain charge, gap potential, and PE escape fraction.

Weak dependence on $H$ is the central validation specific to this coupling. Treat run completion, numerical convergence,
and physical validity as separate conclusions.
