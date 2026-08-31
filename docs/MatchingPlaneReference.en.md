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
| Look up the exact header, 11 columns, units, and Cartesian grid | [Table-backend response CSV v1](#table-backend-response-csv-v1) |
| Build a response table with `beach-zhao-response` | [Build a response table for the table backend](#build-a-response-table-for-the-table-backend) |
| Look up the fixed-point acceptance and relaxation equations | [Fixed-point numerical contract](#fixed-point-numerical-contract) |
| Converge the grid, batch width, and matching-plane height | [Validate convergence and applicability](#validate-convergence-and-applicability) |

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
when they do not bracket a root. Online Zhao uses a guarded secant step with a midpoint fallback.

Online Zhao seeds the solve with the previous outer iteration's endpoint (initially $D_H^n$). For a valid seed it starts
from the explicit endpoint correction, capped by the natural sheath scale

$$
D_{ref}=\sqrt{\epsilon_0 n_i e T_e},
$$

and doubles the width at most 64 times. If the seed lies outside an explicit A, B, or C solution interval, BEACH scans
the branch-compatible sign in $D_{ref}/32$ increments out to $8D_{ref}$. It never forms a bracket across an uncertified
gap. This searches only the current batch endpoint; it does not persistently extend a table. The run stops if the branch
ends before the endpoint, the scan or numeric range is exceeded, or no sign change is found.

The signed scan applies only to an explicit `a`, `b`, or `c` selection. With the default
`zhao_root_selection="require_unique"`, the implicit solver stops instead of choosing a branch when `auto` cannot
certify a unique physical solution at the seed. Implicit integration therefore does not remove strong-PE A/B
coexistence; branchwise solvability still has to be validated.

With `zhao_root_selection="minimum_energy"`, BEACH evaluates every candidate detected by the multistart search over
the full profile from the surface to infinity:

$$
U=-\frac{\epsilon_0}{2}\int_0^\infty E^2\,dx.
$$

It selects the lowest $U$ within an explicit branch, or among numerically certified A, B, and C candidates for `auto`.
The solve stops when a numerical failure prevents certification of the candidate set or when the lowest values are
tied within a relative $10^{-6}$. This energy comparison is a candidate-selection rule; finite multistart search does
not guarantee enumeration of every root or prove time-dependent stability.
The response can be discontinuous where the minimum-energy root switches. If the backward-Euler residual crosses that
discontinuity without an ordinary zero, BEACH reports numerical failure instead of mixing the two roots.

With PEs, the half-Maxwellian reduction gives

$$
\Gamma_{pe}^{escape}(D)=\Gamma_{pe}^{out}
\exp\left[-\frac{\Phi_H(D)-\Phi_{pe,barrier}(D)}
{\langle K_{pe,n}^{out}\rangle}\right].
$$

With $q_{pe}<0$, BEACH solves the endpoint of

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

This normalization keeps `matching_plane_residual <= coupling_rtol` for an accepted trial even when an absolute
tolerance dominates a component.

History response columns contain values evaluated at the accepted trial's $X^m$; feedback columns contain the observed
$X_{raw}^{m+1}$ from that same trial. BEACH does not record an unexecuted relaxation update $X^{m+1}$ after convergence.

BEACH stops rather than using an unconverged value when `coupling_max_iterations` is exhausted, a table query leaves an active range,
or an online solve fails. See the [Output format reference](OutputReference.en.html#matching_plane_quasistatic) for the
accepted-state and residual output contract.

## Validate convergence and applicability

1. Vary `coupling_rtol`, `coupling_atol`, relaxation, and particle count and compare accepted observables.
2. Independently vary table-grid resolution and range, or the explicit online Zhao branch.
3. Converge `batch_duration`, mesh, and periodic-cell resolution.
4. Move $H$ within the overlap region and test invariance of grain charge, gap potential, and PE escape fraction.

Weak dependence on $H$ is the central validation specific to this coupling. Treat run completion, numerical convergence,
and physical validity as separate conclusions.
