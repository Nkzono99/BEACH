title: Use Quasistatic Matching-Plane Coupling

Lang: [English](MatchingPlaneCoupling.en.md) | [日本語](MatchingPlaneCoupling.md)

# Use Quasistatic Matching-Plane Coupling

`surface_current_model.model="matching_plane_quasistatic"` treats the top of the BEACH box,
the z component of `domain.box_max`, as the interface to an outer one-dimensional sheath. BEACH retains every surface charge and the
`k=0` and `k!=0` fields below that interface. The selected response backend supplies the matching potential,
ambient inward fluxes, and access or return barriers. Particle trajectories return outward fluxes to a fixed-point
iteration in every batch, and BEACH commits only a converged trial to the surface charge.

`response_backend="table"` (the default) interpolates a response table produced by an independent Zhao or one-dimensional
PIC parameter sweep. `response_backend="zhao_online"` solves a finite-$H$, charge-driven Zhao A/B/C quasistatic response
from the current $D_H$ and outward PE moments. The latter does not reevaluate a stationary Zhao wall potential and does
not impose zero current. Validate the selected backend over the intended range. `examples/matching_plane_response_synthetic.csv`
is synthetic data for exercising the table path; it is not valid input for a physical study.

## Minimal configuration

**Prerequisite:** Prepare a `periodic2` case that is periodic in x and y, open in z, and contains every mesh vertex
strictly below the top face. Set the external uniform electric and magnetic fields to zero. Define distinct electron and ion
species, plus a photoelectron species only when PE is modeled.

**Action:** Add `[surface_current_model]`.

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "table"
response_table_path = "matching_plane_response.csv"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.0, 0.0, 0.0] # default: relative tolerance only
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

Omitting `response_backend` also selects `"table"`. See `examples/periodic2_matching_plane_quasistatic.toml` for a
complete synthetic smoke case.
BEACH resolves `response_table_path` relative to the directory containing `beach.toml`.

For online Zhao, select the backend and branch instead of supplying a response-table path.

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.05, 0.0, 0.0] # PE-energy tolerance [eV] for finite-ray sampling
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

See `examples/periodic2_matching_plane_zhao_online.toml` for a complete online configuration example.
The second component's `0.05` eV is an example for finite-ray sampling. Select it from a convergence study that varies
the macro-particle count while retaining relative tolerances on the flux components.

Without PE, omit both `photoelectron_species` and the `photo_raycast` species. The response and upstream-reservoir
potential gauge remain active, while PE feedback, return, and escape stay exactly zero.

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
```

**Expected output:** Configuration validation accepts the model and backend wiring. The table backend also validates the
response-table syntax and matching-plane height at startup. `summary.txt` records the latest accepted coupling state,
response backend, two or three species roles, and iteration controls; the table backend also records the response-table content
fingerprint. When
`output.history_stride > 0`, `matching_plane_history.csv` tracks each selected accepted batch's response values, four
feedback moments, return and escape fluxes, iteration count, and normalized residual.

`beachx lint` does not read response-file contents.

**Interpretation:** Successful completion establishes backend evaluation and numerical fixed-point convergence. The table
backend also establishes valid syntax and an in-range interpolation. It does not establish physical validity of the
outer sheath, invariance to matching-plane height, or Monte Carlo convergence.

## Configuration constraints

The model accepts only the following configuration so that no mean field or particle channel is counted twice:

- `field_boundary.mode="periodic2"`
- `periodic2.nonzero_mode_backend="cached_kneq0"` or `"panel_spectral_reference"`
- `periodic2.zero_mode_policy="exclude_k0"`
- `periodic2.lower_boundary_model="symmetric_vacuum"` or `"e_bottom_zero"`
- periodic x/y particle faces and open z-low/z-high faces
- `sim.e0 = [0, 0, 0]` and `sim.b0 = [0, 0, 0]`
- `[reservoir].inflow_model="source_vdf"` (the generic `infinity_barrier` is disabled)
- only the distinct electron, ion, and optional photoelectron roles enabled, each with
  `surface_charge_closure="explicit"`
- ambient electron and ion species with `source_mode="volume_seed"`, `npcls_per_step=0`, and only z-high
  `boundary_inflow="reservoir"`
- negative electron charge and positive ion charge; when specified, the photoelectron is negative, uses
  `photo_raycast`, `inject_face="z_high"`, and `deposit_opposite_charge_on_emit=true`, with an open z-high boundary
- `particle_boundary.ordinary_open_model="escape"`
- `sim.multiple_box_events_policy="abort"` or `"soft_discard"` bounded by count grace, cumulative fraction, and
  absolute charge

`soft_discard` locally removes a rare macro particle whose periodic boundary events cannot be completed. BEACH does not
add events or matching-plane moments from the unresolved step; moments committed by its earlier steps remain. Let $D$ be
the cumulative discard count and $P$ the cumulative number of macro-particles processed in accepted batches. The run
stops when $D$ exceeds `multiple_box_events_soft_discard_count_grace` and $D/P$ exceeds
`multiple_box_events_soft_discard_fraction_limit`. Exceeding the independent
`multiple_box_events_soft_discard_abs_charge_limit` also stops the run; equality at any threshold is allowed. This
absolute-charge limit is the physical error budget. Verify `multiple_box_events_soft_discarded`,
`multiple_box_events_soft_discard_fraction`, and `multiple_box_events_soft_discarded_abs_charge_C` in the
summary/checkpoint. Because a lifetime cumulative fraction can dilute a late burst, also audit the per-batch aggregate
log.

Do not combine this model with `reference_area_m2`, stationary-Zhao source keys, or manual `fixed_current` targets. The
interface area is the x-y area of `domain`, its height is the z component of `domain.box_max`, and its update interval is
one accepted batch. These are derived rather than configurable so that mutually inconsistent values cannot be entered.

The backend-specific constraints are:

- `response_backend="table"` requires `response_table_path` and forbids `zhao_branch`;
- `response_backend="zhao_online"` forbids `response_table_path` and accepts
  `zhao_branch="auto"`, `"a"`, `"b"`, or `"c"`; and
- the online backend requires singly charged species for every role, $T_e>0$, $0\le T_i\le0.1T_e$, positive ion
  number density, and positive inward ambient-electron and ion drift speeds at z-high (`drift_velocity` has a negative
  z component). When PE is specified, it additionally requires equal ambient-electron and PE masses and $T_{pe}>0$.

## Table backend: Response CSV v1

Declare the matching-plane height exactly once before the header. The value must be finite and agree with
the z component of `domain.box_max`.

```csv
# matching_plane_z_m=1.0e-3
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,electron_outward_number_flux_m2_s,ion_outward_number_flux_m2_s,matching_potential_v,electron_inward_number_flux_m2_s,ion_inward_number_flux_m2_s,electron_access_potential_v,ion_access_potential_v,photoelectron_barrier_potential_v
```

The first five columns are input axes, and the final six columns are responses.

Use one common gauge for `matching_potential_v` and all three potential outputs, with the upstream reservoir fixed at
0 V in the outer model. BEACH maps the inward electron and ion VDFs from that 0 V reservoir through the access
potential to the matching potential, so an arbitrary constant cannot be added only to the potential columns.

| Column | Unit | Meaning |
| --- | --- | --- |
| `displacement_c_m2` | C/m2 | Mean $D_z$ immediately below the interface, from BEACH surface charge and the lower closure; positive points in +z |
| `photoelectron_outward_number_flux_m2_s` | 1/(m2 s) | Outward PE number flux that reaches the interface |
| `photoelectron_outward_mean_normal_energy_ev` | eV | Mean normal kinetic energy of the outward PE population |
| `electron_outward_number_flux_m2_s` | 1/(m2 s) | Outward ambient-electron number flux |
| `ion_outward_number_flux_m2_s` | 1/(m2 s) | Outward ion number flux |
| `matching_potential_v` | V | Interface potential $\Phi_H$ returned by the outer sheath |
| `electron_inward_number_flux_m2_s` | 1/(m2 s) | Total electron number flux entering BEACH, including particles returned by the outer sheath |
| `ion_inward_number_flux_m2_s` | 1/(m2 s) | Total ion number flux entering BEACH, including particles returned by the outer sheath |
| `electron_access_potential_v` | V | Access bottleneck for the VDF entering from the electron reservoir |
| `ion_access_potential_v` | V | Access bottleneck for the VDF entering from the ion reservoir |
| `photoelectron_barrier_potential_v` | V | Maximum outer potential barrier faced by outward PEs |

Write one row for every point in the complete Cartesian product of the five input axes. Rows may be in any order, but
duplicate or missing points are invalid. Load-time memory is linear in the row count, and every MPI rank retains the
table; estimate the required memory before adding axis nodes. BEACH applies multilinear interpolation over at most 32 corners and fails on an
out-of-range query along any axis that has more than one node.
Numeric tokens must be decimal reals; Fortran list-directed controls such as `/`, `2*0`, and null fields are rejected.

Fluxes, mean energy, and response fluxes must be nonnegative. Feedback axes 2--5 with two or more nodes must include zero
so that the initial state can be evaluated. Use one node for an input axis that the outer sweep intentionally does not
model. A singleton node may be nonzero; its axis accepts any finite query and explicitly disables that feedback dimension.
This is neither clamping nor implicit
extrapolation.

### Update only the mean charge implicitly

When the explicit mean-current update becomes stiff at a seconds-scale `batch_duration`, the table backend can set
`implicit_zero_mode=true`.

```toml
[periodic2]
lower_boundary_model = "e_bottom_zero"

[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "table"
response_table_path = "matching_plane_type_a_response.csv"
implicit_zero_mode = true
```

The response table for this mode must have at least two `displacement_c_m2` nodes and singleton values on the other four
input axes. With PE, the PE singleton pair must be a positive outward flux and positive mean normal energy; without PE,
both must be zero. The outward ambient-electron and ion singleton values are zero in either case. BEACH solves
$D_H^{n+1}=D_H^n+hJ(D_H^{n+1})$ by bisection within the table's $D_H$ range. It stops instead of extrapolating when the
range does not bracket the root.

Particle tracking and the elementwise `k!=0` distribution still use the batch-start state. BEACH normalizes total
electron/ion absorption to the implicit response. With PE, it also normalizes PE emission to the configured
surface-emission current and PE return to the surface-emission flux minus the outer escape flux. Regression tests cover
endpoint consistency for PE, no PE, and strongly cancelling positive and negative charge. At runtime, BEACH adopts the
finite committed $Q/A$ from the compensated sum of mesh charge as the state for the next batch. Without PE, the
backward-Euler current is only
$q_e\Gamma_e^{in}+q_i\Gamma_i^{in}$ and no PE target is created.
A large value such as 6 s therefore also requires a bracketed response, acceptable
local `k!=0` potential change, adequate macro-particle sampling, and a physically valid response table. This option does
not make an arbitrary `batch_duration` unconditionally stable.

## Online Zhao backend

`response_backend="zhao_online"` converts each query's $D_H$ to $E_H=D_H/\epsilon_0$ and solves a Sagdeev A/B/C root
that connects the matching plane $H$, the z component of `domain.box_max`, to an upstream reservoir at 0 V with zero
field. `zhao_branch="auto"` searches the applicable branches; `"a"`, `"b"`, and `"c"` select one branch explicitly.
`auto` fails closed if it detects more than one physical root or if a numerical failure in a compatible branch prevents a
unique selection. Select `a`, `b`, or `c` explicitly for branch-resolved studies.
The v1 multiple-root check clusters roots found by a finite multistart set; it is not a mathematical root-isolation proof.
For branch-resolved validation, scan the parameters separately with each explicit branch.
This is a finite-$H$, charge-driven boundary response, unlike the wall zero-current root of `zhao_stationary`. It does
not impose $J=0$ while the surface charge is evolving.

$H$ fixes the origin of the outer half-space interface, the `periodic2` zero-mode gauge, and the plane where outward PE
moments are measured. Because online Zhao is planar and translationally symmetric, the absolute coordinate of $H$ is not
a numerical parameter of the Sagdeev equation. This model does not solve a wall-to-$H$ one-dimensional distance
constraint.

The outward PE number flux and mean normal kinetic energy are mapped to the amplitude and energy scale of a
half-Maxwellian that reproduces those two moments. For a zero-PE-flux query, the PE population remains zero. The
configured PE temperature is used as a numerical scale when that role exists; otherwise the ambient-electron temperature
is used. This reduction does not retain energy-tail information beyond the two moments.

The online MVP treats the outward ambient-electron and ion axes as transparent. Their measured values do not alter the
outer profile or return flux and are inactive in the fixed-point residual. The active scales are derived from the Zhao
model's reference fluxes and energy. Do not use this MVP for a case that requires outer return of ambient populations.

The solver is stateless between queries. It does not retain an outer-particle inventory, continuation from the previous
root, an outer flight time, or a delayed-return queue. If the configured branch policy finds no physical root, the
Sagdeev integral is not real, or the nonlinear solve does not converge, BEACH fails closed. It does not silently change
an explicitly selected branch or switch backend.

This backend assumes the planar, collisionless, unmagnetized Zhao VDF family. It is not a full-VDF solver, a 1D PIC
solver, or a time-dependent outer-sheath model. `solar_elevation_deg`, `photoelectron_ref_density_m3`, and
`photoelectron_source_scale` describe the wall source of `zhao_stationary` and must not be set for online matching.

## Fixed-point iteration

Let $X^0$ be the accepted feedback state. The feedback vector and `coupling_atol` use the order

$$
X=(\Gamma_{pe}^{out},\langle K_{z,pe}\rangle^{out},\Gamma_e^{out},\Gamma_i^{out}),
$$

with units [m^-2 s^-1], [eV], [m^-2 s^-1], and [m^-2 s^-1], respectively. `coupling_atol` is a finite,
nonnegative four-vector. Its default `[0, 0, 0, 0]` retains relative-only convergence tests.

Each trial performs these steps:

1. Compute $D_H$ from the current surface charge and evaluate the selected response backend at $(D_H, X^m)$.
2. Apply the response $\Phi_H$ as the `periodic2` zero-mode gauge, the total electron/ion inward fluxes and access maps,
   and the outward PE barrier. Ambient outflow is not reflected locally; it is passed to the next outer response.
3. Replay the particle trial from the same batch-start RNG state and macro-particle residuals, measuring outward moments
   $X_{\mathrm{raw}}^{m+1}$.
4. For every active component $j$, let $s_j$ be the backend scale, $r$ be `coupling_rtol`, and $a_j$ be its absolute
   tolerance. Accept the trial when every component satisfies
   $|X_{raw,j}^{m+1}-X_j^m|\le\max(r s_j,a_j)$.
5. Otherwise, with `coupling_relaxation` $\alpha$, calculate
   $X^{m+1}=(1-\alpha)X^m+\alpha X_{\mathrm{raw}}^{m+1}$ and replay the trial.

Inactive components do not participate in convergence, but their `coupling_atol` entries must be zero. This applies to
singleton feedback axes in a table and to the outward ambient-electron and ion axes in online Zhao. BEACH rejects a
nonzero absolute tolerance on an inactive component no later than provider initialization.

History response columns contain values evaluated at the accepted trial's $X^m$; feedback columns contain the measured
$X_{\mathrm{raw}}^{m+1}$ from that same trial. For each component, `matching_plane_residual` uses
$r|X_{raw,j}^{m+1}-X_j^m|/a_j$ when $a_j>r s_j$, and $|X_{raw,j}^{m+1}-X_j^m|/s_j$ otherwise, then records the
maximum. This conversion preserves `matching_plane_residual <= coupling_rtol` as the acceptance receipt under mixed
tolerances. No unexecuted relaxation step is applied after convergence.

The iteration does not commit surface charge, RNG state, injection residuals, ledger entries, or histories. A rejected
adaptive batch-duration trial also restores the outer state to the batch start. If `coupling_max_iterations` is exhausted
or a table query leaves the interpolation grid, or if the online solve fails, BEACH stops instead of continuing with an
unconverged value.

## Build a response table for the table backend

Sweep all five input axes in the one-dimensional Zhao or PIC model using the same `z=H`, upstream distributions, and
sign conventions. Feed the table the PE flux and normal energy that actually cross the matching plane, not the emission
flux at the material wall. For a nonmonotonic potential, output the maximum barrier over the complete outer profile, not
only a local potential difference.

Cover the expected initial state, stationary state, and intermediate iteration states with margin. Keep axes that the
outer model does not depend on as singletons instead of refining all five dimensions uniformly. A wider numeric grid
does not extend the physical validity of the outer model. Store the table-generation code, upstream conditions, solver
version, and unit conversions with production data.

To evaluate the built-in online Zhao backend over the same five inputs and save a snapshot for the table backend, run:
If `beach-zhao-response` is unavailable in a PyPI distribution, first
[install the current GitHub version](Installation.en.html#install-the-version-described-by-this-site).

```console
beach-zhao-response \
  examples/periodic2_matching_plane_zhao_online.toml \
  examples/matching_plane_zhao_query_grid.csv \
  response.csv
```

This sample grid has PE-flux nodes `0` and `1.0e10` at a fixed PE energy of 3 eV, with only zero on both
ambient-outward axes. It is a minimal wiring and snapshot-generation example, not a production input range.

`beach.toml` must be a complete `model="matching_plane_quasistatic"` configuration with
`response_backend="zhao_online"`. As with a direct online run, do not set `response_table_path` in this configuration.

The query CSV permits blank lines and comments beginning with `#`. Its first noncomment line must be this exact header:

```csv
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,electron_outward_number_flux_m2_s,ion_outward_number_flux_m2_s
```

Every data row contains five finite values; fluxes and PE mean energy are nonnegative. The PE-flux axis must include zero,
and the PE-energy axis must be a singleton. If the PE-flux axis has a positive node, that energy node must also be positive.
The transparent outward ambient axes 4 and 5 must each be a singleton zero node. Enumerate one complete
Cartesian product of the five axes without duplicate or missing points. Row order is arbitrary. The CLI evaluates every query with
the built-in Zhao backend before writing. Only when all queries solve successfully does it write `response.csv` with
`# matching_plane_z_m=...` and the existing 11-column contract. If any query fails, it does not generate the response
table.
This v1 generator produces a flux / $D_H$ response curve at one fixed PE energy. PE-energy dependence requires a future
table format that supports a conditional grid.

Use a separate configuration file for a run that consumes the generated table, and switch the backend as follows:

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "table"
response_table_path = "response.csv"
```

A run that continues to use the built-in online Zhao backend directly does not need a response table.

## Validate the result

**Numerical validation:** Vary these inputs independently and require the reported observables to remain within the
study's error budget:

- `coupling_rtol`, `coupling_atol`, `coupling_relaxation`, and `coupling_max_iterations`
- response-grid resolution and range, or the online Zhao branch selection
- `sim.batch_duration` and macro-particle count
- mesh and periodic-cell resolution

**Coupling validation:** For the table backend, regenerate the outer model for the same physical conditions. For either
backend, move `H` within an overlap region. Grain charge, gap potential, and PE escape fraction should be nearly
invariant; this is the central validation of matching-plane coupling. Also verify `outward = return + escape` for PEs at
the interface.

**Scope:** This quasistatic model does not solve a six-dimensional VDF, an outer flight-time queue, collisions,
magnetized return, outer-sheath transients, or volume plasma charge inside the BEACH region. The online MVP also does
not solve outer return of ambient outward populations. Use one-dimensional--3D kinetic coupling or full PIC for
validation or production when these effects control the result.
