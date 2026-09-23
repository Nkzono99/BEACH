title: Couple an outer sheath at a matching plane

Lang: [日本語](MatchingPlaneCoupling.md) | [English](MatchingPlaneCoupling.en.md)

# Couple an outer sheath at a matching plane

`surface_current_model.model="matching_plane_quasistatic"` makes the top of the BEACH box a matching plane with an
outer one-dimensional sheath. As the surface charge changes, the outer response can update the top potential, ambient
particle inflow, and the photoelectron (PE) return barrier.

This feature does not evolve the outer sheath in time. It is a boundary closure that reconciles the three-dimensional
trajectories inside BEACH with a reduced quasistatic outer response once per accepted batch. See
[Input parameters](Parameters.en.html#surface_current_model-external-sheath-closure) for the complete input contract and
[Output format reference](OutputReference.en.html#matching_plane_quasistatic) for exact output names.

This page covers ordinary selection, configuration, execution, and diagnosis. Use the
[matching-plane numerical and response-table reference](MatchingPlaneReference.en.html) to look up the response CSV,
`implicit_zero_mode`, and fixed-point equations.

## 1. Decide whether to use this model

Choose the model from the outer-plasma behavior that the case needs.

| Required representation | Model |
|---|---|
| Change the outer response from the evolving surface charge each batch | `matching_plane_quasistatic` |
| Hold currents and barriers from a stationary-sheath zero-current root fixed during the run | [`zhao_stationary`](ZhaoStationaryClosure.en.html) |
| Use no outer-sheath closure and model only the field and particles inside BEACH | `none` |

For matching-plane coupling, also choose how to obtain the outer response.

| Backend | Appropriate use | Main constraint |
|---|---|---|
| `response_backend="table"` (default) | Use an independently validated Zhao or 1-D PIC sweep as an auditable fixed snapshot | Queries are cheap; requires a response CSV and finite axis ranges |
| `response_backend="zhao_online"` | Use the built-in Zhao response without preparing a CSV | `implicit_zero_mode` searches within the selected branch; planar, collisionless, unmagnetized, and pays for nonlinear solves |

`examples/matching_plane_response_synthetic.csv` is only a table-path smoke-test fixture. Do not use it for a physical run.
See the [offline kinetic oracle](OuterKineticOracle.en.html) for the separate research workflow that evaluates a 1D1V
Vlasov response.

### With and without PEs

| Configuration | What BEACH represents |
|---|---|
| Without PEs | Ambient-electron and ion inflow plus the outer potential gauge; PE feedback, return, and escape are zero |
| With PEs | The above plus outward PE flux and mean normal energy, and return or escape at the outer barrier |

The outer-sheath connection remains active without PEs. Omit `photoelectron_species` and the `photo_raycast` species;
the response backend still determines the matching potential and ambient inflow.

### Retain the PE distribution

BEACH tracks return and absorption between the emitting surface and H as particle trajectories. The outer sheath is
supplied by the PEs that cross H outward after this transport. Substituting the surface `emit_current_density_a_m2`
directly for the H flux would also supply the exterior with particles that never reached H.

The default `photoelectron_closure="moment_matched_half_maxwellian"` maps the flux and mean normal energy at H
to a Maxwell distribution with the same two moments. To preserve the shape selected by inner transport, add the
following to an online Zhao configuration with PEs:

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
photoelectron_closure = "energy_spectrum"
photoelectron_spectrum_bins_per_decade = 32
# Also specify electron_species / ion_species / photoelectron_species.
```

This option records outward crossings at H before the outer reflection decision and derives exterior density,
return, and escape from the energy-resolved flux. The implicit mean-charge update uses the escape integral of
the same distribution. The PE source remains `photo_raycast`. A complete configuration is provided in
[`periodic2_matching_plane_pe_spectrum.toml`](../examples/periodic2_matching_plane_pe_spectrum.toml).

The existing five-input response table cannot encode a distribution, so this option cannot be used with `table`.
Compare the observed and response-input distributions in `matching_plane_spectrum_history.csv` and check convergence
by increasing bin and ray counts. This remains a planar 1D exterior approximation and does not solve the volume
space charge of PEs flying inside BEACH. See the [numerical reference](MatchingPlaneReference.en.html#photoelectron_closure)
for distribution, iteration, and restart details.

### Distinguish stationary Zhao from online Zhao

The similarly named models solve different constraints.

| | `zhao_stationary` | Matching-plane `zhao_online` |
|---|---|---|
| Constraint | Zero wall current, $J=0$ | Current $D_H/\epsilon_0$ prescribed as the interface field |
| Update | Once at run startup | At every fixed-point query |
| Without PEs | Type C | Not fixed to Type C; $D_H=0$ is the degenerate Type-B state |
| Current treatment | Fixed target for each species | Response fluxes and raw trajectory deposits |

Therefore, omitting PEs does not by itself make an online Zhao result Type C. With the default
`zhao_root_selection="require_unique"`, `zhao_branch="auto"` looks for one physical root compatible with the current
$D_H$ and stops when it cannot establish uniqueness.

## 2. Build a minimal configuration

Start from a `periodic2` case that is periodic in x and y and open in z. Put every mesh vertex strictly below
the z component of `domain.box_max`, set `sim.e0` and `sim.b0` to zero, and inject ambient electrons and ions from the z-high reservoir.

The minimal model selection for `response_backend="zhao_online"` with PEs is:

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.05, 0.0, 0.0]
```

The second component of `coupling_atol` is the PE mean-normal-energy tolerance in eV. The 0.05 eV above is an example;
select it through a convergence study in the ray and macro-particle counts. See
`examples/periodic2_matching_plane_zhao_online.toml` for a complete case.

To use a validated response table, replace the model section with:

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "table"
response_table_path = "matching_plane_response.csv"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.0, 0.0, 0.0]
```

`response_table_path` is resolved relative to the configuration file. The complete wiring example is
`examples/periodic2_matching_plane_quasistatic.toml`.

An online case without PEs does not declare a PE role.

```toml
[periodic2]
lower_boundary_model = "e_bottom_zero"

[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
implicit_zero_mode = true
```

This path needs neither `response_table_path` nor a prebuilt `matching_query.csv`. See
`examples/periodic2_matching_plane_zhao_implicit.toml` for the complete CSV-free case. A query grid and
`beach-zhao-response` are needed only when creating a separate table snapshot. Table generation has no accepted-endpoint
history and therefore rejects `continuation`.

Implicit integration alone does not choose a Zhao branch. By default, it stops if `auto` cannot certify a unique
physical solution at the current seed. With strong PE, scan `a`, `b`, and `c` separately and then select the validated
branch explicitly.

When Type B is the only candidate, BEACH searches the boundary potential directly to reduce missed narrow displacement
intervals. Multiple physical endpoints can remain even when the search follows the spectrum bins. See the
[numerical reference](MatchingPlaneReference.en.html#implicit_zero_mode) for bootstrap uniqueness, nearest-root continuation,
and validity near infinity.

When omitting PEs from a table case, also make the table's PE-flux and PE-energy axes zero-valued singletons. Check the
complete species, boundary, and `periodic2` requirements in
[Input parameters](Parameters.en.html#matching-plane-quasistatic-closure).

First, run the bundled PE-enabled `response_backend="zhao_online"` example without creating another configuration.
It is a four-batch wiring check, not a physically validated research case. From the repository root, run:

```bash
beachx lint examples/periodic2_matching_plane_zhao_online.toml
beach examples/periodic2_matching_plane_zhao_online.toml
beach-inspect outputs/periodic2_matching_plane_zhao_online
```

After a successful run, `outputs/periodic2_matching_plane_zhao_online/` contains at least `summary.txt`, `charges.csv`,
and `matching_plane_history.csv`.

Check `batches=4` and `matching_plane_state_valid=T` in `summary.txt`. They establish four accepted batches and a saved
outer state. Confirm fixed-point convergence from the residuals under [output checks](#4-decide-whether-the-result-succeeded):
a finite state that reaches the iteration limit can also be accepted with a warning.

`beachx lint` checks TOML and known parameter combinations, but it does not read the response CSV. With
`response_backend="table"`, `beach` checks the table header, Cartesian grid, and matching-plane height at startup.

### Choose an online Zhao root policy

After checking the basic example, change the selection policy only when investigating the effect of multiple roots
on your research case.

| `zhao_root_selection` | Selection | What to check |
| --- | --- | --- |
| `require_unique` (default) | Require a unique physical root at each query | Stop if uniqueness cannot be certified |
| `minimum_energy` | Choose the detected candidate with the lowest full-sheath potential energy | A root switch can make the response discontinuous and eliminate the implicit-update endpoint |
| `continuation` | Track from the last accepted root | Requires `implicit_zero_mode=true`; supports `auto` / `a` / `b` / `c` |

To select `continuation`, change these keys in the existing `[surface_current_model]` table:

```toml
response_backend = "zhao_online"
zhao_branch = "auto"
implicit_zero_mode = true
zhao_root_selection = "continuation"
```

The first solve uses multiple guesses derived from the input temperature, flux, and field scales and requires one
physical root. Subsequent solves track from the previous accepted root; when local search cannot reacquire it, full
multistart selects the unique nearest root. Missing roots, ambiguity, or numerical failure stop the solve. An explicit
branch is preserved; `auto` considers certified A / B / C candidates. This does not guarantee retention of the same
physical family or passage through a fold.
Feedback iterations within one batch use the latest valid root; the seed for the next batch is committed only when
the batch is accepted.

No policy proves that a finite set of initial guesses found every mathematical root or establishes time-dependent
stability. The energy definition, distance and ambiguity rules, probe subdivision, and restart behavior are documented
in the [`zhao_root_selection` reference](MatchingPlaneReference.en.html#zhao_root_selection).

## 3. What happens in one accepted batch

1. **Measure the interface state.** Compute the mean displacement $D_H$ immediately below the matching plane from the
   current surface charge and lower boundary.
2. **Obtain the outer response.** From $D_H$ and the outward feedback, the backend returns the matching potential
   $\Phi_H$, electron and ion inward fluxes and access potentials, and the PE barrier. Online implicit mode re-solves
   the backward-Euler endpoint inside each feedback iteration.
3. **Track the same batch.** Use $\Phi_H$ as the `periodic2` zero-mode gauge, track particles, and measure outward moments
   and PE return and escape.
4. **Check the fixed point.** If the measurements disagree with the assumed feedback, relax the feedback and replay.
   Every trial starts from the same batch-start RNG state and macro-particle residuals, and no trial changes state before
   acceptance. If the finite final trial reaches the iteration limit, BEACH records its residual and accepts it with a warning.
5. **Commit once.** Only a converged trial or a finite final trial accepted with a warning updates surface charge, RNG,
   ledger, history, and outer state. The `continuation` root seed belongs to this accepted state as well. If the adaptive
   $k\ne0$ condition rejects the trial, BEACH halves the batch duration and rolls the outer state back.

This replay tests the consistency of the outer response against one particle map instead of comparing different Monte Carlo draws.
Photoelectron emission is sampled before ambient particles, so response-dependent ambient counts do not change the emission
positions and velocities within a batch. This sampling-order correction changes the random sequence relative to earlier
versions even with the same `rng_seed`. Species identifiers and particle assembly order for tracking are unchanged.

### The 0 V reservoir and PE return

The outer model fixes the upstream plasma potential at 0 V. `matching_potential_v` and all three access or barrier
potentials use this same gauge; adding an arbitrary constant only to those potential columns is not equivalent.

When a PE crosses z-high outward, BEACH compares its local potential and normal kinetic energy with the outer PE barrier.
A PE that cannot cross is reflected specularly at z-high and counted as return; one that can cross is counted as escape.
PE return is therefore included, but the external turning-point distance and flight time are reduced to an immediate
boundary reflection.

Outward ambient electrons and ions are not reflected locally. A table can return a total inward flux that includes
populations returned by its outer model. Online Zhao v1 treats ambient-outward feedback as transparent and should not be
used when that outer return controls the result.

For online Zhao, branch and barrier have the following relation.

| Branch | $\Phi_H$ | Electron access / PE barrier |
|---|---:|---:|
| Type A | $\Phi_H>\phi_m$ (may be negative) | $\phi_m<\min(0,\Phi_H)$ |
| Type B | Positive (zero at $D_H=0$) | 0 V |
| Type C | Negative | 0 V |

`surface_current_model_zhao_branch=auto` in `summary.txt` records the selection policy, not the branch chosen by each
query. For an accepted online response, infer the branch from $\Phi_H$ and the access or barrier potential in this table.

$H$ fixes the interface, zero-mode gauge, and PE-moment measurement plane. Because online Zhao is planar and
translationally symmetric, its Sagdeev equation does not use the absolute coordinate of $H$ as a distance parameter and
does not solve a one-dimensional wall-to-$H$ profile.

## 4. Decide whether the result succeeded

With `output.history_stride>0`, BEACH writes `matching_plane_history.csv`. Check these items first.

| Question | Output | Acceptance check |
|---|---|---|
| Is there an accepted state? | `matching_plane_state_valid` | `T` in `summary.txt` |
| Did the fixed point converge? | `matching_plane_residual` | At or below `surface_current_model_coupling_rtol` |
| Is there iteration margin? | `matching_plane_iterations` | Not pinned to the limit; stable when controls change |
| Does PE classification close? | outward / return / escape flux | $\Gamma_{pe}^{out}\simeq\Gamma_{pe}^{return}+\Gamma_{pe}^{escape}$ |
| Is the table provenance identifiable? | response table path / configuration | Matches the production data |
| Are potential and charging stable? | $D_H$, $\Phi_H$, mesh charge / potential | Within tolerance as batch width, particle count, and mesh vary |

See [Output format reference](OutputReference.en.html#matching_plane_quasistatic) for all 17 accepted-state columns,
summary receipts, and the exact time convention.

### Diagnose warnings and stopped runs

| Symptom | Likely cause | Next action |
|---|---|---|
| Response preflight failure | Path, header, $H$, or Cartesian-grid mismatch | Compare the [CSV contract](MatchingPlaneReference.en.html#table-backend-response-csv-v1) with the z component of `domain.box_max` |
| Table query out of range | The active-axis sweep does not cover the transient | Do not extrapolate; regenerate the table over a physically validated range |
| Fixed point reaches the iteration limit and continues with a warning | Particle noise, strong feedback, or overly tight tolerances | Check frequency and residuals in history; adjust ray or macro count, relaxation, or tolerance if needed |
| Online Zhao has no or ambiguous physical solution | Incompatible $D_H$ and branch, multiple roots, or numerical failure | Scan `a`, `b`, and `c` separately; use `minimum_energy` only after validation |
| `continuation` stops | The bootstrap root is not unique, full multistart detects no root, root search or profile certification fails numerically, or nearest-root distances are numerically indistinguishable | Map solvability around the accepted state and reduce `batch_duration`; do not treat this as a fixed-point tolerance miss |
| Table implicit root is not bracketed | The backward-Euler endpoint is absent from the table | Revisit the $D_H$ range under the [`implicit_zero_mode` contract](MatchingPlaneReference.en.html#implicit_zero_mode) or reduce `batch_duration` |
| Online implicit root is not bracketed | The Zhao branch ends, or the finite potential/displacement search detects no sign change | Check the branch and initial charge; reduce `batch_duration` if needed. A missed root alone does not prove physical nonexistence |
| Soft-discard fraction limit or charge warning is reached | Unresolved periodic events are accumulating | Follow the [soft-discard stop conditions](ParticleEvents.en.html#advance-the-time-remaining-after-a-boundary-crossing) and inspect per-batch bursts, cumulative fraction, and absolute charge |

A completed run may contain batches accepted with warnings. Use history residuals to establish fixed-point convergence,
then check physical validity of the outer sheath, dependence on matching-plane height, and Monte Carlo convergence.

## 5. Accepted configuration and model limits

Matching-plane coupling is restricted to the following configuration to prevent double-counting a mean field or particle channel.

| Item | Requirement |
|---|---|
| Box / field | x/y periodic, z open, `field_boundary.mode="periodic2"`, and `sim.e0=sim.b0=[0,0,0]` |
| Periodic split | `cached_kneq0` or `panel_spectral_reference`, `exclude_k0`, and `symmetric_vacuum` or `e_bottom_zero` |
| Reservoir / open faces | `[reservoir].inflow_model="source_vdf"` and `ordinary_open_model="escape"` |
| Ambient species | Only electron and ion roles; `volume_seed`, `npcls_per_step=0`, and z-high reservoir inflow |
| PE species | Optional; negative `photo_raycast` from z-high with opposite charge deposited at emission |
| Surface closure | `explicit` for every role; no manual `fixed_current` target or `neutral_return` |
| Event policy | `abort`, or [`soft_discard` with a fraction limit, count grace, and absolute-charge warning threshold](ParticleEvents.en.html#advance-the-time-remaining-after-a-boundary-crossing) |

Do not specify `reference_area_m2` or stationary-Zhao source keys. The area comes from the domain x-y area, $H$ from
the z component of `domain.box_max`, and the update interval from one accepted batch. See
[Input parameters](Parameters.en.html#matching-plane-quasistatic-closure) for online Zhao's complete charge,
temperature, density, and drift restrictions.

This model does not solve:

- an outer six-dimensional VDF, particle inventory, flight time, or delayed-return queue;
- collisions, magnetized return, or outer-sheath transients;
- volume plasma charge inside the BEACH region; or
- return of outward ambient populations in online Zhao v1.

Online Zhao's default moment closure reduces PEs to a two-moment half-Maxwellian and does not retain the tail shape.
`energy_spectrum` avoids this reduction by using the measured normal-energy distribution, constant within each bin.
Ambient density uses the same orbit mapping as the injected VDF, but upstream neutrality, zero field, and a real
profile remain necessary. For the Type-A / C restriction with positive inward electron drift, see
[density and upstream conditions](MatchingPlaneReference.en.html#consistency-between-the-ambient-vdf-and-sheath-density).

The `auto` multiple-root check compares roots found by a finite multistart set; it is not mathematical root isolation.
Validate branches by scanning explicit `a`, `b`, and `c` selections.

`require_unique` and `minimum_energy` are stateless between queries. `continuation` retains the previous accepted
root as the next search seed. No policy silently switches an explicit branch or backend when a query cannot be solved.

When these effects control the result, validate against an independent one-dimensional--three-dimensional kinetic
coupling or full PIC calculation. The [numerical and response-table reference](MatchingPlaneReference.en.html#validate-convergence-and-applicability)
lists checks that vary the table grid, fixed-point tolerances, and matching-plane height.
