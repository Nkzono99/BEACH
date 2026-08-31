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
| `response_backend="table"` (default) | Production use of an independently validated Zhao or 1-D PIC sweep | Requires a response CSV; supports [`implicit_zero_mode`](MatchingPlaneReference.en.html#implicit_zero_mode) |
| `response_backend="zhao_online"` | Examine the wiring and scope with the built-in Zhao response | Planar, collisionless, unmagnetized; ignores ambient-outward feedback |

`examples/matching_plane_response_synthetic.csv` is only a table-path smoke-test fixture. Do not use it for a physical run.

### With and without PEs

| Configuration | What BEACH represents |
|---|---|
| Without PEs | Ambient-electron and ion inflow plus the outer potential gauge; PE feedback, return, and escape are zero |
| With PEs | The above plus outward PE flux and mean normal energy, and return or escape at the outer barrier |

The outer-sheath connection remains active without PEs. Omit `photoelectron_species` and the `photo_raycast` species;
the response backend still determines the matching potential and ambient inflow.

### Distinguish stationary Zhao from online Zhao

The similarly named models solve different constraints.

| | `zhao_stationary` | Matching-plane `zhao_online` |
|---|---|---|
| Constraint | Zero wall current, $J=0$ | Current $D_H/\epsilon_0$ prescribed as the interface field |
| Update | Once at run startup | At every fixed-point query |
| Without PEs | Type C | Not fixed to Type C; $D_H=0$ is the degenerate Type-B state |
| Current treatment | Fixed target for each species | Response fluxes and raw trajectory deposits |

Therefore, omitting PEs does not by itself make an online Zhao result Type C. With `zhao_branch="auto"`, the solver
looks for one physical root compatible with the current $D_H$ and stops when it cannot establish uniqueness.

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
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
```

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

Check `batches=4` and `matching_plane_state_valid=T` in `summary.txt`. They establish four accepted batches and a solved
fixed point, not physical validity of the outer sheath or particle-count convergence.

`beachx lint` checks TOML and known parameter combinations, but it does not read the response CSV. With
`response_backend="table"`, `beach` checks the table header, Cartesian grid, and matching-plane height at startup.

## 3. What happens in one accepted batch

1. **Measure the interface state.** Compute the mean displacement $D_H$ immediately below the matching plane from the
   current surface charge and lower boundary.
2. **Obtain the outer response.** From $D_H$ and the outward feedback, the backend returns the matching potential
   $\Phi_H$, electron and ion inward fluxes and access potentials, and the PE barrier.
3. **Track the same batch.** Use $\Phi_H$ as the `periodic2` zero-mode gauge, track particles, and measure outward moments
   and PE return and escape.
4. **Check the fixed point.** If the measurements disagree with the assumed feedback, relax the feedback and replay.
   Every trial starts from the same batch-start RNG state and macro-particle residuals, and an unconverged trial changes no state.
5. **Commit once.** Only the converged trial updates surface charge, RNG, ledger, history, and outer state. If the adaptive
   $k\ne0$ condition rejects it, BEACH halves the batch duration and rolls the outer state back as well.

This replay tests the consistency of the outer response against one particle map instead of comparing different Monte Carlo draws.

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
| Type A | Positive | $\phi_m<0$ |
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
| Is the table provenance identifiable? | response path / content fingerprint | Matches the production data |
| Are potential and charging stable? | $D_H$, $\Phi_H$, mesh charge / potential | Within tolerance as batch width, particle count, and mesh vary |

See [Output format reference](OutputReference.en.html#matching_plane_quasistatic) for all 17 accepted-state columns,
summary receipts, and the exact time convention.

### Diagnose a stopped run

| Symptom | Likely cause | Next action |
|---|---|---|
| Response preflight failure | Path, header, $H$, or Cartesian-grid mismatch | Compare the [CSV contract](MatchingPlaneReference.en.html#table-backend-response-csv-v1) with the z component of `domain.box_max` |
| Table query out of range | The active-axis sweep does not cover the transient | Do not extrapolate; regenerate the table over a physically validated range |
| Fixed point reaches the iteration limit | Particle noise, strong feedback, or overly tight tolerances | Increase ray or macro counts; if residual decreases, adjust relaxation or the limit |
| Online Zhao has no or ambiguous physical solution | Incompatible $D_H$ and branch, multiple roots, or numerical failure | Scan `a`, `b`, and `c` separately and validate the physical branch |
| Implicit root is not bracketed | The backward-Euler endpoint is absent from the table | Revisit the $D_H$ range under the [`implicit_zero_mode` contract](MatchingPlaneReference.en.html#implicit_zero_mode) or reduce `batch_duration` |
| Soft-discard limit is reached | Unresolved periodic events exceed the error budget | Follow the [soft-discard stop conditions](ParticleEvents.en.html#advance-the-time-remaining-after-a-boundary-crossing) and inspect per-batch bursts, cumulative fraction, and absolute charge |

A completed run establishes backend evaluation and numerical fixed-point convergence. It does not establish physical
validity of the outer sheath, invariance to matching-plane height, or Monte Carlo convergence.

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
| Event policy | `abort`, or [`soft_discard` bounded by fraction, count grace, and absolute charge](ParticleEvents.en.html#advance-the-time-remaining-after-a-boundary-crossing) |

Do not specify `reference_area_m2` or stationary-Zhao source keys. The area comes from the domain x-y area, $H$ from
the z component of `domain.box_max`, and the update interval from one accepted batch. See
[Input parameters](Parameters.en.html#matching-plane-quasistatic-closure) for online Zhao's complete charge,
temperature, density, and drift restrictions.

This model does not solve:

- an outer six-dimensional VDF, particle inventory, flight time, or delayed-return queue;
- collisions, magnetized return, or outer-sheath transients;
- volume plasma charge inside the BEACH region; or
- return of outward ambient populations in online Zhao v1.

Online Zhao reduces the PE flux and mean normal energy to a half-Maxwellian that reproduces those two moments. It does
not retain the high-energy tail.

The `auto` multiple-root check compares roots found by a finite multistart set; it is not mathematical root isolation.
Validate branches by scanning explicit `a`, `b`, and `c` selections.

Every query is stateless and does not continue from the previous root. If it cannot solve a query, it does not silently
switch an explicit branch or backend.

When these effects control the result, validate against an independent one-dimensional--three-dimensional kinetic
coupling or full PIC calculation. The [numerical and response-table reference](MatchingPlaneReference.en.html#validate-convergence-and-applicability)
lists checks that vary the table grid, fixed-point tolerances, and matching-plane height.
