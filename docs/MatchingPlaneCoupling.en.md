title: Use Quasistatic Matching-Plane Coupling

Lang: [English](MatchingPlaneCoupling.en.md) | [日本語](MatchingPlaneCoupling.md)

# Use Quasistatic Matching-Plane Coupling

`surface_current_model.model="matching_plane_quasistatic"` treats the top of the BEACH box,
the z component of `domain.box_max`, as the interface to an outer one-dimensional sheath. BEACH retains every surface charge and the
`k=0` and `k!=0` fields below that interface. A precomputed outer-model response table supplies the matching potential,
ambient inward fluxes, and access or return barriers. Particle trajectories return outward fluxes to a fixed-point
iteration in every batch, and BEACH commits only a converged trial to the surface charge.

This feature does not reevaluate the Zhao stationary wall potential during a run. Build the response table from an
independent Zhao or one-dimensional PIC parameter sweep and validate that outer model over the intended range.
`examples/matching_plane_response_synthetic.csv` is synthetic data for exercising the execution path; it is not valid
input for a physical study.

## Minimal configuration

**Prerequisite:** Prepare a `periodic2` case that is periodic in x and y, open in z, and contains every mesh vertex
strictly below the top face. Set the external uniform electric and magnetic fields to zero. Define distinct electron, ion, and
photoelectron species.

**Action:** Add `[surface_current_model]`.

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_table_path = "matching_plane_response.csv"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

See `examples/periodic2_matching_plane_quasistatic.toml` for a complete synthetic smoke case.
BEACH resolves `response_table_path` relative to the directory containing `beach.toml`.

**Expected output:** Configuration validation accepts the model wiring, and startup validates the response-table syntax
and matching-plane height. `summary.txt` records the latest accepted coupling state together with the response-table
content fingerprint, three species roles, and iteration controls. When
`output.history_stride > 0`, `matching_plane_history.csv` tracks each selected accepted batch's response values, four
feedback moments, return and escape fluxes, iteration count, and normalized residual.

`beachx lint` does not read response-file contents.

**Interpretation:** Successful completion establishes valid table syntax, an in-range interpolation, and numerical
fixed-point convergence. It does not establish physical validity of the outer sheath, invariance to matching-plane
height, or Monte Carlo convergence.

## Configuration constraints

The model accepts only the following configuration so that no mean field or particle channel is counted twice:

- `field_boundary.mode="periodic2"`
- `periodic2.nonzero_mode_backend="cached_kneq0"` or `"panel_spectral_reference"`
- `periodic2.zero_mode_policy="exclude_k0"`
- `periodic2.lower_boundary_model="symmetric_vacuum"` or `"e_bottom_zero"`
- periodic x/y particle faces and open z-low/z-high faces
- `sim.e0 = [0, 0, 0]` and `sim.b0 = [0, 0, 0]`
- `[reservoir].inflow_model="source_vdf"` (the generic `infinity_barrier` is disabled)
- exactly three enabled species, one for each distinct role, each with `surface_charge_closure="explicit"`
- ambient electron and ion species with `source_mode="volume_seed"`, `npcls_per_step=0`, and only z-high
  `boundary_inflow="reservoir"`
- a negative-charge `photo_raycast` species
- an open outward z-high boundary for the photoelectrons
- `particle_boundary.ordinary_open_model="escape"`
- `sim.multiple_box_events_policy="abort"`

Do not combine this model with `reference_area_m2`, Zhao-only keys, or manual `fixed_current` targets. The interface area
is the x-y area of `domain`, its height is the z component of `domain.box_max`, and its update interval is one accepted batch. These are
derived rather than configurable so that mutually inconsistent values cannot be entered.

## Response CSV v1

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

Fluxes, mean energy, and response fluxes must be nonnegative. Feedback axes 2--5 must include zero so that the initial
state can be evaluated. Use one node for an input axis that the outer sweep intentionally does not model. A singleton
axis accepts any finite query and explicitly disables that feedback dimension; it is neither clamping nor implicit
extrapolation.

## Fixed-point iteration

Starting from the accepted feedback state $X^0$, each trial performs these steps:

1. Compute $D_H$ from the current surface charge and interpolate the table at $(D_H, X^m)$.
2. Apply the response $\Phi_H$ as the `periodic2` zero-mode gauge, the total electron/ion inward fluxes and access maps,
   and the outward PE barrier. Ambient outflow is not reflected locally; it is passed to the next outer response.
3. Replay the particle trial from the same batch-start RNG state and macro-particle residuals, measuring outward moments
   $X_{\mathrm{raw}}^{m+1}$.
4. Accept this trial as converged when the maximum difference between $X^m$ and $X_{\mathrm{raw}}^{m+1}$, normalized
   by each response-table feedback-axis span, is no greater than `coupling_rtol`.
5. Otherwise, with `coupling_relaxation` $\alpha$, calculate
   $X^{m+1}=(1-\alpha)X^m+\alpha X_{\mathrm{raw}}^{m+1}$ and replay the trial.

History response columns contain values evaluated at the accepted trial's $X^m$; feedback columns contain the measured
$X_{\mathrm{raw}}^{m+1}$ from that same trial. Their difference is bounded by the recorded residual, and no unexecuted
relaxation step is applied after convergence.

The iteration does not commit surface charge, RNG state, injection residuals, ledger entries, or histories. A rejected
adaptive batch-duration trial also restores the outer state to the batch start. If `coupling_max_iterations` is exhausted
or a query leaves the interpolation grid, BEACH fails instead of continuing with an unconverged value.

## Build a response table

Sweep all five input axes in the one-dimensional Zhao or PIC model using the same `z=H`, upstream distributions, and
sign conventions. Feed the table the PE flux and normal energy that actually cross the matching plane, not the emission
flux at the material wall. For a nonmonotonic potential, output the maximum barrier over the complete outer profile, not
only a local potential difference.

Cover the expected initial state, stationary state, and intermediate iteration states with margin. Keep axes that the
outer model does not depend on as singletons instead of refining all five dimensions uniformly. A wider numeric grid
does not extend the physical validity of the outer model. Store the table-generation code, upstream conditions, solver
version, and unit conversions with production data.

## Validate the result

**Numerical validation:** Vary these inputs independently and require the reported observables to remain within the
study's error budget:

- `coupling_rtol`, `coupling_relaxation`, and `coupling_max_iterations`
- response-grid resolution and range
- `sim.batch_duration` and macro-particle count
- mesh and periodic-cell resolution

**Coupling validation:** Regenerate the outer model for the same physical conditions while moving `H` within an overlap
region. Grain charge, gap potential, and PE escape fraction should be nearly invariant; this is the central validation
of matching-plane coupling. Also verify `outward = return + escape` for PEs at the interface.

**Scope:** This quasistatic model does not solve a six-dimensional VDF, an outer flight-time queue, collisions,
magnetized return, outer-sheath transients, or volume plasma charge inside the BEACH region. Use one-dimensional--3D
kinetic coupling or full PIC for validation or production when these effects control the result.
