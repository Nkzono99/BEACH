title: BEACH Algorithm Guide

Lang: [English](Algorithms.en.md) | [日本語](Algorithms.md)

# BEACH Algorithm Guide

This document describes the numerical algorithms and execution order used by the current Fortran implementation of BEACH.
BEACH is not a grid PIC code. It is a surface charging simulator that couples a BEM-like Coulomb field evaluation, with charges on triangular boundary elements as sources, to test-particle tracking in batches.

Implementation links point to the current file and the main symbol. If line numbers drift in later edits, prefer the linked file and symbol name.

---

## Table of Contents

### Part I: Overview
1. [BEACH computation model](#1-beach-computation-model)
2. [Execution entry point and initialization](#2-execution-entry-point-and-initialization)
3. [Batch loop](#3-batch-loop)

### Part II: Field Calculation
4. [Coulomb field from boundary-element charge](#4-coulomb-field-from-boundary-element-charge)
5. [Switching between direct / treecode / FMM](#5-switching-between-direct--treecode--fmm)
6. [periodic2 field boundary](#6-periodic2-field-boundary)

### Part III: Particle Calculation
7. [Particle generation and injection state](#7-particle-generation-and-injection-state)
8. [Boris pusher](#8-boris-pusher)
9. [Collision detection](#9-collision-detection)
10. [Box boundary conditions](#10-box-boundary-conditions)

### Part IV: Surface Charge and Output
11. [Charge deposition and surface models](#11-charge-deposition-and-surface-models)
12. [Statistics, history, and restart](#12-statistics-history-and-restart)
13. [Parallelization and performance profiling](#13-parallelization-and-performance-profiling)

### Part V: FMM and Batch Stability
14. [Coulomb FMM core details](#14-coulomb-fmm-core-details)
15. [`batch_duration` stability and steady value](#15-batch_duration-stability-and-steady-value)

---

## 1. BEACH computation model

**Source**:
[`bem_simulator`](../src/runtime/simulator/bem_simulator.f90),
[`bem_simulator_loop`](../src/runtime/simulator/bem_simulator_loop.f90),
[`bem_field_solver`](../src/physics/field_solver/bem_field_solver.f90),
[`bem_injection`](../src/particles/bem_injection.f90)

The main state in BEACH is the charge `q_elem(i)` on triangular mesh elements plus the particles generated for each batch.
The mesh geometry is fixed. Particles are absorbed when they hit a surface. The absorbed particle charge is accumulated on the hit element and affects the field calculation in the next batch.

### 1.1 State variables

The main data types are defined in [`bem_types`](../src/core/bem_types.f90).

| Type | Main contents | Role |
| --- | --- | --- |
| `sim_config` | `dt`, `batch_count`, `max_step`, `field_solver`, `field_bc_mode`, `box_min/max`, `bc_low/high` | Time advance, field calculation, boundary conditions |
| `mesh_type` | `v0/v1/v2`, `centers`, `normals`, `bb_min/max`, `q_elem`, `elem_surface_model` | Triangular boundary elements and charge |
| `particles_soa` | `x`, `v`, `q`, `m`, `w`, `species_id`, `alive` | SoA representation of particles in a batch |
| `injection_state` | `macro_residual(:)` | Carries fractional `reservoir_face` macro-particle counts across batches and restarts |
| `sim_stats` | `processed_particles`, `absorbed`, `escaped`, `survived_max_step`, `batches`, `last_rel_change` | Run statistics |

### 1.2 Physical approximation

Each triangular element contributes to the field as a point charge `q_i` at its centroid `c_i`.
The current implementation uses a centroid point-charge approximation, not an exact BEM integration over a continuous surface charge distribution.

Softened Coulomb kernel:

$$
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\lVert\mathbf{r}\rVert^2 + \epsilon^2}}
$$

Potential:

$$
\phi(\mathbf{x}) =
k_\mathrm{c}
\sum_j q_j G_\epsilon(\mathbf{x} - \mathbf{c}_j)
$$

Electric field:

$$
\mathbf{E}(\mathbf{x}) =
k_\mathrm{c}
\sum_j q_j
\frac{\mathbf{x} - \mathbf{c}_j}
{\left(\lVert\mathbf{x} - \mathbf{c}_j\rVert^2 + \epsilon^2\right)^{3/2}}
$$

Here `k_coulomb` is the Coulomb constant from [`bem_constants`](../src/core/bem_constants.f90).
The field applied to particles is the boundary-element field plus the uniform external field `sim.e0`.

### 1.3 Execution unit

BEACH advances batches up to `sim.batch_count`. Each batch means:

1. Generate particles for each species from the current surface charge state.
2. Refresh the field solver with the current `q_elem`.
3. Track each particle for at most `sim.max_step` steps.
4. Accumulate charge from collided particles into element-wise deltas `dq`.
5. Add `dq` to `q_elem` and apply surface-model relaxation if needed.
6. Update statistics and history.

`sim.tol_rel` is an output and monitoring value. In the current Fortran implementation it is not an early-stop condition.

---

## 2. Execution entry point and initialization

**Source**:
[`app/main.f90`](../app/main.f90),
[`bem_app_config_runtime`](../src/config/bem_app_config_runtime.f90),
[`bem_mesh`](../src/mesh/bem_mesh.f90),
[`bem_restart`](../src/runtime/bem_restart.f90)

### 2.1 Main program order

`app/main.f90` is the CLI entry point. Its high-level order is:

| Order | Processing | Main implementation |
| --- | --- | --- |
| 1 | Handle CLI options that finish before config loading, such as `--version` | `handle_early_cli` |
| 2 | Initialize MPI and the performance profile | `mpi_initialize`, `perf_configure_from_env` |
| 3 | Load config, build mesh, read restart or initialize RNG seed | `load_or_init_run_state` |
| 4 | Open history CSV writers | `open_history_writer`, `open_potential_history_writer` |
| 5 | Run the batch simulation | `run_absorption_insulator` |
| 6 | Write the summary and final CSV files | `print_run_summary`, `write_result_files` |
| 7 | Write RNG state and macro residuals as checkpoints | `write_rng_state_file`, `write_macro_residuals_file` |
| 8 | Write performance profile and shut down MPI | `perf_write_outputs`, `mpi_shutdown` |

If a config path is passed explicitly, BEACH uses it. Otherwise it uses `beach.toml` in the current directory.
If no config file exists, BEACH runs with `default_app_config` defaults.

### 2.2 Mesh construction

The mesh is built from templates or an OBJ file according to `mesh.mode`. For template meshes,
[`build_template_mesh`](../src/config/bem_app_config_runtime.f90) does the following:

1. Iterate over `mesh.templates`.
2. Skip templates with `enabled=false`.
3. Dispatch by `kind` to `make_plane`, `make_box`, `make_cylinder`, `make_sphere`, and related builders.
4. Assign a `mesh_id` per template.
5. Expand `surface_model` and `epsilon_r` to element arrays.
6. Concatenate all template triangle arrays and pass them to `init_mesh`.

`init_mesh` precomputes:

- element centroids `centers(:, i)`
- element normals `normals(:, i)`
- element AABBs `bb_min/max(:, i)`
- representative length `h_elem(i) = sqrt(area_i)`
- initial charge `q_elem(i)`
- collision grid

### 2.3 periodic2 collision mesh

When `sim.field_bc_mode="periodic2"`, `prepare_periodic2_collision_mesh` translates each triangle to a canonical position along the two periodic axes.
This is separate from the periodic image sum used in the field calculation. It stabilizes the primitive-cell mesh for collision tests.
Element indices remain those of the base element, so a hit on a periodic image is still deposited on the base element.

### 2.4 Restart

When `output.resume=true`, `load_restart_checkpoint` reads:

- `summary.txt`: completed batch count and statistics
- `charges.csv`: element charge
- `rng_state.txt` or `rng_state_rankNNNNN.txt`: RNG state
- `macro_residuals.csv` or rank-local residual files: reservoir fractional particle counts

`sim.batch_count` is the cumulative target batch count. If the checkpoint has `batches=100` and `batch_count=150`, the resumed run executes only 50 batches.

---

## 3. Batch loop

**Source**:
[`run_absorption_insulator`](../src/runtime/simulator/bem_simulator_loop.f90),
[`prepare_batch_state`](../src/runtime/simulator/bem_simulator_loop.f90),
[`process_particle_batch`](../src/runtime/simulator/bem_simulator_loop.f90),
[`commit_batch_charge`](../src/runtime/simulator/bem_simulator_loop.f90)

### 3.1 Loop skeleton

When `run_absorption_insulator` receives `initial_stats`, it resumes from `initial_stats%batches`.

```text
final_batch_idx = sim.batch_count
batch_count_this_run = final_batch_idx - stats.batches

field_solver.init(mesh, sim)

for local_batch_idx = 1..batch_count_this_run:
    prepare_batch_state(...)
    field_solver.refresh(mesh)
    process_particle_batch(...)
    commit_batch_charge(...)
    count_batch_outcomes(...)
    MPI allreduce statistics
    accumulate_batch_stats(...)
    write charge/potential history when stride matches
```

### 3.2 Batch state

`prepare_batch_state` prepares:

- current batch number `batch_idx = stats%batches + 1`
- particles from `init_particle_batch_from_config`
- `dq_thread(nelem, nth)`: charge delta per OpenMP thread
- `photo_emission_dq(nelem)`: source-side opposite charge from `photo_raycast`
- `escaped_boundary_flag(:)` and `absorbed_flag(:)`

Splitting `dq_thread` by thread lets collision deposition be accumulated without atomics.

### 3.3 Particle processing

`process_particle_batch` advances each particle for at most `sim.max_step` time steps.

| Order | Processing |
| --- | --- |
| 1 | Read current position `x0` and velocity `v0` |
| 2 | Evaluate the electric field from boundary-element charge with `field_solver%eval_e(mesh, x0, e)` |
| 3 | Add the uniform external electric field `sim.e0` |
| 4 | Compute candidate position `x1` and velocity `v1` with `boris_push` |
| 5 | Test segment collision with `find_first_hit(mesh, x0, x1, hit, sim=sim)` |
| 6 | If there is a hit, add `q * w` to `dq_thread(elem, tid)` for the hit element and absorb the particle |
| 7 | If there is no hit, apply open / reflect / periodic with `apply_box_boundary` |
| 8 | If the particle remains alive, update `x` and `v` and continue to the next step |

If `BEACH_WARN_LONG_PARTICLE_STEPS` is set to a positive integer, BEACH prints diagnostics for long-lived particles at that step interval.

### 3.4 Charge commit

`commit_batch_charge` merges thread-local deltas and photo-emission deltas.

$$
\Delta q_i =
\sum_{\mathrm{thread}} \Delta q_{i,\mathrm{thread}}
{}+ \Delta q_{i,\mathrm{photo}}
$$

In MPI runs, `mpi_allreduce_sum_real_dp_array` sums `dq` across ranks.
Then BEACH applies:

$$
q_i^{new} = q_i^{old} + \Delta q_i
$$

and runs surface-model relaxation. The monitoring value is computed from the actually changed charge:

$$
\operatorname{rel\_change} =
\frac{\lVert\mathbf{q}^{\mathrm{new}} - \mathbf{q}^{\mathrm{old}}\rVert_2}
{\max\left(\lVert\mathbf{q}^{\mathrm{new}}\rVert_2, q_\mathrm{floor}\right)}
$$

This value is stored in `stats%last_rel_change` and written to history output.

---

## 4. Coulomb field from boundary-element charge

**Source**:
[`bem_field_solver`](../src/physics/field_solver/bem_field_solver.f90),
[`bem_field_solver_config`](../src/physics/field_solver/bem_field_solver_config.f90),
[`bem_field_solver_eval`](../src/physics/field_solver/bem_field_solver_eval.f90)

### 4.1 Direct evaluation

In direct mode, every element contributes directly to the field at evaluation point `r`.

$$
\mathbf{E}(\mathbf{r}) =
k_c \sum_{i=1}^{N}
q_i
\frac{\mathbf{r} - \mathbf{c}_i}
{\left(\lVert\mathbf{r} - \mathbf{c}_i\rVert^2 + \epsilon^2\right)^{3/2}}
$$

The cost is `O(nelem)` per evaluation point.
As the particle step count and particle count grow, this dominates, so large-element cases use treecode or FMM.

### 4.2 Length normalization

`sim.field_normalization` selects the internal length scale `L0`.

| Value | `L0` |
| --- | --- |
| `si` | `1 m` |
| `length` | `sim.field_length_scale` |
| `box` | `max(box_max - box_min)` |
| `mesh` | maximum width of the mesh bounding box |

Internally the solver evaluates with:

$$
\mathbf{x}' = \frac{\mathbf{x} - \mathbf{x}_0}{L_0},
\quad
\epsilon' = \frac{\epsilon}{L_0}
$$

The field is converted back to SI by multiplying by `k_c / L0^2`; the potential is converted by `k_c / L0`.
Input settings and output CSV files remain in SI units.

---

## 5. Switching between direct / treecode / FMM

**Source**:
[`init_field_solver`](../src/physics/field_solver/bem_field_solver_config.f90),
[`refresh_field_solver`](../src/physics/field_solver/bem_field_solver_tree.f90),
[`eval_e_field_solver`](../src/physics/field_solver/bem_field_solver_eval.f90),
[14. Coulomb FMM core details](#14-coulomb-fmm-core-details)

### 5.1 Mode selection

`sim.field_solver` accepts:

| Value | Behavior |
| --- | --- |
| `direct` | Always use direct sum |
| `treecode` | Use octree plus monopole approximation |
| `fmm` | Use the Coulomb FMM core |
| `auto` | Use treecode if `nelem >= tree_min_nelem`, otherwise direct |

In the current implementation, `periodic2` requires `field_solver="fmm"`.

### 5.2 Treecode

The treecode partitions element centroids into an octree.

1. Put all element indices in `elem_order`.
2. Compute the AABB of element centroids in the node.
3. Stop at a leaf if the element count is `leaf_max` or smaller, or if the node cannot be split.
4. Otherwise classify by the node center into eight octants and build child nodes recursively.

`refresh_field_solver` recomputes node monopoles bottom-up.

$$
Q_n = \sum_{i \in n} q_i
$$

$$
\mathbf{c}_{Q,n} =
\begin{cases}
Q_n^{-1}\sum_{i \in n} q_i \mathbf{c}_i, & |Q_n| > 0, \\
\mathbf{c}_{\mathrm{node}}, & Q_n \approx 0
\end{cases}
$$

During evaluation, a node is accepted as far field from its radius `R` and distance `d` to the target point.
Accepted nodes are evaluated as monopoles. Rejected nodes are descended, and leaves use direct sum.
For nodes with strong charge cancellation, monopole error can be large; if `abs(Q) < charge_cancellation_tol * sum(abs(q_i))`, far-field acceptance is suppressed.

### 5.3 FMM

FMM mode calls a simulator-independent Coulomb FMM core.
The field-solver adapter handles:

1. Convert mesh centroids to source coordinates `src_pos(3, nelem)`.
2. Build the geometry-dependent plan with `build_plan(plan, src_pos, options)`.
3. Update the charge-dependent state with `update_state(plan, state, q_elem)`.
4. Call `eval_point(plan, state, r, e)` for each particle position.

For P2M / M2M / M2L / L2L / L2P details, see
[14. Coulomb FMM core details](#14-coulomb-fmm-core-details).

---

## 6. periodic2 field boundary

**Source**:
[`bem_field_solver_config`](../src/physics/field_solver/bem_field_solver_config.f90),
[`bem_coulomb_fmm_periodic`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90),
[`bem_collision`](../src/physics/bem_collision.f90)

`sim.field_bc_mode="periodic2"` treats exactly two of the three axes as periodic.
An axis is periodic when `bc_low(axis) == bc_high(axis) == periodic`.
The third axis is the open direction.

### 6.1 Validation

`periodic2` requires:

- `sim.use_box = true`
- exactly two periodic axes
- `box_max - box_min > 0` for each periodic axis
- `sim.field_solver = "fmm"`

`field_periodic_far_correction` accepts `auto`, `none`, and `m2l_root_oracle`.
The current implementation normalizes `auto` to `none` for compatibility.
`m2l_root_oracle` is a diagnostic far correction enabled only by explicit request.

### 6.2 Near images and far correction

`field_periodic_image_layers = N` enumerates near images along the two periodic axes as:

$$
(i, j) \in [-N, N]^2
$$

The FMM core combines primary-cell sources and image sources for near contributions.
With `m2l_root_oracle`, the build stage fits a root-local correction from the Ewald residual and injects it into the root local expansion at runtime.

### 6.3 Collision side

The collision-side `periodic2` logic is separate from the field-calculation FMM.
`find_first_hit_periodic2` computes the needed image-shift range from the particle segment and the canonical mesh AABB, then tests the shifted segment against the base mesh.
If multiple candidates have the same `t`, ties are broken deterministically by element index and image shift.

---

## 7. Particle generation and injection state

**Source**:
[`init_particle_batch_from_config`](../src/config/bem_app_config_runtime.f90),
[`bem_injection`](../src/particles/bem_injection.f90),
[`bem_sheath_runtime`](../src/physics/sheath/bem_sheath_runtime.f90)

### 7.1 Species processing

`init_particle_batch_from_config` traverses all species, determines the number generated on each rank, and interleaves them into the SoA particle array.
In MPI runs, `mpi_split_count` distributes the global count across ranks.

Supported `source_mode` values:

| source_mode | Count | Position | Velocity |
| --- | --- | --- | --- |
| `volume_seed` | `npcls_per_step` | Uniform random in `pos_low` to `pos_high` | Shifted Maxwellian |
| `reservoir_face` | Determined dynamically from flux and batch duration | Injection-face rectangle | Maxwellian weighted by incoming flux, or velocity grid |
| `photo_raycast` | Number of hits from `rays_per_batch` rays | First surface hit by each ray | Flux-weighted Maxwellian along the surface normal |

### 7.2 reservoir_face macro count

The incoming flux of a drifting Maxwellian is computed from the normal velocity component against the injection-face inward normal `n`:

$$
u_n = \mathbf{u} \cdot \mathbf{n}
$$

and the thermal speed:

$$
\sigma = \sqrt{\frac{k_B T}{m}}
$$

The implementation uses `flux_weighted_normal_tail(vmin, u_n, sigma)` and counts only particles whose normal velocity is at least `vmin_normal`.

For area `A`, batch duration `T_b`, and macro-particle weight `w`:

$$
N_\mathrm{phys} = \Gamma_\mathrm{in} A T_b
$$

$$
N_\mathrm{macro,expected} = \frac{N_\mathrm{phys}}{w}
$$

The fractional part is carried in `injection_state%macro_residual(species)`.

$$
B = r_\mathrm{old} + N_\mathrm{macro,expected}
$$

$$
N_\mathrm{macro} = \lfloor B \rfloor,\quad
r_\mathrm{new} = B - \lfloor B \rfloor
$$

In MPI runs, `batch_duration_scale = 1 / nrank`, so each rank generates its share of the global flux.

### 7.3 Reservoir sampling

For `reservoir_face`, positions are sampled uniformly from the injection-face rectangle and, if requested, shifted slightly away from the face by `position_jitter_dt=sim.dt`.
Velocity is generated by one of the following methods:

- Generate a shifted Maxwellian and resample only the normal component from a flux-weighted distribution.
- If `velocity_distribution="grid"`, read a CSV velocity grid. With `phase_space`, use `max(v_n,0) f(v)`; with `flux_weighted`, treat the input values as the incoming distribution.

With `reservoir_potential_model="infinity_barrier"`, the normal-velocity lower bound is corrected from the difference between the injection-face mean potential and `phi_infty`.

### 7.4 photo_raycast

`photo_raycast` launches rays from the injection face and emits particles from the first mesh element hit by each ray.

1. Sample the ray origin uniformly on the injection-face rectangle.
2. Normalize `ray_direction` and confirm that it points inward.
3. Extend the ray to the box boundary and run `find_first_hit` on that segment.
4. If it hits the mesh, use the element normal to construct emission position and velocity.
5. If it hits a box boundary and the boundary condition reflects or periodically wraps it, continue up to `raycast_max_bounce`.

The weight per hit is:

$$
w_\mathrm{hit} =
\frac{J_\perp A_\perp T_b}
{|q| N_\mathrm{ray,global}}
$$

where `A_perp = A * abs(dot(ray_direction, inward_normal))`.

If `deposit_opposite_charge_on_emit=true`, the source element receives:

$$
\Delta q_\mathrm{emit} = -q_\mathrm{particle} w_\mathrm{eff}
$$

as `photo_emission_dq`.

### 7.5 Photo escape closure

With `photo_escape_model="boltzmann_cutoff"`, BEACH uses the centroid potential excluding the contribution from the emitting element itself.

$$
\mathrm{barrier} =
\max(\phi_\mathrm{emit} - \phi_\infty, 0)
$$

$$
f_\mathrm{escape} =
\exp\left(
{}-\frac{|q|\,\mathrm{barrier}}{k_B T_\mathrm{PE}}
\right)
$$

The effective weight is:

$$
w_\mathrm{eff} = w_\mathrm{hit} f_\mathrm{escape}
$$

This is a reduced closure: returning photoelectrons are not tracked individually and are treated as immediate neutralization.

---

## 8. Boris pusher

**Source**:
[`bem_pusher`](../src/physics/bem_pusher.f90)

Particle motion is advanced with the Boris method using the uniform magnetic field `sim.b0` and the electric field `E` evaluated at the particle position.

Inputs:

- position `x`
- velocity `v`
- charge `q`
- mass `m`
- time step `dt`
- electric field `E`
- magnetic flux density `B`

Update equations:

$$
\mathbf{v}^- =
\mathbf{v}^n
{}+ \frac{q}{m}\mathbf{E}\frac{\Delta t}{2}
$$

$$
\mathbf{t} =
\frac{q}{m}\mathbf{B}\frac{\Delta t}{2}
,\quad
\mathbf{s} =
\frac{2\mathbf{t}}{1 + \lVert\mathbf{t}\rVert^2}
$$

$$
\mathbf{v}' =
\mathbf{v}^- + \mathbf{v}^- \times \mathbf{t}
$$

$$
\mathbf{v}^+ =
\mathbf{v}^- + \mathbf{v}' \times \mathbf{s}
$$

$$
\mathbf{v}^{n+1} =
\mathbf{v}^+ + \frac{q}{m}\mathbf{E}\frac{\Delta t}{2}
$$

$$
\mathbf{x}^{n+1} =
\mathbf{x}^{n} + \mathbf{v}^{n+1}\Delta t
$$

BEACH runs collision detection on the segment `x^n -> x^{n+1}`.
If there is a collision, the particle is absorbed and `x^{n+1}` is not saved back into particle state.

---

## 9. Collision detection

**Source**:
[`bem_collision`](../src/physics/bem_collision.f90),
[`bem_mesh`](../src/mesh/bem_mesh.f90)

### 9.1 Broad phase

`init_mesh` builds each element AABB and the collision grid.
Small meshes use linear search; larger meshes use a uniform grid plus 3D-DDA.

Collision grid:

1. Build the bounding box of all element AABBs.
2. Estimate cell width from `target_elems_per_cell`.
3. Register each element index, in CSR form, into cells overlapped by its AABB.

The particle segment `p0 -> p1` is first intersected with the grid AABB, and only traversed cells are enumerated by 3D-DDA.
Only elements registered in those cells are passed to the narrow phase.

### 9.2 Narrow phase: Moller-Trumbore

The segment

$$
\mathbf{p}(t) = \mathbf{p}_0 + t(\mathbf{p}_1 - \mathbf{p}_0), \quad 0 \le t \le 1
$$

and the triangle

$$
\mathbf{v}(u,v) =
\mathbf{v}_0 + u(\mathbf{v}_1-\mathbf{v}_0) + v(\mathbf{v}_2-\mathbf{v}_0)
$$

are tested with the Moller-Trumbore method.
The conditions are:

- the triangle is not degenerate
- the segment direction is not nearly parallel to the triangle plane
- `0 <= u <= 1`
- `0 <= v`
- `u + v <= 1`
- `0 <= t <= 1`

If multiple elements are hit, the one with the smallest `t` is selected.

### 9.3 periodic2 collision

In `periodic2`, the mesh itself contains only base elements.
`find_first_hit_periodic2` computes the required image-shift range from the segment and mesh AABB.

$$
n_\mathrm{min} =
\left\lceil
\frac{\operatorname{min}(p_0, p_1) - \operatorname{max}(\mathrm{mesh}) - \mathrm{tol}}{L}
\right\rceil
$$

$$
n_\mathrm{max} =
\left\lfloor
\frac{\operatorname{max}(p_0, p_1) - \operatorname{min}(\mathrm{mesh}) + \mathrm{tol}}{L}
\right\rfloor
$$

For each image, the segment is shifted by `-shift` and intersected against the base mesh.
The hit record stores both the physical image-coordinate hit position `hit%pos` and the primary-cell wrapped position `hit%pos_wrapped`.

---

## 10. Box boundary conditions

**Source**:
[`bem_boundary`](../src/physics/bem_boundary.f90)

If a particle does not hit the mesh and the candidate updated position leaves the simulation box, axis-wise boundary conditions are applied.

| Boundary condition | Processing |
| --- | --- |
| `open` | Remove the particle and count it as `escaped_boundary` |
| `reflect` | Mirror the position at the boundary plane and reverse the normal velocity component |
| `periodic` | Wrap to the opposite side |

`apply_box_boundary` checks the three axes in order.
If one step crosses multiple periodic lengths, `periodic` still returns the particle to the box using `modulo`.
For `reflect` and `periodic`, the result is clamped slightly inside the box, around `1e-12`, to avoid numerical instability from landing exactly on a boundary.

---

## 11. Charge deposition and surface models

**Source**:
[`commit_batch_charge`](../src/runtime/simulator/bem_simulator_loop.f90),
[`bem_surface_models`](../src/physics/bem_surface_models.f90)

### 11.1 Insulator accumulation

With the default `surface_model="insulator"`, absorbed particle charge is accumulated directly on the hit element.

When particle `p` hits element `i`:

$$
\Delta q_i \mathrel{+}= q_p w_p
$$

This accumulated charge is used as source charge in the next batch's field-solver refresh.

### 11.2 Photo emission bookkeeping

When `photo_raycast` uses `deposit_opposite_charge_on_emit=true`, the opposite-sign charge is added to the emitting element.
This is accumulated as `photo_emission_dq`, separately from charge deposited by later particle collisions, and merged at batch commit.

### 11.3 Floating conductor relaxation

`surface_model="conductor"` is available only when `field_bc_mode="free"`.
Conductor elements are grouped by `mesh_id` as floating conductor groups.
The goal is to equalize element potentials within each conductor object while conserving the total charge of that object.

Unknowns:

- conductor element charges `q_i`
- equipotential value `V_g` for each conductor group

For an element `i` in group `g(i)`:

$$
\sum_j A_{ij} q_j - V_{g(i)} =
{}-\phi_i^\mathrm{fixed}
$$

where

$$
A_{ij} =
\begin{cases}
1/\epsilon, & i=j,\ \epsilon>0, \\
2\sqrt{\pi}/h_i, & i=j,\ \epsilon=0, \\
1/\sqrt{\lVert\mathbf{c}_i-\mathbf{c}_j\rVert^2+\epsilon^2}, & i\ne j
\end{cases}
$$

`phi_fixed` is the potential from non-conductor charge and the uniform external electric field, divided by `k_coulomb`.

Each group adds a total-charge conservation constraint:

$$
\sum_{i \in g} q_i = Q_g^\mathrm{before}
$$

The resulting square linear system is solved by Gaussian elimination with partial pivoting, and the conductor elements' `q_elem` values are replaced.

### 11.4 Dielectric metadata

`surface_model="dielectric"` and `epsilon_r` are metadata in the current version.
Dielectric polarization and dielectric boundary conditions are not yet reflected in the field calculation.

---

## 12. Statistics, history, and restart

**Source**:
[`bem_simulator_stats`](../src/runtime/simulator/bem_simulator_stats.f90),
[`bem_simulator_io`](../src/runtime/simulator/bem_simulator_io.f90),
[`bem_output_writer`](../src/runtime/bem_output_writer.f90),
[`bem_restart`](../src/runtime/bem_restart.f90)

### 12.1 Batch outcomes

`count_batch_outcomes` counts local-rank particles in five categories.

| Index | Meaning |
| --- | --- |
| 1 | Total particles in the batch |
| 2 | Particles absorbed by the mesh |
| 3 | Particles counted as escaped |
| 4 | Particles that left through an open box boundary |
| 5 | Particles that survived to `max_step` |

Among particles without `absorbed_flag`, those removed by an open boundary are `escaped_boundary`; those still `alive` at the end are `survived_max_step`.

In MPI runs, the batch-count array is allreduced so particles on non-root ranks are included in statistics.

### 12.2 History output

When `history_stride > 0`, BEACH writes `charge_history.csv`.
The output condition is:

$$
(stats.batches - 1) \bmod history\_stride = 0
$$

so batch 1 is always written.

When `output.write_potential_history=true`, `potential_history.csv` is written at the same stride.
For potential history, the field solver is refreshed with the current `q_elem`, and centroid potentials are computed and written.

### 12.3 Final output

When `output.write_files=true`, the root rank writes the main final outputs:

- `summary.txt`
- `charges.csv`
- `mesh_potential.csv` when enabled
- `mesh_triangles.csv`
- `mesh_sources.csv`

All ranks save checkpoint RNG state and macro residuals.
In MPI runs, file names are rank-local.

### 12.4 Restart consistency

On restart, BEACH validates:

- checkpoint `mesh_nelem` matches the current mesh element count
- MPI world size matches the previous run
- statistics in `summary.txt` are finite and non-negative
- `charges.csv` element count and charge values are valid
- RNG state and macro residuals can be loaded

If a required checkpoint is missing, `output.resume=true` stops instead of falling back to a new run.

---

## 13. Parallelization and performance profiling

**Source**:
[`bem_mpi`](../src/core/bem_mpi.F90),
[`bem_performance_profile`](../src/runtime/bem_performance_profile.f90),
[`bem_simulator_loop`](../src/runtime/simulator/bem_simulator_loop.f90)

### 13.1 OpenMP

Particle tracking is parallelized over particle indices with OpenMP.

- `dq_thread(nelem, nth)` accumulates collision charge thread-locally.
- The schedule is `dynamic, 1` to reduce load imbalance from different particle lifetimes.
- Some OpenMP loops are also used in field-solver refresh and treecode node accumulation.

### 13.2 MPI

MPI parallelism splits particle generation and tracking by rank.
Each rank holds the mesh and `q_elem`; at batch commit, `dq` is allreduced so every rank has the same charge state.

Main allreduces:

- sum of `dq(nelem)`
- sum of batch outcome counts

Only the root rank writes human-readable final CSV and history files.
RNG state and macro residuals are saved per rank.

### 13.3 Performance profile

Set `BEACH_PROFILE=1` to write times for major regions to `performance_profile.csv`.
Main regions:

- `load_or_init`
- `field_solver_init`
- `prepare_batch`
- `field_refresh`
- `particle_batch`
- `commit_charge`
- `mpi_reduce`
- `stats_update`
- `history_write`
- `write_results`
- `write_checkpoint`

For MPI scaling evaluation, it is usually best to inspect the maximum time across ranks, `rank_max_s`.

---

## 14. Coulomb FMM core details

This section summarizes the specification and algorithms of the current Fortran Coulomb FMM core,
[`bem_coulomb_fmm_core` module page](../module/bem_coulomb_fmm_core.html),
and its split implementation files.

- Public API / boundary: `src/physics/field_solver/fmm/api/`
- Internal shared implementation: `src/physics/field_solver/fmm/internal/common/`
- Tree / plan implementation: `src/physics/field_solver/fmm/internal/tree/`
- State / eval implementation: `src/physics/field_solver/fmm/internal/runtime/`
- periodic2 implementation: `src/physics/field_solver/fmm/internal/periodic/`

The target is a simulator-independent internal API. It does not directly `use` `mesh_type` or `sim_config`.
On the BEACH side, the field-solver adapter calls this core.

### 1. Purpose

The FMM core returns Coulomb electric fields quickly at many evaluation points for a fixed source point set `src_pos(3,n)` and variable charges `src_q(n)`.

Current design goals:

- kernel is only 3D Coulomb
- source geometry and charge updates are separated
- only `free` and `periodic2` are supported
- near direct sum is also handled inside the core
- simulator code sees only array APIs

### 2. Public API

The core provides four main procedures:

```fortran
call build_plan(plan, src_pos, options)
call update_state(plan, state, src_q)
call eval_points(plan, state, target_pos, e)
call eval_point(plan, state, r, e)
```

Input and output meanings:

- `src_pos(3,n)`: source point coordinates, fixed after `build_plan`
- `src_q(n)`: source charges, updateable at each `update_state`
- `target_pos(3,m)` or `r(3)`: evaluation points
- `e(3,m)` or `e(3)`: electric field vectors

Notes:

- The returned field does not include `k_coulomb`; the BEACH adapter multiplies it at the end.
- `build_plan` is geometry-dependent processing, and `update_state` is charge-dependent processing.
- `eval_point(s)` assumes `plan` and `state` are ready.

#### 2.2 C ABI / Python integration

`src/physics/field_solver/bem_field_kernel_c.f90` exposes this Fortran API as an `iso_c_binding` opaque-handle API.
`make build-kernel` builds the shared library as `build/libbeach_field_kernel.so`.

Main C ABI:

```text
beach_kernel_create(handle)
beach_kernel_destroy(handle)
beach_kernel_build(handle, src_pos, options...)
beach_kernel_update_charges(handle, src_q)
beach_kernel_eval_e(handle, target_pos, e)
beach_kernel_eval_phi(handle, target_pos, phi)
beach_kernel_force_on_charges(handle, target_pos, target_q, origin, force, torque)
```

The Python side calls this ABI with `ctypes` through `beach.fortran_results.kernel.FieldKernel`.
`calc_object_forces_kernel` evaluates `sum(q_i E_not_self(r_i))` by zeroing the object's own source charge, avoiding self-force contamination while using the same field kernel, including `periodic2 + m2l_root_oracle`.
`Beach.scene()` / `BeachScene` temporarily apply rigid translations and rotations of objects on the Python side and pass the edited centroid array to the same ABI.
The rigid-transform helper path uses NumPy by default and can use an optional Numba backend, but field evaluation itself is done by the Fortran kernel.

#### 2.3 BEACH adapter usage

The BEACH field-solver adapter passes mesh element centroids to this core as `src_pos`.

- During initialization, it calls `update_state` immediately after `build_plan`.
- During later refreshes, normal operation assumes mesh geometry is unchanged, so the existing `plan` is reused and only `update_state` is called with updated `src_q`.
- `build_plan` and legacy tree metadata are synchronized again only when the plan is missing, the source count changes, or zero elements caused plan/state disposal.

### 3. Data structures

#### 3.1 `fmm_options_type`

Main internal options:

- `theta`: parameter for well-separated tests
- `leaf_max`: maximum source count in a source-octree leaf
- `order`: Cartesian expansion order
- `softening`: `epsilon` of the softened Coulomb kernel
- `use_periodic2`: enable two-periodic-axis mode
- `periodic_axes(2)`, `periodic_len(2)`: periodic axes and lengths
- `periodic_image_layers`: near image-sum layer count `N`
- `periodic_far_correction`: core values are `auto`, `none`, and `m2l_root_oracle`; with `periodic2`, `auto` is normalized to `none`, and `m2l_root_oracle` is enabled only when explicit
- `periodic_ewald_alpha`, `periodic_ewald_layers`: decomposition parameter and cutoff depth used by the build-time Ewald fit for `m2l_root_oracle`
- `target_box_min/max`: box used for a dual-target tree

The BEACH adapter currently uses `order = 4`, while the core itself accepts variable order.
For `periodic2`, `auto` is normalized to `none`; `m2l_root_oracle` explicitly enables far correction.

#### 3.2 `fmm_plan_type`

This is geometry-dependent immutable data:

- multi-index tables `alpha`, `deriv_alpha`
- source octree
- optional target tree
- source leaf list `source_leaf_nodes`
- target leaf list `leaf_nodes`
- near lists `near_start/near_nodes`
- far node lists `far_start/far_nodes`
- M2L pair cache `m2l_target_nodes/m2l_source_nodes`
- periodic image-shift arrays
- M2L derivative table `m2l_deriv`
- P2M basis table `source_p2m_basis`
- compressed translation tables for M2M/L2L

#### 3.3 `fmm_state_type`

This is charge-dependent data updated on each refresh:

- `src_q(n)`
- `multipole(ncoef, nnode)`
- `local(ncoef, n_target_nodes)`
- `multipole_active(nnode)`
- `local_active(n_target_nodes)`

`multipole` stores multipole coefficients per source-tree node, and `local` stores local expansion coefficients per target-tree node.
`*_active` flags are 0/1 flags used to skip zero nodes quickly.

### 4. Mathematical definitions

#### 4.1 Kernel

The current core uses the softened Coulomb kernel:

$$
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\lVert\mathbf{r}\rVert^2 + \epsilon^2}}
$$

$$
\phi(\mathbf{x}) = \sum_j q_j \, G_\epsilon(\mathbf{x} - \mathbf{x}_j)
$$

$$
\mathbf{E}(\mathbf{x}) = - \nabla \phi(\mathbf{x})
$$

Both the near direct sum and far-field expansions use the same $G_\epsilon$.

#### 4.2 Multi-index

The core uses a multi-index $\alpha = (\alpha_x, \alpha_y, \alpha_z)$.

$$
|\alpha| = \alpha_x + \alpha_y + \alpha_z
$$

$$
\alpha! = \alpha_x! \, \alpha_y! \, \alpha_z!
$$

$$
\mathbf{r}^\alpha = r_x^{\alpha_x} r_y^{\alpha_y} r_z^{\alpha_z}
$$

#### 4.3 P2M

For node center $c$, leaf-node multipole coefficients are:

$$
M_\alpha(c) = \sum_{j \in \text{leaf}} q_j
\frac{(\mathbf{x}_j - \mathbf{c})^\alpha}{\alpha!}
$$

#### 4.4 M2M

Child-node coefficients are translated to the parent center and accumulated.
With $\mathbf{d} = c_{\mathrm{child}} - c_{\mathrm{parent}}$:

$$
M_\beta(c_{\mathrm{parent}}) =
\sum_{\alpha \le \beta}
M_\alpha(c_{\mathrm{child}})
\frac{\mathbf{d}^{\beta-\alpha}}{(\beta-\alpha)!}
$$

The current implementation precomputes, during `build_plan`, the index for $\beta - \alpha$ and the value
$\mathbf{d}^{\beta-\alpha} / (\beta-\alpha)!$.

#### 4.5 M2L

For source-node center $c_s$ and target-node center $c_t$, let $R = c_t - c_s$.

Local expansion coefficients are updated as:

$$
L_\alpha(c_t) \mathrel{+}=
\sum_\beta (-1)^{|\beta|}
M_\beta(c_s)
D^{\alpha+\beta} G_\epsilon(R)
$$

Here $D^\gamma$ is a multi-index derivative.
The current implementation precomputes $D^{\alpha+\beta} G_\epsilon(R)$ per pair as `m2l_deriv(:, pair)`.

#### 4.6 L2L

The local expansion at parent center $c_{\mathrm{parent}}$ is translated to child center $c_{\mathrm{child}}$.
With $\mathbf{d} = c_{\mathrm{child}} - c_{\mathrm{parent}}$:

$$
L_\alpha(c_{\mathrm{child}}) \mathrel{+}=
\sum_{\gamma \ge \alpha}
L_\gamma(c_{\mathrm{parent}})
\frac{\mathbf{d}^{\gamma-\alpha}}{(\gamma-\alpha)!}
$$

The shift monomials are also precomputed during `build_plan`.

#### 4.7 L2P

Let $c_{\mathrm{leaf}}$ be the center of the target leaf that contains evaluation point $x$, and let
$\mathbf{dr} = x - c_{\mathrm{leaf}}$.

$$
E_k(x) = - \sum_{|\alpha| < p} L_{\alpha + e_k}(c_{\mathrm{leaf}}) \frac{\mathbf{dr}^\alpha}{\alpha!}
$$

Here $e_k$ is the unit multi-index for axis $k$.

### 5. `build_plan` algorithm

`build_plan` performs only geometry-dependent work.

#### 5.1 Source tree

The source-coordinate bounding box is recursively split into eight octants to build the octree.
The stopping condition is either:

- source count `<= leaf_max`
- the bounding box is small enough that further subdivision is not useful

#### 5.2 Target topology

There are two target-side modes:

- `target_box` disabled:
  reuse source-tree leaves as target leaves
- `target_box` enabled:
  build a separate target tree that covers the whole box

In `periodic2`, target points are wrapped into the box before target leaf lookup.

#### 5.3 Near/far lists and M2L pair cache

For each target leaf, the source tree is traversed recursively to build near nodes and far nodes.

The well-separated test is:

$$
(r_s + r_t)^2 < \theta_{\mathrm{eff}}^2 \, \lVert\mathbf{d}\rVert^2
$$

where:

- $r_s$ is the source-node radius
- $r_t$ is the target-node radius
- $\mathbf{d}$ is the vector between node centers
- $\theta_{\mathrm{eff}} = \theta$ for both `free` and `periodic2`

In `periodic2`, a minimum-image correction is applied to $\mathbf{d}$.

Then a dual-tree recursion builds the M2L pair cache and prepares index arrays per target node.

#### 5.4 Build-time precomputation

At the end of `build_plan`, quantities that do not change between refreshes are precomputed:

- `source_parent_of`
- `parent_of`
- `source_p2m_basis`
- `m2m_term_count`, `m2m_alpha_list`, `m2m_delta_list`
- `l2l_term_count`, `l2l_gamma_list`, `l2l_delta_list`
- `source_shift_monomial`
- `target_shift_monomial`
- `shift_axis1`, `shift_axis2`
- `periodic_ewald`
- `periodic_root_operator`
- `m2l_deriv`

This makes `update_state` close to charge-dependent accumulation only.

#### 5.5 Pseudocode

```text
build_plan(src_pos, options):
  initialize_basis_tables(order)
  build_source_tree(src_pos)
  precompute_source_p2m_basis()
  build_target_topology(target_box)
  build_interactions()
  precompute_translation_operators()
  precompute_periodic2_ewald_data()
  precompute_periodic_root_operator()
  precompute_m2l_derivatives()
```

### 6. `update_state` algorithm

`update_state` corresponds to refresh in the legacy implementation.
Source coordinates are fixed; only `src_q` changes.

#### 6.1 Processing order

```text
update_state(plan, state, src_q):
  ensure_state_capacity()
  copy src_q
  clear active flags
  clear multipole/local only when the tree has no source leaves or no M2L pairs
  P2M on source leaves
  M2M bottom-up
  M2L on cached pairs
  L2L top-down
  mark state ready
```

#### 6.2 OpenMP parallelization

OpenMP is currently used in:

- one parallel region around the full `update_state`, including `src_q` copy and active-flag initialization
- `P2M`: loop over source leaves
- `M2M`: loop over nodes at the same depth
- `M2L`: loop over target nodes
- `L2L`: loop over nodes at the same depth
- translation and M2L derivative precomputation during `build_plan`

The loops are written to map roughly one node to one thread, and shared-array updates are independent at node granularity.

#### 6.3 Implementation optimizations

`update_state` avoids unnecessary work by:

- not recomputing the multi-index difference $\beta - \alpha$
- not rebuilding powers of parent-child center shifts
- precomputing the `P2M` monomial basis per source during build
- storing only valid compressed `(alpha, delta)` terms for `M2M/L2L`
- using source-node active flags to skip zero nodes in `M2L` per pair
- accumulating `M2L` contributions in thread-local `local_acc` before writing back to target-node columns
- using source-leaf-specific indices in `P2M`, not target-leaf indices

### 7. `eval_point(s)` algorithm

Evaluation proceeds as:

```text
eval_point(r):
  if plan is not built or state is not ready:
    return zero vector

  if periodic2:
    wrap r into target box

  leaf = locate_target_leaf(r)
  if leaf not found or leaf is not mapped to a leaf slot:
    use direct sum over all sources
    if periodic2 and far correction is m2l_root_oracle:
      add exact periodic Ewald correction
    return

  evaluate local expansion at leaf center
  add near direct interactions
  root local already carries periodic root correction when enabled
```

#### 7.1 Leaf lookup

- In `periodic2`, the evaluation point is wrapped into the target box before lookup.
- If a target tree exists, its leaves are used.
- If no target tree exists, source-tree leaves are used.
- If lookup fails, or the leaf cannot map to a tree leaf slot, evaluation falls back to direct sum.

#### 7.2 Near direct

Source indices in the near list are evaluated by direct sum.
In `periodic2`, image shifts in `[-N, N] x [-N, N]` are handled explicitly.
Fallback uses the same direct kernel; if explicit `m2l_root_oracle` is enabled in `periodic2`, the oracle correction is added separately.

#### 7.3 Out-of-box fallback

When a dual-target tree is used, evaluation points can leave the target box.
Then there is no target leaf, so evaluation falls back to direct sum over all sources.
With explicit `m2l_root_oracle`, the same exact periodic correction used as the build-time Ewald-fit teacher is added to direct fallback.

#### 7.4 Location of root correction

The `m2l_root_oracle` root correction is injected into `state%local(:, root)` during `update_state`.
Therefore normal leaf evaluation in `eval_point(s)` does not recompute the root correction; it just uses the local expansion carried by `state`.

### 8. `periodic2` and far correction

#### 8.1 `periodic2`

`periodic2` means exactly two axes are periodic and the remaining axis is open.

The near image sum explicitly adds the finite images:

$$
i, j \in [-N, N]
$$

M2L uses the same image-shift set and precomputes each pair derivative as an image sum.

#### 8.2 `periodic2` Ewald (Ewald2P) correction

`bem_coulomb_fmm_periodic_ewald.f90` implements an Ewald-form correction for the two-periodic, one-open Coulomb field.
Here `exact` means the finite sum actually evaluated by the code. It is not the theoretical infinite sum; it is a build-time oracle whose real-space and reciprocal-space cutoffs are controlled by `field_periodic_image_layers = N` and `field_periodic_ewald_layers = L`.

##### 8.2.1 Notation

Let the periodic axes be `a_1, a_2` and the open axis be `f`.
Define periodic lengths, cell area, image set, and reciprocal-lattice set as:

$$
L_1 = \operatorname{periodic\_len}(1),\qquad
L_2 = \operatorname{periodic\_len}(2),\qquad
A = L_1 L_2
$$

$$
\mathcal I_N = \{(i,j)\in\mathbb Z^2 \mid |i|,|j|\le N\},\qquad
\mathcal K_L = \{(m,n)\in\mathbb Z^2 \mid |m|,|n|\le L,\ (m,n)\neq(0,0)\}
$$

Image shifts and reciprocal-lattice vectors are:

$$
\mathbf L_{ij} = iL_1\,\mathbf e_{a_1} + jL_2\,\mathbf e_{a_2},\qquad
\mathbf k_{mn} = 2\pi\left(\frac{m}{L_1}\mathbf e_{a_1} + \frac{n}{L_2}\mathbf e_{a_2}\right)
$$

For source position \(\mathbf s\) and evaluation point \(\mathbf r\), define:

$$
\mathbf R_{ij} = \mathbf r - \mathbf s - \mathbf L_{ij},\qquad
R_{ij} = \lVert\mathbf R_{ij}\rVert,\qquad
z = (\mathbf r - \mathbf s)\cdot \mathbf e_f
$$

Below, \(\alpha =\) `field_periodic_ewald_alpha` and \(\epsilon =\) `softening`.

##### 8.2.2 Real-space term

The screened Coulomb field implemented by `add_screened_point_charge` is:

$$
\mathbf E_\alpha(\mathbf R) =
q\left(
\frac{\operatorname{erfc}(\alpha R)}{R^3}
{}+\frac{2\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 R^2}}{R^2}
\right)\mathbf R
$$

It is the gradient of the potential:

$$
\Phi_\alpha(\mathbf R) = q\,\frac{\operatorname{erfc}(\alpha R)}{R}
$$

The direct kernel used by `add_softened_point_charge` is:

$$
\mathbf E_\epsilon(\mathbf R) =
q\,\frac{\mathbf R}{(R^2+\epsilon^2)^{3/2}}
$$

It uses the same softening as the normal runtime direct path.

The implemented real-space correction is:

$$
\mathbf E_{\mathrm{real}} =
\sum_{(i,j)\in\mathcal I_{N+L}} \mathbf E_\alpha(\mathbf R_{ij})
{}- \sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
$$

Terms with `r2 <= tiny(1.0d0)` are skipped, so self-interaction is excluded.
If the direct fallback contribution \(\sum_{(i,j)\in\mathcal I_N}\mathbf E_\epsilon\) is added to `add_periodic2_exact_ewald_correction_single_source`, the softened inner-image part cancels and the outer shell is replaced by the screened form.

##### 8.2.3 Reciprocal-space term

For \((m,n)\neq(0,0)\), `add_exact_periodic2_reciprocal_space_correction` defines:

$$
\theta_{mn} = \mathbf k_{mn}\cdot(\mathbf r-\mathbf s),\qquad
k_{mn} = \lVert\mathbf k_{mn}\rVert
$$

$$
G^\pm_{mn}(z) =
e^{\pm k_{mn} z}\operatorname{erfc}\!\left(\frac{k_{mn}}{2\alpha}\pm \alpha z\right)
$$

and uses:

$$
\mathbf E_{\mathrm{rec}} =
q \sum_{(m,n)\in\mathcal K_L}
\frac{\pi}{A}
\begin{pmatrix}
\frac{(\mathbf k_{mn})_{a_1}}{k_{mn}}\sin\theta_{mn}\,\bigl(G^+_{mn}(z)+G^-_{mn}(z)\bigr) \\
\frac{(\mathbf k_{mn})_{a_2}}{k_{mn}}\sin\theta_{mn}\,\bigl(G^+_{mn}(z)+G^-_{mn}(z)\bigr) \\
\cos\theta_{mn}\,\bigl(G^-_{mn}(z)-G^+_{mn}(z)\bigr)
\end{pmatrix}
$$

In code these correspond to `term_p`, `term_m`, and `pair_sum`.
This term represents the high-frequency reciprocal-lattice components excluding `k=0`.

##### 8.2.4 `k=0` term

The zero-mode correction implemented by `add_exact_periodic2_k0_correction` is:

$$
\mathbf E_0 =
q\,\frac{2\pi}{A}\operatorname{erf}(\alpha z)\,\mathbf e_f
$$

The single-source oracle keeps this form as the `k=0` electric-field contribution.

##### 8.2.5 Implemented correction

Together, the correction added by `add_periodic2_exact_ewald_correction_single_source` for one source is:

$$
\mathbf E_{\mathrm{corr}} =
\mathbf E_{\mathrm{real}}
{}+ \mathbf E_{\mathrm{rec}}
{}+ \mathbf E_0
$$

`add_periodic2_exact_ewald_correction_all_sources` first sums this over all sources.

##### 8.2.6 `charged_walls` total-charge correction

For the non-neutral slab `charged_walls` closure, `add_periodic2_exact_ewald_correction_all_sources` adds this total-charge correction after summing all sources:

$$
\mathbf E_{\mathrm{walls}}(z) =
\begin{cases}
\frac{2\pi Q_{\mathrm{tot}}}{A}\,\mathbf e_f, & z < z_{\mathrm{low}}, \\
0, & z_{\mathrm{low}} \le z \le z_{\mathrm{high}}, \\
{}-\frac{2\pi Q_{\mathrm{tot}}}{A}\,\mathbf e_f, & z > z_{\mathrm{high}}
\end{cases}
$$

Here `A = L_1 L_2` is the periodic-cell area, `Q_tot = \sum_j q_j`, and `z_low/high` are the nonperiodic-axis bounds of `target_box_min/max`.
This term corresponds to the field from two compensation walls. It cancels exactly inside the slab, so it does not affect a root oracle built inside the target box or normal particle advancement. It affects only direct fallback evaluations outside the target box.

If `field_periodic_ewald_alpha <= 0`, `resolve_periodic2_ewald_alpha` selects:

$$
\alpha = \frac{1.2}{(N+1)\min(L_1,L_2)}
$$

automatically. If `min(L_1,L_2) <= 0`, it sets `alpha = 0` and disables the oracle.
Internally, `kmax = max(1, field_periodic_ewald_layers)` defines the reciprocal-space finite sum.

The actual runtime direct fallback is:

$$
\mathbf E_{\mathrm{fallback}} =
\sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
{}+
\mathbf E_{\mathrm{corr}}
{}+
\mathbf E_{\mathrm{walls}}
$$

In the build-time fit for `m2l_root_oracle`, check points are inside the target box, so `\mathbf E_{\mathrm{walls}} = 0`; the teacher uses only the single-source `\mathbf E_{\mathrm{corr}}`.
Because the `periodic_root_operator` does not use the constant potential mode, the monopole column is fixed to zero.

##### 8.2.7 `m2l_root_oracle`

`m2l_root_oracle` is an explicit opt-in, high-cost diagnostic mode that uses this Ewald2P correction as the teacher and fits an operator from root multipole to root local at proxy/check points.
Normal production runs use `none`.

- `periodic_image_layers = N`: near image shells left explicit at runtime
- `periodic_ewald_layers = L`: build-time oracle real-space outer shell `N < max(|i|,|j|) <= N+L` and reciprocal cutoff `|m|, |n| <= L`
- `periodic_ewald_alpha = alpha`: Ewald decomposition parameter, auto-selected when `<= 0`
- during build, exact periodic Ewald correction is evaluated at check points, and the field residual is fit by least squares to create the root-local operator
- at runtime, only `local(:, root) += T_root_oracle * multipole(:, root)` is added, so the eval path does not contain the Ewald implementation
- outside-tree fallback directly adds exact periodic correction to direct sum, reducing periodic residual outside the target box
- the fit uses field, not potential, and fixes the constant potential mode of the local expansion to zero

### 9. Interpreting computational cost

With fixed order $p$ and bounded interaction lists, practical costs are approximately:

- `build_plan`: close to $O(n \log n)$
- `update_state`: close to $O(n)$
- `eval_point`: close to $O(\log n + n_{\mathrm{near}} \, n_{\mathrm{img}}^2)$
- `eval_points`: parallel execution of the above point evaluation for each target

The constant factors depend strongly on:

- `order`
- `theta`
- `leaf_max`
- `periodic_image_layers`
- `periodic_ewald_layers`
- whether a target tree exists

### 10. Current implementation limits

This FMM core is not a generic kernel FMM.

- kernel is fixed to Coulomb
- the simulator adapter default order is `order = 4`
- source coordinates are considered immutable after `build_plan`
- supported boundaries are `free` and `periodic2`
- `periodic2` requires exactly two periodic axes
- far correction modes are `none` by default, `auto`, and `m2l_root_oracle`; `periodic2` `auto` normalizes to `none`, and `m2l_root_oracle` is explicit opt-in
- `eval_point(s)` return values do not include `k_coulomb`

### 11. Implementation mapping

Main implementation locations:

- Public API / wrapper:
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_build.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_state.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90`
- Shared type definitions:
  `fmm_options_type`, `fmm_plan_type`, `fmm_state_type`
  in `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90`
- Plan construction:
  `build_plan`
  in `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90`
- Charge refresh:
  `update_state`, `p2m_leaf_moments`, `m2m_upward_pass`, `m2l_accumulate`, `l2l_downward_pass`
  in `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90`
- Evaluation:
  `eval_point`, `eval_points`
  in `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90`
- periodic2 helpers:
  `has_valid_target_box`, `use_periodic2_m2l_root_oracle`,
  `use_periodic2_root_operator`, `build_periodic_shift_values`, `add_point_charge_images_field`,
  `wrap_periodic2_point`, `apply_periodic2_minimum_image`, `distance_to_source_bbox`,
  `distance_to_source_bbox_periodic`
  in `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90`
- periodic2 Ewald/oracle:
  `resolve_periodic2_ewald_alpha`, `precompute_periodic2_ewald_data`,
  `add_periodic2_exact_ewald_correction_single_source`, `add_periodic2_exact_ewald_correction_all_sources`
  in `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90`
- periodic2 root operator:
  `precompute_periodic_root_operator`
  in `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`
- BEACH adapter:
  `src/physics/field_solver/bem_field_solver_config.f90`,
  `src/physics/field_solver/bem_field_solver_tree.f90`,
  `src/physics/field_solver/bem_field_solver_eval.f90`

Design responsibilities:

- Core:
  geometry preprocessing, expansion-coefficient updates, near direct, point evaluation
- BEACH adapter:
  build `src_pos` from `mesh_type`, pass `q_elem` into `src_q`, and multiply by `k_coulomb` at the end

---

## 15. `batch_duration` stability and steady value

This section organizes the theoretical relationship between `sim.batch_duration` in the BEACH batch loop, or the one-batch physical time determined from `sim.batch_duration_step`, and the validity and stability of the converged wall-charge distribution.
In the current implementation, when `batch_duration_step` is used, `sim.batch_duration = sim.dt * sim.batch_duration_step`; in `reservoir_face` injection, the physical inflow per batch determines the macro-particle count or weight.

Implementation entry points:

- Batch procedure: [`SPEC.md` §4](https://github.com/Nkzono99/BEACH/blob/main/SPEC.md)
- Parameter definitions: [Parameters](Parameters.en.html) for `sim.batch_duration` / `sim.batch_duration_step`
- Injection usage: `src/particles/bem_injection.f90` (`reservoir_face` / `photo_raycast`)
- Batch generation and weight resolution: `src/config/bem_app_config_runtime.f90`

### 1. Reduction to a continuous-time model

Let $q_j(t)$ be the accumulated charge of insulator wall element `j`, and let $J_j(\mathbf q)$ be the incident charge flux per unit wall area at that charge state.
The absorption-only model becomes:

$$
\frac{dq_j}{dt} \;=\; J_j(\mathbf q)\, A_j
$$

where $A_j$ is the element area. Since $J$ depends on the field created by wall charge, it is generally **nonlinear**.

One BEACH batch can be viewed as an **explicit update that freezes the field at the start of the batch**. In expectation:

$$
\mathbf q^{n+1} \;=\; \mathbf q^n \;+\; \Delta t_b \cdot \mathbf J(\mathbf q^n)\,\mathbf A \;+\; \boldsymbol\eta^n
$$

where:

- $\Delta t_b = $ `sim.batch_duration`
- $\mathbf A$ is the element-area vector
- $\boldsymbol\eta^n$ is Monte Carlo sampling error within the batch

The implementation follows this picture: particles in a batch see the same field $E(\mathbf q^n)$, and the charge delta is applied to the wall at the end of the batch.
Thus `batch_duration` is the time step of this explicit update.

### 2. Validity of the steady value

Write the mean update map as:

$$
F_{\Delta t_b}(\mathbf q) \;=\; \mathbf q \;+\; \Delta t_b\, \mathbf J(\mathbf q)\,\mathbf A
$$

Its fixed point $\mathbf q^\*$ satisfies:

$$
F_{\Delta t_b}(\mathbf q^\*) = \mathbf q^\*
\quad\Longleftrightarrow\quad
\mathbf J(\mathbf q^\*) = 0
$$

Therefore, **the fixed point of the mean model itself does not depend on $\Delta t_b$**.

In this sense, if the iteration converges stably and Monte Carlo error is sufficiently averaged, changing `batch_duration` does not change the targeted continuous-time steady solution.

However, this statement applies only to the **fixed point of the mean model**. Actual runs include:

- finite-sample error per batch
- fluctuation in monitoring quantities used to judge convergence
- residual error from stopping at a finite batch count

So the observed converged value can retain weak `batch_duration` dependence. The safe statement is:

> The mean fixed point of the iteration is independent of `batch_duration`, but finite-sample and finite-time calculations can show small step-size dependence.

### 3. Linear stability

Near a fixed point $\mathbf q^\*$, define perturbations $\delta\mathbf q^n = \mathbf q^n - \mathbf q^\*$.
The linearized mean update is:

$$
\delta \mathbf q^{n+1} \;=\; \bigl(I + \Delta t_b\, M\bigr)\,\delta \mathbf q^n,
\qquad
M_{ij} \;\equiv\; \frac{\partial (J_i A_i)}{\partial q_j}\bigg|_{\mathbf q^\*}
$$

The stability condition for a general multi-degree-of-freedom system is the spectral-radius condition:

$$
\rho\!\left(I + \Delta t_b\, M\right) < 1
$$

For each eigenvalue $\lambda_k$:

$$
|1 + \Delta t_b\, \lambda_k| < 1
$$

This is the essential BEACH stability condition.

As an insulator wall accumulates charge, it tends to attract fewer particles of the same sign, so the dominant eigenvalues of $M$ are expected to be real negative, $\mathrm{Re}(\lambda_k) < 0$.
Only under this **real-negative dominant mode assumption**, using response time scale $\tau_k \equiv 1/|\lambda_k|$, the fastest mode avoids divergence when:

$$
0 \;<\; \Delta t_b \;<\; \frac{2}{|\lambda_{\max}|} \;=\; 2\,\tau_{\min}
$$

and converges monotonically, or overdamped, when:

$$
0 \;<\; \Delta t_b \;<\; \frac{1}{|\lambda_{\max}|} \;=\; \tau_{\min}
$$

Practical interpretation:

- $\Delta t_b < 2\,\tau_{\min}$: non-divergence under the real-negative dominant mode assumption
- $\Delta t_b < \tau_{\min}$: monotone convergence under the same assumption
- in a general coupled system, the precise condition is $\rho(I + \Delta t_b\, M) < 1$

Thus the $2\tau$ / $\tau$ rule is better described as an **explicit-Euler stability guide under a real-negative dominant mode assumption**, not as a strict BEACH CFL condition.

### 4. Relation to Monte Carlo noise

For a one-mode approximation with noise:

$$
\delta q^{n+1} \;=\; \left(1 - \frac{\Delta t_b}{\tau}\right)\,\delta q^n \;+\; \xi^n
$$

the steady variance depends on the variance of $\xi^n$.
The key point is that **the $\Delta t_b$ dependence of $\mathrm{Var}(\xi^n)$ depends on injection normalization**.
BEACH `reservoir_face` has two modes.

#### 4.1 Fixed `w_particle`

When `w_particle` is specified directly, the physical inflow count changes in proportion to $\Delta t_b$, so the expected macro-particle count per batch follows:

$$
N_\text{macro} \;\propto\; \Delta t_b
$$

The shot-noise variance of the batch charge increment can be regarded roughly as proportional to $\Delta t_b$:

$$
\mathrm{Var}(\xi^n) \;\approx\; \alpha\, \Delta t_b
$$

In the limit $\Delta t_b \ll \tau$, the steady variance does not depend strongly on `batch_duration`.

#### 4.2 Fixed `target_macro_particles_per_batch`

When `w_particle` is solved from `target_macro_particles_per_batch`, the weight is determined as in `src/config/bem_app_config_runtime.f90:644`:

$$
w_\text{particle} \;\propto\; \frac{\Gamma\, A\, \Delta t_b}{N_\text{target}}
$$

so the noise dependence on $\Delta t_b$ differs from the simple $\mathrm{Var}(\xi^n) \propto \Delta t_b$ of §4.1.
The macro-particle count is fixed, while the contribution per particle is proportional to $\Delta t_b$.

#### 4.3 Practical interpretation

The useful separation is:

- `batch_duration` mainly controls **deterministic stability**
- the main knobs for statistical noise are **`w_particle` or `target_macro_particles_per_batch`**

In particular, neither of these is generally true:

> Making `batch_duration` smaller always lowers noise.
> Making `batch_duration` larger leaves noise almost unchanged.

The answer depends on injection normalization.

### 5. Physical estimate of $\tau_{\min}$

$\tau_{\min}$ is the fastest effective response time that controls numerical stability.
It is hard to express with one general physical formula because it depends on geometry, potential distribution, upstream distribution function, and injection model.
In practice, estimate two different quantities.

#### 5.1 Charging / sheath relaxation time

A natural estimate uses an effective capacitance $C_\text{eff}$ and effective conductance $G_\text{eff}$:

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}}{G_\text{eff}}
$$

or a typical potential change $\Delta\phi$ and effective current $I_\text{eff}$:

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}\,\Delta\phi}{I_\text{eff}}
$$

This is a relatively slow charging timescale affected by geometry and shielding.

#### 5.2 Inverse plasma frequency

Another fast reference is:

$$
\tau_{pe} \;=\; \omega_{pe}^{-1} \;=\; \sqrt{\frac{\varepsilon_0\, m_e}{n_e\, e^2}}
$$

This is the microscopic fast timescale of an electron plasma and is useful as a reference for how sharply the system can respond.

However, treating $\omega_{pe}^{-1}$ directly as an upper bound on $\tau_{\min}$ is too strong.
It is better viewed as a **fast-side physical reference**. The effective time constant that limits `batch_duration` often comes from $\tau_\text{charge}`, including geometry and incoming-flux limits.

#### 5.3 Practical choice

For $\tau_{\min}$, estimate:

- $\omega_{pe}^{-1}$: microscopic fast reference
- $\tau_\text{charge}$: system-specific charging / sheath relaxation timescale

Then refine with numerical experiments.

> $\omega_{pe}^{-1}$ is only a fast reference; the actual stability limit is set by an effective response time that often includes $\tau_\text{charge}$.

### 6. Practical usage

1. Estimate both $\omega_{pe}^{-1}$ and $\tau_\text{charge}$ as physical scales.
2. Start with a conservative, smaller `batch_duration` if oscillation should be avoided.
3. Compare charge history and monitoring quantities with `batch_duration` multiplied by 1/2 and 2, as a **step-size sensitivity check**.
4. If the converged values nearly agree and no oscillation or divergence is visible, the `batch_duration` is practically adequate.
5. If noise is large, first adjust `w_particle` or `target_macro_particles_per_batch`. Do not try to solve noise only by changing `batch_duration`.
6. Oscillation in `charge_history.csv` `last_rel_change`, or jitter in element charge time series, is a useful diagnostic. This is better called a **step-size sensitivity check** than strict Richardson extrapolation, because it does not assume a power law of the error.

### 7. Summary

| Item | Conclusion |
|---|---|
| Validity of steady value | The fixed point of the mean update does not depend on `batch_duration` |
| Exact stability condition | $\rho(I + \Delta t_b\, M) < 1$ |
| $2\tau$, $\tau$ rules | Explicit-Euler approximate guide under real-negative dominant modes |
| Role of $\omega_{pe}^{-1}$ | Microscopic fast reference, not generally a direct stability upper bound |
| Noise and `batch_duration` | Dependence is set by injection normalization |
| Main noise-reduction knobs | Adjust `w_particle` or `target_macro_particles_per_batch` |
| Practical check | Step-size sensitivity check by varying `batch_duration` |

It is theoretically clean to say that the steady value of the mean model does not depend on how `batch_duration` is chosen.
The general stability condition $\rho(I + \Delta t_b M) < 1$ follows directly from classical stability analysis.
The remaining uncertainty is the value of $\tau_{\min}$ itself, which must be narrowed down with both case-specific physical estimates and numerical experiments.

### Related documents

- [Fortran parameter file specification](Parameters.en.html) — how to set `sim.batch_duration` / `sim.batch_duration_step`
- [Fortran-centered workflow](Workflow.en.html) — batch-loop execution control
- `SPEC.md` — one-batch computation procedure and stop conditions
