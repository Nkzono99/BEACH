title: Field Solvers and Boundary Conditions

Lang: [English](FieldSolvers.en.md) | [日本語](FieldSolvers.md)

# Field Solvers and Boundary Conditions

## Choosing a Solver First

| Goal | Recommended setting | Notes |
| --- | --- | --- |
| Small mesh smoke test | `field_solver = "auto"` or `"direct"` | `direct` is exact but costs `O(nelem)` |
| Larger production run | `field_solver = "fmm"` | See [Coulomb FMM Core Details](FMMCore.en.html) |
| Two-periodic-axis boundary | `field_bc_mode = "periodic2"` and `field_solver = "fmm"` | Exactly two axes must be periodic |
| Accuracy checks and debugging | Compare `direct` and `fmm` on a small case | Keep mesh and particle settings identical |

In the current implementation, `periodic2` is FMM-only. Treat `auto` as the convenient choice for small and medium `field_bc_mode="free"` cases.

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
[Coulomb FMM core details](FMMCore.en.html)

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
In Stage 1, let $A_n = \sum_{i \in n}|q_i|$. Even when an internal node passes the geometric criterion, its monopole is accepted only if
`abs(A_n - abs(Q_n)) <= 64 * epsilon(1.0d0) * max(A_n, abs(Q_n))`.
This machine-epsilon tolerance permits only rounding differences in charge aggregation.
Mixed-sign internal nodes descend to leaves for direct interactions, while same-sign nodes retain the existing monopole path and performance.

A future monopole-plus-dipole error criterion is expected to recover mixed-sign performance without restoring the unstable signed-charge-centroid approximation.

### 5.3 FMM

FMM mode calls a simulator-independent Coulomb FMM core.
The field-solver adapter handles:

1. Convert mesh centroids to source coordinates `src_pos(3, nelem)`.
2. Build the geometry-dependent plan with `build_plan(plan, src_pos, options)`.
3. Update the charge-dependent state with `update_state(plan, state, q_elem)`.
4. Call `eval_point(plan, state, r, e)` for each particle position.

For P2M / M2M / M2L / L2L / L2P details, see
[Coulomb FMM core details](FMMCore.en.html).

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

`field_periodic_far_correction` accepts `auto`, `none`, `m2l_root_oracle`, and `cached_kneq0`.
The current implementation normalizes `auto` to `none` for compatibility.
`m2l_root_oracle` is a diagnostic far correction enabled only by explicit request.
`cached_kneq0` is the production backend and requires x/y periodic, z open,
`exclude_k0`, and `e_bottom_zero`.

### 6.2 Near images and far correction

`field_periodic_image_layers = N` enumerates near images along the two periodic axes as:

$$
(i, j) \in [-N, N]^2
$$

The FMM core combines primary-cell sources and image sources for near contributions.
With `m2l_root_oracle`, the build stage fits a root-local correction from the Ewald residual and injects it into the root local expansion at runtime.
`cached_kneq0` stores the smooth full-periodic Ewald residual as a versioned
operator and subtracts a symmetric `k=0` state built during charge refresh.
On a cache miss only MPI rank 0 generates the operator under a filesystem lock, publishes
it with an atomic rename, and broadcasts it. Warm field evaluation and charge
refresh contain no all-source Ewald sum. The electrostatic snapshot then adds
the boundary-condition-specific `k=0` exactly once.

### 6.3 Collision side

The collision-side `periodic2` logic is separate from the field-calculation FMM.
`find_first_hit_periodic2` computes the needed image-shift range from the particle segment and the canonical mesh AABB, then tests the shifted segment against the base mesh.
If multiple candidates have the same `t`, ties are broken deterministically by element index and image shift.

---
