title: FMM

Lang: [日本語](FMM.md) | [English](FMM.en.md)

# FMM

The Fast Multipole Method (FMM) reuses a plan built from fixed geometry and compresses distant sources into expansion
coefficients. Each batch refreshes charge-dependent coefficients; each particle target combines the far expansion with a near
Direct sum.

| Suitable calculation | Check first |
| --- | --- |
| Many elements and many particle steps per batch | Difference from Direct and release-build runtime |
| Large triangle P0 mesh | Difference from panel Direct and mesh refinement |
| `periodic2` field | Image sum, nonzero/zero modes, and zero-mode boundary condition |

```toml
[sim]
field_solver = "fmm"

[domain]
box_min = [-0.5, -0.5, -0.1]
box_max = [ 0.5,  0.5,  1.0]

[field_boundary]
mode = "free"
```

The `[domain]` box also defines the target tree covering the particle region. For a free field, a target outside that tree
falls back to a Direct sum over every source.

## Combine far expansions with near Direct sums

FMM partitions sources and targets into octrees and separates near and far interactions.

```text
source charge
   │ P2M: form leaf multipoles
   ▼
source multipoles
   │ M2M: aggregate children into parents
   ▼
far interactions
   │ M2L: translate into target-node local expansions
   ▼
target locals
   │ L2L: propagate parents into children
   ▼
L2P at target + near Direct interactions
```

Far interactions use Cartesian multipole and local expansions. Sources in the near list use the analytic panel kernel
directly. The current simulator adapter fixes the expansion order at `order = 4`; it is not a user configuration value.

### Cartesian-expansion notation

For multi-index $\alpha=(\alpha_x,\alpha_y,\alpha_z)$, define

$$
|\alpha|=\alpha_x+\alpha_y+\alpha_z,qquad
\alpha!=\alpha_x!\alpha_y!\alpha_z!,qquad
\mathbf r^\alpha=r_x^{\alpha_x}r_y^{\alpha_y}r_z^{\alpha_z}.
$$

The expansion uses the Laplace kernel

$$
G(\mathbf r)=\frac{1}{\lVert\mathbf r\rVert},\qquad
\mathbf E=-\nabla\phi.
$$

The following coefficients omit the Coulomb constant. The BEACH adapter applies it together with unit scaling after evaluation.

### Aggregate sources with P2M and M2M

For leaf center $\mathbf c$, each triangle-P0 source contributes the exact surface average

$$
M_{i,\alpha}=\frac{q_i}{A_i}\int_{T_i}
\frac{(\mathbf y-\mathbf c)^\alpha}{\alpha!}\,dA_{\mathbf y}.
$$

To translate child center $\mathbf c_c$ to parent center $\mathbf c_p$, let $\mathbf d=\mathbf c_c-\mathbf c_p$ and sum

$$
M_\beta(\mathbf c_p)=
\sum_{\alpha\le\beta}M_\alpha(\mathbf c_c)
\frac{\mathbf d^{\beta-\alpha}}{(\beta-\alpha)!}
$$

over every child. P2M bases and M2M shift monomials depend only on geometry and are precomputed with the plan.

### Convert multipoles into target locals with M2L

For source center $\mathbf c_s$, target center $\mathbf c_t$, and $\mathbf R=\mathbf c_t-\mathbf c_s$, a
well-separated node pair contributes

$$
L_\alpha(\mathbf c_t)\mathrel{+}=
\sum_\beta(-1)^{|\beta|}M_\beta(\mathbf c_s)
D^{\alpha+\beta}G(\mathbf R).
$$

$D^\gamma$ is a multi-index derivative. The plan stores source/target combinations in the M2L pair cache and stores
$D^{\alpha+\beta}G(\mathbf R)$ in `m2l_deriv`. A batch `update_state` therefore performs mainly products and sums
with current multipole coefficients.

### Propagate and evaluate with L2L and L2P

For parent-to-child displacement $\mathbf d=\mathbf c_c-\mathbf c_p$,

$$
L_\alpha(\mathbf c_c)\mathrel{+}=
\sum_{\gamma\ge\alpha}L_\gamma(\mathbf c_p)
\frac{\mathbf d^{\gamma-\alpha}}{(\gamma-\alpha)!}.
$$

For a target $\mathbf x$ in a leaf centered at $\mathbf c_l$, let $\mathbf r=\mathbf x-\mathbf c_l$. Then

$$
\phi_\mathrm{far}(\mathbf x)=
\sum_\alpha L_\alpha(\mathbf c_l)\frac{\mathbf r^\alpha}{\alpha!},
$$

$$
E_{\mathrm{far},k}(\mathbf x)=
-\sum_\alpha L_{\alpha+\mathbf e_k}(\mathbf c_l)
\frac{\mathbf r^\alpha}{\alpha!}.
$$

The fixed P0 triangle kernel then evaluates the leaf near list directly. M2L does not represent every interaction by
itself: an FMM result is always the sum of the far local expansion and near Direct interactions.

## Reuse geometry separately from charge updates

To reuse a fixed mesh over many batches, FMM separates geometry-dependent `plan` data from charge-dependent `state` data.

| Data | Main contents | Update time |
| --- | --- | --- |
| plan | Source tree, target tree, near/far lists, P2M basis, translation operators | Initialization; rebuild if geometry, source count, or major options change |
| state | Current `q_elem`, multipole coefficients, local coefficients | Whenever the field snapshot is refreshed |

`build_plan` precomputes geometry-only quantities. `update_state` performs P2M, M2M, M2L, and L2L for current element charges.
Deposited charge is committed at batch end and enters `state` at the next field snapshot.

## Cover the target tree with `[domain]`

In public configuration, `[domain]` defines the simulation and target-tree extent through `box_min` / `box_max` or
`box_origin` / `box_size`. Cover every position particles can reach. `domain.periodic_axes` owns topology, while
`field_boundary.mode` selects the `free` or `periodic2` field closure.

For a free field, FMM sums every source directly when no target leaf is found; frequent fallback reduces performance to
Direct-like cost. `cached_kneq0` requires a fixed target topology and rejects out-of-box targets.

<a id="source-kernel"></a>

## Treat triangle P0 consistently in near and far fields

Triangle P0 retains each complete triangle as source geometry. Source-tree boxes contain all three vertices rather than only
the centroid, and P2M uses exact surface averages of triangle monomials. Near interactions and out-of-tree fallback use the same
analytic panel integral as [Direct](DirectSolver.en.html#triangle-p0).

Triangle P0 requires non-degenerate triangles and resolved vacuum sides. If a target lies exactly on a panel, the current FMM core uses the
principal-value trace for the panel field. Particle collision events normally stop trajectories before such an evaluation, but
surface-field verification must account for the difference from the Direct vacuum-side trace.

See [`examples/panel_fmm.toml`](../examples/panel_fmm.toml) for a complete configuration.

## Tune accuracy through the near/far split

For source-node radius $r_s$, target-node radius $r_t$, and center distance $d$, a node pair is conceptually well separated when

$$
(r_s+r_t)^2 < \theta^2 d^2.
$$

Well-separated pairs enter far expansions. Other pairs are refined and ultimately evaluated with near Direct interactions.

| Key | Effect |
| --- | --- |
| `tree_theta` | Smaller values make far acceptance stricter and are generally more accurate and slower |
| `tree_leaf_max` | Changes near Direct work, tree depth, and memory balance |
| `field_normalization` | Changes numerical coordinate scale, not physical units or the model |

If `tree_theta` and `tree_leaf_max` are absent from the input, even explicit `fmm` uses element-count settings
`(0.40, 12)`, `(0.50, 16)`, `(0.58, 20)`, and `(0.65, 24)`, with boundaries at 1,500, 10,000, and 50,000 elements.
They are starting values, not error guarantees.

Because the expansion order is currently fixed at four, users cannot perform an order-convergence sweep. Measure solver
sensitivity with `tree_theta`, `tree_leaf_max`, and Direct comparisons; measure source-discretization convergence separately by
refining the mesh.

## Evaluate fields, potentials, and history output

At an ordinary target, both field and potential use a local expansion plus near Direct interactions. For element-center output
in `potential_history.csv`, triangle P0 includes the analytic panel self integral in its near kernel. See
[Direct](DirectSolver.en.html#potential-at-element-centers) for the definition.

When potential history is written, state is refreshed using the latest element charge. Enabling it can therefore add a separate
refresh and evaluation over every element target beyond the normal batch field update.

## Connect periodic2 far correction to local expansions

For `periodic2`, the smooth field beyond the finite images is inserted as a root-multipole-to-target-local operator after M2L
and before L2L. See [periodic2 Far Correction](PeriodicFarCorrection.en.html) for `none`, production `cached_kneq0`, the
Ewald2P teacher, and cache behavior; see [periodic2 electrostatics](PeriodicElectrostatics.en.html) for ownership of the complete
field including zero mode. Small reference cases can combine `field_solver="direct"` with
`panel_spectral_reference`.

## Measure setup cost separately from particle evaluation

For fixed order and bounded interaction lists, useful estimates are:

| Operation | Estimate |
| --- | --- |
| Plan construction | Near $O(N\log N)$; one-time, but includes translation and periodic-operator preparation |
| State update | Near $O(N)$ |
| One target | Near $O(\log N+N_{\mathrm{near}})$ |
| Many targets | Target evaluations are parallelized with OpenMP |

Actual constants and memory depend on the number of expansion coefficients, source and target nodes, M2L pairs, near-list size,
and periodic images. For a small problem, plan and state fixed costs can exceed Direct. Compare performance in a release profile,
separating initialization, batch refresh, and particle-evaluation time.

## Converge from Direct comparisons through charging results

Use the following order:

1. Run a reduced Direct and FMM case with identical mesh and normalization.
2. Compare field and potential near surfaces, far away, across strong charge gradients, and in cancellation regions.
3. Reduce `tree_theta`, vary `tree_leaf_max`, and check stability of the observables.
4. Refine the source mesh to measure kernel-discretization convergence independently from FMM approximation.
5. Compare impact elements, escape counts, post-batch `q_elem`, and the charge ledger.
6. For periodic2, separately check finite image layers, far correction, zero mode, and cold/warm cache agreement.

Relative error is unstable where the reference field is nearly zero, so report an absolute error or an error normalized by a
representative field as well. A small pointwise FMM error does not remove errors from a coarse panel mesh or a large particle time
step.

## Choices supported by the current implementation

- Coulomb kernel only
- expansion order fixed at four in the simulator adapter
- source geometry fixed after plan construction
- free and periodic2 field boundaries
- analytic triangle-P0 near interactions and exact P2M
- Direct fallback outside the target tree for a free field
- the physical periodic zero mode is added exactly once outside the FMM core
- the FMM core itself does not solve the outer plasma/sheath; matching-plane closure is composed outside the core

Coefficient arrays, interaction caches, translation precomputation, and parallel loops are documented in
[Coulomb FMM internals](FMMCore.en.html). Follow [periodic2 Far Correction](PeriodicFarCorrection.en.html) for periodic-operator
internals.

## Code reference

- BEACH adapter and fixed order: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- Plan/state refresh: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- FMM public API: [`bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90)
- Plan and interaction construction: [`bem_coulomb_fmm_plan_ops.f90`](../src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90)
- State upward/downward passes: [`bem_coulomb_fmm_state_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90)
- Target evaluation and fallback: [`bem_coulomb_fmm_eval_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90)
- Main solver regressions: [`test_dynamics_fmm.f90`](../tests/fortran/test_dynamics_fmm.f90), [`test_dynamics_panel_fmm.f90`](../tests/fortran/test_dynamics_panel_fmm.f90)
