title: FMM

Lang: [日本語](FMM.md) | [English](FMM.en.md)

# FMM

The Fast Multipole Method (FMM) reuses a computation plan built from fixed geometry and compresses distant interactions into
expansion coefficients. Each batch refreshes charge-dependent coefficients; each particle target combines the far expansion
with a near Direct sum. It reduces work for many boundary elements and many particle targets while its error is measured against
Direct.

| Suitable calculation | Check first |
| --- | --- |
| Many elements and many particle steps per batch | Difference from Direct and release-build runtime |
| Large triangle P0 mesh | Difference from panel Direct and mesh refinement |
| Legacy `periodic2` field | Image sum, nonzero/zero modes, and outer-model composition |

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "free"
use_box = true
box_min = [-0.5, -0.5, -0.1]
box_max = [ 0.5,  0.5,  1.0]
```

`use_box=true` creates a target tree covering the particle domain. A target outside that tree is evaluated by a Direct sum over
every source.

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

Far interactions use Cartesian multipole and local expansions. Sources in the near list use the selected point or panel kernel
directly. The current simulator adapter fixes the expansion order at `order = 4`; it is not a user configuration value.

## Reuse geometry separately from charge updates

To reuse a fixed mesh over many batches, FMM separates geometry-dependent `plan` data from charge-dependent `state` data.

| Data | Main contents | Update time |
| --- | --- | --- |
| plan | Source tree, target tree, near/far lists, P2M basis, translation operators | Initialization; rebuild if geometry, source count, or major options change |
| state | Current `q_elem`, multipole coefficients, local coefficients | Whenever the field snapshot is refreshed |

`build_plan` constructs trees and interaction lists and precomputes geometry-only quantities. `update_state` performs P2M, M2M,
M2L, and L2L for the current element charges. At each particle position, evaluation needs only the local expansion of its target
leaf plus near Direct interactions.

State is not refreshed particle by particle within a batch. Deposited charge is committed at the end of the batch and enters the
state used by the next field snapshot.

## Cover the particle region with the target tree

The target side of an FMM plan is formed in one of two ways:

- `use_box=false`: source-tree leaves are also used as target leaves
- `use_box=true`: an independent target tree covers `box_min` through `box_max`

If the particle evaluation domain is wider than the source bounding box, the first form makes targets outside the tree likely.
Normally, define a box covering every position particles can reach. A free-boundary calculation may still use a box, allowing
particle boundary geometry and the target-tree extent to share the same definition.

When no target leaf can be found, FMM does not extrapolate a local expansion. It sums every source directly. This fallback favors
continuity and correctness, but frequent fallback reduces performance to Direct-like cost. If many particles move outside the
box, verify both that boundary events are processed as intended and that the box covers the intended field-evaluation domain.

<a id="source-kernel"></a>

## Keep source discretization consistent in near and far fields

### Centroid point charge

Point mode uses element centroids as source positions. Near interactions use the point kernel with `sim.softening`; far
interactions are multipole expansions of the same charges. Softening is divided by $L_0$ together with coordinates in the
normalized internal calculation.

### Constant surface charge over a triangle

Triangle P0 mode retains each complete triangle as source geometry. Source-tree boxes contain all three vertices rather than only
the centroid, and P2M uses exact surface averages of triangle monomials. Near interactions and out-of-tree fallback use the same
analytic panel integral as [Direct](DirectSolver.en.html#triangle-p0).

There is therefore no near-field fallback to a centroid point source. Triangle P0 requires `softening=0`, non-degenerate
triangles, insulator surfaces, and resolved vacuum sides. If a target lies exactly on a panel, the current FMM core uses the
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
| `softening` | Applies only to the point kernel; triangle P0 requires zero |

If `tree_theta` and `tree_leaf_max` are absent from the input, even explicit `fmm` uses element-count settings
`(0.40, 12)`, `(0.50, 16)`, `(0.58, 20)`, and `(0.65, 24)`, with boundaries at 1,500, 10,000, and 50,000 elements.
They are starting values, not error guarantees.

Because the expansion order is currently fixed at four, users cannot perform an order-convergence sweep. Measure solver
sensitivity with `tree_theta`, `tree_leaf_max`, and Direct comparisons; measure source-discretization convergence separately by
refining the mesh.

## Evaluate fields, potentials, and history output

At an ordinary target, both field and potential use a local expansion plus near Direct interactions. For element-center output in
`potential_history.csv`, the point kernel adds its finite self term after FMM evaluation. Triangle P0 includes the analytic panel
self integral in its near kernel. See [Direct](DirectSolver.en.html#potential-at-element-centers) for the definitions.

When potential history is written, state is refreshed using the latest element charge. Enabling it can therefore add a separate
refresh and evaluation over every element target beyond the normal batch field update.

## Own the nonzero modes in periodic2

In the legacy `sim.field_bc_mode="periodic2"` path, targets are wrapped into the primary box and
`field_periodic_image_layers=N` explicitly covers the $[-N,N]^2$ near images. Far correction is one of:

| `field_periodic_far_correction` | Role |
| --- | --- |
| `none` | Finite image sum; default |
| `auto` | Currently normalized to `none` |
| `m2l_root_oracle` | High-cost diagnostic fitting an Ewald residual into root locals |
| `cached_kneq0` | Production path adding the infinite-periodic nonzero mode with a versioned operator |

With `cached_kneq0`, the FMM core computes the nonzero mode; the physical zero mode and outer response are composed once by the
field snapshot. The point-only `m2l_root_oracle` is unavailable for panel sources. These are not merely FMM accuracy controls:
they form a coupled configuration with particle boundaries and the outer plasma. See
[periodic2 electrostatics](PeriodicElectrostatics.en.html) and
[outer-plasma models](OuterPlasmaModels.en.html).

Small split-reference cases combine `field_solver="direct"` with a panel spectral backend. Their configuration is described in
[Infinite-periodic periodic2 with outer plasma](InfinitePeriodicOuterConfiguration.en.html).

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

1. Run a reduced Direct and FMM case with identical mesh, kernel, softening, and normalization.
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
- triangle P0 without softening, with insulator-only Phase 1 surfaces
- Direct fallback for targets outside the tree
- periodic zero mode and outer response are not completed by the FMM core alone

Expansion equations, multi-indices, translations, parallel loops, and periodic-cache generation are documented separately in
[Coulomb FMM internals](FMMCore.en.html).

## Code reference

- BEACH adapter and fixed order: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- Plan/state refresh: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- FMM public API: [`bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90)
- Plan and interaction construction: [`bem_coulomb_fmm_plan_ops.f90`](../src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90)
- State upward/downward passes: [`bem_coulomb_fmm_state_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90)
- Target evaluation and fallback: [`bem_coulomb_fmm_eval_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90)
- Main solver regressions: [`test_dynamics_fmm.f90`](../tests/fortran/test_dynamics_fmm.f90), [`test_dynamics_panel_fmm.f90`](../tests/fortran/test_dynamics_panel_fmm.f90)
