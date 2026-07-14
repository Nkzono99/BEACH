title: Treecode

Lang: [日本語](Treecode.md) | [English](Treecode.en.md)

# Treecode

Treecode groups sources in an octree and evaluates a distant group as one monopole. Near nodes and nodes with cancelling
positive and negative charge are expanded to a Direct sum, so the fraction of accepted far nodes determines both accuracy and
speed. It targets fewer source evaluations than Direct for medium free-boundary problems.

| Property | Description |
| --- | --- |
| Kernel | point, triangle P0 |
| Field boundary | free only |
| Far approximation | Node total charge at its center of charge |
| Near evaluation | Direct leaf sum of softened point charges or the analytic triangle-P0 panel kernel |
| Potential | Uses the same tree traversal and near/far classification as the field |

```toml
[sim]
field_solver = "treecode"
field_bc_mode = "free"
softening = 1.0e-6
```

For triangle P0, also select the element kernel and set `softening=0`.

```toml
[field]
element_kernel = "triangle_p0"
```

## Build the octree once from fixed geometry

At initialization, triangle centroids are used as source positions and an octree is built as follows:

1. Find the axis-aligned bounding box of the source centroids in a node.
2. Split at the box center into eight octants.
3. Recurse until a node contains at most `tree_leaf_max` sources.

If the centroids are numerically coincident, or splitting leaves every source in one octant, the node remains a leaf even when
it exceeds the nominal limit. Tree topology is reused while the mesh geometry remains fixed.

Triangle P0 still uses centroids for partitioning, but expands each MAC node radius to contain every triangle vertex in the
node. This prevents a panel that reaches near a target from being misclassified as a far monopole merely because its centroid
is distant.

## Refresh node charge moments for each batch

Changing surface charge does not rebuild the tree. Instead, node charge moments are accumulated from leaves to the root. For
node $n$, BEACH updates

$$
Q_n=\sum_{i\in n}q_i,
\qquad
A_n=\sum_{i\in n}|q_i|,
$$

$$
\mathbf{c}_{Q,n}=
\begin{cases}
Q_n^{-1}\sum_{i\in n}q_i\mathbf{c}_i, & |Q_n|>\mathrm{tiny},\\
\mathbf{c}_{n}, & \text{otherwise},
\end{cases}
$$

where $\mathbf{c}_n$ is the geometric node center. This refresh forms the field snapshot at the start of a batch. Element
charge committed at the end of a batch is used by the next refresh.

## Group nodes according to distance

All sources in a leaf are evaluated directly. For an internal node, the node radius $R$, distance $d$ from its center to the
target, and `tree_theta` define the far candidate condition

$$
R < \theta(d-R).
$$

A target inside the node's bounding sphere always causes traversal into the children. A smaller $\theta$ opens more nodes and
is slower and more accurate; a larger $\theta$ groups more interactions and is faster and coarser.

An accepted far node is evaluated as a monopole $Q_n$ at $\mathbf{c}_{Q,n}$ for both field and potential. Point sources use
the configured softening; triangle P0 is unsoftened. Leaves use softened point sums or analytic panel sums, so triangle-P0
near fields, on-surface jumps, and principal-value self potentials never fall back to centroid point charges. Treecode
traverses the source tree for every target; unlike FMM, it does not construct target-side local expansions.

## Expand cancelling nodes to a Direct sum

In a node containing positive and negative charge, cancellation can make $Q_n$ small and move $\mathbf{c}_{Q,n}$ far outside
the node. A single-monopole approximation is then unstable. In addition to the geometric condition, BEACH accepts a node only if

$$
|A_n-|Q_n||
\le 64\,\epsilon_{\mathrm{mach}}\max(A_n,|Q_n|).
$$

This tolerance permits only accumulation roundoff for same-sign charge. An effectively mixed-sign node is opened even when it
is geometrically distant, eventually reaching leaf Direct sums.

The guard preserves accuracy under cancellation but can reduce acceleration when signs are mixed at fine scales. For such a
distribution, measure runtime as well as error against Direct, and compare FMM when appropriate.

## Settings that control accuracy and speed

| Key | Role | Constraint |
| --- | --- | --- |
| `tree_theta` | Geometric acceptance of a far node | $0 < \theta \le 1$ |
| `tree_leaf_max` | Nominal number of sources in a leaf | at least 1 |
| `tree_min_nelem` | Switching threshold for `field_solver="auto"` | at least 1 |
| `softening` | Softening for point leaf sums and monopoles [m] | non-negative; zero for triangle P0 |

When `tree_theta` and `tree_leaf_max` are not explicitly present in the input, element-count values are selected even for an
explicit `treecode` mode.

| `nelem` | `tree_theta` | `tree_leaf_max` |
| ---: | ---: | ---: |
| `< 1500` | 0.40 | 12 |
| `1500`–`9999` | 0.50 | 16 |
| `10000`–`49999` | 0.58 | 20 |
| `>= 50000` | 0.65 | 24 |

These are starting values for speed and accuracy, not case-specific error bounds. If only one of the two keys is explicitly
set, that value overrides the table while the other still comes from the table.

## Accelerate field and potential with the same tree

The Treecode path evaluates electric fields at particle positions, `eval_potential` at arbitrary points, and
`potential_history.csv` at all element centers with the same source tree. Mesh potential traverses the tree once per element
center and reduces source evaluations relative to the $O(N^2)$ Direct sum when enough far nodes are accepted. It preserves the
area-equivalent point self term for zero softening and the analytic triangle-P0 centroid self potential.

The tree is built during initialization, node moments are refreshed after charge changes, and the tree is traversed for every
target. Well-behaved distributions require far fewer interactions per target than Direct, but the exact cost depends on tree
imbalance, $\theta$, leaf size, and the fraction of mixed-sign nodes.

## Measure approximation error against Direct

Compare against Direct with the same source kernel, softening, normalization, and mesh.

1. Sample strong-field regions, near-surface points, far points, and charge-cancellation regions.
2. Reduce `tree_theta` and verify convergence toward Direct.
3. Vary `tree_leaf_max` and measure both result and runtime sensitivity.
4. Compare particle impact elements and post-batch `q_elem`, not only fields at isolated points.
5. Measure performance in a release build with representative particle and step counts.

Regression tests separately verify that same-sign nodes retain the field and potential monopole paths, strongly cancelling
mixed-sign nodes preserve Direct accuracy, triangle-P0 near and self terms match analytic panel sums, and length normalization
returns consistent SI values.

## Code reference

- Octree construction and moment refresh: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- MAC, mixed-sign guard, and traversal: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- Parameter selection: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- Accuracy, cancellation, and normalization tests: [`test_dynamics_field_solver.f90`](../tests/fortran/test_dynamics_field_solver.f90)
- Direct-equivalence test for batch results: [`test_simulator.f90`](../tests/fortran/test_simulator.f90)
