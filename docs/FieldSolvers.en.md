title: Field evaluation

Lang: [日本語](FieldSolvers.md) | [English](FieldSolvers.en.md)

# Field evaluation

At the start of each batch, BEACH constructs a field from the current element charges, `q_elem`, and uses that same field
for all particles tracked in the batch. Charge carried to a surface is committed at the end of the batch, so it changes the
field starting with the next batch.

Field evaluation has two main choices: how one element's charge is represented and how interactions from many elements are
evaluated.

| Choice | Role | Values |
| --- | --- | --- |
| source kernel | Charge distribution on one element | `point` / `triangle_p0` |
| solver | How interactions from many sources are accumulated | `direct` / `treecode` / `fmm` / `auto` |
| field boundary | How periodic images and the far field are included | `free` / `periodic2` |

## Select a solver from problem size and source kernel

| Solver | Main use | Source kernel | Field boundary | Approximation |
| --- | --- | --- | --- | --- |
| [Direct](DirectSolver.en.html) | Small problems, reference results, and split references | point, triangle P0 | free, constrained periodic2 | Evaluates the selected kernel for every element |
| [Treecode](Treecode.en.html) | Medium free-space problems | point, triangle P0 | free | Uses monopoles for far nodes and the selected kernel in near leaves |
| [FMM](FMM.en.html) | Large problems and many targets | point, triangle P0 | free, periodic2 | Approximates far interactions with multipole and local expansions |
| `auto` | Select by element count for a free boundary | point, triangle P0 | free | Direct/Treecode for point; Direct/FMM for triangle P0 |

With `auto`, Direct is used when `nelem < tree_min_nelem`. At and above the threshold, point sources use Treecode and
triangle P0 sources use FMM. The default threshold is `256`. Runtime depends not only on the element count, but also on the
number of particles, number of steps, and target distribution, so measure a reduced case representative of the production run.

## Solver and Field-Boundary Compatibility

This table is the canonical compatibility reference for solvers, source kernels, and field boundaries.

| Solver | `free` | `periodic2` |
| --- | --- | --- |
| `direct` | `point` / `triangle_p0` | Split reference only. Requires `triangle_p0`, `periodic2.nonzero_mode_backend="panel_spectral_reference"`, `zero_mode_policy="exclude_k0"`, and a supported lower-boundary model |
| `treecode` | `point` / `triangle_p0` | Unsupported |
| `fmm` | `point` / `triangle_p0` | Supported. Infinite-periodic production runs use `cached_kneq0` |
| `auto` | `point` / `triangle_p0` | Unsupported |

`periodic2` additionally requires `sim.use_box=true`, exactly two periodic axes, and one open axis.
The Direct split reference is intended for reduced reference and validation cases; the normal periodic2 production path uses FMM.
The [official Direct split-reference configuration](../examples/periodic2_linear_outer_reference.toml) provides a complete example.

## Choose element-charge discretization with the source kernel

### Point

`field.element_kernel="point"` places each element's total charge at its triangle centroid. `sim.softening` can regularize
the singularity near that centroid. This is the compatibility default.

### Triangle P0

`field.element_kernel="triangle_p0"` distributes each element's total charge with constant density over its triangle.
Its analytic triangle kernel handles the near field and self potential, making it a different discretization from a centroid
point charge.

Triangle P0 requires `sim.softening=0`, finite non-degenerate triangles, and a resolved vacuum side for every element.
Current Phase 1 supports insulator surfaces only. Treecode evaluates near leaves with the analytic panel kernel and distant
nodes as monopoles for both field and potential. See [Direct](DirectSolver.en.html#triangle-p0),
[Treecode](Treecode.en.html), and [FMM](FMM.en.html#source-kernel) for details.

## Normalize length to control numerical scale

`sim.field_normalization` selects the representative internal length $L_0$. It does not change input or output units.

| Value | $L_0$ | Origin $\mathbf{x}_0$ |
| --- | --- | --- |
| `si` | 1 m | 0 |
| `length` | `field_length_scale` | 0 |
| `box` | Largest of the three box extents | `box_min` |
| `mesh` | Largest mesh-bounding-box extent | Mesh-bounding-box minimum |

The internal calculation uses

$$
\mathbf{x}'=\frac{\mathbf{x}-\mathbf{x}_0}{L_0},
\qquad
\epsilon'=\frac{\epsilon}{L_0},
$$

then multiplies electric field by $k_c/L_0^2$ and potential by $k_c/L_0$ to restore SI units. `box` requires
`use_box=true` and positive box extents. `mesh` falls back to `field_length_scale` only when the mesh is empty.

## Combine the solver with boundary components in periodic2

`periodic2` is not defined by one solver choice. It is a coupled configuration combining a finite or infinite image sum,
nonzero and zero modes, an outer sheath, reservoir acceleration or deceleration, and photoelectron escape or return.

The legacy `sim.field_bc_mode="periodic2"` path uses FMM, while small-system split-reference cases use a Direct panel spectral
backend. Their component compositions are described in
[periodic2 electrostatics](PeriodicElectrostatics.en.html) for the field decomposition and
[outer-plasma models](OuterPlasmaModels.en.html) for the outer response.

## Measure solver and source-discretization errors separately

For a new mesh or kernel, first compare a reduced identical case against Direct. Then vary the source-mesh resolution to
measure discretization error and vary solver controls to measure approximation error. Compare not only fields at a few points,
but also the observables used in the study, such as impact locations, accumulated charge, and conserved quantities. See
[Validating results](ValidationGuide.en.html) for the complete workflow.

## Code reference

- Initialization and automatic selection: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- Field and potential evaluation: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- Treecode/FMM charge refresh: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- Configuration values: [Configuration parameters](Parameters.en.html#field-solver)
