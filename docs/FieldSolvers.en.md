title: Field evaluation

Lang: [日本語](FieldSolvers.md) | [English](FieldSolvers.en.md)

# Field evaluation

At the start of each batch, BEACH constructs a field from the current element charges, `q_elem`, and uses that same field
for all particles tracked in the batch. Charge carried to a surface is committed at the end of the batch, so it changes the
field starting with the next batch.

Element charge is always represented as a constant density over its triangle. Users select how interactions are evaluated
and which field boundary is applied.

| Choice | Role | Values |
| --- | --- | --- |
| element kernel | Charge distribution on one element | `triangle_p0` (fixed) |
| solver | How interactions from many sources are accumulated | `direct` / `treecode` / `fmm` / `auto` |
| field boundary | How periodic images and the far field are included | `free` / `periodic2` |

## Select a solver from problem size

| Solver | Main use | Field boundary | Approximation |
| --- | --- | --- | --- |
| [Direct](DirectSolver.en.html) | Small problems, reference results, and split references | free, constrained periodic2 | Evaluates the analytic panel kernel for every element |
| [Treecode](Treecode.en.html) | Medium free-space problems | free | Uses monopoles for far nodes and the analytic panel kernel in near leaves |
| [FMM](FMM.en.html) | Large problems and many targets | free, periodic2 | Approximates far interactions with multipole and local expansions |
| `auto` | Select by element count for a free boundary | free | Direct or FMM |

With `auto`, Direct is used when `nelem < tree_min_nelem` and FMM is used at and above the threshold. The default threshold
is `256`. Runtime depends not only on the element count, but also on the
number of particles, number of steps, and target distribution, so measure a reduced case representative of the production run.

## Solver and Field-Boundary Compatibility

This table is the canonical compatibility reference for solvers and field boundaries.

| Solver | `free` | `periodic2` |
| --- | --- | --- |
| `direct` | Supported | Split reference only. Requires `periodic2.nonzero_mode_backend="panel_spectral_reference"`, `zero_mode_policy="exclude_k0"`, and a supported lower-boundary model |
| `treecode` | Supported | Unsupported |
| `fmm` | Supported | Supported. Infinite-periodic production runs use `cached_kneq0` |
| `auto` | Supported | Unsupported |

`periodic2` requires a `[domain]` box, `periodic_axes=["x", "y"]`, and a nonperiodic z axis. The Direct split reference is
for reduced validation cases; the normal periodic2 production path uses FMM.

## Separate domain topology from the field boundary

```toml
[domain]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0, 1.0, 1.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"
```

`[domain]` owns the simulation cell and periodic topology; `[field_boundary]` owns field closure. `mode="free"` selects a
free-space field, while `mode="periodic2"` uses the domain's x/y periodicity. The current periodic2 implementation accepts only
x/y periodic and z nonperiodic. Species boundary settings may select `open`, `reflect`, or `redistributed_reflect` on
nonperiodic faces, but cannot add, remove, or override periodic axes.

## Discretize element charge with triangle P0

BEACH distributes each element's total charge with constant density over its triangle. `triangle_p0` is the only element
kernel and is implicit rather than selected in a `[field]` table. It requires finite non-degenerate triangles and a resolved
vacuum side for every element. Treecode evaluates near leaves with the analytic panel kernel and distant
nodes as monopoles for both field and potential. See [Direct](DirectSolver.en.html#triangle-p0),
[Treecode](Treecode.en.html), and [FMM](FMM.en.html#source-kernel) for details.

The former `[field]` table and `sim.softening` key have been removed. Leaving either in an input fails as an unknown table or
key; BEACH does not reinterpret it as another element model.

## Normalize length to control numerical scale

`sim.field_normalization` selects the representative internal length $L_0$. It does not change input or output units.

| Value | $L_0$ | Origin $\mathbf{x}_0$ |
| --- | --- | --- |
| `si` | 1 m | 0 |
| `length` | `field_length_scale` | 0 |
| `box` | Largest extent in `[domain]` | `domain.box_min` |
| `mesh` | Largest mesh-bounding-box extent | Mesh-bounding-box minimum |

The internal calculation uses

$$
\mathbf{x}'=\frac{\mathbf{x}-\mathbf{x}_0}{L_0},
$$

then multiplies electric field by $k_c/L_0^2$ and potential by $k_c/L_0$ to restore SI units. `box` requires a `[domain]` with
positive extents. `mesh` falls back to `field_length_scale` only when the mesh is empty.

## Combine the solver with boundary components in periodic2

`periodic2` combines a finite- or infinite-image nonzero mode with a physical zero mode. Production uses FMM; reduced split
references use a Direct panel spectral backend. See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for component
ownership.

## Measure solver and mesh-discretization errors separately

For a new mesh, first compare a reduced identical case against Direct. Then vary the source-mesh resolution to
measure discretization error and vary solver controls to measure approximation error. Compare not only fields at a few points,
but also the observables used in the study, such as impact locations, accumulated charge, and conserved quantities. See
[Validating results](ValidationGuide.en.html) for the complete workflow.

## Code reference

- Initialization and automatic selection: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- Field and potential evaluation: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- Treecode/FMM charge refresh: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- Configuration values: [Configuration parameters](Parameters.en.html#field-solver)
