title: Direct solver

Lang: [日本語](DirectSolver.md) | [English](DirectSolver.en.md)

# Direct solver

The Direct solver evaluates every source-element contribution at every target. It evaluates the selected source discretization
directly, serving both as a solver for small problems and as the reference used to measure Treecode and FMM approximation error.
The central choice is whether each element charge is represented at its centroid or integrated over the triangle surface.

| Property | Description |
| --- | --- |
| Cost | $O(MN)$ for $N$ sources and $M$ targets |
| Kernels | point, triangle P0 |
| Normal field boundary | free |
| Geometry plan | None; evaluation reads the current `q_elem` |

```toml
[sim]
field_solver = "direct"
field_bc_mode = "free"
```

## Represent element charge at the centroid

With `field.element_kernel="point"`, the total charge $q_i$ of element $i$ is placed at its triangle centroid
$\mathbf{c}_i$. The field and potential at $\mathbf{r}$ are

$$
\mathbf{E}(\mathbf{r})=
k_c\sum_{i=1}^{N}q_i
\frac{\mathbf{r}-\mathbf{c}_i}
{\left(\lVert\mathbf{r}-\mathbf{c}_i\rVert^2+\epsilon^2\right)^{3/2}},
$$

$$
\phi(\mathbf{r})=
k_c\sum_{i=1}^{N}
\frac{q_i}{\sqrt{\lVert\mathbf{r}-\mathbf{c}_i\rVert^2+\epsilon^2}},
$$

where $\epsilon$ is `sim.softening`. Softening regularizes the point-charge singularity. Use triangle P0 to resolve charge over
the triangle surface; when using the point kernel, measure sensitivity to softening.

When `softening=0`, an ordinary point evaluation skips a contribution whose source centroid exactly equals the target to avoid
division by zero. This does not define a physical point-kernel self field at the centroid.

<a id="triangle-p0"></a>

## Integrate surface charge over each triangle

With `field.element_kernel="triangle_p0"`, a constant density $\sigma_i=q_i/A_i$ is placed on the triangle of area $A_i$:

$$
\phi_i(\mathbf{r})=
k_c\int_{T_i}\frac{\sigma_i}{\lVert\mathbf{r}-\mathbf{r}'\rVert}\,dS',
\qquad
\mathbf{E}_i(\mathbf{r})=
k_c\int_{T_i}\sigma_i
\frac{\mathbf{r}-\mathbf{r}'}{\lVert\mathbf{r}-\mathbf{r}'\rVert^3}\,dS'.
$$

BEACH evaluates these integrals with an analytic kernel based on edge logarithms and the solid angle. The same panel kernel is
used for far and near interactions; there is no centroid-point fallback.

Potential is continuous on the panel, and self potential uses its principal value. The normal field has the surface-charge jump

$$
\mathbf{E}^{\pm}=\mathbf{E}^{\mathrm{PV}}
\pm\frac{\sigma_i}{2\epsilon_0}\mathbf{n}_i.
$$

Therefore, each element needs a `normal_plus` or `normal_minus` vacuum side. Set `[mesh].surface_side` for OBJ input or
`surface_side` on every `[[mesh.templates]]` entry. `outward_closed` is valid only for consistently oriented closed
two-manifolds.

Current triangle P0 requirements are:

- `sim.softening = 0`
- finite triangles with positive area
- a resolved vacuum side for every element
- insulator surfaces only in Phase 1
- `direct`, `fmm`, or `auto`, which selects between those two

See [`examples/panel_direct.toml`](../examples/panel_direct.toml) for a complete configuration.

<a id="potential-at-element-centers"></a>

## Add self terms to element-center potentials

Element-center potential written to outputs such as `potential_history.csv` treats the self term differently from an arbitrary
point evaluation.

Triangle P0 uses the analytic panel integral including the source element. For the point kernel, BEACH adds a finite
representative self term. With $h_i=\sqrt{A_i}$, its coefficient is

$$
C_{ii}=\begin{cases}
1/\epsilon, & \epsilon>0,\\
2\sqrt{\pi}/h_i, & \epsilon=0,
\end{cases}
$$

and the self potential is $k_c C_{ii}q_i$. The zero-softening expression is an output self term for the point kernel, not the
exact self integral of a triangle P0 element. Preserve that kernel distinction when interpreting element-center potentials.

Computing the potential at every element center with Direct costs $O(N^2)$. A small history stride can therefore add substantial
cost separately from field evaluations at particle positions.

## Evaluate committed charge directly

Coordinates are normalized according to `field_normalization`, and fields and potentials are converted back to SI. See
[Field evaluation](FieldSolvers.en.html#length-normalization) for the equations and options.

Direct holds no geometry tree or expansion coefficients and reads the current `mesh%q_elem` at evaluation. The snapshot remains
fixed within a batch; deposited charge is committed at batch end and enters the field of the next batch, as for the other solvers.

## Separate errors using Direct as the reference

Direct has no solver approximation error, but it still contains errors from:

- the point-versus-triangle-P0 source discretization
- surface-mesh resolution and geometry
- particle time stepping and collision detection
- point-kernel softening
- periodic and outer-model assumptions when those models are present

Agreement with Direct therefore validates the Treecode or FMM approximation, not convergence of the full physical solution.
Refine the source mesh and particle time step independently.

## Code reference

- Direct field, potential, and self terms: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- Analytic triangle P0 kernel: [`bem_panel_kernel.f90`](../src/physics/panel/bem_panel_kernel.f90)
- Panel geometry: [`bem_panel_geometry.f90`](../src/physics/panel/bem_panel_geometry.f90)
- Configuration validation: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)
- Main regression tests: [`test_dynamics_field_solver.f90`](../tests/fortran/test_dynamics_field_solver.f90), [`test_panel_kernel.f90`](../tests/fortran/test_panel_kernel.f90)
