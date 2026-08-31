title: Surface-charge update numerics

Lang: [日本語](SurfaceChargeNumerics.md) | [English](SurfaceChargeNumerics.en.md)

# Surface-charge update numerics

This numerical and implementation reference documents conserved surface state, batch-end ordering, the floating-conductor
linear system, and parallel reduction. See [how surfaces charge](SurfaceModels.en.html) if you only need to choose a model.

## Stored quantity and sign

`q_elem(i)` is total charge [C] on element $i$. Only field evaluation divides by element area $A_i$ to obtain P0 surface
charge density:

$$
\sigma_i=\frac{q_i}{A_i}.
$$

When macro particle $p$ is absorbed on element $i$, the pending change is

$$
\Delta q_i\mathrel{+}=q_p w_p.
$$

Electrons deposit negative charge and positive ions deposit positive charge; the absorbed particle is removed. Ordered
triangles for collision and one-sided field-trace signs come from the same mesh geometry. A surface model does not rewrite
triangle winding.

## Ordering at batch end

1. Build the field held fixed during particle tracking from batch-start `q_elem`.
2. Accumulate charge from each particle's first mesh hit in thread-local storage.
3. Add reaction charge from surface emission sources.
4. Apply global `neutral_return` reweighting when enabled.
5. Sum threads and reduce the global charge change across MPI ranks.
6. Add the change to `q_elem` exactly once.
7. If conductor objects exist, conserve their charge and equalize potential.
8. Calculate the net pre/post difference, statistics, and histories.

Particles in one batch share the field from stage 1. If an incomplete collision or photo-ray query invalidates a batch,
BEACH does not apply partial particle arrays or emission changes.

## Floating-conductor linear system

Elements with `surface_model="conductor"` form groups by `mesh_id`. For element $i$ in group $g(i)$, unknown element
charges $q_j$ and group potentials $V_g$ satisfy

$$
\sum_j A_{ij}q_j-V_{g(i)}=-\phi_i^\mathrm{fixed}.
$$

For source triangle $T_j$ carrying unit total charge, the P0 potential coefficient is

$$
A_{ij}=\frac{1}{A_j}\int_{T_j}
\frac{1}{|\mathbf c_i-\mathbf y|}\,dA_{\mathbf y}.
$$

Here $\mathbf c_i$ is the target-element centroid and $A_j$ is source-element area. The analytic panel potential includes
the self term under the principal-value side convention. $\phi_i^\mathrm{fixed}$ is potential from nonconductor charge and
the uniform external field divided by `k_coulomb`.

Every group also conserves total charge:

$$
\sum_{i\in g}q_i=Q_g^\mathrm{before}.
$$

The dense square system is solved by Gaussian elimination with partial pivoting, and only conductor elements are replaced.
Charge does not move between objects and nonconductor elements are unchanged. The current model uses centroid collocation
and P0-triangle influence and accepts only `field_boundary.mode="free"`.

## OpenMP, MPI, and restart

The particle loop gathers absorbed charge in thread-local arrays and sums them after the loop. Photoelectron-closure emission
and return quantities and local charge changes are reduced across MPI, so every rank holds the same global `q_elem`.
Conductor relaxation is then applied deterministically to that same state on every rank.

Check [run and resume](Execution.en.html) and [input parameters](Parameters.en.html) for checkpoint compatibility and
required files for restoring committed `q_elem` and related state.

## Definition of `tol_rel`

Including conductor relaxation, define the pre/post update difference as

$$
\Delta\mathbf q=\mathbf q^\mathrm{after}-\mathbf q^\mathrm{before}.
$$

The monitor is

$$
\mathrm{tol\_rel\ metric}
=\frac{\|\Delta\mathbf q\|_2}{\max(\|\mathbf q^\mathrm{after}\|_2,q_\mathrm{floor})}.
$$

`tol_rel` is an output metric, not an early-stop condition in the current implementation.

## Implementation and tests

- Particle absorption and batch commit: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Surface-model facade: [`bem_surface_models.f90`](../src/physics/bem_surface_models.f90)
- Floating-conductor solver: [`bem_surface_models_conductor.f90`](../src/physics/bem_surface_models_conductor.f90)
- Batch statistics: [`bem_simulator_stats.f90`](../src/runtime/simulator/bem_simulator_stats.f90)
- Model regression: [`test_surface_models.f90`](../tests/fortran/test_surface_models.f90)
