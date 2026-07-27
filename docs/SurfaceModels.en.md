title: Surface charge update

Lang: [日本語](SurfaceModels.md) | [English](SurfaceModels.en.md)

# Surface charge update

Charge absorbed by the surface is held as a difference during a batch. When a surface particle source is active, its reaction
charge joins the same difference. The differences are committed to `q_elem` at batch end, and surface-model processing after
commit produces the source charge for the next batch field.

## Update order within a batch

1. Build an immutable field snapshot from batch-start `q_elem`.
2. Add $q_pw_p$ from each particle's first mesh hit to thread-local differences.
3. Add surface-source reaction charge to `photo_emission_dq`.
4. Sum OpenMP-thread differences and MPI all-reduce global `dq`.
5. Apply `q_elem <- q_elem + dq`.
6. If conductors are present, equalize potential while conserving object charge.
7. Compute net pre/post-commit difference and the `tol_rel` metric.

Later particles in the same batch do not see charge from steps 2 or 3. The change enters the field at the next snapshot refresh.
This lag produces dependence on `batch_duration`; see [Batch duration and stability](BatchDurationStability.en.html).

## Conserved quantity and sign

`q_elem(i)` is total charge [C] on element $i$. Even with the implicit P0 panel discretization, stored state is not surface density; field
evaluation alone converts it to

$$
\sigma_i=\frac{q_i}{A_i}.
$$

When macro particle $p$ is absorbed on element $i$,

$$
\Delta q_i\mathrel{+}=q_pw_p.
$$

An electron deposits negative and a positive ion deposits positive charge. The particle is removed after absorption.

Ordered triangles for collision and `elem_vacuum_sign` for one-sided field traces are derived from the same mesh geometry. Surface
models do not rewrite triangle winding.

## Surface models

| `surface_model` | Post-commit treatment | Current role |
| --- | --- | --- |
| `insulator` | Retain charge on the hit element | Main v1.0 model |
| `conductor` | Equalize each `mesh_id` object while conserving total charge | `field_bc_mode="free"` only |
| `dielectric` | Retain charge and output `epsilon_r` metadata | Polarization not implemented |

## Insulator accumulation

`insulator` performs no post-commit redistribution. Charge changed by absorption or emission remains on that element.

The model does not include lateral surface conduction, bulk leakage, finite-resistivity relaxation, secondary-electron emission,
or specular/diffuse reflection. v1.0 interaction is primarily absorption and must not be interpreted as containing these effects.

## Floating conductor

Elements with `surface_model="conductor"` form one floating-conductor object per `mesh_id`. This is not a grounded fixed-potential
boundary. After particle differences are committed, potentials at element centroids are equalized while preserving each object's
total charge $Q_g^\mathrm{before}$.

Unknowns are all conductor charges $q_j$ and one scaled equipotential $V_g$ per group. For element $i$ in group $g(i)$,

$$
\sum_jA_{ij}q_j-V_{g(i)}=-\phi_i^\mathrm{fixed},
$$

The potential coefficient for source element $j$ carrying unit total charge as
a P0 triangle is

$$
A_{ij}=\frac{1}{A_j}\int_{T_j}
\frac{1}{|\mathbf c_i-\mathbf y|}\,dA_{\mathbf y}.
$$

$\mathbf c_i$ is the target-element centroid, while $A_j$ and $T_j$ are the
source-element area and triangle. The analytic P0 panel potential is used,
including its self term, under the principal-value side convention.
$\phi_i^\mathrm{fixed}$ is potential from nonconductor charge and uniform external field divided by `k_coulomb`.

Every group adds charge conservation

$$
\sum_{i\in g}q_i=Q_g^\mathrm{before}.
$$

The dense square system of size $N_\mathrm{cond}+N_\mathrm{group}$ is solved by Gaussian elimination with partial pivoting. Only
conductor $q_i$ values are replaced, so no charge moves between objects and nonconductor elements are unchanged.

This floating-conductor model uses centroid collocation with a P0 triangle
influence matrix. It cannot be combined with periodic or outer fields; current
validation accepts only `field_bc_mode="free"`. Check convergence of object
potential and charge distribution with mesh refinement.

## Dielectric metadata

`surface_model="dielectric"` and `epsilon_r` currently preserve geometry and material identity in output. BEACH does not solve
permittivity interface conditions, normal-$\mathbf D$ jumps, polarization charge, or internal field. Supplying `epsilon_r` does
not change field or charge update from the insulator behavior and must not be interpreted as dielectric polarization.

## OpenMP and MPI commit

The particle loop accumulates absorbed charge in `dq_thread(nelem,nthreads)`, avoiding an atomic update on every hit. After the
loop, thread columns are summed. Local differences plus `photo_emission_dq` are MPI all-reduced so all ranks hold identical global
`q_elem`.

Conductor relaxation runs deterministically on that same post-allreduce mesh state on every rank. A batch with an incomplete
collision or photo-ray query never reaches commit and does not use partial particle arrays or emission differences.

## `tol_rel`

With net difference after conductor relaxation

$$
\Delta\mathbf q=\mathbf q^\mathrm{after}-\mathbf q^\mathrm{before},
$$

the monitor is

$$
\mathrm{tol\_rel\ metric}
=\frac{\|\Delta\mathbf q\|_2}{\max(\|\mathbf q^\mathrm{after}\|_2,q_\mathrm{floor})}.
$$

Under the current contract, `tol_rel` is a monitoring and output metric, not an early-stop condition.

See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for reaction-charge signs and emitted-particle tracking.
See [Inspect Output Files](OutputGuide.en.html) for species-resolved charge balance, history, final output, and restart files.

## Code reference

- Particle absorption and batch commit: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Conductor relaxation facade: [`bem_surface_models.f90`](../src/physics/bem_surface_models.f90)
- Floating-conductor solver: [`bem_surface_models_conductor.f90`](../src/physics/bem_surface_models_conductor.f90)
- Batch statistics: [`bem_simulator_stats.f90`](../src/runtime/simulator/bem_simulator_stats.f90)
