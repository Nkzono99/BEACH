title: Surface charge update

Lang: [日本語](SurfaceModels.md) | [English](SurfaceModels.en.md)

# Surface charge update

Charge absorbed by the surface and reaction charge from surface emission are held as differences during a batch and committed to
`q_elem` together at batch end. Surface-model processing after commit produces the source charge for the next batch field.

## Update order within a batch

1. Build an immutable field snapshot from batch-start `q_elem`.
2. Add $q_pw_p$ from each particle's first mesh hit to thread-local differences.
3. Add photoelectron source reaction charge to `photo_emission_dq`.
4. Sum OpenMP-thread differences and MPI all-reduce global `dq`.
5. Apply `q_elem <- q_elem + dq`.
6. If conductors are present, equalize potential while conserving object charge.
7. Compute net pre/post-commit difference and the `tol_rel` metric.

Later particles in the same batch do not see charge from steps 2 or 3. The change enters the field at the next snapshot refresh.
This lag produces dependence on `batch_duration`; see [Batch duration and stability](BatchDurationStability.en.html).

## Conserved quantity and sign

`q_elem(i)` is total charge [C] on element $i$. Even for the `triangle_p0` kernel, stored state is not surface density; field
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

## Photoelectron emission

With `photo_raycast` and `deposit_opposite_charge_on_emit=true`, source element $i$ receives

$$
\Delta q_{i,\mathrm{emit}}=-q_pw_p.
$$

This is accumulated separately in `photo_emission_dq` and merged into the same commit as later collision deposition.

If the emitted particle returns to element $j$, ordinary absorption adds $+q_pw_p$ there. Return to the same element cancels source
charge; return elsewhere transfers net charge between surface elements. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for ray weight, reduced escape, and tracked outer return.

## Floating conductor

Elements with `surface_model="conductor"` form one floating-conductor object per `mesh_id`. This is not a grounded fixed-potential
boundary. After particle differences are committed, potentials at element centroids are equalized while preserving each object's
total charge $Q_g^\mathrm{before}$.

Unknowns are all conductor charges $q_j$ and one scaled equipotential $V_g$ per group. For element $i$ in group $g(i)$,

$$
\sum_jA_{ij}q_j-V_{g(i)}=-\phi_i^\mathrm{fixed},
$$

with

$$
A_{ij}=
\begin{cases}
1/\epsilon, & i=j,\ \epsilon>0,\\
2\sqrt\pi/h_i, & i=j,\ \epsilon=0,\\
1/\sqrt{|\mathbf c_i-\mathbf c_j|^2+\epsilon^2}, & i\ne j.
\end{cases}
$$

$\mathbf c_i$ is element centroid, $h_i$ element length scale, and $\epsilon$ field softening.
$\phi_i^\mathrm{fixed}$ is potential from nonconductor charge and uniform external field divided by `k_coulomb`.

Every group adds charge conservation

$$
\sum_{i\in g}q_i=Q_g^\mathrm{before}.
$$

The dense square system of size $N_\mathrm{cond}+N_\mathrm{group}$ is solved by Gaussian elimination with partial pivoting. Only
conductor $q_i$ values are replaced, so no charge moves between objects and nonconductor elements are unchanged.

This is a simplified centroid point-potential conductor, not an exact triangle-P0 conductor BEM boundary integral. It cannot be
combined with periodic or outer fields; current validation accepts only `field_bc_mode="free"`. Check convergence of object
potential and charge distribution with mesh refinement and softening.

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

## Charge ledger

`charge_ledger.csv` separates species-resolved signed charge and counts.

| Flux | Meaning |
| --- | --- |
| `injected_from_remote` | Charge entering from `volume_seed` or `reservoir_face` |
| `emitted_from_surface` | Tracked charge leaving through `photo_raycast` |
| `absorbed_on_surface` | Charge absorbed by the mesh |
| `escaped_to_infinity` | Charge escaping through open or outer models |
| `discarded_unresolved` | Charge discarded alive at `max_step` |
| `interface_outward_gross` / `returned_gross` | Gross transfer across the outer ownership interface |

The transactional residual combines before/after changes of surface, local-flight, outer-flight, and unresolved stocks with remote
injection, escape, and discard. Surface emission and absorption are internal transfers between surface and flight stocks and are
not counted twice as independent external sources.

Separately from a residual where opposite species can cancel, `discarded_unresolved_abs` checks
$\sum_s|Q_{s,\mathrm{discard}}|$. A small residual does not validate a run with large unresolved discard.

## History, final output, and restart

With `history_stride>0`, `charge_history.csv` is written when

$$
(\texttt{stats.batches}-1)\bmod\texttt{history\_stride}=0,
$$

so batch 1 is always included. With `write_potential_history=true`, the current `q_elem` refreshes the field snapshot at the same
stride and writes element-centroid `potential_history.csv`.

Main files are `charges.csv`, `charge_history.csv`, `charge_ledger.csv`, `summary.txt`, and `mesh_sources.csv`. Restart validates
mesh element count, MPI world size, finite statistics and charge, ordered mesh/species/model fingerprints, ledger stocks, and any
required outer profile. Missing required checkpoints with `output.resume=true` stop rather than start a new run.

## Code reference

- Particle absorption and batch commit: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Conductor relaxation: [`bem_surface_models.f90`](../src/physics/bem_surface_models.f90)
- Transactional charge ledger: [`bem_charge_ledger.f90`](../src/runtime/coupling/bem_charge_ledger.f90)
- Batch statistics: [`bem_simulator_stats.f90`](../src/runtime/simulator/bem_simulator_stats.f90)
- History output: [`bem_simulator_io.f90`](../src/runtime/simulator/bem_simulator_io.f90)
- Final output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- Restart validation: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
