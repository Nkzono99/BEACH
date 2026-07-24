title: Particle collision and boundary events

Lang: [日本語](ParticleEvents.md) | [English](ParticleEvents.en.md)

# Particle collision and boundary events

BEACH accepts only the first collision or boundary crossing on the current trajectory segment. A triangle hit absorbs the
particle at that position. A reflect or periodic crossing starts a newly integrated trajectory for the remaining time. Repeating
this operation completes one particle step.

## Accept only the first collision or crossing

Given current state $(\mathbf{x}_0,\mathbf{v}_0)$ and candidate $(\mathbf{x}_1,\mathbf{v}_1)$ from the
[Boris particle update](BorisPusher.en.html), BEACH proceeds as follows:

1. Query the first mesh hit on $\mathbf{x}(t)=\mathbf{x}_0+t(\mathbf{x}_1-\mathbf{x}_0)$.
2. If the endpoint is inside the box, commit either the mesh hit or the endpoint.
3. If it is outside, compare the first box-face fraction against the mesh-hit $t$.
4. Absorb when the mesh is simultaneous or earlier; otherwise apply the box-face action.
5. If reflect or periodic keeps the particle alive, rebuild a candidate for the remaining time.

If mesh-hit and box-face fractions differ by at most $64\epsilon_\mathrm{mach}\max(1,|t|)$, the mesh is treated as first.
Even if one outer step contains several events, BEACH always advances to the first event of the current segment before
re-integrating the remainder.

| Earliest event | Committed action |
| --- | --- |
| mesh | Return hit position and element index; absorb the particle |
| open face | Escape at the event position or pass to an outer interface |
| reflect face | Reverse normal velocity and advance the remainder from just inside the box |
| periodic face | Move just inside the opposite face, retain velocity, and advance the remainder |

## Intersect the trajectory segment with triangles

### Reduce the triangle candidates

Mesh initialization builds an axis-aligned bounding box (AABB) for every triangle.

| Element count | Candidate search |
| ---: | --- |
| below 64 | Linearly inspect every element AABB |
| 64 and above | Use a uniform grid and 3D DDA to inspect only cells crossed by the trajectory segment |

The uniform grid targets eight elements per cell and caps each axis at 128 cells. Triangle indices are stored in CSR form for
every cell overlapped by the triangle AABB. These are current fixed implementation values, not input parameters.

At query time, BEACH intersects the trajectory segment with the grid AABB and visits crossed cells with 3D DDA. A triangle may be registered
in more than one cell, but only the smallest intersection parameter is retained. The DDA iteration bound is
`nx + ny + nz + 3`. If cell indices, increments, or parameters do not make bounded finite progress, the result is not treated as
"no hit"; the query returns `collision_query_grid_stalled`.

### Select the first intersection

Each candidate triangle uses Möller–Trumbore intersection. Writing the triangle as

$$
\mathbf{r}(u,v)=\mathbf{v}_0
+u(\mathbf{v}_1-\mathbf{v}_0)
+v(\mathbf{v}_2-\mathbf{v}_0),
$$

an intersection is accepted only when

$$
0\le u\le1,
\qquad
0\le v,
\qquad
u+v\le1,
\qquad
0\le t\le1.
$$

If the determinant magnitude is at most $64\epsilon_\mathrm{mach}$ times the product of segment length and the two edge lengths,
the geometry is considered degenerate or nearly parallel. The determinant sign is not culled, so collisions are detected from
both sides of a triangle. Triangle winding does not determine collision sidedness; it is used separately for quantities such as
the field vacuum-side trace.

The smallest $t$ is selected among multiple triangle hits, with
$\mathbf{h}=\mathbf{x}_0+t(\mathbf{x}_1-\mathbf{x}_0)$.

## Map periodic-image collisions to the primary cell

With periodic2, the mesh still stores only base elements in the primary cell. BEACH enumerates only image shifts whose canonical
mesh AABB can overlap the trajectory-segment AABB on the two periodic axes. For period $L$, segment range $[p_{\min},p_{\max}]$, and mesh range
$[m_{\min},m_{\max}]$ on one axis,

$$
n_{\min}=\left\lceil
\frac{p_{\min}-m_{\max}-\mathrm{tol}}{L}
\right\rceil,
\qquad
n_{\max}=\left\lfloor
\frac{p_{\max}-m_{\min}+\mathrm{tol}}{L}
\right\rfloor.
$$

For each image, the trajectory segment is shifted by $-nL$ into the base-mesh frame and passed to the ordinary intersection query. A hit stores:

| Value | Meaning |
| --- | --- |
| `hit%pos` | Physical coordinate on the periodic image that was hit |
| `hit%pos_wrapped` | Coordinate wrapped into the primary cell |
| `hit%image_shift` | Image indices on the two periodic axes |
| `hit%elem_idx` | Base-mesh element index |

If candidate $t$ values agree within a relative tolerance of `1e-12`, selection is deterministic by element index, first image
index, then second image index. If the image count on one axis or its Cartesian product exceeds 4,096, the query does not skip
work and report no hit; it returns `collision_query_image_limit`.

This image range is determined by which mesh images a particle trajectory segment can hit. It is unrelated to
`field_periodic_image_layers`, which controls the field sum. See
[periodic2 electrostatics](PeriodicElectrostatics.en.html) for field images.

## Apply the action of each box face

With `use_box=true`, every low and high axis face is configured as `open`, `reflect`, or `periodic`. When a trajectory segment reaches several
faces simultaneously at an edge or corner, faces within a machine-epsilon tolerance of the minimum fraction are combined into
one mask.

### Open, reflect, and periodic faces

With `external_boundary.ordinary_open.model="escape"`, any open face in the mask removes the particle and increments `escaped_boundary`.
A reflect-only face reverses the corresponding velocity component. A periodic face moves the particle to the opposite face
without changing velocity. A surviving event uses `nearest` so its new coordinate is one floating-point value inside the box,
not exactly on the face.

Reflect and periodic actions at a corner are applied from one face mask, making the result independent of axis traversal order.

### Potential barrier

`external_boundary.ordinary_open.model="potential_barrier"` compares outward normal kinetic energy at a single open face,

$$
K_n=\frac{1}{2}m v_\mathrm{out}^2,
$$

against

$$
\Delta U=q\left(\phi_\infty-\phi_\mathrm{boundary}\right).
$$

The particle reflects when $v_\mathrm{out}>0$, $\Delta U>0$, and $K_n<\Delta U$; otherwise it escapes. This is a reduced
single-face model. A corner involving multiple simultaneous open faces is not generalized and returns
`particle_step_ambiguous_open_corner`.

Reservoir acceleration, the outer sheath, and photoelectron return are physical boundary models built around this event
mechanism. They are documented separately in [Particle sources and boundaries](ParticleSourcesBoundaries.en.html) and
[Outer-plasma models](OuterPlasmaModels.en.html).

## Advance the time remaining after a boundary crossing

Position and velocity at a box event are interpolated between candidate input and output states, and only the active face
coordinate is set exactly to the box value. After applying the surviving action, BEACH uses
a guard derived from the box coordinates and span to place reflected or periodic particles inside the face. This avoids a
subnormal one-ULP offset at a zero-valued face and prevents the next event fraction from underflowing to zero.

$$
\Delta t_\mathrm{remain}=(1-t_\mathrm{event})\Delta t_\mathrm{segment}
$$

to build a new predicted-midpoint field and Boris candidate. The earliest mesh or box event is then queried on the remainder.
The pre-reflection field sample is not reused for that remaining time.

At most eight box events are processed in one local continuation. All eight are supported; if a ninth is required,
`particle_step_multiple_box_events` is returned and the incomplete state is not committed. This also acts as a safety signal for
an excessively large `dt`, narrow box, or very fast particle.
The scale-aware guard does not change this existing eight-event safety limit.

The default `multiple_box_events_policy="abort"` fails the run closed at this point. Only an explicit
`"soft_discard"` removes the affected macro-particle and records batch, rank, particle, species, step, macro charge,
position, and velocity on standard error. Count and absolute macro charge are also retained in `summary.txt`, restart,
and the charge ledger. The run aborts when either configured cumulative limit is exceeded. This is a bounded numerical
workaround for qualitative comparisons, not a replacement for a physical boundary model.

## Transfer particles through z-high to an outer model

When particle transfer is active, the z-high open face is returned to the caller as an interface crossing instead of immediately
applying escape. For this face only, crossing time is refined from a quadratic trajectory consistent with the Boris endpoint
positions and velocities rather than using the straight-segment fraction unchanged. The payload contains face, fraction of the full step,
position, velocity, and remaining time.

For a step whose candidate endpoint is outside z-high, BEACH evaluates candidate x/y crossing times on the same quadratic
trajectory and selects the event that actually occurs first. An earlier lateral periodic or reflecting face is applied before
reintegrating the remainder; an earlier z-high event is transferred to the outer model. Simultaneous lateral actions are composed
before transfer.

This refinement covers crossings detected because the candidate endpoint is outside the box. It does not search for a temporary
excursion that returns to an inside endpoint; mesh hits retain the existing chord test.

If the outer model returns the particle locally, ordinary particle stepping resumes from the returned position and velocity for
the remaining time. If it escapes to infinity, it is removed. See [Outer-plasma models](OuterPlasmaModels.en.html) for outer
acceleration and return conditions.
Across the original local step, BEACH processes at most eight external events; a ninth returns
`particle_step_multiple_external_events`. This budget is distinct from the box-event budget of each local continuation.

## Stop when the query cannot be completed

A collision query is `ok` only when all required candidates were examined.

| Status | Code | Meaning |
| --- | ---: | --- |
| `collision_query_ok` | 0 | Query completed |
| `collision_query_image_limit` | 1 | Periodic image enumeration exceeded 4,096 |
| `collision_query_index_range` | 2 | Image bounds were non-finite or outside integer range |
| `collision_query_invalid_segment` | 3 | A trajectory-segment endpoint was non-finite |
| `collision_query_grid_stalled` | 4 | Invalid grid geometry or DDA failed to progress |
| `particle_step_invalid_boundary` | 1001 | Invalid particle, box, or event geometry |
| `particle_step_multiple_box_events` | 1002 | A ninth box event was needed in one step |
| `particle_step_ambiguous_open_corner` | 1003 | Multiple open faces occurred at an outer-owned or potential-barrier event |
| `particle_step_multiple_external_events` | 1004 | A ninth external event was needed in the original local step |

The former name `particle_step_unsupported_barrier_corner` remains as a compatibility alias for code 1003.

Treating these states as "no hit" could let particles pass through a surface and corrupt the charge ledger, so normal tracking
fails closed. Within OpenMP, the earliest particle and step failure is selected. In MPI, failure rank and particle state are
shared so every rank reports the same batch, rank, particle, step, and status before stopping. Photo raycasts apply the same rule
to species, ray, and bounce.

To inspect the internal cause of `grid_stalled`, set `BEACH_COLLISION_DIAGNOSTICS=1`. On the failing path, BEACH then writes
the DDA branch name, `p0` / `p1`, collision-grid bounds, cell indices, and values such as `t_cur`, `t_next`, and `t_delta` to
standard error. This environment variable enables diagnostics only and does not change collision physics.

## Converge collision positions and charging results

1. Halve `dt` and verify stability of the segment-hit element and position.
2. Build reduced cases with thin surfaces, edge and corner hits, and nearly parallel trajectory segments.
3. Test a mesh hit during the remainder after reflection and after periodic wrapping.
4. At a periodic seam, verify that `pos` and `pos_wrapped` identify the intended base element.
5. If `multiple_box_events` occurs, reduce `dt` first. If a qualitative finite-image comparison uses soft discard,
   verify that its rate and absolute charge cannot affect the conclusion.
6. Refine the mesh and check convergence of the first hit and final surface-charge distribution.

## Code reference

- Box/mesh event ordering and remainder: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- Triangles, grid DDA, and periodic images: [`bem_collision.f90`](../src/physics/bem_collision.f90)
- Box-face detection and actions: [`bem_boundary.f90`](../src/physics/bem_boundary.f90)
- Collision-grid construction: [`bem_mesh.f90`](../src/mesh/bem_mesh.f90)
- Event ordering and eight-event-bound tests: [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)
- Boundary unit tests: [`test_boundary.f90`](../tests/fortran/test_boundary.f90)
- Collision, periodic, and fail-closed tests: [`test_dynamics_basic.f90`](../tests/fortran/test_dynamics_basic.f90)
