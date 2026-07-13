title: Boris particle update

Lang: [日本語](BorisPusher.md) | [English](BorisPusher.en.md)

# Boris particle update

BEACH advances one particle step in two stages: it constructs a candidate trajectory, then resolves the first collision or
boundary crossing on that trajectory. The candidate combines a predicted-midpoint field sample, a Boris velocity update, and
a trapezoidal position update. Public state and checkpoints keep $(\mathbf{x},\mathbf{v})$ at the same time level before and
after the step.

| State | Time level |
| --- | --- |
| Input | $\mathbf{x}^n,\mathbf{v}^n$ |
| Field sample | Predicted position $\mathbf{x}_\mathrm{mid}$ |
| Output candidate | $\mathbf{x}^{n+1},\mathbf{v}^{n+1}$ |
| Magnetic field | Uniform `sim.b0` |

## Sample the field at a predicted midpoint

For step interval $\Delta t$, current velocity predicts the midpoint

$$
\mathbf{x}_\mathrm{mid}=
\mathbf{x}^n+\frac{1}{2}\mathbf{v}^n\Delta t.
$$

The field snapshot is evaluated once there:

$$
\mathbf{E}_\mathrm{mid}=
\mathbf{E}_\mathrm{snapshot}(\mathbf{x}_\mathrm{mid}).
$$

The snapshot composes the element-charge field, `sim.e0`, and any selected periodic zero mode or outer profile exactly once.
Particles in one batch share the snapshot formed from element charge at the start of that batch. See
[Field evaluation](FieldSolvers.en.html) for its construction.

With `use_box=true`, only the field-sample position is mapped into the solver's valid region.

| Axis | Treatment of $\mathbf{x}_\mathrm{mid}$ |
| --- | --- |
| Both low and high faces periodic | Wrap into the primary box with `modulo` |
| Otherwise | Clamp between `box_min` and `box_max` |

Candidate trajectory endpoints remain in physical coordinates so that the next stage can order a triangle collision against a
box-boundary crossing.

## Advance velocity with Boris

For particle charge $q$, mass $m$, and $\mathbf{B}=\texttt{sim.b0}$, first apply a half electric acceleration:

$$
\mathbf{v}^-=
\mathbf{v}^n+\frac{q}{m}\mathbf{E}_\mathrm{mid}\frac{\Delta t}{2}.
$$

The magnetic rotation is

$$
\mathbf{t}=\frac{q}{m}\mathbf{B}\frac{\Delta t}{2},
\qquad
\mathbf{s}=\frac{2\mathbf{t}}{1+\lVert\mathbf{t}\rVert^2},
$$

$$
\mathbf{v}'=\mathbf{v}^-+\mathbf{v}^-\times\mathbf{t},
\qquad
\mathbf{v}^+=\mathbf{v}^-+\mathbf{v}'\times\mathbf{s},
$$

followed by the other half electric acceleration:

$$
\mathbf{v}^{n+1}=
\mathbf{v}^++\frac{q}{m}\mathbf{E}_\mathrm{mid}\frac{\Delta t}{2}.
$$

The implementation uses the standard Boris rotation with $t^2=\mathbf{t}\cdot\mathbf{t}$.

## Form the candidate position with the trapezoidal rule

After velocity is updated, the candidate position uses same-time input and output velocities in the trapezoidal rule:

$$
\mathbf{x}^{n+1}=
\mathbf{x}^n+\frac{1}{2}
\left(\mathbf{v}^n+\mathbf{v}^{n+1}\right)\Delta t.
$$

For a uniform electric field alone, this matches the constant-acceleration velocity and displacement. Regression tests for a
smooth spatial electrostatic field verify second-order convergence of candidate position and velocity when combined with the
predicted-midpoint field sample.

## Split the step at collisions and boundary crossings

The equations above produce a full-step candidate $(\mathbf{x}^{n+1},\mathbf{v}^{n+1})$. It becomes the accepted state when
the candidate trajectory has no triangle collision or box-boundary crossing.

Triangle collisions and ordinary box crossings are ordered along the **straight trajectory segment within the step** that
connects $\mathbf{x}^n$ and $\mathbf{x}^{n+1}$. A triangle collision absorbs the particle at the hit point. If a reflect or
periodic face is crossed first, BEACH interpolates position and velocity to the crossing time and constructs a new Boris
candidate for the remaining interval. At the z-high outer interface, the applicable coupling path refines crossing time from a
quadratic trajectory consistent with the Boris endpoints. See
[Particle collision and boundary events](ParticleEvents.en.html).

## Properties of the Boris update

| Condition | Property |
| --- | --- |
| Pure B | Preserves $\lVert\mathbf{v}\rVert$ to roundoff |
| Prescribed uniform E/B | Velocity update is reversible under a sign change of $\Delta t$ |
| Smooth prescribed field | Current same-time update is second order in position and velocity |
| Mesh collision or open-boundary escape | Irreversible process that ends tracking at the collision or boundary point |
| Self-consistent field changing between batches | Snapshot fixed within a batch and refreshed after charge commit |

## Choose `dt` from trajectories and charging results

`dt` must resolve particle rotation, spatial field variation, and collision geometry:

- cyclotron angle $|q|\lVert\mathbf{B}\rVert\Delta t/m$
- transit time over the spatial electric-field variation scale
- transit time near triangles and box faces
- transit time through a rapidly changing outer-interface profile

Triangle intersections use the straight trajectory segment within a step. When an orbit bends strongly in one step or passes
fine geometry, reduce the time step until collision positions converge. Compare `dt`, `dt/2`, and when needed `dt/4` using
trajectories, hit elements, absorption and escape statistics, and post-batch element charge. Keeping `max_step` fixed halves the
maximum physical tracking time when `dt` is halved, so adjust `max_step` when comparing the same time horizon.

## Code reference

- Boris velocity rotation and trapezoidal position: [`bem_pusher.f90`](../src/physics/bem_pusher.f90)
- Predicted-midpoint field and remainder reintegration: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- E/B, pure-B conservation, reversibility, and convergence tests: [`test_dynamics_basic.f90`](../tests/fortran/test_dynamics_basic.f90)
- Charged-mesh and field-sample tests: [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)
