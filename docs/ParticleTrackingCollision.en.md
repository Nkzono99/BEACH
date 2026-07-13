title: Particle tracking and collision

Lang: [English](ParticleTrackingCollision.en.md) | [日本語](ParticleTrackingCollision.md)

# Particle tracking and collision

BEACH stores position and velocity at the same time and treats the Boris velocity update, trapezoidal position
update, and box/mesh event selection as one step.

## Same-time Boris update

For input $(\mathbf{x}^n,\mathbf{v}^n)$, the electric field is evaluated once at a predicted midpoint.

$$
\mathbf{x}_\mathrm{mid}=\mathbf{x}^n+\frac12\mathbf{v}^n\Delta t,
\qquad
\mathbf{E}_\mathrm{mid}=\mathbf{E}(\mathbf{x}_\mathrm{mid})+\mathbf{E}_0
$$

A Boris rotation with this field and uniform `sim.b0` gives $\mathbf{v}^{n+1}$. Position uses a trapezoidal
update.

$$
\mathbf{x}^{n+1}=\mathbf{x}^n+
\frac12(\mathbf{v}^n+\mathbf{v}^{n+1})\Delta t
$$

Public state does not store half-step velocity. The contract is second-order position and velocity in a smooth
prescribed field, and speed-norm preservation to roundoff in a pure magnetic field.

## Event ordering

The earliest triangle collision or box crossing is selected along the candidate motion.

| First event | Behavior |
| --- | --- |
| Mesh hit | Absorb at the surface; discard the candidate endpoint |
| Open face | Escape or transfer to an outer model |
| Reflecting face | Reverse normal velocity and reintegrate remaining time |
| Periodic face | Wrap and reintegrate remaining time |

A safety limit bounds box events in one step. Event processing that cannot make finite progress fails closed.

## Triangle collision

Small meshes use a linear element scan. Larger meshes use a uniform AABB grid and 3D-DDA traversal to select
candidates. Möller–Trumbore intersection is applied to candidate triangles, and the smallest segment parameter
defines the first hit.

Degenerate triangles, nonfinite coordinates, stalled DDA state, and excessive periodic-image ranges do not
continue as “no hit”; they return status and stop coherently.

## periodic2 collision

The mesh stores primary-cell elements only. BEACH computes the periodic images that can intersect a particle
segment and shifts the segment back to the base mesh for each test. Both physical-image and wrapped hit positions
are retained.

Field images and collision images serve different purposes. See
[periodic2 field evaluation](PeriodicElectrostatics.en.html) for the field treatment.

## Unresolved particles

A particle that reaches `sim.max_step` without absorption or escape is counted as `survived_max_step`. It is not
silently reclassified as absorbed or escaped.

## Numerical checks

- Halve `sim.dt` and compare trajectory, hit position, and outcomes.
- Refine the mesh and check first-hit stability.
- Verify that `survived_max_step` does not affect the conclusion.
- Test trajectories crossing periodic seams.

Full equations and status tables remain in the legacy combined
[Particle tracking and surface charge accumulation](ParticleChargeLoop.en.html#8-particle-time-integration-boris-velocity-update-and-same-time-state).
