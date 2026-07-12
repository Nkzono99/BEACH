# Kinetic Sheath Particle Coupling Design

- Date: 2026-07-12
- Status: approved
- Related decision: `docs/adr/0001-kinetic-outer-plasma.md`

## Objective

Connect the converged `outer_plasma.model="kinetic_1d"` state to particle
inflow and outward transfer so that the field, boundary VDF, escape decision,
turning point, and return time all use the same potential profile referenced to
`outer_plasma.infinity_potential`.

The supported physics remains electrostatic, collisionless, unmagnetized,
one-dimensional, and monotonic. Non-monotonic virtual cathodes, trapped
populations, magnetic outer orbits, and persistent outer queues remain outside
this change and fail closed.

## Existing Defect

The field snapshot already solves and broadcasts the kinetic profile, including
`interface_potential` and `infinity_potential`. Particle injection does not own
that state. With `electrostatic_1d_instant_return`, it reconstructs an interface
barrier from `lambda_D * Q_surface / (epsilon_0 A)`, which is a linear-Debye
approximation and need not equal the nonlinear kinetic solution.

Outward transfer uses the solved endpoint potential for escape, but computes
return time with an analytic exponential Debye profile. The configuration also
forbids `photoelectron_closure="kinetic_mean"` together with tracked-particle
return even though ADR 0001 assigns mean closure only to outer space charge and
tracked particles to the surface charge ledger.

## State Ownership

`electrostatic_snapshot_type%outer` is the only runtime source of the active
outer potential. Batch injection receives a read-only outer state after the
snapshot has been refreshed for that batch. No injection routine recomputes a
sheath potential from total surface charge or a free-space face average when
`kinetic_1d` coupling is active.

The kinetic state must be ready and MPI-identical before injection. A missing,
stale, non-finite, or incompatible state is a fatal model error; there is no
fallback to `reservoir_potential_model`, `sheath_injection_model`, or a linear
Debye estimate.

## Inflow Mapping

The configured negative and positive `z_high` `reservoir_face` species remain
infinity VDFs. For each species, define

```text
delta_phi = phi_interface - phi_infinity
v_interface^2 = v_infinity^2 - 2 q delta_phi / m.
```

Only infinity normal velocities with a non-negative interface radicand are
accessible. Macro-particle flux and residual counting use the accessible
infinity flux. Sampled normal velocities are then mapped to the interface by
the same energy equation. Tangential velocity is unchanged. This mapping is
applied every batch from the refreshed snapshot, including after restart.

## Outward Mapping

An outward particle crossing `z_high` is evaluated against the same kinetic
profile. Conservation of normal energy gives

```text
K_normal(z) = K_normal(interface) + q [phi_interface - phi(z)].
```

If `K_normal(infinity) >= 0`, the particle escapes. Otherwise the first root of
`K_normal(z)=0` on the monotonic profile is the turning point. The return time
is twice the outward integral of `dz / v_normal(z)`. The implementation treats
the stored potential as piecewise linear and integrates each segment
analytically, including the integrable square-root endpoint at the turning
point. It rejects non-monotonic, non-finite, or sign-inconsistent profiles.

The returned particle keeps its tangential velocity, advances tangentially by
`v_t * tau_outer`, wraps x/y periodically, reverses its normal velocity at the
local interface, and resumes the remaining local step. Existing frozen-field
and queue limits remain mandatory.

## Photoelectron Contract

`photoelectron_closure="kinetic_mean"` supplies only the mean outer charge
density and current diagnostics used by the kinetic Poisson solve. It does not
deposit return charge and does not inject a second statistical particle.

Tracked `photo_raycast` particles may therefore use the same kinetic-profile
return mapper. Surface charge changes only through the existing tracked
emission and later absorption events. The outgoing histogram remains an audit
diagnostic and must not create charge. `deposit_opposite_charge_on_emit=true`
continues to be required for charge-conserving tracked photoelectrons.

## Configuration Contract

For `kinetic_1d` particle transfer:

- `coupling.particle_transfer_mode="electrostatic_1d_instant_return"`
- `outer_plasma.return_model="kinetic_1d_profile_return"`
- open `z_high`, periodic x/y, and `b0=0`
- positive `field_evolution_timescale` and `max_frozen_field_ratio`
- `outer_queue_enabled=false`
- `sim.reservoir_potential_model="none"`
- `sim.sheath_injection_model="none"`

The existing `electrostatic_1d_instant_return` return-model identifier remains
available only for `linear_debye`. Using that identifier with `kinetic_1d` is
rejected so an exponential flight-time approximation cannot be selected by
accident.

## Restart And MPI

The existing checkpointed outer profile remains the restart source. A resumed
run refreshes or restores the snapshot before the next batch injection. Root
solves and broadcasts the kinetic state; all ranks derive the same accessibility
threshold and effective global reservoir count before distributing particles.
No new rank-local physical state is introduced.

The return-model identifier participates in the existing model fingerprint.
Changing from a linear return model to `kinetic_1d_profile_return` is therefore
restart-incompatible unless the checkpoint already declares that model.

## Verification

Tests must demonstrate:

1. Inflow kinetic energy changes by `-q delta_phi` and inaccessible infinity
   velocities are excluded.
2. The injection barrier uses the snapshot interface potential, not total
   surface charge or `debye_length`.
3. Piecewise-linear profile integration agrees with constant-field analytic
   flight time and converges for a smooth monotonic profile.
4. Escape and return conserve normal energy; x/y displacement and wrapping use
   the computed return time.
5. MPI-global reservoir counts and sampled accessibility are rank-count
   invariant.
6. Restart produces the same first resumed batch as an uninterrupted run.
7. `kinetic_mean` plus tracked return changes the surface ledger only through
   tracked emission/absorption.
8. Linear Debye, Zhao, and legacy barrier configurations retain their existing
   behavior, while mixed or unsupported kinetic configurations fail closed.

Heavy MPI and full FMM verification remains opt-in under the repository test
tiers and runs only on a KUDPC compute node.
