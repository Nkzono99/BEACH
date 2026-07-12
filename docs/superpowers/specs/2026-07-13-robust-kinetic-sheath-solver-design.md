# Robust Kinetic Sheath Solver Design

## Purpose

The `kinetic_1d` outer-plasma model currently stops when a finite-difference
Newton step cannot find a residual-reducing point on the monotonic physical
branch.  A 100,000-batch lunar-regolith validation reached batch 14,400 and
failed at `E_interface = -0.7289832458 V/m`, even though the failure did not
establish that the boundary-value problem had no physical solution.

This change makes the nonlinear solve robust without changing the kinetic
density closures, the Bohm entry condition, the monotonic branch definition,
or the accepted residual tolerance.

## Physical Invariants

- The solved equations remain the one-dimensional Poisson equation with the
  absorbing ambient-electron, cold drifting-ion, and emitted-photoelectron
  closures already implemented in `bem_outer_plasma_kinetic`.
- The interface Neumann condition and far-side Robin condition are unchanged.
- A result is accepted only when the original residual infinity norm is no
  larger than `residual_tolerance`.
- The electron-repelling monotonic branch remains a hard feasibility
  constraint.  Numerical globalization may not cross into an inaccessible
  closure state.
- A violated Bohm condition or inaccessible density closure is classified as
  `outer_plasma_no_physical_solution`.  Exhausting numerical recovery methods
  is classified as `outer_plasma_numerical_failure`.
- No linear-Debye, previous-profile, or other physically different fallback is
  returned as a converged kinetic solution.

## Analytic Jacobian

Each interior Poisson row depends locally on `phi(j-1:j+1)` and nonlocally on
`phi(1)` through the absorbing-electron and emitted-photoelectron closures.
The analytic Jacobian is therefore represented as a tridiagonal matrix plus a
first-column border.  It is not stored as a dense matrix.

The closure routines return both density and derivatives with respect to local
potential and interface potential.  Existing electron and ion susceptibility
expressions are retained but made explicit.  The emitted-photoelectron
derivatives are added from the same piecewise formulas as its density.  At the
piecewise square-root endpoints, the solver uses the analytic limiting branch
where finite and a one-sided, scale-aware derivative only for the singular
endpoint itself.

The bordered system is solved in linear time by solving the tridiagonal base
system for the right-hand side and the border vector, followed by the scalar
rank-one correction.  Pivot and denominator checks must reject non-finite or
numerically singular systems.

## Nonlinear Globalization

The primary path is analytic Newton with branch-preserving backtracking.  A
trial is accepted when it is physically feasible and reduces the original
residual infinity norm.

If Newton cannot produce an accepted step, the solver enters pseudo-transient
continuation.  The diagonal is shifted by a positive inverse pseudo-time term,
and the pseudo-time is increased after successful residual reduction and
decreased after rejection.  This changes only the iteration path: convergence
is still checked against the unshifted original residual.  Once a sufficiently
Newton-like step is available, the solver returns to the primary path.

The nonlinear iteration budget counts all Newton and pseudo-transient steps.
Diagnostics record the total iterations, final residual, and whether
pseudo-transient recovery was used.

## Interface-Field Continuation

When an initial profile is supplied, its interface field is reconstructed and
the requested field is reached through adaptive parameter continuation.  A
successful subsolve permits a larger next field increment; a numerical failure
halves the increment and retries from the last converged profile.  Only the
target-field profile is returned.  Closure-inaccessibility at a trial field is
not reclassified as numerical success.

The zero-initial-profile path uses the existing analytic seed and nonlinear
globalization directly.  Continuation has bounded subdivision and terminates
with `outer_plasma_numerical_failure` if the minimum progress threshold is
exhausted.

## Performance

For `N` grid points, a nonlinear step changes from `N` residual evaluations
plus an `O(N^3)` dense solve to one residual/Jacobian evaluation plus an `O(N)`
bordered-tridiagonal solve.  The default `N=128` solve should therefore be
cheap relative to a particle batch and suitable for the existing periodic
outer-plasma update cadence.

## Verification

1. Compare every analytic closure derivative and the assembled Jacobian with
   scale-aware finite differences away from piecewise singular endpoints.
2. Retain vacuum analytic, grid refinement, Gauss closure, photoelectron,
   warm-start, and sub-Bohm regression tests.
3. Add a regression that crosses the previously failing lunar field
   `-0.7289832458 V/m` while maintaining the monotonic branch and residual
   tolerance.
4. Run L1 and L2 on a SysA compute node.
5. Run the 3,000-batch smoke case and the 100,000-batch lunar-regolith case;
   the latter must pass the old batch-14,400 failure point and complete without
   accepting a residual larger than `1e-8`.
