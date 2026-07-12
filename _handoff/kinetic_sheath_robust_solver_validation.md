# Kinetic Sheath Robust Solver Validation

## Scope

This report records validation of the analytic bordered-tridiagonal Newton,
pseudo-transient recovery, and adaptive interface-field continuation added to
the `kinetic_1d` outer-plasma solver.

The physical equations, monotonic branch, Bohm condition, Neumann/Robin
boundaries, and final residual tolerance are unchanged. No alternate physical
model is used as a fallback.

## Implementation

- Branch: `feature/robust-kinetic-sheath-solver`
- Design commit: `1cf584b`
- Closure derivative commit: `a8221c6`
- Structured Jacobian and recovery commit: `e673aac`
- Adaptive continuation commit: `f56122f`

## Focused SysA Evidence

| Job | Target | Result |
|---|---|---|
| `8113020` | pre-change kinetic core + solver baseline | 7/7 core and 9/9 solver tests passed |
| `8113022` | closure derivative RED | expected compile failure on missing derivative API |
| `8113028` | closure derivative GREEN | 10 tests, 24 assertions passed |
| `8113031` | analytic Jacobian RED | expected compile failure on missing Jacobian action API |
| `8113034` | analytic Jacobian without recovery | Jacobian action passed; lunar `-0.70 V/m` source reproduced Newton line-search failure |
| `8113040` | pseudo-transient GREEN | 10 core and 11 solver tests passed; `-0.70` and `-0.7289832458 V/m` converged |
| `8113041` | continuation RED | expected compile failure on missing continuation diagnostic |
| `8113047` | continuation GREEN | 10 core and 11 solver tests passed, 41 solver assertions |
| `8113048` | `make check` | debug OpenMP build passed on SysA |

The difficult-field regression first solves the lunar ambient state at zero
field, reaches `-0.70 V/m` through multiple continuation steps, then reaches
the former runtime failure field `-0.72898324579369622 V/m`. The final state is
monotonic and its original residual is at most `1e-8`.

## Pending Tiered and Runtime Evidence

- L1 result: pending
- L2 result: pending
- 3,000-batch smoke result: pending
- 100,000-batch lunar-regolith result: pending
