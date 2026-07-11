# P0 Panel Side Contract

## Scope

This contract defines the free-space, direct-solver `triangle_p0` source model introduced in Phase 1. It fixes geometry ownership, source normalization, one-sided limits, self terms, and unsupported combinations before FMM or periodic panel support is added.

The existing `point` source model remains the default and retains its current numerical behavior.

## Geometry And Charge

For an ordered triangle `(v0, v1, v2)`, the winding normal and area are

```text
n_w = cross(v1 - v0, v2 - v0) / |cross(v1 - v0, v2 - v0)|
A   = |cross(v1 - v0, v2 - v0)| / 2
```

Degenerate or non-finite triangles are invalid. `q_elem(i)` is the total charge on the panel and `sigma_i = q_elem(i) / A_i` is its uniform P0 surface-charge density.

The mesh owns `elem_vacuum_sign(i)` with value `+1` or `-1`, and derives

```text
n_vac(i) = elem_vacuum_sign(i) * n_w(i)
```

Triangle winding is preserved. Physical side metadata never rewrites vertex order.

## Side Policies

- `normal_plus`: vacuum is the `+n_w` side; `elem_vacuum_sign=+1`.
- `normal_minus`: vacuum is the `-n_w` side; `elem_vacuum_sign=-1`.
- `outward_closed`: every connected mesh group must be closed and consistently orientable. Positive signed volume gives `+1`; negative signed volume gives `-1`. Open, non-manifold, or inconsistently wound groups are rejected.

An open surface in `triangle_p0` mode requires `normal_plus` or `normal_minus`. The legacy `point` mode does not interpret surface-side metadata.

## Kernel

For target `x` and triangle `T_i`,

```text
phi_i(x) = k_c * q_i / A_i * integral_Ti 1 / |x-y| dS_y
E_i(x)   = k_c * q_i / A_i * integral_Ti (x-y) / |x-y|^3 dS_y
```

Potential is continuous across the panel. Electric field is evaluated using an explicit side:

- `principal_value`: finite-part tangential field with zero self normal contribution.
- `normal_plus`: `E_PV + sigma/(2*eps0) * n_w`.
- `normal_minus`: `E_PV - sigma/(2*eps0) * n_w`.
- `vacuum`: `E_PV + sigma/(2*eps0) * n_vac`.

Therefore

```text
E_normal_plus - E_normal_minus = sigma / eps0 * n_w
```

Off-surface evaluation has no side ambiguity and converges to these limits.

## Analytic Direct Evaluator

The production direct evaluator uses the polygon edge/solid-angle representation. For projection `p`, signed height `h=dot(x-v0,n_w)`, each edge outward co-normal `m_e`, edge distance `d_e=dot(v_e-p,m_e)`, edge line integral `L_e`, and signed solid angle `Omega` chosen so `sign(Omega)=sign(h)`:

```text
I_phi = sum_e d_e * L_e - h * Omega
I_E   = sum_e m_e * L_e + n_w * Omega
```

At an on-panel target, `h=0` and `Omega=0` define the principal value. The centroid self potential uses the same edge formula, not `1/softening`.

A separate adaptive/Duffy quadrature implementation is the correctness oracle. It is not used silently as a point-source fallback.

## Moments

Geometry initialization precomputes area, centroid, winding normal, and normalized raw moments required by later panel-integrated multipoles. Phase 1 tests exact constant, first, and second Cartesian moments under rigid transforms, vertex permutations, and scale changes.

## Configuration And Validation

`field.element_kernel` accepts `point` and `triangle_p0`. `triangle_p0` is available only with:

- `field_solver="direct"`
- `field_bc_mode="free"`
- insulator elements only
- `softening=0`
- a resolved vacuum side for every element

`auto`, treecode, FMM, periodic2, conductor, dielectric, nonzero softening, and unresolved sides are errors. They never fall back to point sources.

Python reconstruction and the point-only C kernel reject output with `field_source_model=triangle_p0` until they implement this same kernel contract.

## Verification Contract

Phase 1 requires:

- analytic direct results agree with high-order quadrature off surface;
- potential limits from both sides agree;
- the normal-field jump equals `sigma/eps0`;
- centroid self potential agrees with singular quadrature;
- rigid transforms and cyclic vertex permutations preserve results;
- charge scaling is linear;
- legacy point tests remain unchanged;
- every unsupported solver, surface, and Python path fails closed.

