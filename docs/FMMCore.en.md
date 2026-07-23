title: Coulomb FMM Core Details

Lang: [English](FMMCore.en.md) | [日本語](FMMCore.md)

# Coulomb FMM Core Details

This section summarizes the specification and algorithms of the current Fortran Coulomb FMM core,
[`bem_coulomb_fmm_core` module page](../module/bem_coulomb_fmm_core.html),
and its split implementation files.

See [FMM](FMM.en.html) for the user-facing equations and computation flow, and
[periodic2 Far Correction](PeriodicFarCorrection.en.html) for root-operator construction and operation. This page focuses on
Fortran internal arrays and implementation steps.

- Public API / boundary: `src/physics/field_solver/fmm/api/`
- Internal shared implementation: `src/physics/field_solver/fmm/internal/common/`
- Tree / plan implementation: `src/physics/field_solver/fmm/internal/tree/`
- State / eval implementation: `src/physics/field_solver/fmm/internal/runtime/`
- periodic2 implementation: `src/physics/field_solver/fmm/internal/periodic/`

The target is a simulator-independent internal API. It does not directly `use` `mesh_type` or `sim_config`.
On the BEACH side, the field-solver adapter calls this core.

## 1. Purpose

The FMM core returns Coulomb electric fields quickly at many evaluation points for a fixed source point set `src_pos(3,n)` and variable charges `src_q(n)`.

Current design goals:

- kernel is only 3D Coulomb
- source geometry and charge updates are separated
- only `free` and `periodic2` are supported
- near direct sum is also handled inside the core
- simulator code sees only array APIs

## 2. Public API

The core provides four main procedures:

```fortran
call build_plan(plan, src_pos, options)
call update_state(plan, state, src_q)
call eval_points(plan, state, target_pos, e)
call eval_point(plan, state, r, e)
```

Input and output meanings:

- `src_pos(3,n)`: source point coordinates, fixed after `build_plan`
- `src_q(n)`: source charges, updateable at each `update_state`
- `target_pos(3,m)` or `r(3)`: evaluation points
- `e(3,m)` or `e(3)`: electric field vectors

Notes:

- The returned field does not include `k_coulomb`; the BEACH adapter multiplies it at the end.
- `build_plan` is geometry-dependent processing, and `update_state` is charge-dependent processing.
- `eval_point(s)` assumes `plan` and `state` are ready.

#### 2.2 C ABI / Python integration

`src/physics/field_solver/bem_field_kernel_c.f90` exposes this Fortran API as an `iso_c_binding` opaque-handle API.
`make build-kernel` builds the shared library as `build/libbeach_field_kernel.so`.

Main C ABI:

```text
beach_kernel_get_abi_version(major, minor)
beach_kernel_get_build_info(buffer, capacity, length)
beach_kernel_create(handle)
beach_kernel_destroy(handle)
beach_kernel_build(handle, src_pos, options...)
beach_kernel_update_charges(handle, src_q)
beach_kernel_eval_e(handle, target_pos, e)
beach_kernel_eval_phi(handle, target_pos, phi)
beach_kernel_force_on_charges(handle, target_pos, target_q, origin, force, torque)
```

The public header is `beach/include/beach_field_kernel.h` in the Python package.
The current ABI is `1.0`. Before calling the other functions, a C caller should
call `beach_kernel_get_abi_version` and require an equal major version and a
library minor version greater than or equal to its required minor version.
Coordinate and vector arrays use `values[3 * point_index + component]` storage.
The public header defines status codes, periodic far-correction codes, and handle
ownership.

The Python side calls this ABI with `ctypes` through `beach.fortran_results.kernel.FieldKernel`.
The Python wrapper checks compatible libraries that provide the version-query
symbol when loading them. It accepts older libraries without that symbol for
transition compatibility, while newly built libraries are required to provide it.
`calc_object_forces_kernel` evaluates `sum(q_i E_not_self(r_i))` by zeroing the object's own source charge, avoiding self-force contamination while using the same field kernel, including `periodic2 + cached_kneq0`.
`Beach.scene()` / `BeachScene` temporarily apply rigid translations and rotations of objects on the Python side and pass the edited centroid array to the same ABI.
The rigid-transform helper path uses NumPy by default and can use an optional Numba backend, but field evaluation itself is done by the Fortran kernel.

#### 2.3 BEACH adapter usage

The BEACH field-solver adapter passes different geometry to the core according
to the source model.

| source model | plan construction | meaning of `src_q(i)` |
|---|---|---|
| `point` | pass element centroids as `src_pos` to `build_plan` | point charge located at the centroid |
| `triangle_p0` | pass all three vertices to `build_panel_plan` | total charge on the triangle; surface density is `src_q(i)/area(i)` |

- During initialization, it calls `update_state` immediately after `build_plan` or `build_panel_plan`.
- During later refreshes, normal operation assumes mesh geometry is unchanged, so the existing `plan` is reused and only `update_state` is called with updated `src_q`.
- `build_plan` and legacy tree metadata are synchronized again only when the plan is missing, the source count changes, or zero elements caused plan/state disposal.

## 3. Data structures

#### 3.1 `fmm_options_type`

Main internal options:

- `theta`: parameter for well-separated tests
- `leaf_max`: maximum source count in a source-octree leaf
- `order`: Cartesian expansion order
- `softening`: `epsilon` of the softened Coulomb kernel for `point` sources; it must be zero for `triangle_p0`
- `use_periodic2`: enable two-periodic-axis mode
- `periodic_axes(2)`, `periodic_len(2)`: periodic axes and lengths
- `periodic_image_layers`: near image-sum layer count `N`
- `periodic_far_correction`: core values are `auto`, `none`, and `cached_kneq0`; with `periodic2`, `auto` is normalized to `none`
- `periodic_ewald_alpha`, `periodic_ewald_layers`: decomposition parameter and cutoff depth used by the build-time Ewald fit for `cached_kneq0`
- `target_box_min/max`: box used for a dual-target tree

The BEACH adapter currently uses `order = 4`, while the core itself accepts variable order.
For `periodic2`, `auto` is normalized to `none`; `cached_kneq0` explicitly enables far correction.

#### 3.2 `fmm_plan_type`

This is geometry-dependent immutable data:

- multi-index tables `alpha`, `deriv_alpha`
- source octree
- optional target tree
- source leaf list `source_leaf_nodes`
- target leaf list `leaf_nodes`
- near lists `near_start/near_nodes`
- far node lists `far_start/far_nodes`
- M2L pair cache `m2l_target_nodes/m2l_source_nodes`
- periodic image-shift arrays
- M2L derivative table `m2l_deriv`
- P2M basis table `source_p2m_basis`
- compressed translation tables for M2M/L2L

#### 3.3 `fmm_state_type`

This is charge-dependent data updated on each refresh:

- `src_q(n)`
- `multipole(ncoef, nnode)`
- `local(ncoef, n_target_nodes)`
- `multipole_active(nnode)`
- `local_active(n_target_nodes)`

`multipole` stores multipole coefficients per source-tree node, and `local` stores local expansion coefficients per target-tree node.
`*_active` flags are 0/1 flags used to skip zero nodes quickly.

## 4. Mathematical definitions

#### 4.1 Source kernels

The core supports two source kernels: `point` and `triangle_p0`.

##### Point sources

`point` evaluates charges located at element centroids with the softened
Coulomb kernel:

$$
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\lVert\mathbf{r}\rVert^2 + \epsilon^2}}
$$

$$
\phi(\mathbf{x}) = \sum_j q_j \, G_\epsilon(\mathbf{x} - \mathbf{x}_j)
$$

$$
\mathbf{E}(\mathbf{x}) = - \nabla \phi(\mathbf{x})
$$

The near direct sum and far-field expansion are both based on this
$G_\epsilon$.

##### P0 triangle sources

`triangle_p0` treats `q_i` as the total charge on triangle $T_i$ of area $A_i$,
with constant surface density $\sigma_i=q_i/A_i$:

$$
\phi_i(\mathbf{x}) =
\frac{q_i}{A_i}\int_{T_i}
\frac{1}{\lVert\mathbf{x}-\mathbf{y}\rVert}\,dA_{\mathbf{y}}
$$

Near direct evaluation uses the analytic P0 panel kernel based on logarithmic
edge terms and the solid angle. Far-field P2M uses exact area-averaged
monomials over the triangle relative to tree-node center $\mathbf{c}$:

$$
M_{i,\alpha} =
q_i\left\langle(\mathbf{y}-\mathbf{c})^\alpha\right\rangle_{T_i}
= \frac{q_i}{A_i}\int_{T_i}
(\mathbf{y}-\mathbf{c})^\alpha\,dA_{\mathbf{y}}
$$

Area weighting is therefore contained in the panel integral and the P2M basis.
Because `q_i` is already the total element charge, it is not multiplied by
$A_i$ again. M2M/M2L/L2L then use the unsoftened Coulomb/Laplace expansion for
these panel moments. `triangle_p0` enforces `softening=0` and never falls back
to softened point sources.

#### 4.2 Multi-index

The core uses a multi-index $\alpha = (\alpha_x, \alpha_y, \alpha_z)$.

$$
|\alpha| = \alpha_x + \alpha_y + \alpha_z
$$

$$
\alpha! = \alpha_x! \, \alpha_y! \, \alpha_z!
$$

$$
\mathbf{r}^\alpha = r_x^{\alpha_x} r_y^{\alpha_y} r_z^{\alpha_z}
$$

#### 4.3 P2M

For node center $c$, leaf-node multipole coefficients are:

$$
M_\alpha(c) = \sum_{j \in \text{leaf}} q_j
\frac{(\mathbf{x}_j - \mathbf{c})^\alpha}{\alpha!}
$$

#### 4.4 M2M

Child-node coefficients are translated to the parent center and accumulated.
With $\mathbf{d} = c_{\mathrm{child}} - c_{\mathrm{parent}}$:

$$
M_\beta(c_{\mathrm{parent}}) =
\sum_{\alpha \le \beta}
M_\alpha(c_{\mathrm{child}})
\frac{\mathbf{d}^{\beta-\alpha}}{(\beta-\alpha)!}
$$

The current implementation precomputes, during `build_plan`, the index for $\beta - \alpha$ and the value
$\mathbf{d}^{\beta-\alpha} / (\beta-\alpha)!$.

#### 4.5 M2L

For source-node center $c_s$ and target-node center $c_t$, let $R = c_t - c_s$.

Local expansion coefficients are updated as:

$$
L_\alpha(c_t) \mathrel{+}=
\sum_\beta (-1)^{|\beta|}
M_\beta(c_s)
D^{\alpha+\beta} G_\epsilon(R)
$$

Here $D^\gamma$ is a multi-index derivative.
The current implementation precomputes $D^{\alpha+\beta} G_\epsilon(R)$ per pair as `m2l_deriv(:, pair)`.

#### 4.6 L2L

The local expansion at parent center $c_{\mathrm{parent}}$ is translated to child center $c_{\mathrm{child}}$.
With $\mathbf{d} = c_{\mathrm{child}} - c_{\mathrm{parent}}$:

$$
L_\alpha(c_{\mathrm{child}}) \mathrel{+}=
\sum_{\gamma \ge \alpha}
L_\gamma(c_{\mathrm{parent}})
\frac{\mathbf{d}^{\gamma-\alpha}}{(\gamma-\alpha)!}
$$

The shift monomials are also precomputed during `build_plan`.

#### 4.7 L2P

Let $c_{\mathrm{leaf}}$ be the center of the target leaf that contains evaluation point $x$, and let
$\mathbf{dr} = x - c_{\mathrm{leaf}}$.

$$
E_k(x) = - \sum_{|\alpha| < p} L_{\alpha + e_k}(c_{\mathrm{leaf}}) \frac{\mathbf{dr}^\alpha}{\alpha!}
$$

Here $e_k$ is the unit multi-index for axis $k$.

## 5. `build_plan` algorithm

`build_plan` performs only geometry-dependent work.

#### 5.1 Source tree

The source-coordinate bounding box is recursively split into eight octants to build the octree.
The stopping condition is either:

- source count `<= leaf_max`
- the bounding box is small enough that further subdivision is not useful

#### 5.2 Target topology

There are two target-side modes:

- `target_box` disabled:
  reuse source-tree leaves as target leaves
- `target_box` enabled:
  build a separate target tree that covers the whole box

In `periodic2`, target points are wrapped into the box before target leaf lookup.

#### 5.3 Near/far lists and M2L pair cache

For each target leaf, the source tree is traversed recursively to build near nodes and far nodes.

The well-separated test is:

$$
(r_s + r_t)^2 < \theta_{\mathrm{eff}}^2 \, \lVert\mathbf{d}\rVert^2
$$

where:

- $r_s$ is the source-node radius
- $r_t$ is the target-node radius
- $\mathbf{d}$ is the vector between node centers
- $\theta_{\mathrm{eff}} = \theta$ for both `free` and `periodic2`

In `periodic2`, a minimum-image correction is applied to $\mathbf{d}$.

Then a dual-tree recursion builds the M2L pair cache and prepares index arrays per target node.

#### 5.4 Build-time precomputation

At the end of `build_plan`, quantities that do not change between refreshes are precomputed:

- `source_parent_of`
- `parent_of`
- `source_p2m_basis`
- `m2m_term_count`, `m2m_alpha_list`, `m2m_delta_list`
- `l2l_term_count`, `l2l_gamma_list`, `l2l_delta_list`
- `source_shift_monomial`
- `target_shift_monomial`
- `shift_axis1`, `shift_axis2`
- `periodic_ewald`
- `periodic_root_operator`
- `m2l_deriv`

This makes `update_state` close to charge-dependent accumulation only.

#### 5.5 Pseudocode

```text
build_plan(src_pos, options):
  initialize_basis_tables(order)
  build_source_tree(src_pos)
  precompute_source_p2m_basis()
  build_target_topology(target_box)
  build_interactions()
  precompute_translation_operators()
  precompute_periodic2_ewald_data()
  precompute_periodic_root_operator()
  precompute_m2l_derivatives()
```

## 6. `update_state` algorithm

`update_state` corresponds to refresh in the legacy implementation.
Source coordinates are fixed; only `src_q` changes.

#### 6.1 Processing order

```text
update_state(plan, state, src_q):
  ensure_state_capacity()
  copy src_q
  clear active flags
  clear multipole/local only when the tree has no source leaves or no M2L pairs
  P2M on source leaves
  M2M bottom-up
  M2L on cached pairs
  L2L top-down
  mark state ready
```

#### 6.2 OpenMP parallelization

OpenMP is currently used in:

- one parallel region around the full `update_state`, including `src_q` copy and active-flag initialization
- `P2M`: loop over source leaves
- `M2M`: loop over nodes at the same depth
- `M2L`: loop over target nodes
- `L2L`: loop over nodes at the same depth
- translation and M2L derivative precomputation during `build_plan`

The loops are written to map roughly one node to one thread, and shared-array updates are independent at node granularity.

#### 6.3 Implementation optimizations

`update_state` avoids unnecessary work by:

- not recomputing the multi-index difference $\beta - \alpha$
- not rebuilding powers of parent-child center shifts
- precomputing the `P2M` monomial basis per source during build
- storing only valid compressed `(alpha, delta)` terms for `M2M/L2L`
- using source-node active flags to skip zero nodes in `M2L` per pair
- accumulating `M2L` contributions in thread-local `local_acc` before writing back to target-node columns
- using source-leaf-specific indices in `P2M`, not target-leaf indices

## 7. `eval_point(s)` algorithm

Evaluation proceeds as:

```text
eval_point(r):
  if plan is not built or state is not ready:
    return zero vector

  if periodic2:
    wrap r into target box

  leaf = locate_target_leaf(r)
  if leaf not found or leaf is not mapped to a leaf slot:
    use direct sum over all sources
    return

  evaluate local expansion at leaf center
  add near direct interactions
  root local already carries periodic root correction when enabled
```

#### 7.1 Leaf lookup

- In `periodic2`, the evaluation point is wrapped into the target box before lookup.
- If a target tree exists, its leaves are used.
- If no target tree exists, source-tree leaves are used.
- If lookup fails, or the leaf cannot map to a tree leaf slot, evaluation falls back to direct sum.

#### 7.2 Near direct

Source indices in the near list are evaluated by direct sum.
In `periodic2`, image shifts in `[-N, N] x [-N, N]` are handled explicitly.
Fallback uses the same direct kernel.

#### 7.3 Out-of-box fallback

When a dual-target tree is used, evaluation points can leave the target box.
Then there is no target leaf, so evaluation falls back to direct sum over all sources.
`cached_kneq0` assumes a fixed target topology and rejects evaluation outside the target box.

#### 7.4 Location of root correction

The `cached_kneq0` root correction is injected into target-anchor local expansions during `update_state`.
Therefore normal leaf evaluation in `eval_point(s)` does not recompute the root correction; it just uses the local expansion carried by `state`.

## 8. `periodic2` and far correction

This section retains formulas and fallback details from the FMM-core viewpoint. Configuration selection, operator fitting,
cache lifecycle, and `k=0` ownership are separated into [periodic2 Far Correction](PeriodicFarCorrection.en.html).

#### 8.1 `periodic2`

`periodic2` means exactly two axes are periodic and the remaining axis is open.

The near image sum explicitly adds the finite images:

$$
i, j \in [-N, N]
$$

M2L uses the same image-shift set and precomputes each pair derivative as an image sum.

#### 8.2 `periodic2` Ewald (Ewald2P) correction

`bem_coulomb_fmm_periodic_ewald.f90` implements an Ewald-form correction for the two-periodic, one-open Coulomb field.
Here `exact` means the finite sum actually evaluated by the code. It is not the theoretical infinite sum; it is a build-time oracle whose real-space and reciprocal-space cutoffs are controlled by `field_periodic_image_layers = N` and `field_periodic_ewald_layers = L`.

Ewald2P is a build-time teacher for `cached_kneq0`, not the
runtime particle kernel. In the cached `triangle_p0` path, the teacher is applied
to proxy point charges and fitted as a root-multipole-to-local operator. Real
triangles still use the analytic panel kernel in the near field and
triangle-averaged P2M in the far source representation.

$\alpha$ is a numerical parameter balancing real- and reciprocal-space
convergence; it is not Debye screening. See
[the Ewald2P teacher in periodic2 electrostatics](PeriodicElectrostatics.en.md#ewald2p-teacher)
for the intuitive split and its relation to runtime evaluation.

##### 8.2.1 Notation

Let the periodic axes be `a_1, a_2` and the open axis be `f`.
Define periodic lengths, cell area, image set, and reciprocal-lattice set as:

$$
L_1 = \operatorname{periodic\_len}(1),\qquad
L_2 = \operatorname{periodic\_len}(2),\qquad
A = L_1 L_2
$$

$$
\mathcal I_N = \{(i,j)\in\mathbb Z^2 \mid |i|,|j|\le N\},\qquad
\mathcal K_L = \{(m,n)\in\mathbb Z^2 \mid |m|,|n|\le L,\ (m,n)\neq(0,0)\}
$$

Image shifts and reciprocal-lattice vectors are:

$$
\mathbf L_{ij} = iL_1\,\mathbf e_{a_1} + jL_2\,\mathbf e_{a_2},\qquad
\mathbf k_{mn} = 2\pi\left(\frac{m}{L_1}\mathbf e_{a_1} + \frac{n}{L_2}\mathbf e_{a_2}\right)
$$

For source position $\mathbf s$ and evaluation point $\mathbf r$, define:

$$
\mathbf R_{ij} = \mathbf r - \mathbf s - \mathbf L_{ij},\qquad
R_{ij} = \lVert\mathbf R_{ij}\rVert,\qquad
z = (\mathbf r - \mathbf s)\cdot \mathbf e_f
$$

Below, $\alpha =$ `field_periodic_ewald_alpha` and $\epsilon =$ `softening`.

##### 8.2.2 Real-space term

The screened Coulomb field implemented by `add_screened_point_charge` is:

$$
\mathbf E_\alpha(\mathbf R) =
q\left(
\frac{\operatorname{erfc}(\alpha R)}{R^3}
{}+\frac{2\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 R^2}}{R^2}
\right)\mathbf R
$$

It is the gradient of the potential:

$$
\Phi_\alpha(\mathbf R) = q\,\frac{\operatorname{erfc}(\alpha R)}{R}
$$

The direct kernel used by `add_softened_point_charge` is:

$$
\mathbf E_\epsilon(\mathbf R) =
q\,\frac{\mathbf R}{(R^2+\epsilon^2)^{3/2}}
$$

It uses the same softening as the normal runtime direct path.

The implemented real-space correction is:

$$
\mathbf E_{\mathrm{real}} =
\sum_{(i,j)\in\mathcal I_{N+L}} \mathbf E_\alpha(\mathbf R_{ij})
{}- \sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
$$

Terms with `r2 <= tiny(1.0d0)` are skipped, so self-interaction is excluded.
If the direct fallback contribution $\sum_{(i,j)\in\mathcal I_N}\mathbf E_\epsilon$ is added to `add_periodic2_exact_ewald_correction_single_source`, the softened inner-image part cancels and the outer shell is replaced by the screened form.

##### 8.2.3 Reciprocal-space term

For $(m,n)\neq(0,0)$, `add_exact_periodic2_reciprocal_space_correction` defines:

$$
\theta_{mn} = \mathbf k_{mn}\cdot(\mathbf r-\mathbf s),\qquad
k_{mn} = \lVert\mathbf k_{mn}\rVert
$$

$$
G^\pm_{mn}(z) =
e^{\pm k_{mn} z}\operatorname{erfc}\!\left(\frac{k_{mn}}{2\alpha}\pm \alpha z\right)
$$

and uses:

$$
\mathbf E_{\mathrm{rec}} =
q \sum_{(m,n)\in\mathcal K_L}
\frac{\pi}{A}
\begin{pmatrix}
\frac{(\mathbf k_{mn})_{a_1}}{k_{mn}}\sin\theta_{mn}\,\bigl(G^+_{mn}(z)+G^-_{mn}(z)\bigr) \\
\frac{(\mathbf k_{mn})_{a_2}}{k_{mn}}\sin\theta_{mn}\,\bigl(G^+_{mn}(z)+G^-_{mn}(z)\bigr) \\
\cos\theta_{mn}\,\bigl(G^-_{mn}(z)-G^+_{mn}(z)\bigr)
\end{pmatrix}
$$

In code these correspond to `term_p`, `term_m`, and `pair_sum`.
This term represents the high-frequency reciprocal-lattice components excluding `k=0`.

##### 8.2.4 `k=0` term

The zero-mode correction implemented by `add_exact_periodic2_k0_correction` is:

$$
\mathbf E_0 =
q\,\frac{2\pi}{A}\operatorname{erf}(\alpha z)\,\mathbf e_f
$$

The single-source Ewald teacher keeps this form as the `k=0` electric-field contribution.

##### 8.2.5 Implemented correction

Together, the correction added by `add_periodic2_exact_ewald_correction_single_source` for one source is:

$$
\mathbf E_{\mathrm{corr}} =
\mathbf E_{\mathrm{real}}
{}+ \mathbf E_{\mathrm{rec}}
{}+ \mathbf E_0
$$

This single-source correction is evaluated at proxy/check points to generate the cached operator.

If `field_periodic_ewald_alpha <= 0`, `resolve_periodic2_ewald_alpha` selects:

$$
\alpha = \frac{1.2}{(N+1)\min(L_1,L_2)}
$$

automatically. If `min(L_1,L_2) <= 0`, it sets `alpha = 0` and disables operator generation.
Internally, `kmax = max(1, field_periodic_ewald_layers)` defines the reciprocal-space finite sum.

The `cached_kneq0` cold build evaluates Ewald-residual field and potential at check points and fits a
root-multipole-to-target-local operator. The constant-potential coefficient, which is not determined by the field fit, is fitted
separately from the potential residual and stored in the versioned cache. `m2l_root_oracle` has been removed and is rejected.
Infinite-periodic production runs use `cached_kneq0`.

## 9. Interpreting computational cost

With fixed order $p$ and bounded interaction lists, practical costs are approximately:

- `build_plan`: close to $O(n \log n)$
- `update_state`: close to $O(n)$
- `eval_point`: close to $O(\log n + n_{\mathrm{near}} \, n_{\mathrm{img}}^2)$
- `eval_points`: parallel execution of the above point evaluation for each target

The constant factors depend strongly on:

- `order`
- `theta`
- `leaf_max`
- `periodic_image_layers`
- `periodic_ewald_layers`
- whether a target tree exists

## 10. Current implementation limits

This FMM core is not a generic kernel FMM.

- kernel is fixed to Coulomb
- the simulator adapter default order is `order = 4`
- source coordinates are considered immutable after `build_plan`
- supported boundaries are `free` and `periodic2`
- `periodic2` requires exactly two periodic axes
- far correction modes are `none` by default, `auto`, and `cached_kneq0`; `periodic2` `auto` normalizes to `none`
- `eval_point(s)` return values do not include `k_coulomb`

### 10.1 Cached periodic nonzero operator

#### What the operator accelerates

In `periodic2`, copies of the primary-cell charge distribution extend infinitely
along x/y. The ordinary FMM efficiently evaluates the near image shell
`[-N,N]^2`, but a smooth contribution from all images outside that shell remains.
Instead of running an Ewald sum for every particle evaluation, `cached_kneq0`
precomputes only this far difference as a linear map into FMM local expansions.

| Input | Cached operator | Output |
| --- | --- | --- |
| root multipole built from current charges | apply a geometry-fixed matrix | far local expansion for each target anchor |

The cache does not contain field samples, particle positions, or charge history.
It contains a geometry-specific matrix mapping source multipoles to far local
expansions, so it remains reusable when charges change between batches.

#### What one field evaluation adds

| Order | Component | Purpose |
| ---: | --- | --- |
| 1 | primary cell plus finite near images | evaluate singular/near interactions with ordinary FMM and direct kernels |
| 2 | cached Ewald residual | restore the smooth infinite-periodic field outside the finite shell |
| 3 | subtract the symmetric `k=0` carried by the cached teacher | leave only `k!=0` in the nonzero backend |
| 4 | snapshot adds the physical `k=0` | apply `symmetric_vacuum`, `e_bottom_zero`, or outer-plasma physics |

`cached_kneq0` alone is therefore not the complete field. Steps 1--3 belong to
the nonzero backend; step 4 belongs to `electrostatic_snapshot`. `exclude_k0`
does not discard the mean field. It prevents double counting because a separate
boundary provider adds that field exactly once.

#### Relation to the formula

Steps 1--3 give the runtime nonzero-mode kernel:

$$
K_{k\ne0}
= K_\mathrm{shell}(N)
+ R_\mathrm{Ewald}^{\mathrm{full}}
- K_0^\mathrm{sym}
$$

| Term | Built when | Role |
| --- | --- | --- |
| $K_\mathrm{shell}(N)$ | ordinary FMM plan/runtime | primary cell and finite near images |
| $R_\mathrm{Ewald}^{\mathrm{full}}$ | cold cache build | smooth difference between full-periodic Ewald and the finite image shell |
| $K_0^\mathrm{sym}$ | charge-state refresh | remove the symmetric $k=0$ part from the cached full-periodic kernel |

$K_0^\mathrm{sym}$ evaluates a piecewise-polynomial source-height prefix state by binary search in $O(\log n)$ per target.
The final surface field is

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}.
$$

The triangle-height integral, lower-boundary closure, and outer-plasma coupling
for $K_0^\mathrm{physical}$ are described in
[periodic2 electrostatics](PeriodicElectrostatics.en.md) and
[Outer plasma models](OuterPlasmaModels.en.md).

Field-fit columns and the constant potential mode have different units and are not mixed in one least-squares system.
The potential gauge is fixed separately from the mean residual.

#### Cache lifecycle

| Stage | Action |
| --- | --- |
| Build identity | fingerprint geometry, target topology, order, periods, image layers, and generator/build versions |
| Warm read | accept only matching version, fingerprint, shape, and checksum |
| Miss or corruption | acquire the filesystem lock and regenerate the operator |
| Publish | close a same-directory `.tmp` file, then atomically rename it |
| Checkpoint | omit the regenerable operator payload |

Independent jobs therefore serialize on one lock file, and a reader cannot accept a partially written operator.

#### MPI/OpenMP cold build

| Unit of work | Owner |
| --- | --- |
| Cache I/O and lock | MPI rank 0 only |
| Target operator slices | distributed across MPI ranks with at most one-target imbalance |
| Proxy columns within one target | evaluated with OpenMP |
| Regularized QR | built once per target and reused for every proxy RHS |
| Complete operator | assembled on every rank with `MPI_Allreduce(SUM)` |

Warm field evaluation and charge refresh contain neither an all-source Ewald sum nor an operator refit.

#### Cold versus warm execution

| Path | Ewald teacher | QR fit | Particle evaluation |
| --- | --- | --- | --- |
| first cache miss | run | fit and publish the operator | starts after the build |
| cache hit | skipped | skipped | uses the loaded operator |
| batch charge refresh | skipped | skipped | applies the same operator to new multipoles |

A cold build is expensive because it generates a reusable matrix once, not
because the infinite-periodic field is recomputed every batch. The warm hot path
contains no all-source Ewald sum.

#### SysA measurements

The fixture was the archived 2026-07-12 regolith input with order 4, 64 targets, 280 proxy points, and 840 check points.
The timing scope is stated explicitly because the rows do not all measure the same interval.

| Layout | Measured time | Scope |
| --- | ---: | --- |
| Former root-only, 1 rank x 1 thread | 31 min 24 s | cold operator build |
| Reusable QR, 1 rank x 1 thread | about 25 min 45 s | through operator publication |
| 1 rank x 112 threads | 47.0 s | cache prime plus batch 1 |
| 2 ranks x 112 threads | 36.7 s | cache prime plus batch 1 |
| 4 ranks x 112 threads | 31.5 s | cache prime plus batch 1 |
| 6 ranks x 112 threads | 30.3 s | cache prime plus batch 1 |

All parallel layouts produced the same cache checksum. Their Frobenius relative difference from the former operator was `1.73e-15`.
These are measurements for this fixture, not timing guarantees for arbitrary geometry.

#### Operating guidance

| Situation | Recommendation |
| --- | --- |
| Dedicated cache prime | use 1 rank x 112 threads as the core-efficiency and queue-footprint baseline |
| Production allocation already exists | generate within that allocation; 6 x 112 took about 30 s for this fixture |
| Add ranks only for cold build | marginal benefit from 4--6 ranks is small |
| 1 core | not operational; the measured job ran out of memory in the particle batch after publication |
| Warm cache | validate fingerprint and checksum, then reuse it |

## 11. Implementation mapping

Main implementation locations:

- Public API / wrapper:
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_build.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_state.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90`
- Shared type definitions:
  `fmm_options_type`, `fmm_plan_type`, `fmm_state_type`
  in `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90`
- Plan construction:
  `build_plan`
  in `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90`
- Charge refresh:
  `update_state`, `p2m_leaf_moments`, `m2m_upward_pass`, `m2l_accumulate`, `l2l_downward_pass`
  in `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90`
- Evaluation:
  `eval_point`, `eval_points`
  in `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90`
- periodic2 helpers:
  `has_valid_target_box`, `use_periodic2_cached_kneq0`,
  `use_periodic2_root_operator`, `build_periodic_shift_values`, `add_point_charge_images_field`,
  `wrap_periodic2_point`, `apply_periodic2_minimum_image`, `distance_to_source_bbox`,
  `distance_to_source_bbox_periodic`
  in `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90`
- periodic2 Ewald/oracle:
  `resolve_periodic2_ewald_alpha`, `precompute_periodic2_ewald_data`,
  `add_periodic2_exact_ewald_correction_single_source`
  in `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90`
- periodic2 root operator:
  `precompute_periodic_root_operator`
  in `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`
- BEACH adapter:
  `src/physics/field_solver/bem_field_solver_config.f90`,
  `src/physics/field_solver/bem_field_solver_tree.f90`,
  `src/physics/field_solver/bem_field_solver_eval.f90`

Design responsibilities:

- Core:
  geometry preprocessing, expansion-coefficient updates, near direct, point evaluation
- BEACH adapter:
  build `src_pos` from `mesh_type`, pass `q_elem` into `src_q`, and multiply by `k_coulomb` at the end

---
