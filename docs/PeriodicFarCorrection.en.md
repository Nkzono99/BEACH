title: periodic2 Far Correction

Lang: [English](PeriodicFarCorrection.en.md) | [日本語](PeriodicFarCorrection.md)

# periodic2 Far Correction

The `periodic2` FMM evaluates the primary cell and finite images with ordinary tree FMM, then corrects the smooth
infinite-periodic field outside that range with an additional operator. This page is the canonical reference for the Ewald2P
teacher, `cached_kneq0`, cache, and connection to FMM state.

## Select a finite-image model or an infinite-periodic approximation

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0e-4, 1.0e-4, 3.2e-4]
periodic_axes = ["x", "y"]

[fields]
boundary = "periodic2"

[fields.solver]
method = "fmm"

[fields.periodic]
backend = "cached_kneq0"
image_layers = 1
cache_dir = ".beach_cache/periodic2"
generation_tolerance = 1.0e-8
```

`[domain]` owns the fixed target topology and periodic axes; `fields.boundary` selects the periodic2 field. `cached_kneq0`
accepts only targets inside an x/y-periodic, z-nonperiodic domain box.

| `fields.periodic.backend` | Computation | Use |
| --- | --- | --- |
| `finite_images` | Primary cell and finite images only | Finite-image model and convergence comparisons |
| `cached_kneq0` | Build and reuse a versioned operator | Production infinite-periodic nonzero modes |
| `panel_spectral_reference` | Nonzero modes by a direct Fourier sum (direct solver) | Reference for small validation cases |

`cached_kneq0` does not determine the z-directed plane-average `k=0` field. Field composition adds the physical `k=0` exactly
once as described in [periodic2 electrostatics](PeriodicElectrostatics.en.html). The operator itself does not solve the
outer plasma/sheath; the outer-sheath closure supplies its potential reference in a layer separate from zero-mode construction.

## Evaluate the finite image shell with ordinary FMM

Let axes 1 and 2 be periodic with lengths $L_1,L_2$. Image layer $N$ contains

$$
\mathcal I_N=\{(i,j)\mid -N\le i,j\le N\}
$$

with source shifts

$$
\mathbf s_{ij}=iL_1\mathbf e_1+jL_2\mathbf e_2.
$$

Near lists use the fixed P0 triangle kernel directly. Well-separated node pairs use M2L derivatives summed over
the same image shifts. Denote this part by $K_\mathrm{shell}(N)$.

With `fields.periodic.backend="finite_images"`, $K_\mathrm{shell}(N)$ is the complete result. Images outside the range are omitted,
not approximated implicitly.

## Use Ewald2P as a teacher, not a runtime kernel

For a Coulomb sum periodic in two directions and open in one, split the kernel with numerical parameter $\alpha$:

$$
\frac{1}{r}=
\frac{\operatorname{erfc}(\alpha r)}{r}
+\frac{\operatorname{erf}(\alpha r)}{r}.
$$

The first term is evaluated as a real-space image sum. The second is represented by reciprocal modes in the two periodic
directions while retaining the open coordinate in real space.

| Component | Numerical role |
| --- | --- |
| Real space | Converge the short-range part as a screened Coulomb sum |
| Reciprocal `k\ne0` | Add smooth laterally varying modes |
| `k=0` | Plane-average field, treated separately from the physical z-boundary condition |

$\alpha$ distributes work between real and reciprocal sums; it is not Debye screening. The implemented teacher is a
high-accuracy finite sum truncated at `fields.periodic.ewald_layers`. It does not directly evaluate the theoretical infinite sum
for every target.

Define the smooth residual after subtracting the finite shell as

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N).
$$

Ewald evaluation occurs only on the cold path that constructs the operator, never during ordinary particle evaluation.

## Fit a source-node-multipole-to-target-local operator

With fixed geometry, the far residual is linear in source charge. For multipole coefficients $\mathbf M_s$ of source-tree nodes
$s$ at depth 2 and far local coefficients $\mathbf L_t^\mathrm{far}$ at target anchor $t$ (target-tree depth 2),

$$
\mathbf L_t^\mathrm{far}=\sum_s\mathbf A_{t,s}\mathbf M_s.
$$

Seen through a single root multipole, a particle layer whose height approaches the distance to the excluded images (about $2L$)
barely converges, and in a charged three-layer regolith the nonzero modes leaked a plane mean of several volts. Splitting the
sources into depth-2 nodes, and the targets into smaller anchors, makes the expansions converge at the same order.
The cold build constructs $\mathbf A_{t,s}$ as follows:

1. For each source node, place proxy points around the node and form the proxy-to-multipole matrix and its minimum-norm
   pseudoinverse.
2. Place check points for each target anchor and factor the local-field matrix with regularized QR (shared by all source nodes).
3. For each source node, evaluate each proxy source's Ewald residual at the check points (including the plane-mean correction
   below) and solve for the local coefficients.
4. Compose proxy-to-local with the pseudoinverse.
5. Build a fingerprint from geometry, source nodes, target anchors, and configuration, then publish the operator with a checksum.

There are (target anchors)×(source nodes) blocks, up to 64×64 in the three-layer regolith, but the run-time addition is only
small matrix-vector products. The first build is longer (about 27 minutes for the three-layer regolith on one 28-thread process;
MPI distributes the target anchors over ranks) and is reused from the cache.

`cached_kneq0` stores this linear map in a versioned cache and reuses it only when fingerprint, shape, and checksum all match.
It applies the proxy-point operator to the triangle-averaged P2M coefficients
of the fixed P0 triangle source discretization.

A field-only fit cannot determine the constant-potential coefficient of a local expansion, so `cached_kneq0` fits it separately
from potential residuals at the same check points. This avoids mixing field and potential, which have different units, in one
least-squares column set.

Multipole rows have different powers of length. Before the pseudoinverse QR factorization, each row and its
constraint right-hand side are scaled together. This prevents an absolute rank threshold from discarding higher
moments in micrometre SI geometry. The cache fingerprint includes the generator version, so a version with a changed
construction does not reuse older caches.

The physical mean field retains its height distribution inside the source slab. Above all sources, cumulative charge
is set exactly to the total charge $Q$, preventing cancellation residue in inclined-panel polynomials from being
extrapolated into the vacuum. This does not remove the height-dependent mean field.

## Apply correction after ordinary M2L and before L2L

Each batch `update_state` follows this sequence:

```text
P2M
  -> M2M
  -> ordinary tree M2L
  -> L_far(t) += sum_s A_ts M_s  periodic far correction
  -> L2L
  -> L2P + near Direct
```

For every target anchor node, the implementation performs

$$
\mathbf L_t\leftarrow\mathbf L_t+
\sum_s\mathbf A_{t,s}\mathbf M_s.
$$

Ordinary L2L then propagates the corrected locals into leaves, and L2P evaluates them at particle positions. The warm-path cost
is therefore a matrix-vector product and ordinary local evaluation, not an all-source Ewald sum.

The tree's M2L pair cache and `m2l_deriv` remain unchanged; correction is an additional stage between ordinary M2L and L2L.

## Remove symmetric `k=0` from the cached operator

The full Ewald residual contains a symmetric-vacuum `k=0` term needed to define the fit. The cell average of the finite image
shell $K_\mathrm{shell}(N)$ is the infinite-sheet mean field $K_0^\mathrm{sym}$ ($-2\pi|h|/A$ per unit charge at height $h$) plus
a smooth remainder $s(h)$. $s$ is the closed-form integral of $1/r$ over the image block $[c-(N+\tfrac12)L,\,c+(N+\tfrac12)L]^2$
minus $K_0^\mathrm{sym}$, analytic even at $h=0$ ($c$ is the centre of the source-tree root). The plane mean of the residual is
exactly $-s$, so the fit receives only the mean-free $R+s$, and run time subtracts $S(z)=\sum_i q_i s_i(z)$ from a per-source
height table. The FMM-side nonzero kernel is

$$
K_{k\ne0}=K_\mathrm{shell}(N)
+\mathcal F\left[R_\mathrm{Ewald}^{\mathrm{full}}+s\right]-K_0^\mathrm{sym}-S,
$$

where $\mathcal F$ is the operator fit above. The height table of $S$ uses cubic Hermite interpolation at spacing $L/32$; the
field is the derivative of the same polynomial, so it is the gradient of the potential. Field composition adds the physical mean
field exactly once:

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}.
$$

`zero_mode_policy="exclude_k0"` prevents double counting between the FMM backend and physical zero mode; it does not discard
the mean field. See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for the boundary models and Gauss-law construction.

## Understand cache and parallelization boundaries

| Phase | Main work |
| --- | --- |
| Cold miss | Ewald teacher, QR fit, operator construction, atomic publication with checksum |
| Warm hit | Validate fingerprint, shape, and checksum; load operator |
| Charge refresh | P2M/M2M/M2L, apply stored operator, update zero-mode state |
| Particle evaluation | Local expansion plus near Direct; no Ewald reevaluation |

MPI root owns cache I/O and locking. During a cold build, target-operator slices are distributed across MPI ranks, proxy columns
within a target are evaluated with OpenMP, and the completed operator is collected across ranks. Charge history and particle
positions are not stored in the cache.

`cached_kneq0` requires targets to remain inside the configured target box and does not provide an out-of-box Direct fallback
because the cached operator is a production kernel for a fixed target topology.

## Check convergence and cache diagnostics

1. Increase `fields.periodic.image_layers` and verify that the split between near shell and operator does not change results.
2. Vary `fields.periodic.ewald_layers` and, when needed, `fields.periodic.ewald_alpha` to check teacher truncation.
3. Require cold-build and warm-cache fields, potentials, and primary observables to agree to roundoff.
4. Inspect `periodic2_cache_fingerprint`, `periodic2_cache_hit`, and `periodic2_operator_build_count`.
5. Check the selected physical `k=0` and Gauss residual separately from the nonzero mode.
6. Check that the nonzero modes have zero plane mean. With an outer-sheath closure, the mean top-face potential in
   `top_reference_history.csv` equals the outer-root $\phi_0$. `make test-fortran-far-correction` checks the same property for a
   uniform sheet.

Length-scale invariance and cache reuse alone do not establish an accuracy that can be compared with a weak external sheath
field. With the charges of a charged three-layer regolith (one image layer, expansion order 4), a direct Fourier sum in the vacuum
above the layers gives a plane mean of the nonzero modes of 0.02 V or less, a lateral potential error of 0.01–0.03 V, and a
field error of about 2×10³ V/m.

See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for the meaning of finite images
and for the infinite-periodic and zero-mode configuration.

## Code reference

- Mode selection and periodic helpers: [`bem_coulomb_fmm_periodic.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90)
- Ewald2P teacher: [`bem_coulomb_fmm_periodic_ewald.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90)
- Operator construction (source-node and target-anchor pairs): [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
- Smooth remainder $s$ of the image-shell plane mean and its height table: [`bem_periodic_shell_mean.f90`](../src/physics/field_solver/periodic/bem_periodic_shell_mean.f90)
- Cache I/O: [`bem_coulomb_fmm_periodic_cache.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_cache.f90)
- Operator application to state: [`bem_coulomb_fmm_state_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90)
