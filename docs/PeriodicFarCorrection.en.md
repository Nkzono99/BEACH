title: periodic2 Far Correction

Lang: [English](PeriodicFarCorrection.en.md) | [日本語](PeriodicFarCorrection.md)

# periodic2 Far Correction

The `periodic2` FMM evaluates the primary cell and finite images with ordinary tree FMM, then corrects the smooth
infinite-periodic field outside that range with an additional operator. This page is the canonical reference for the Ewald2P
teacher, `cached_kneq0`, cache, and connection to FMM state.

## Select a finite-image model or an infinite-periodic approximation

```toml
[sim]
field_solver = "fmm"
field_periodic_image_layers = 1
field_periodic_far_correction = "cached_kneq0"
field_periodic_cache_dir = ".beach_cache/periodic2"
field_periodic_generation_tolerance = 1.0e-8

[domain]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0, 1.0, 1.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"
```

`[domain]` owns the fixed target topology and periodic axes; `[field_boundary]` selects the periodic2 field. The current
`cached_kneq0` path accepts only targets inside an x/y-periodic, z-nonperiodic domain box.

| `field_periodic_far_correction` | Computation | Use |
| --- | --- | --- |
| `none` | Primary cell and finite images only | Finite-image model and convergence comparisons |
| `auto` | Currently normalized to `none` | Compatibility only |
| `cached_kneq0` | Build and reuse a versioned operator | Production infinite-periodic nonzero modes |

`cached_kneq0` does not determine the z-directed plane-average `k=0` field. Field composition adds the physical `k=0` exactly
once as described in [periodic2 electrostatics](PeriodicElectrostatics.en.html). Self-consistent outer-plasma/sheath models are
unsupported.
`m2l_root_oracle` has been removed and is rejected at startup with guidance to use `cached_kneq0`.

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

With `field_periodic_far_correction="none"`, $K_\mathrm{shell}(N)$ is the complete result. Images outside the range are omitted,
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
high-accuracy finite sum truncated at `field_periodic_ewald_layers`. It does not directly evaluate the theoretical infinite sum
for every target.

Define the smooth residual after subtracting the finite shell as

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N).
$$

Ewald evaluation occurs only on the cold path that constructs the operator, never during ordinary particle evaluation.

## Fit a root-multipole-to-target-local operator

With fixed geometry, the far residual is linear in source charge. For root multipole $\mathbf M_\mathrm{root}$ and far local
coefficients $\mathbf L_t^\mathrm{far}$ at target anchor $t$,

$$
\mathbf L_t^\mathrm{far}=\mathbf A_t\mathbf M_\mathrm{root}.
$$

The cold build constructs $\mathbf A_t$ as follows:

1. Place proxy points around the source root and form the proxy-to-root-multipole matrix.
2. Place check points for each target anchor and evaluate each proxy source's Ewald residual.
3. Solve a regularized QR problem for local coefficients reproducing the check-point field.
4. Compose proxy-to-local with the minimum-norm pseudoinverse of proxy-to-multipole.
5. Build a fingerprint from geometry and configuration, then publish the operator with a checksum.

`cached_kneq0` stores this linear map in a versioned cache and reuses it only when fingerprint, shape, and checksum all match.
It applies the proxy-point operator to the triangle-averaged P2M coefficients
of the fixed P0 triangle source discretization.

A field-only fit cannot determine the constant-potential coefficient of a local expansion, so `cached_kneq0` fits it separately
from potential residuals at the same check points. This avoids mixing field and potential, which have different units, in one
least-squares column set.

## Apply correction after ordinary M2L and before L2L

Each batch `update_state` follows this sequence:

```text
P2M
  -> M2M
  -> ordinary tree M2L
  -> L_far(t) += A_t M_root      periodic far correction
  -> L2L
  -> L2P + near Direct
```

For every target anchor node, the implementation performs

$$
\mathbf L_t\leftarrow\mathbf L_t+
\mathbf A_t\mathbf M_\mathrm{root}.
$$

Ordinary L2L then propagates the corrected locals into leaves, and L2P evaluates them at particle positions. The warm-path cost
is therefore a matrix-vector product and ordinary local evaluation, not an all-source Ewald sum.

The tree's M2L pair cache and `m2l_deriv` remain unchanged; correction is an additional stage between ordinary M2L and L2L.

## Remove symmetric `k=0` from the cached operator

The full Ewald residual contains a symmetric-vacuum `k=0` term needed to define the fit. Because `lower_boundary_model` selects
the actual mean field, the FMM-side nonzero kernel is

$$
K_{k\ne0}=K_\mathrm{shell}(N)
+R_\mathrm{Ewald}^{\mathrm{full}}-K_0^\mathrm{sym},
$$

and field composition adds the physical mean field exactly once:

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

1. Increase `field_periodic_image_layers` and verify that the split between near shell and operator does not change results.
2. Vary `field_periodic_ewald_layers` and, when needed, `field_periodic_ewald_alpha` to check teacher truncation.
3. Require cold-build and warm-cache fields, potentials, and primary observables to agree to roundoff.
4. Inspect `periodic2_cache_fingerprint`, `periodic2_cache_hit`, and `periodic2_operator_build_count`.
5. Check the selected physical `k=0` and Gauss residual separately from the nonzero mode.

See [Finite Periodic Configuration](FinitePeriodicConfiguration.en.html) for
finite-image use and [periodic2 electrostatics](PeriodicElectrostatics.en.html)
for the complete infinite-periodic and zero-mode configuration.

## Code reference

- Mode selection and periodic helpers: [`bem_coulomb_fmm_periodic.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90)
- Ewald2P teacher: [`bem_coulomb_fmm_periodic_ewald.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90)
- Root-operator construction: [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
- Cache I/O: [`bem_coulomb_fmm_periodic_cache.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_cache.f90)
- Operator application to state: [`bem_coulomb_fmm_state_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90)
