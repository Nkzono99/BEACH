title: Validating Simulation Results

Lang: [English](ValidationGuide.en.md) | [日本語](ValidationGuide.md)

# Validating Simulation Results

A zero exit status does not establish physical or numerical validity.

## Level 1: execution completed

- exit status is zero;
- `summary.txt`, `charges.csv`, and required histories exist;
- `batches == sim.batch_count`;
- `beachx inspect` can read the output;
- restart model, mesh, and species fingerprints match.

## Level 2: numerical qualification is established

- absorbed, escaped, and max-step populations are physically explainable;
- charge-ledger residual and unresolved discard are checked separately;
- history length supports the claimed steady value; `tol_rel` is monitoring only;
- key observables are stable under smaller `sim.dt`;
- mesh, FMM order/tolerance, outer grid, and sampling are refined as applicable;
- `batch_duration` 0.5x/2x checks do not change the conclusion;
- stochastic conclusions include seed or ensemble sensitivity.

Qualification here means satisfying declared discretization and convergence
criteria. Process completion, generated CSV files, or one
`status="converged"` value is not sufficient by itself.

## Level 3: a physical conclusion is supported

- all input differences between compared cases, except the intended physical
  model change, are enumerated;
- boundary conditions, self interaction, surface trace, and source/target
  motion policy are explicit;
- finite-box, finite-time, and finite-shell results are not extrapolated to
  infinity, steady state, or infinite periodicity without evidence;
- numerical, stochastic, and model uncertainty are reflected in the claim.

## Model-Specific Checks

| Model | Required diagnostics |
| --- | --- |
| periodic2 cached | cache fingerprint, cold/warm agreement, zero-mode/Gauss residual |
| unified outer | accessibility refinement, linearity, outer energy/frozen-field error |
| kinetic outer | solver status, Poisson residual, Bohm/branch applicability |
| photoelectron | emission/return ledger, ambient charge ratio, histogram range |
| object detachment | primary-only self exclusion, PV trace, work/potential agreement, quadrature, finite-shell/cache, from-rest barrier |

## Additional Gates for Periodic Object Detachment

1. Record whether the conclusion uses `configured` or `infinite_physical`.
   `configured` reproduces the run's finite or cached policy;
   `infinite_physical` combines cached `k != 0` with the physical
   `E_bottom=0` zero mode.
2. Confirm `exclude_primary_keep_images`. Do not confuse it with the legacy
   `kernel-forces` policy `exclude_target_lattice` or the potential
   reconstruction self term `area_equivalent`.
3. For triangle sources, check order-3/order-7 quadrature or mesh refinement.
   Object mechanics uses the PV trace, while the particle pusher uses the
   one-sided plus trace. Cached `cached_kneq0_trace_correction` metadata has
   already been applied to `periodic_kneq0`; never add it again.
4. On a frozen-source path, require `integral(F_z dh)` and
   `U_env(0)-U_env(h)` to agree within the declared tolerance and require
   `path.status="converged"`. The frozen external energy
   `U_env=sum(q phi_env)` has no factor of `1/2`.
5. Use the native canonical-unwrapped source representation for finite shells
   and retain raw symmetric and `E_bottom=0`-corrected results. Two small raw
   adjacent-shell increments produced false convergence, so selection now
   combines `force_tail_proxy_N` / `work_tail_proxy_J` with
   `reference_force_error_N` / `reference_work_error_J` for an
   `infinite_physical` reference. `increment_converged` is this combined gate
   and must be true twice consecutively. `status="not_converged"` correctly
   carries `selected_image_layers=None` and `selected_path=None`.
6. For `evaluate_release()`, check the full `barrier_free_from_rest` path with
   finite-range adhesion and gravity, not endpoint energy alone.
7. A non-neutral x/y-periodic cell can retain a constant far field and linear
   potential. Do not report a finite-height work or speed as escape energy or
   escape speed at infinity.
8. Also qualify the cached infinite-periodic implementation with opt-in
   analytic oracles. For a uniform non-neutral triangle plane
   (`sigma=Q/A`) under `E_bottom=0`, require `E_z=0` below,
   `E_z=sigma/epsilon0` above, the surface PV
   `E_z=sigma/(2 epsilon0)`, surface pressure
   `sigma^2/(2 epsilon0)`, and total object force
   `Q^2/(2 epsilon0 A)`. For a neutral `sigma_0 cos(kx)` sheet, require the
   `k != 0` field and potential amplitudes to decay as
   `exp(-k |z-z0|)` and the analytic error to decrease under triangle-mesh
   refinement.

A zero CLI exit code and four output artifacts establish Level 1. Path and
shell convergence are only parts of Level 2. The model-selection and
sensitivity checks above are needed before considering a Level-3 physical
claim.

See [Physics release verification](PhysicsReleaseVerification.en.html) for small reference contracts.
