# Periodic Object Detachment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reusable BEACH Python API that computes the instantaneous 3D wrench and frozen-charge detachment path of one central object while retaining its periodic replicas, then validate the result against the archived finite-image simulation and paired SysA finite/infinite simulations.

**Architecture:** A native C ABI exposes the simulator's physical periodic zero mode and the existing FMM cache controls without changing existing ABI entry points. An immutable Python `ObjectInteractionSnapshot` composes the same Fortran FMM nonzero field, the same native physical `k=0` provider, an exact free-primary evaluation through the same Fortran core, and the configured uniform field; a separate `ObjectProbe` moves target quadrature only. Mechanics, finite-shell convergence, CLI/report generation, and SysA simulation staging remain separate modules so electrostatic policy cannot be changed by plotting or release assumptions.

**Tech Stack:** Fortran 2008, ISO C binding, fpm, OpenMP/MPI, Python 3.10+, ctypes, NumPy, pytest, Ruff, Markdown, KUDPC Slurm/SysA/SysB.

## Global Constraints

- Run test controllers and analysis payloads only on KUDPC compute nodes through `tssrun` or `sbatch`; never run them directly on `laurel31`.
- Do not run multiple `fpm test` commands concurrently in the shared worktree.
- The primary electric-field result must call the same `bem_coulomb_fmm_core` used by the simulator; Python direct sums are only small-fixture or legacy-compatibility oracles.
- Infinite periodic results are complete only as `cached_kneq0 + physical k0(E_bottom=0) + sim.e0`; reject an unlabeled `k!=0`-only result.
- Exclude only the selected primary-cell object and retain all of its nonzero periodic replicas.
- Keep the original source snapshot immutable; move target probe geometry only.
- Preserve the existing `calc_object_forces_kernel` public behavior as the explicit `exclude_target_lattice` compatibility policy.
- Never silently substitute point-centroid and triangle-P0 source or target models.
- Use principal-value surface traces for mechanical target force; simulator pusher `plus` traces remain unchanged.
- First release supports `outer_plasma.model="none"`; reject active outer closures rather than dropping their field.
- Return SI units and immutable/read-only public result records.
- Preserve test-tier policy: analytic tests in L1/L2, heavy FMM in L3, cached/far-correction diagnostics opt-in only.
- Store caches, staged simulations, and generated validation artifacts under `/LARGE1/gr20001/b36291/codex-tmp/`; do not commit them.
- Use TDD: every behavior change starts with a focused failing test and ends with a focused green test before commit.
- Full paired simulations retain `rng_seed=12345`, six MPI ranks, 112 OpenMP threads per rank, and all archived physical parameters except the explicit periodic far-correction/cache difference.

---

### Task 1: Native Physical Zero-Mode C ABI

**Files:**
- Create: `src/physics/periodic_zero_mode/bem_periodic_zero_mode_c.f90`
- Create: `tests/fortran/test_periodic_zero_mode_c.f90`
- Modify: `fpm.toml`
- Modify: `Makefile`

**Interfaces:**
- Produces C symbols `beach_zero_mode_create`, `beach_zero_mode_destroy`, `beach_zero_mode_build`, `beach_zero_mode_update`, and `beach_zero_mode_eval`.
- `beach_zero_mode_build(handle, nsrc, source_heights_3xn, area_xy)` accepts three vertex heights per source; a point source repeats its centroid height three times.
- `beach_zero_mode_update(handle, nsrc, charge, e_bottom, z_gauge, phi_gauge)` refreshes the existing simulator state.
- `beach_zero_mode_eval(handle, ntarget, z, trace, phi, ez)` accepts trace `-1`, `0`, or `1` using the existing minus/principal-value/plus constants.
- Status codes mirror the field-kernel convention: `0=ok`, `1=invalid handle`, `2=invalid argument`, `3=not ready`.

- [ ] **Step 1: Write the failing C-ABI test**

  Add a two-sheet non-neutral fixture that calls the five new symbols and compares every target against direct calls to `build_periodic_zero_mode_height_plan`, `refresh_periodic_zero_mode_state`, and `eval_periodic_zero_mode`:

  ```fortran
  source_heights(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  source_heights(:, 2) = [1.0d0, 1.0d0, 1.0d0]
  charge = [2.0d-9, -0.5d-9]
  z = [-0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0]
  status = beach_zero_mode_build(handle, 2_c_int, c_loc(source_heights), 4.0d0)
  status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), 0.0d0, -0.5d0, 0.0d0)
  status = beach_zero_mode_eval(handle, 5_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
  ```

  Add an inclined source with three distinct vertex heights and compare its
  exact height projection as well. Assert invalid area, invalid trace,
  update-before-build, size mismatch, eval-before-update, NULL pointers,
  overflowing counts, and every NaN/Infinity in source heights, charges,
  `e_bottom`, gauges, and target heights return status rather than terminating.

- [ ] **Step 2: Run RED on an available SysB compute node**

  First inspect `hostname`, `module list`, `spartition`, and `qgroup`, then choose
  a visible SysB partition. For example, when `eb` is visible, run:

  ```bash
  tssrun -p eb -t 0:20:0 --rsc p=1:t=2:c=2 \
    bash -lc 'cd /LARGE1/gr20001/b36291/codex-tmp/BEACH-periodic-object-force && \
    FPM_ACTION=test ./build.sh --target test_periodic_zero_mode_c'
  ```

  Expected: fpm reports the missing test target/module or undefined C symbols.

- [ ] **Step 3: Implement the native handle without duplicating zero-mode physics**

  Define a private handle that owns `periodic_zero_mode_plan_type`,
  `periodic_zero_mode_state_type`, `built`, and `charged`. Before delegating,
  validate counts against the addressable Fortran shape, pointer presence,
  finite `area_xy`, every source height/charge, `e_bottom`, both gauges, every
  target height, and trace membership. This prevents the delegated routines'
  `error stop` paths from crossing the ABI. Delegate all numerical work:

  ```fortran
  call build_periodic_zero_mode_height_plan(real(source_heights, dp), real(area_xy, dp), &
                                             zero%plan, plan_status, message)
  call refresh_periodic_zero_mode_state(zero%plan, real(charge, dp), real(e_bottom, dp), &
                                        real(z_gauge, dp), real(phi_gauge, dp), zero%state)
  call eval_periodic_zero_mode(zero%plan, zero%state, real(z(i), dp), int(trace, i32), &
                               phi_value, field_value)
  ```

  Do not reimplement cumulative-charge polynomials in the C wrapper.

- [ ] **Step 4: Run GREEN and existing zero-mode regressions**

  Run the new target and `test_periodic_zero_mode` sequentially in one `tssrun`
  controller. Add `test_periodic_zero_mode_c` explicitly to
  `FORTRAN_L2_TARGETS`; expected: both pass with no `error stop` on invalid ABI
  inputs.

- [ ] **Step 5: Commit**

  ```bash
  git add src/physics/periodic_zero_mode/bem_periodic_zero_mode_c.f90 \
          tests/fortran/test_periodic_zero_mode_c.f90 fpm.toml Makefile
  git commit -m "feat: expose physical periodic zero mode"
  ```

### Task 2: Field-Kernel Cache Controls And Diagnostics

**Files:**
- Modify: `src/physics/field_solver/bem_field_kernel_c.f90`
- Modify: `tests/fortran/test_field_kernel_c.f90`
- Create: `tests/fortran/test_field_kernel_cache_c.f90`
- Modify: `fpm.toml`
- Modify: `Makefile`
- Modify: `beach/fortran_results/kernel.py`
- Modify: `beach/fortran_results/potential.py`
- Modify: `tests/python/test_field_kernel.py`
- Create: `tests/python/test_field_kernel_cache.py`
- Modify: `tests/python/test_far_correction_contract.py`

**Interfaces:**
- Adds ABI-safe setter `beach_kernel_set_periodic_cache(handle, path_ptr, path_len, tolerance)`; existing build function signatures remain unchanged.
- Adds `beach_kernel_get_periodic_cache_info(handle, hit_ptr, build_count_ptr,
  fingerprint_ptr, fingerprint_capacity, fingerprint_length_ptr, path_ptr,
  path_capacity, path_length_ptr)`. Text is UTF-8 in caller-owned buffers;
  returned lengths exclude the NUL byte and always report required length.
  Capacities must include the NUL byte; short buffers return invalid-argument
  without truncation. Not-yet-built handles return not-ready. A successfully
  built non-cached plan returns hit/count zero, lengths zero, and NUL-only text.
- Extends `FieldKernelOptions` with `periodic_cache_dir: str = ".beach_cache/periodic2"` and `periodic_generation_tolerance: float = 1.0e-8`.
- Adds immutable
  `FieldKernelDiagnostics(periodic_cache_hit: bool | None,
  periodic_operator_build_count: int,
  periodic_cache_fingerprint: str | None,
  periodic_cache_path: Path | None)` and `FieldKernel.diagnostics()`.
- `_options_from_result` reads `sim.field_periodic_cache_dir` and `sim.field_periodic_generation_tolerance`.
- `_coerce_periodic2` gains keyword `allow_cached_kneq0=False`, validates both
  mapping and existing 7-tuple inputs, and propagates the keyword through auto
  result/config resolution. Pure-Python potential/Coulomb paths stay fail-closed;
  only the native field-kernel resolver passes `allow_cached_kneq0=True`.
- New symbols are configured lazily. Loading an older shared library for legacy free/finite calls remains supported; selecting `cached_kneq0`, native zero mode, direct-eval methods, or diagnostics without the required symbol raises a model-specific `FieldKernelError`.

- [ ] **Step 1: Write failing Fortran and Python tests**

  Lightweight Fortran ABI tests in `test_field_kernel_c` cover setter/getter
  argument validation and non-cached/not-ready behavior without generating an
  operator. The separate opt-in `test_field_kernel_cache_c` configures a
  temporary cache directory before `beach_kernel_build`, builds the same cached
  fixture twice, and requires first build count `1`, second cache hit true,
  identical fingerprint/path, and invalid empty path/tolerance rejection.

  The default-tier `test_field_kernel.py` uses a fake ABI to verify option
  marshaling, lazy-symbol errors, tuple validation, and diagnostics decoding
  without constructing an operator. The real cache build is isolated in
  opt-in `test_field_kernel_cache.py`:

  ```python
  options = FieldKernelOptions(
      periodic2=((0, 1), (2.0, 2.0), (0.0, 0.0), 0, "cached_kneq0", 0.0, 4),
      box_min=(0.0, 0.0, -1.0),
      box_max=(2.0, 2.0, 1.0),
      periodic_cache_dir=str(tmp_path / "cache"),
      periodic_generation_tolerance=2.5e-9,
  )
  with FieldKernel(source_pos, source_q, options=options, library_path=lib) as kernel:
      diagnostics = kernel.diagnostics()
  assert diagnostics.periodic_cache_path.parent == tmp_path / "cache"
  assert diagnostics.periodic_cache_fingerprint
  ```

- [ ] **Step 2: Run RED**

  Run lightweight `test_field_kernel_c` and
  `pytest -q tests/python/test_field_kernel.py` through one compute-node
  controller after `make build-kernel`; run `test_field_kernel_cache_c` and
  `tests/python/test_field_kernel_cache.py` only through a dedicated
  far-correction opt-in Make target excluded from default pytest/L1 collection.
  Expected: missing fields/symbols.

- [ ] **Step 3: Implement cache configuration and safe C-string copies**

  Store cache settings in `field_kernel_handle`; initialize defaults in create;
  copy them into each `fmm_options_type` before point/panel build. Setter strings
  are length-delimited UTF-8 and reject embedded NUL, empty, invalid UTF-8, and
  over-256-byte paths. Getter strings follow the capacity/required-length/NUL
  contract above. In Python, call the setter before geometry build and configure
  ctypes signatures once; convert empty diagnostics to `None`, never `Path(".")`.

- [ ] **Step 4: Run GREEN and ABI compatibility tests**

  Run the lightweight focused Fortran/Python tests, then the dedicated opt-in
  cache target, and rebuild/load the shared library using unchanged legacy
  calls. Expected: old tests pass without supplying cache options and no real
  operator is generated by default L1/L2.

- [ ] **Step 5: Commit**

  ```bash
  git add src/physics/field_solver/bem_field_kernel_c.f90 \
          tests/fortran/test_field_kernel_c.f90 tests/fortran/test_field_kernel_cache_c.f90 \
          fpm.toml Makefile beach/fortran_results/kernel.py \
          beach/fortran_results/potential.py tests/python/test_field_kernel.py \
          tests/python/test_field_kernel_cache.py \
          tests/python/test_far_correction_contract.py
  git commit -m "feat: expose periodic kernel cache controls"
  ```

### Task 3: Exact Direct Evaluation Through The Fortran Core

**Files:**
- Modify: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`
- Modify: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90`
- Modify: `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90`
- Modify: `src/physics/field_solver/bem_field_kernel_c.f90`
- Modify: `tests/fortran/test_field_kernel_c.f90`
- Modify: `beach/fortran_results/kernel.py`
- Modify: `tests/python/test_field_kernel.py`

**Interfaces:**
- Adds public core procedures `eval_direct_points` and `eval_direct_potential_points` that bypass tree/FMM far approximations and sum every configured source/image with the existing point or analytic panel kernels.
- Adds C symbols `beach_kernel_eval_e_direct` and `beach_kernel_eval_phi_direct` without changing existing eval symbols.
- Adds `FieldKernel.eval_e_direct(points)` and `FieldKernel.eval_phi_direct(points)`.
- Point evaluation uses the plan's configured softening and coincident-source convention. Panel evaluation calls `panel_potential_field` with `panel_side_principal_value`.
- Primary subtraction builds a non-periodic target-only plan and calls these exact-direct methods; it never assumes a large `leaf_max` makes an ordinary FMM evaluation exact.
- Exact-direct public methods reject every periodic plan. The existing internal
  direct-all-source routines enumerate only configured finite shifts and do not
  implement target wrapping or cached far correction, so these methods are
  deliberately restricted to the non-periodic primary evaluator.

- [ ] **Step 1: Write failing exact-direct tests**

  Add a point fixture deliberately configured with an aggressive FMM
  approximation and prove direct entry points match an explicit all-source
  softened sum while ordinary FMM may differ. Add an off-surface and on-surface
  panel fixture and compare field/potential against `panel_potential_field` with
  principal-value side. Python calls both new methods, checks charge refresh,
  and requires a periodic plan to fail closed.

- [ ] **Step 2: Run RED**

  Run `test_field_kernel_c` and the focused Python field-kernel tests on `eb`. Expected: undefined direct-eval procedures/symbols.

- [ ] **Step 3: Expose existing direct-all-source internals, do not duplicate kernels**

  Promote batched wrappers around the existing
  `eval_direct_all_sources_scalar` and
  `eval_direct_all_sources_potential_scalar`:

  ```fortran
  !$omp parallel do default(none) schedule(static) &
  !$omp   shared(plan, state, target_pos, field, ntarget) private(i)
  do i = 1_i32, ntarget
    call eval_direct_all_sources_scalar( &
      plan, state, target_pos(1, i), target_pos(2, i), target_pos(3, i), &
      plan%options%softening**2, field(1, i), field(2, i), field(3, i) &
    )
  end do
  !$omp end parallel do
  ```

  Keep Coulomb scaling in the same C ABI layer as ordinary `eval_e`/`eval_phi`.

- [ ] **Step 4: Run GREEN and core regressions**

  Run `test_field_kernel_c`, `test_coulomb_fmm_core_basic`,
  `test_coulomb_fmm_core_panel`, and `test_field_kernel.py` sequentially on a
  compute node. Expected: direct values match exact oracles and existing FMM
  results are unchanged.

- [ ] **Step 5: Commit**

  ```bash
  git add src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90 \
          src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90 \
          src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90 \
          src/physics/field_solver/bem_field_kernel_c.f90 \
          tests/fortran/test_field_kernel_c.f90 beach/fortran_results/kernel.py \
          tests/python/test_field_kernel.py
  git commit -m "feat: expose exact direct field evaluation"
  ```

### Task 4: Python Zero-Mode Binding And Mechanical Records

**Files:**
- Create: `beach/fortran_results/periodic_zero_mode.py`
- Create: `beach/fortran_results/detachment.py`
- Create: `tests/python/test_periodic_zero_mode.py`
- Create: `tests/python/test_detachment.py`
- Modify: `beach/fortran_results/__init__.py`
- Modify: `beach/__init__.py`

**Interfaces:**
- Produces `PeriodicZeroMode(source_heights_m, source_charges_C, area_xy_m2, *, e_bottom_V_m=0.0, z_gauge_m=None, phi_gauge_V=0.0, library_path=None)` with `update_charges`, `eval(z_m, trace="principal_value")`, `close`, and context-manager support.
- Produces immutable/read-only `WrenchComponent`, `ObjectWrench`, `ObjectForcePath`, `DetachmentResult`, and `AdhesionProfile`.
- `AdhesionProfile.none()`, `.finite_range_constant(force_N, range_m)`, and `.tabulated(displacement_m, force_N)` expose vectorized `force_N(h)` and `work_J(h)`.
- `ObjectForcePath.from_samples(displacement_m, force_N, torque_Nm, potential_energy_J=None)` constructs deterministic test/adapter paths.
- `ObjectForcePath.evaluate_release(mass_kg, gravity_m_s2, adhesion,
  eta_translation=1.0, energy_tolerance_J=1.0e-18,
  dissipation_work_J=None)` returns all work/speed/barrier arrays and checks the
  continuous piecewise model, not only the caller's sample points.

- [ ] **Step 1: Write failing wrapper and mechanics tests**

  Verify the wrapper matches the analytic two-sheet field for principal/plus traces and that `update_charges` changes the result without rebuilding. Verify mechanics with a constant `2 N` force over `[0, 1, 2] m`, `m=1 kg`, `g=1 m/s2`, and finite adhesion `0.25 N` through `0.5 m`:

  ```python
  path = ObjectForcePath.from_samples(
      displacement_m=np.array([0.0, 1.0, 2.0]),
      force_N=np.array([[0.0, 0.0, 2.0]] * 3),
      torque_Nm=np.zeros((3, 3)),
  )
  result = path.evaluate_release(
      mass_kg=1.0,
      gravity_m_s2=1.0,
      adhesion=AdhesionProfile.finite_range_constant(0.25, 0.5),
  )
  np.testing.assert_allclose(result.electrostatic_work_J, [0.0, 2.0, 4.0])
  assert result.barrier_free_from_rest
  ```

  Add negative-initial-work, positive-endpoint-with-barrier,
  `eta_translation=0.5`, invalid tables, and read-only array tests. Include a
  finite-range adhesion endpoint and a tabulated-force knot strictly between
  path samples where endpoint energies are positive but the intervening energy
  minimum is negative; both must report a barrier.

- [ ] **Step 2: Run RED**

  Build the shared kernel, then run the two Python files on `eb`. Expected: import failures for new classes.

- [ ] **Step 3: Implement thin binding and mechanics-only calculations**

  `PeriodicZeroMode` must only marshal arrays/status; all zero-mode values come
  from Task 1. Use cumulative trapezoids for sampled work and exact integration
  of the piecewise-linear tabulated adhesion force. Form a mechanics grid from
  the union of force-path samples, finite-range endpoints, adhesion table knots,
  and dissipation-work knots. On each resulting interval the interpolated net
  force is linear; evaluate interval endpoints and every interior derivative
  root to obtain the true piecewise-model energy minimum and first barrier.
  Never infer barrier freedom from `np.min` on the original displacement samples
  alone. At the public sample points compute:

  ```python
  available = (electric_work - mass_kg * gravity_m_s2 * displacement
               - adhesion.work_J(displacement) - dissipation_work)
  speed = np.sqrt(2.0 * eta_translation * np.maximum(available, 0.0) / mass_kg)
  barrier_free = bool(minimum_continuous_available >= -energy_tolerance_J)
  ```

  Make defensive copies and mark returned arrays non-writeable. Compute maximum
  reachable speed only over the contiguous prefix before the first energy
  barrier; later positive energy does not make a from-rest endpoint reachable.

- [ ] **Step 4: Run GREEN and Ruff**

  Run both focused Python tests and `ruff check beach/fortran_results/periodic_zero_mode.py beach/fortran_results/detachment.py tests/python/test_periodic_zero_mode.py tests/python/test_detachment.py` on a compute node. Expected: all pass.

- [ ] **Step 5: Commit**

  ```bash
  git add beach/fortran_results/periodic_zero_mode.py beach/fortran_results/detachment.py \
          beach/fortran_results/__init__.py beach/__init__.py \
          tests/python/test_periodic_zero_mode.py tests/python/test_detachment.py
  git commit -m "feat: add detachment mechanics records"
  ```

### Task 5: Point-Source Object Interaction Snapshot

**Files:**
- Create: `beach/fortran_results/object_interaction.py`
- Create: `tests/python/test_object_interaction.py`
- Modify: `beach/fortran_results/kernel.py`
- Modify: `beach/fortran_results/__init__.py`
- Modify: `beach/__init__.py`
- Modify: `beach/fortran_results/facade.py`
- Modify: `tests/fortran/test_periodic2_cached_snapshot.f90`

**Interfaces:**
- Produces `ObjectInteractionSnapshot.from_result(result, *, step=-1, config_path=None, periodic_model="configured", cache_dir=None, generation_tolerance=None, library_path=None)`.
- `periodic_model` is exactly `"configured"` or `"infinite_physical"`; the
  latter overrides only far correction to `cached_kneq0`. Any resolved
  `cached_kneq0` field, including a configured one, is completed with the native
  physical zero mode. A configured finite shell (`far_correction="none"`) has no
  separately added zero mode.
- Produces `ObjectInteractionSnapshot.object_probe(mesh_id, *, self_policy="exclude_primary_keep_images", target_integration="auto", quadrature_order=7)`.
- Produces `ObjectProbe.wrench(transform=None, transform_origin="geometric_area_centroid", torque_origin="geometric_area_centroid", components=True) -> ObjectWrench`.
- Physical component keys are exactly `other_objects_all_images`, `target_periodic_images`, `external_uniform`, and `total_external`; numerical metadata also exposes `periodic_kneq0`, `physical_k0`, and `primary_free_subtraction`.
- Snapshot/probe are closeable context managers. A probe transform never mutates source geometry or source charges.
- `field_source_model` must be exactly point-compatible or `triangle_p0`; unknown values raise instead of falling back to point centroids.

- [ ] **Step 1: Write failing analytic point tests**

  Add fixtures for:

  1. one free point object: primary exclusion gives zero force/torque;
  2. two free objects: exact Coulomb action/reaction and torque-origin shift;
  3. one net-neutral target in a symmetric x/y periodic lattice: own-image `Fz=0` at the original position;
  4. a non-neutral two-height periodic fixture: both an
     `infinite_physical` override and `periodic_model="configured"` with an
     already configured `cached_kneq0` field equal the same explicit
     `periodic k!=0 + native physical k0 + uniform` composition;
  5. the old `calc_object_forces_kernel` still removes the full target lattice;
  6. active outer-plasma config is rejected.

  Keep component algebra and immutable-record tests runnable without a shared
  library by injecting internal fake periodic/zero/direct evaluators. Native
  equivalence cases may use `_kernel_lib()` and skip only when the library is
  genuinely absent.

  Extend `test_periodic2_cached_snapshot` with off-surface target points and
  compare the production `electrostatic_snapshot` total field against C-ABI
  `cached_kneq0 + beach_zero_mode_eval + prescribed_e`. This is the binding
  acceptance test required to claim the same Fortran calculation method.

  The key policy assertion is:

  ```python
  with ObjectInteractionSnapshot.from_result(
      result,
      periodic_model="infinite_physical",
      config_path=config,
      cache_dir=tmp_path / "cache",
      library_path=lib,
  ) as snapshot:
      wrench = snapshot.object_probe(1).wrench(components=True)
  np.testing.assert_allclose(
      wrench.components["other_objects_all_images"].force_N
      + wrench.components["target_periodic_images"].force_N
      + wrench.components["external_uniform"].force_N,
      wrench.force_N,
  )
  ```

- [ ] **Step 2: Run RED**

  Run `pytest -q tests/python/test_object_interaction.py` after `make build-kernel` on `eb`. Expected: missing API.

- [ ] **Step 3: Implement source-state composition using one Fortran geometry plan**

  Resolve original centers/triangles/charges/mesh IDs once. Build the periodic `FieldKernel` with `external_e0=(0,0,0)` so uniform field stays separable. Evaluate two periodic charge states without rebuilding geometry:

  ```python
  q_target = np.where(target_mask, charges, 0.0)
  q_other = np.where(target_mask, 0.0, charges)

  periodic.update_charges(q_other)
  e_other_kneq0 = periodic.eval_e(target_points)
  periodic.update_charges(q_target)
  e_target_lattice_kneq0 = periodic.eval_e(target_points)
  ```

  Build the primary free kernel from target sources with periodicity disabled
  and evaluate it through Task 3's explicit exact-direct core entry points.
  Compose:

  ```python
  e_own_images = (
      e_target_lattice_kneq0 + e_target_k0 - e_primary_free
  )
  e_total = e_other_kneq0 + e_other_k0 + e_own_images + external_e0
  ```

  Compose uniform potential with the same explicit gauge as
  `phi_uniform(r) = -external_e0 dot (r - gauge_origin)` so force-integrated and
  potential-difference work remain consistent. Resolve the effective far
  correction before composition: `none` means finite configured images with no
  separate zero mode, whereas every `cached_kneq0` result, whether configured or
  overridden by `infinite_physical`, must add the native physical zero mode.

  Aggregate all forces/torques through one shared NumPy helper that exactly
  implements `sum(qE)` and `sum(r-origin) cross qE`; compare it against
  `beach_kernel_force_on_charges` in tests. Prevalidate finite target points and
  their full box/interface membership before any native call. Serialize charge
  substitutions under an internal lock and restore the full charge state in
  `finally`, including when an evaluator raises.

- [ ] **Step 4: Run GREEN, compatibility, and leak tests**

  Run `test_object_interaction.py`, `test_field_kernel.py`,
  `test_periodic2_cached_snapshot`, and repeated context-manager
  construction/destruction. Expected: simulator/C-ABI total fields agree,
  component sums close, legacy records unchanged, and no use-after-close.

- [ ] **Step 5: Commit**

  ```bash
  git add beach/fortran_results/object_interaction.py beach/fortran_results/kernel.py \
          beach/fortran_results/__init__.py beach/fortran_results/facade.py beach/__init__.py \
          tests/python/test_object_interaction.py tests/fortran/test_periodic2_cached_snapshot.f90
  git commit -m "feat: add periodic object interaction snapshot"
  ```

### Task 6: Triangle-P0 Target Integration And Rigid Wrench Covariance

**Files:**
- Create: `beach/fortran_results/panel_quadrature.py`
- Modify: `beach/fortran_results/object_interaction.py`
- Create: `tests/python/test_object_panel_interaction.py`
- Modify: `tests/python/test_field_kernel.py`

**Interfaces:**
- Produces internal `panel_target_quadrature(triangles_m, element_charges_C, order)` for Gauss-Duffy orders `3` and `7`; one-point behavior is named `target_integration="centroid_compatibility"` rather than called quadrature.
- Returned `(points_m, charge_weights_C, element_index)` satisfies exact per-element charge conservation.
- `target_integration="auto"` selects saved centroid points for point source and Gauss-Duffy order 7 for triangle P0; it never changes the source model.
- Triangle primary free evaluation uses a non-periodic panel `FieldKernel` and Task 3's exact-direct methods with principal-value panel kernels.
- `ObjectProbe.wrench` accepts any existing `RigidTransform`; transform and torque origins accept `"geometric_area_centroid"`, `"origin"`, or a finite absolute 3-vector.

- [ ] **Step 1: Write failing triangle/quadrature tests**

  Add tests that:

  ```python
  points, q_weight, elem = panel_target_quadrature(triangles, q_elem, order=7)
  for i, q in enumerate(q_elem):
      assert np.sum(q_weight[elem == i]) == pytest.approx(q)
  ```

  Verify a uniform external field gives `F=Q E`, a single isolated triangle object gives zero external wrench after primary exclusion, orders 3 and 7 converge on a separated two-panel fixture, and rigid translation/rotation obey:

  ```text
  F' = R F
  tau'_(R o + a) = R tau_o
  tau_b = tau_a - (b-a) cross F
  ```

  Add a near-surface fixture demonstrating centroid compatibility is labeled and differs from converged quadrature.

- [ ] **Step 2: Run RED**

  Run the new Python test on `eb`; expected missing quadrature module and current point-source rejection.

- [ ] **Step 3: Implement Gauss-Duffy target rules and panel dispatch**

  Use `numpy.polynomial.legendre.leggauss(order)` on the Duffy square-to-triangle
  map, normalize charge weights to sum exactly to each element charge, then
  apply the probe transform around `transform_origin`. Pass original source
  triangles to both periodic and primary kernels. At `h=0`, call the native
  physical zero mode with `trace="principal_value"`; do not change simulator
  `zero_mode_trace_plus`.

- [ ] **Step 4: Run GREEN and panel regressions**

  Run `test_object_panel_interaction.py`, `test_field_kernel.py`, Fortran `test_panel_kernel`, and `test_coulomb_fmm_core_panel` sequentially on a compute node. Expected: all pass; the heavy panel core remains in its opt-in tier.

- [ ] **Step 5: Commit**

  ```bash
  git add beach/fortran_results/panel_quadrature.py \
          beach/fortran_results/object_interaction.py \
          tests/python/test_object_panel_interaction.py tests/python/test_field_kernel.py
  git commit -m "feat: integrate triangle object wrenches"
  ```

### Task 7: Fixed-Source Paths, Work Consistency, And Finite-Shell Oracle

**Files:**
- Create: `beach/fortran_results/periodic_force_oracle.py`
- Modify: `beach/fortran_results/object_interaction.py`
- Modify: `beach/fortran_results/detachment.py`
- Create: `tests/python/test_object_detachment_path.py`
- Create: `tests/python/test_periodic_force_oracle.py`

**Interfaces:**
- `ObjectProbe.vertical_path(displacement_m, *, adaptive=True,
  relative_tolerance=5.0e-3, force_absolute_tolerance_N=1.0e-12,
  work_absolute_tolerance_J=1.0e-18, max_refinement=8,
  torque_origin="geometric_area_centroid", components=True) -> ObjectForcePath`.
- The source snapshot remains exactly at `G_0`; only target quadrature receives translation `(0,0,h)`.
- Every transformed target point must remain at or below configured `box_max[2]` or the resolved outer interface; the first release has no extrapolation beyond that surface.
- `ObjectForcePath` includes line-integrated work, potential-difference work, relative mismatch, refinement diagnostics, and `status in {"converged", "not_converged"}`.
- Produces `finite_shell_wrench(snapshot, probe, transform, image_layers, closure)` where closure is `"symmetric"` or `"e_bottom_zero"`.
- Produces `finite_shell_convergence(snapshot, probe, displacement_m, max_layers=12, relative_tolerance=1.0e-2, force_floor_N=1.0e-12, work_floor_J=1.0e-18)`.
- Large old-run shell rows use Fortran `periodic_far_correction="none"`; exact explicitly replicated direct sums are limited to small fixtures and legacy `M=1` comparison.

- [ ] **Step 1: Write failing moved-probe and work tests**

  Use one net-neutral periodic target and assert original wrench own-image
  `Fz=0`, but a positive z displacement yields the independently summed image
  force while source coordinates remain byte-identical. Use a separate
  non-neutral target to verify the physical closure shift. For a two-charge
  conservative fixture, require:

  ```python
  np.testing.assert_allclose(
      path.electrostatic_work_J,
      path.potential_difference_work_J,
      rtol=5.0e-3,
      atol=1.0e-18,
  )
  ```

  Add a sharply varying near-contact curve that triggers midpoint refinement, a
  force-zero-crossing case that converges through absolute tolerances, and a
  max-refinement case that returns `not_converged`.

- [ ] **Step 2: Write failing shell-closure tests and run RED**

  For a non-neutral sheet fixture, show raw symmetric shells differ from the physical boundary by `Q_cell/(2 epsilon0 A) ez` and corrected shells converge. Run both new test files; expected missing path/oracle APIs.

- [ ] **Step 3: Implement batched path evaluation and closure-visible shells**

  Validate every transformed point for finiteness and box/interface membership
  before entering native code. Batch all points for a candidate displacement
  grid, evaluate each source component without rebuilding plans, integrate
  cumulative trapezoids, and recursively insert interval midpoints where the
  embedded coarse/fine force or work estimate exceeds
  `absolute_tolerance + relative_tolerance * scale`. Preserve both raw and
  corrected shell values:

  ```python
  closure_shift = np.array([0.0, 0.0, total_charge_C / (2.0 * EPS0 * area_xy_m2)])
  e_corrected = e_symmetric + closure_shift
  ```

  Never silently select the final shell when two successive increments fail convergence.

- [ ] **Step 4: Run GREEN and performance guard**

  Run both focused files plus a benchmark assertion that a 65-height path builds the periodic geometry plan once. Expected: work/potential agreement and no per-height plan rebuild.

- [ ] **Step 5: Commit**

  ```bash
  git add beach/fortran_results/periodic_force_oracle.py \
          beach/fortran_results/object_interaction.py beach/fortran_results/detachment.py \
          tests/python/test_object_detachment_path.py tests/python/test_periodic_force_oracle.py
  git commit -m "feat: add fixed-source detachment paths"
  ```

### Task 8: Public CLI, Example, And Bilingual Documentation

**Files:**
- Create: `beach/cli/object_detachment.py`
- Modify: `beach/cli/main.py`
- Modify: `tests/python/test_beachx_cli.py`
- Create: `tests/python/test_object_detachment_cli.py`
- Create: `examples/analyze_periodic_object_detachment.py`
- Modify: `docs/PythonPostprocessAPI.md`
- Modify: `docs/PythonPostprocessAPI.en.md`
- Modify: `docs/PostprocessTutorial.md`
- Modify: `docs/PostprocessTutorial.en.md`
- Modify: `docs/ValidationGuide.md`
- Modify: `docs/ValidationGuide.en.md`
- Modify: `docs/Workflow.md`
- Modify: `docs/Workflow.en.md`
- Modify: `plugins/beach-context/references/python_postprocess_api.md`

**Interfaces:**
- Adds `beachx object-detachment OUTPUT_DIR --config CONFIG --target-mesh-id ID --periodic-model {configured,infinite-physical} --z-max-m Z --z-points N --mass-kg M`.
- Optional flags include `--gravity-m-s2`, `--adhesion-force-n`, `--adhesion-range-m`, `--eta-translation`, `--torque-origin X,Y,Z`, `--cache-dir`, `--library`, and `--output-dir`.
- Writes `instantaneous_wrench.csv`, `path.csv`, `summary.json`, and `report.md` with physics policy, original/override config, cache diagnostics, component sums, convergence status, and assumptions.
- Existing `beachx kernel-forces` stays unchanged except documentation labels its self policy `exclude_target_lattice`.

- [ ] **Step 1: Write failing parser/output tests**

  Require root help to list `object-detachment`, validate positive mass/z-point constraints, and use a fake snapshot evaluator to assert deterministic CSV/JSON columns without requiring a native library. Add one native integration test behind `_kernel_lib()` that exercises a two-object fixture.

- [ ] **Step 2: Run RED**

  Run `pytest -q tests/python/test_beachx_cli.py tests/python/test_object_detachment_cli.py` on `eb`. Expected: missing subcommand/module.

- [ ] **Step 3: Implement CLI as orchestration only**

  Parse inputs, construct `ObjectInteractionSnapshot`, evaluate wrench/path/release, and serialize public records. Do not reproduce force, work, or adhesion formulas in the CLI. JSON must reject NaN/Infinity and encode unavailable potential work as `null` plus a diagnostic reason.

- [ ] **Step 4: Document exact semantics and run docs checks**

  Document `configured` versus `infinite_physical`, primary-only exclusion, point/panel target integration, PV versus pusher trace, frozen-source path, finite-range adhesion, from-rest barrier, non-neutral finite-height restriction, and the old `area_equivalent` distinction. Run focused tests, `ruff check .`, docs sync check, and Markdown link checks on a compute node.

- [ ] **Step 5: Commit**

  ```bash
  git add beach/cli/object_detachment.py beach/cli/main.py \
          tests/python/test_beachx_cli.py tests/python/test_object_detachment_cli.py \
          examples/analyze_periodic_object_detachment.py \
          docs/PythonPostprocessAPI.md docs/PythonPostprocessAPI.en.md \
          docs/PostprocessTutorial.md docs/PostprocessTutorial.en.md \
          docs/ValidationGuide.md docs/ValidationGuide.en.md \
          docs/Workflow.md docs/Workflow.en.md \
          plugins/beach-context/references/python_postprocess_api.md
  git commit -m "docs: add periodic object detachment workflow"
  ```

### Task 9: Reproducible Archived-Run And Simulation-Pair Tooling

**Files:**
- Create: `tools/periodic_object_validation.py`
- Create: `tests/python/test_periodic_object_validation_tool.py`
- Create: `examples/job_scripts/periodic_object_validation_sysa.sh`
- Modify: `docs/ValidationGuide.md`
- Modify: `docs/ValidationGuide.en.md`

**Interfaces:**
- `python tools/periodic_object_validation.py stage --archive-run RUN --validation-root ROOT --binary BIN` creates cache/input/run/submit/analysis directories and structurally verified TOML files.
- `python tools/periodic_object_validation.py verify-inputs --validation-root ROOT`
  permits only output/cache paths, far correction, cache controls,
  `batch_count`, and `history_stride` to differ in their declared stages. It
  normalizes an absent Ewald-layer key to the simulator default before comparing;
  writing the explicit default `4` is therefore metadata normalization, not an
  allowed physics difference.
- `python tools/periodic_object_validation.py verify-run --case-dir DIR --expected-batches N` checks summary, charge/particle ledger identities, finite data, mesh/species/model fingerprints, cache expectations, and executable/input hashes.
- `python tools/periodic_object_validation.py analyze --archive-run RUN --validation-root ROOT --library LIB` invokes public BEACH APIs and writes the comparison artifacts; it does not import scripts from `project_dust_release` as physics evaluators.
- The committed SysA script is a template consumed by `stage`; generated scripts contain resolved absolute paths and metadata commands.

- [ ] **Step 1: Write failing staging/verification tests**

  Use a minimal archived TOML fixture and assert exact patches:

  ```python
  finite = load_toml(root / "input/full/finite_configured.toml")
  infinite = load_toml(root / "input/full/infinite_physical.toml")
  assert finite["sim"]["field_periodic_far_correction"] == "none"
  assert infinite["sim"]["field_periodic_far_correction"] == "cached_kneq0"
  assert infinite["sim"]["field_periodic_ewald_layers"] == 4
  assert infinite["sim"]["field_periodic_generation_tolerance"] == 1.0e-8
  assert finite["sim"]["rng_seed"] == infinite["sim"]["rng_seed"] == 12345
  ```

  Mutate `dt`, mesh radius, MPI count metadata, Ewald layers away from the
  normalized default, and seed in separate tests and require fail-closed
  verification. Test cache-prime/100-batch/full expected diagnostics with
  synthetic summary files.

- [ ] **Step 2: Run RED**

  Run `pytest -q tests/python/test_periodic_object_validation_tool.py` on `eb`. Expected: missing tool.

- [ ] **Step 3: Implement deterministic staging and artifact schemas**

  Create the fixed validation layout from the design. Use structured TOML parsing/writing, SHA-256 every input/binary, and write a `manifest.json` containing commit, dirty flag, module list, Slurm resources, source run/version, and allowed config differences. Generated jobs use KUDPC `#SBATCH --rsc`, plain `srun`, `/usr/bin/time -v`, and no generic `--ntasks`, `--cpus-per-task`, or nested MPI launcher.

- [ ] **Step 4: Run GREEN and dry-stage the real archive**

  Run focused tests, then stage the real archive under a temporary pytest path and call `verify-inputs`; do not submit. Expected: only declared finite/infinite and smoke/full differences.

- [ ] **Step 5: Commit**

  ```bash
  git add tools/periodic_object_validation.py \
          tests/python/test_periodic_object_validation_tool.py \
          examples/job_scripts/periodic_object_validation_sysa.sh \
          docs/ValidationGuide.md docs/ValidationGuide.en.md
  git commit -m "feat: add periodic validation harness"
  ```

### Task 10: Archived R20260625-0002 Post-Processing Validation

**Files:**
- Generate only under: `/LARGE1/gr20001/b36291/codex-tmp/beach-periodic-object-force-validation/analysis/`

**Interfaces:**
- Uses archived version `1.3.0-v1.3.0`, config, mesh, final charges, and charge-history rows at batches `149001`, `180001`, `279001`, and `280000` when available.
- Evaluates target meshes 6 and 7 for configured finite shell, infinite physical field, own-image/other/k0 components, `0..2R` paths, work/barrier/speed, and shell convergence.
- Imports old generated CSVs as labeled compatibility evidence only.

- [ ] **Step 1: Build the final native shared library on a compute node**

  Require a clean worktree and record the exact commit first. Inspect `hostname`,
  `module list`, `spartition`, and `qgroup`, then choose a visible SysB compute
  partition. For example, when `eb` is visible, run:

  ```bash
  tssrun -p eb -t 0:45:0 --rsc p=1:t=8:c=8 \
    bash -lc 'cd /LARGE1/gr20001/b36291/codex-tmp/BEACH-periodic-object-force && \
    make build-kernel && sha256sum build/libbeach_field_kernel.so'
  ```

  Record library hash and branch commit in the validation manifest.

- [ ] **Step 2: Run archived analysis through Slurm**

  Submit one SysB batch job with `p=1:t=28:c=28`, four-hour limit, repository on `PYTHONPATH`, cache/output under the validation root, and:

  ```bash
  srun python tools/periodic_object_validation.py analyze \
    --archive-run /LARGE1/gr20001/b36291/project_dust_release/runs/local_charging_release_baseline/R20260625-0002 \
    --validation-root /LARGE1/gr20001/b36291/codex-tmp/beach-periodic-object-force-validation \
    --library /LARGE1/gr20001/b36291/codex-tmp/BEACH-periodic-object-force/build/libbeach_field_kernel.so
  ```

- [ ] **Step 3: Verify numerical gates**

  Require component sums within 1%, work/potential mismatch within 0.5%, two-step shell convergence or explicit `not_converged`, no path point above the box/interface, and no unlabeled incomplete cached field. Confirm `area_equivalent` appears only as potential context.

- [ ] **Step 4: Write the Japanese review from machine-readable results**

  `review_ja.md` must state separately whether positive instantaneous `Fz`, adhesion/gravity margin, endpoint work, and barrier-free-from-rest release survive each boundary/self policy. It must identify the share from target replicas and physical `k=0` and label frozen-charge/contact limitations.

- [ ] **Step 5: Preserve artifact manifest**

  Hash every CSV/JSON/figure/report, record the analysis command/job ID/runtime, and ensure the git worktree remains clean because generated data is outside the repository.

### Task 11: Paired SysA Smoke And Full Simulations

**Files:**
- Generate only under: `/LARGE1/gr20001/b36291/codex-tmp/beach-periodic-object-force-validation/`

**Interfaces:**
- Runs `cache_prime` for 1 batch, paired smoke cases for 100 batches/history stride 10, then paired full cases for 280,000 batches/history stride 1000.
- `finite_configured` changes only output path from archived input.
- `infinite_physical` additionally changes far correction to `cached_kneq0`,
  writes the simulator-default Ewald layers `4` explicitly, and sets cache
  directory and generation tolerance `1.0e-8`.
- Uses one clean feature-branch executable for both cases, `rng_seed=12345`, `p=6:t=112:c=112`, Intel 2023.2/Intel MPI 2023.2, and direct `srun`.

- [ ] **Step 1: Recheck SysA and build one MPI executable**

  On the login control plane inspect `hostname` and `module list`, detect the
  currently loaded `Sys*` module, switch that module to `SysA` only when needed,
  then run `module list`, `spartition`, and `qgroup`; require visible and allocated
  `gr20001a` immediately before submission. Require a clean worktree and record
  the exact commit before building/installing on SysA compute resources:

  ```bash
  tssrun -p gr20001a -t 1:00:0 --rsc p=1:t=28:c=28 \
    bash -lc 'cd /LARGE1/gr20001/b36291/codex-tmp/BEACH-periodic-object-force && \
    PREFIX=/LARGE1/gr20001/b36291/codex-tmp/beach-periodic-object-force-validation/install \
    BEACH_VERSION_MODE=git FPM_ACTION=install FPM_PROFILE=release \
    FPM_FC=mpiifort FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" \
    ./build.sh --compiler mpiifort'
  ```

  Copy/hash the installed binary into `bin/<commit>/beach`; both cases use that
  immutable path. The hash is invalid if source dirtiness appears between the
  clean check and build completion.

- [ ] **Step 2: Stage and verify all inputs**

  Run `stage` and `verify-inputs` on a compute node. Require an empty output root, clean git commit, matching binary/input hashes, and no undeclared physical differences.

- [ ] **Step 3: Submit and gate cache-prime/smoke job**

  Re-run `module list`, `spartition`, and `qgroup`, then submit generated
  `smoke_sysa.sh` to visible `gr20001a` for one hour with
  `p=6:t=112:c=112`. The generated script exports
  `OMP_NUM_THREADS=112`, `OMP_PROC_BIND=spread`, and `OMP_PLACES=cores`. It runs
  cache prime, finite 100, and infinite 100 sequentially. Require prime
  `cache_hit=false/build_count=1`, later infinite
  `cache_hit=true/build_count=0`, identical fingerprint, 100 completed batches,
  3392 elements, seven meshes, finite values, ledger identities, and no fallback.

- [ ] **Step 4: Submit the full pair after smoke passes**

  Submit `full_finite_sysa.sh` with 12 hours and `full_infinite_sysa.sh` with 14
  hours, each `p=6:t=112:c=112` and the same three OpenMP exports. They may run
  concurrently because outputs are separate and the cached operator is already
  atomically published. Monitor with `squeue` at reasonable intervals; if a job
  reaches wall time, resume only with matching model/cache fingerprints.

- [ ] **Step 5: Verify and compare full outputs**

  Require exit code zero, `batches=280000`, 280 saved history rows, matching mesh/species fingerprints, and only the declared model fingerprint difference. Run `verify-run` then `analyze` on the pair. Produce `run_summary.csv`, `charge_history_pair.csv`, `particle_ledger_pair.csv`, `mesh_potential_pair.csv`, `object_wrench.csv`, `object_path_curves.csv`, `object_path_summary.csv`, `finite_shell_convergence.csv`, `comparison_matrix.csv`, figures, and updated `review_ja.md`.

  Interpret `archived_v1.3 -> new finite` as version/reproducibility drift and `new finite -> new infinite` as the boundary-model effect. A single paired seed is not a probability estimate; effects below block/history variability are labeled stochastic-sensitive.

### Task 12: Full Verification, Review, Merge, And Push

**Files:**
- Verify: all changed source, tests, docs, and validation manifests

- [ ] **Step 1: Run static and tiered suites without concurrent fpm**

  On SysB compute nodes run, sequentially:

  ```bash
  make test-l0
  make test-l1
  make test-l2
  make test-l3
  make test-fortran-far-correction
  ruff check .
  ```

  Use one `sbatch` job with enough wall time for L3/far correction. Record each target runtime and peak RSS. Do not add new cached tests to default L1.

- [ ] **Step 2: Run MPI cache and simulator smoke on SysA**

  Build MPI tests once, then launch MPI payloads directly with KUDPC `tssrun`/`srun`, not nested `mpiexec`. Require cold/warm cache fingerprint agreement across ranks and clean two-rank completion.

- [ ] **Step 3: Review complete branch and validation conclusions**

  Run a task-scoped review after each implementation task and a final whole-branch physics/code review against the design. Resolve every Critical/Important issue, confirm public API/docs parity, verify generated artifacts from exact hashes, and ensure no claim exceeds frozen-charge/contact/one-seed evidence.

- [ ] **Step 4: Commit final evidence references and push feature branch**

  Commit only documentation links/summary values that are stable and useful in the repository; keep large artifacts in `codex-tmp`. Push `feat/periodic-object-detachment` and verify its remote hash.

- [ ] **Step 5: Merge main and push after all gates pass**

  Fast-forward `main`, push `origin/main`, and verify local `main`, remote-tracking `origin/main`, and `git ls-remote origin refs/heads/main` all equal the reviewed feature hash. Do not merge while a SysA full job or required test is still running.
