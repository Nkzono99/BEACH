# Periodic object force and detachment-path design

Date: 2026-07-11

Status: written specification for user review

## 1. Purpose

BEACH output must support two related frozen-charge diagnostics for one selected
object in an x/y-periodic, z-open electrostatic system:

1. Evaluate the instantaneous three-dimensional force and torque at the saved
   position, including the selected object's own nonzero periodic replicas.
2. Hold the original periodic source lattice fixed, move only the selected
   central object, integrate the force along the path, and convert the net work
   into a conditional translational speed.

The initial validation target is
`/LARGE1/gr20001/b36291/project_dust_release/runs/local_charging_release_baseline/R20260625-0002`.
Generated validation data belongs under
`/LARGE1/gr20001/b36291/codex-tmp/` and is not committed to BEACH.

This feature is a frozen-charge post-process. It does not claim that the saved
charge distribution, plasma, polarization, contact geometry, or neighboring
objects remain unchanged during a real release event.

## 2. Existing-result audit

The old run currently has several different quantities that must not be mixed:

- The mesh-center model replaces each object by its net charge at an
  area-weighted center. It has no periodic images and gives no force crossing in
  this run; its largest direct `Fz` is approximately `-5.91e-12 N`.
- `area_equivalent` is a panel self-potential approximation. It is not an
  area-averaged object force and is not evidence for detachment.
- The element-resolved beachx-compatible force uses a finite 3x3 x/y image
  shell of non-target objects. It removes the selected object before expanding
  images, so the selected object's own periodic replicas are also missing. Its
  saved maxima are approximately `4.55e-6 N` for mesh 6 and `6.96e-6 N` for
  mesh 7.
- The existing moving-sphere model moves the selected panel charges through the
  central-cell field of non-target panels. It contains no periodic images and
  subtracts a constant adhesion-plus-gravity force over the entire `0..2R`
  path.
- `FieldKernel` can aggregate 3D force and torque, but the object helper zeros
  the target source charges. Under a periodic kernel this removes both the
  primary target and every target replica.
- The Python `cached_kneq0` kernel returns the nonzero periodic modes only. The
  simulator separately adds the physical `k=0` provider with `E_bottom=0`.
  Therefore the current Python cached force is not the simulator's total field
  for a non-neutral cell.

The new result set must reproduce these definitions as compatibility rows, then
identify which changed assumption causes any sign or magnitude change.

## 3. Physical contract

### 3.1 Frozen source snapshot and probe

Let `Q_all` and `G_0` be the saved charges and original geometry. Let `T` be the
selected central object. The source snapshot is immutable. A probe contains the
selected object's original charge distribution and geometry and may be rigidly
transformed without modifying the source snapshot.

The external field used for the selected object is

```text
E_env(r) = E_periodic_physical[Q_all, G_0](r)
         - E_free[Q_T, G_T,0](r).
```

This is the single-grain defect construction. It removes only the selected
object in the primary cell. It retains:

- all non-target objects in the primary cell;
- all periodic replicas of non-target objects; and
- every nonzero periodic replica of the selected object.

It excludes the selected physical body's internal electrostatic self-force and
self-torque. Moving the selected object in the source snapshot, or zeroing its
source charges before building the periodic field, is incorrect for this
contract.

### 3.2 Infinite periodic field

For `field_periodic_far_correction="cached_kneq0"`, the physical field is

```text
E_periodic_physical = E_cached,k!=0 + ez * E_physical,k=0 + E_uniform.
```

The physical zero mode uses the same source-height projection and
`E_bottom=0` boundary condition as the simulator. Point sources are horizontal
sheets in the zero-mode projection. Triangle P0 sources use the exact triangle
height projection already implemented by the simulator.

The force API must not silently call a `k!=0` result a total periodic field.
It must either compose the physical zero mode or label the result as incomplete.
For this feature, incomplete cached results are rejected.

The first release supports simulator snapshots with no active outer-plasma
closure. A configured outer model is rejected with an actionable error until a
matching frozen outer-field snapshot is available. A configured uniform
`sim.e0` is retained and reported as a separate component.

For a non-neutral cell, the physical zero mode approaches
`Q_cell/(epsilon0 A)` above all sources. Work to infinite height therefore does
not generally converge. The path must end at or below the configured box or
outer-interface height unless the caller supplies a physically closed outer
field. Reported speeds are local frozen-field scales at the chosen endpoint,
not escape speeds to infinity.

### 3.3 Force and torque

For a point-source result, the target approximation is

```text
F = sum_i q_i E_env(r_i)
tau_o = sum_i (r_i - o) cross q_i E_env(r_i).
```

For a triangle P0 result, the target is integrated over each uniformly charged
panel:

```text
F = sum_i integral_Ai sigma_i E_env(r) dA
tau_o = sum_i integral_Ai (r - o) cross sigma_i E_env(r) dA.
```

The source model is never silently downgraded. Triangle target integration uses
explicit quadrature and reports its order. At the original surface, full-field
and primary-subtraction evaluations use the same principal-value convention.
Near-contact conclusions require quadrature and mesh-refinement convergence;
one centroid sample per panel is retained only as a labeled compatibility
approximation.

The torque origin is explicit. The API accepts a user vector and named
geometric references. A geometric area centroid is not called a center of mass.
The validation reports torque about the fitted sphere center and, when a contact
point is supplied, about that pivot.

### 3.4 Physical source decomposition

Every wrench can be returned as the sum of these physical groups:

- `other_objects_all_images`;
- `target_periodic_images`;
- `external_uniform`;
- `total_external`.

For diagnosis, the periodic groups also expose `periodic_kneq0` and
`physical_k0`. The negative free-primary term is available as implementation
metadata but is not presented as an additional physical object.

This decomposition determines whether a positive `Fz` is due to nearby objects,
the selected object's replicas, the non-neutral zero mode, or a prescribed
uniform field.

### 3.5 Single-object path

For the requested vertical path,

```text
r_i(h) = r_i(0) + h ez,
```

while every source in `G_0`, including the target replicas, remains at its
original position. Each path sample returns the full force and torque, not only
`Fz`.

The electrostatic work is

```text
W_elec(h) = integral_0^h Fz(s) ds.
```

When a consistent potential is available, the implementation also computes

```text
W_phi(h) = U_env(0) - U_env(h),
U_env(h) = sum_i q_i phi_env(r_i(h)),
```

and treats agreement between `W_elec` and `W_phi`, and between `Fz` and
`-dU/dh`, as validation checks. The path sampler refines intervals until the
requested force/work tolerance is met or reports non-convergence. It does not
rebuild the full periodic source plan at every height.

### 3.6 Release and speed semantics

Electrostatics and mechanics are separate steps. For a body initially at rest,

```text
K(h) = W_elec(h) - m g h - W_adh(h) - W_diss(h),
v(h) = sqrt(2 eta_translation max(K(h), 0) / m).
```

The result contains at least:

- electrostatic work and electrostatic-only speed;
- gravity-corrected work and speed;
- adhesion-corrected net work and speed;
- the minimum cumulative net work and a `barrier_free_from_rest` flag;
- the first height at which a from-rest path becomes energetically inaccessible;
- the endpoint and maximum reachable speed under the supplied model.

Positive endpoint work alone is not called detachment from rest. A from-rest
path requires nonnegative available kinetic energy over the whole sampled path,
within numerical tolerance.

Adhesion is represented by an explicit finite-range force/work profile. The
compatibility report may include the old constant `5.46875e-7 N` resistance
over the full `0..2R` interval, but labels it `legacy_constant_path_proxy`.
It is not the default physical adhesion law. The translational partition
defaults to `eta_translation=1` for a constrained pure-translation path; the
old value `0.5` is a labeled sensitivity row. Rotational motion requires a
generalized path and is outside the first vertical-path implementation.

The instantaneous detachment table separately reports

```text
Fz_electric,
Fz_electric - m g,
Fz_electric - m g - F_adhesion(0),
```

so raw force, gravity, and contact assumptions remain distinguishable.

## 4. Public Python design

The existing APIs remain source compatible. Their current target-zeroing
behavior is documented as `exclude_target_lattice` and is not silently changed.
New behavior uses an immutable snapshot and a separate probe:

```python
snapshot = beach.ObjectInteractionSnapshot.from_result(
    result,
    step=-1,
    config_path="input/beach.toml",
    periodic_model="infinite_physical",
    cache_dir=(
        "/LARGE1/gr20001/b36291/codex-tmp/"
        "beach-periodic-object-force-cache"
    ),
)

probe = snapshot.object_probe(
    mesh_id=6,
    self_policy="exclude_primary_keep_images",
    target_integration="auto",
    quadrature_order=7,
)

wrench = probe.wrench(
    transform=None,
    torque_origin=fitted_sphere_center,
    components=True,
)

path = probe.vertical_path(
    displacement_m=np.linspace(0.0, 2.0 * radius_m, 65),
    adaptive=True,
    components=True,
)

release = path.evaluate_release(
    mass_kg=mass_kg,
    gravity_m_s2=1.62,
    adhesion=adhesion_profile,
    eta_translation=1.0,
)
```

The concrete public records are immutable dataclasses:

- `ObjectWrench`: identifiers, transform, torque origin, total charge, force,
  torque, component wrenches, and numerical metadata.
- `ObjectForcePath`: transforms/displacements, force and torque arrays,
  component arrays, cumulative electrostatic work, potential-difference work,
  and convergence diagnostics.
- `DetachmentResult`: mass and resistance assumptions, cumulative net energy,
  speeds, barrier status, instantaneous margin, and endpoint metadata.
- `AdhesionProfile`: `none`, finite-range constant traction, or tabulated
  traction. All models expose both force and cumulative work.

Arrays are SI-valued NumPy arrays with documented shapes and are read-only in
returned records. Invalid mesh IDs, unsupported boundary closures, missing
native kernels, non-finite inputs, and paths beyond the supported interface
raise explicit exceptions. Adaptive or shell non-convergence returns a result
whose diagnostic status is `not_converged`; it never returns an unlabeled last
iterate. There is no fallback to a different physics model.

`periodic_model="configured"` reproduces the run configuration, including its
finite shell when far correction is `none`. `periodic_model="infinite_physical"`
is an explicit analysis override to `cached_kneq0 + physical_k0`; output
metadata records both the original configuration and the override. The old run
comparison evaluates both modes.

## 5. Internal architecture

### 5.1 Kernel composition

The implementation has four independent evaluators:

1. A full periodic nonzero-mode `FieldKernel` built once from all original
   sources.
2. A native physical zero-mode evaluator using the simulator's existing
   height-plan implementation and `E_bottom=0` state.
3. A free-space primary-target evaluator built once from the original selected
   object. Point sources use the exact blocked direct kernel with the configured
   softening. Triangle sources use the exact direct panel kernel with the same
   principal-value trace convention as the full field.
4. The existing uniform-field term.

An object probe batches all quadrature points for one or more transforms,
evaluates each component, subtracts the free primary field, and integrates
target charge weights into force, torque, and potential energy.

The native zero-mode binding is a separate ABI handle rather than adding
arguments to the existing field-kernel build symbols. It accepts point heights
or triangle vertex heights, periodic area, source charges, `E_bottom`, and a
potential gauge, and evaluates field/potential with an explicit trace policy.
This keeps the existing `FieldKernel` ABI stable.

`FieldKernelOptions` is extended to carry periodic cache directory and
generation tolerance from the config. The C ABI gains versioned build entry
points for those options; the existing entry points remain available for ABI
compatibility.

### 5.2 Target integration

- Point source: target charge points are the saved element centroids.
- Triangle P0: quadrature points and weights are generated once on the original
  target triangles, then rigidly transformed. Charge weights sum to each
  element's saved total charge.
- The free primary evaluator uses the same point softening or triangle panel
  kernel and principal-value convention as the full source evaluator.
- Compatibility centroid integration is selectable, labeled, and never chosen
  by `auto` for triangle P0.

### 5.3 Performance

- Source plans and cached periodic operators are built once per source geometry.
- Charge updates do not rebuild geometry.
- Target points from multiple path heights are evaluated in bounded batches.
- Component decomposition reuses linearity and cached source states.
- Direct finite-shell validation is blocked/vectorized and limited to requested
  targets and saved batches.
- Cache files and large validation products use
  `/LARGE1/gr20001/b36291/codex-tmp/`.

The validation report records plan-build time, cache hit/miss, per-height
evaluation time, peak target batch size, FMM order, image layers, quadrature
order, and tolerances. Performance changes never alter the physical policy.

## 6. Finite-shell oracle

The independent oracle explicitly constructs the `(2M+1)^2` x/y source cells,
removes only the central-cell target, and evaluates the transformed target
probe. It reports separate contributions from other objects and target replicas.

For a non-neutral cell, a symmetric finite shell converges to the symmetric
zero-mode closure, not automatically to `E_bottom=0`. Comparison with the
configured physical infinite-periodic field therefore adds

```text
+ Q_cell / (2 epsilon0 A) ez
```

to the symmetric-shell field. The report contains both the raw symmetric-shell
value and the corrected `E_bottom=0` value so the closure is visible.

Shells are increased over integer `M` from zero through `M_max`. The validation
default is `M_max=12`. Convergence requires two successive shell increments to
change both force and path work by at most 1% relative, with absolute floors of
`1e-12 N` and `1e-18 J`. Lack of convergence is a result, not a reason to
select the last shell silently.

## 7. Test and validation matrix

### 7.1 Automated tests

- One isolated point charge under free boundaries has zero external force and
  torque after primary exclusion.
- One charge in a symmetric periodic lattice has zero whole-object `Fz` at
  `h=0`; moving only the probe produces the expected image contribution.
- Two-object free-space forces match direct Coulomb values and action-reaction
  symmetry.
- `exclude_target_lattice` and `exclude_primary_keep_images` differ exactly by
  the target-replica contribution.
- A non-neutral periodic fixture verifies
  `cached k!=0 + physical k0` against the simulator snapshot oracle.
- Raw symmetric shells plus the analytic closure shift converge to the
  `E_bottom=0` result.
- Physical group components and numerical kernel components sum to the total
  within tolerance.
- Force and torque obey translation/rotation covariance and the torque-origin
  shift identity.
- Triangle quadrature converges with order, and a mesh-refinement fixture checks
  force/torque convergence near a surface.
- Path work converges with sampling refinement and agrees with potential-energy
  difference.
- Speed conversion, negative-work clamping, finite-range adhesion work, and the
  barrier-free-from-rest rule have unit tests.
- Existing `calc_object_forces_kernel` results remain unchanged unless a new
  policy is explicitly selected.

Small analytic fixtures use test-specific tolerances derived from their exact
scale. The old-run study uses these declared defaults unless a stricter value
is recorded: 0.5% relative work-integration tolerance, 1% component-sum and
shell tolerance, and 1% change between target quadrature orders 3 and 7. A
failed tolerance remains visible in the report and blocks an unqualified
physical conclusion.

Slow cached/operator and large-shell tests are assigned to the existing
L3/far-correction opt-in tiers. Small analytic tests remain in L1/L2. No heavy
FMM target is added to the default `make test` path.

### 7.2 R20260625-0002 validation

The saved source model in this run is point-centroid legacy. For meshes 6 and
7, evaluate the final state, each mesh's first beachx-compatible adhesion-only
crossing batch (`180001` for mesh 6 and `149001` for mesh 7), and the latest
available charge-history row. Produce:

- instantaneous 3D force/torque and detachment margins;
- vertical force/torque curves from `0` to `2R`;
- cumulative electrostatic, gravity, adhesion, and net work;
- electrostatic-only, gravity-corrected, and adhesion-bracket speed curves;
- barrier-free and endpoint-positive classifications;
- component decomposition for other objects, target replicas, `k!=0`, physical
  `k=0`, and uniform field;
- finite-shell convergence for integer `M` from zero through `M_max`;
- point-centroid compatibility and triangle/target-quadrature sensitivity where
  the saved source model permits it;
- numerical timings and cache diagnostics.

The comparison table includes the old mesh-center, beachx-compatible 3x3,
moving-sphere, configured finite-shell, corrected finite-shell convergence, and
infinite-periodic defect definitions. It states explicitly that
`area_equivalent` is a potential self-term, not a force model.

Validation artifacts include machine-readable CSV/JSON, force/work/speed plots,
and a Japanese Markdown review. The review labels each conclusion as one of:

- numerically verified under the selected model;
- sensitive to periodic closure or discretization;
- conditional on frozen charge;
- conditional on adhesion/contact assumptions; or
- not supported by the available output.

### 7.3 Paired SysA simulations

Post-processing the archived charge state is the primary boundary-model
comparison, but it does not show how the boundary model changes charge
accumulation. Run two new simulations from the same initial state with the
feature-branch executable:

- `finite_configured`: the archived input unchanged except for output/cache
  paths; it keeps `field_periodic_far_correction="none"` and one image layer.
- `infinite_physical`: identical input except for
  `field_periodic_far_correction="cached_kneq0"`, an explicit cache directory,
  and the physical `E_bottom=0` zero mode selected by the simulator.

Both cases retain `rng_seed=12345`, six MPI ranks, 112 OpenMP threads per rank,
the point-centroid source model, mesh, particle injection, time step, batch
duration, and `batch_count=280000`. They use the same newly built executable so
that `finite_configured - archived_v1.3` measures version/reproducibility effects
and `infinite_physical - finite_configured` isolates the boundary-model effect.
Bitwise identity is not required across compiler/runtime versions; charge,
particle-ledger, field, force, work, and convergence differences are reported.

The durable validation root is
`/LARGE1/gr20001/b36291/codex-tmp/beach-periodic-object-force-validation/`.
Before the full pair, run matched 100-batch smoke cases with `history_stride=10`.
The smoke gate requires successful cache generation/reuse, finite diagnostics,
matching mesh/particle counts, no model fallback, and valid restart/output
metadata. Only after both smoke cases pass are the full cases submitted to the
visible SysA group queue. Each batch script records the executable hash, git
commit, Modules, input hash, Slurm resources, cache fingerprint, start/end time,
and exit status.

The paired comparison includes saved batches common to both cases and the final
state. It compares total and per-object charge, element-charge spatial norms,
absorbed/escaped ledgers, `rel_change`, `k!=0`/physical-`k=0` field components,
instantaneous wrench, path work, barrier status, and speed brackets for meshes 6
and 7. Conclusions distinguish archived-version drift, stochastic sampling,
boundary-model response, and post-processing-only changes.

## 8. Success criteria

The feature is complete when:

1. Both requested diagnostics use the same documented primary-only exclusion
   policy and include target periodic replicas.
2. A cached infinite-periodic result includes the simulator's physical zero
   mode and exposes the component decomposition.
3. The direct-shell oracle and analytic closure reproduce the infinite result
   within stated convergence tolerances on controlled fixtures.
4. Point and triangle source models never silently substitute for one another.
5. Work, potential difference, and speed conversion pass their consistency
   checks, including the from-rest barrier rule.
6. The old run has a reproducible comparison report under `codex-tmp` that says
   whether its positive `Fz` and positive work survive the corrected model.
7. Public documentation and examples explain assumptions, limitations, and SI
   units, and all relevant L0-L3 tests pass in their intended tiers.

## 9. Explicit non-goals

- Recharging, conductor-like redistribution, dielectric polarization, or plasma
  response while the object moves.
- Contact-mechanics inference from the BEACH surface mesh alone.
- A claim of natural release probability from saved batches or multiplier
  sweeps.
- Escape speed to infinity for a non-neutral `E_bottom=0` cell without an outer
  closure.
- Coupled translation and rotation dynamics in the first implementation.
- Changing the simulator's production boundary condition or the legacy object
  force API by default.
