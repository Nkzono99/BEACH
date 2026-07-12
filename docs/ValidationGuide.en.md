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
A path/work/shell convergence check on a fixed saved mesh and source
discretization covers only a subset of the numerical axes required for Level 2.
Passing that subset does not establish Level 2 as a whole when mesh or source
refinement has not been performed.

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
   reconstruction self term `area_equivalent`. For an object crossing a
   periodic seam, keep the all-source/cache input in its saved representation
   and unwrap only the target probe into a mesh-connected branch. Require
   `target_geometry_representation="periodic2_mesh_connected"` and confirm that
   primary subtraction, target integration, torque origin, and geometry radius
   use that same branch. Production mass, adhesion, and the `0..2R` path use the
   explicitly declared model radius; compare the connected-geometry bounding
   radius against it within a declared tolerance. Record both the origin policy
   and origin coordinates for every 3D torque.
3. For triangle target integration, compare Gauss-Duffy order 3 and order 7.
   This checks target-side area integration; it is not source-mesh refinement.
   Check source discretization independently with source meshes at different
   resolutions.
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
8. Qualify the cached infinite-periodic implementation with analytic oracles.
   The lightweight pytest oracle is opt-in through
   `BEACH_RUN_FIELD_KERNEL_CACHE_TESTS=1`. In the production strict chain,
   however, the analysis job must run `probe-periodic-oracles` before
   `analyze --require-complete` and provide a write-once receipt bound to the
   same staged library. The primary oracle uses the production `point` source
   and `softened_point` kernel; a `triangle_p0` oracle is evaluated with the
   same library as an auxiliary check.
9. For a uniform non-neutral plane (`sigma=Q/A`) under the `E_bottom=0`
   closure, require `E_z=0` below, the surface PV
   `E_z=sigma/(2 epsilon0)`, `E_z=sigma/epsilon0` above, surface pressure
   `sigma^2/(2 epsilon0)`, and total object force
   `Q^2/(2 epsilon0 A)`. Fix the potential gauge at sheet height `z0` so that
   `phi=0` on and below the sheet and
   `phi=-sigma (z-z0)/epsilon0` above it. The surface PV applies to the field
   trace; potential is continuous in this gauge. This force is
   closure-dependent Maxwell traction, not a universal free-space self force.
10. The production `point` oracle spatially averages a 24x24 x/y sample at
   each of `z=-0.25, 0, +0.25` and evaluates 4x4 and 8x8 source grids.
   Off-surface modulation RMS must decrease from 4 to 8 and be at most
   `0.12 V/m` on the fine grid. Point-centroid wrenches on both grids must be
   within 12% of the analytic result, and the inter-grid force difference must
   be within 1%. Force, transverse-force, and torque errors must each decrease
   from 4 to 8. The primary-free subtraction residual must also decrease and
   be within 1% on the fine grid. The decomposition requires other force,
   external force, and `total_external-target_periodic_images` to each be at
   most `1e-12`, and confirms that the reported component residual equals the
   maximum of those three terms and the primary-free residual. A
   separate single-primary check verifies
   removal of primary self force and subtraction of softened self potential
   energy `-K q^2/a` to relative `1e-11`. This is the
   `exclude_primary_keep_images` contract: exclude only the primary while
   retaining periodic images. The auxiliary `triangle_p0` oracle requires the
   Gauss-Duffy order-3/order-7 wrench difference to be within 1%.
11. For a neutral `sigma_0 cos(kx)` sheet, the same 4x4/8x8 snapshots and
   operator caches evaluate analytic `Ex`, `Ez`, and `phi` at
   `|z-z0|=0.25 m, 0.50 m`. Field and potential errors must each decrease from
   the 4x4 to the 8x8 grid and be within 8% on the fine grid. The field-vector
   L2 amplitude ratio and potential L2 amplitude ratio from 0.25 m to 0.50 m
   are compared with `exp(-k*0.25 m)`. Both-grid ratios and relative errors are
   stored in the receipt, and each fine-grid ratio error must be within 18%.
   Both heights are evaluated together in each snapshot, so this adds no
   operator-cache build. The production point source requires a
   charge-neutrality ratio `<=1e-12` on both grids and paired `+z/-z` samples at
   both heights. Tangential field is even, normal field is odd, and
   potential is even; separate field and potential parity errors must each be
   within 8%. A softened-point micro-oracle with `a/L=1e-6` compares analytic
   softened field/potential with the direct evaluator and ordinary with direct
   at four points `r/a=0,1,2,3`, all to relative `1e-11`; normalized self field
   must be at most `32 epsilon_machine`. This is a local kernel-contract check,
   not a substitute for the periodic closure. The 12% uniform-plane, 8%
   cosine analytic, and 18% decay-ratio thresholds are smoke gates for the
   production ABI/cache path. They
   neither replace the 0.5% object-path and 1% finite-shell convergence
   criteria nor establish 8%/12%/18% physical accuracy. Do not reuse plane-oracle
   errors as force or torque error bars for the saved sphere mesh or its source
   discretization. Sphere and source refinement remain `not_evaluated` until
   they are performed separately.

## SysA Comparison of Archived, Finite, and Infinite Runs

`tools/periodic_object_validation.py` is a fail-closed harness for comparing an
existing archive with `finite_configured` and `infinite_physical` runs that
retain the same physical inputs. Put the validation root outside the repository
and require it to be empty before staging.
The four `stage` path arguments (archive, validation root, binary, and library)
accept only `[A-Za-z0-9_./:+-]+`; whitespace, `@`, `$`, quotes, and newlines are
rejected. The validation root must be outside both the repository and archive,
and staging rejects any existing symlink ancestor in that write destination.
Archive, binary, and library are read-only inputs and are resolved to their
canonical targets after the safe-character check. After staging, verify and
strict paths bind validation-root descendants and archive-analysis metadata to
their exact canonical strings and reject existing symlinks from root to leaf.

```bash
current_sys="$(module -t list 2>&1 | grep -E '^Sys(A|B|C|CL|G)/' | head -n 1 || true)"
if [ -n "${current_sys}" ] && [ "${current_sys}" != "SysA/2022" ]; then
  module switch "${current_sys}" SysA/2022
elif [ -z "${current_sys}" ]; then
  module load SysA/2022
fi
module load intel/2023.2 intelmpi/2023.2

python3.11 tools/periodic_object_validation.py stage \
  --archive-run /path/to/RUN \
  --validation-root /LARGE1/.../beach-periodic-object-force-validation \
  --binary /path/to/clean-build/beach \
  --library /path/to/clean-build/libbeach_field_kernel.so \
  --require-clean-source
python3.11 tools/periodic_object_validation.py verify-inputs \
  --validation-root /LARGE1/.../beach-periodic-object-force-validation
bash /LARGE1/.../beach-periodic-object-force-validation/submit/submit_chain.sh
```

`stage` snapshots the executable, Python source, analysis tool, and native
kernel library with hashes. Operationally, build the executable and library
from the same clean commit under SysA/Intel before passing them to `stage`.
`stage --require-clean-source` reads executable `--build-info` and the library
C ABI build info. It accepts them only when version, mode, full source SHA, and
`SHA:clean` agree with each other and with the staged source commit.
`verify-inputs` rereads the staged artifacts, while each simulation
`summary.txt` and the plane-oracle receipt are bound to the same build origin.
`analyze --require-complete` requires this build origin unconditionally, so a
production manifest staged without `--require-clean-source` cannot pass strict
qualification.
Production staging and strict analysis require both
`input/release_kernel_base.toml` and
`analysis/local_release/release_model_summary.json`. Their canonical paths,
hashes, schemas, finite numeric values, and physical ranges are rechecked, so
work and speed calculations cannot silently fall back to defaults. Non-strict
archive-only analysis retains the prior default-compatible behavior.
Compiler identity is not embedded in the artifacts; retain the SysA/Intel
module and build logs as separate execution evidence.

Generated jobs unset `PYTHONHOME`, set `PYTHONNOUSERSITE=1`, and set
`PYTHONPATH="${SOURCE_ROOT}"` exactly. They do not inherit the submitting
environment's user site or `PYTHONPATH`.

The generated DAG has six jobs: smoke (including cache prime), finite 140000,
finite 280000, infinite 140000, infinite 280000, and analysis. Each model
starts fresh through batch 140000 and restarts from the verified checkpoint
into a separate 280000-segment directory. The final analysis job depends on
both 280000 jobs with `afterok`, runs the mandatory production plane oracle,
and then invokes `analyze --require-complete`. `submit_chain.sh` atomically
updates a persistent `job_ids.json` journal after each submission and uses an
exclusive lock. It refuses to submit the same DAG again even after a partial
submission. Strict analysis rechecks six unique job IDs, SysA/Intel
module/hash logs for every job, source commit, resources, and zero exit status
for the five predecessors, and the ID of the running analysis job.

`verify-run` checks the schema, geometry, charge and particle ledgers,
`mesh_sources.csv`, `mesh_potential.csv`, every MPI-rank checkpoint, history,
restart fingerprints, cache fingerprint, and cache-file hash. It writes an
execution receipt once under `provenance/verified/`, including every regular
file in the output tree. Later verification never replaces that receipt; it
compares current files against it and also fixes the restart-parent and
cache-prime receipt dependencies. Cached runs and post-processing evaluators
are bound to the verified cache-prime receipt hash as well as the cache
fingerprint, path, and file hash; they cannot silently reuse a different
cache. Strict final analysis rechecks the staged
archive input, manifest, source snapshot, case graph, binary, and library, and
requires an existing execution receipt for all seven cases. Archive-only
preflight may omit `--require-complete`, but its `missing` and `not_evaluated`
rows cannot support a conclusion.

The comparison semantics are fixed. `archive -> new finite` includes version,
restart, and runtime drift and is not a boundary-only comparison. Only
`new finite -> new infinite`, produced by the same new executable and inputs,
isolates the boundary-model effect. The `vdw_work` speed bracket uses an
equivalent finite-range constant force that preserves initial adhesion force
and total adhesion work; it does not reproduce the original `1/s^2` barrier
shape. The primary constrained-translation result uses `eta_translation=1`;
the archived value `0.5` is retained separately as a sensitivity value.

The harness separately audits the archived `Fz>0` estimators. It binds
`force_timeseries.csv`, `moving_sphere_force_curves.csv`,
`moving_sphere_release_summary.csv`, and `moving_sphere_model_summary.json` by
hash and exact schema, then compares only the common batch/mesh keys
`(149001,7)`, `(180001,6)`, `(279001,6)`, and `(279001,7)`. The archived
center-charge direct estimate, local-pair proxy, and moving-sphere force/work at
`z=0` and `z=2R` are recorded beside current native-finite other-object,
target-periodic-image, and total components in
`legacy_estimator_comparison.csv`. Batch 280000 is not interpolated because it
is absent from the legacy tables. Numerical closeness is not a gate: the old
estimator excludes the entire target from its self sources, while the current
evaluator excludes only the primary and retains target periodic images. The
strict gate covers input coverage, finite values, charge/radius/batch mapping,
and internal endpoint consistency between the legacy curve and summary.

The analysis compares per-object charge and element-charge L1/L2/Linf
differences at every saved batch common to the runs. Expensive force/path
evaluation is restricted to mesh 7 at batch 149001, mesh 6 at batch 180001,
both at 279001, and final batch 280000. `snapshot_manifest.csv` distinguishes
history from final charge and hashes both the charge vector and source file.
`comparison_matrix.csv` assigns separate `comparison_kind` values to archived
version drift, field-closure changes on the same charge, charging-history
response under a common infinite evaluator, and the end-to-end difference.
The strict comparison artifact contract allows exactly those four kinds and
requires force, endpoint work, minimum available energy, from-rest barrier,
and endpoint reachability for each. Every referenced snapshot must resolve;
the frozen-field comparison must use one charge snapshot; effective
far-correction pairs must match their declarations; and no structural
`missing`, `invalid`, or `not_evaluated` row is allowed.

Mesh IDs and ordering must match exactly across archived, new-finite, and
new-infinite cases. Triangle coordinates use the declared tolerance
`max(1e-18 m, 64 epsilon Lbox)`. In addition, the new-finite and new-infinite
`field_source_model`, `field_kernel_id`, and the `mesh_sources.csv`
source/template kind, surface model, `epsilon_r`, and element count must match
the staged inputs and each other exactly. Missing equivalent metadata in the
old archive does not weaken the new-run contract. For cached evaluators,
`object_wrench.csv` also records cache hit, build count, fingerprint, path,
file hash, and cache-prime receipt hash.

`analyze --require-complete` performs strict input, receipt, oracle, and
geometry checks and generates artifacts in a temporary directory. A failure
confined to that analysis temporary directory removes it only when all
write-once state, including the qualified oracle receipt, is complete and
valid, and `analysis/` remains unpublished (absent or empty). The same
validation root may then be retried only inside the original allocation while
preserving the original analysis job ID, or through a site-permitted same-ID
requeue. An ordinary new `sbatch`, any source/tool/library/input change, or the
absence of a same-ID retry path requires a new validation root. A partial
oracle-generation failure, oracle configuration or cache files without a
receipt, and any retry after atomic publication of `analysis/` also require a
new validation root.

`numerical_qualification_for_local_frozen_model` is a subset gate evaluated
with the saved sphere mesh and source discretization held fixed. It covers the
exact 30 path/wrench keys, six shell groups, and path/work convergence. The
finite-shell relative tolerance is 1%; the adaptive path tolerance is 0.5%,
with a `1e-12 N` absolute force floor and a `1e-18 J` absolute work floor. It
also requires the force floor to be no more than 0.5% of peak force, the
instantaneous detachment-force margin to exceed `1e-12 N`, and the energy margin
from the barrier decision boundary to exceed 10% of the energy tolerance.

For this gate, `status="qualified"` means only that path integration,
work/potential consistency, decision resolution, and finite-shell convergence
passed on the fixed saved discretization. It does not evaluate saved-sphere mesh
refinement, source-discretization refinement, or an absolute sphere force/torque
error bound; those remain `not_evaluated`. The 8%/12%/18% plane-oracle thresholds do
not fill in a sphere error bar. A barrier or speed is not claimed when these
resolution gates fail, even if integration itself converged. If a path or shell
is unconverged, barrier and speed are `not_claimed_unqualified`. In a non-neutral
system with a remaining constant upper field, `0..2R` work and speed are local
frozen-field diagnostics, not escape energy or speed at infinity.
In `object_path_summary.csv`, `path_start_m` and `path_end_m` record the
validated final grid after adaptive refinement, rather than the initial grid,
and prove that it retained the `0..2R` interval.

Strict analysis creates exactly 14 artifacts in a temporary directory:
`run_summary.csv`, `charge_history_pair.csv`, `particle_ledger_pair.csv`,
`mesh_potential_pair.csv`, `snapshot_manifest.csv`, `object_wrench.csv`,
`object_path_curves.csv`, `object_path_summary.csv`,
`finite_shell_convergence.csv`, `comparison_matrix.csv`,
`legacy_estimator_comparison.csv`, `analysis_summary.json`, `review_ja.md`, and
`artifacts.json`. It checks the
exact file set, nonempty files, and size/hash manifest, fsyncs them, and then
atomically publishes `analysis/`. Strict success then requires both
`numerical_qualification_for_local_frozen_model.status="qualified"` and
`comparison_artifact_contract.status="complete"`. If either gate fails, the
CLI exits nonzero but never overwrites the published diagnostics. Rerunning
strict analysis after publication requires a new validation root.

An exactly complete and readable 14-artifact set is Level-1 execution
evidence. A zero strict-CLI exit code also establishes the structural, oracle,
comparison, and fixed-saved-discretization path/work/shell subset gates above.
It does not establish Level 2 as a whole. While saved-sphere mesh and source
refinement remain `not_evaluated`, Level 2 as a whole remains incomplete. Those
refinements, model-selection checks, and sensitivity checks are required before
considering a Level-3 physical claim.

See [Physics release verification](PhysicsReleaseVerification.en.html) for small reference contracts.
