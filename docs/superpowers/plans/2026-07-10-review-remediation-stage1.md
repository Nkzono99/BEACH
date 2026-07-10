# Review Remediation Stage 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the highest-risk silent failures and data corruption behaviors before the larger particle, restart, periodic-kernel, and potential-closure redesigns.

**Architecture:** Introduce a shell-free filesystem boundary, make periodic collision incompleteness explicit, force unsafe mixed-sign tree nodes to descend, and make Python output readers reject corrupted dense histories. Keep existing public entry points compatible where possible, but change silent failure into an explicit error. Align inspect and documentation with the actual cost and far-correction behavior.

**Tech Stack:** Fortran 2008 with `iso_c_binding`, fpm, OpenMP, Python 3, NumPy, pytest, Ruff, JSON/Markdown documentation, KUDPC Slurm.

## Global Constraints

- Always review and report in Japanese.
- The final architecture remains the design in `docs/superpowers/specs/2026-07-10-review-remediation-design.md`; Stage 1 guards are temporary where the design says so.
- Do not change public APIs without updating docs and examples in the same stage.
- Use TDD: add a focused failing test, verify the expected failure, implement the minimum correction, then verify focused and broader tests.
- Do not run `make test*`, `fpm test`, Fortran executables, or substantial Python test payloads directly on KUDPC login nodes.
- Run only one fpm test controller at a time because targets share `build/`.
- Route runtime tests through `tssrun` or `sbatch` + `srun` with `p=1` for controller commands.
- Keep long FMM diagnostics out of the L1 default path.
- Do not add Python dependencies.
- Ignore `*.i90` files.
- Do not modify or squash the pre-existing commit `2e47428`.
- Each task ends in a focused commit after its task review is clean.

---

### Task 1: Shell-Free Output Directory Creation

**Files:**
- Create: `src/runtime/bem_filesystem.f90`
- Modify: `src/runtime/bem_output_writer.f90`
- Test: `tests/fortran/test_output_writer_io.f90`

**Interfaces:**
- Preserve `bem_output_writer::ensure_output_dir(out_dir)`.
- Add an internal `bem_filesystem::create_directories(path, status)` implementation.
- Bind directly to POSIX `mkdir`, `opendir`, and `closedir`; never construct a shell command from `out_dir`.

- [ ] **Step 1: Add a failing literal-path regression test**

  Add a test path containing spaces, `$()`, semicolon, double quote, and single quote. Verify that `ensure_output_dir` creates that literal directory and does not create the marker path named inside the shell syntax.

- [ ] **Step 2: Run the focused target on a compute node and verify RED**

  Run:

  ```bash
  tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
    bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && FPM_ACTION=test ./build.sh --target test_output_writer_io'
  ```

  Expected: the current quoted shell command either fails to create the literal path or performs an unintended expansion; the marker assertion demonstrates the defect without deleting unrelated paths.

- [ ] **Step 3: Implement recursive POSIX directory creation**

  Convert a Fortran string to a null-terminated `c_char` array. Walk path prefixes separated by `/`, call `mkdir(mode=0777)`, and treat an existing prefix as success only when `opendir` confirms it is a directory. Support absolute paths, repeated separators, `.`, `..`, spaces, quotes, and shell metacharacters as literal path content.

- [ ] **Step 4: Remove `execute_command_line` from production output directory creation**

  Keep the public wrapper in `bem_output_writer`, delegate to `bem_filesystem`, and return a clear error for an empty path, an existing non-directory component, or an OS failure.

- [ ] **Step 5: Verify GREEN and static safety**

  Re-run `test_output_writer_io`, then run `rg -n "execute_command_line" src/runtime` and confirm no user output path reaches a shell command.

- [ ] **Step 6: Commit**

  Commit message: `fix: create output directories without a shell`

---

### Task 2: Fail Closed On Incomplete Periodic Collision Queries

**Files:**
- Modify: `src/physics/bem_collision.f90`
- Modify: `src/physics/bem_collision_types.f90` if the hit type is defined separately
- Modify: `src/runtime/simulator/bem_simulator_loop.f90`
- Modify: `src/runtime/simulator/bem_simulator.f90`
- Test: `tests/fortran/test_dynamics_basic.f90`
- Test: `tests/fortran/test_simulator.f90`

**Interfaces:**
- Add public collision status constants: `collision_query_ok`, `collision_query_image_limit`, `collision_query_index_range`.
- Add an optional status output to `find_first_hit`, or a status field to `hit_info`, while preserving existing positional arguments.
- Existing callers that do not request status must receive an explicit error on incomplete work rather than `has_hit=.false.`.

- [ ] **Step 1: Replace the existing fail-open assertion with a failing status assertion**

  The 5000-cell test must require `collision_query_image_limit`, not no-hit success. Add a second test with coordinates outside safe `i32` conversion range and require `collision_query_index_range`.

- [ ] **Step 2: Run `test_dynamics_basic` on a compute node and verify RED**

  Expected: no status exists or the query still reports ordinary no-hit.

- [ ] **Step 3: Implement range-safe image enumeration**

  Compute cell bounds in `i64` or real range-checked intermediates. Check conversion bounds before converting. Return an explicit status when the configured safety limit would be exceeded.

- [ ] **Step 4: Propagate incomplete status through the simulator**

  Do not invoke `error stop` inside the OpenMP particle loop. Store a thread-safe failure flag/status, finish the parallel region, and terminate with a message containing batch, particle, and status context.

- [ ] **Step 5: Verify focused tests**

  Run `test_dynamics_basic` and `test_simulator` sequentially on a compute node.

- [ ] **Step 6: Commit**

  Commit message: `fix: fail closed on incomplete collision queries`

---

### Task 3: Disable Unsafe Mixed-Sign Tree Monopoles

**Files:**
- Modify: `src/physics/field_solver/bem_field_solver_eval.f90`
- Modify: `src/physics/field_solver/bem_field_solver_tree.f90` only if additional node state is required
- Test: `tests/fortran/test_dynamics_field_solver.f90`
- Document: `docs/FieldSolvers.md`
- Document: `docs/FieldSolvers.en.md`

**Interfaces:**
- Preserve field solver configuration and method signatures.
- Mixed-sign internal nodes must descend until a leaf/direct interaction in Stage 1.
- Same-sign nodes retain existing monopole behavior and performance.

- [ ] **Step 1: Add a failing cancellation sweep**

  Construct positive/negative source pairs with net ratios around the old `1e-10` threshold. Choose a query where the exact dipole field is nonzero and assert treecode matches direct within the existing tree accuracy contract.

- [ ] **Step 2: Add a same-sign characterization test**

  Confirm existing same-sign accuracy remains unchanged; this test should pass before implementation.

- [ ] **Step 3: Run `test_dynamics_field_solver` and verify the mixed-sign test fails**

- [ ] **Step 4: Implement correctness-first node acceptance**

  Accept a monopole only when `abs(node_q)` equals `node_abs_q` within a rounding tolerance derived from machine epsilon. Any physically mixed-sign node descends. Do not use the unstable signed charge centroid for an accepted mixed-sign approximation.

- [ ] **Step 5: Update solver documentation**

  Document the Stage 1 rule and note that a later monopole+dipole criterion will recover mixed-sign performance without reintroducing the error.

- [ ] **Step 6: Verify and commit**

  Run `test_dynamics_field_solver`; commit as `fix: descend mixed-sign treecode nodes`.

---

### Task 4: Reject Corrupted Dense History Snapshots

**Files:**
- Modify: `beach/fortran_results/history.py`
- Test: `tests/python/test_fortran_results.py`
- Document: `docs/PythonPostprocessAPI.md`
- Document: `docs/PythonPostprocessAPI.en.md`

**Interfaces:**
- Preserve `FortranChargeHistory` public construction and snapshot access.
- Replace implicit zero-fill and last-row-wins behavior with `ValueError` for malformed snapshots.

- [ ] **Step 1: Add failing tests for malformed snapshots**

  Cover missing element, duplicate element, out-of-range index, non-finite charge, inconsistent `processed_particles`, inconsistent `rel_change`, and decreasing batch order. Assert error messages identify the batch and defect.

- [ ] **Step 2: Run focused pytest on a compute node and verify RED**

  Run only the new test node IDs with one pytest controller.

- [ ] **Step 3: Implement strict batch indexing**

  Validate rows before allocating the dense snapshot. Require exactly one row for every `elem_idx=1..mesh_nelem`, finite values, and constant batch metadata. Keep file indexing lazy where currently lazy.

- [ ] **Step 4: Update Python API documentation**

  State that Fortran history is dense and corruption raises `ValueError`; missing elements are not physical zero.

- [ ] **Step 5: Verify focused and full Python tests**

  Run the focused tests first, then `pytest -q` on a compute node.

- [ ] **Step 6: Commit**

  Commit message: `fix: reject incomplete history snapshots`

---

### Task 5: Make Inspect Potential Recalculation Explicit

**Files:**
- Modify: `beach/cli/inspect_fortran_output.py`
- Test: `tests/python/test_cli_error_messages.py` or a focused inspect CLI test module
- Document: `docs/PythonPostprocessAPI.md`
- Document: `docs/PythonPostprocessAPI.en.md`

**Interfaces:**
- Add `--recompute-potential` as an explicit command flag.
- Without the flag, print min/max only from `mesh_potential.csv` when present.
- Do not call `Beach.compute_potential` during ordinary inspection.

- [ ] **Step 1: Add a failing no-recompute test**

  Patch or instrument `Beach.compute_potential` to fail if called. Run ordinary inspect on a result with triangles but no precomputed potential and assert inspection succeeds without potential min/max.

- [ ] **Step 2: Add precomputed and explicit-recompute tests**

  Assert `mesh_potential.csv` min/max is reported without recomputation, and `--recompute-potential` invokes the existing calculation.

- [ ] **Step 3: Run focused pytest and verify RED**

- [ ] **Step 4: Implement the explicit cost boundary**

  Reuse `result.mesh_potential_v` when available. Preserve plotting behavior: a requested potential plot is an explicit expensive action and may compute potential.

- [ ] **Step 5: Update help and docs, verify, and commit**

  Run focused CLI tests and full Python tests; commit as `perf: avoid implicit quadratic inspect work`.

---

### Task 6: Align Far-Correction Contract And Clear Lint Debt

**Files:**
- Modify: `SPEC.md`
- Modify: `tools/sync_starlight_docs.py`
- Verify: `docs/Parameters.md`
- Verify: `docs/FMMCore.md`
- Verify: `beach/config/schemas/beach.schema.json`
- Test: relevant schema/config tests if text changes reveal a mismatch

**Interfaces:**
- During the staged migration, `field_periodic_far_correction="auto"` remains an alias of `none`.
- `m2l_root_oracle` remains explicit opt-in and is not described as production exact periodic physics.
- The future analytic M2L default is deferred to a versioned Stage 5 migration.

- [ ] **Step 1: Add a static contract check if one is absent**

  Ensure tests or a lightweight script assert that defaults and `auto` semantics agree across Fortran config, FMM core, schema, and SPEC-facing generated metadata.

- [ ] **Step 2: Correct `SPEC.md`**

  Remove the statement that `auto` defaults to and normalizes as `m2l_root_oracle`. Describe `none` as a finite-image approximation and the oracle as an expensive diagnostic backend pending the periodic-kernel redesign.

- [ ] **Step 3: Remove the unused `re` import**

  Keep this isolated lint cleanup in the same documentation/tooling commit because it is the only existing Ruff failure reported by the review baseline.

- [ ] **Step 4: Verify static and lint gates**

  Run `git diff --check`, JSON schema parse checks, and `ruff check .`. Run any focused config tests on a compute node when behavior assertions change.

- [ ] **Step 5: Commit**

  Commit message: `docs: align periodic far-correction contract`

---

## Stage 1 Completion Gate

- [ ] All six task commits exist and have clean task reviews.
- [ ] `git diff --check` passes.
- [ ] `ruff check .` passes.
- [ ] Full Python tests pass on a compute node.
- [ ] `make check` passes in the supported build environment.
- [ ] L1 Fortran targets pass sequentially on a compute node.
- [ ] The branch remains based on the user's existing `2e47428` commit and contains no unrelated changes.
- [ ] A Stage 1 whole-branch review has no unresolved Critical or Important findings.
