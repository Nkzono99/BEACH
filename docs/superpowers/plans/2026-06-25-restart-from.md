# restart_from Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `[output].restart_from` so resume checkpoint input can differ from the new output directory.

**Architecture:** Add a string field to `app_config`, parse it from `[output]`, validate it only applies with `resume = true`, and pass a derived checkpoint directory into the existing restart loader. Keep every writer on `output.dir`. Mirror the config change in JSON schema, Python workload estimation, and docs.

**Tech Stack:** Fortran fpm runtime, toml-f parser, Python CLI utilities, JSON Schema Draft 7, pytest/fpm tests.

## Global Constraints

- Always review in Japanese.
- On KUDPC login nodes, do not run Fortran tests directly; use Slurm when runtime tests are needed.
- Do not change public APIs without updating docs and examples.
- Add or extend tests when modifying resume logic.
- Keep Python side lightweight; avoid new dependencies.
- Ignore `*.i90` files.

---

### Task 1: Config And Runtime Restart Directory

**Files:**
- Modify: `src/config/bem_app_config_types.f90`
- Modify: `src/config/app_config_parser/bem_app_config_parser.f90`
- Modify: `app/main.f90`
- Test: `tests/fortran/test_app_config_parser.f90`

**Interfaces:**
- Consumes: `[output].restart_from` string from TOML.
- Produces: `app_config%output_restart_from` and runtime checkpoint source selection.

- [ ] **Step 1: Write the failing Fortran parser test**

Add `restart_from = "outputs/parent"` to the existing `[output]` fixture and assert:

```fortran
call assert_true(trim(cfg%output_restart_from) == 'outputs/parent', 'output.restart_from mismatch')
```

- [ ] **Step 2: Run parser test on a compute node and verify it fails**

Run: `tssrun --rsc p=1:t=00:10:00 -q gr20001g FPM_ACTION=test ./build.sh --target test_app_config_parser`

Expected: FAIL because `output_restart_from` is not defined or is empty.

- [ ] **Step 3: Implement minimal Fortran support**

Add `character(len=256) :: output_restart_from = ''`, parse `restart_from`, reject it unless `resume_output` is true, and use it for `load_restart_checkpoint` input only.

- [ ] **Step 4: Run parser test again**

Run: `tssrun --rsc p=1:t=00:10:00 -q gr20001g FPM_ACTION=test ./build.sh --target test_app_config_parser`

Expected: PASS.

### Task 2: Schema, Lint, And Workload Estimator

**Files:**
- Modify: `schemas/beach.schema.json`
- Modify: `beach/config/schemas/beach.schema.json`
- Modify: `beach/cli/estimate_fortran_workload.py`
- Test: `tests/python/test_config_cli.py`
- Test: `tests/python/test_estimate_fortran_workload.py`

**Interfaces:**
- Consumes: `output.restart_from` in Python-loaded config dictionaries.
- Produces: schema validation and workload resume batch counts that use `restart_from` when present.

- [ ] **Step 1: Write failing Python tests**

Add tests proving `completed_batches_from_resume_config` reads `restart_from/summary.txt` and that lint accepts `restart_from` with `resume = true` but rejects it without resume.

- [ ] **Step 2: Run focused pytest and verify failure**

Run: `pytest -q tests/python/test_estimate_fortran_workload.py::test_completed_batches_from_resume_config_reads_restart_from tests/python/test_config_cli.py::test_lint_accepts_output_restart_from tests/python/test_config_cli.py::test_lint_rejects_restart_from_without_resume`

Expected: FAIL before implementation.

- [ ] **Step 3: Implement Python/schema support**

Add the schema property and conditional validation. Change workload estimator to choose `restart_from` over `dir` for checkpoint summaries.

- [ ] **Step 4: Run focused pytest again**

Run: same focused `pytest -q ...`

Expected: PASS.

### Task 3: Documentation And Verification

**Files:**
- Modify: `SPEC.md`
- Modify: `docs/fortran_workflow.md`
- Modify: `docs/fortran_parameter_file.md`
- Modify: `docs/agent-user-guide.md`

**Interfaces:**
- Consumes: implemented `[output].restart_from` semantics.
- Produces: user-facing docs matching runtime/schema behavior.

- [ ] **Step 1: Update docs**

Document that `restart_from` changes checkpoint input only, while all new files are written to `output.dir`.

- [ ] **Step 2: Run safe verification**

Run local Python tests with `pytest -q` if acceptable. Run Fortran build/test through `tssrun` or `sbatch` on `gr20001g`.

Expected: focused tests pass; broader verification passes or any blocker is reported.
