# Release Test Design Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce duplicated HPC release testing and make expensive Fortran targets individually observable without weakening correctness coverage.

**Architecture:** Keep cumulative development tiers intact, add a dedicated release-correctness target, and route target execution through a timing-aware sequential shell runner. Split assertion-free diagnostics and release-profile benchmarks from far-correction correctness.

**Tech Stack:** GNU Make, Bash, fpm, Fortran 2008, pytest.

## Global Constraints

- Do not run Fortran tests directly on KUDPC login nodes.
- Do not execute multiple fpm tests concurrently because they share `build/`.
- Keep existing cumulative L0-L3 commands backward compatible.
- Do not introduce machine-dependent hard runtime thresholds.

---

### Task 1: Lock the release-gate contract

**Files:**
- Modify: `tests/python/test_physics_release_gate_contract.py`
- Modify: `Makefile`
- Modify: `tools/run_physics_release_gate.sh`

**Interfaces:**
- Produces: `test-fortran-release-correctness`, diagnostic and benchmark targets, and manifest timing CSV keys.

- [x] Add failing source-contract tests for the dedicated release target, diagnostic separation, benchmark profile, and timing artifacts.
- [x] Run the focused pytest target through a KUDPC compute-node allocation and confirm the expected failures.
- [x] Add the Makefile target lists and update the release gate commands.
- [x] Re-run the focused tests and confirm they pass.

### Task 2: Add sequential per-target timing

**Files:**
- Create: `tools/run_fortran_targets.sh`
- Create: `tests/python/test_run_fortran_targets.py`
- Modify: `Makefile`

**Interfaces:**
- Consumes: target names as positional arguments and existing build variables through the environment.
- Produces: optional CSV columns `target,profile,status,elapsed_seconds`.

- [x] Add failing tests using a temporary fake build command for success, failure, CSV creation, and sequential order.
- [x] Run the focused pytest target through a compute node and confirm the expected failures.
- [x] Implement the runner and replace the inline Make recipe.
- [x] Confirm focused and release-contract Python tests pass.

### Task 3: Separate correctness, diagnostics, and benchmark code

**Files:**
- Modify: `tests/fortran/test_periodic2_infinite_operator.f90`
- Create: `benchmarks/fortran/benchmark_periodic2_runtime.f90`
- Modify: `fpm.toml`
- Modify: `Makefile`

**Interfaces:**
- Produces: debug correctness test without timing assertions and an opt-in release-profile runtime benchmark.

- [x] Add source-contract tests that fail while the runtime block remains in the correctness test.
- [x] Move the runtime loop into a non-installable fpm example benchmark.
- [x] Verify the Python contract and Fortran build on a compute node.

### Task 4: Remove redundant cache generation

**Files:**
- Modify: `tests/fortran/test_periodic2_operator_cache.f90`

**Interfaces:**
- Produces: a unique cache directory and exactly one intentional cold build before the warm load.

- [x] Add a source-contract test rejecting the probe-plan build pattern.
- [x] Replace the probe plan with a unique temporary cache directory and existing cleanup helpers.
- [x] Run the focused cache target on a compute node and confirm cold/warm/corruption behavior.

### Task 5: Update documentation and verify tiers

**Files:**
- Modify: `docs/PhysicsReleaseVerification.md`
- Modify: `docs/PhysicsReleaseVerification.en.md`
- Modify: `docs/Workflow.md`
- Modify: `docs/Workflow.en.md`
- Modify: `docs/agent-user-guide.md`
- Modify: `docs/agent-user-guide.en.md`

**Interfaces:**
- Documents: release correctness subset, diagnostic opt-in, benchmark opt-in, and per-target timing CSV.

- [x] Update Japanese and English documentation together.
- [x] Run Python tests and L0/L1 checks through a compute node.
- [x] Run the new release-correctness and far-correction correctness gates sequentially on a compute node.
- [x] Inspect timing CSVs and compare the new release-gate composition with the previous 940-second baseline.
