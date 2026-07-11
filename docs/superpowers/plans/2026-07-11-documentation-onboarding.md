# Documentation Onboarding Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Provide one verified beginner path from installation through physical and numerical validity checks.

**Architecture:** Keep `docs/` as the canonical bilingual source and let `tools/sync_starlight_docs.py` generate Starlight content. Preserve existing slugs while adding four focused onboarding pages and making `beachx config init` identical to one committed runnable example.

**Tech Stack:** Markdown, Astro Starlight, Python 3.10+, pytest, TOML/JSON Schema, Fortran/fpm.

## Global Constraints

- Preserve every existing public documentation slug.
- Keep Japanese and English page structure synchronized.
- Keep `docs/` canonical; do not hand-edit generated Starlight content.
- Do not add dependencies.
- Run Fortran smoke tests only on a KUDPC compute node.

---

### Task 1: Official Beginner Configuration

**Files:**
- Create: `examples/tutorial_insulator.toml`
- Modify: `beach/config/core.py`
- Modify: `tests/python/test_config_cli.py`

**Interfaces:**
- Produces: `default_config() -> dict[str, Any]` semantically equal to `examples/tutorial_insulator.toml`.

- [ ] **Step 1: Write the failing semantic-equivalence test**

```python
def test_default_config_matches_official_tutorial_case() -> None:
    assert normalize_config_document(default_config()) == load_config_file(
        Path("examples/tutorial_insulator.toml")
    )
```

- [ ] **Step 2: Run the test on a compute node and verify it fails because the example is missing or differs**

Run: `python3.11 -m pytest -q tests/python/test_config_cli.py::test_default_config_matches_official_tutorial_case`

- [ ] **Step 3: Add the direct/free-space one-batch insulator case and make `default_config()` return the same data**

The case contains one `volume_seed` electron species, one plane template, `field_solver="direct"`, `field_bc_mode="free"`, and file output under `outputs/latest`.

- [ ] **Step 4: Run the test and lint the example**

Run: `python3.11 -m pytest -q tests/python/test_config_cli.py`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add beach/config/core.py examples/tutorial_insulator.toml tests/python/test_config_cli.py
git commit -m "docs: establish one official beginner case"
```

### Task 2: Bilingual Onboarding Pages

**Files:**
- Create: `docs/Installation.md`, `docs/Installation.en.md`
- Create: `docs/Tutorial.md`, `docs/Tutorial.en.md`
- Create: `docs/ValidationGuide.md`, `docs/ValidationGuide.en.md`
- Create: `docs/Troubleshooting.md`, `docs/Troubleshooting.en.md`
- Modify: `docs/index.md`, `docs/index.en.md`
- Modify: `docs/ConfigurationRecipes.md`, `docs/ConfigurationRecipes.en.md`

**Interfaces:**
- Consumes: `examples/tutorial_insulator.toml` from Task 1.
- Produces: four stable page stems consumed by the sync script.

- [ ] **Step 1: Add a failing sync-contract test for the four bilingual page pairs and Japanese descriptions**

Add assertions in `tests/python/test_docs_sync.py` that every page stem appears for `root` and `en`, and that Japanese pages have Japanese descriptions.

- [ ] **Step 2: Run the test and verify missing page registrations fail**

Run: `python3.11 -m pytest -q tests/python/test_docs_sync.py`

- [ ] **Step 3: Write the four page pairs**

Installation covers prerequisites/PATH/update/removal. Tutorial embeds the complete official case and expected output checks. Validation separates execution completion from physical/numerical checks. Troubleshooting maps common symptoms to concrete commands and unsupported combinations.

- [ ] **Step 4: Rewrite the indexes and mark recipe snippets as diffs from the official case**

The index contains the supported-scope table, command/data flow, first-use sequence, and four goal links; remove duplicate documentation lists.

- [ ] **Step 5: Run Markdown/static tests**

Run: `python3.11 -m pytest -q tests/python/test_docs_sync.py`
Expected: pass.

- [ ] **Step 6: Commit**

```bash
git add docs tests/python/test_docs_sync.py
git commit -m "docs: add installation tutorial and validity path"
```

### Task 3: Workflow-Ordered Starlight Navigation

**Files:**
- Modify: `tools/sync_starlight_docs.py`
- Modify: `docs-site/astro.config.mjs`
- Modify: `tests/python/test_docs_sync.py`

**Interfaces:**
- Consumes: bilingual page stems from Task 2.
- Produces: generated slugs `installation`, `tutorial`, `validation-guide`, and `troubleshooting`.

- [ ] **Step 1: Extend the failing test with exact sidebar group and item order**

Expected groups: Start, Usage, Reference, Numerics, Developers, AI Agents.

- [ ] **Step 2: Register pages and link mappings in `sync_starlight_docs.py`**

Use Japanese descriptions for root locale and English descriptions for `en`.

- [ ] **Step 3: Reorder the Starlight sidebar without changing old slugs**

Move `agent-user-guide` into its own AI group and `workflow` into Developers.

- [ ] **Step 4: Generate and check content**

Run: `python3.11 tools/sync_starlight_docs.py && python3.11 tools/sync_starlight_docs.py --check`
Expected: no out-of-date pages.

- [ ] **Step 5: Commit**

```bash
git add tools/sync_starlight_docs.py docs-site/astro.config.mjs tests/python/test_docs_sync.py
git commit -m "docs: order site navigation by user workflow"
```

### Task 4: End-to-End Verification

**Files:**
- Verify only; fix scoped defects in files from Tasks 1-3.

**Interfaces:**
- Consumes: official case and generated docs.
- Produces: compute-node smoke evidence and a production Starlight build.

- [ ] **Step 1: Run Python and style checks on a compute node**

Run: `python3.11 -m pytest -q && python3.11 -m ruff check .`
Expected: all pass.

- [ ] **Step 2: Run the official beginner Fortran case on a compute node**

Run: built `beach examples/tutorial_insulator.toml` with a temporary output directory.
Expected: exit 0, `summary.txt` exists, `batches=1`, charge ledger residual is finite.

- [ ] **Step 3: Build the Starlight production site**

Run: `make docs-starlight`
Expected: Astro build succeeds without broken internal links.

- [ ] **Step 4: Run the physics release gate after all docs commits**

Run inside Slurm: `make test-physics-release`
Expected: manifest schema 2 status passed, every gate within memory budget, convergence categories complete.

- [ ] **Step 5: Confirm the worktree is clean**

Run: `git status --short`
Expected: no uncommitted onboarding files remain.
