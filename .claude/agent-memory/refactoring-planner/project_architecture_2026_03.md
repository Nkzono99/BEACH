---
name: BEACH architecture analysis (2026-03)
description: Comprehensive codebase structure analysis identifying coupling issues, duplicated code, god-types, and phased refactoring plan
type: project
---

Full codebase audit findings as of 2026-03-25 (code at commit 09ab222):

**Scale:** Fortran ~13,700 LOC (48 files), Python ~10,900 LOC (41 files), 20 Fortran + 8 Python test suites.

**Fortran structural issues (prioritized):**
1. `resolve_periodic2_collision_config` duplicated identically in `bem_collision.f90` (L520-561), `bem_mesh.f90` (L288-329), and partially in `bem_field_solver_config.f90` (L38-84).
2. `cross()` function duplicated in `bem_pusher.f90`, `bem_collision.f90`, `bem_mesh.f90`.
3. `field_solver_type` is a god-type (78 fields, ~190 LOC type def alone) mixing treecode, FMM, periodic state. Submodules: tree (514 LOC), eval (447 LOC), config (180 LOC).
4. `sim_config` has 51 fields including 5 `has_*` booleans; `particle_species_spec` has 10 `has_*` booleans -- parse concern leaking into runtime types.
5. `default_app_config()` overrides `sim_config()` defaults with different values (e.g., dt=1e-9 vs 1e-12) -- dual-default drift risk.
6. `app/main.f90` L49-102: 8-branch conditional for optional parameter passing to `run_absorption_insulator`.
7. `bem_app_config_runtime.f90` (900 LOC) mixes mesh building, particle seeding, batch init, reservoir face, MPI split, sheath injection.
8. 349 `error stop` across 27 files with no structured error reporting.
9. 24 `.gcov` files in repo root (should be gitignored).

**Python structural issues (prioritized):**
1. `_coerce_periodic2` returns a 7-element tuple (axes, lengths, origins, image_layers, far_correction, ewald_alpha, ewald_layers) -- should be a dataclass.
2. `_auto_periodic2_from_result` called from 3 files (animation, mesh, potential) with repeated TOML probing.
3. `potential.py` (953 LOC) mixes computation, periodic coercion, precomputed shortcircuit, and TOML probing.
4. `__init__.py.__all__` exports private names: `_select_frame_columns`, `_surface_charge_density`.
5. `facade.py` (701 LOC) duplicates full parameter signatures from delegated functions.
6. 6 legacy `cli_*.py` shim files in `beach/` root alongside `beach/cli/legacy.py`.
7. Fortran-Python output contract (CSV headers, summary.txt keys) is implicit -- column changes silently break Python readers.

**Phased plan (in priority order):**
- Phase 1: Extract shared `resolve_periodic2_collision_config` + `cross()` to utility modules
- Phase 2: Python periodic2 dataclass + TOML probe caching
- Phase 3: Split `field_solver_type` into treecode_state/fmm_state/periodic2_config sub-types
- Phase 4: Separate `has_*` booleans from `sim_config`/`particle_species_spec` into parse_state
- Phase 5: Eliminate main.f90 8-branch via sentinel values or run_options struct
- Phase 6: Split `bem_app_config_runtime.f90` into mesh_builder + batch_builder
- Phase 7: Remove private names from `__all__`, delete legacy CLI shims
- Phase 8: Explicit Fortran-Python output contract constants + smoke tests

**Why:** These findings drive the refactoring roadmap. Phase 1-2 are non-breaking and unblock later phases.
**How to apply:** Use as the reference checklist when planning or reviewing refactoring PRs. Phase order is important -- later phases depend on earlier ones.
