# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- Added a config-driven field-kernel runtime benchmark that separates mesh construction, solver initialization, charge refresh, and volume/near-panel field and potential evaluation for `point` and `triangle_p0` comparisons.

### Fixed
- The Camphor `ifx` install profile no longer passes unsupported inline/report options or enables the unstable IPO link path, allowing release executables to link reliably with Intel 2023.2.
- Reflected and periodic particle events now place surviving particles a scale-aware distance inside the box, preventing zero-time boundary chatter caused by a subnormal one-ULP offset at zero-valued faces. The existing eight-event safety limit is unchanged.

### Changed
- Documented `kinetic_1d` as the standard outer-sheath model and reclassified `unified_linear_response` as advanced rough-surface linear screening. Model identifiers and runtime defaults are unchanged.
- **BREAKING**: Split `outer_plasma.photoelectron_closure` into `photoelectron_density_model` (`none` / `kinetic_mean`) and the independent `photoelectron_histogram_enabled` diagnostic switch. Particle return remains controlled only by `outer_plasma.return_model` and `coupling.particle_transfer_mode`.
- Renamed the tracked photoelectron-return example to `examples/periodic2_photoelectron_return.toml`.

### Removed
- **BREAKING**: Removed `photo_escape_model` and its `boltzmann_cutoff` reduced-emission closure. Photoelectrons now keep their full emitted weight and use ordinary tracking plus `potential_barrier` or an outer-plasma return model.
- **BREAKING**: Removed the `individual_return` photoelectron-closure setting; outgoing histograms and applicability checks are enabled with `outer_plasma.photoelectron_histogram_enabled=true`.

## [1.6.2] - 2026-07-15

### Added
- Opt-in `BEACH_COLLISION_DIAGNOSTICS=1` output for collision-grid traversal failures.

### Fixed
- Collision-grid DDA now classifies subnormal periodic-seam displacements as stationary using a mesh-scale tolerance, preventing false `grid_stalled` failures caused by an infinite axis increment.

## [1.6.1] - 2026-07-15

### Changed
- Package metadata now uses an SPDX license expression and declares the license file through the current `pyproject.toml` fields.

### Fixed
- PyPI source installations now include Fortran benchmark sources referenced by `fpm.toml`.
- Pull-request and release workflows now validate every explicit fpm target in the built source distribution and build a wheel from that artifact before publication.

## [1.6.0] - 2026-07-14

### Added
- Free-space `triangle_p0` Treecode evaluation with exact analytic panel near interactions, all-vertex MAC radii, and monopole far interactions for both electric field and potential.
- A machine-readable output manifest and documentation-contract tests for configuration, output, and restart guidance.

### Changed
- Treecode now accelerates arbitrary-point and mesh-centroid potential evaluation in addition to electric-field evaluation.
- Field-solver documentation now provides a canonical compatibility matrix for solver, kernel, boundary condition, and periodic far-correction combinations.
- Parameter, schema, and agent guidance now document the supported Direct `periodic2` split-reference exception without repeating cross-reference prose.
- Output and restart documentation now describes conditional outer-plasma and photoelectron checkpoints, MPI-global residuals, and OBJ mesh-source paths consistently.
- Starlight development pages now identify main-branch documentation, expose last-updated and edit links, and receive content/type checks in CI.

### Fixed
- Repaired the malformed `[coupling]` parameter table and documented every top-level configuration table, including `coupling.update_mode`.

## [1.5.0] - 2026-07-11

### Added
- `sim.open_boundary_model="potential_barrier"` for open-boundary return/escape decisions based on normal kinetic energy and the local potential barrier.
- Direct continuous triangle-P0 panel potential/field kernels with explicit vacuum-side, self-potential, principal-value, and field-jump contracts.
- Panel-aware FMM topology, near correction, exact triangle multipole moments, and a matching public C/Python kernel path.
- Explicit periodic2 zero-mode projection and panel-spectral nonzero-mode reference evaluation for non-neutral surface charge.
- Linear Debye, nonlinear 1D kinetic, and unified linear-response outer-plasma models with fail-closed applicability and branch diagnostics.
- Conservative particle transfer across the outer interface, including infinity-to-interface reservoir mapping, individual photoelectron return accounting, and explicit 3D electrostatic outer orbits.
- A versioned cached infinite-periodic nonzero-mode operator with checksum validation, corruption recovery, MPI root-only generation, and warm-cache reuse.
- Checkpoint schema 3, which restores complete outer-plasma potential, electric-field, charge-density, residual, and current state while retaining schema-2 read-only migration.
- Physics release gates covering L3, far-correction, two-rank MPI, cache concurrency, memory budgets, and six required convergence categories.
- Bilingual installation, tutorial, validation, troubleshooting, and physical-redesign completion-audit documentation.

### Changed
- `beach-inspect` now uses only precomputed mesh potential for its ordinary summary; `--recompute-potential` explicitly requests the quadratic reconstruction path.
- Periodic far-correction documentation now consistently describes `none` as the default, `auto` as its compatibility alias, and `m2l_root_oracle` as an explicit diagnostic mode.
- `outer_plasma.interface_z` is an ownership/handoff plane in the unified model rather than a truncation of the electrostatic field domain.
- The Starlight documentation navigation now follows installation, first run, output inspection, and numerical-validation workflows.

### Fixed
- Field-solver submodules no longer re-import host-associated FMM types, restoring editable builds with GFortran 13.
- Boris particle candidates now use same-time trapezoidal position updates and sample the boundary-element electric field at the predicted midpoint. The existing `boris_push` signature is unchanged, but corrected trajectories differ from previous releases.
- Particle tracking now orders the first mesh hit against the first box-face event, advances one reflected/periodic remainder, and fails closed when one outer step would cross a second box face.
- Reservoir injection now rounds one global expected macro-particle count before distributing particles across MPI ranks, making count and residual sequences independent of MPI world size.
- MPI checkpoints now store one root-owned global reservoir residual while retaining rank-local RNG state; ambiguous legacy rank-local residual files are rejected on resume.
- The workload estimator now rounds reservoir counts globally before MPI rank splitting and reports both global and selected-rank totals.
- The production particle loop now keeps box-event state off the no-crossing hot path; a representative 800-triangle direct-field benchmark measured 1.13% median overhead versus the pre-event baseline.
- Output-directory creation now uses the POSIX filesystem binding without shell evaluation, including paths with literal shell metacharacters.
- Periodic collision queries now fail closed on image/range limits and propagate deterministic OpenMP/MPI failure context for both particle tracking and `photo_raycast` injection.
- Treecode now descends mixed-sign charge nodes instead of accepting an unreliable near-cancelled monopole approximation.
- Charge-history readers now reject incomplete dense batch snapshots before exposing the batch through indexed or per-step access.

### Removed
- Removed the user-facing output-generation subcommand from `beachx config`; high-level BEACH authoring keys are now accepted directly by the Fortran config loader.

## [1.4.0] - 2026-07-07

### Added
- `output.restart_from` for reading resume checkpoints from a separate directory while writing new outputs to `output.dir`
- `photo_raycast` example species in the default `beachx config init` template
- `photo_escape_model="boltzmann_cutoff"` for `photo_raycast` species to approximate photoelectron return/escape with a reduced Boltzmann closure

### Changed
- README now focuses on quick start and links detailed usage to `docs/`
- `periodic2` far correction now defaults to `none`; `auto` remains a compatibility alias for `none`, and `m2l_root_oracle` is an explicit opt-in diagnostic mode
- Moved `m2l_root_oracle` far-correction diagnostics out of the default heavy Fortran target list and behind `make test-fortran-far-correction`

### Fixed
- Maxwellian velocity sampling now bounds extreme tails to avoid runaway particle steps
- `periodic2` collision search now guards against runaway periodic image enumeration

## [1.3.0] - 2026-06-22

### Added
- `beachx lint <beach.toml>` for TOML parsing, packaged JSON Schema validation, high-level config normalization, and semantic BEACH config checks

### Removed
- **BREAKING**: Removed the preset/case config layer, including `beachx preset`, `case.toml`, preset package data, and case/preset schemas

### Changed
- Fortran config loading now uses `toml-f` for standard TOML syntax, including multiline arrays, dotted keys, inline tables, and literal `#` characters in strings
- Python config normalization now uses `tomli-w` instead of maintaining duplicate in-project TOML writer code
- `beachx config` now creates, validates, normalizes, and diffs direct `beach.toml` files
- BEACH context plugin references now document the direct `beach.toml` workflow

### Fixed
- Serial builds now reject multi-rank launcher configuration before entering MPI initialization

## [1.2.0] - 2026-06-14

### Added
- Codex context plugin for BEACH repository workflows
- Tiered developer test targets (`test-l0`, `test-l1`, `test-l2`, `test-l3`, `test-heavy`, `test-full`)
- Object surface material models, including limited free-space floating conductor support
- Stable development version metadata mode for faster incremental fpm builds

### Changed
- **BREAKING**: `sim.batch_count` is now the cumulative target batch count when `output.resume=true`; resuming from `batches=100` with `batch_count=150` runs 50 more batches
- Resume progress, history, and workload estimates now use remaining batches derived from checkpoint state
- Surface material validation and documentation are aligned with current conductor/dielectric support

### Fixed
- Resume now fails when required checkpoint files are missing instead of silently starting a new run
- Restart, config, CLI, and Python result readers now reject non-finite or invalid numeric values more consistently
- Zero-softening self-field singularities are skipped in direct and tree field paths
- Workload estimator accepts the current normalization and sheath configuration keys

## [1.1.0] - 2026-05-27

### Added
- Internal field normalization controls (`field_normalization`, `field_length_scale`) for direct, treecode, and FMM Coulomb kernels
- External electric-field and velocity-grid injection support
- Field-kernel scene analysis utilities
- Batch-duration stability theory documentation
- Agent user guide for AI-assisted simulation workflows

### Changed
- Periodic mesh visualization can preserve whole mesh objects while wrapping periodic domains
- Coulomb force matrix plotting supports parallel workers for larger mesh-group analyses
- JSON Schema numeric constraints now consistently reject nonzero-only invalid values

### Fixed
- Small-scale collision tolerances for regolith-scale geometries

## [1.0.0] - 2026-03-26

### Added
- `beachx` unified CLI with subcommands: `inspect`, `animate`, `workload`, `slices`, `coulomb`, `mobility`, `profile`, `config`, `preset`, `model`
- Preset-based config workflow with validation and preset management commands
- JSON Schema for `beach.toml` config validation
- OBJ mesh transform support (scale, rotation, offset) with CRLF line ending handling
- `periodic2` support for Coulomb force calculation and 3D electric field line plotting
- Periodic mesh replication for visualization (`periodic2_repeat`)
- Fortran potential history output (`write_potential_history`, `potential_history.csv`)
- OpenMP parallelization for mesh potential computation and FMM evaluation
- `beachx model close-pack` for generating close-packed sphere models
- Sheath injection model (`zhao_auto`) with configurable parameters
- 3D electric field line tracing and plotting (`trace_field_lines`, `plot_field_lines_3d`)
- Coulomb mobility analysis (`analyze_coulomb_mobility`)
- `Beach` facade class for high-level Python API
- Performance profiling output and `beachx profile` visualization
- Comprehensive Fortran test suite (20 test targets)
- Python post-processing API documentation (`docs/python_postprocess_api.md`)
- SPEC.md behavioral specification

### Changed
- **BREAKING**: Old CLI entry points (`beach-inspect`, `beach-animate-history`, etc.) are deprecated in favor of `beachx` subcommands
- `sim_stats` particle counters upgraded from 32-bit to 64-bit integers to prevent overflow in large simulations
- `total_particles_from_config` now validates against integer overflow
- `plot_potential()` default colormap unified to `"viridis"` (was `"jet"` in facade)
- Internal `_select_frame_columns` and `_surface_charge_density` removed from public `__all__`

### Fixed
- Empty CSV handling in `load_fortran_result` for edge cases with zero elements
- Empty `mesh_triangles.csv` no longer causes `IndexError`
- Invalid `m2l_root_trunc` reference removed from example config comments

## [0.8.0] - 2026-03-08

### Added
- Preset-based `beachx config` workflow with high-level spatial config normalization
- JSON Schema for case and preset TOML files
- `beachx preset save` and `beachx preset edit` commands
- High-level config validation

### Changed
- Removed unused reserved sim keys

## [0.7.0] - 2026-02-15

### Added
- `periodic2` m2l_root_oracle far correction (Ewald residual build-time fit)
- Exact periodic2 Ewald far correction
- Mesh view angle controls for Python plotting

### Changed
- Reorganized FMM source layout into subdirectories
- Switched Fortran formatter to 2-space indent

## [0.6.0] - 2026-01-25

### Added
- Plate-hole template and per-cap cylinder controls
- Simulation-box potential slice plotting with configurable vmin/vmax
- Expanded Fortran test cases

### Changed
- Split large Fortran modules into focused submodules (simulator, field_solver, config_parser)
- Organized split modules into subdirectories

## [0.5.0] - 2026-01-12

### Added
- Treecode field solver with auto-tuning heuristics
- Treecode sim config controls (`tree_theta`, `tree_leaf_max`, `tree_min_nelem`)
- Treecode parameter support in workload estimate CLI

### Fixed
- Photo raycast now ignores out-of-box mesh hits

## [0.4.0] - 2025-12-28

### Added
- Mesh-group charge APIs and Coulomb interaction calculation
- `rel_change` displayed in batch progress output

### Changed
- **BREAKING**: Simplified Fortran history API
- Renamed `calc_coulomb` args to `target`/`source`
- Default mesh step selection to latest history

## [0.3.0] - 2025-12-15

### Changed
- **BREAKING**: Removed deprecated config keys

## [0.2.0] - 2025-12-01

### Added
- pip-installable packaging with bundled Fortran binary
- MPI-rank-aware workload estimation
- Generic HPC job script template
