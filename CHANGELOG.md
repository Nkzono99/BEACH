# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [1.0.0] - 2026-03-26

### Added
- `beachx` unified CLI with subcommands: `inspect`, `animate`, `workload`, `slices`, `coulomb`, `mobility`, `profile`, `config`, `preset`, `model`
- Preset-based config workflow (`beachx config render`, `beachx config validate`, `beachx preset list/show/save/edit`)
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
- Preset-based `beachx config` workflow with high-level spatial config rendering
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
