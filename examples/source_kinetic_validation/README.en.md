# Reproducing source-connected validation

Lang: [日本語](README.md) | [English](README.en.md)

This example compares Fortran static roots with the independent #38 distribution and probes time dependence with the
existing kinetic oracle. It separates root agreement, finite search results, and dynamic convergence.
See the [model reference](../../docs/SourceKineticSheath.en.md) for physics and CSV contracts.

## Running

Use a BEACH development environment, GNU Fortran, NumPy, and Matplotlib. Running the independent Python code also
requires SciPy and the external `lunar-sheath-solver-v0.1.0.zip` input. These validation dependencies are not added to
BEACH runtime requirements. Run from the repository root on a compute node; on KUDPC, enter a `tssrun` allocation first.

```bash
mkdir -p build/source-validation
gfortran -O2 -J build/source-validation -I build/source-validation \
  src/core/bem_kinds.f90 src/core/bem_constants.f90 \
  src/physics/sheath/bem_source_kinetic_sheath.f90 \
  examples/source_kinetic_validation/scan.f90 -o build/source-validation/scan
python examples/source_kinetic_validation/compare_handoff.py \
  _handoff/lunar-sheath-solver-v0.1.0.zip build/source-validation/scan \
  build/source-validation/archive
unzip _handoff/lunar-sheath-solver-v0.1.0.zip -d build/source-validation/reference
PYTHONPATH=build/source-validation/reference/lunar-sheath-solver/src \
  python -m unittest discover -s build/source-validation/reference/lunar-sheath-solver/tests
PYTHONPATH=build/source-validation/reference/lunar-sheath-solver/src \
  python examples/source_kinetic_validation/multiaxis.py \
  build/source-validation/scan build/source-validation/multiaxis
PYTHONPATH=. python examples/source_kinetic_validation/oracle_probe.py build/source-validation/oracle
```

`compare_handoff.py` writes `native_roots.csv`, a `comparison.json` with differences and checksums, and the E–G
`existence.png` for 1,211 conditions. `multiaxis.py` compares 181 M–E, tau–G, and Type A fixture queries against
`solve_extended`, retaining 481/961-point scans, every root, a comparison JSON, and two figures. Plot cells are samples,
not precise bifurcation curves. Both scripts exit with code 1 on disagreement.

`oracle_probe.py` matches the electron source to each root's $A_e$ and separately changes spatial resolution, velocity
resolution, outer length, ion temperature, duration, and velocity range. It retains full oracle outputs for all fourteen
conditions and an aggregate `comparison.json`. The existing oracle has warm ions; comparing $T_i=0.12$ and $0.03$ eV
probes the cold limit. Exit code 0 means all probes ran, not that every probe was classified `steady`.

## Results on 2026-09-22

Independent ZIP SHA-256:
`7c3a2e885c5df9bf816c0a54c4facafcb35d76381be12b6a966d1200bbfd11f4`.

- Independent Python: 86 tests passed.
- Distribution: status, count, and Type agreed for 739 roots across 1,211 conditions. Maximum scaled root difference: $3.6\times10^{-13}$.
- The 980-point E–G slice had 600 accepted, 255 analytically excluded, and 125 unresolved queries, with 624 detected roots.
- All 181 additional M–E, tau–G, and A fixture queries agreed with the extended Python solver (maximum difference $2.1\times10^{-13}$).
  Refining from 481 to 961 points did not change root count, Type, or classification.
- Final L1 checks based on BEACH's public main (`c3fc093`): all Fortran targets passed; Python had 843 passed and 49 skipped.
  These ran on a KUDPC compute node with `CC=gcc`, including documentation navigation and bundled-reference consistency.
- Seven source orbit/Poisson regression cases with sixty assertions, dark online/table charging, and a two-batch example passed.

Time-dependent comparisons used $M=10,\tau=0.2,G=0.3,E_H=0.1,T_e=12$ eV. Static potentials are 0.274641 V for B and
−1.595019 V for N.

| Changed axis | B potential [V] | B status | N potential [V] | N status |
|---|---:|---|---:|---|
| Base: nz=32, nv=128, L=3λ, 6 ion transits | 0.303067 | far boundary unconverged | −1.385123 | far boundary unconverged |
| nz=64 | 0.305228 | far boundary unconverged | −1.390159 | far boundary unconverged |
| nv=256 | 0.282565 | steady | −1.388728 | far boundary unconverged |
| L=6λ, nz=64 | 0.298639 | far boundary unconverged | −1.574772 | steady |
| Ti=0.03 eV, nv=256 | 0.282565 | steady | −1.388725 | far boundary unconverged |
| 12 ion transits | 0.306552 | far boundary unconverged | −1.385124 | far boundary unconverged |
| vmax=9ve, same Δv | 0.303135 | far boundary unconverged | −1.385129 | far boundary unconverged |

Stability of the static roots, simultaneous refinement of all axes, and table/online convergence of BEACH charging with
PE remain unestablished. Across fourteen probes, the maximum charge-budget residual was $2.01\times10^{-14}$,
Gauss residual $1.01\times10^{-27}$ C/m², and velocity-boundary loss fraction $5.46\times10^{-4}$.
Conservation does not establish far-boundary or grid convergence. Root order is not used to select a stable branch, and unconverged oracle results are not
converted to tables. This implementation does not complete #38's dynamic validation.
