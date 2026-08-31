title: Use FMM

Lang: [日本語](FMM.md) | [English](FMM.en.md)

# Use FMM

This how-to helps BEACH users select FMM, configure its computational domain, and check both error against Direct and
measured runtime. FMM accelerates repeated evaluation at many particle positions by approximating the distant field from
large triangle meshes. Its fixed costs can dominate small cases, so adopt it only after checking accuracy and runtime.

After this page, you can configure free-space FMM, tune its accuracy, compare it with Direct, and profile one case.
[Coulomb FMM internals](FMMCore.en.html) is the source of truth for expansion formulas, internal data structures, update
algorithms, and parallel implementation.

## Decide whether to use FMM

| Condition | Choice |
| --- | --- |
| Many elements and many particle steps per batch | Consider FMM and benchmark it against Direct in a release build |
| Production infinite two-axis periodic field | Use FMM with `cached_kneq0` |
| Small case, reference result, or approximation check | Use Direct |
| Automatic selection by element count in free space | Use `field_solver="auto"` and check the resolved solver in output |

The crossover depends on particle count, tracking steps, and target distribution as well as element count.
See [Field evaluation](FieldSolvers.en.html) for the full choice including Treecode.

## Use a minimal free-space FMM configuration

Start from a case that already passes lint, then set `[sim]`, `[domain]`, `[field_boundary]`, and `[output]` as below.
The box values are examples; use the physical extent of the actual case.

```toml
[sim]
field_solver = "fmm"

[domain]
box_min = [-0.5, -0.5, -0.1]
box_max = [ 0.5,  0.5,  1.0]
periodic_axes = []

[field_boundary]
mode = "free"

[output]
write_files = true
dir = "outputs/fmm"
```

See [`examples/panel_fmm.toml`](../examples/panel_fmm.toml) for a complete case. That short-check example has
`write_files=false`; change it to the `[output]` settings above before inspecting results. Then validate the input, run it,
and verify the resolved solver:

```bash
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/fmm
grep '^field_reconstruction_resolved_field_solver=' outputs/fmm/summary.txt
```

The lint output must end with `status=ok`, the simulator must exit with code `0`, and the final command must print
`field_reconstruction_resolved_field_solver=fmm`. This establishes completion through the FMM path; it does not establish
acceptable approximation error or physical validity.

## Cover the evaluation range with the domain

`[domain]` defines both the particle region and the target tree where FMM provides fast evaluation. Enclose every position
that particles can reach. When writing element-center potentials or history, keeping the mesh inside the box also avoids
out-of-range evaluations.

Specify either `box_min` with `box_max`, or `box_origin` with a positive `box_size`. In a free field, BEACH falls back to a
Direct sum over every element for a point outside the box. The result remains available, but frequent fallback removes the
performance benefit. A periodic2 run with `cached_kneq0` rejects an out-of-box target.

## Tune accuracy

First omit `tree_theta` and `tree_leaf_max` and use the element-count-dependent resolved values as a baseline.
This automatic tuning also applies when `field_solver="fmm"` is explicit. Inspect the resolved values with:

```bash
grep -E '^field_reconstruction_(tree_theta|tree_leaf_max|fmm_expansion_order)=' \
  outputs/fmm/summary.txt
```

| Setting | How to use it |
| --- | --- |
| `tree_theta` | Smaller values make far acceptance stricter and are generally more accurate and slower; `0 < theta <= 1` |
| `tree_leaf_max` | Changes the balance between tree depth and near Direct work; compare multiple values of at least `1` |
| `field_normalization` | Changes internal coordinate scaling only; do not use it as a physical or accuracy acceptance control |

The current simulator fixes the expansion order at four; it is not configurable. Tune in this order:

1. Run with the automatically resolved values.
2. Add a run with a smaller `tree_theta`.
3. Vary `tree_leaf_max` and find where the observables are stable.
4. Accept a combination only when it meets the Direct error tolerance and reduces measured runtime.

See [Input parameters](Parameters.en.html#field-solver) for types, defaults, and the element-count tuning table.

<a id="source-kernel"></a>

## Check accuracy against Direct

Direct and FMM both represent element charge with the same triangle-P0 source model. Their difference on the same mesh
therefore measures mainly the FMM approximation. Check mesh-discretization error separately by refining the mesh.

Duplicate a reduced free-space case and change only the solver and output directory. For the first solver-error comparison,
set `sim.batch_count=1` so both runs produce element charge from the same initial surface-charge state.

| Configuration | `sim.field_solver` | `output.dir` |
| --- | --- | --- |
| `direct.toml` | `"direct"` | `"outputs/direct"` |
| `fmm.toml` | `"fmm"` | `"outputs/fmm"` |

Keep the mesh, particles, random seed, thread count, `dt`, `field_normalization`, and output controls identical.

```bash
beachx lint direct.toml
beachx lint fmm.toml
OMP_NUM_THREADS=1 beach direct.toml
OMP_NUM_THREADS=1 beach fmm.toml
beachx inspect outputs/direct
beachx inspect outputs/fmm
```

Use [Electric Field Calculation](PythonPostprocessAPI.en.html#6-electric-field-calculation) and
[Potential Reconstruction](PythonPostprocessAPI.en.html#4-potential-reconstruction) to compare the same evaluation points.
Relative error becomes unstable where the reference field is nearly zero, so also use absolute error or normalize by a
representative field. After the solver error passes, extend to a representative multi-batch run and compare absorption and
escape counts, element-wise final charge, and `charge_ledger.csv`. Set tolerances from the research objective before comparing
and follow [Validating simulation results](ValidationGuide.en.html) for acceptance.

Do not use a target exactly on a panel for an ordinary pointwise field agreement check because FMM and Direct use different
trace conventions there. A periodic2 Direct reference cannot be created by changing only the solver name; use the split
reference described in [periodic2 electrostatics](PeriodicElectrostatics.en.html).

## Measure performance

After accuracy passes, measure a release build with a representative particle count. Keep MPI ranks, OpenMP threads, mesh,
particles, and output controls identical between Direct and FMM, and run locally or inside a compute-node allocation.

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach direct.toml
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach fmm.toml
beachx profile outputs/direct/performance_profile.csv \
  --save outputs/direct/performance_profile.png
beachx profile outputs/fmm/performance_profile.csv \
  --save outputs/fmm/performance_profile.png
```

Use `rank_max_s` on the `simulation_total` row for the total comparison. FMM setup appears in `field_solver_init`, per-batch
field refresh in `field_refresh`, and particle tracking including target evaluation in `particle_batch`. Setup cost can make
FMM slower on a small case.

Output controls such as `write_potential_history` also change workload and must match. For periodic2 `cached_kneq0`, measure
cold-cache and warm-cache runs separately. See [Runtime profiling](Execution.en.html#profiling) for phase definitions.

## Use FMM with periodic2

FMM is the normal production path for `periodic2`. Use `cached_kneq0` for infinite-periodic nonzero modes and specify the
physical zero mode under `[periodic2]`. `field_periodic_far_correction="none"` is a finite-image approximation, not an
infinite-periodic solution.

See [`examples/periodic2_cached_panel.toml`](../examples/periodic2_cached_panel.toml) for a complete configuration.
[periodic2 electrostatics](PeriodicElectrostatics.en.html) is canonical for component selection and validation, while
[periodic2 Far Correction](PeriodicFarCorrection.en.html) owns cache generation and reuse.

## Check current limits

- FMM sources are Coulomb triangle-P0 elements; no alternative kernel can be selected.
- The simulator expansion order is fixed at four and cannot be swept from configuration.
- Supported field boundaries are `free` and `periodic2`; follow the complete
  [compatibility table](FieldSolvers.en.html#solver-and-field-boundary-compatibility).
- Mesh geometry remains fixed during a run; only charge changes between batches.
- Changing `field_normalization` does not change the SI units of output fields and potentials.
- FMM does not solve an outer plasma or sheath; matching-plane response is composed outside FMM.

Developers looking for formulas, internal APIs, geometry/charge-state separation, out-of-range fallback, or OpenMP details
should continue with [Coulomb FMM internals](FMMCore.en.html).
