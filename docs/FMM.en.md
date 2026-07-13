title: Field evaluation with FMM

Lang: [English](FMM.en.md) | [日本語](FMM.md)

# Field evaluation with FMM

The Fast Multipole Method accelerates Coulomb-field evaluation from many boundary elements at many particle
targets. This page contains what a user needs to select the solver and verify its accuracy.

## When to use it

| Case | Recommendation |
| --- | --- |
| Small mesh or reference solution | `direct` |
| Medium free-boundary case | `auto` or `treecode` |
| Large mesh with many targets | `fmm` |
| `periodic2` | `fmm` is required |

FMM is approximate. For a new mesh, kernel, or periodic configuration, first compare a small case with direct
evaluation or an independent oracle.

## Computation

FMM partitions sources in an octree and summarizes well-separated cell interactions with expansions.

```text
source charge
   ↓ P2M
multipole expansion
   ↓ M2M
upward tree
   ↓ M2L
local expansion
   ↓ L2L / L2P
target field
```

Near cells remain direct, preserving local accuracy while aggregating far interactions.

## Plan and state

BEACH separates fixed mesh geometry from element charge that changes each batch.

- Plan: source positions, tree, interaction lists, and other geometry-dependent data
- State: multipole and local coefficients built from current `q_elem`

The plan is normally built once during initialization, while only state is updated each batch. A geometry or
element-count change requires rebuilding the plan.

## Source kernels

| `field.element_kernel` | FMM source |
| --- | --- |
| `point` | Softened point charge at the element centroid |
| `triangle_p0` | Constant density over a triangle |

These are different discretizations. Check source-mesh and kernel refinement independently of FMM order.

## Accuracy controls

- Expansion order
- Leaf size and well-separated criterion
- Point-kernel softening
- `triangle_p0` near correction
- Periodic image layers and far correction
- Source-mesh resolution

A small FMM truncation error at one order does not establish a small source-discretization error.

## periodic2

Under `periodic2`, finite near images from FMM can be combined with the `cached_kneq0` far operator. The zero mode
and outer response are composed by the electrostatic snapshot rather than completed inside the FMM core alone.

See [periodic2 field evaluation](PeriodicElectrostatics.en.html).

## Recommended verification

1. Compare FMM with direct evaluation on a small mesh.
2. Increase FMM order and check field and potential convergence.
3. Refine the source mesh and check conclusion stability.
4. Under periodic2, check cold/warm cache equality and zero-mode diagnostics.
5. Measure performance with a release build.

See [Coulomb FMM core implementation](FMMCore.en.html) for internal APIs, tree construction, Cartesian
expansions, M2L, and cached operators.
