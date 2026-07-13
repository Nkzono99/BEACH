title: periodic2 field evaluation

Lang: [English](PeriodicElectrostatics.en.md) | [日本語](PeriodicElectrostatics.md)

# periodic2 field evaluation

`periodic2` makes x/y periodic while leaving z open. The field is not defined by adding only a finite number of
images; BEACH separates lateral Fourier modes from open-direction boundary closure.

## Field decomposition

| Component | Meaning | Evaluation |
| --- | --- | --- |
| Near images | Strong local field around the primary cell | Finite-image direct/tree/FMM sum |
| `k != 0` far field | Infinite-periodic field varying in x/y | `cached_kneq0` operator |
| Surface `k = 0` | Plane-averaged surface-charge field | Cumulative charge in height |
| Outer response | Plasma beyond the interface | Selected outer model |

The electrostatic snapshot combines each required component exactly once.

## Finite images and cached operator

With `field_periodic_far_correction="none"`, only configured `periodic_images` are evaluated. This is a finite-image
model and must not be interpreted as an infinite-periodic field without image refinement.

`cached_kneq0` precomputes a root-level linear operator for the smooth far field obtained from an Ewald2P
reference after subtracting near images and the surface zero mode. A cache can be reused only for the same
geometry, kernel, softening, periodic cell, and FMM settings.

The cache stores `k != 0` only. The snapshot builder adds the physical zero mode once, avoiding double counting.

## Surface zero mode

Triangle area is projected in height, and cumulative charge below each height produces the plane-average field.
A non-neutral cell can retain a constant far field and linear potential in z, so lower and upper boundary closure
must be explicit.

The zero mode is a physical Gauss-law component, not a numerical term that can simply be discarded.

## Relation to particles and collision

Field targets are wrapped to the periodic cell, while physical trajectory-event positions are retained. Collision
searches the geometric periodic images required by a segment. Field-image and collision-image counts are not the
same numerical concept.

## Selection guide

| Purpose | Configuration |
| --- | --- |
| Small finite-image comparison | `periodic_images` with no far correction |
| Infinite-periodic production | `cached_kneq0` plus explicit zero/outer closure |
| Solver verification | Image refinement, cold/warm cache equality, and Ewald/oracle comparison |

## Diagnostics

- Cache fingerprint and cold/warm equality
- No duplicate primary, near, far, or zero component
- Gauss residual and lower/upper closure
- Convergence against image count or a reference operator
- No interpretation of finite-height potential in a non-neutral cell as escape energy at infinity

See [Field solvers and boundary conditions](FieldSolvers.en.html) for configuration and
[Outer plasma models](OuterPlasmaModels.en.html) for outer closure. Ewald and operator implementation details
remain in the legacy combined [periodic2 zero mode and outer plasma](PeriodicZeroModeOuterPlasma.en.html).
