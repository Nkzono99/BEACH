title: periodic2 electrostatics

Lang: [English](PeriodicElectrostatics.en.md) | [日本語](PeriodicElectrostatics.md)

# periodic2 electrostatics

`fields.boundary="periodic2"` computes the field of a surface that repeats infinitely in x/y and does not repeat in z. Use it
when one cell represents part of a wide surface such as regolith. BEACH splits this field into three components and adds each
once. This page explains what each component represents and which settings compute it.

## What is solved

The field is that of the surface charge in the cell and all its copies in x/y. The sum over copies is split into three components
with different character.

| Component | Physical meaning | How it is computed |
|---|---|---|
| Near copies | Strongly varying local field of the cell and the $N$ layers of copies around it | Direct sum over finite images (Direct or FMM) |
| Nonzero mode of far copies ($k\ne0$) | Smooth field varying in x/y from copies beyond the $N$ layers | Far correction (`cached_kneq0`); not included with finite images only |
| Zero mode ($k=0$) | Field averaged over x/y, set by the total charge below each height | Computed analytically from the height distribution of the triangles |

The zero mode is not a component that may be removed numerically; it is a physical component set by Gauss's law and the lower
boundary condition.

## Configuration

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0e-4, 1.0e-4, 1.0e-3]
periodic_axes = ["x", "y"]

[fields]
boundary = "periodic2"

[fields.solver]
method = "fmm"

[fields.periodic]
backend = "cached_kneq0"
image_layers = 1
lower_boundary_model = "symmetric_vacuum"
```

The periodic axes are set by `domain.periodic_axes` and are shared by the field, the particles, and the illumination rays.
`periodic2` accepts only the combination of periodic x/y and non-periodic z.

### Choose a backend

| `fields.periodic.backend` | Components included | Use |
|---|---|---|
| `finite_images` | Near copies only; the zero mode is not separated | Finite-image models and small comparisons |
| `cached_kneq0` | All three components | Regular calculations with FMM |
| `panel_spectral_reference` | All three components (the nonzero mode by a direct Fourier sum) | Small reference solutions with Direct |

Results with `finite_images` cannot be treated as the solution for an infinitely repeating surface until convergence in the
number of image layers is confirmed. Combinations of solvers and field boundaries are in [Choose a field solver](FieldSolvers.en.html).

## Finite images

`fields.periodic.image_layers` $=N$ adds the $N$ layers of copies $(i,j)\in[-N,N]^2$ around the cell directly.

| $N$ | Cells added |
|---|---|
| 0 | The cell only ($1\times1$) |
| 1 | One surrounding layer ($3\times3$) |
| 2 | Two surrounding layers ($5\times5$) |

The far correction adds only the copies outside these $N$ layers. The layers summed directly by FMM must match the layers the far
correction subtracts, and the far-correction operator is built for each number of layers.

## Nonzero mode of far copies

`cached_kneq0` takes the sum over infinitely repeating copies (the Ewald sum), removes the contribution of the $N$ directly summed
layers and the zero mode, and adds the result to FMM with a precomputed operator. The zero mode is removed so that it is not
counted twice with the analytic zero mode below. How the operator is built, the cache, and the accuracy are in
[periodic2 far correction](PeriodicFarCorrection.en.html).

## Zero mode ($k=0$)

With total charge $q_i$ of triangle $i$ and the fraction $F_i(z)$ of its area at or below height $z$, the total charge below $z$ is

$$
C(z)=\sum_iq_iF_i(z)
$$

$F_i$ is piecewise quadratic between the heights of the three vertices, so $C(z)$ is also quadratic in each interval. With cell
area $A=L_xL_y$ and far field below the surface $E_\mathrm{bottom}$, Gauss's law gives

$$
E_0(z)=E_\mathrm{bottom}+\frac{C(z)}{\epsilon_0A}
$$

The potential from a reference point $(z_g,\phi_g)$ is

$$
\phi_0(z)=\phi_g-E_\mathrm{bottom}(z-z_g)
-\frac1{\epsilon_0A}\int_{z_g}^zC(\zeta)\,d\zeta
$$

In a cell whose total charge is not zero, a constant field and a linear potential remain above the surface.

### Lower boundary condition

| `lower_boundary_model` | $E_\mathrm{bottom}$ | $E_\mathrm{top}$ | Meaning |
|---|---:|---:|---|
| `symmetric_vacuum` | $-Q/(2\epsilon_0A)$ | $+Q/(2\epsilon_0A)$ | The same vacuum extends above and below, with no external field |
| `e_bottom_zero` | $0$ | $Q/(\epsilon_0A)$ | The field is zero below the lowest charge |

$Q=\sum_iq_i$ is the total charge of the cell. If the total charge is zero, the two conditions coincide. Neither solves dielectric
polarization or shielding inside objects.

### Potential reference

With the outer-sheath connection, the reference point of the zero mode is the top face $(z_\mathrm{high},\phi_0)$, and potentials
are referenced to the upstream plasma at 0 V ([Connecting to the outer sheath](ZhaoStationaryClosure.en.html#match-the-top-face-potential-to-the-wall-potential)).
Otherwise, read potentials relative to the top-face mean.

### Why the mean field of the outer sheath cannot replace it

Inside a particle layer, $C(z)$ changes strongly when positive and negative charge separate by height, even if the total charge is
small. This mean field inside the layer is also part of the zero mode. Replacing it with the planar field of the outer sheath
drops this component created by the surface charge. With the outer-sheath connection, the mean field inside the layer is computed
from the surface charge, and only the potential reference and the barriers are taken from the outer sheath.

## Particle collisions and periodic images

The field is evaluated at the position mapped back into the cell, while particle orbits are tracked in physical coordinates.
Collisions with surfaces are found by searching geometrically the periodic images an orbit can reach
([Collision and boundary events](ParticleEvents.en.html)).

## Check convergence

- With `finite_images`, increase `image_layers` and confirm that the quantity of interest does not change.
- With `cached_kneq0`, confirm that rebuilding and reusing the cache, and changing the thread and MPI layout, give the same result ([periodic2 far correction](PeriodicFarCorrection.en.html)).
- In a cell whose total charge is not zero, do not read a potential difference at finite height directly as the energy needed to escape to infinity.

## Scope

- Only the two axes x/y are periodic. The z direction does not repeat.
- Dielectric polarization, resistance, and shielding inside objects are not solved.
