title: reservoir_face inflow and velocity sampling

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# `reservoir_face` inflow and velocity sampling

This page describes the injection algorithm used after selecting `source_mode="reservoir_face"` on the
[particle-source overview](ParticleSourcesBoundaries.en.html). It answers two questions:

1. How does an upstream velocity distribution function (VDF) determine physical and macro-particle inflow per batch?
2. How are initial positions derived from the injection aperture and face-arrival velocities generated for the selected particles?

```text
Upstream VDF and aperture
          ↓ determine accessible flux
Inflow during one batch
          ↓ convert to a macro-particle count
Initial position from the aperture and face-arrival velocity
          ↓
Common particle tracker
```

`reservoir_face` itself does not select the source mode, choose a sheath model, or process escape and return after a particle
leaves the box. Those responsibilities are separated into the related pages listed below.

## Inputs and outputs of this page

| Kind | Content |
| --- | --- |
| Inputs | Upstream VDF, injection aperture, `batch_duration`, macro-particle weight, and optionally the potential drop from upstream to the face |
| Outputs | Global macro-particle count, initial positions derived from the injection aperture, face-arrival velocities, and the remainder carried into the next batch |

The upstream VDF can be a Maxwell distribution or a velocity grid. Both are interpreted as surface-crossing inflow; a
potential-drop correction is added only when the selected model requires it.

## Define the aperture and inflow area

`inject_face` is one of `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, or `z_high`. `pos_low` and `pos_high` define a
rectangular aperture on that box face. Let $\mathbf n$ be its inward unit normal and $A$ its area.

`reservoir_face` requires `sim.use_box=true` and `sim.batch_duration>0`. The normal coordinates of the aperture must coincide
with the box face, and both tangential extents must be positive.

## Weight a Maxwell distribution by inflow flux

For an upstream drifting Maxwell distribution, define normal drift $u_n=\mathbf u\cdot\mathbf n$ and thermal standard
deviation $\sigma=\sqrt{k_\mathrm{B}T/m}$. With upstream normal-speed floor $v_{\min}$, the one-way inflow flux is

$$
\Gamma_\mathrm{in}
=n\int_{v_n\ge v_{\min}}v_n f_n(v_n)\,dv_n
=n\left[u_n\{1-\Phi(a)\}+\sigma\varphi(a)\right],
\qquad
a=\frac{v_{\min}-u_n}{\sigma},
$$

where $\Phi$ and $\varphi$ are the standard-normal CDF and density. At zero temperature,
$\Gamma_\mathrm{in}=nu_n$ when $u_n\ge v_{\min}$ and zero otherwise.

Tangential velocities are ordinary Gaussian samples. Normal velocity is sampled from the flux-weighted half-range distribution
proportional to $v_n f_n(v_n)$. This accounts for the probability of crossing the face being proportional to normal velocity,
rather than placing a volume Maxwell distribution directly on the face.

## Convert velocity-grid values into flux

With `velocity_distribution="grid"`, velocities are sampled from CSV points and nonnegative values $f$.

| `velocity_grid_pdf_kind` | Meaning of CSV $f$ | Sampling weight |
| --- | --- | --- |
| `phase_space` | Phase-space distribution | $\max(v_n,0)f$ |
| `flux_weighted` | Already weighted for inflowing particles | $f$ |

Both retain only entries with $v_n\ge v_{\min}$ and $v_n>0$. `velocity_grid_sampling="rectilinear"` interpolates a complete
rectilinear grid, `discrete` samples input points directly, and `auto` uses `rectilinear` when possible.

For a velocity-grid distribution, count comes from `particle_flux_m2_s` or `current_density_a_m2`; density and temperature do not
set the count. The current implementation rejects combination with Zhao injection correction.

## Carry residuals forward to preserve long-time flux

Let $\Gamma_\mathrm{in}$ be either the flux derived from a Maxwell distribution or the flux supplied for a velocity grid. The
expected physical and macro-particle inflows are

$$
N_\mathrm{phys}=\Gamma_\mathrm{in}A\,\Delta t_\mathrm{batch},
\qquad
N_\mathrm{macro,expected}=\frac{N_\mathrm{phys}}{w}.
$$

If $r$ is the fractional remainder carried from the preceding batch,

$$
N_\mathrm{macro}=\left\lfloor r+N_\mathrm{macro,expected}\right\rfloor,
\qquad
r\leftarrow r+N_\mathrm{macro,expected}-N_\mathrm{macro}.
$$

Individual batch counts can vary, while the total over many batches approaches the expected flux.

Supply either a positive `w_particle` or derive it from `target_macro_particles_per_batch`. The latter adjusts the number of
tracked samples without changing physical flux. A value of `-1` shares the weight resolved for species 1. The two forms are
mutually exclusive.

## Use one potential drop for accessibility and face velocity

Let $\phi_\infty$ be the potential at infinity or upstream and $\phi_f$ the injection-face or interface potential. The 1-D
electrostatic map preserves tangential velocity and normal energy:

$$
\frac12m v_{n,f}^2+q\phi_f
=\frac12m v_{n,\infty}^2+q\phi_\infty,
$$

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}.
$$

| $B$ | Upstream accessibility and face velocity |
| ---: | --- |
| $B>0$ | Only $v_{n,\infty}\ge\sqrt B$ reaches the face and decelerates on the way |
| $B=0$ | Normal velocity is unchanged |
| $B<0$ | Every inflowing particle is accessible and accelerates toward the face |

The upstream VDF is truncated to accessible particles with $v_{\min}=\sqrt{\max(B,0)}$, and accepted velocities are mapped to
the face with the same $B$. One potential drop therefore determines count and face-arrival velocity consistently.

## Receive the potential drop or upstream-VDF correction from another model

This page converts a supplied upstream VDF and potential drop into particles. Use
[Selecting boundary and outer-domain models](OuterPlasmaModels.en.html) as the canonical guide for choosing the external model.

| Configuration | Value supplied to `reservoir_face` |
| --- | --- |
| No correction | Use $B=0$ and interpret the configured VDF as the distribution at the face |
| `reservoir_potential_model="infinity_barrier"` | One potential drop derived from the mean aperture potential and `sim.phi_infty` |
| Split `linear_debye` | Interface potential drop derived from the surface zero mode |
| Split `kinetic_1d` | $\phi_I-\phi_\infty$ from the converged outer profile |

`infinity_barrier` evaluates the batch-start potential on an `injection_face_phi_grid_n` by
`injection_face_phi_grid_n` cell-centered grid in the aperture and uses the mean $\bar\phi_f$. It follows the same potential
convention as the field snapshot, including the point or triangle kernel, periodic field, zero mode, outer state, and `sim.e0`,
but does not solve intermediate $E(z)$, turning position, flight time, or space charge. The same evaluation accumulates the
population standard deviation, minimum, and maximum. For a Maxwell reservoir, a large in-face variation relative to its
characteristic energy produces a warning in the first and final batch.

In a split outer model, z-high `reservoir_face` species are interpreted as VDFs at infinity. For `kinetic_1d`, the interface field
is a Poisson boundary condition, while particle velocity changes through $\phi_I-\phi_\infty$.

Zhao-family models and `floating_no_photo` do not supply a generic potential drop $B$; they correct upstream-VDF density, drift,
and speed floor. See [Inflow-VDF sheath closures](SheathInjectionClosures.en.html) for those rules.

## Disperse initial positions with a virtual flight interval

Creating every particle at the same face-time position introduces an artificial alignment. BEACH assigns each particle a uniform
virtual interval $\tau\in[0,\texttt{sim.dt})$ and shifts only its initial position before tracking:
$\mathbf x\leftarrow\mathbf x+\mathbf v\tau$.

Internally, the runtime passes this interval as `position_jitter_dt=sim.dt` to the sampler. This is an internal identifier, not a
user-facing configuration key. The position shift does not add random noise to velocity. It preserves global simulation time and each particle's tracking horizon. Axes periodic on both faces are
wrapped into the primary box; other axes are clamped inside the box.

## Preserve the global remainder across MPI and restart

With MPI, the root rank determines the global count and remainder before splitting particles across ranks. The remainder is stored
in `macro_residuals.csv` and broadcast to all ranks on restart. Expected inflow and the remainder sequence therefore do not depend
on MPI world size.

When checking results, distinguish species-resolved injected count and charge, resolved macro-particle weight, remainder, applied
$v_{\min}$ and potential drop, and batch-averaged absorption and escape currents.

Changing `batch_duration` changes both expected particles per batch and the field-update interval. See
[Batch duration and stability](BatchDurationStability.en.html) for convergence of physical steady results.

## Related pages

- Choose a particle source: [Particle-source overview](ParticleSourcesBoundaries.en.html)
- Combine inflow correction, outflow boundaries, and outer fields: [Selecting boundary and outer-domain models](OuterPlasmaModels.en.html)
- Correct the VDF with Zhao-family models or `floating_no_photo`: [Inflow-VDF sheath closures](SheathInjectionClosures.en.html)
- Process particles that leave through an open face: [Particle escape and return](ParticleEscapeReturn.en.html)

## Code reference

- Flux integration, macro-particle count, and Maxwell or velocity-grid sampling: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Velocity corrections from barriers and outer profiles: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Aperture geometry and input-combination validation: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- Zhao species correction: [`bem_sheath_runtime.f90`](../src/physics/sheath/bem_sheath_runtime.f90)
- MPI-global macro-particle remainder: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
