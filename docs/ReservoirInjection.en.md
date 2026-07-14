title: Reservoir injection

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# Reservoir injection

`source_mode="reservoir_face"` converts a distribution outside the box into a macro-particle flux through a selected face. It
integrates the flux of upstream particles that can reach the face, maps their normal velocity with the same potential drop, and
then samples position and velocity.

## Define the aperture and inflow area

`inject_face` is one of `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, or `z_high`. `pos_low` and `pos_high` define a
rectangular aperture on that box face. Let $\mathbf n$ be its inward unit normal and $A$ its area.

`reservoir_face` requires `sim.use_box=true` and `sim.batch_duration>0`. The normal coordinates of the aperture must coincide
with the box face, and both tangential extents must give a positive area.

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
proportional to $v_n f_n(v_n)$. This is the essential difference from placing a volume Maxwell distribution directly on a face.

## Carry residuals forward to preserve long-time flux

The expected physical and macro-particle inflows are

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

Individual batch counts may vary, but the total over many batches approaches the expected flux.

Supply either a positive `w_particle` or derive it from `target_macro_particles_per_batch`. The latter controls sampling count
without changing physical flux. A value of `-1` shares the weight resolved for species 1. The two forms are mutually exclusive.

## Convert velocity-grid values into flux

With `velocity_distribution="grid"`, velocities are sampled from CSV points and nonnegative values $f$.

| `velocity_grid_pdf_kind` | Meaning of CSV $f$ | Sampling weight |
| --- | --- | --- |
| `phase_space` | Phase-space distribution | $\max(v_n,0)f$ |
| `flux_weighted` | Already weighted for inflowing particles | $f$ |

Both retain only entries with $v_n\ge v_{\min}$ and $v_n>0$. `velocity_grid_sampling` selects interpolation on a complete
rectilinear grid, discrete sampling of input points, or `auto`, which uses interpolation when possible.

For a grid distribution, count comes from `particle_flux_m2_s` or `current_density_a_m2`; density and temperature are unused.
The current implementation rejects combination with Zhao injection correction.

## Use one potential drop for accessibility and face velocity

Let $\phi_\infty$ be the potential at infinity or upstream and $\phi_f$ the face or interface potential. The 1-D
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

| $B$ | Accessibility and face velocity |
| ---: | --- |
| $B>0$ | Only $v_{n,\infty}\ge\sqrt B$ reaches the face and decelerates on the way |
| $B=0$ | Normal velocity is unchanged |
| $B<0$ | Every inflowing particle is accessible and accelerates toward the face |

The upstream VDF is truncated to accessible trajectories with $v_{\min}=\sqrt{\max(B,0)}$, and accepted velocities are mapped
to the face with the same $B$. Particle count and velocity at face arrival are therefore determined by one potential drop.

## Select the outer model that supplies the potential drop

| Configuration | Source of $\phi_f-\phi_\infty$ | Spatial outer orbit |
| --- | --- | --- |
| No correction | Zero | None |
| `reservoir_potential_model="infinity_barrier"` | Mean aperture potential and `sim.phi_infty` | Not solved |
| Split `linear_debye` | Interface drop derived from the surface zero mode | Analytic profile for return |
| Split `kinetic_1d` | Interface potential of the converged outer profile | Discrete profile for return |
| Zhao injection closure | Branch-specific local VDF and cutoff | Does not add $E(z)$ to the pusher |

Legacy `infinity_barrier` evaluates the batch-start field snapshot on an
`injection_face_phi_grid_n` by `injection_face_phi_grid_n` cell-centered grid in the aperture and uses the scalar mean
$\bar\phi_f$. It follows the same potential convention as the snapshot, including the point or triangle kernel, periodic field,
zero mode, outer state, and `sim.e0`, but it does not solve intermediate $E(z)$, turning position, flight time, or space charge.

In a split outer model, z-high `reservoir_face` species are interpreted as distributions at infinity. For `kinetic_1d`, the
interface field is a Poisson boundary condition, while the particle velocity changes through the solved
$\phi_I-\phi_\infty$. See [Outer plasma models](OuterPlasmaModels.en.html) for model selection and return through the same
profile.

See [Sheath injection closures](SheathInjectionClosures.en.html) for Zhao branches and `floating_no_photo` rules that override
source density, cutoff, and drift.

## Spread the initial position phase with jitter

Position is sampled uniformly in the aperture and moved $10^{-12}$ m inward along the face normal. If
`position_jitter_dt>0`, each particle receives a uniform virtual interval
$\tau\in[0,\texttt{position\_jitter\_dt})$ and $\mathbf x\leftarrow\mathbf x+\mathbf v\tau$. Axes periodic on both faces are
then wrapped into the primary box; other axes are clamped inside the box.

The jitter spreads the artificial phase in which all particles in a batch start exactly on one face. It changes only creation
position while preserving global simulation time and each particle's tracking horizon.

## Preserve global residuals across MPI and restart

With MPI, the root determines the global count and remainder before splitting particles across ranks. The global remainder is
stored in `macro_residuals.csv` and broadcast to all ranks on restart. Check species-resolved injected count and charge, the
resolved macro weight, residual, applied $v_{\min}$ and potential drop, and batch-averaged absorption and escape currents.

Changing `batch_duration` changes both expected particles per batch and the field-update interval. See
[Batch duration and stability](BatchDurationStability.en.html) for convergence of physical steady results.

## Code reference

- Flux integral, macro count, and Maxwell/grid sampling: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Barrier and outer-profile velocity correction: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Geometry and input-combination validation: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- Zhao species correction: [`bem_sheath_runtime.f90`](../src/physics/sheath/bem_sheath_runtime.f90)
- MPI-global macro residual: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
