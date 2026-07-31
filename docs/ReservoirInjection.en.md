title: Reservoir inflow through simulation boundaries

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# Reservoir inflow through simulation boundaries

`[particles.species.boundary_inflow]` defines particles crossing a simulation boundary from a plasma reservoir outside the
box. The upstream VDF, box-face area, `batch_duration`, and macro-particle weight determine the per-face injection count,
initial positions, face-arrival velocities, and remainder carried to the next batch.

Boundary inflow has a separate responsibility from `source_mode`. In the initial implementation, it can share a species only
with `source_mode="volume_seed"`. Use `source_mode="plane_source"` for an explicit internal emission plane, but combining it
with boundary inflow for the same species fails closed. See [Particle sources and boundary inflow](ParticleSourcesBoundaries.en.html) to choose a source and
[particle escape and local return](ParticleEscapeReturn.en.html) for outward open-face processing.

## Select inflow faces

Set a face to `"reservoir"` to inject through the complete box face.

```toml
[[particles.species]]
species_key = "electron"
source_mode = "volume_seed"
npcls_per_step = 100
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
w_particle = 1.0e5
drift_velocity = [0.0, 0.0, -4.0e5]

[particles.species.boundary_inflow]
z_high = "reservoir"
```

The available keys are `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, and `z_high`; more than one can be enabled.
Reservoir inflow cannot be assigned to faces on `domain.periodic_axes`. Positions are uniform over each selected box face,
so `pos_low`, `pos_high`, and `inject_face` do not define an aperture for boundary inflow.

Outward particle treatment is configured separately. Use `[particle_boundary]` for global actions and
`[particles.species.boundary]` for species overrides. The effective outward action on a reservoir inflow face must be
`open`; `boundary_inflow` does not override it.

## Weight a Maxwell distribution by inflow flux

Let $\mathbf n$ be the selected face's inward unit normal and $A$ its area. For an upstream drifting Maxwell distribution,
define normal drift $u_n=\mathbf u\cdot\mathbf n$ and thermal standard deviation
$\sigma=\sqrt{k_\mathrm{B}T/m}$. With upstream normal-speed floor $v_{\min}$, the one-way inflow flux is

$$
\Gamma_\mathrm{in}
=n\int_{v_n\ge v_{\min}}v_n f_n(v_n)\,dv_n
=n\left[u_n\{1-\Phi(a)\}+\sigma\varphi(a)\right],
\qquad
a=\frac{v_{\min}-u_n}{\sigma},
$$

where $\Phi$ and $\varphi$ are the standard-normal CDF and density. At zero temperature,
$\Gamma_\mathrm{in}=nu_n$ when $u_n\ge v_{\min}$ and zero otherwise.

Tangential velocities are Gaussian samples. Normal velocity follows the half-range distribution proportional to
$v_n f_n(v_n)$ because surface-crossing probability is proportional to normal velocity; it is not a volume Maxwell
distribution placed directly on the boundary.

## Convert velocity-grid values into flux

With `velocity_distribution="grid"`, velocities are sampled from CSV points and nonnegative values $f$.

| `velocity_grid_pdf_kind` | Meaning of CSV $f$ | Sampling weight |
| --- | --- | --- |
| `phase_space` | Phase-space distribution | $\max(v_n,0)f$ |
| `flux_weighted` | Already weighted for inflowing particles | $f$ |

Both retain only entries with $v_n\ge v_{\min}$ and $v_n>0$. `velocity_grid_sampling="rectilinear"` interpolates a complete
rectilinear grid, `discrete` samples input points directly, and `auto` uses `rectilinear` when possible.

For a velocity-grid distribution, count comes from `particle_flux_m2_s` or `current_density_a_m2`; density and temperature do
not set the count.

## Carry one remainder per face to preserve long-time flux

Let $\Gamma_\mathrm{in}$ be either the flux derived from a Maxwell distribution or the flux supplied for a velocity grid.
Expected physical and macro-particle inflows are

$$
N_\mathrm{phys}=\Gamma_\mathrm{in}A\,\Delta t_\mathrm{batch},
\qquad
N_\mathrm{macro,expected}=\frac{N_\mathrm{phys}}{w}.
$$

If $r$ is the fractional remainder for one species-face pair,

$$
N_\mathrm{macro}=\left\lfloor r+N_\mathrm{macro,expected}\right\rfloor,
\qquad
r\leftarrow r+N_\mathrm{macro,expected}-N_\mathrm{macro}.
$$

Individual batch counts can vary, while the long-time total approaches the expected flux. Supply a positive `w_particle` or
derive it from `target_macro_particles_per_batch`. The latter controls the total sample count across all enabled inflow
faces, not physical flux.

## Select inflow correction as a boundary condition

`[reservoir]` supplies the external-plasma condition shared by boundary inflow and the scalar barrier for particles leaving
through open faces.
For an arbitrary internal plane, neither which side represents the external plasma nor how its plane potential maps to the
upstream reference is unique. At a simulation boundary, the reservoir is naturally defined outside the box, so
`infinity_barrier` can map `phi_infty` to the boundary-face potential as an outside-to-boundary correction. For this reason,
new configurations apply inflow-side potential/barrier correction to `boundary_inflow`, not to `plane_source`.
Deprecated `reservoir_face` retains its previous correction as compatibility behavior for existing cases.

```toml
[particle_boundary]
z_high = "open"
ordinary_open_model = "potential_barrier"

[reservoir]
inflow_model = "infinity_barrier"
phi_infty = 0.0
face_potential_grid_n = 5
```

| Configuration | Effect on boundary inflow |
| --- | --- |
| `reservoir.inflow_model="source_vdf"` | Interpret the configured VDF at the boundary and apply no potential correction |
| `reservoir.inflow_model="infinity_barrier"` | Correct accessible flux and normal velocity using the mean face potential and `reservoir.phi_infty` |

Let $\phi_\infty$ be the potential at infinity or upstream and $\phi_f$ the boundary-face potential. The 1-D electrostatic
map preserves tangential velocity and normal energy:

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}.
$$

| $B$ | Upstream accessibility and boundary-face velocity |
| ---: | --- |
| $B>0$ | Only $v_{n,\infty}\ge\sqrt B$ reaches the face and decelerates on the way |
| $B=0$ | Normal velocity is unchanged |
| $B<0$ | Every inflowing particle is accessible and accelerates toward the face |

`infinity_barrier` evaluates every inflow face on an $N\times N$ cell-centered grid set by
`reservoir.face_potential_grid_n` and uses its batch-start mean $\bar\phi_f$. It follows the field snapshot's potential
convention, including the fixed P0 triangle kernel, periodic field, zero mode, and `sim.e0`.

`particle_boundary.ordinary_open_model="potential_barrier"` uses the same `phi_infty` and the local potential at the
particle's actual crossing point to decide outward return or escape. Thus, inflow `infinity_barrier` uses a face mean,
whereas outward `potential_barrier` uses the event position.

## Model scope

Defining this model on a boundary gives density, temperature, VDF, and `phi_infty` a consistent interpretation as conditions
outside the box. BEACH nevertheless solves only a scalar velocity map to the boundary and particle motion after crossing it.
It does not solve outside-box trajectories, intermediate $E(z)$, turning positions, flight times, space charge, or a
self-consistent external sheath. This setting does not restore the removed `[outer_plasma]` or `[coupling]` models.

A uniform external electric field has no finite potential at infinity. When it is combined with `infinity_barrier`, define
`phi_infty` as an effective reservoir reference and validate the physical interpretation separately.

## Initial position, MPI, and restart

To avoid artificial alignment on one face-time plane, BEACH assigns each particle a uniform virtual interval
$\tau\in[0,\texttt{sim.dt})$ and shifts only its initial position before tracking:
$\mathbf x\leftarrow\mathbf x+\mathbf v\tau$. Global simulation time and each particle's tracking horizon are unchanged.

With MPI, the root rank determines global counts and remainders before splitting particles across ranks. Remainders are stored
in checkpoints and restored on resume. When checking results, distinguish species- and face-resolved injected count and charge,
resolved macro-particle weight, remainder, applied $v_{\min}$ and potential drop, absorbed current, and escape current.

Checkpoint schema v6 writes `macro_residuals.csv` as `species_idx,face,residual`: `face=0` denotes the legacy source and
`1..6` denote boundary faces. The older two-column form remains readable.

Changing `batch_duration` changes both expected particles per batch and the field-update interval. See
[Batch duration and stability](BatchDurationStability.en.html) for convergence of steady results.

## Migrate from legacy `reservoir_face`

`source_mode="reservoir_face"` remains accepted for compatibility, but new configurations should use
`[particles.species.boundary_inflow]`. The legacy form has an aperture selected by `inject_face` and
`pos_low` / `pos_high`; when moving to a complete-face boundary inflow, verify that the physical inflow does not change.
Move explicit internal rectangular emission surfaces to `source_mode="plane_source"`. BEACH does not silently convert the
legacy form to either setting.

## Code reference

- Flux integration, macro-particle count, and Maxwell or velocity-grid sampling: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Boundary inflow and scalar-barrier velocity correction: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Input-combination validation: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- MPI-global macro-particle remainders: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
