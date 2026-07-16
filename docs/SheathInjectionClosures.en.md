title: Sheath injection closures

Lang: [日本語](SheathInjectionClosures.md) | [English](SheathInjectionClosures.en.md)

# Sheath injection closures

`sim.sheath_injection_model` derives density, drift, and cutoff values for `reservoir_face` and `photo_raycast` from an analytic
sheath closure. These corrections apply to source sampling; generated particles advance in the separately composed field
snapshot.

| Model | Quantity solved | Values applied to particle sources |
| --- | --- | --- |
| `none` | None | Use configured VDFs unchanged |
| `floating_no_photo` | Negative floating potential without photoelectrons | Electron normal cutoff |
| `zhao_a/b/c` | Analytic closure on the selected Zhao branch | Electron/ion/photoelectron density, drift, cutoff, and emission current |
| `zhao_auto` | Zhao branch search ordered by solar elevation | Same corrections from the converged branch |

See [Kinetic 1-D outer plasma](KineticOuterPlasma.en.html) for a self-consistent outer Poisson profile and
[Unified linear response](UnifiedLinearResponse.en.html) for a linear field including rough surfaces.

## Select species by source role

Enabled species are scanned in configuration order:

- first negative `reservoir_face`: solar-wind or ambient electron;
- first positive `reservoir_face`: ion;
- for Zhao models, first negative `photo_raycast`: photoelectron.

Electron and ion reservoirs must share `inject_face`. Inward drift comes from each source `drift_velocity` or is resolved by
`sheath_electron_drift_mode` and `sheath_ion_drift_mode`; both must point inward. Zhao requires positive electron and
photoelectron temperatures, ion density, reference photoelectron density, ion drift, and particle masses.

## Balance electron and ion inflow with `floating_no_photo`

Ion inflow $\Gamma_i$ is computed from the ordinary flux-weighted drifting Maxwell distribution. For trial negative potential
$\phi_0\le0$, electron cutoff is

$$
v_{e,\min}(\phi_0)
=\sqrt{\frac{2e\max(0,-\phi_0)}{m_e}}.
$$

Bisection solves

$$
F(\phi_0)=\Gamma_e(v_n\ge v_{e,\min})-\Gamma_i=0.
$$

The initial lower bound is minus 128 times electron temperature in eV, expanded until a bracket is found, with at most 80
iterations. A nonpositive ion flux or failure to bracket a negative root stops.

The resulting $\phi_0$ becomes $v_{\min}$ for the first electron reservoir. Ion VDF is unchanged and no photoelectron is included.
This is a reduced current-balance closure with no spatial $\phi(z)$, $E(z)$, turning point, or flight time.

## Form the dimensionless Zhao variables

For solar elevation $\alpha$ and reference photoelectron density $n_{\mathrm{phe,ref}}$, surface source density is

$$
n_{\mathrm{phe},0}=n_{\mathrm{phe,ref}}\sin\alpha.
$$

The implementation also forms

$$
v_{\mathrm{swe,th}}=\sqrt{\frac{2eT_\mathrm{swe}}{m_e}},
\qquad
v_{\mathrm{phe,th}}=\sqrt{\frac{2eT_\mathrm{phe}}{m_e}},
$$

$$
c_s=\sqrt{\frac{eT_\mathrm{swe}}{m_i}},
\quad
M=\frac{u_i}{c_s},
\quad
u=\frac{u_e}{v_\mathrm{swe,th}},
\quad
\tau=\frac{T_\mathrm{swe}}{T_\mathrm{phe}}.
$$

A nonlinear root of analytic density and current conditions determines surface potential $\phi_0$, potential minimum $\phi_m$
where applicable, and effective solar-wind electron density $n_{\mathrm{swe},\infty}$.

## Select a Zhao branch from solar elevation and roots

| Branch | Potential and population structure |
| --- | --- |
| A | Nonmonotone $\phi_0>0$, $\phi_m<0$; captured photoelectrons below the minimum and reflected solar-wind electrons above it |
| B | Monotone $\phi_0>0$; photoelectron capture without reflected solar-wind electrons |
| C | Monotone $\phi_0<0$; reflected solar-wind electrons and no photoelectron cutoff |

`zhao_auto` tries C→A→B for $\alpha<20^\circ$ and A→B→C otherwise. This is branch selection within the Zhao family. Failure of
all branches stops without switching to `floating_no_photo` or an outer Poisson model. Explicit `zhao_a/b/c` solves only the
requested branch and does not move to another branch on failure.

## Apply closure results to each source

### Ambient electron

Effective density becomes $n_{\mathrm{swe},\infty}$, with branch-dependent cutoff

$$
v_{e,\min}=\sqrt{\frac{2e\max(0,\Delta\phi_e)}{m_e}},
\qquad
\Delta\phi_e=
\begin{cases}
-\phi_m & A,\\
0 & B,\\
-\phi_0 & C.
\end{cases}
$$

### Photoelectron

Emission normal cutoff is

$$
v_{pe,\min}=\sqrt{\frac{2e\max(0,\Delta\phi_{pe})}{m_{pe}}},
\qquad
\Delta\phi_{pe}=
\begin{cases}
\phi_0-\phi_m & A,\\
\phi_0 & B,\\
0 & C.
\end{cases}
$$

Normal drift is overridden to zero, and free-photoelectron current density becomes

$$
J_{pe}=\frac{|q_{pe}|n_{\mathrm{phe},0}v_{\mathrm{phe,th}}}{2\sqrt\pi}
\begin{cases}
\exp[(\phi_m-\phi_0)/T_\mathrm{phe}] & A,\\
\exp[-\phi_0/T_\mathrm{phe}] & B,\\
1 & C.
\end{cases}
$$

$T$ and $\phi$ use the corresponding eV/V convention in these exponentials. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for converting current density into ray macro weights.

### Ion and local reference plane

With `sheath_reference_coordinate`, that plane is sheath coordinate $z_s=0$. The Zhao 1-D profile is sampled at the distance to
the common reservoir face, reconstructing free/reflected electrons, free/captured photoelectrons, and ion speed and density. These
values override electron density/cutoff and ion cold-beam normal speed.

Because this path has already built a local VDF, it does not also apply a generic `reservoir_potential_model` barrier shift.

## Scope of the Zhao closure

The Zhao result is a literature-based source-VDF precorrection. The Boris field snapshot and Zhao root are not updated from
batch-varying `q_elem`, so this closure does not represent a self-consistent external field for arbitrary 3-D surface geometry.

Use `kinetic_1d + kinetic_1d_profile_return` or a validated unified configuration when incoming and outgoing particles share one
self-consistent potential profile.

## Apply only one source correction

- `sheath_injection_model` is rejected with `reservoir_potential_model`.
- A `velocity_distribution="grid"` reservoir cannot currently use Zhao or floating correction.
- `kinetic_1d_profile_return` rejects Zhao and `reservoir_potential_model`.
- Zhao requires negative electron, positive ion, and negative photoelectron species.
- `floating_no_photo` uses only negative electron and positive ion species and does not modify photoemission.

Photoelectrons are tracked as ordinary particles after emission. Finite boxes use
`open_boundary_model="potential_barrier"`; external sheaths use `outer_plasma.return_model` and
`coupling.particle_transfer_mode` to decide return or escape. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for charge balance.

## Code reference

- Floating root and Zhao density/current/root solve: [`bem_sheath_model_core.f90`](../src/physics/sheath/bem_sheath_model_core.f90)
- Species detection and runtime overrides: [`bem_sheath_runtime.f90`](../src/physics/sheath/bem_sheath_runtime.f90)
- Apply corrections to source sampling: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Input and combination validation: [`bem_app_config_parser.f90`](../src/config/app_config_parser/bem_app_config_parser.f90)
