title: Connecting to the outer sheath (Zhao stationary sheath)

Lang: [English](ZhaoStationaryClosure.en.md) | [日本語](ZhaoStationaryClosure.md)

# Connecting to the outer sheath (Zhao stationary sheath)

A BEACH cell is much shorter than the Debye length, so the outer sheath sees the top face (z-high) as the surface itself.
This model assumes a planar stationary sheath outside the cell, solves its zero-current root, and gives BEACH the currents,
barriers, and potential reference at the top face. It does not evolve the sheath in time.

By default the root is solved once at the start of the run. With outer-root refresh
(`sheath.coupling.outflow_refresh_batches`), the root is re-solved periodically from the photoelectrons observed to leave
through the top face, so that re-absorption of photoelectrons by the structure inside the cell is fed back to the outer sheath.
How to assemble a whole case is described in [Set up a periodic surface in plasma](PeriodicPlasmaSurface.en.html).

## What is solved

The model is a one-dimensional sheath formed by the solar wind upstream (potential 0 V) and photoelectrons emitted from the
surface. It finds the stationary solution with zero net current to the surface (the zero-current root).

$$
J_e + J_i + J_{escape} = 0
$$

$J_e$ is the absorbed electron current, $J_i$ the absorbed ion current, and $J_{escape}$ the current of photoelectrons that
escape; signs are contributions to surface charging. Without photoelectrons the condition is $J_e + J_i = 0$. Roots fall into
three types by the shape of the potential. Without photoelectrons only Type C exists.

| Type | Shape of the potential |
|---|---|
| Type A | Has a potential minimum $\phi_m$ between the surface and the upstream plasma |
| Type B | Decreases monotonically from the surface to upstream (positive surface) |
| Type C | Increases monotonically from the surface to upstream (negative surface) |

The root gives the type (branch), the wall potential $\phi_0$, the potential minimum $\phi_m$, the upstream electron density
$n_{e,\infty}$, and the current of each species. BEACH solves the root with the external library
[sheath-model](https://github.com/Nkzono99/sheath-model) (the commit is pinned in `fpm.toml`). sheath-model also checks that an
algebraic root satisfies the following conditions and rejects roots that do not.

- $E^2\ge0$ over the whole potential profile
- The profile approaches a neutral, field-free state upstream
- The ion flow is not blocked on its way

`sheath.zhao.branch="auto"` tries C → A → B when the solar elevation is below 20 degrees and A → B → C otherwise, and returns
the first root that holds. If no root holds, the run stops at startup with the reason.

## Assumptions

- The sheath is planar, collisionless, unmagnetized, and stationary.
- Solar-wind electrons are a Maxwellian drifting at the same normal solar-wind velocity as the ions. Slow electrons that turn
  back upstream are reflected as a mirror image of the incoming distribution.
- Ions are a cold beam, and the ion temperature is not used for the root.
- Photoelectrons are emitted from the surface as a half-Maxwellian with density $n_{pe,0}=s_{UV}\,n_{pe,ref}\sin\alpha$
  ($\alpha$ is the solar elevation, $n_{pe,ref}$ the reference density, and $s_{UV}$ a scale factor).
- The cell is a wall of zero thickness for the outer sheath.

### Electron drift and the upstream band

Give the electrons the same normal solar-wind drift as the ions and set `sheath.zhao.upstream_band_tolerance` to 0.1. The
examples use this setting.

**The upstream band:** in Type A and Type C with inward-drifting electrons and reflected slow electrons, the upstream
distribution of reflected electrons is the mirror image of the drifting incoming distribution and has a kink at zero velocity.
Because of this kink, a term $u\,h\log(1/h)$ appears in the electron density close to upstream; as $h\to0$ it exceeds the other
terms and makes $E^2<0$ ($u$ is the ratio of the drift to the electron thermal speed and $h$ the depth below the upstream
potential divided by $T_e$). No solution connects exactly to the neutral, field-free upstream state. $E^2<0$ occurs only in a
narrow potential band next to the upstream potential. `upstream_band_tolerance` accepts roots whose band is no wider than that
fraction of $\lvert\phi_m\rvert$ (A) or $\lvert\phi_0\rvert$ (C) and records the width as
`surface_current_model_upstream_negative_band_V` in `summary.txt`. Roots with a wider band are not used.

**Physical meaning of the band:** the band is where the idealization that slow upstream electrons are an exact mirror image of
the drifting incoming distribution fails. That distribution has a cusp-shaped dip at $v_z=0$; the integral of the Penrose
stability criterion diverges logarithmically, so it is kinetically unstable. In reality, relaxation by the instability and
scattering by surface roughness are expected to smooth the slow-electron distribution.

**Error estimate:** removing the drift only from incoming electrons slower than $a_s v_{th}$, so that the distribution is smooth,
and solving again with orbit-consistent densities removes the band and gives $E^2\ge0$ everywhere. For $n=5$ cm⁻³, $T_e=10$ eV,
a 400 km/s solar wind at normal incidence, and photoelectrons of 4.5 µA/m² at 2.2 eV, the differences from the roots with an
accepted band were as follows ([validation record](https://github.com/Nkzono99/sheath-model/blob/main/outputs/upstream_relaxation_20261009/REPORT.md)).

| Type | Band | $a_s$ that removes the band | $\phi_0$ difference for $a_s=0.07$–$0.2$ | $\phi_m$ difference | Absorbed electron flux difference |
|---|---:|---:|---:|---:|---:|
| A | 52 mV | 0.07 | −23 to −81 mV | −24 to −68 mV | below 0.3% |
| C (no photoelectrons) | 2.9 mV | 0.02 | −1 to −81 mV ($a_s=0.02$–$0.2$) | — | 0 |

The spread over $a_s$ measures the uncertainty from the treatment of slow electrons. It is much smaller than the effect of
neglecting the drift (about 1 V for A and 2.7 V for C, [Known limitations](#known-limitations)), so BEACH treats the drift by
accepting the band. A tolerance of 0.1 accepts both A (band 6.6% of $\lvert\phi_m\rvert$) and C (0.04%) above.

The default 0 is the exact condition, which rejects drifting A/C (for example, at a solar elevation of 60 degrees B is selected,
and at 10 degrees no root holds). Then set the electron drift to zero. The ion drift is required in both cases.

## How BEACH uses the root

### Fix the currents to the root

The three role species (electrons, ions, photoelectrons) use the fixed-current charge closure (`fixed_current`).
In each batch, the tracked distribution is scaled uniformly to the target charge, which is the root current multiplied by the
batch duration.

| Channel | Sign | Treatment |
|---|---:|---|
| Electron absorption $J_e$ | negative | Target for absorption on the surface |
| Ion absorption $J_i$ | positive | Target for absorption on the surface |
| Photoelectron emission $J_{emit}$ | positive | Target for the reaction charge left on the emitter |
| Photoelectron return $J_{return}$ | ≤ 0 | Target for re-absorption on the surface |
| Photoelectron escape $J_{escape}$ | positive | Current that leaves the box and is not deposited |

The root satisfies $J_{return}=J_{escape}-J_{emit}$ and $J_e+J_i+J_{escape}=0$, so the surface currents close as
$J_e+J_i+J_{return}+J_{emit}=0$. Emission and return are scaled separately, and the net photoelectron current, a small
difference of two large currents, is never used as the denominator of a scale. The total surface charge is therefore
constrained by the floating condition, and only its distribution over the elements is set by tracking.

### Match the top-face potential to the wall potential

With a backend that treats the zero mode of the field separately (`fields.periodic.backend` set to `cached_kneq0` or
`panel_spectral_reference`), the horizontal mean potential of the top face is matched to the wall potential $\phi_0$.
Adding a constant in vacuum does not change potential differences in the cell or particle orbits; it makes the potentials
values referenced to the upstream plasma at 0 V. With a setting that does not separate the zero mode, the potential cannot be
matched to $\phi_0$ and the run warns at startup.

### Particles entering and leaving through the top face

**Inflow:** from the upstream distribution at 0 V, only particles that can pass the outer-sheath barrier and reach the top face
are selected, and their normal velocity is mapped to the top face by energy conservation.

$$
\frac12 m v_{n,f}^{2}=\frac12 m v_{n,\infty}^{2}-q\phi_f
$$

$\phi_f$ is the mean top-face potential at the start of the batch and $v_{n,\infty}$ the upstream normal velocity.

**Outflow:** electrons and photoelectrons that cross the top face outward are tested against the outer barrier with the
potential $\phi_c$ at the crossing point and the normal kinetic energy. Particles that cannot pass are returned with their normal
velocity reversed; particles that can pass escape.

| Type | Barrier for incoming electrons | Barrier for outgoing electrons and photoelectrons | Ions |
|---|---:|---:|---:|
| Type A | $\phi_m$ | $\phi_m$ | 0 V |
| Type B / C | 0 V | 0 V | 0 V |

For example, with Type A, $\phi_0\approx6.9$ V, and $\phi_m\approx-0.07$ V, an electron leaving the top face is returned if its
normal energy is below about 7 eV. In Type C the wall is negative, so every outgoing electron escapes.

In an x/y-periodic cell, the return position is chosen again uniformly on the top face. The lateral distance a particle travels
in the outer sheath, the flight time times the tangential velocity, is of order meters and much larger than the cell width.
The normal kinetic energy is corrected with the potential $\phi_r$ at the return position so that the total energy is conserved.

$$
\frac12 m v_{n,r}^{2}=\frac12 m v_{n,c}^{2}+q(\phi_c-\phi_r)
$$

A slow particle for which the right-hand side is zero or negative stays only briefly outside and hardly moves sideways, so it is
returned at the crossing point. The tangential velocity is unchanged.

### Inject electrons at the upstream density

Electrons that reach the wall are absorbed and do not come back, so the upstream electron distribution lacks its fast outward
part. The root solves the density $n_{e,\infty}$ of the electron Maxwellian so that the plasma is quasi-neutral upstream
including this missing part. For Type C without photoelectrons,

$$
n_{e,\infty}=\frac{n_i}{1-\frac12\operatorname{erfc}\sqrt{-\phi_0/T_e}}
$$

which is larger than the ion density $n_i$. BEACH injects electrons through the top face at this $n_{e,\infty}$ and ions at the
configured $n_i$. The injected electron flux then matches the electron current of the root, and the electron scale is 1 within
statistical error. $n_{e,\infty}$ is updated each time the outer root is refreshed.

### Refresh the outer root

Structure inside the cell changes the emission seen by the outer sheath. Part of the photoelectrons emitted by the surface is
re-absorbed in the cell, so fewer photoelectrons leave through the top face than are emitted, with a different energy
distribution. Outer-root refresh feeds this change back to the outer sheath. Every `outflow_refresh_batches` accepted batches:

1. Sum the number $N_{out}$ and the normal kinetic energy of photoelectrons crossing the top face outward, and the number
   $N_{emit}$ emitted by the surface. Compute the transmission $\eta=N_{out}/N_{emit}$ and the mean normal energy $\bar K$.
2. Replace the photoelectron source of the outer sheath with a half-Maxwellian of flux $\eta\Gamma_{emit}$ and temperature
   $\bar K$, and solve the zero-current root. $\Gamma_{emit}$ is the surface emission flux set by the configuration; the surface
   emission target does not change.
3. Use the currents, barriers, wall potential, and upstream electron density of the new root from the next batch.

| Situation | Behavior |
|---|---|
| `sheath.zhao.branch` is explicit | Only that type is solved |
| `sheath.zhao.branch="auto"` | The previous type is tried first, then the usual order. The previous root is the initial guess |
| No root is found, or the outflow is zero | Warn and keep the previous root |
| Restart | Rebuild the root from the saved photoelectron source. A partial sum is discarded and counting restarts at the restart point |

Choose the interval so that the photoelectron samples within it fix the mean energy to a few percent and the interval is much
shorter than the change of surface charging. If only the steady state matters, the state reached does not depend on the interval.

## Configuration

```toml
[sheath]
closure = "zero_current"

[sheath.zhao]
branch = "auto"
upstream_band_tolerance = 0.1

[sheath.species]
electron = "solar_wind_electron"
ion = "solar_wind_ion"
photoelectron = "photoelectron"

[sheath.photoelectrons]
solar_elevation_deg = 60.0
ref_density_m3 = 6.4e7

[sheath.coupling]
outflow_refresh_batches = 50
```

- Role species use `charging.closure="fixed_current"`. Electrons and ions enter by boundary inflow through the top face;
  photoelectrons use `photo_raycast` from the top face with reaction charge. The top face is open.
- Write the same solar-wind density for electrons and ions (a mismatch is a configuration error). Set both the electron and
  the ion drift to the normal solar-wind velocity (negative z). With zero electron drift no tolerance is needed
  ([Electron drift and the upstream band](#electron-drift-and-the-upstream-band)).
- Without photoelectrons (Type C), set `sheath.photoelectrons.source_scale=0.0` and omit the photoelectron species and its keys.
  Outer-root refresh is not available.

All keys and combination constraints are in [Input parameters](Parameters.en.html#surface_current_model-external-sheath-closure). Complete examples are
[`examples/grouped/zero_current.toml`](../examples/grouped/zero_current.toml),
[`zero_current_refresh.toml`](../examples/grouped/zero_current_refresh.toml), and
[`zero_current_no_photo.toml`](../examples/grouped/zero_current_no_photo.toml).

## Outputs to check

| Output | What to check |
|---|---|
| `summary.txt` | `surface_current_model_zhao_branch`, `surface_current_model_phi0_V`, `surface_current_model_phi_m_V` (initial root). The two residuals `surface_current_model_pe_budget_residual_current_density_A_m2` and `surface_current_model_surface_budget_residual_current_density_A_m2` are near zero |
| `fixed_current_history.csv` | `target_over_tracked`. If the outer root and tracking in the cell agree, every species stays near 1 |
| `charge_ledger.csv` | Tracked, target, and applied charge per species, `fixed_*_weight_scale`, and counts |
| `matching_plane_history.csv` (outer-root refresh) | `phi_H_V` (the current $\phi_0$), the photoelectron outflow flux and mean normal energy, and `residual` |
| `top_reference_history.csv` | The mean top-face potential is close to $\phi_0$ (see the known limitations below) |

If the scales differ strongly between species, the balance of currents on each element departs from the state set by
tracking alone. Column definitions are in [Output formats](OutputReference.en.html#zhao_stationary).

A small zero-current residual does not mean that the tracked spatial distribution of absorption and return has converged.
Vary the ray count, batch duration, and random seed and confirm that the per-element distribution does not change.

## Scope

- The field, space charge, Debye shielding, distance to turning points, and flight times outside the box are not solved.
  Reflection at the top face is an approximation that shortens the round trip outside to a reversal at the boundary.
- The uniform magnetic field must be zero.
- Because the solution is planar and stationary, it does not follow curvature, collisions, transients of the outer sheath, or
  illumination and plasma conditions that change during the run. Outer-root refresh follows only the change of photoelectrons
  leaving through the top face.

### Known limitations

- **Accepting the upstream band:** a root accepted with the tolerance does not follow the idealized electron distribution inside
  the upstream band. The effect on the wall quantities is tens of mV ([Electron drift and the upstream band](#electron-drift-and-the-upstream-band)).
- **Error of a zero electron drift:** at the same potential, the drift increases the electron flux at first order in the drift
  (the one-way flux of a drifting Maxwellian through a plane is $e^{-u^2}+\sqrt\pi\,u\,(1+\operatorname{erf}u)$ times the flux
  without drift). At a zero-current root the absorbed electron flux is constrained by the ions and the photoelectron escape, so
  the difference appears mainly in the potentials. Under the same conditions as above ($u=0.21$, flux factor 1.42), the
  zero-drift roots differed from the drifting roots as follows.

  | Type | $\phi_0$ (drift 400 km/s → 0) | $\phi_m$ | Absorbed electron flux |
  |---|---|---|---|
  | A | 6.04 → 7.05 V | −0.79 → −0.13 V | −6% |
  | C (no photoelectrons) | −7.34 → −4.63 V | — | Unchanged (equal to the ion flux) |

  In a three-layer particle bed (4000 s), zero drift raised all potentials in the cell by about 1 V (A) and 3–4 V (C), while the
  potential differences between layers changed by 1% or less.
- **Reduced photoelectron source:** outer-root refresh replaces the photoelectrons leaving through the top face with a
  half-Maxwellian defined by only two quantities, flux and mean energy.
- **Escaping electrons are not fed back:** solar-wind electrons that are turned back in the cell and leave through the top face are
  not included in the outer root, which treats every electron reaching the wall as absorbed. The fixed-current closure deposits
  this part as absorption, so the balance of total currents is kept. The escaping fraction can be read from the departure of the
  electron scale from 1 in `fixed_current_history.csv`.
- **Offset of the mean top-face potential:** as the surface charge grows, the mean top-face potential departs from $\phi_0$. The
  far correction (`cached_kneq0`) does not remove the zero mode completely; the offset is proportional to the charge and creates
  a spurious vertical field above the particles. Check it by comparing `potential_mean_V` in `top_reference_history.csv` with
  $\phi_0$ ([accuracy of the far correction](PeriodicFarCorrection.en.html)).
