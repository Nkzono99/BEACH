title: Electron drift and the upstream band

Lang: [English](ZhaoUpstreamBand.en.md) | [日本語](ZhaoUpstreamBand.md)

# Electron drift and the upstream band

When solar-wind electrons are given the same drift as the ions, the Type A / C roots of the outer sheath are rejected by
default. This page explains why, and what to write in `sheath.zhao.upstream_field_floor_v_m`. The recommended value is the
solar-wind motional field $\lvert\mathbf{u}\times\mathbf{B}\rvert$, about 2 mV/m for 400 km/s and 5 nT. After reading you know
where an accepted root is an approximation and how much it affects the wall-side quantities. The outer-sheath model as a whole
is described in [Connecting to the outer sheath](ZhaoStationaryClosure.en.html).

## Why the band appears

Write the upstream distribution of electrons moving toward the surface as $g(a)$, with $a$ the inward normal speed in the
surface frame. Slow electrons (those that do not reach $\phi_m$ in Type A or the wall in Type C) are reflected by the potential
and return upstream with the same speed. Mirror reflection and Liouville's theorem give, at a depth $h=-\phi/T_e$ below the
upstream potential $0$, the electron density

$$
n_e(h)-n_e(0)=\int_0^\infty\bigl[g(\sqrt{w^2+h})-g(w)\bigr]\,dw
\approx g(0)\,O(h)+g'(0)\,\tfrac{h}{2}\log\tfrac{1}{h}
$$

The coefficient $g'(0)$ is the slope of the incoming distribution at zero speed in the surface frame. A Maxwellian without drift
has $g'(0)=0$ and no such term, but a distribution with drift $u$ has $g'(0)=2u\,g(0)>0$, so an electron excess grows as
$h\log(1/h)$ toward upstream. The ion and photoelectron densities are regular in $\phi$ ($O(h)$) and cannot cancel it. Including
the ion temperature does not change this.

The band requires all three of the following assumptions at once.

1. $g'(0)\ne0$ in the surface frame (an incoming distribution with drift)
2. Slow electrons are mirror-reflected exactly by a stationary one-dimensional potential
3. The profile connects exactly to a neutral, field-free state at upstream infinity

## What happens to the potential profile

sheath-model does not integrate the profile in $z$; it treats the first integral of Poisson's equation
$E^2(\phi)=(2/\varepsilon_0)\int_0^\phi\rho\,d\phi'$ as a function of $\phi$. The root conditions (neutrality at 0 V upstream, zero
current, and the Type A connection $E^2(0)=0$) are all integral conditions in $\phi$ space, so they remain solvable when $E^2$
becomes negative on the way.

For an accepted root, $E^2(\phi)$ crosses zero at $\phi_b$ just before upstream (about −42 mV for Type A in the example),
is negative beyond it, and returns to zero at $\phi=0$. There is no real $z(\phi)$ where $E^2<0$, so the profile ends at $\phi_b$
with a maximum of zero slope. $\phi_b$ lies at a finite distance of about $3\lambda_D$ (about 30 m) from the top face, and the
wall-side profile $z<z_b$ is built normally. The profile generation in sheath-model also stops before the band. The quantities
BEACH uses ($\phi_0$, $\phi_m$, $n_{e,\infty}$, the currents, and the barriers) come from the algebraic root of the integral
conditions, so the existence of the band itself does not affect them.

The field missing inside the band, $\sqrt{-\min E^2}$, is the electric field that has to be supplied from outside for the
one-dimensional electrostatic sheath to connect upstream. It is the quantity of the acceptance criterion below.

## Acceptance criteria

| Key | Default | Accepted roots |
|---|---:|---|
| `sheath.zhao.upstream_field_floor_v_m` | 0 | The field deficit $\sqrt{-\min E^2}$ inside the band is at most this value [V/m] |
| `sheath.zhao.upstream_band_tolerance` | 0 | The potential width of the band is at most this fraction of $\lvert\phi_m\rvert$ (A) or $\lvert\phi_0\rvert$ (C) |

When both are given, a root must satisfy both; with both 0, drifting A/C are rejected. Either criterion accepts the same
algebraic root, so the wall-side quantities do not depend on the choice.

**Field floor (recommended).** This replaces the upstream boundary "infinity with $E=0$" by "the point where the field has
fallen to the external level". The solar wind carries the motional field $\lvert\mathbf{u}\times\mathbf{B}\rvert$ (about 2 mV/m
for 400 km/s and 5 nT), and its magnetic field and density fluctuations act at the same level, so the one-dimensional
electrostatic sheath cannot claim fields below it. The measured deficits are all below the motional field.

| Conditions | Type | Band width | Field deficit |
|---|---|---:|---:|
| Example ($n=8.7$ cm⁻³, $T_e=12$ eV, 405 km/s, solar elevation 60 degrees, photoelectrons) | A | 42 mV | 0.63 mV/m |
| Same, without photoelectrons | C | 1.8 mV | 0.03 mV/m |
| $n=5$ cm⁻³, $T_e=10$ eV, 400 km/s, photoelectrons 4.5 µA/m² at 2.2 eV | A | 52 mV | 0.67 mV/m |

**Width tolerance.** The same band judged by its potential width. Under the conditions above the A band is 4–7% of
$\lvert\phi_m\rvert$ and the C band 0.04%.

**Zero electron drift.** The band disappears, but the electron flux changes by $e^{-u^2}+\sqrt\pi\,u\,(1+\operatorname{erf}u)$
(1.42 for $u=0.21$), and at a zero-current root the difference appears in the potentials. Under the conditions of the third row,
the zero-drift root changed $\phi_0$ from 6.04 to 7.05 V (A) and from −7.34 to −4.63 V (C). In a three-layer particle bed all
potentials in the cell rose by about 1 V (A) and 3–4 V (C), while the potential differences between layers changed by 1% or
less. This is much larger than the uncertainty of an accepted band (tens of mV), so BEACH treats the drift by accepting the
band.

## What fills the band in a real plasma

The band concerns electrons with $\lvert v_z\rvert\lesssim0.07\,v_{th}$ (about 150 km/s) in the surface frame, below 1 eV. The
physics left out of the one-dimensional electrostatic model dominates for exactly this population.

| Quantity ($n=8.7$ cm⁻³, $T_e=12$ eV, $B=5$ nT) | Value |
|---|---:|
| Debye length $\lambda_D$ | about 9 m |
| Band region ($\lvert\phi\rvert<50$ mV, $3$–$7\lambda_D$) | 30–60 m |
| Electron cyclotron period | 7 ms |
| Transit time $\tau$ of slow electrons through the band region | 0.2–0.4 ms |
| Gyration angle $\Omega\tau$ during the transit | 0.2–0.35 rad |

- **Gyration mixing by the magnetic field.** Thermal electrons cross the band in 20–30 µs, so the unmagnetized approximation
  holds for them, but slow electrons gyrate by 0.2–0.35 rad during the crossing. With a field component perpendicular to the
  normal, part of $v_\perp\approx v_{th}$ turns into $v_z$ and changes it by 400–700 km/s, more than the width of the cusp
  (150 km/s), so the distribution near $v_z\approx0$ is smoothed.
- **Solar-wind fluctuations.** Density fluctuations of a few to 10% correspond to potential fluctuations
  $T_e\,\delta n/n\approx0.1$–$1$ V, below which the 50 mV band is buried. A stationary connection with $E\to0$ never forms.
- **Kinetic instability.** A distribution cusped at zero velocity violates the Penrose condition and flattens by quasilinear
  diffusion. Slow electrons stay for 5–10 plasma periods (38 µs).
- **Finite body and convection.** A field line stays connected to the dayside surface for $R_{\rm Moon}/u\approx4$ s, and
  reflected electrons are carried away along it, so the semi-infinite one-dimensional stationary mirror population never
  accumulates.

Which mechanism dominates has not been confirmed directly, because this population is buried in the spacecraft potential of
a few volts and is not observed. Read the field floor $E_{floor}\approx\lvert\mathbf{u}\times\mathbf{B}\rvert$ as the field level
at which this outside physics takes over.

## Outputs to check

| Output | What to check |
|---|---|
| `surface_current_model_upstream_field_deficit_V_m` in `summary.txt` | The field deficit inside the band; the root was accepted with it at most `upstream_field_floor_v_m` |
| `surface_current_model_upstream_negative_band_V` in `summary.txt` | The potential width of the band, a measure of the uncertainty of the wall quantities |
| `surface_current_model_zhao_branch` in `summary.txt` | Whether the expected type was selected. After a rejection, `auto` moves to the next type (B at 60 degrees) |

## Scope

- An accepted root does not follow the idealized electron distribution inside the band. Compared with orbit-consistent
  solutions in which only the slow incoming electrons ($a\lesssim a_s$) lose their drift, the wall potential differed by
  20–80 mV for $a_s=0.07$–$0.2$ and the absorbed electron flux by less than 0.3%
  ([validation record](https://github.com/Nkzono99/sheath-model/blob/main/outputs/upstream_relaxation_20261009/REPORT.md)).
  The uncertainty of the wall quantities is these tens of mV.
- The field floor does not solve the physics that fills the band; it is a boundary condition that hands the responsibility to
  the outside of the one-dimensional electrostatic model. Under conditions whose magnetic-field direction or fluctuation level
  departs from typical solar wind, set $E_{floor}$ to the external field of those conditions.
- The ion temperature (fluid model) does not change the band (for $T_i\le12$ eV the band width and field deficit change by
  less than 1%).
