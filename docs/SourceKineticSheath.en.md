title: Source-connected sheath response

Lang: [日本語](SourceKineticSheath.md) | [English](SourceKineticSheath.en.md)

# Source-connected sheath response

Select `density_model="source_kinetic"` to compute the outer matching-plane density from boundary-connected particle
orbits. This reference defines the model, root classification, and diagnostic CSV. It does not select a dynamically stable branch.

## Configuration and compatibility

Use the following settings with an otherwise valid [matching-plane case](MatchingPlaneCoupling.en.html).
The runnable dark example is `examples/periodic2_matching_plane_source_kinetic.toml`.
`beach examples/periodic2_matching_plane_source_kinetic.toml` runs two batches and writes the summary and histories
to `outputs/periodic2_matching_plane_source_kinetic/`.

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
density_model = "source_kinetic"
zhao_branch = "auto"
zhao_root_selection = "require_unique"
```

The default `density_model="zhao_legacy"` preserves the released BEACH density implementation and root policies. Source-connected
kinetics is the candidate for the primary physical model; existing inputs retain their original model while validation
continues. There is no scheduled removal of this compatibility behavior. Stationary Zhao and table configurations reject `density_model`.

The initial model is one-dimensional, collisionless, and unmagnetized, with zero-drift ambient Maxwellian electrons,
half-Maxwellian emitted electrons, and cold nonreflecting ions. Set ambient electron `drift_velocity=[0,0,0]`, ion
`temperature_ev=0`, and inward ion drift. Nonzero electron drift is unsupported; no legacy drift factor is substituted.
Trapped orbits disconnected from either source are empty. Ambient outward feedback remains transparent.

`zhao_branch` accepts `auto`, `a`, `n`, `b`, `c`, and `b0`. The solver searches and validates all candidates before filtering
by Type, then accepts exactly one detected root. Multiplicity and unresolved queries cannot define a response and stop
the run. No alternative Type, legacy density, or finite-boundary model is substituted. `minimum_energy` and `continuation`
are unsupported for this model. A single detected root does not establish mathematical uniqueness.

## Relation to Zhao's original formulation

It is inaccurate to label Zhao's entire electron model a Boltzmann density model.
[Zhao et al. (2020), §III.A, equations (3)–(6)](https://scholarworks.indianapolis.iu.edu/server/api/core/bitstreams/e20ede43-d66d-4b73-b89c-3d3192672188/content)
integrate Maxwellian distributions over distinct passing, reflected, and returning velocity populations.
At zero drift, equation (3) is $n_{e,f}=(A_e/2)e^\phi\operatorname{erfc}(\sqrt{\phi-\phi_m})$.
Applying it to Type B with the path minimum $\phi_m=0$ gives the source-connected density below.
This applies the original population formula to zero-drift Type B; it does not replace an entirely Boltzmann Zhao theory.

The existing BEACH `evaluate_zhao_density_hat` Type B has no position-dependent velocity cutoff and reduces to
$(A_e/2)e^\phi$ at zero drift. Source mode corrects that difference. Existing A/C densities already retain velocity
cutoffs and are not simply Boltzmann. The name `zhao_legacy` denotes reproduction of the previous BEACH implementation,
not fidelity to the original paper.

[Zhao's 2022 dissertation](https://scholarsmine.mst.edu/doctoral_dissertations/3176/), §4.3.1, equation (4.29), and
Appendix equation (1) retain the same velocity cutoff and explicitly define the minimum over the entire profile.
The full text of the 2021 IEEE TPS version
[DOI:10.1109/TPS.2021.3110946](https://doi.org/10.1109/TPS.2021.3110946) was not retrieved, so its Type B-specific
equations are not claimed to have been checked. Exact transport of a nonzero-drift boundary VDF still requires separate validation.

## Normalization and densities

Use $\phi=e(\Phi-\Phi_{out})/(k_BT_e)$, $x=(z-H)/\lambda_{De}$, $E=-\phi'$, $\tau=T_{ph}/T_e$, and
$G=\Gamma_{ph,H}^{out}/[n_i\sqrt{k_BT_e/m_e}]$. The legacy Zhao internal temperature ratio is reciprocal.
Normalize the physical field $D_H/\epsilon_0$ by $T_e[\mathrm{eV}]/\lambda_{De}$. The emitted temperature is the measured
mean normal energy at the interface, and $G$ is the outward flux crossing that interface.

Type B retains the velocity cutoff of accelerated incoming electrons:

$$
n_e(\phi)=\frac{A_e}{2}\operatorname{erfcx}(\sqrt\phi).
$$

Let $\phi_*$ be the path minimum and $q=G\exp[-(\phi_H-\phi_*)/\tau]$. Inside and outside an A/N minimum,

$$
n_{e,in}=\frac{A_e}{2}e^{\phi_*}\operatorname{erfcx}(\sqrt{\phi-\phi_*}),\qquad
n_{e,out}=A_e e^\phi-n_{e,in},
$$
$$
n_{ph,out}=\sqrt{\frac\pi{2\tau}}q\operatorname{erfcx}(\sqrt{(\phi-\phi_*)/\tau}),\qquad
n_{ph,in}=2\sqrt{\frac\pi{2\tau}}G e^{(\phi-\phi_H)/\tau}-n_{ph,out}.
$$

B uses the inner formula with $\phi_*=0$; C uses the outer formula with $\phi_*=\phi_H$.
Cold ions have $n_i=(1-2\phi/M^2)^{-1/2}$, requiring $\max\phi<M^2/2$ over the whole path. The turning limit is never clipped.

The solver imposes $\phi(\infty)=E(\infty)=\rho(\infty)=0$ and the specified $E_H$, and solves the ambient amplitude
$A_e$ from neutrality. This differs from prescribing an upstream VDF amplitude. Current
$J=M/\sqrt{m_i/m_e}+q-A_e e^{\phi_*}/\sqrt{2\pi}$ is an output; $J=0$ is not imposed. Validation checks neutrality,
first integrals, sampled path $E^2$ and charge extrema, minimum curvature, ion turning margin, emission barriers, and edge accessibility.

| Type | Potential profile | Path minimum |
|---|---|---|
| B | Positive interface, monotonic decrease to zero | 0 |
| C | Negative interface, monotonic increase to zero | $\phi_H$ |
| A | $\phi_H>0>\phi_m$, interior minimum | $\phi_m$ |
| N | $\phi_m<\phi_H\le0$, interior minimum | $\phi_m$ |
| B0 | Homogeneous $\phi=E=0$ | 0 |

Zero-field queries search both B0 and nonuniform B/C boundary extrema. B0 requires
$A_e=2[1-\sqrt{\pi/(2\tau)}G]>0$. Electron access and PE barrier potentials are $\Phi_m$ for A/N and the existing 0 V gauge for B/C/B0.

## Candidate diagnostics

Passing the source model to `beach-zhao-atlas` writes one row per detected candidate instead of the legacy one-row-per-A/B/C
format. The input uses the same [three-column atlas query contract](MatchingPlaneReference.en.html).

```bash
beach-zhao-atlas source.toml queries.csv source-atlas.csv
```

`classification` is `solutions`, `analytically_excluded`, `search_unresolved`, `numerical_exception`, or `invalid_input`.
`numerical_failures` counts nonfinite candidate evaluations, including when another valid root was found.
Operational exclusions such as unsupported drift or selection policies fail configuration validation; they are not analytical nonexistence claims.
`root_count` includes accepted
roots of all Types. Rejected candidates have `accepted=F` and a reason. Queries without candidates retain a row with NaN
potentials. `phi_h_te` and `phi_min_te` are in $T_e/e$ units; the latter is the physical path minimum, not a legacy monotonic
placeholder. The CSV also records amplitude, escaping/incoming fluxes, current, residuals, bounds, and barrier roundoff.

The finite search spans depths $10^{-24}$ to $10^8$ with a base of 481 samples and additional shallow/deep samples. It refines
crossings, sampled extrema, and physical-domain boundaries. Only analytical exclusions set `search_complete=T`.
For $M\ge1,E_H>0$, $E_H^2<2\sqrt{2\pi\tau}G$ is necessary. Type B additionally requires
$q\ge\sqrt{2/\pi}\tau/(1+\sqrt\tau)$; this is not a nonexistence certificate for other Types.

`deep_root=T` marks depth above $10^4$. `barrier_resolution_limited=T` marks
$\epsilon_{mach}(|\phi_H|+|\phi_m|)/(\phi_H-\phi_m)>10^{-7}$. Deep roots are not counted as practical extensions of lunar
applicability. Multiply the potential span by $T_e$ to assess the nonrelativistic approximation separately.
The atlas is not a response table. `beach-zhao-response` writes only a complete Cartesian product after root selection;
missing or ambiguous responses are not filled. Generated tables record `# density_model=...`; online summaries record the density model too.

## Validation scope

The independent `lunar-sheath-solver-v0.1.0.zip` distribution passed 86 tests separately from BEACH's Fortran regressions.
All 739 roots across its 1,211 conditions agreed in count and Type; the maximum scaled difference in potentials, amplitude,
fluxes, and current was $3.6\times10^{-13}$. Regressions cover independent B/A/N velocity integrals, A/N Poisson first
integrals, B/N coexistence, B0 and a nonuniform zero-field C endpoint, depth truncation, and SI conversion.
See the [validation example README](../examples/source_kinetic_validation/README.en.md) for reproduction and scope.

For a small dark BEACH charging fixture, source online and a 129-row table generated with the same model agreed within
0.2% in matching potential and displacement. This coupling regression does not establish spatial or temporal convergence of PE charging.

Fourteen existing kinetic-oracle probes with each coexistence root's own $A_e$ produced three `steady` and eleven
`far_boundary_not_converged` results. Refining B's velocity grid and extending N's outer length brought potentials closer
to the static roots, but these probes do not establish a stable branch, the semi-infinite limit, or grid convergence.
Unconverged results were not promoted to response tables.

These checks validate static roots. Comparisons with the [time-dependent kinetic oracle](OuterKineticOracle.en.html) must
match the upstream electron source amplitude to each root's $A_e$. Spatial/velocity grids, the cold-ion limit, outer length,
time convergence, and BEACH charging and PE return/escape convergence require separate checks. Root order and small static
residuals do not establish time-dependent stability.
