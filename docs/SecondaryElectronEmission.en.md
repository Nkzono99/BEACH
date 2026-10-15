title: Secondary electron emission

Lang: [日本語](SecondaryElectronEmission.md) | [English](SecondaryElectronEmission.en.md)

# Secondary electron emission

Secondary electron emission (SEE) creates electrons from the energy and angle of a particle absorbed on a surface.
Set the emission model on the incident species and the temperature on the target species; BEACH tracks the secondaries
in the frozen field of the same batch. Reaction charge, return, and escape can be checked in the charge budget.

## What it solves

For incident mass $m$ and velocity $\mathbf v$, the energy is $E=m|\mathbf v|^2/(2e)$ in eV.
With unit normal $\hat{\mathbf n}$ pointing toward the incident side, use $c=\max(0.1,-\mathbf v\cdot\hat{\mathbf n}/|\mathbf v|)$.
The `whipple` yield is

$$
\delta = \frac{1.114\delta_{max}}{c}\left(\frac{E_{max}}{E}\right)^{0.35}
\left[1-\exp\left(-2.28c\left(\frac{E}{E_{max}}\right)^{1.35}\right)\right]
$$

(Whipple 1981; Marchand et al. 2014, Eq. (2)). The yield is zero for $E\le0$.
`constant` uses `const_yield` independently of incident energy and angle; it can also be assigned to ions.

## Assumptions

The emitted macro count is the integer part of $\delta|q_{in}/q_{see}|$ plus a Bernoulli draw with the fractional probability.
Each secondary inherits the incident weight $w_{in}$. The velocity sampler is shared with photoelectrons:
the normal component follows a flux-weighted half-Maxwellian; the tangential components are zero-mean Gaussians
truncated at $6\sigma$. Normal drift is zero and the temperature comes from the target species.

## Treatment in BEACH

After absorbing the incident particle and depositing $q_{in}w_{in}$ on the hit element, generate secondaries
$10^{-12}$ m from the collision point along the emission normal. Canonicalize periodic coordinates.
For each secondary, always deposit $-q_{see}w_{in}$ on the same element, balancing emitted particle charge.
This uses the photoelectron surface-emission charge channel.

Append and track each generation in the same frozen field, including absorption, escape, and the step limit.
An absorbed secondary may itself have emission enabled. Primaries are generation 0; a particle at the configured
maximum generation does not emit further. Record the would-be count in diagnostics without depositing its reaction charge.
Random draws depend on the seed, batch, rank, parent particle, step, generation, and child ordinal, independently of
thread count and execution order.

For an external sheath, assign the SEE species to `sheath.species.secondary_electron`. It receives the same top barrier
as photoelectrons ($\phi_m$ in Type A, 0 V in B/C) and the same uniform return over a periodic cell.
No SEE-specific term is added to the outer root: sum the PE and SEE top-outflow spectra.
For normalization and fixed-current allocation, see [the stationary outer sheath](ZhaoStationaryClosure.en.html#refresh-the-outer-root).
The initial root uses only the PE surface half-Maxwellian; without PE it starts from zero-source Type C.
SEE enters at the first refresh.

## Configuration

The generation limit is `particles.secondary_emission.max_generations` (default `3`, nonnegative); primaries are generation 0.

Add this difference to an existing incident species:

```toml
[particles.species.surface_interaction.secondary_emission]
enabled = true
model = "whipple"
yield_max = 1.5
energy_max_ev = 500.0
target_species = "secondary_electron"
```

The target is a species with no independent injection source:

```toml
[[particles.species]]
species_key = "secondary_electron"
charge_c = -1.602176634e-19
mass_kg = 9.1093837139e-31

[particles.species.source]
mode = "secondary_emission"

[particles.species.distribution]
temperature_ev = 2.0

[particles.species.boundary]
z_high = "open"
```

Omit `sampling.*`, `inflow.*`, emission current density, ray count, and `charging.closure`.
See [parameters](Parameters.en.html) for all defaults and constraints.
Complete configurations are [`zero_current_secondary.toml`](../examples/zero_current_secondary.toml) and
the night-side [`zero_current_no_photo_secondary.toml`](../examples/zero_current_no_photo_secondary.toml).

## Outputs to check

`secondary_emission_history.csv` records hits, emitted count and charge, mean yield, and the generation-capped count
for each incident species in every accepted batch. `summary.txt` contains cumulative values;
`charge_ledger.csv` records SEE emission, reabsorption, and escape. With outer-root refresh, check
`secondary_share` and `secondary_budget_residual_a_m2` in `matching_plane_history.csv`.
See [the output reference](OutputReference.en.html) for units and column definitions.

## Scope

There are no material or `mesh_id` dependent parameters, elastic/inelastic reflection, or surface conduction.
The cosine floor of 0.1 is a provisional divergence safeguard and does not establish material behavior at grazing incidence.
Phase 1 fixed-current allocation is provisional. Only the PE escape target is allocated by the PE fraction in the latest window;
SEE remains fully tracked, so its budget residual is not constrained to zero.
Assigning `fixed_current` / `neutral_return` to a SEE species is a configuration error.
Particles finish within the batch; their flight state is not carried into the next batch.

Workload estimates count only independent injection sources and exclude collision-generated SEE.
With SEE enabled they print `secondary_emission_scope=excluded_collision_generated_particles`.
