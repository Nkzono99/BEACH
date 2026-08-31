title: Use the outer 1D1V kinetic oracle for validation

Lang: [日本語](OuterKineticOracle.md) | [English](OuterKineticOracle.en.md)

# Use the outer 1D1V kinetic oracle for validation

`beach-kinetic-response` is an offline research tool that fixes a matching-plane input and tests whether the outer
collisionless 1D1V Vlasov--Poisson system reaches a steady or statistically stationary response. It does not evolve
the outer system with the BEACH runtime. Only a certified, complete Cartesian response subset can be converted for the
existing `response_backend="table"` path.

A successful process exit does not establish physical validity. Check each query's `classification`, grid convergence,
far boundary, and conservation diagnostics before conversion.

## 1. Run the reference sweep

Run this command from the repository root. The output directory must be new or empty.

```bash
beach-kinetic-response \
  examples/outer_kinetic_reference.toml \
  examples/outer_kinetic_k1_k6.csv \
  outputs/outer-kinetic-atlas
```

The example is a diagnostic sweep containing unique Zhao roots, multiple roots, and no-root regions. It is not a
production fixture whose every point is expected to certify. The command writes:

| Output | Purpose |
|---|---|
| `kinetic_response_raw.csv` | Per-query response, classification, stationarity, far-boundary, and conservation data |
| `kinetic_response_manifest.json` | Grid, velocity ranges, CFL, tolerances, commit, source hash, and raw-CSV hash |
| `kinetic_response_profiles/` | Potential, field, charge-density profiles, and time histories |

## 2. Interpret classifications

| `classification` | Table use |
|---|---|
| `steady` | Eligible after grid and domain convergence |
| `stationary_average` | Eligible only after checking autocorrelation and SEM and explicitly allowing it in the converter |
| `unresolved_transient` | Reject; revisit run time or averaging windows |
| `secular` | Reject; do not replace it with an arbitrary-time average |
| `far_boundary_not_converged` | Reject; validate the outer length and far-boundary physics first |
| `numerical_failure` | Reject; inspect `failure_reason`, the velocity domain, and conservation |

The solver imposes $E(H)=D_H/\epsilon_0$ and $\Phi(L)=0$. It uses `E(L)=0` as a convergence diagnostic rather than a
third boundary constraint. Version 1 accepts only zero outward ambient-population inputs.
In the raw CSV, `far_field_v_m` and `far_charge_imbalance` are maximum absolute values over the final averaging window.
Use the corresponding NPZ for signed values and profiles.

`ambient_electron.number_density_m3` and `ambient_ion.number_density_m3` are amplitudes of the complete drifting
Maxwellians feeding the far boundary. Only their $v<0$ half spaces are imposed there, so these values need not equal the
computed far-field densities. Version 1 does not shoot the electron-source amplitude to `E(L)=0`. For a Zhao overlap
test, use the ambient electron source density solved by that Zhao root. Where Zhao has no solution, scan the source
density as an independent convergence axis. Do not interpret a far-boundary failure caused by an inconsistent fixed
amplitude as insufficient outer length.

## 3. Convert a certified subset

The converter checks the raw-CSV hash, every selected classification, duplicates, and the complete Cartesian product.
It does not fill missing points and accepts only `steady` rows by default.

```bash
beach-kinetic-table \
  outputs/outer-kinetic-atlas/kinetic_response_raw.csv \
  outputs/matching_plane_response.csv \
  --range displacement_c_m2=-1e-11:1e-11
```

By default it finds the manifest beside the raw CSV. Use `--manifest` when it is elsewhere. Add
`--allow-stationary-average` only when the study explicitly accepts statistically stationary rows.

The resulting CSV is consumed by the existing
[matching-plane table backend](MatchingPlaneCoupling.en.html#2-build-a-minimal-configuration).

## 4. Compare before production use

Vary at least these controls independently:

- `nz` and each species' `nv`;
- both velocity limits;
- `z_length_m`;
- `cfl`; and
- `max_time_s` and the averaging window.

At points with a unique Zhao root, compare $\Phi_H$, ambient inward fluxes, and PE return and escape fluxes. Then run the
same BEACH case with `zhao_online` and the kinetic-derived table and compare surface charge, gap potential, and absorbed
currents. Finally vary the matching-plane height at three or more points. If the principal observables do not converge,
classify the setup as `no_matching_overlap`.

This oracle is collisionless, unmagnetized, and one-dimensional, and reconstructs PEs as a half-Maxwellian from flux and
mean normal energy. If a strong PE beam prevents far-field quasineutrality, do not pass certification by extending the
domain alone. Test whether collisions, magnetic fields, geometric dilution, or an energy-resolved VDF require a different
outer model.
