title: Validating Simulation Results

Lang: [English](ValidationGuide.en.md) | [日本語](ValidationGuide.md)

# Validating Simulation Results

A successful process exit does not by itself establish numerical or physical validity. Check three questions in order:

1. Did the configured run finish, with explainable particle and charge balances?
2. Are the observables of interest stable when the numerical resolution is changed?
3. Does the selected physical model apply, and does the conclusion remain within its scope?

Use [Inspect Output Files](OutputGuide.en.html) to locate each value. This page uses those values to decide whether a run
is acceptable.

## Define acceptance criteria first

Before running convergence cases, define what “unchanged” means:

- primary observables: total surface charge, potential distribution, absorbed or escaping flux, or the quantity used directly in the scientific conclusion;
- tolerance: an acceptable absolute or relative difference in each primary observable;
- evaluation interval: a final value, a post-transient mean, a maximum, or another declared statistic;
- statistical uncertainty: variation across random seeds or an ensemble of otherwise identical runs.

Do not adjust the criterion after seeing a convenient result. BEACH has no universal convergence tolerance for all cases;
the appropriate criterion depends on the purpose of the simulation.

## 1. Confirm that the run completed

```bash
beachx inspect outputs/latest
```

At minimum, check the following:

| Location | Acceptance condition |
| --- | --- |
| process exit status | `0` |
| `summary.txt` | `batches == sim.batch_count` |
| `charges.csv` | final element charges are present |
| required histories | records are complete at the configured stride |
| restarted run | mesh fingerprint matches; any model or species fingerprint warning is recorded |

These checks establish only that the requested computation completed. They do not establish convergence of a physical quantity.

## 2. Check particle and charge balances

Confirm that particle counts in `summary.txt` are consistent with the configured sources and boundaries.

| Field | Interpretation |
| --- | --- |
| `absorbed` | particles absorbed at a surface |
| `escaped_boundary` | particles leaving through the box boundary or outer model |
| `survived_max_step` | particles whose absorption or escape was unresolved at `sim.max_step` |

`survived_max_step` is neither absorption nor escape. If it is large enough to affect the conclusion, revise
`sim.max_step`, `sim.dt`, or the box dimensions and show that the unresolved population becomes negligible.

Then inspect the charge balance in `charge_ledger.csv` and `summary.txt`:

- `charge_ledger_residual_C`: conservation residual including injection, emission, surface absorption, and escape to infinity;
- `charge_ledger_discarded_unresolved_abs_C`: absolute unresolved charge summed without cancellation between species.

A small conservation residual does not compensate for a large unresolved charge. Assess the two values separately.

## 3. Check time evolution and statistical variation

Use `charge_history.csv` and, when needed, `potential_history.csv` to assess whether the system has reached the regime
used in the conclusion. Inspect an evaluation interval rather than only its last sample:

- no remaining systematic increase or decrease;
- no unexplained batch-to-batch oscillation or excursion;
- scatter around the interval mean is below the declared tolerance;
- increasing `sim.batch_count` does not change the mean or conclusion.

`last_rel_change` and `sim.tol_rel` are monitoring values. In the current implementation they do not stop a run early,
and `last_rel_change < tol_rel` alone does not establish a steady state.

Separate Monte Carlo variation from time-discretization effects.

| Quantity varied | Main effect tested |
| --- | --- |
| `sim.rng_seed` | random-sampling variation |
| macro-particle count or particle weight | Monte Carlo noise |
| `sim.batch_count` | adequacy of the simulated duration |
| `sim.batch_duration` or `sim.batch_duration_step` | stability of the end-of-batch surface-charge update |

As a practical check, compare `batch_duration` at 0.5x and 2x the baseline. See
[`batch_duration` Stability and Steady Values](BatchDurationStability.en.html) for the underlying interpretation.

## 4. Check numerical-resolution dependence

Copy the baseline case and normally vary one setting at a time. Compare the same primary observables in every case,
and require the changes to fall below the previously declared tolerances.

| Convergence axis | Example comparison | Error being tested |
| --- | --- | --- |
| particle time step | halve `sim.dt` | orbit, collision location, absorption, and escape |
| tracking length | increase `sim.max_step` | truncation from unresolved particles |
| surface mesh | refine the triangles | spatial discretization of charge and potential |
| field solver | compare with Direct on a small case | Treecode or FMM approximation |
| finite periodic images | increase `field_periodic_image_layers` | image-sum truncation |
| outer model | increase grid points or surface samples | outer profile and geometry sampling |

Converging a path integral or image shell while holding the mesh fixed does not quantify mesh-discretization error.
Record which convergence axes were tested and which remain unevaluated.

## 5. Check diagnostics for the selected physical model

In addition to the common checks, inspect diagnostics specific to the active model.

| Configuration | Additional checks |
| --- | --- |
| finite-image `periodic2` | dependence on `field_periodic_image_layers`; do not extrapolate a finite-image result to infinite periodicity |
| `cached_kneq0` | cache fingerprint, cold/warm agreement, zero mode, and Gauss residual |
| `infinity_barrier` | mean reservoir-face potential, image-layer dependence, and warnings for in-plane potential variation |
| `potential_barrier` | open-face crossing potential, `reservoir.phi_infty`, and the normal-kinetic-energy reflection/escape decision |
| photoelectrons | emission, return, and escape charge balance; convergence with ray sampling and time step |
| `matching_plane_quasistatic` | common: fixed-point residual and iterations, PE `outward = return + escape`, relaxation, and matching-height dependence; table: response-grid dependence; online: branch policy, root solve, and moment reduction |

Definitions and applicability limits are documented in [Finite Periodic Configuration](FinitePeriodicConfiguration.en.html)
and [Particle Escape and Return](ParticleEscapeReturn.en.html). See
[Quasistatic Matching-Plane Coupling](MatchingPlaneCoupling.en.html) for model selection and scope, and the
[matching-plane numerical and response-table reference](MatchingPlaneReference.en.html) for table generation,
fixed-point checks, and the height sweep.

## 6. Limit the physical conclusion to what was tested

Finally, decide whether the numerically stable result supports the intended physical claim.

- Enumerate every input difference between compared cases beyond the intended model change.
- State the boundary conditions, particle sources, surface model, and field-evaluation method.
- Do not treat a finite box as infinity, finite simulated time as a steady state, or a finite image sum as an infinite periodic system.
- Reflect numerical error, Monte Carlo variation, and model applicability in the reported precision or uncertainty.

At minimum, retain the input file, output directory, primary observables and acceptance criteria, and the differences
observed in each convergence case. Small fixtures and HPC gates used to qualify a BEACH code release are separate from
user-case validation and are documented in [Physics Release Verification](PhysicsReleaseVerification.en.html).
