title: Validating Simulation Results

Lang: [English](ValidationGuide.en.md) | [日本語](ValidationGuide.md)

# Validating Simulation Results

A zero exit status does not establish physical or numerical validity.

## Level 1: execution completed

- exit status is zero;
- `summary.txt`, `charges.csv`, and required histories exist;
- `batches == sim.batch_count`;
- `beachx inspect` can read the output;
- restart model, mesh, and species fingerprints match.

## Level 2: quantitative use is justified

- absorbed, escaped, and max-step populations are physically explainable;
- charge-ledger residual and unresolved discard are checked separately;
- history length supports the claimed steady value; `tol_rel` is monitoring only;
- key observables are stable under smaller `sim.dt`;
- mesh, FMM order/tolerance, outer grid, and sampling are refined as applicable;
- `batch_duration` 0.5x/2x checks do not change the conclusion;
- stochastic conclusions include seed or ensemble sensitivity.

Model-specific diagnostics include periodic cache/zero-mode closure, unified accessibility and orbit errors,
kinetic-solver branch status, and photoelectron ledger/histogram limits.
See [Physics release verification](PhysicsReleaseVerification.en.html) for small reference contracts.
