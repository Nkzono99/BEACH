# BEACH Plugin Skill Guide

For BEACH-specific user support, choose skills as follows.

| Skill | Use Case |
| --- | --- |
| `beach-config-review` | Review `beach.toml`, `case.toml`, presets, schema validation, and input consistency. |
| `beach-run-diagnose` | Diagnose failed install/build/run attempts, missing outputs, restart issues, or abnormal results. |
| `beach-case-design` | Design new cases, presets, parameter sweeps, and settings from a research objective. |
| `beach-output-analysis` | Read `outputs/latest`, use `beachx`, use Python APIs, and analyze histories. |
| `beach-simulator-guide` | Explain BEACH concepts, learning paths, docs, and user tutorials. |
| `beach-method-summary` | Draft method descriptions for papers, presentations, README, and design notes. |
| `beach-issue-report` | Turn bugs, improvements, features, and documentation gaps into GitHub Issue drafts. |

In a development checkout where the full repository is available, prefer root docs and source files. Treat plugin `references/` as snapshots for users working outside the repository.

Shared rules:

- Respond in the user's main language. Default to Japanese when Japanese appears.
- Do not translate code identifiers, TOML keys, file paths, commands, or CSV columns.
- Do not invent defaults or undocumented behavior.
- Treat the Fortran implementation and `SPEC.md` as the source of truth for simulation behavior.
- `tol_rel` is an output monitoring metric, not an early-stop condition in the current implementation.
