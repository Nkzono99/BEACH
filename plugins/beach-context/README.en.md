Lang: [日本語](README.md) | [English](README.en.md)

# BEACH Context Plugin

This repo-local plugin gives Codex BEACH-specific context for configuration files, run workflows, output analysis, known failure modes, simulator guidance, and issue reporting. Skills are expected to respond in the user's language, including Japanese and English, while keeping code identifiers, commands, file names, and parameter names unchanged.

The plugin assumes that many users installed BEACH only through `pip install beach-bem` and cannot read the full repository. For that reason, `references/` bundles snapshots of `SPEC.md`, the config workflow, the final `beach.toml` specification, the Python post-processing API, schemas, representative examples, and built-in presets. In a development checkout where the full repository is available, prefer the latest root docs with the same names when needed.

## Installation

To sparse-install only the marketplace metadata and plugin from GitHub:

```bash
codex plugin marketplace add Nkzono99/BEACH \
  --ref main \
  --sparse .agents/plugins \
  --sparse plugins/beach-context
```

At this point the marketplace is registered, but the `BEACH Context` plugin is not enabled yet. Start Codex, open `/plugins`, and install `BEACH Context`.

```bash
codex
# Open /plugins inside Codex
```

After installing the plugin, restart Codex. The plugin skills will then be available even when Codex is started outside the repository.

To update an already registered marketplace:

```bash
codex plugin marketplace upgrade beach
```

To use a local checkout as the marketplace:

```bash
codex plugin marketplace add /path/to/BEACH
```

In this case too, install `BEACH Context` from `/plugins` after registering the marketplace.

## Bundled Skills

- `beach-config-review`: Review `beach.toml`, `case.toml`, and presets.
- `beach-run-diagnose`: Diagnose install, build, run, abnormal-exit, missing-output, and restart problems.
- `beach-case-design`: Design BEACH cases, presets, and sweeps from a physical objective.
- `beach-output-analysis`: Guide analysis of `outputs/latest`, CSV files, histories, Python APIs, and `beachx`.
- `beach-simulator-guide`: Provide learning and usage guides for BEACH users.
- `beach-method-summary`: Draft method descriptions for papers, presentations, and README-style docs.
- `beach-issue-report`: Turn bug reports, improvement requests, and feature requests into GitHub Issue drafts.

See [docs/skills-guide.en.md](docs/skills-guide.en.md) for detailed skill selection guidance.

## Bundled References

`references/` includes the following files so the plugin can support users without reading the full repository:

- `README.md`, `SPEC.md`
- `agent-user-guide.md`, `fortran_workflow.md`
- `config_workflow.md`, `fortran_parameter_file.md`
- `python_postprocess_api.md`
- `fortran_fmm_core.md`, `batch_duration_stability.md`
- `schemas/beach.schema.json`, `schemas/beach.case.schema.json`, `schemas/beach.preset.schema.json`
- `examples/beach.toml`, `examples/periodic2_basic/case.toml`
- built-in preset snapshots under `presets/...`

For pre-run checks, prefer `beachx config validate`; for rendered-config estimates, use `beachx estimate-workload`; for output inspection, use `beachx inspect outputs/latest`.

## Distribution Policy

BEACH-specific knowledge lives in this plugin. Cross-cutting Fortran, fpm, Slurm, and Codex operation workflows should live in shared plugins or the repo-root `AGENTS.md`.
