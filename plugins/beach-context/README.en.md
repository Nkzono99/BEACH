Lang: [日本語](README.md) | [English](README.en.md)

# BEACH Context Plugin (Codex / Claude Code)

This self-contained plugin provides BEACH configuration, case design, run diagnosis, output analysis, simulator guidance, and feedback workflows for both Codex and Claude Code. Codex loads `.codex-plugin/plugin.json`; Claude Code loads `.claude-plugin/plugin.json`; both use the shared `skills/` and `references/` trees.

The plugin assumes that many users installed BEACH only through `pip install beach-bem` and cannot read the full repository. For that reason, `references/` bundles snapshots of `SPEC.md`, the config workflow, the `beach.toml` specification, the Python post-processing API, schemas, and representative examples. In a development checkout where the full repository is available, prefer the latest root docs with the same names when needed.

## Install in Codex

To sparse-install only the marketplace metadata and plugin from GitHub:

```bash
codex plugin marketplace add Nkzono99/BEACH \
  --ref main \
  --sparse .agents/plugins \
  --sparse plugins/beach-context
```

At this point only the marketplace is registered, so install the plugin itself:

```bash
codex plugin add beach-context@beach
```

After installation, start a new Codex thread. The plugin skills are available even when Codex is started outside the repository.

To update an already registered marketplace:

```bash
codex plugin marketplace upgrade beach
codex plugin add beach-context@beach
```

To use a local checkout as the marketplace:

```bash
codex plugin marketplace add /path/to/BEACH
codex plugin add beach-context@beach
```

## Install in Claude Code

Install from the GitHub marketplace:

```bash
claude plugin marketplace add Nkzono99/BEACH
claude plugin install beach-context@beach-claude
```

Update an existing installation:

```bash
claude plugin marketplace update beach-claude
claude plugin update beach-context@beach-claude
```

Test the checkout without installing it:

```bash
claude --plugin-dir ./plugins/beach-context
```

The Claude manifest intentionally omits a fixed version so Git-backed installs use the commit SHA for updates.

## Bundled Skills

- `beach-config-review`: Review `beach.toml` configuration files.
- `beach-run-diagnose`: Diagnose install, build, run, abnormal-exit, missing-output, and restart problems.
- `beach-case-design`: Design BEACH configs and sweeps from a physical objective.
- `beach-output-analysis`: Guide analysis of `outputs/latest`, CSV files, histories, Python APIs, and `beachx`.
- `beach-simulator-guide`: Provide learning and usage guides for BEACH users.
- `beach-method-summary`: Draft method descriptions for papers, presentations, and README-style docs.
- `feedback-beach`: Organize bug reports, improvement requests, feature requests, and documentation feedback into GitHub Issue drafts.

Claude Code also discovers the specialized subagents under `agents/`. They inherit the active model instead of pinning an older model alias.

See [docs/skills-guide.en.md](docs/skills-guide.en.md) for detailed skill selection guidance.

## Bundled References

`references/` includes the following files so the plugin can support users without reading the full repository:

- `README.md`, `SPEC.md`
- `agent-user-guide.md`, `fortran_workflow.md`
- `config_workflow.md`, `fortran_parameter_file.md`
- `python_postprocess_api.md`
- `fortran_fmm_core.md`, `batch_duration_stability.md`
- `schemas/beach.schema.json`
- `examples/tutorial_insulator.toml`, `examples/beach.toml`, `examples/periodic2_basic/beach.toml`

For pre-run checks, prefer `beachx lint beach.toml`; for pre-run estimates, use `beachx estimate-workload`; for output inspection, use `beachx inspect outputs/latest`.

## Distribution Policy

BEACH-specific knowledge lives in this plugin. Cross-cutting Fortran, fpm, Slurm, and Codex operation workflows should live in shared plugins or the repo-root `AGENTS.md`.
