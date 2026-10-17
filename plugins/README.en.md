Lang: [日本語](README.md) | [English](README.en.md)

# BEACH Codex Plugins

This directory contains Codex plugins that package the context needed to use, analyze, and maintain BEACH.

## Available Plugins

| Plugin | Contents |
| --- | --- |
| [beach-context](beach-context/README.en.md) | Skills and bundled references for config review, run diagnosis, case design, output analysis, method summaries, simulator guidance, and feedback reporting |

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

After installing the plugin, restart Codex. The `beach-context` skills will then be available even when Codex is started outside the repository.

To update an already registered marketplace:

```bash
codex plugin marketplace upgrade beach
```

To use a local checkout as the marketplace:

```bash
codex plugin marketplace add /path/to/BEACH
```

In this case too, install `BEACH Context` from `/plugins` after registering the marketplace.

## Skill Visibility

The repository root `AGENTS.md` contains project-local instructions for BEACH developers. It is loaded only when Codex is started inside the BEACH repository.

The `plugins/beach-context/skills/` directory contains user-facing plugin skills. After installing the plugin from `/plugins` and restarting Codex, these skills are available from other working directories such as `~`.

## runops Integration

`beach-bem` publishes the `beach` adapter for [runops](https://github.com/Nkzono99/runops) 2.x
through the `runops.adapters` entry point (`beach.runops_adapter:Adapter`). When `beach-bem` is
installed in a runops project's `.venv`, `runo init beach` and `runo update` find the adapter and
copy the skills and references from `plugins/beach-context` (shipped in the wheel as
`beach.agent_context`) into the project's `.agents/skills/` and `.claude/skills/`. No plugin
installation or restart is needed.

For each run the adapter applies the `--set` dotted keys (e.g. `run.batch_count`) to
`input/beach.toml` and sets the output directory to `work/latest`. Keys the schema does not define
produce a warning at creation. A run is completed when `summary.txt` reports `batches` equal to
`run.batch_count`. `runo new --from RUN --continue --set run.batch_count=<cumulative target>`
resumes from the source run's checkpoint. `runo tidy --purge-work` keeps the final checkpoint, so a
purged run can still be continued.

## Placement Policy

BEACH-specific physics, configuration specifications, output specifications, known failure modes, and learning paths live in this simulator repo plugin. Cross-cutting Fortran, fpm, Slurm, and Codex operation workflows should live in repo-root `AGENTS.md` or shared plugins.
