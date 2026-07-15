Lang: [日本語](README.md) | [English](README.en.md)

# BEACH Plugins for Codex and Claude Code

This directory contains plugins that package the context needed to use, analyze, and maintain BEACH in Codex and Claude Code.

## Available Plugins

| Plugin | Contents |
| --- | --- |
| [beach-context](beach-context/README.en.md) | The shared Codex/Claude Code plugin with skills, agents, and bundled references for config review, run diagnosis, case design, output analysis, method summaries, simulator guidance, and feedback reporting |
| [beach-context-claude](beach-context-claude/README.md) | Compatibility package for manually copied legacy Claude Code project agents and commands |

## Install in Codex

Sparse-install the marketplace metadata and plugin from GitHub:

```bash
codex plugin marketplace add Nkzono99/BEACH \
  --ref main \
  --sparse .agents/plugins \
  --sparse plugins/beach-context
codex plugin add beach-context@beach
```

Start a new Codex thread after installation. The skills are available outside the BEACH repository too.

Update an existing installation:

```bash
codex plugin marketplace upgrade beach
codex plugin add beach-context@beach
```

Use a local checkout as the marketplace:

```bash
codex plugin marketplace add /path/to/BEACH
codex plugin add beach-context@beach
```

## Install in Claude Code

The formal Claude Code distribution uses the same [beach-context](beach-context/README.en.md) plugin:

```bash
claude plugin marketplace add Nkzono99/BEACH
claude plugin install beach-context@beach-claude
```

Test the checkout directly with:

```bash
claude --plugin-dir ./plugins/beach-context
```

The [beach-context-claude](beach-context-claude/README.md) directory remains for workflows that manually copy legacy project-local agents and commands.

## Instruction and Skill Visibility

The repository-root `AGENTS.md` and `CLAUDE.md` contain project-local instructions for BEACH developers.
The formal plugin shares `skills/` and `references/` across both products, while Claude Code also discovers `agents/`.

## Placement Policy

BEACH-specific physics, configuration specifications, output specifications, known failure modes, and learning paths live in this simulator plugin. Cross-cutting Fortran, fpm, and Slurm workflows belong in repository-root guidance or shared plugins.
