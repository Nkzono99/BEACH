---
name: beach-issue-report
description: Help BEACH users prepare GitHub issues for bugs, improvement requests, feature requests, documentation gaps, reproducible simulation problems, and confusing analysis workflows.
---

Use this skill when a user wants to report a bug, request an improvement or feature, turn a failed run into an issue, or prepare a GitHub issue for BEACH.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, TOML keys, file paths, commands, labels, and GitHub fields unchanged.
- If language is unclear, default to Japanese.

## Sources

Assume ordinary users may only have the pip-installed simulator and this plugin, not the full repo. Prefer:

- User-provided command, version, config snippet, log, output listing, expected behavior, and actual behavior.
- Bundled references: `../../references/SPEC.md`, `../../references/README.md`, `../../references/fortran_parameter_file.md`, `../../references/config_workflow.md`, `../../references/python_postprocess_api.md`.
- Bundled docs: `../../docs/known-failure-modes.md`, `../../docs/usage-workflows.md`, `../../docs/skills-guide.md`.

Use repo root docs only when the full checkout is available and you need to confirm current behavior.

## Issue Types

Classify the issue before writing:

- Bug report: reproducible wrong result, crash, parser error, invalid output, documentation mismatch, analysis failure.
- Improvement request: existing workflow is confusing, slow, hard to inspect, or missing diagnostics.
- Feature request: new solver behavior, input option, CLI, output field, surface model, or example.
- Documentation request: unclear parameter, missing example, outdated command, confusing install or analysis guide.

## Required Information

Collect or ask for missing critical fields when possible:

- BEACH version: `beachx --version` if available, package version, or git commit/tag.
- Install method: `pip install beach-bem`, Git URL, editable source build, local binary.
- Execution command: `beach ...`, `fpm run ...`, MPI/Slurm command if used.
- Config type: `case.toml`, preset, rendered `beach.toml`.
- Minimal config or redacted snippet.
- First meaningful error log, not only the final wrapper abort line.
- Expected behavior and actual behavior.
- Output directory listing and relevant `summary.txt` lines.
- Whether the case uses `reservoir_face`, `photo_raycast`, `periodic2`, FMM/treecode, resume, MPI, or Python post-processing.

## Privacy and Size

- Remove secrets, tokens, personal paths, private hostnames if unnecessary, and large generated CSV payloads.
- Prefer short log excerpts plus file listings over full logs.
- For large configs/outputs, include the smallest reproducer or attach files through GitHub rather than pasting huge text.

## Output

If GitHub tooling is available and the user explicitly asks to create the issue, create it in `Nkzono99/BEACH` after showing the sanitized title/body. Otherwise, prepare a copy-ready draft:

```text
## Issue Draft

### Title
...

### Type
Bug report / Improvement request / Feature request / Documentation request

### Body
...

### Suggested Labels
- ...

### Missing Information
- ...
```

For bug reports, use this body shape:

```markdown
## Summary

## Environment
- BEACH version:
- Install method:
- Platform / cluster:
- Compiler / fpm:

## Reproduction
1.

## Expected Behavior

## Actual Behavior

## Input / Configuration

## Output / Logs

## Additional Context
```
