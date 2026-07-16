---
name: beach-documentation-writing
description: Write, revise, and review repository documentation for BEACH, including README files, tutorials, user guides, reference pages, and numerical-method pages. Use for documentation changes that improve clarity, navigation, expected results, canonical ownership, or Japanese/English consistency; prefer a dedicated method-summary skill when available if the task only produces a standalone paper or presentation summary.
---

# BEACH Documentation Writing

Write for the reader's task while preserving BEACH's documented behavior.

If the requested output is a standalone paper paragraph, presentation summary, or method overview rather than a repository documentation change, use `beach-method-summary` when that skill is available. Otherwise, continue with this skill's sources and constraints without editing repository documentation.

## Sources

Inspect the target document, its Japanese or English counterpart, nearby navigation, `SPEC.md`, and the relevant Fortran implementation before changing behavioral claims.

Read `references/style-guide.md` before drafting or reviewing prose. Apply its document-type rules instead of imposing one narrative style on every page.

## Workflow

1. Classify the page as an entry page, tutorial, task guide, reference, explanation, or contributor guide.
2. Identify the reader's immediate question and the observable result that answers it.
3. Find duplicated instructions and choose one canonical page. Keep summaries elsewhere short and link to the canonical page.
4. Revise claims before wording. Preserve uncertainty, model scope, identifiers, commands, and expected outputs.
5. Match the structure to the page type. Keep reference pages searchable; reserve cognitive rhythm for explanations grounded in physical behavior, numerical evidence, or tradeoffs.
6. Review the paired `.md` and `.en.md` files. Keep their behavior, commands, warnings, and navigation equivalent even when sentence structure differs.
7. Inspect the final diff for broken links, unsupported claims, hidden prerequisites, and prose-only narration such as "next we will see."

## BEACH Constraints

- Treat the Fortran implementation and `SPEC.md` as the source of truth for simulation behavior.
- Keep code identifiers, TOML keys, commands, file paths, CSV columns, and algorithm names unchanged.
- State expected outputs concretely when documenting a runnable example.
- Distinguish a completed run from a numerically or physically validated result.
- Do not describe `tol_rel` as an early-stop condition in the current implementation.
- Do not imply that the default insulator-accumulation workflow validates conductor, resistive, or other extension models.
- Present `kinetic_1d` as the standard self-consistent outer-sheath model. Present `unified_linear_response` separately as an advanced rough-surface linear-screening model with conditional applicability, not as a higher-accuracy replacement.

## Review Output

Follow the repository language policy and review in Japanese. Lead with correctness or navigation problems, cite file and line locations, and separate required fixes from optional polish. When editing files, summarize the reader-visible outcome and the checks performed.
