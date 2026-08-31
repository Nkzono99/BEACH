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

1. Classify the page as exactly one of: entry page, tutorial, task guide, reference, explanation, or contributor guide.
   If it mixes page types, split the material or move details to an existing canonical page before polishing prose.
2. Identify the reader's immediate question and the observable result that answers it.
3. For a human-facing page, state its question, one-sentence answer, and reader outcome near the opening. Explain the
   ordinary path first. Put options, deprecated behavior, implementation details, and exhaustive constraints later or
   move them to reference and internals pages.
4. Classify existing sections internally as KEEP, SHORTEN, MOVE, or DELETE. Preserve information by linking to its
   canonical owner; do not keep a dense copy merely to make every page self-contained.
5. Find duplicated instructions and choose one canonical page. Keep summaries elsewhere short and link to the canonical page.
6. Revise claims before wording. Preserve uncertainty, model scope, identifiers, commands, and expected outputs.
7. Match the structure to the page type. Keep reference pages searchable; reserve cognitive rhythm for explanations grounded in physical behavior, numerical evidence, or tradeoffs.
8. Review the paired `.md` and `.en.md` files. Keep their behavior, commands, warnings, and navigation equivalent even when sentence structure differs.
9. Inspect the final diff for broken links, unsupported claims, hidden prerequisites, and prose-only narration such as "next we will see."

## Scope Boundary

Entry pages, tutorials, task guides, and explanations help a researcher build a correct mental model or complete one
task. They are not exhaustive specifications. Do not organize them around Fortran modules, internal arrays, private
APIs, MPI/OpenMP reduction steps, RNG rollback, or deprecated inputs. Mention those only when they change the user's
decision, then link to a reference or contributor page.

Reference pages, schemas, `SPEC.md`, migration notes, and API documentation may enumerate every identifier, default,
unit, constraint, equation, and unsupported combination. Do not copy that enumeration into a human-facing page.

Prefer Japanese concept names in Japanese prose and add exact identifiers only where readers configure, inspect, or
search for them. For example, write 電荷差分の確定反映（commit） and バッチ中に固定する場（field snapshot）
at first use instead of using internal English terms as the explanation itself.

## BEACH Constraints

- Treat the Fortran implementation and `SPEC.md` as the source of truth for simulation behavior.
- Keep code identifiers, TOML keys, commands, file paths, CSV columns, and algorithm names unchanged.
- State expected outputs concretely when documenting a runnable example.
- Distinguish a completed run from a numerically or physically validated result.
- Do not describe `tol_rel` as an early-stop condition in the current implementation.
- Do not imply that the default insulator-accumulation workflow validates conductor, resistive, or other extension models.
- Do not present a self-consistent outer-plasma or outer-sheath model as supported. Document the current local `reservoir_face` plus closed-photoelectron workflow with its model limits; removed `[outer_plasma]`, `[coupling]`, and legacy selectors are rejected input.

## Review Output

Follow the repository language policy and review in Japanese. Lead with correctness or navigation problems, cite file and line locations, and separate required fixes from optional polish. When editing files, summarize the reader-visible outcome and the checks performed.
