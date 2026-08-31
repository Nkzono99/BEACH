# BEACH documentation style guide

This guide adapts Japanese technical-writing and cognitive-rhythm principles to BEACH documentation. It is a decision guide, not a list of words to replace mechanically.

## Choose the page structure

### Entry pages and README

- State what BEACH computes, its current scope, and the result a new user can obtain.
- Keep one runnable path near the top and state its expected result.
- Link to canonical installation, configuration, validation, and development pages instead of reproducing them.
- Separate user tasks from contributor tasks.

### Tutorials and task guides

- Use the order: prerequisite, action, expected output, interpretation, next choices.
- Keep the main path uninterrupted. Move complete configurations and optional variants to links, tables, or collapsible details.
- Explain what a successful command proves and what it does not prove.
- Show only the smallest configuration difference needed for the task. Link to a complete file under `examples/`
  instead of duplicating it in prose.

### Reference pages

- Optimize headings, tables, identifiers, defaults, constraints, and cross-links for lookup.
- Put the definition or decision rule before background discussion.
- Avoid rhetorical questions, suspense, and narrative detours.

### Explanations and numerical-method pages

- Start from a real physical or numerical question, competing expectation, measurement, or tradeoff.
- Alternate concrete evidence with interpretation when that makes the reasoning easier to follow.
- Close questions that the page opens. Do not create tension by narrating the document itself.
- Put practical selection or validation criteria before long derivations when readers need them to act.
- If a page must teach both model meaning and implementation internals, split it. A reader should not fall from an
  overview into data structures, parallel reductions, rollback state, or private APIs without choosing an internals page.

## Keep human-facing pages at one abstraction level

- Introduce no more than three major new concepts in one H2 section.
- Keep overview diagrams to seven nodes or fewer.
- Use H3 as the usual maximum depth for entry, tutorial, task-guide, and explanation pages.
- If one prose paragraph contains five or more code identifiers, move the identifiers to a table, a minimal example,
  or the canonical reference.
- Keep important limits visible, but put the complete limit matrix in reference documentation.
- In Japanese prose, explain an implementation term in Japanese before giving its searchable identifier.

## Write Japanese technical prose

- Keep one topic per paragraph and state the paragraph's claim early.
- Prefer active constructions that identify the actor or calculation.
- Use Japanese terms in prose when established translations exist. Preserve exact identifiers in code formatting.
- Add a half-width space between Japanese text and Latin letters or digits, except next to Japanese punctuation.
- Keep commands, logs, configuration, and output in fenced code blocks.
- Remove empty previews and summaries such as `重要なのは`, `ここでは`, and `次に〜を見る` when the sentence adds no information about BEACH.
- Replace vague intensifiers and borrowed business terms with conditions, mechanisms, values, or named components.
- Do not shorten sentences at the cost of missing scope, assumptions, or uncertainty.

## Apply cognitive rhythm selectively

Use rhythm to clarify explanatory reasoning, not to decorate procedures. A useful transition comes from BEACH itself: an assumption fails under a boundary condition, an approximation changes error behavior, or a diagnostic separates two possible causes. Do not insert invented hesitation, rhetorical questions, or dramatic one-line conclusions into installation steps, command references, or parameter tables.

## Final checks

- The opening identifies the reader and the page's useful outcome.
- Commands have prerequisites and observable expected results.
- Terms keep one meaning across the page and its linked pages.
- Behavioral claims match `SPEC.md` and the current implementation.
- Links point to the canonical page and do not create competing instructions.
- Japanese and English counterparts contain the same commands, constraints, and warnings.
- The page distinguishes execution success, numerical convergence, and physical validity.

## Influences

- k16shikano, cognitive-rhythm-writing: https://gist.github.com/k16shikano/eb2929f13ed19c97188393d297be8432
- hikimay, japanese-tech-writing: https://github.com/hikimay/japanese-tech-writing
