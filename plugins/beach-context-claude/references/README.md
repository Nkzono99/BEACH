# BEACH Claude Context References

この directory は、Claude Code 版 context package が参照すべき一次情報の地図です。
現行の正式 document は repo root の docs にあります。

## 優先参照

| 用途 | 参照先 |
| --- | --- |
| quick start / install / command | `../../../README.md` |
| simulation behavior の source of truth | `../../../SPEC.md` |
| algorithm / FMM / periodic2 / batch duration | `../../../docs/Algorithms.md` |
| parameter reference | `../../../docs/Parameters.md` |
| configuration workflow | `../../../docs/Configuration.md` |
| developer workflow | `../../../docs/Workflow.md` |
| Python post-processing API | `../../../docs/PythonPostprocessAPI.md` |
| user guide | `../../../docs/agent-user-guide.md` |
| schema | `../../../schemas/beach.schema.json` |
| examples | `../../../examples/beach.toml`, `../../../examples/periodic2_basic/beach.toml` |

## 注意

- FMM と batch duration stability の詳細は `docs/Algorithms.md` に統合済みです。
- 正式 docs は PascalCase 名に統一されています。
- `plugins/beach-context/` の `references/` は Codex plugin 用 snapshot です。Claude Code 版では、可能な限り
  repo root の現行 docs を優先してください。
- KUDPC login node 上では、build/test/simulation/大規模解析を直接実行しないでください。
