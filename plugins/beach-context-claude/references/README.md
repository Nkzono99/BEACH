# BEACH Claude Context References

このdirectoryは、互換用project-local agentが参照する一次情報の地図です。正式Claude Code pluginは
隣接する`plugins/beach-context/`の同梱`references/`を使います。

## 優先参照

| 用途 | repo内の参照先 | 正式plugin単体の参照先 |
| --- | --- | --- |
| quick start / install | `../../../README.md` | `../../beach-context/references/README.md` |
| simulation behavior | `../../../SPEC.md` | `../../beach-context/references/SPEC.md` |
| beginner config | `../../../docs/Tutorial.md` | `../../beach-context/references/examples/tutorial_insulator.toml` |
| parameter reference | `../../../docs/Parameters.md` | `../../beach-context/references/fortran_parameter_file.md` |
| configuration workflow | `../../../docs/Configuration.md` | `../../beach-context/references/config_workflow.md` |
| execution / developer workflow | `../../../docs/Execution.md`, `../../../docs/Workflow.md` | `../../beach-context/references/fortran_workflow.md` |
| algorithm overview | `../../../docs/Algorithms.md` | `../../beach-context/references/SPEC.md` |
| field solver / FMM | `../../../docs/FieldSolvers.md`, `../../../docs/FMM.md`, `../../../docs/FMMCore.md` | `../../beach-context/references/fortran_fmm_core.md` |
| periodic2 | `../../../docs/PeriodicElectrostatics.md`, `../../../docs/FinitePeriodicConfiguration.md` | `../../beach-context/references/periodic_zero_mode_outer_plasma.md` |
| output | `../../../docs/OutputGuide.md` | `../../beach-context/references/SPEC.md` |
| Python post-processing | `../../../docs/PythonPostprocessAPI.md` | `../../beach-context/references/python_postprocess_api.md` |
| schema | `../../../schemas/beach.schema.json` | `../../beach-context/references/schemas/beach.schema.json` |

## 注意

- `docs/Algorithms.md`は全体像と索引です。FMM、periodic2、batch durationの詳細は専用文書で確認します。
- repoを読める場合はrootの現行docsとFortran実装を優先します。
- KUDPC login nodeではbuild/test/simulation/大規模解析を直接実行しません。
