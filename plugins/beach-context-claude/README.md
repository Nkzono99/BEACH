# BEACH Claude Code Compatibility Package

このdirectoryは、旧来のClaude Code project-local agent / commandを手動配置するworkflowとの互換用です。
新規導入では、Codexと共通の正式plugin `plugins/beach-context/`を使ってください。正式pluginは
`.claude-plugin/plugin.json`、repo rootの`.claude-plugin/marketplace.json`、自己完結したskills・agents・
referencesを備えています。

## 推奨: 正式pluginをinstallする

```bash
claude plugin marketplace add Nkzono99/BEACH
claude plugin install beach-context@beach-claude
```

更新:

```bash
claude plugin marketplace update beach-claude
claude plugin update beach-context@beach-claude
```

checkoutをinstallせず検証:

```bash
claude plugin validate ./plugins/beach-context
claude --plugin-dir ./plugins/beach-context
```

正式pluginはversionをmanifestへ固定せず、Git sourceのcommit SHAを更新単位として使います。

## 互換用の手動配置

既存workflowを維持する必要がある場合だけ、`agents/*.md`を`.claude/agents/`、
`commands/*.md`を`.claude/commands/`へ配置してください。`commands/`はClaude Codeで引き続き
読み込めますが、新しい正式pluginでは`skills/`を優先します。

## 含まれるagent / command

| 名前 | 用途 |
| --- | --- |
| `beach-simulator-guide` | BEACHの概念、学習順、公式docsへの案内 |
| `beach-config-review` | `beach.toml`のschema、物理設定、mesh、solver、outputのレビュー |
| `beach-case-design` | 物理目的からconfigとparameter sweepを設計 |
| `beach-run-diagnose` | install/build/run失敗、異常終了、出力欠損、restart問題の診断 |
| `beach-output-analysis` | output、CSV、history、Python API、`beachx`による解析支援 |
| `beach-method-summary` | 論文・発表・README向けの手法説明、algorithm要約 |
| `feedback-beach` | feedback整理とGitHub Issue草案作成 |

## 参照方針

- 正式plugin内では`plugins/beach-context/references/`の同梱snapshotを使えます。
- BEACH repo内ではrootの`README.md`, `SPEC.md`, `docs/`, `schemas/`, `src/`, `app/`を優先します。
- algorithm全体の入口は`docs/Algorithms.md`、詳細は`docs/FMM.md`, `docs/FMMCore.md`,
  `docs/PeriodicElectrostatics.md`, `docs/ParticleEvents.md`などの専用文書です。
- 設定生成の現行既定は`examples/tutorial_insulator.toml`と一致し、x/y periodic、
  `field_solver="fmm"`、`field_periodic_image_layers=1`、far correctionなしの$3\times3$有限画像和です。
- KUDPC login nodeではbuild、test、simulation、大きなoutput解析を直接実行しません。
