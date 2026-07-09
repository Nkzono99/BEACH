# BEACH Claude Code Context Package

BEACH を Claude Code で扱うための repo-local context package です。Codex plugin
`plugins/beach-context/` の skill を、Claude Code の subagent / slash command
として使いやすい形に移植しています。

この package は BEACH repo と一緒に使う前提です。正式な仕様確認では repo root の
`README.md`, `SPEC.md`, `docs/Algorithms.md`, `docs/Parameters.md` を優先してください。

## 導入

BEACH repo 内で Claude Code を起動する場合は、この directory を参照しながら使えます。

```bash
cd /path/to/BEACH
claude
```

Claude Code の project-local agent として常用する場合は、必要な agent 定義だけを
`.claude/agents/` に配置してください。既存の `.claude/agents/` には開発者向け agent
があるため、この package 側の agent は利用者支援・設定診断・出力解析に寄せています。

## 含まれる agent

| Agent | 用途 |
| --- | --- |
| `beach-simulator-guide` | BEACH の概念、学習順、公式 docs への案内 |
| `beach-config-review` | `beach.toml` の schema、物理設定、mesh、solver、output のレビュー |
| `beach-run-diagnose` | install/build/run 失敗、異常終了、出力欠損、restart 問題の診断 |
| `beach-output-analysis` | `outputs/latest`、CSV、history、Python API、`beachx` による解析支援 |
| `beach-method-summary` | 論文・発表・README 向けの手法説明、アルゴリズム要約 |
| `feedback-beach` | バグ報告、改善要望、ドキュメント要望の feedback 整理と GitHub Issue 草案作成 |

各 agent は日本語を既定とし、TOML key、command、file path、CSV column は翻訳しません。

## 含まれる command

`commands/` には Claude Code の custom slash command 風に使える短い prompt を置いています。
必要に応じて `.claude/commands/` へ配置してください。

| Command | 用途 |
| --- | --- |
| `beach-config-review.md` | 指定した `beach.toml` を設定レビューする |
| `beach-run-diagnose.md` | 失敗ログや出力 directory から原因を切り分ける |
| `beach-output-analysis.md` | 出力の読み方・解析手順を組む |
| `beach-method-summary.md` | 手法説明文やアルゴリズム要約を作る |
| `feedback-beach.md` | feedback を Issue 草案に整理する |
| `beach-simulator-guide.md` | 利用者向けの学習導線を作る |

## 参照 document

参照マップは `references/README.md` にあります。現行 docs は PascalCase 名に統一されています。
FMM core と batch duration stability の詳細は `docs/Algorithms.md` に統合済みです。

## 運用方針

- Fortran 実装と `SPEC.md` を simulation behavior の source of truth とします。
- 設定値の網羅確認は `docs/Parameters.md` と `schemas/beach.schema.json` を見ます。
- アルゴリズム説明は `docs/Algorithms.md` を起点にします。
- `plugins/beach-context/` は Codex plugin として残し、この Claude 版からは破壊的に変更しません。
- KUDPC login node 上では、ビルド、テスト、simulation 実行、大きな出力解析を直接実行しません。
