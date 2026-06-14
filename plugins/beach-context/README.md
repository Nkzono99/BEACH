Lang: [English](README.en.md) | [日本語](README.md)

# BEACH Context Plugin

BEACH の設定ファイル、実行手順、出力解析、既知の失敗モード、学習導線、issue 作成補助を Codex に渡すための repo-local plugin です。各 skill はユーザーの言語に合わせて、日本語・英語などで応答する前提です。

利用者が `pip install beach-bem` だけで BEACH を入れており、repo 全体を読めない状況を前提にしています。そのため、この plugin は `references/` に `SPEC.md`、設定 workflow、`beach.toml` 仕様、Python 後処理 API、schema、代表 example などのスナップショットを同梱しています。repo 全体を読める開発環境では、同名の最新 docs を優先して確認してください。

## 導入

GitHub から marketplace と plugin だけを読む場合:

```bash
codex plugin marketplace add Nkzono99/BEACH \
  --ref main \
  --sparse .agents/plugins \
  --sparse plugins/beach-context
```

この時点では marketplace が登録されるだけで、`BEACH Context` plugin はまだ有効化されていません。Codex を起動し、`/plugins` から `BEACH Context` を install してください。

```bash
codex
# Codex 内で /plugins を開く
```

install 後に Codex を再起動すると、repo 外の作業ディレクトリでもこの plugin の skill が利用できます。

登録済み marketplace を更新する場合:

```bash
codex plugin marketplace upgrade beach
```

ローカル checkout を marketplace として使う場合:

```bash
codex plugin marketplace add /path/to/BEACH
```

この場合も、登録後に `/plugins` から `BEACH Context` を install します。

## 同梱 skill

- `beach-config-review`: `beach.toml` の設定レビューと検証方針
- `beach-run-diagnose`: install/build/実行失敗、異常終了、出力欠損、restart 問題の診断
- `beach-case-design`: 物理意図から BEACH 設定、parameter sweep を設計
- `beach-output-analysis`: `outputs/latest`、CSV、履歴、Python API、`beachx` による解析導線
- `beach-simulator-guide`: 利用者向けの学習・運用ガイド
- `beach-method-summary`: 論文・発表・README 用の手法説明
- `beach-issue-report`: バグ報告、改善要望、追加機能要望の GitHub Issue 化

各 skill の詳しい使い分けは [docs/skills-guide.md](docs/skills-guide.md) を参照してください。

## 同梱 reference

`references/` には、plugin 単体でも利用者支援ができるように以下を同梱しています。

- `README.md`, `SPEC.md`
- `agent-user-guide.md`, `fortran_workflow.md`
- `config_workflow.md`, `fortran_parameter_file.md`
- `python_postprocess_api.md`
- `fortran_fmm_core.md`, `batch_duration_stability.md`
- `schemas/beach.schema.json`
- `examples/beach.toml`, `examples/periodic2_basic/beach.toml`

設定ファイルの事前検査には `beachx config validate`、render 後の軽量確認には `beachx estimate-workload`、出力確認には `beachx inspect outputs/latest` を優先します。

## 配布方針

BEACH 固有の知識はこの plugin に置き、Fortran/fpm/Slurm/Codex 運用の横断的な手順は共通 plugin や repo root の `AGENTS.md` に寄せます。
