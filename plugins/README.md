Lang: [English](README.en.md) | [日本語](README.md)

# BEACH Codex Plugins

このディレクトリには、BEACH の利用・解析・保守に必要な文脈を Codex に配布するための plugin を置きます。

## 利用できる plugin

| Plugin | 内容 |
| --- | --- |
| [beach-context](beach-context/README.md) | 設定レビュー、実行診断、ケース設計、出力解析、手法説明、学習ガイド、feedback 整理用の skill と同梱 reference |
| [beach-context-claude](beach-context-claude/README.md) | Claude Code 向けの BEACH agent / slash command 風 prompt / 参照マップ |

## 導入

GitHub から marketplace と plugin だけを sparse install する場合:

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

install 後に Codex を再起動すると、repo 外の作業ディレクトリでも `beach-context` の skill が利用できます。

登録済み marketplace を更新する場合:

```bash
codex plugin marketplace upgrade beach
```

ローカル checkout を marketplace として使う場合:

```bash
codex plugin marketplace add /path/to/BEACH
```

この場合も、登録後に `/plugins` から `BEACH Context` を install します。

## Claude Code 版

[beach-context-claude](beach-context-claude/README.md) は Codex plugin ではなく、Claude Code の
agent / command 定義として使うための context package です。BEACH repo 内で参照するか、必要な
`agents/*.md` と `commands/*.md` だけを Claude Code の project-local 設定に配置してください。

## Skill の見え方

repo root の `AGENTS.md` は BEACH 開発者向けの project-local 指示です。BEACH repo 配下で Codex を起動した場合だけ読み込まれます。

一方、`plugins/beach-context/skills/` は利用者向け plugin skill です。`/plugins` で install して Codex を再起動すると、`~` など repo 外で Codex を起動しても利用できます。

Claude Code では [beach-context-claude](beach-context-claude/README.md) の `agents/` と `commands/` を
参照してください。必要なものだけを `.claude/agents/` または `.claude/commands/` に配置して使う想定です。

## runops 連携

`beach-bem` は [runops](https://github.com/Nkzono99/runops) 2.x 向けの `beach` adapter を
`runops.adapters` entry point（`beach.tooling.runops_adapter:Adapter`）で公開します。runops project の
`.venv` に `beach-bem` が入っていれば、`runo init beach` と `runo update` が adapter を見つけ、
`plugins/beach-context` の skill と reference（wheel では `beach.agent_context`）を project の
`.agents/skills/` と `.claude/skills/` へ展開します。plugin の install や再起動は不要です。

adapter は run ごとに `input/beach.toml` へ `--set` の dotted key（例: `run.batch_count`）を適用し、
出力先を `work/latest` にします。schema に無い key は作成時に warning を出します。
`summary.txt` の `batches` が `run.batch_count` に達した run を completed とします。
`runo new --from RUN --continue --set run.batch_count=<累積の目標>` は元 run の checkpoint から再開します。
`runo tidy --purge-work` は最終 checkpoint を残すので、整理した run からも再開できます。

## 配置方針

BEACH 固有の物理、設定仕様、出力仕様、既知の失敗モード、学習導線は simulator repo 内の plugin に置きます。Fortran/fpm/Slurm/Codex 運用の横断的な手順は、repo root の `AGENTS.md` や共通 plugin に寄せます。
