Lang: [English](README.en.md) | [日本語](README.md)

# BEACH Codex Plugins

このディレクトリには、BEACH の利用・解析・保守に必要な文脈を Codex に配布するための plugin を置きます。

## 利用できる plugin

| Plugin | 内容 |
| --- | --- |
| [beach-context](beach-context/README.md) | 設定レビュー、実行診断、ケース設計、出力解析、手法説明、学習ガイド、issue 作成補助用の skill と同梱 reference |

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

## Skill の見え方

repo root の `AGENTS.md` は BEACH 開発者向けの project-local 指示です。BEACH repo 配下で Codex を起動した場合だけ読み込まれます。

一方、`plugins/beach-context/skills/` は利用者向け plugin skill です。`/plugins` で install して Codex を再起動すると、`~` など repo 外で Codex を起動しても利用できます。

## 配置方針

BEACH 固有の物理、設定仕様、出力仕様、既知の失敗モード、学習導線は simulator repo 内の plugin に置きます。Fortran/fpm/Slurm/Codex 運用の横断的な手順は、repo root の `AGENTS.md` や共通 plugin に寄せます。
