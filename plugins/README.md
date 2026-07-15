Lang: [English](README.en.md) | [日本語](README.md)

# BEACH Plugins for Codex and Claude Code

このdirectoryには、BEACHの利用・解析・保守に必要な文脈をCodexとClaude Codeへ配布するpluginを置きます。

## 利用できるplugin

| Plugin | 内容 |
| --- | --- |
| [beach-context](beach-context/README.md) | Codex / Claude Code共通の正式plugin。設定レビュー、実行診断、case設計、出力解析、手法説明、学習ガイド、feedback整理用skill・agent・reference |
| [beach-context-claude](beach-context-claude/README.md) | 旧来のClaude Code project-local agent / commandを手動配置するための互換package |

## Codexへの導入

GitHubからmarketplaceとpluginだけをsparse installする場合:

```bash
codex plugin marketplace add Nkzono99/BEACH \
  --ref main \
  --sparse .agents/plugins \
  --sparse plugins/beach-context
codex plugin add beach-context@beach
```

install後は新しいCodex threadで利用してください。repo外の作業directoryでもskillが利用できます。

登録済みmarketplaceを更新する場合:

```bash
codex plugin marketplace upgrade beach
codex plugin add beach-context@beach
```

ローカルcheckoutをmarketplaceとして使う場合:

```bash
codex plugin marketplace add /path/to/BEACH
codex plugin add beach-context@beach
```

## Claude Codeへの導入

正式Claude Code pluginも同じ[beach-context](beach-context/README.md)を使います。

```bash
claude plugin marketplace add Nkzono99/BEACH
claude plugin install beach-context@beach-claude
```

checkoutを直接試す場合:

```bash
claude --plugin-dir ./plugins/beach-context
```

[beach-context-claude](beach-context-claude/README.md)は、既存のproject-local agent / commandを
手動配置するworkflow向けに残した互換packageです。

## 指示とskillの見え方

repo rootの`AGENTS.md`と`CLAUDE.md`はBEACH開発者向けのproject-local指示です。
正式pluginの`skills/`と`references/`は両製品で共有し、Claude Codeは`agents/`も読み込みます。

## 配置方針

BEACH固有の物理、設定仕様、出力仕様、既知の失敗mode、学習導線はsimulator repo内のpluginに置きます。
Fortran/fpm/Slurm運用の横断的な手順はrepo rootの`AGENTS.md`や共通pluginに寄せます。
