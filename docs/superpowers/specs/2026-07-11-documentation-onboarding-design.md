# Documentation onboarding redesign

## Goal

初見利用者が、対応範囲を理解し、前提ツールを導入し、一つの公式ケースを実行し、正常終了と
物理・数値妥当性を区別して判断できる導線を作る。既存の技術リファレンスとURLは維持する。

## Chosen approach

大規模な文書分割ではなく、4つの入門ページとsidebar再編を行う。

- `Installation`: Python 3.10+、Fortran compiler、fpm、PATH、更新・削除
- `Tutorial`: 公式入門ケースの目的、完全config、実行、期待出力、次に変える値
- `ValidationGuide`: 正常終了と物理・数値妥当性を分離したchecklist
- `Troubleshooting`: installation、runtime、output、restart、unsupported modelの代表障害

日本語と英語を同じ情報構造で提供する。既存ページのslugは変更しない。

## Official beginner case

`beachx config init` と `examples/tutorial_insulator.toml` を同じ小ケースに固定する。
direct/free-space、単一 `volume_seed` electron、絶縁体平面、1 batch とし、FMM、periodic2、
photoelectron、conductor、dielectricを含めない。多機能例は設定レシピから発展例として残す。

Python test は既定configとexampleのsemantic一致を固定する。Fortran smoke testは計算ノードで
exampleを実行し、`summary.txt`、batch完了、吸収/escape内訳を確認する。

## Navigation

sidebar順序:

1. はじめに: overview、Installation、Tutorial、OutputGuide、Troubleshooting
2. 使い方: recipes、configuration、postprocess、ValidationGuide
3. リファレンス: Parameters、Python API、Physics release verification
4. 数値アルゴリズム
5. 開発者向け
6. AIエージェント・自動化向け

トップページは対応範囲表、command/data flow、初回4ステップ、目的別リンクだけに整理し、
重複する文書一覧を削除する。

## Scope and constraints

- v1.0主対象のinsulator accumulationを前面に出す。
- conductorは条件付き、dielectric polarizationは未実装、`tol_rel`はmonitoring-onlyと明記する。
- `beach` binary、`beachx` CLI、Python `beach` packageの役割を区別する。
- 新しい依存は追加しない。
- Starlight syncのcanonical sourceは引き続き`docs/`とする。

## Verification

- Python: default/example semantic契約、docs sync unit tests
- Static: Markdown link target、schema、`git diff --check`
- Site: `tools/sync_starlight_docs.py --check` と Starlight production build
- Fortran: official beginner caseのcompute-node smoke run

## Non-goals

- 既存の長いParameters/Algorithmsを全面分割しない。
- 新しいsimulation physicsを導入しない。
- 公開済みslugを削除・変更しない。
