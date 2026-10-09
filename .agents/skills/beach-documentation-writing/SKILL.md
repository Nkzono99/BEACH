---
name: beach-documentation-writing
description: Write, restructure, or review BEACH documentation (docs/*.md and *.en.md, README, docs-site/navigation.json, plugin reference copies). Use whenever a code change alters user-visible behavior, configuration keys, outputs, or physics models, or when a page has become long, duplicated, or hard to read.
---

# BEACH ドキュメントの書き方

BEACH のドキュメントは、機能を足すたびに段落を継ぎ足すと、同じ説明が複数のページに散り、
削除した機能の経緯が本文に残って読めなくなる。このスキルは、**どのページに何を書くか**を先に決め、
**正本を書き換える**ことで、ページ数と長さを保ったまま内容を最新にするための手順である。

## 手順

1. **変更を分類し、直すページを決める。**
   [references/page-map.md](references/page-map.md) の「変更から直すページを引く」表で、変更の種類から
   正本のページを決める。正本以外のページは、要約（2 文以内）とリンクが古くなっていないかだけを確認する。
   どのページにも収まらない内容なら、新しいページを作る前に page-map.md の構成を見直す。
2. **正本の該当節を書き換える。** 末尾への追記ではなく、該当節を現在の挙動の説明に置き換える。
   - 経緯（以前はこうだった、何を削除した）は書かない。CHANGELOG と移行ページ（`GroupedConfiguration.md`）へ。
   - 設定例は新形式（7 グループ）のキーだけを使う。旧形式のキーは移行ページの対応表にだけ置く。
   - 他のページの内容を説明し直さない。必要なら 2 文以内で要約してリンクする。
3. **ページの型に合わせる。** [references/page-types.md](references/page-types.md) の骨組みに沿う。
   手順のページに内部実装を、モデルのページに全キーの一覧を入れない。
4. **文体と用語を揃える。** [references/style.md](references/style.md) に従う。用語は
   [`docs/Glossary.md`](../../../docs/Glossary.md) を正本とする。
5. **英語版を対訳で揃える。** 日本語版と同じ見出し構成、同じコマンド・設定・警告・リンク先にする。
6. **確認する。** [references/checklist.md](references/checklist.md) の項目を通す。

## 判断の原則

- **1 ページ 1 問い、1 概念 1 正本。** 同じ説明が 2 か所にあれば、正本でない方を要約とリンクに置き換える。
- **現行の挙動だけを書く。** 実装（Fortran）と `SPEC.md` を正とし、ドキュメントの記述と食い違えば実装に合わせて直す。
- **読む人が使う順に書く。** 通常の使い方を先に、選択肢・制約・内部実装を後に書く。
- **長さの上限を守る。** 手順とモデルのページは 200 行程度まで。超えたら、内容を正本のページへ移すか、
  page-map.md を更新してページを分ける。
- **結論の強さを区別する。** 実行の完了、数値的な収束、物理的な妥当性を別々に書く。

## レビューするとき

日本語でレビューする。読む人の誤解や操作の失敗につながる問題（誤った記述、壊れたリンク、正本の重複、
旧形式のキー）を先に、文体の改善を後に挙げ、ファイルと行を示す。
