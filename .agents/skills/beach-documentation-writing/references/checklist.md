# 完了前の確認

## 内容

- [ ] 書いた挙動を実装（Fortran）または `SPEC.md` で確かめた。設定キーは schema、出力名は
      `schemas/beach.output-manifest.json` と書き出し処理で確かめた。
- [ ] 同じ説明が正本以外のページに残っていない（`grep` で主要な用語・キーを検索する）。
- [ ] 設定例が新形式のキーだけで書かれている。旧形式のキーは `GroupedConfiguration.md` にだけある。
- [ ] 削除済み・非推奨の入力への言及が、移行ページ以外にない。
- [ ] 実行の完了、数値的な収束、物理的な妥当性を区別して書いた。
- [ ] 日本語版と英語版の見出し構成、設定、コマンド、警告、リンク先が一致する。

## 構成

- [ ] ページを足した・消した・改題した場合、`docs-site/navigation.json` と
      [page-map.md](page-map.md) を更新した。
- [ ] 消したページへのリンクが残っていない（`grep -rn "OldPage.html" docs README*.md`）。
- [ ] plugin のコピーを正本と一致させた（page-map.md の「一緒に更新するもの」）。

## テスト

KUDPC ではログインノードで実行せず、`tssrun` で計算ノードに投入する。

```bash
python -m pytest -q tests/python/test_docs_sync.py \
  tests/python/test_documentation_contracts.py \
  tests/python/test_far_correction_contract.py
python tools/sync_starlight_docs.py --check   # docs-site を生成済みの場合
```

テストがドキュメントの古い文言を固定している場合は、文言ではなく契約（キーの網羅、リンクの解決、
コピーの一致）を検査するようにテストを直す。
