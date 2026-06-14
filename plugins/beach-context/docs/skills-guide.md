# BEACH plugin skill guide

BEACH 固有の利用者支援では、次の順で skill を選びます。

| Skill | 使う場面 |
| --- | --- |
| `beach-config-review` | `beach.toml`、schema validation、入力設定の整合性確認 |
| `beach-run-diagnose` | install/build/run が失敗した、出力がない、restart が壊れた、実行結果が異常 |
| `beach-case-design` | 新しい計算設定、parameter sweep、研究目的から設定値を決める |
| `beach-output-analysis` | `outputs/latest` の読み方、`beachx` 可視化、Python API、履歴解析 |
| `beach-simulator-guide` | BEACH の全体像、学習順、ドキュメント案内、ユーザー向け tutorial |
| `beach-method-summary` | 論文、発表、README、設計メモ向けの手法説明 |
| `beach-issue-report` | バグ、改善要望、機能追加、ドキュメント不足を GitHub Issue 草案にする |

repo 全体を読める開発 checkout では root docs と source を優先し、plugin の `references/` は repo 外利用者向け snapshot として扱います。

共通ルール:

- ユーザーの主言語で応答する。日本語が混じる場合は原則日本語。
- code identifier、TOML key、file path、command、CSV column は翻訳しない。
- 不明な既定値や未文書化挙動は推測で断定しない。
- Fortran 実装と `SPEC.md` を simulation behavior の source of truth とする。
- `tol_rel` は監視出力であり、現行実装では早期終了条件ではない。
