title: BEACH ドキュメント
ordered_subpage: Installation.md
ordered_subpage: Tutorial.md
ordered_subpage: OutputGuide.md
ordered_subpage: Troubleshooting.md
ordered_subpage: ConfigurationRecipes.md
ordered_subpage: Configuration.md
ordered_subpage: PostprocessTutorial.md
ordered_subpage: ValidationGuide.md
ordered_subpage: Parameters.md
ordered_subpage: PythonPostprocessAPI.md
ordered_subpage: Algorithms.md
ordered_subpage: FieldSolvers.md
ordered_subpage: ParticleChargeLoop.md
ordered_subpage: FMMCore.md
ordered_subpage: BatchDurationStability.md
ordered_subpage: PhysicsReleaseVerification.md
ordered_subpage: Workflow.md
ordered_subpage: FortranDependencyMap.md
ordered_subpage: agent-user-guide.md

Lang: [日本語](index.md) | [English](index.en.md)

# BEACH

BEACHは、絶縁体表面に蓄積する電荷と、その電荷が作る電場中のテスト粒子軌道をbatch単位で
計算するシミュレータです。

## 初めて使う方

1. [動作環境を確認してインストールする](Installation.html)
2. [公式入門ケースを実行する](Tutorial.html)
3. [出力を確認する](OutputGuide.html)
4. [物理・数値的な妥当性を確認する](ValidationGuide.html)

環境構築済みなら、`beachx config init`から始められます。

```bash
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

## コマンドとデータの流れ

```text
beach.toml
    ├─ beachx lint ── 設定検査
    ▼
beach ────────────── Fortranシミュレーション
    ▼
outputs/<case>/
    ├─ beachx inspect / animate ── 確認・可視化
    └─ Python package beach ───── 独自解析
```

## 対応範囲

| 機能 | 状況 | 注記 |
| --- | --- | --- |
| 絶縁体表面への電荷蓄積 | 対応 | 現行版の主対象 |
| 浮遊導体 | 条件付き | 利用可能な場境界とsurface modelの組合せを確認 |
| 誘電分極 | 未対応 | `epsilon_r`は現行ではmetadataで、独立した分極境界条件ではない |
| 2軸周期境界 | 対応 | finite image、cached infinite operator、zero-modeの意味を区別する |
| outer plasma | 条件付き | model applicabilityとerror contractを満たす場合のみ |
| 自動収束停止 | 未対応 | `tol_rel`はmonitoring metric |

非対応の組合せは別modelへsilent fallbackせず停止します。

## 目的別の入口

| 目的 | ページ |
| --- | --- |
| 自分のcaseを作る | [設定レシピ](ConfigurationRecipes.html) |
| 後処理と図を作る | [後処理チュートリアル](PostprocessTutorial.html) |
| 全入力keyを調べる | [入力パラメータ](Parameters.html) |
| 数値モデルを確認する | [アルゴリズム概要](Algorithms.html) |
| 問題を解決する | [トラブルシューティング](Troubleshooting.html) |

Fortran APIは[FORD API](https://nkzono99.github.io/BEACH/fortran/)から参照できます。
