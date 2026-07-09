title: BEACH ドキュメント
ordered_subpage: OutputGuide.md
ordered_subpage: OutputGuide.en.md
ordered_subpage: ConfigurationRecipes.md
ordered_subpage: ConfigurationRecipes.en.md
ordered_subpage: Parameters.md
ordered_subpage: Parameters.en.md
ordered_subpage: Configuration.md
ordered_subpage: Configuration.en.md
ordered_subpage: PostprocessTutorial.md
ordered_subpage: PostprocessTutorial.en.md
ordered_subpage: PythonPostprocessAPI.md
ordered_subpage: PythonPostprocessAPI.en.md
ordered_subpage: Algorithms.md
ordered_subpage: Algorithms.en.md
ordered_subpage: FieldSolvers.md
ordered_subpage: FieldSolvers.en.md
ordered_subpage: ParticleChargeLoop.md
ordered_subpage: ParticleChargeLoop.en.md
ordered_subpage: FMMCore.md
ordered_subpage: FMMCore.en.md
ordered_subpage: BatchDurationStability.md
ordered_subpage: BatchDurationStability.en.md
ordered_subpage: Workflow.md
ordered_subpage: Workflow.en.md
ordered_subpage: FortranDependencyMap.md
ordered_subpage: FortranDependencyMap.en.md
ordered_subpage: agent-user-guide.md
ordered_subpage: agent-user-guide.en.md

Lang: [日本語](index.md) | [English](index.en.md)

# BEACH ドキュメント

BEACH (BEM + Accumulated CHarge) は、三角形境界要素上の表面電荷が作る Coulomb 場と、テスト粒子追跡を batch 単位で結合する表面帯電シミュレータです。
現行版の主対象は、絶縁体表面への電荷蓄積です。Fortran 実行系 `beach` が計算を行い、Python CLI / API の `beachx` と `beach` パッケージが設定検査、後処理、可視化を担当します。

## 3分クイックスタート

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem

mkdir run_periodic2
cd run_periodic2

beachx config init
beachx lint beach.toml
beachx config render beach.toml --output beach.rendered.toml

beach beach.rendered.toml
beachx inspect outputs/latest
```

`beach` が読むのは最終キーに展開済みの TOML です。高水準記法を使う場合は、実行前に `beachx config render` を実行してください。
出力先を指定しない `render` は入力ファイルを上書きするため、初回は `--output beach.rendered.toml` を推奨します。

## 成功確認

通常完了では、`outputs/latest/summary.txt` が生成され、`batches` が設定した `sim.batch_count` に到達します。
次に `beachx inspect outputs/latest` で吸収数、脱出数、要素電荷、mesh 情報を確認してください。
出力ファイルの意味は [出力の読み方](OutputGuide.html) にまとめています。

## 次に読むページ

| 目的 | ページ |
| --- | --- |
| まず動かした後に出力を確認したい | [出力の読み方](OutputGuide.html) |
| 設定例から自分のケースを作りたい | [設定レシピ](ConfigurationRecipes.html) |
| `beach.toml` の全キーを調べたい | [入力パラメータリファレンス](Parameters.html) |
| `beachx config render` の挙動を知りたい | [beachx config / 高水準記法ガイド](Configuration.html) |
| 図を作る最短手順を知りたい | [後処理チュートリアル](PostprocessTutorial.html) |
| Python API を網羅的に確認したい | [Python 後処理 API リファレンス](PythonPostprocessAPI.html) |
| 数値モデルを確認したい | [BEACH アルゴリズム概要](Algorithms.html) |

## ユースケース

| やりたいこと | 入口 |
| --- | --- |
| 小さな template mesh で動作確認する | [設定レシピ](ConfigurationRecipes.html) の「平面メッシュで最小実行」 |
| 粒子注入方式を選ぶ | [設定レシピ](ConfigurationRecipes.html) と [入力パラメータリファレンス](Parameters.html) の `particles` |
| 2軸周期境界を使う | [設定レシピ](ConfigurationRecipes.html) の `periodic2` と [場ソルバーと境界条件](FieldSolvers.html) |
| `batch_duration` を調整する | [`batch_duration` の安定性](BatchDurationStability.html) |
| 実装 API を追う | [Fortran API](https://nkzono99.github.io/BEACH/fortran/) と [Fortran 依存関係マップ](FortranDependencyMap.html) |

## ドキュメント一覧

- [出力の読み方](OutputGuide.html)
- [設定レシピ](ConfigurationRecipes.html)
- [Fortran パラメータファイル仕様](Parameters.html)
- [beachx config / 高水準記法ガイド](Configuration.html)
- [後処理チュートリアル](PostprocessTutorial.html)
- [Python 後処理 API リファレンス](PythonPostprocessAPI.html)
- [BEACH アルゴリズム概要](Algorithms.html)
- [場ソルバーと境界条件](FieldSolvers.html)
- [粒子追跡と表面電荷蓄積](ParticleChargeLoop.html)
- [Coulomb FMM コア詳細](FMMCore.html)
- [`batch_duration` の安定性](BatchDurationStability.html)
- [実行・開発ワークフロー](Workflow.html)
- [Fortran 依存関係マップ](FortranDependencyMap.html)
- [BEACH Agent User Guide](agent-user-guide.html)

API 単位の詳細は、FORD が生成する [Fortran API](https://nkzono99.github.io/BEACH/fortran/) から辿れます。
