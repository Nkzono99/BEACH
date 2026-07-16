title: 10 分チュートリアル

Lang: [日本語](Tutorial.md) | [English](Tutorial.en.md)

# 10 分チュートリアル

このチュートリアルでは、1 個の電子を絶縁体平面へ向けて追跡し、吸収と表面電荷の堆積までを確認します。
公式入門ケースは非周期の最小構成です。FMM、周期境界、光電子、導体・誘電体モデルは使いません。

## 0. インストールを確認する

```bash
beach --version
beachx --help
```

両方のコマンドが実行できれば準備完了です。コマンドが見つからない場合は、先に
[インストール](Installation.html)を確認してください。

## 1. 設定を作る

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
```

`beachx config init` が作るのは、非周期の公式入門ケースです。生成内容はリポジトリの
[`examples/tutorial_insulator.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/tutorial_insulator.toml)
と同一です。設定全文を読む必要はありません。まずは次のキーだけ押さえてください。

| キー | 値 | このケースでの意味 |
| --- | --- | --- |
| `batch_count` | `1` | 1 バッチだけ実行する |
| `npcls_per_step` | `1` | 1 個の電子を生成する |
| `drift_velocity` | `[0.0, 0.0, -1.0e6]` | 電子を平面へ向ける |
| `surface_model` | `"insulator"` | 衝突した電子の電荷を表面に蓄積する |
| `field_solver` | `"direct"` | 直接計算で電場を求める |
| `field_bc_mode` | `"free"` | 周期境界を使わない |
| `dir` | `"outputs/latest"` | 結果の出力先 |

## 2. 実行して成功を確認する

```bash
beach beach.toml
beachx inspect outputs/latest
```

正常終了すると、`outputs/latest/summary.txt` と `outputs/latest/charges.csv` が生成されます。
この決定論的なケースでは、次の件数を確認できます。

```text
processed_particles=1
absorbed=1
batches=1
survived_max_step=0
```

これは、1 個の電子を 1 バッチで処理し、最大ステップ数に達する前に平面へ吸収できたことを表します。

## 3. 堆積した電荷を見る

```bash
beachx inspect outputs/latest --save-mesh outputs/latest/charge.png
```

`charges.csv` では、電子が衝突した要素だけに負電荷が堆積します。この結果は、粒子の生成から
衝突判定、吸収、表面電荷の堆積までの実行経路を確認するものです。定常帯電状態を表すものではありません。

![公式入門ケースの表面電荷密度](media/tutorial_insulator_charge.png)

## 4. 目的に合わせて次へ進む

| 目的 | 最初に試すこと | 関連ページ |
| --- | --- | --- |
| 粒子数、発射位置、速度、表面解像度を変える | `npcls_per_step`、`pos_low` / `pos_high`、`drift_velocity`、`nx` / `ny` を編集する | [設定レシピ](ConfigurationRecipes.html) |
| 出力ファイルを読み、履歴や分布を調べる | `summary.txt`、`charges.csv`、履歴ファイルの役割を確認する | [出力ガイド](OutputGuide.html)、[後処理チュートリアル](PostprocessTutorial.html) |
| 粒子更新、衝突、表面電荷更新の仕組みを理解する | 1 バッチ内の計算順序を追う | [アルゴリズム](Algorithms.html) |
| 結果を研究へ使う | 時間刻み、メッシュ、粒子数への依存性を確認する | [計算結果の妥当性確認](ValidationGuide.html) |
