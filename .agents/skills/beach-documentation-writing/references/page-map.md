# ページの地図

各ページが答える問いと、正本として持つ内容。新しい内容はまずこの表で書く場所を決める。
ページを足す・統合するときは、この表と `docs-site/navigation.json` を同じ変更で更新する。

型: 入口 / チュートリアル / 手順 / モデル / リファレンス / 内部実装 / 開発（[page-types.md](page-types.md)）

## はじめる

| ファイル | slug | 型 | 答える問い | 正本として持つもの |
|---|---|---|---|---|
| `index.md` | `index` | 入口 | BEACH は何を計算し、何を計算しないか。どこから読むか | 適用範囲の表、目的別の読み始め方 |
| `Installation.md` | `installation` | 手順 | 使える状態にするには | 動作要件、インストール・更新 |
| `Tutorial.md` | `tutorial` | チュートリアル | 1 ケースを動かし、バッチ間の帰還を見るには | 公式入門ケースと参照出力 |

## ケースを作る

| ファイル | slug | 型 | 答える問い | 正本として持つもの |
|---|---|---|---|---|
| `ConfigurationRecipes.md` | `configuration-recipes` | 手順 | 研究ケースを何の順で決めるか | 判断の順序、組み込み形状の選び方 |
| `Configuration.md` | `configuration` | 手順 | 設定ファイルを作り、検査するには | 7 グループの構成、`beachx config init/lint/validate/diff`、`beach --check-config` |
| `ParticleSourcesBoundaries.md` | `particle-sources-boundaries` | 手順 | 粒子をどこから入れるか | 粒子源の選び方（`volume_seed`、`plane_source`、境界流入、`photo_raycast`） |
| `ReservoirInjection.md` | `reservoir-injection` | 手順 | 外部プラズマの分布を境界から入れるには | 流入分布（Maxwell、速度 grid）、流入の写像（`source_vdf`、`infinity_barrier`） |
| `PeriodicPlasmaSurface.md` | `periodic-plasma-surface` | 手順 | 太陽風と光電子の中の周期表面で、上端をどう閉じるか | 4 方式（開放、閉じた光電子、スカラー障壁、外部シース）の選び方と設定の差分 |
| `FieldSolvers.md` | `field-solvers` | 手順 | 場ソルバと場境界をどう選ぶか | ソルバと場境界の対応表 |
| `BatchDurationStability.md` | `batch-duration-stability` | 手順 | バッチ幅をどう決めるか | 固定幅の比較手順、適応進行の使い方 |
| `Execution.md` | `execution` | 手順 | 実行・並列化・再開するには | 実行コマンド、checkpoint からの再開手順 |

## 結果を調べる

| ファイル | slug | 型 | 答える問い | 正本として持つもの |
|---|---|---|---|---|
| `OutputGuide.md` | `output-guide` | 手順 | 実行は正常に終わったか。最初にどのファイルを見るか | 最初に見るファイル、入門ケースでの照合 |
| `PostprocessTutorial.md` | `postprocess-tutorial` | チュートリアル | 結果を図にするには | `beachx inspect/animate` と Python API の最初の使い方 |
| `ObjectForcesDetachment.md` | `object-forces-detachment` | 手順 | 物体に働く力と離脱を調べるには | 力・離脱の解析手順 |
| `ValidationGuide.md` | `validation-guide` | 手順 | 結果を研究に使ってよいか | 収支・収束・モデル診断の確認手順、モデル別の診断表 |
| `Troubleshooting.md` | `troubleshooting` | 手順 | 失敗の原因は何か | 症状から原因への索引 |

## モデルを理解する

| ファイル | slug | 型 | 答える問い | 正本として持つもの |
|---|---|---|---|---|
| `Algorithms.md` | `algorithms` | モデル | 表面電荷・場・粒子をどの順に更新するか | バッチの 6 段階、`dt` / バッチ幅 / バッチ数の区別 |
| `SurfaceModels.md` | `surface-models` | モデル | 吸収した電荷を表面でどう保持するか | 絶縁体、浮遊導体 |
| `PhotoelectronEmission.md` | `photoelectron-emission` | モデル | 光電子をどう放出し、その電荷をどう閉じるか | 光線追跡による放出、反作用電荷、`neutral_return`、`fixed_current` |
| `ParticleEscapeReturn.md` | `particle-escape-return` | モデル | box の面に達した粒子をどう扱うか | `open`（脱出）、`reflect`、`redistributed_reflect`、`potential_barrier` |
| `ZhaoStationaryClosure.md` | `zhao-stationary-closure` | モデル | 外部シースの定常解を上端へどう接続するか | 零電流根、上端への写像、固定電流、外部根の更新、既知の制約 |
| `PeriodicElectrostatics.md` | `periodic-electrostatics` | モデル | x/y 周期の場をどう組み立てるか | 3 成分（有限画像、面内変動、面平均）、下側境界、有限画像の意味 |

## リファレンス

| ファイル | slug | 型 | 正本として持つもの |
|---|---|---|---|
| `Parameters.md` | `parameters` | リファレンス | 新形式の全キー（型、単位、既定値、制約） |
| `GroupedConfiguration.md` | `configuration-groups` | リファレンス | 旧形式からの移行: 旧キーと新キーの対応表、非推奨・削除済みの入力とサポート期間 |
| `OutputReference.md` | `output-reference` | リファレンス | 全出力ファイル、列、`summary.txt` の key、checkpoint の契約 |
| `Glossary.md` | `glossary` | リファレンス | 用語（日本語、英語、識別子） |
| `PythonPostprocessAPI.md` | `python-postprocess-api` | リファレンス | Python API |
| `SurfaceChargeNumerics.md` | `surface-charge-numerics` | 内部実装 | 確定反映の順序、浮遊導体の連立方程式、並列集約 |
| `BatchDurationTheory.md` | `batch-duration-theory` | 内部実装 | バッチ幅依存性と適応制御の理論 |
| `DirectSolver.md` | `direct-solver` | 内部実装 | 三角形要素の直接評価 |
| `Treecode.md` | `treecode` | 内部実装 | Treecode |
| `FMM.md` | `fmm` | 内部実装 | FMM の設定、精度と性能の測り方 |
| `ParticleTrackingCollision.md` | `particle-tracking-collision` | 内部実装 | 粒子 1 step の流れ |
| `BorisPusher.md` | `boris-pusher` | 内部実装 | Boris 法 |
| `ParticleEvents.md` | `particle-events` | 内部実装 | 衝突・境界イベントの判定順、同時イベント |
| `PeriodicFarCorrection.md` | `periodic-far-correction` | 内部実装 | periodic2 の遠方補正（Ewald、演算子、cache）とその精度 |

## 開発者向け

| ファイル | slug | 型 | 正本として持つもの |
|---|---|---|---|
| `Architecture.md` | `architecture` | 開発 | 実行の流れ、状態の持ち主、機能から実装・テストへの対応（Code reference はここに集約） |
| `Workflow.md` | `workflow` | 開発 | 開発環境、変更ごとのテスト |
| `PhysicsReleaseVerification.md` | `physics-release-verification` | 開発 | 物理リリースの判定 |
| `FMMCore.md` | `fmm-core` | 内部実装 | FMM core の数式と内部 API |
| `FortranDependencyMap.md` | `fortran-dependency-map` | 開発 | 生成物（`make` で再生成。手で編集しない） |
| `agent-user-guide.md` | `agent-user-guide` | 入口 | AI エージェント向けの入口（サイドバー非表示） |

## 変更から直すページを引く

| 変更 | 正本 | あわせて確認するもの |
|---|---|---|
| 設定キーの追加・変更 | `Parameters.md` | 挙動を説明するモデルのページ、旧キーがあれば `GroupedConfiguration.md` の対応表、schema |
| 出力ファイル・列・`summary.txt` の key | `OutputReference.md` | `schemas/beach.output-manifest.json`、最初に見るべきものなら `OutputGuide.md` |
| 物理モデルの追加・変更 | そのモデルのページ | `Parameters.md`、`ValidationGuide.md` のモデル別診断表、関係する手順のページ |
| 周期表面の上端の扱い（光電子、障壁、外部シース） | 該当するモデルのページ | `PeriodicPlasmaSurface.md` の比較表と設定差分 |
| 粒子源・境界流入 | `ParticleSourcesBoundaries.md` / `ReservoirInjection.md` | `Parameters.md` |
| 場ソルバ・周期場 | `FieldSolvers.md` / `PeriodicElectrostatics.md` | 内部実装のページ |
| 実行・再開の手順 | `Execution.md` | `OutputReference.md` の checkpoint、`Troubleshooting.md` |
| 新しい失敗の形 | `Troubleshooting.md` | 原因を説明する正本 |
| 入力の削除・非推奨化 | `GroupedConfiguration.md` | CHANGELOG。他のページからは記述を消す |
| ソースの配置・責務 | `Architecture.md` | `Workflow.md` のテスト対応 |
| 用語の追加・変更 | `Glossary.md` | 使っている全ページ |

## 一緒に更新するもの

- `docs-site/navigation.json`: ページの追加・削除・改題。`docs/*.md` は全て載せる（テストで検査）。
- plugin のコピー（バイト一致をテストで検査）:
  `plugins/beach-context/references/fortran_parameter_file.md` ← `docs/Parameters.md`、
  `fortran_fmm_core.md` ← `docs/FMMCore.md`、`agent-user-guide.md` ← `docs/agent-user-guide.md`、
  `python_postprocess_api.md` ← `docs/PythonPostprocessAPI.md`、`SPEC.md` ← `SPEC.md`。
  `GroupedConfiguration*.md` も同じ内容に保つ。
- `README.md` / `README.en.md`: 入口から外れたページへのリンク。
