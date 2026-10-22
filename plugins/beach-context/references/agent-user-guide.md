title: BEACH Agent User Guide

# BEACH Agent User Guide

AI エージェントが BEACH の設定、実行、検証を始めるための運用ガイドです。
パラメータや API の網羅的な説明は再掲せず、それぞれの正本へ案内します。

## 最初に確認すること

BEACH は、境界要素法による電場評価と粒子追跡を組み合わせ、主に絶縁体表面の
電荷蓄積を計算します。現行の基準構成は次のとおりです。

- 外部領域のプラズマ場や粒子輸送は解かない
- 周囲粒子は `[particles.species.inflow]` で非周期box面から注入する
- 光電子は `photo_raycast`、`source.inject_face` と同じ species 面の反射、`neutral_return` で閉じる
- mesh hit は吸収し、表面電荷差分を batch ごとに commit する
- `run.batch_count` まで実行し、`particles.tracking.max_steps_per_particle` を粒子ごとの上限とする
- `output.diagnostics.charge_rel_change_threshold` は監視値であり、早期停止条件ではない

実装挙動の正本は [`SPEC.md`](../SPEC.md)、設定キーの正本は
[入力パラメータ](reference/Parameters.md) です。削除済みの `[outer_plasma]`、
`[coupling]`、旧 selector は設定として受理されません。

## 最短の実行手順

```bash
pip install beach-bem
beachx config init
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

リポジトリから開発する場合:

```bash
python -m pip install -e . --no-build-isolation
make check
make run CONFIG=examples/beach.toml
```

`beachx lint` が通ることは設定契約を満たすことを示しますが、物理解の収束や妥当性を
保証しません。実行後は出力と目的量を別途確認してください。

## ケースの選び方

| 目的 | 出発点 | 次に読むページ |
|---|---|---|
| 入門の帯電更新 | `examples/tutorial_electron_beam.toml` | [入門チュートリアル](start/Tutorial.md) |
| 1 粒子の開発用 smoke | `examples/smoke_single_electron.toml` | [設定レシピ](guide/ConfigurationRecipes.md) |
| 外部 plasma reservoir | `[particles.species.inflow]` | [Reservoir 注入](guide/ReservoirInjection.md) |
| 閉じた光電子再分配 | `examples/closed_photoelectron.toml` | [光電子放出](models/PhotoelectronEmission.md) |
| free-space の場評価 | `fields.boundary="free"` | [場の評価](guide/FieldSolvers.md) |
| 2 軸周期の場評価 | `fields.boundary="periodic2"` | [周期静電場](models/PeriodicElectrostatics.md) |
| 出力の解析 | `beachx inspect` / `Beach(...)` | [後処理チュートリアル](guide/PostprocessTutorial.md) |
| checkpoint から再開 | `[run.restart]` | [実行と再開](guide/Execution.md) |

境界 reservoir + closed PE の例は外部シースを自己無撞着に解くモデルではありません。
box 高さ、周期画像、`particles.tracking.dt_s`、`particles.tracking.max_steps_per_particle`、`run.batch.duration_s`、粒子数、ray 数は目的量に
対して収束確認が必要です。

## 設定を編集するとき

通常の流れは次の 4 段階です。

1. 近い `examples/*.toml` をコピーする。
2. [設定レシピ](guide/ConfigurationRecipes.md) で必要な table を選ぶ。
3. [入力パラメータ](reference/Parameters.md) で単位、既定値、排他条件を確認する。
4. `beachx lint beach.toml` を実行する。lint は schema 検査の後に Fortran の preflight (`beach --check-config`) を呼ぶ。

特に次の制約を先に確認してください。

| 機能 | 必須条件 |
|---|---|
| `inflow` | `[domain]`のbox、解決後の`run.batch.duration_s>0`、非周期面、外部VDF |
| `photo_raycast` | `[domain]` の box、解決後の `run.batch.duration_s>0`、正の電流密度、`sampling.rays_per_batch>=1` |
| closed PE | 負電荷、`source.deposit_opposite_charge_on_emit=true`、`[particles.species.boundary]` で `source.inject_face` と同じ面が `reflect` または `redistributed_reflect`、`charging.closure="neutral_return"` |
| `periodic2` | `[domain]` の box、`domain.periodic_axes` がちょうど 2 軸、`[fields.periodic]` の zero-mode 設定 |
| resume | `output.enabled=true`、必要な checkpoint、同じ MPI size |

`periodic2` の通常経路は FMM です。Direct は
`triangle_p0 + panel_spectral_reference + exclude_k0` の検証用 split reference に限られます。
有限画像近似の `fields.periodic.backend="finite_images"` と、無限周期の非零モードを扱う
`"cached_kneq0"` を区別してください。詳細な互換表は
[場の評価](guide/FieldSolvers.md#solverと場境界の互換表)を正本とします。

## 実行と検証

開発時の代表コマンド:

```bash
make check       # Fortran build check
make test-l0     # static/schema/build
make test        # L1: Python + quick Fortran
pytest -q        # Python only
ruff check .     # Python lint
```

長時間の FMM、MPI、release 検証は通常ループから分離されています。

```bash
make test-l2
make test-l3
make test-heavy
make test-fortran-far-correction
make test-field-kernel-cache
make test-mpi
```

個別 Fortran target は
`FPM_ACTION=test ./build.sh --target <name>` で実行できます。KUDPC ではホスト役割を
確認し、ログインノード上で test payload を直接実行せず、`tssrun` または
`sbatch` + `srun` を使ってください。詳細は[開発ワークフロー](developer/Workflow.md)にあります。

テスト通過と物理検証は分けて報告します。物理ケースでは少なくとも、時間刻み、
空間離散、box/周期画像、粒子統計、未解決率に対する目的量の変化を確認してください。

## 出力を確認するとき

まず次を確認します。

```bash
beachx inspect outputs/latest
```

主要ファイルは `summary.txt`、`charges.csv`、`mesh_triangles.csv`、
`mesh_sources.csv`、`charge_ledger.csv`、`checkpoint_complete.txt` です。履歴、potential、
reservoir residual、性能 profile は設定に応じて追加されます。生成条件と restart 契約は
[出力ガイド](guide/OutputGuide.md)を正本とします。

Python では次の入口を使います。

```python
from beach import Beach

run = Beach("outputs/latest")
fig, ax = run.plot_mesh()
fig, ax = run.plot_charges(step=-1)
```

CLI の全オプションは `beachx <command> --help`、Python API は
[Python 後処理 API](reference/PythonPostprocessAPI.md)を参照してください。

## 変更時のチェックリスト

- Fortran の挙動と `SPEC.md` の記述が一致している
- 公開設定を変えた場合、schema、parser、例、日英ドキュメント、テストを更新した
- field、collision、boundary、injection、resume の変更には対応する回帰テストがある
- `output.diagnostics.charge_rel_change_threshold` を早期停止として説明していない
- 実行完了と数値・物理妥当性を混同していない
- 日英ページでコマンド、識別子、警告、リンク先の責務が対応している

## ドキュメント案内

| 質問 | 正本 |
|---|---|
| シミュレータは何を計算するか | [`SPEC.md`](../SPEC.md)、[アルゴリズム](models/Algorithms.md) |
| ケースをどう組み立てるか | [設定レシピ](guide/ConfigurationRecipes.md) |
| キー、型、単位、制約は何か | [入力パラメータ](reference/Parameters.md) |
| 実行、負荷見積もり、再開はどうするか | [実行と再開](guide/Execution.md) |
| 粒子源と境界をどう選ぶか | [粒子源と境界](guide/ParticleSourcesBoundaries.md) |
| 衝突と box event はどう処理されるか | [粒子イベント](developer/ParticleEvents.md) |
| 場 solver と周期設定をどう選ぶか | [場の評価](guide/FieldSolvers.md)、[FMM](guide/FMM.md) |
| 出力をどう読み、可視化するか | [出力ガイド](guide/OutputGuide.md)、[後処理チュートリアル](guide/PostprocessTutorial.md) |
| 開発テストをどう選ぶか | [開発ワークフロー](developer/Workflow.md)、[検証ガイド](guide/ValidationGuide.md) |
