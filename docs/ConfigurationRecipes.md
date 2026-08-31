title: ケース設計の流れ

Lang: [日本語](ConfigurationRecipes.md) | [English](ConfigurationRecipes.en.md)

# ケース設計の流れ

このページは、動作確認済みの公式チュートリアルを研究用ケースへ変更するときの判断順を示します。
設定全体をここへ複製せず、各段階で必要な最小差分だけを示します。全キーと組合せ制約は
[入力パラメータリファレンス](Parameters.html)を正本とします。

**出発点:** [10 分チュートリアル](Tutorial.html)を完了し、`beach.toml` と
`outputs/tutorial` がある作業ディレクトリで、基準設定を複製します。

```bash
cp beach.toml case.toml
beachx lint case.toml
```

以降は一つの判断だけを変更し、そのたびに `beachx lint case.toml` を実行します。チュートリアルの
出力は基準結果として残し、新しい出力先は手順 7 で指定します。

## 1. 目的と合格条件を決める

設定を編集する前に、次を一文ずつ書き出します。

- どの表面へ、どの環境から粒子が到達するか
- 観測したい量は電荷、電位、吸収率、escape、または計算時間のどれか
- どの保存則、基準解、収束幅を満たせば結果を採用するか

公式チュートリアルは batch 間の帯電 feedback を学ぶ教材です。そこで使う粒子重みや
`batch_duration=0` を、そのまま特定の実環境の時間発展とは解釈できません。先に合格条件を決めると、
後の mesh 分割、粒子数、solver 精度を必要な範囲へ絞れます。

## 2. 形状を選ぶ

簡単な形状は `mode="template"` の組み込み形状を使います。実測形状や CAD 由来形状は三角形面を持つ
OBJ を使います。複数の独立した物体、特に独立した導体を区別したい場合は、別々の template として
`mesh_id` を分けます。

| `kind` | 形状 | 主な寸法 | 主な分割数 |
| --- | --- | --- | --- |
| `plane` | XY 長方形平面 | `size_x`, `size_y` | `nx`, `ny` |
| `plate_hole`, `plane_hole` | 円形穴付き長方形平面（後者は別名） | `size_x`, `size_y`, `radius` | `n_theta`, `n_r` |
| `disk` | 円板 | `radius` | `n_theta`, `n_r` |
| `annulus` | 同心リング | `radius`, `inner_radius` | `n_theta`, `n_r` |
| `box` | 閉じた直方体表面 | `size` | `nx`, `ny`, `nz` |
| `cylinder` | z 軸方向の円柱 | `radius`, `height`, `cap` | `n_theta`, `n_z` |
| `sphere` | 球面 | `radius` | `n_lon`, `n_lat` |

たとえば、チュートリアルの平面を球へ置き換える場合は、既存の `[[mesh.templates]]` block だけを
次のように置き換えます。

```toml
[[mesh.templates]]
kind = "sphere"
enabled = true
surface_side = "outward_closed"
center = [0.5, 0.5, 0.5]
radius = 0.2
n_lon = 24
n_lat = 12
```

OBJ へ切り替える場合は既存の `[[mesh.templates]]` block を削除し、`[mesh]` を次へ置き換えます。

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj" # 実際の OBJ path に置き換える
surface_side = "outward_closed"
```

`outward_closed` は、面の向きが整合した閉じた two-manifold にだけ使えます。開いた面では、OBJ の面法線と
真空側に応じて `normal_plus` または `normal_minus` を選びます。配置変換、要素数、各形状の制約は
[組み込み形状と OBJ のリファレンス](Parameters.html#mesh-メッシュ入力)、複数形状の実例は
[`examples/beach.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/beach.toml)を参照してください。

## 3. 表面 model を選ぶ

| 目的 | `surface_model` | 現行モデルの範囲 |
| --- | --- | --- |
| 命中位置へ電荷を蓄積する | `insulator` | 表面内伝導、bulk 漏洩、誘電分極を解かない |
| 物体の総電荷を保って等電位化する | `conductor` | 自由空間中の浮遊導体。`field_boundary.mode="free"` のみ |

チュートリアルは `insulator` です。導体へ変える場合は対象の template、または OBJ 用の `[mesh]` で
次の 1 行を変更します。

```toml
surface_model = "conductor"
```

`dielectric`、`epsilon_r`、抵抗性表面は現行入力では利用できません。導体へ変えたケースは、チュートリアルの
結果を基準にせず、物体電位と総電荷を別に検証します。モデルの意味と制約は
[表面はどう帯電するか](SurfaceModels.html)を参照してください。

## 4. 粒子源を選ぶ

粒子をどこから供給するかで経路を選びます。`source_mode` を変更するときは、旧 mode 専用キーを残さず、
その species entry を専用ページの最小設定へ置き換えます。

| `source_mode` | 新しいケースでの用途 | 物理的な供給量 | macro sampling の調整 |
| --- | --- | --- | --- |
| `volume_seed` | box 内へ指定個数を置く軌道試験や初期粒子 | 面流束には対応しない | `npcls_per_step` |
| `plane_source` | box 内部の矩形面から一方向に供給 | 流束 × 面積 × `batch_duration` | `w_particle` または `target_macro_particles_per_batch` |
| `photo_raycast` | 照射 ray が命中した表面から放出 | 電流密度 × 投影面積 × `batch_duration` | `rays_per_batch` と ray の命中率 |
| `reservoir_face` | deprecated な既存設定を読むための互換入力のみ | 新規ケースでは使用しない | 新規ケースでは使用しない |

チュートリアルの `volume_seed` を維持して供給量や領域だけ変えるなら、既存 species entry の次の値だけを
調整します。

```toml
npcls_per_step = 500
pos_low = [0.2, 0.2, 0.8]
pos_high = [0.8, 0.8, 0.8]
drift_velocity = [0.0, 0.0, -1.0e6]
```

外部 plasma から非周期の box 面全体を通す流入は `source_mode` ではなく
`[particles.species.boundary_inflow]` です。現行 schema では、境界流入だけの species にも
`source_mode="volume_seed"` と `npcls_per_step=0` を指定します。選択の詳細は
[粒子をどこから入れるか](ParticleSourcesBoundaries.html)、境界流入は
[境界から粒子を流入させる](ReservoirInjection.html)、光電子は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)を参照してください。

## 5. box・場境界・粒子境界を決める

まず `[domain]` の box が mesh と到達しうる全粒子軌道を含むようにします。次に場の closure を
`field_boundary.mode`、非周期面を出る粒子の作用を `[particle_boundary]` で選びます。

| 判断 | 選択 |
| --- | --- |
| 有限物体の自由空間場 | `periodic_axes=[]`, `field_boundary.mode="free"` |
| x/y に無限反復する場 | `periodic_axes=["x", "y"]`, `field_boundary.mode="periodic2"` |
| 非周期面を出た粒子 | `open`, `reflect`, `redistributed_reflect` |

周期 topology を指定する公開キーは `domain.periodic_axes` だけです。`particle_boundary` や species 別境界で
周期軸を追加、削除、上書きはできません。`periodic2` は x/y 周期・z 非周期に限定されるため、新しく使う場合は
完全例 [`examples/periodic2_basic/beach.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/periodic2_basic/beach.toml)
と [periodic2 静電場](PeriodicElectrostatics.html)から始めてください。粒子の escape と return は
[粒子の escape と return](ParticleEscapeReturn.html)にまとめています。

## 6. solver を選ぶ

| `sim.field_solver` | 選ぶ場面 | 確認方法 |
| --- | --- | --- |
| `direct` | 小規模計算と基準解 | 近似 solver の比較基準にする |
| `treecode` | 中規模の自由空間計算 | opening 条件を変えて観測量を比較する |
| `fmm` | 大規模計算、多数の評価点、通常の `periodic2` | `tree_theta` と `tree_leaf_max` を変えて Direct と比較する |
| `auto` | 自由空間で要素数に応じ自動選択 | 出力で実際に選ばれた solver を確認する |

チュートリアルの基準計算は `field_solver="direct"` のまま残します。高速 solver を採用する前に、同じ mesh と
粒子条件を縮小したケースで Direct との差を測ります。互換性と選択基準は
[場の評価方法](FieldSolvers.html)、FMM の設定と精度調整は [FMM を使う](FMM.html)を参照してください。

## 7. 出力先を分けて実行する

最後に、既存の `[output]` で `dir` だけを変更し、基準結果を上書きしない出力先へ分けます。

```toml
dir = "outputs/case"
```

HPC では login node で simulation を直接実行せず、計算 node の割当内で `beach case.toml` を実行します。

```bash
beachx lint case.toml
beach case.toml
beachx inspect outputs/case
```

ローカル実行、MPI、checkpoint と再開は[実行・再開する](Execution.html)を参照してください。

最小の合格条件は、`lint` が `status=ok`、実行が終了 code 0、`inspect` の `batches` が
`sim.batch_count` と一致し、`summary.txt` と `charges.csv` が存在することです。これは完走した証拠であり、
物理的に正しい証拠ではありません。各出力の読み方は[出力ファイルを調べる](OutputGuide.html)を参照してください。

## 8. 妥当性を確認する

目的に対応する観測量について、少なくとも次を確認します。

- `processed_particles = absorbed + escaped` が成り立ち、`escaped` の内訳である `escaped_boundary`、
  `survived_max_step`、`multiple_box_events_soft_discarded` を説明でき、表面・放出を含む電荷 ledger が閉じる
- mesh を細分化しても結論が変わらない
- `dt`、`max_step`、`batch_duration`、batch 数、粒子数を変えても結論が許容幅内にある
- Treecode や FMM の小規模結果が Direct の基準結果と一致する
- 乱数 seed または反復実行を変えた統計変動が許容幅内にある
- 採用する結論が、選んだ表面 model、境界、粒子源の適用範囲を超えない

変更前後を同じ物理量・同じ時刻で比較し、一度に複数の数値条件を変えません。具体的な確認コマンドと
採用基準は[計算結果の妥当性確認](ValidationGuide.html)に従ってください。
