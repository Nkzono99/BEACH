title: BEACH 入力パラメータリファレンス

Lang: [日本語](Parameters.md) | [English](Parameters.en.md)

# 入力パラメータリファレンス

本文書は、Fortran 実行系が読む最終的な `beach.toml` のパラメータリファレンスです。
単位は、特に断りがない限り SI 単位です。

初めて設定を組む場合は、先に [設定レシピ](ConfigurationRecipes.html) を読むと全体像を掴みやすいです。

Fortran parser が解決する高水準記法は [Configuration](Configuration.html) にまとめています。
本書では、正規化後に実行時設定として使われるキーを中心に説明します。

| 関連ドキュメント | 内容 |
|---|---|
| [設定レシピ](ConfigurationRecipes.html) | 目的別の設定手順と調整ポイント |
| [Configuration](Configuration.html) | `beachx config`、高水準記法、schema/lint |
| [Algorithms](Algorithms.html) | BEM 場計算、粒子 push、衝突、蓄積電荷の計算手順への導線 |
| [Workflow](Workflow.html) | 実行、開発、テスト、KUDPC での注意 |
| [FMMCore](FMMCore.html) | `field_solver="fmm"` の数値アルゴリズム |

---

## 目次

- [入力パラメータリファレンス](#入力パラメータリファレンス)
  - [目次](#目次)
  - [読み込みルール](#読み込みルール)
  - [単位と座標](#単位と座標)
  - [最小例](#最小例)
  - [セクション一覧](#セクション一覧)
  - [パラメータ詳細リファレンス](#パラメータ詳細リファレンス)
    - [`[sim]`: 実行制御と場計算](#sim-実行制御と場計算)
    - [`[[particles.species]]`: 粒子種](#particlesspecies-粒子種)
    - [`sim.sheath_injection_model`: シース流入補正](#simsheath_injection_model-シース流入補正)
    - [`[mesh]`: メッシュ入力](#mesh-メッシュ入力)
    - [`[[mesh.templates]]`: 組み込み形状](#meshtemplates-組み込み形状)
    - [`[output]`: 出力と再開](#output-出力と再開)
  - [高水準記法との関係](#高水準記法との関係)
  - [検証ルール](#検証ルール)

---

## 読み込みルール

| 項目 | 仕様 |
|---|---|
| 明示指定 | `beach path/to/beach.toml` |
| 既定入力 | 引数なしではカレントディレクトリの `beach.toml` |
| 開発実行 | `fpm run -- path/to/beach.toml` でも同じ |
| 形式 | TOML。複数行配列も利用可能 |
| 未知キー | 未知のセクション名・キー名はエラー |
| schema | `schemas/beach.schema.json` |
| lint | `beachx lint beach.toml` |

Editor schema を使う場合は、`beach.toml` の先頭にコメント directive を置きます。
Fortran パーサは最初のセクションより前の通常キーを受け付けないため、`"$schema" = "..."` は使いません。

```toml
#:schema ../schemas/beach.schema.json
```

GitHub Raw URL を指定することもできます。

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

相対パスは、その `beach.toml` 自身から見た相対パスです。
`outputs/.../beach.toml` など深い場所に置く場合は `../../schemas/beach.schema.json` のように調整してください。

---

## 単位と座標

| 種類 | 代表キー | 単位・向き |
|---|---|---|
| 時間 | `dt`, `batch_duration` | 秒 |
| 長さ | `box_min`, `box_max`, `pos_low`, `pos_high` | m |
| 電荷 | `q_particle`, 要素電荷出力 | C |
| 質量 | `m_particle` | kg |
| 速度 | `drift_velocity`, `ray_direction` | m/s。ただし `ray_direction` は方向ベクトル |
| 電場 | `e0`, `e0_abs` | V/m |
| 磁場 | `b0` | T |
| 密度 | `number_density_cm3`, `number_density_m3` | cm^-3 または m^-3 |
| 温度 | `temperature_k`, `temperature_ev` | K または eV。両方の同時指定は不可 |
| 角度 | `e0_phi_xy_deg`, `e0_phi_z_deg`, `sheath_alpha_deg` | degree |

`*_low` / `*_high` は各軸の下限・上限です。
`inject_face` は `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, `z_high` のいずれかを指定します。

---

## 最小例

初見では、物理流入を直接指定できる `source_mode = "reservoir_face"` を推奨します。

```toml
[sim]
dt = 2.0e-8
batch_duration_step = 60000.0
batch_count = 100
max_step = 10000
softening = 1.0e-6
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_far_correction = "none"

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 5000
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
target_macro_particles_per_batch = -1
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]

[output]
write_files = true
write_mesh_potential = false
dir = "outputs/latest"
history_stride = 1
```

---

## セクション一覧

| セクション | 必須 | 内容 |
|---|---:|---|
| `[sim]` | 条件付き | 時間刻み、バッチ数、場ソルバ、境界、外部場、シース補正 |
| `[particles]` | yes | `[[particles.species]]` のコンテナ |
| `[[particles.species]]` | yes | 粒子種、注入方式、速度分布、マクロ粒子重み |
| `[mesh]` | no | OBJ または組み込み template の選択 |
| `[[mesh.templates]]` | no | `mode="template"` で使う組み込み形状 |
| `[output]` | no | 出力先、履歴、checkpoint 再開 |

`reservoir_face` または `photo_raycast` を使う場合、`[sim]` は必須です。
`[[particles.species]]` は 1 件以上必要です。

---

## パラメータ詳細リファレンス

### `[sim]`: 実行制御と場計算

#### 実行制御

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `dt` | float | `1.0e-9` | 時間刻み [s] |
| `rng_seed` | int | `12345` | 乱数シード |
| `batch_count` | int | `1` | 通常実行では処理するバッチ数。`output.resume=true` では累積の到達バッチ数 |
| `batch_duration` | float | `0.0` | 1 バッチの物理時間 [s] |
| `batch_duration_step` | float | `0.0` | `batch_duration = dt * batch_duration_step` として解決 |
| `max_step` | int | `400` | 粒子 1 個あたりの最大 push 回数 |
| `tol_rel` | float | `1.0e-8` | 相対変化量の監視値。停止条件には未使用 |
| `q_floor` | float | `1.0e-30` | `rel_change` 計算時の分母下限 |

`batch_duration` と `batch_duration_step` の同時指定はエラーです。
`reservoir_face` / `photo_raycast` では、解決後の `batch_duration > 0` が必須です。

#### 場ソルバ

`field_solver` は、境界要素電荷から評価点の Coulomb 電場を計算する方式です。
選択肢ごとの対応パラメータは下表のとおりです。

| `field_solver` | 用途 | 対応する場境界 |
|---|---|---|
| `direct` | 要素数が小さい場合の厳密な全対全評価 | `field_bc_mode="free"` |
| `treecode` | 中規模以上の近似評価 | `field_bc_mode="free"` |
| `fmm` | 大規模評価、`periodic2`、FMM コア検証 | `field_bc_mode="free"` / `"periodic2"` |
| `auto` | 要素数に応じて direct / treecode を自動選択 | `field_bc_mode="free"` |

共通キー:

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `softening` | float | `1.0e-6` | Coulomb 場計算の softening 長さ [m] |
| `field_solver` | string | `"auto"` | `direct` / `treecode` / `fmm` / `auto` |
| `field_normalization` | string | `"si"` | `si` / `box` / `mesh` / `length` |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` で使う長さスケール [m] |

`field_normalization` は場計算内部の座標・softening・周期 cell の正規化だけを変えます。
出力される電場・電位は SI に戻されます。

| `field_normalization` | 長さスケール |
|---|---|
| `si` | 入力 SI 座標をそのまま使う |
| `box` | `sim.box_max - sim.box_min` の最大幅。`sim.use_box=true` が必須 |
| `mesh` | mesh bounding box の最大幅。mesh が空なら `field_length_scale` |
| `length` | `field_length_scale` |

##### `field_solver = "direct"`

全 source 要素を直接足し上げます。
近似誤差はなく、計算量は評価点数を `M`、要素数を `N` として `O(MN)` です。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"direct"` を指定 |
| `softening` | float | `1.0e-6` | source と評価点が近いときの特異性を緩和 |
| `field_normalization` | string | `"si"` | direct 評価前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `field_bc_mode` | string | `"free"` | `direct` では `"free"` のみ |

`tree_theta`、`tree_leaf_max`、`tree_min_nelem` は `direct` では使いません。

##### `field_solver = "treecode"`

source octree を作り、遠方 node は multipole 近似、近傍 node は direct 和で評価します。
FMM のような local expansion は使わず、評価点ごとに木を走査します。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"treecode"` を指定 |
| `softening` | float | `1.0e-6` | 近傍 direct 和と multipole 評価の softening |
| `field_normalization` | string | `"si"` | tree 構築前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_theta` | float | `0.5` | MAC パラメータ。`0 < theta <= 1`。大きいほど速く粗い |
| `tree_leaf_max` | int | `16` | 葉 node あたり最大 source 数。`>= 1` |
| `field_bc_mode` | string | `"free"` | `treecode` では `"free"` のみ |

`tree_min_nelem` は `field_solver="auto"` 用のしきい値なので、明示 `treecode` では切替には使いません。

##### `field_solver = "fmm"`

simulator 非依存の Coulomb FMM コアを使います。
source 幾何の plan と、電荷更新ごとの state を分け、P2M/M2M/M2L/L2L/L2P と近傍 direct 和で評価します。
詳細は [Coulomb FMM コア詳細](FMMCore.html) を参照してください。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"fmm"` を指定 |
| `softening` | float | `1.0e-6` | 近傍 direct 和と FMM 評価の softening |
| `field_normalization` | string | `"si"` | FMM plan 構築前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_theta` | float | `0.5` | near/far 判定の MAC パラメータ。`0 < theta <= 1` |
| `tree_leaf_max` | int | `16` | source tree の葉 node あたり最大 source 数。`>= 1` |
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | `periodic2` の近傍画像層数 |
| `field_periodic_far_correction` | string | `"none"` | `periodic2` の遠方補正。`none` / `auto` / `m2l_root_oracle` |
| `field_periodic_ewald_alpha` | float | `0.0` | `m2l_root_oracle` 用 Ewald 分解パラメータ |
| `field_periodic_ewald_layers` | int | `4` | Ewald oracle の real-space outer shell / reciprocal cutoff 深さ |

`field_periodic_*` は `field_bc_mode="periodic2"` のときだけ意味を持ちます。
`tree_min_nelem` は明示 `fmm` では使いません。

##### `field_solver = "auto"`

要素数が `tree_min_nelem` 未満なら direct、以上なら treecode を使います。
`auto` は FMM へは切り替えません。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"auto"` を指定 |
| `softening` | float | `1.0e-6` | direct / treecode の softening |
| `field_normalization` | string | `"si"` | 自動選択前に共通で使う正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_min_nelem` | int | `256` | treecode へ切り替える要素数しきい値。`>= 1` |
| `tree_theta` | float | `0.5` | treecode 選択時の MAC パラメータ |
| `tree_leaf_max` | int | `16` | treecode 選択時の葉 node あたり最大 source 数 |
| `field_bc_mode` | string | `"free"` | `auto` では `"free"` のみ |

`tree_theta` と `tree_leaf_max` は、明示指定がなければ要素数から次の値を使います。

| 要素数 `nelem` | `tree_theta` | `tree_leaf_max` |
|---:|---:|---:|
| `< 1500` | `0.40` | `12` |
| `1500 <= nelem < 10000` | `0.50` | `16` |
| `10000 <= nelem < 50000` | `0.58` | `20` |
| `50000 <= nelem` | `0.65` | `24` |

#### 場境界と periodic2

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | `periodic2` の近傍画像層数 |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `m2l_root_oracle` |
| `field_periodic_ewald_alpha` | float | `0.0` | `m2l_root_oracle` 用 Ewald 分解パラメータ |
| `field_periodic_ewald_layers` | int | `4` | Ewald oracle の outer shell / reciprocal cutoff 深さ |

legacy `periodic2` は `field_solver="fmm"` を使います。小規模検証用のsplit referenceだけは、`field_solver="direct"`と以下の3 tableを明示します。

| table.key | 既定 | 意味 |
|---|---:|---|
| `periodic2.nonzero_mode_backend` | 必須 | `panel_spectral_reference` |
| `periodic2.zero_mode_policy` | 必須 | `exclude_k0` |
| `periodic2.lower_boundary_model` | 必須 | `e_bottom_zero` |
| `periodic2.reference_mode_layers` | `4` | Fourier mode cutoff |
| `periodic2.panel_quadrature_order` | `12` | panel面積積分次数 |
| `periodic2.interface_sample_n` | `5` | interface各軸の診断点数 |
| `periodic2.interface_phi_tolerance` | `1e-3` | 非零モード電位比上限 |
| `periodic2.interface_field_tolerance` | `1e-3` | 非零モード電場比上限 |
| `outer_plasma.interface_z` | 必須 | z上側interface。初期モデルではbox上面 |
| `outer_plasma.debye_length` | 必須 | 線形Debye長 |
| `outer_plasma.thermal_voltage` | 必須 | 線形性・診断の電位scale |
| `outer_plasma.max_linearity_ratio` | `0.25` | `abs(phi_t-phi_inf)/thermal_voltage`上限 |
| `outer_plasma.max_gap_ratio` | `5` | `(z_t-z_mesh,max)/lambda`上限 |
| `outer_plasma.max_local_charge_ratio` | `50` | 局所平均plasma電荷推定比上限 |
| `coupling.outer_update_stride` | `1` | outer profile更新batch間隔 |
| `outer_plasma.return_model` | `none` | `electrostatic_1d_instant_return`で個別粒子を返却 |
| `coupling.particle_transfer_mode` | `none` | return modelと同じIDを指定 |
| `coupling.field_evolution_timescale` | `0` | frozen-field比較時間 [s]。instant returnでは正値必須 |
| `coupling.max_frozen_field_ratio` | `0.1` | `tau_outer/field_evolution_timescale`上限 |
| `coupling.outer_queue_enabled` | `false` | 現在は`true`を拒否 |

完全な例は`examples/periodic2_linear_outer_reference.toml`です。閾値違反時にlegacy modelへfallbackしません。
粒子移送を含む例は`examples/periodic2_outer_particle_transfer.toml`です。instant returnはz-high、`b0=0`、x/y periodicだけに対応します。
`sim.use_box=true`、かつ 2 軸だけが `periodic`、残り 1 軸が開放のときに有効です。
場評価だけでなく、collision と `photo_raycast` の raycast でも periodic image を考慮します。

`field_periodic_far_correction="auto"` は互換用に受理され、現在は `none` と同じ扱いです。
`m2l_root_oracle` は build 時の診断用で、exact periodic Ewald residual を root operator に fit する高コストモードです。

#### 外部場

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | 一様外部電場 [V/m] |
| `e0_abs` | float | 未指定 | 一様外部電場の大きさ [V/m] |
| `e0_phi_xy_deg` | float | `0.0` | `e0_abs` 指定時の xy 面内方位角 [deg] |
| `e0_phi_z_deg` | float | `0.0` | `e0_abs` 指定時の xy 面からの仰角 [deg] |
| `b0` | float[3] | `[0,0,0]` | 一様磁場 [T] |

一様外部電場は、`e0 = [Ex, Ey, Ez]` で直接指定するか、`e0_abs` と角度で指定します。
両形式の混在はエラーです。

#### 流入補助とシース補正

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `reservoir_potential_model` | string | `"none"` | `none` / `infinity_barrier` |
| `phi_infty` | float | `0.0` | 無限遠基準電位 [V] |
| `open_boundary_model` | string | `"escape"` | `escape` / `potential_barrier` |
| `injection_face_phi_grid_n` | int | `3` | 注入面平均電位の `N x N` 評価格子 |
| `raycast_max_bounce` | int | `16` | `photo_raycast` の最大反射回数 |
| `sheath_injection_model` | string | `"none"` | `none` / `zhao_auto` / `zhao_a` / `zhao_b` / `zhao_c` / `floating_no_photo` |
| `sheath_alpha_deg` | float | `60.0` | Zhao シースの太陽高度角 [deg] |
| `sheath_photoelectron_ref_density_cm3` | float | `64.0` | Zhao シースの基準光電子密度 [cm^-3] |
| `sheath_reference_coordinate` | float | 未指定 | シース 1D 座標の基準平面位置 [m] |
| `sheath_electron_drift_mode` | string | `"normal"` | `normal` / `full` |
| `sheath_ion_drift_mode` | string | `"normal"` | `normal` / `full` |

`sheath_injection_model != "none"` は、現状 `reservoir_potential_model="none"` と組み合わせて使います。
詳細は [`sim.sheath_injection_model`](#simsheath_injection_model-シース流入補正) を参照してください。

#### 計算領域と粒子境界

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `use_box` | bool | `false` | ボックス境界を有効化 |
| `box_min` | float[3] | `[-1,-1,-1]` | ボックス下限 [m] |
| `box_max` | float[3] | `[1,1,1]` | ボックス上限 [m] |
| `bc_x_low`, `bc_x_high` | string | `"open"` | x 方向下限・上限の粒子境界 |
| `bc_y_low`, `bc_y_high` | string | `"open"` | y 方向下限・上限の粒子境界 |
| `bc_z_low`, `bc_z_high` | string | `"open"` | z 方向下限・上限の粒子境界 |

粒子境界は `open`, `reflect`, `periodic` を指定します。
`open` は `outflow`, `escape` も同義語として受理されます。

`open_boundary_model="potential_barrier"` では、`open` 面を越えた粒子について、
境界通過点の BEM 電位 `phi_boundary` と `phi_infty` から電位障壁
`q_particle * (phi_infty - phi_boundary)` を評価します。障壁が正で、開境界法線方向の運動エネルギー
`0.5 * m_particle * v_normal^2` より大きい場合は法線速度を反転して反射し、それ以外は脱出として扱います。
一様外部電場 `e0` は無限遠基準の電位を定義しないため、この障壁評価には含めません。

`periodic2` の mesh は、runtime で collision 用 canonical unwrapped 表現へ平行移動してから ray-triangle 判定します。
raw 頂点は periodic 軸で box 外を含んでも構いませんが、triangle を頂点ごとに mod 折り返すことはしません。

---

### `[[particles.species]]`: 粒子種

`[[particles.species]]` は 1 件以上必須です。
`source_mode` によって、使うキーと制約が変わります。

#### 共通キー

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `enabled` | bool | `true` | 種を有効化 |
| `source_mode` | string | `"volume_seed"` | `volume_seed` / `reservoir_face` / `photo_raycast` |
| `q_particle` | float | `-1.602176634e-19` | 粒子電荷 [C] |
| `m_particle` | float | `9.10938356e-31` | 粒子質量 [kg] |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | 位置下限 [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | 位置上限 [m] |
| `drift_velocity` | float[3] | `[0,0,-8e5]` | ドリフト速度 [m/s] |
| `temperature_k` | float | `2.0e4` | 温度 [K] |
| `temperature_ev` | float | 未指定 | 温度 [eV]。`temperature_k` と排他 |
| `velocity_distribution` | string | `"maxwellian"` | `maxwellian` / `grid` |
| `inject_face` | string | 未指定 | 注入面。`reservoir_face` / `photo_raycast` で必須 |

#### `source_mode = "volume_seed"`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `npcls_per_step` | int | `0` | 1 バッチに生成するマクロ粒子数 |
| `w_particle` | float | `1.0` | マクロ粒子重み |

制約:

| 条件 | 内容 |
|---|---|
| 粒子数 | 有効 species 全体で `npcls_per_step` 合計が 1 以上必要 |
| 重み自動解決 | `target_macro_particles_per_batch` は使用不可 |

#### `source_mode = "reservoir_face"`

| キー | 型 | 説明 |
|---|---|---|
| `number_density_cm3`, `number_density_m3` | float | 上流密度。どちらか片方を指定 |
| `w_particle` | float | マクロ粒子重み。正値 |
| `target_macro_particles_per_batch` | int | `w_particle` 自動解決用。`>0` または `-1` |
| `velocity_grid_path` | string | `velocity_distribution="grid"` の CSV パス |
| `velocity_grid_pdf_kind` | string | `phase_space` / `flux_weighted` |
| `velocity_grid_sampling` | string | `auto` / `rectilinear` / `discrete` |
| `particle_flux_m2_s`, `current_density_a_m2` | float | `grid` 分布の流入量。どちらか片方を指定 |

基本制約:

| 条件 | 内容 |
|---|---|
| 領域 | `sim.use_box=true` が必須 |
| 時間 | `sim.batch_duration > 0` が必須 |
| 注入面 | `inject_face` が必須 |
| 注入範囲 | `pos_low` / `pos_high` は指定 face 上にある必要がある |
| 重み | `w_particle` と `target_macro_particles_per_batch` は同時指定不可 |
| 重み共有 | `target_macro_particles_per_batch=-1` は species 2 以降だけ可。species 1 の `w_particle` を共有 |

Maxwellian 分布では、`number_density_*` と温度から drifting Maxwellian の片側流束を計算します。

Grid 分布では、`velocity_grid_path` の CSV を読み込みます。
必要列は `vx_m_s, vy_m_s, vz_m_s, f` です。
`f` は内部で `sum f = 1` に正規化されます。
この場合、`number_density_*` / `temperature_*` は使わず、`particle_flux_m2_s` または `current_density_a_m2` で流入量を決めます。
`current_density_a_m2` は `abs(J / q_particle)` として粒子 flux に変換します。

| `velocity_grid_sampling` | 挙動 |
|---|---|
| `auto` | 完全な直交格子なら三線形補間。不完全格子や散布点なら離散サンプル |
| `rectilinear` | 直交格子補間を強制。直交格子でなければエラー |
| `discrete` | CSV 行を直接サンプル |

| `velocity_grid_pdf_kind` | 挙動 |
|---|---|
| `phase_space` | 流入面の内向き法線速度 `v_n` を掛けた `v_n f(v)` でサンプル |
| `flux_weighted` | CSV の `f` を流束重み済み分布として扱う |

どちらの PDF でも、`v_n > 0` の速度だけを使います。
`velocity_grid_path` の相対パスは実行時のカレントディレクトリ基準です。
現状、`velocity_distribution="grid"` は `sim.sheath_injection_model="none"` のときのみ有効です。

粒子数は次のように決まります。

```text
n_macro_expected = gamma_in * A * batch_duration / w_particle
n_injected = floor(residual + n_macro_expected)
```

残差は次バッチへ繰り越されます。
`target_macro_particles_per_batch > 0` のときは、その値に近づくよう `w_particle` を自動計算します。

#### `source_mode = "photo_raycast"`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `emit_current_density_a_m2` | float | `0.0` | レイ垂直面基準の放出電流密度 [A/m^2] |
| `rays_per_batch` | int | `0` | 1 バッチの発射レイ数 |
| `deposit_opposite_charge_on_emit` | bool | `false` | 放出元要素に逆符号電荷を堆積 |
| `photo_escape_model` | string | `"none"` | `none` / `boltzmann_cutoff` |
| `normal_drift_speed` | float | `0.0` | 放出法線方向ドリフト [m/s] |
| `ray_direction` | float[3] | 注入面内向き法線 | レイ方向 |

制約:

| 条件 | 内容 |
|---|---|
| 領域 | `sim.use_box=true` が必須 |
| 時間 | `sim.batch_duration > 0` が必須 |
| 放出量 | `emit_current_density_a_m2 > 0`, `rays_per_batch > 0` が必須 |
| 注入面 | `inject_face` が必須 |
| 粒子属性 | `q_particle` は非ゼロ、`m_particle > 0` |
| レイ方向 | 正規化可能で、注入面内向き法線との内積が正 |
| 使用不可 | `npcls_per_step`, `number_density_*`, `w_particle`, `target_macro_particles_per_batch` |

レイ 1 本が命中したときの重みは次の式です。

```text
w_hit = J_perp * A_perp * batch_duration / (|q_particle| * rays_per_batch)
```

実際の生成粒子数はレイの命中率で決まるため、バッチごとの生成数は `rays_per_batch` 以下です。
`field_bc_mode="periodic2"` では、periodic image に命中しても primary cell に wrap した hit 座標から放出します。

`photo_escape_model="boltzmann_cutoff"` では、放出元要素の自己寄与を除いた中心電位で障壁を評価します。

```text
barrier = max(phi_emit - phi_infty, 0)
escape_factor = exp(-|q_particle| * barrier / (k_B * T_PE))
```

`deposit_opposite_charge_on_emit=true` の場合、放出元要素へ残す逆符号電荷にも同じ実効重みを使います。
抑制された光電子電流は即時 return として扱います。

---

### `sim.sheath_injection_model`: シース流入補正

`sim.sheath_injection_model` は、既存の `reservoir_face` / `photo_raycast` species を束ね、シースに対応する流束や法線速度 cutoff を上書きする共有設定です。

| 値 | 内容 |
|---|---|
| `none` | 補正なし |
| `zhao_auto` | 太陽高度角に応じて Zhao の Type A/B/C 分枝を自動選択 |
| `zhao_a`, `zhao_b`, `zhao_c` | Zhao の 1D 光電子シース条件を指定分枝で使用 |
| `floating_no_photo` | 光電子を含まない簡易 floating sheath |

Zhao 系モデルでは、次の species が自動検出されます。

| 対象 | 検出条件 |
|---|---|
| solar-wind electron | 最初の負電荷 `reservoir_face` species |
| ion | 最初の正電荷 `reservoir_face` species |
| photoelectron | 最初の負電荷 `photo_raycast` species |

Zhao 系モデルの効果:

| 対象 | 上書き内容 |
|---|---|
| electron reservoir | 有効密度を Zhao 解の `n_swe_inf` に置換し、障壁に応じた `vmin_normal` を適用 |
| ion reservoir | `sheath_reference_coordinate` 指定時に局所密度・局所法線速度・冷たいビーム近似へ更新 |
| photoelectron | `emit_current_density_a_m2` を Zhao の自由光電子電流へ置換し、`normal_drift_speed=0` として扱う |

`floating_no_photo` では、最初の負電荷 / 正電荷 `reservoir_face` species の電流釣り合いから負の浮遊電位を解きます。
electron reservoir species には cutoff を掛け、`photo_raycast` species があっても放出電流は 0 とみなします。

補足:

- Zhao 系モデルは `temperature_*`, `number_density_*`, `drift_velocity`, `m_particle`, `q_particle` を背景プラズマ条件として再利用します。
- `sheath_reference_coordinate` は、共有 `inject_face` の法線軸に沿った基準平面位置です。
- 例: `inject_face="z_high"` かつ `sheath_reference_coordinate=0.02` なら、平面 `z=0.02` を `z_sheath=0` とみなします。
- `sheath_reference_coordinate` 未指定時は、共有 cutoff ベースの補正だけを適用します。
- Fortran 実装では、Type A は 1 次積分、Type B/C は単調分枝の 1 次積分で局所プロファイルを再構成します。
- `zhao_auto` は `alpha < 20 deg` で `C -> A -> B`、それ以外で `A -> B -> C` の順に分枝解を試みます。

---

### `[mesh]`: メッシュ入力

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `mode` | string | `"auto"` | `auto` / `obj` / `template` |
| `obj_path` | string | `"examples/simple_plate.obj"` | OBJ ファイルパス |
| `surface_model` | string | `"insulator"` | OBJ 全体の表面モデル |
| `epsilon_r` | float | `1.0` | OBJ 全体の相対誘電率。`>= 1` |
| `obj_scale` | float | `1.0` | OBJ 読み込み後の一様スケール |
| `obj_rotation` | float[3] | `[0,0,0]` | OBJ 読み込み後の回転角 [deg] |
| `obj_offset` | float[3] | `[0,0,0]` | OBJ 読み込み後の平行移動 [m] |

`mode="auto"` では、`obj_path` が存在すれば OBJ、なければ template を使います。
OBJ 変換順序は `scale -> rotate -> offset` です。

```text
v_new = R(rotation) * (v_old * obj_scale) + obj_offset
```

OBJ 入力では、ファイル全体を `mesh_id=1` として読みます。
1 つの OBJ 内に離れた `conductor` 部品があっても同じ浮遊導体として扱われます。
独立導体として扱う場合は、template 入力などで `mesh_id` を分けてください。

OBJ の対応範囲:

| 項目 | 対応 |
|---|---|
| 改行 | LF / CRLF |
| 面行 | `f v`, `f v/vt`, `f v/vt/vn`, `f v//vn` |
| 多角形 | 四角形以上はファン三角形分割 |

---

### `[[mesh.templates]]`: 組み込み形状

共通キー:

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `enabled` | bool | `true` | template を有効化 |
| `kind` | string | `"plane"` | `plane` / `plate_hole` / `plane_hole` / `disk` / `annulus` / `box` / `cylinder` / `sphere` |
| `surface_model` | string | `"insulator"` | `insulator` / `conductor` / `dielectric` |
| `surface_side` | string | 未指定 | `triangle_p0` の真空側: `normal_plus` / `normal_minus` / `outward_closed` |
| `epsilon_r` | float | `1.0` | 相対誘電率。`>= 1` |
| `center` | float[3] | `[0,0,0]` | 形状中心 [m] |

`[[mesh.templates]]` を書いた場合、実際に使うテンプレート数は定義件数で決まります。

### `[field]`: 要素核

`element_kernel="point"` が互換既定です。`element_kernel="triangle_p0"` は各要素の `q_elem` を三角形上の一定面密度として扱う厳密 direct 核です。Phase 1 の制約は `sim.field_solver="direct"`、`sim.field_bc_mode="free"`、`sim.softening=0`、全表面 `insulator` です。OBJ では `[mesh].surface_side`、template では各 `[[mesh.templates]].surface_side` を明示してください。`outward_closed` は閉じた向き整合 two-manifold にだけ使えます。
無効化された template は mesh に追加されず、`mesh_id` も消費しません。

`kind` の概要:

| `kind` | 生成形状 | 基準面・軸 |
|---|---|---|
| `plane` | 長方形平面 | XY 平面、`z=center[3]` |
| `plate_hole`, `plane_hole` | 中央に円形穴を持つ長方形平面 | XY 平面、穴中心は `center` |
| `disk` | 円板 | XY 平面、中心は `center` |
| `annulus` | 同心リング | XY 平面、中心は `center` |
| `box` | 閉じた直方体表面 | 各軸に平行な 6 面 |
| `cylinder` | 円柱側面と任意の上下キャップ | z 軸方向 |
| `sphere` | 球面 | 中心は `center` |

#### `kind = "plane"`

XY 平面上の長方形を `nx * ny` 個の矩形セルに分け、各セルを 2 三角形へ分割します。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | 平面中心 `[x, y, z]` [m] |
| `size_x` | float | `1.0` | x 方向サイズ [m]。`> 0` |
| `size_y` | float | `1.0` | y 方向サイズ [m]。`> 0` |
| `nx` | int | `1` | x 方向分割数。`>= 1` |
| `ny` | int | `1` | y 方向分割数。`>= 1` |

要素数は `2 * nx * ny` です。

#### `kind = "plate_hole"` / `"plane_hole"`

XY 平面上の長方形プレートから、中心の円形穴を除いた形状です。
`plane_hole` は `plate_hole` の別名です。
穴境界は `n_theta` 分割の多角形で近似し、穴縁から外周までを `n_r` 層に分けます。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | プレート中心および穴中心 `[x, y, z]` [m] |
| `size_x` | float | `1.0` | x 方向サイズ [m]。`> 0` |
| `size_y` | float | `1.0` | y 方向サイズ [m]。`> 0` |
| `radius` | float | `0.5` | 穴半径 [m]。実行時には `0 < radius < min(size_x, size_y) / 2` |
| `n_theta` | int | `24` | 穴境界の周方向分割数。`>= 3` |
| `n_r` | int | `4` | 穴縁から外周までの半径方向分割数。`>= 1` |

外周は長方形境界に一致します。
円形穴の半径が半幅または半高以上になる設定はエラーです。
共通 default の `radius=0.5` は、既定の `size_x=size_y=1.0` ではこの制約に当たるため、`plate_hole` では `radius` を明示指定してください。

#### `kind = "disk"`

XY 平面上の円板です。
内部は極座標で分割され、中心から外周へ向かって三角形化されます。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | 円板中心 `[x, y, z]` [m] |
| `radius` | float | `0.5` | 円板半径 [m]。`> 0` |
| `n_theta` | int | `24` | 周方向分割数。`>= 3` |
| `n_r` | int | `4` | 半径方向分割数。`>= 1` |

内部的には `inner_radius=0` の `annulus` と同じ生成経路を使います。

#### `kind = "annulus"`

XY 平面上の同心リングです。
内半径から外半径までを `n_r` 層に分けます。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | リング中心 `[x, y, z]` [m] |
| `radius` | float | `0.5` | 外半径 [m]。`> 0` |
| `inner_radius` | float | `0.25` | 内半径 [m]。`0 <= inner_radius < radius` |
| `n_theta` | int | `24` | 周方向分割数。`>= 3` |
| `n_r` | int | `4` | 半径方向分割数。`>= 1` |

`inner_radius=0` も受理されますが、円板を作る場合は `kind="disk"` の方が意図が明確です。

#### `kind = "box"`

閉じた直方体表面です。
6 面すべてを三角形化し、法線は外向きになるように頂点順序を設定します。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | 直方体中心 `[x, y, z]` [m] |
| `size` | float[3] | `[1,1,1]` | x, y, z 方向サイズ [m]。各成分 `> 0` |
| `nx` | int | `1` | x 方向分割数。`>= 1` |
| `ny` | int | `1` | y 方向分割数。`>= 1` |
| `nz` | int | `1` | z 方向分割数。`>= 1` |

要素数は `4 * (nx * ny + ny * nz + nx * nz)` です。
これは、各面の矩形セルを 2 三角形へ分け、対向する 2 面分を数えたものです。

#### `kind = "cylinder"`

z 軸方向の円柱です。
側面を `n_theta * n_z` の矩形セルに分け、必要に応じて上下キャップを追加します。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | 円柱中心 `[x, y, z]` [m] |
| `radius` | float | `0.5` | 円柱半径 [m]。`> 0` |
| `height` | float | `1.0` | z 方向高さ [m]。`> 0` |
| `n_theta` | int | `24` | 周方向分割数。`>= 3` |
| `n_z` | int | `1` | 軸方向分割数。`>= 1` |
| `cap` | bool | `true` | 上下キャップをまとめて有効化 |
| `cap_top` | bool | `cap` の値 | 上面キャップ。指定時は `cap` より優先 |
| `cap_bottom` | bool | `cap` の値 | 下面キャップ。指定時は `cap` より優先 |

円柱は `z = center[3] - height/2` から `z = center[3] + height/2` まで伸びます。
側面の要素数は `2 * n_theta * n_z` です。
各キャップを有効化すると、それぞれ `n_theta` 個の三角形が追加されます。

#### `kind = "sphere"`

経度・緯度分割に基づく球面です。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | 球中心 `[x, y, z]` [m] |
| `radius` | float | `0.5` | 球半径 [m]。`> 0` |
| `n_lon` | int | `24` | 経度方向分割数。`>= 3` |
| `n_lat` | int | `12` | 緯度方向分割数。`>= 2` |

要素数は `2 * n_lon * (n_lat - 1)` です。
極付近は 1 三角形、その他の緯度帯は 2 三角形で構成されます。

表面モデル:

| `surface_model` | 挙動 |
|---|---|
| `insulator` | 衝突粒子の電荷を要素へ蓄積 |
| `conductor` | `mesh_id` ごとの浮遊導体として、総電荷を保存しながら等電位になるよう要素電荷を再配分 |
| `dielectric` | `epsilon_r` をメタデータとして保存。現行の場計算・電荷蓄積では誘電体分極をまだ分岐しない |

`conductor` は現時点で `field_bc_mode="free"` の直接 Coulomb 係数を使って再配分します。
`field_bc_mode="periodic2"` とは併用できません。
導体要素数が大きいケースでは、dense solve によりバッチごとの追加コストが増えます。

---

### `[output]`: 出力と再開

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `write_files` | bool | `true` | ファイル出力の有効/無効 |
| `write_mesh_potential` | bool | `false` | `mesh_potential.csv` を出力 |
| `write_potential_history` | bool | `false` | `potential_history.csv` を出力 |
| `dir` | string | `"outputs/latest"` | 出力先ディレクトリ |
| `history_stride` | int | `1` | 履歴 CSV の出力間隔 [batch] |
| `resume` | bool | `false` | 既存 checkpoint から再開 |
| `restart_from` | string | なし | `resume=true` 時の checkpoint 読み込み元 |

出力ファイル:

| ファイル | 条件・内容 |
|---|---|
| `summary.txt` | 実行統計と設定概要 |
| `charges.csv` | 最終要素電荷 |
| `mesh_triangles.csv` | 要素 geometry。`mesh_id` 列を含む |
| `mesh_sources.csv` | `mesh_id` ごとの元メッシュ種別、表面モデル、`epsilon_r`、要素数 |
| `mesh_potential.csv` | `write_mesh_potential=true` のとき |
| `charge_history.csv` | `history_stride > 0` のとき |
| `potential_history.csv` | `write_potential_history=true` かつ `history_stride > 0` のとき |
| `rng_state.txt` | 乱数状態 |
| `macro_residuals.csv` | マクロ粒子数の残差繰越 |
| `charge_ledger.csv` | species 別の signed charge flux、粒子数、再開用累積値 |

`mesh_potential.csv` は要素重心での電位 [V] を記録します。
自己項は `softening > 0` なら `1/softening`、そうでなければ面積等価半径近似を使います。
`periodic2` では explicit image shell を加えます。
`field_periodic_far_correction="m2l_root_oracle"` のときだけ exact Ewald residual も加えます。

`potential_history.csv` は `charge_history.csv` と同じ `history_stride` で要素ごとの電位を記録します。
形式は `batch, elem_idx, potential_V` です。
履歴ごとに `field_solver%refresh` と `compute_mesh_potential` が走るため、有効化すると計算コストが増えます。

`resume=true` の要件:

| 条件 | 内容 |
|---|---|
| 出力 | `write_files=true` が必須 |
| 読み込み元 | `restart_from` 未指定なら `output.dir`、指定時は `restart_from` |
| 必須ファイル | `summary.txt`, `charges.csv`, `rng_state.txt` |
| 任意ファイル | `macro_residuals.csv`。schema v2 で台帳 metadata がある場合は `charge_ledger.csv` も必須 |
| 挙動 | 必須 checkpoint がなければ新規実行にフォールバックせず停止 |

`restart_from` は checkpoint の読み込み元だけを変更します。
新しい出力は常に `output.dir` に書かれます。

MPI 実行時:

| ファイル | 内容 |
|---|---|
| `rng_state_rankNNNNN.txt` | rank 別乱数状態 |
| `macro_residuals_rankNNNNN.csv` | rank 別残差 |

`summary.txt` の `mpi_world_size` は現在の rank 数と一致している必要があります。schema v2 では model / ordered mesh / ordered species fingerprint も一致する必要があります。

`[[particles.species]].species_key` は restart fingerprint 用の安定 ID です。省略時は `species_<1-based index>` を割り当てます。明示する場合は粒子種間で一意にしてください。

---

## 高水準記法との関係

次のキーは schema には含まれますが、Fortran parser が実行時キーへ解決する高水準記法です。
Fortran parser は読み込み時に右列のキーとして扱います。

| 高水準キー | 解決先・用途 |
|---|---|
| `sim.box_origin`, `sim.box_size` | `sim.box_min`, `sim.box_max` |
| `inject_region_mode`, `uv_low`, `uv_high` | `inject_face` 上の `pos_low`, `pos_high` |
| `mesh.groups`, template の `group`, `center_local` | template ごとの実座標 |
| template の `placement_mode`, `anchor`, `offset`, `offset_frac` | `center` |
| template の `size_mode`, `size_frac` | `size_x`, `size_y`, `size`, `radius` など |

高水準記法の詳細、例、lint 時の扱いは [Configuration](Configuration.html) を参照してください。

---

## 検証ルール

| 項目 | ルール |
|---|---|
| 未知キー | すべてエラー |
| `[particles]` | `[[particles.species]]` のコンテナとしてのみ使用。直下に `key = value` は書かない |
| 旧キー | 旧名は未知キーとして扱う |
| 型 | schema と Fortran パーサの両方で検証 |
| 値域 | `beachx lint` と実行時 parser が既知制約を検証 |

実行前には次を推奨します。

```bash
beachx lint beach.toml
```
