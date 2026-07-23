title: BEACH 入力パラメータリファレンス

Lang: [日本語](Parameters.md) | [English](Parameters.en.md)

# 入力パラメータリファレンス

本文書は、Fortran実行系が読む`beach.toml`のパラメータリファレンスです。
単位は、特に断りがない限り SI 単位です。

初めて設定を組む場合は、先に[シミュレーションケースを設計する](ConfigurationRecipes.html)を読むと全体像を掴みやすいです。

box基準の座標・配置を指定する補助パラメータも通常のinput keyとして掲載し、どの値を計算または上書きするかを
[座標・配置の補助パラメータ](#座標配置の補助パラメータ)に明記しています。

| 関連ドキュメント | 内容 |
|---|---|
| [シミュレーションケースを設計する](ConfigurationRecipes.html) | 目的別の設定手順と調整ポイント |
| [`beach.toml`を作成・検証する](Configuration.html) | `beachx config`、schema、lint |
| [Algorithms](Algorithms.html) | BEM 場計算、粒子 push、衝突、蓄積電荷の計算手順への導線 |
| [Workflow](Workflow.html) | 実行、開発、テスト、KUDPC での注意 |
| [FMM](FMM.html) | `field_solver="fmm"`の選択と精度確認 |

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

Editor schema を使う場合は、`beach.toml` の先頭に GitHub Raw URL のコメント directive を置きます。
Fortran パーサは最初のセクションより前の通常キーを受け付けないため、`"$schema" = "..."` は使いません。

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

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

## 公式入門ケース

最初の実行には [10 分チュートリアル](Tutorial.html) と
`examples/tutorial_insulator.toml` を使ってください。`beachx config init`
も同一の設定を生成します。この入門ケースは x/y 軸の `periodic2` を
`field_periodic_image_layers=1`、`field_periodic_far_correction="none"` の有限画像和で扱います。

`reservoir_face`、無限周期補正、外部シースは、そのケースが実行できた後に
[シミュレーションケースを設計する](ConfigurationRecipes.html)から追加する応用設定です。

---

## TOML の階層とセクション一覧

`[sim]`、`[field]`、`[particles]`、`[mesh]`、`[periodic2]`、`[external_boundary]`、
`[output]` が通常の公開構成です。旧 `[outer_plasma]` と `[coupling]` は互換入力として残りますが、
`[external_boundary]` とは混在できません。ネスト関係は次のとおりです。

```text
beach.toml
├── [sim]
├── [field]
├── [particles]
│   └── [[particles.species]]       # 1 件以上の array-of-tables
├── [mesh]
│   ├── [mesh.groups.<name>]        # 名前付きの子 table
│   └── [[mesh.templates]]          # 0 件以上の array-of-tables
├── [periodic2]
├── [external_boundary]
│   ├── [external_boundary.field]
│   ├── [external_boundary.particles]
│   └── [external_boundary.ordinary_open]
├── [outer_plasma]                 # 旧互換・正規化後の実行時表現
├── [coupling]                     # 旧互換・正規化後の実行時表現
└── [output]
```

本文中の `sim.dt` や `external_boundary.field.model` は「table 名.key」の参照表記です。TOML ではそれぞれ
`[sim]` の下へ `dt = ...`、`[external_boundary.field]` の下へ `model = ...` と書きます。

| TOML table | 親 | 件数・必須条件 | 内容 |
|---|---|---|---|
| `[sim]` | root | 条件付き | 時間刻み、バッチ数、場ソルバ、境界、外部場、シース補正 |
| `[field]` | root | 任意 | 要素電荷の離散化 kernel |
| `[particles]` | root | 必須 | `[[particles.species]]` のコンテナ。直下に通常 key は置かない |
| `[[particles.species]]` | `[particles]` | 1 件以上 | 粒子種、注入方式、速度分布、マクロ粒子重み |
| `[mesh]` | root | 任意 | OBJ または組み込み template の選択 |
| `[mesh.groups.<name>]` | `[mesh]` | 0 件以上 | 複数 template で共有する配置と scale |
| `[[mesh.templates]]` | `[mesh]` | 0 件以上 | `mode="template"` で使う組み込み形状 |
| `[periodic2]` | root | 条件付き | split periodic2 の非零モード・零モード・下側境界モデル |
| `[external_boundary.field]` | `[external_boundary]` | facade 使用時は必須 | 外部場モデルと物理・診断パラメータ |
| `[external_boundary.particles]` | `[external_boundary]` | facade 使用時は必須 | z-high 粒子、流入、時間・orbit guard |
| `[external_boundary.ordinary_open]` | `[external_boundary]` | 任意 | outer が所有しない open 面。既定は `escape` |
| `[outer_plasma]` | root | 旧互換のみ | 正規化後の外部場・return設定 |
| `[coupling]` | root | 旧互換のみ | 正規化後のouter更新・粒子移送設定 |
| `[output]` | root | 任意 | 出力先、履歴、checkpoint 再開 |

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
| `direct` | 要素数が小さい場合の厳密な全対全評価、split reference | `free`、または条件を満たす`periodic2` split reference |
| `treecode` | 中規模以上の近似評価 | `field_bc_mode="free"` |
| `fmm` | 大規模評価、`periodic2`、FMM コア検証 | `field_bc_mode="free"` / `"periodic2"` |
| `auto` | point は direct / treecode、triangle P0 は direct / FMM を要素数で自動選択 | `field_bc_mode="free"` |

solver、kernel、場境界の組合せは[場の評価の互換表](FieldSolvers.html#solverと場境界の互換表)で確認できます。

共通キー:

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `softening` | float | `1.0e-6` | point kernel の softening 長さ [m]。`triangle_p0` では `0` が必須 |
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
| `softening` | float | `1.0e-6` | point source の特異性を緩和。`triangle_p0` は `0` |
| `field_normalization` | string | `"si"` | direct 評価前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `field_bc_mode` | string | `"free"` | 通常は`free`。`triangle_p0`のsplit referenceだけ`periodic2`を使用可能 |

`tree_theta`、`tree_leaf_max`、`tree_min_nelem` は `direct` では使いません。

##### `field_solver = "treecode"`

source octree を作り、遠方 node は monopole 近似、近傍 node は選択した source kernel の direct 和で評価します。
FMM のような local expansion は使わず、評価点ごとに木を走査します。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"treecode"` を指定 |
| `softening` | float | `1.0e-6` | point の近傍 direct 和と monopole 評価の softening。`triangle_p0` は `0` |
| `field_normalization` | string | `"si"` | tree 構築前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_theta` | float | `0.5` | MAC パラメータ。`0 < theta <= 1`。大きいほど速く粗い |
| `tree_leaf_max` | int | `16` | 葉 node あたり最大 source 数。`>= 1` |
| `field_bc_mode` | string | `"free"` | `treecode` では `"free"` のみ |

`tree_min_nelem` は `field_solver="auto"` 用のしきい値なので、明示 `treecode` では切替には使いません。

##### `field_solver = "fmm"`

simulator 非依存の Coulomb FMM コアを使います。
source 幾何の plan と、電荷更新ごとの state を分け、P2M/M2M/M2L/L2L/L2P と近傍 direct 和で評価します。
選択と精度確認は[FMM](FMM.html)、内部実装は
[Coulomb FMMコア詳細](FMMCore.html)にまとめています。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"fmm"` を指定 |
| `softening` | float | `1.0e-6` | point source の近傍 direct 和と FMM 評価に使用。`triangle_p0` は `0` |
| `field_normalization` | string | `"si"` | FMM plan 構築前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_theta` | float | `0.5` | near/far 判定の MAC パラメータ。`0 < theta <= 1` |
| `tree_leaf_max` | int | `16` | source tree の葉 node あたり最大 source 数。`>= 1` |
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | `periodic2` の近傍画像層数 |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / production `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | operator generator の Ewald 分解パラメータ |
| `field_periodic_ewald_layers` | int | `4` | Ewald generator の real-space / reciprocal cutoff 深さ |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | versioned periodic operator cache directory |
| `field_periodic_generation_tolerance` | float | `1e-8` | cache fingerprint に含める generation tolerance |

`field_periodic_*` は `field_bc_mode="periodic2"` のときだけ意味を持ちます。
`tree_min_nelem` は明示 `fmm` では使いません。

##### `field_solver = "auto"`

`element_kernel="point"` は要素数が `tree_min_nelem` 未満なら direct、以上なら treecode を使います。
`element_kernel="triangle_p0"` は同じしきい値で direct / FMM を選びます。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"auto"` を指定 |
| `softening` | float | `1.0e-6` | point の direct / treecode に使用。`triangle_p0` は `0` |
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

#### 場境界

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | `periodic2` の近傍画像層数 |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | operator generator の Ewald 分解パラメータ |
| `field_periodic_ewald_layers` | int | `4` | Ewald generator の outer shell / reciprocal cutoff 深さ |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | versioned periodic operator cache directory |
| `field_periodic_generation_tolerance` | float | `1e-8` | cache fingerprint に含める generation tolerance |

### `[external_boundary]`: 外部条件の公開設定

新しい設定ファイルで外部条件を指定する正本は、この`[external_boundary]` facadeです。
後段の`[outer_plasma]`、`[coupling]`、旧`[sim]` selectorは、既存設定の読解と
正規化後のruntime表現を確認するための互換リファレンスです。

外部条件は、内部実装の selector を組み合わせるのではなく、次の3つの責務で指定します。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 2.0e-3
thermal_voltage = 3.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
```

- `field` は外部 plasma 応答の電位・電場モデルを選ぶ。linear/kinetic は
  z-high interface 外、unified は rough surface から far 領域までを扱う。
- `particles` は z-high を通過する粒子の扱いと、流入 VDF の所有者を選ぶ。
- `ordinary_open` は outer model が所有しない open 面だけを扱う。

`field` と `particles` は facade を使うときの必須 table です。`ordinary_open` は省略でき、
その場合は `model="escape"` です。外部条件を使わない通常ケースでは `[external_boundary]`
自体を省略できます。

#### `[external_boundary.field]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `model` | string | 必須 | `none` / `linear_debye` / `kinetic_1d` / `unified_linear_response` |
| `kinetic_closure` | string | `absorbing_maxwellian` | `kinetic_1d` のみ。`absorbing_maxwellian` / `zhao_charge_driven` |
| `zhao_branch` | string | `auto` | `zhao_charge_driven` のみ。`auto` / `a` / `b` / `c` |
| `photoelectron_source_scale` | float | `1` | `zhao_charge_driven` の解析光電子 source 倍率。UV なしは `0` |
| `photoelectron_density_model` | string | `none` | `kinetic_1d + absorbing_maxwellian` の任意平均密度。`none` / `kinetic_mean` |
| `photoelectron_histogram_enabled` | bool | `false` | `linear_debye + same_batch` で histogram を有効にするときだけ `true` を指定。無効時は key 自体を省略 |
| `infinity_potential` | float | `0` | `linear_debye` だけで指定できる無限遠基準電位 [V] |
| `debye_length` | float | active model で必須 | `linear_debye` / `kinetic_1d` / `unified_linear_response` の長さ scale [m] |
| `thermal_voltage` | float | active model で必須 | `linear_debye` / `kinetic_1d` / `unified_linear_response` の電位 scale [V] |
| `unified_grid_points` | int | `129` | `unified_linear_response` の zero-mode Poisson grid 点数 |
| `accessible_fraction_tolerance` | float | `0.1` | `unified_linear_response` の rough-surface accessible fraction 収束許容値 |
| `max_linearity_ratio` | float | `0.25` | `linear_debye` / `unified_linear_response` の線形性診断上限 |
| `max_gap_ratio` | float | `5` | `linear_debye` / `kinetic_1d` の interface–mesh gap 診断上限 |
| `max_local_charge_ratio` | float | `50` | `linear_debye` / `kinetic_1d` の局所 plasma 電荷推定比上限 |
| `photoelectron_histogram_bins` | int | `32` | `linear_debye + same_batch` で histogram を有効にしたときの bin 数 |
| `photoelectron_histogram_energy_max` | float | histogram 時に必須 | 有効にした histogram の上端 [J] |
| `photoelectron_ambient_charge_scale` | float | histogram 時に必須 | 有効にした histogram の ambient signed-charge scale [C] |
| `max_photoelectron_charge_ratio` | float | `0.1` | 有効にした histogram の光電子電荷比上限 |

`interface_z` と粒子 return model は指定しません。`interface_z` は
`sim.box_max[2]` から、return model は `field.model` と `particles.mode` から導出されます。
`model="none"` では `model` 以外の field key を書けません。

まず model/closure に対応する行だけを使い、その他の条件付き key は必要な診断や機能を有効にするときだけ追加します。
選んだ行にない key は、既定値と同じ値を書いても no-op として拒否されます。

| 選択 | 必須・通常設定する key | 条件付きで追加する key |
|---|---|---|
| `none` | `model` のみ | 追加不可 |
| `linear_debye` | `debye_length`, `thermal_voltage` | `infinity_potential`、linearity / gap / local-charge 上限、`same_batch` histogram |
| `kinetic_1d + absorbing_maxwellian` | `debye_length`, `thermal_voltage` | `kinetic_closure`、`photoelectron_density_model`、gap / local-charge 上限 |
| `kinetic_1d + zhao_charge_driven` | `debye_length`, `thermal_voltage`, `kinetic_closure`、必要な `sim.sheath_*` 物理値 | source scale / branch、gap / local-charge 上限 |
| `unified_linear_response` | `debye_length`, `thermal_voltage` | grid、accessible fraction、linearity 上限 |

`kinetic_closure="zhao_charge_driven"` でも `debye_length` と `thermal_voltage` は schema 上必須ですが、
Zhao root/profile の物理 scale ではありません。光電子ありでは $T_{pe}$ と基準密度、光電子なしでは
ambient $T_e$ と $n_\infty$ から profile scale を導出します。

この 2 値は split-interface の gap、lateral field、local-charge 診断の基準量です。
`zhao_branch="auto"`、zero gauge、density model なしなど、
既定値または model が固定する値は通常は書きません。

#### `[external_boundary.particles]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `mode` | string | 必須 | `local_source` / `same_batch` / `zhao_queue` |
| `inflow_model` | string | `auto` | `auto` / `source_vdf` / `infinity_barrier` / `legacy_sheath` |
| `legacy_sheath_model` | string | 条件付き | `legacy_sheath` 時の `floating_no_photo` / `zhao_auto` / `zhao_a` / `zhao_b` / `zhao_c` |
| `steady_start_mode` | string | `none` | `kinetic_1d + zhao_charge_driven + same_batch` で `zhao_floating` を有効にするときだけ指定 |
| `steady_start_mesh_id` | int | `1` | `steady_start_mode="zhao_floating"` で使う mesh ID |
| `outer_update_stride` | int | `1` | `local_source` / `same_batch` の `linear_debye` / `kinetic_1d` だけで指定できる更新間隔 [batch] |
| `field_evolution_timescale` | float | `0` | `same_batch` / `zhao_queue` の frozen-field 診断時間 scale [s] |
| `max_frozen_field_ratio` | float | `0.1` | `same_batch` / `zhao_queue` の frozen-field 適用性比上限 |
| `outer_orbit_dt` | float | `0` | unified 3D `same_batch` の外部軌道刻み。正値を明示 |
| `outer_orbit_max_steps` | int | `100000` | unified 3D `same_batch` の外部軌道最大 step 数 |
| `outer_orbit_energy_tolerance` | float | `1e-4` | unified 3D `same_batch` の相対 energy 誤差上限 |

`mode` は z-high を出た粒子を外部領域へ移送するか、その結果をいつ反映するかだけを選びます。
reservoir 流入の分布や補正は `inflow_model` で独立に選びます。

| `mode` | 意味 | 対応する field |
|---|---|---|
| `local_source` | 外部軌道へ移送せず、z-high 通過を `ordinary_open` で処理 | すべて |
| `same_batch` | z-high 通過粒子を同じ batch 内で return / escape 判定 | `linear_debye` / `kinetic_1d` / `unified_linear_response` |
| `zhao_queue` | Zhao reservoir queue に保持し、後続 batch で再注入 | `kinetic_1d` + `zhao_charge_driven` + `zhao_branch="auto"` |

mode ごとの追加必須値は次のとおりです。表にない内部 return / transfer ID は書きません。

| `particles.mode` と field | 追加必須・制約 |
|---|---|
| `local_source` + `none` / unified | transport / time / orbit key は追加しない |
| `local_source` + linear / kinetic | 必要なら `outer_update_stride` だけを追加 |
| `same_batch` + linear / kinetic | `field_evolution_timescale > 0`、`inflow_model="auto"`。必要なら update stride / time guard |
| `same_batch` + unified | `field_evolution_timescale > 0`、`outer_orbit_dt > 0`。必要なら time / orbit guard |
| `zhao_queue` | `sim.batch_duration > 0`、`field_evolution_timescale > 0`、正の photoelectron source。update stride は内部で 1 に固定 |

`steady_start_*` は Zhao `same_batch`、`outer_orbit_*` は unified `same_batch`、
time guard は `same_batch` / `zhao_queue` にだけ書けます。mode に効果のない key はエラーです。

`inflow_model="auto"` は、tracked `linear_debye` / `kinetic_1d` では同じ1D profileへ
流入を委ね、それ以外では `source_vdf` として解決します。この1D tracked構成では
別の inflow correction を重ねられません。`unified_linear_response` は外部軌道を所有しますが
流入 VDF は所有しないため、`source_vdf`、`infinity_barrier`、`legacy_sheath`を選べます。

#### `[external_boundary.ordinary_open]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `model` | string | `escape` | `escape` / `potential_barrier` |

`local_source` では z-high も `ordinary_open` が処理します。`same_batch` / `zhao_queue` では
outer model が z-high を所有し、`ordinary_open` は残りの open 面だけを処理します。流入の
`particles.inflow_model="infinity_barrier"` と、流出の
`ordinary_open.model="potential_barrier"` は別々に選べます。

facade は読み込み時に旧 `[sim]` / `[outer_plasma]` / `[coupling]` の実行時表現へ正規化されます。
両方の書式を同じファイルに混在させるとエラーです。選択手順と制約は
[外部境界を構成する](OuterPlasmaModels.html)、旧書式との対応は
[境界設定の移行ガイド](BoundaryConfigurationMigration.html)にまとめています。

### `[periodic2]`: 非零モード・零モード・下側境界

`[periodic2]` は `[sim]` の子ではなく、トップレベル table です。
legacy `periodic2` では `field_solver="fmm"` を使います。小規模検証用の split reference に限り、
`field_solver="direct"` と `[periodic2]`、`[external_boundary]` を明示します。

| キー | 既定値 | 意味 |
|---|---:|---|
| `nonzero_mode_backend` | 必須 | `panel_spectral_reference` / `cached_kneq0` |
| `zero_mode_policy` | 必須 | `exclude_k0` |
| `lower_boundary_model` | 必須 | `symmetric_vacuum` / `e_bottom_zero` |
| `reference_mode_layers` | `4` | Fourier mode cutoff |
| `panel_quadrature_order` | `12` | panel 面積積分次数 |
| `interface_sample_n` | `5` | interface 各軸の診断点数 |
| `interface_phi_tolerance` | `1e-3` | 非零モード電位比上限 |
| `interface_field_tolerance` | `1e-3` | 非零モード電場比上限 |

### `[outer_plasma]`: 互換・正規化後 runtime リファレンス

> **新規設定用の API ではありません。** 新規ファイルでは `[external_boundary.field]` を使ってください。
> 以下の raw 名は旧ファイルと正規化後の runtime contract を読むためのもので、
> `[external_boundary]` と同じファイルへ追加できません。

facade では `interface_z` と `return_model` を自動導出します。旧 `[outer_plasma]` はトップレベル table で、
`[periodic2]` の子ではありません。

| キー | 既定値 | 意味 |
|---|---:|---|
| `model` | 必須 | `linear_debye` / `kinetic_1d` / `unified_linear_response` |
| `interface_z` | 必須 | z 上側 interface。初期モデルでは box 上面 |
| `infinity_potential` | `0` | 無限遠基準電位 [V] |
| `debye_length` | 必須 | 線形/`absorbing_maxwellian` tailとsplit診断の長さscale。Zhao profileの物理scaleではない |
| `thermal_voltage` | 必須 | 線形性とsplit診断の電位scale。Zhao profileの物理scaleではない |
| `unified_grid_points` | `129` | unified zero-mode Poisson grid 点数（17 以上） |
| `accessible_fraction_tolerance` | `0.1` | rough surface 高さ標本を各軸 2 倍にしたときの accessible fraction 最大差 |
| `max_linearity_ratio` | `0.25` | `abs(phi_t-phi_inf)/thermal_voltage` 上限 |
| `max_gap_ratio` | `5` | `(z_t-z_mesh,max)/lambda` 上限 |
| `max_local_charge_ratio` | `50` | 局所平均 plasma 電荷推定比上限 |
| `kinetic_closure` | `absorbing_maxwellian` | `kinetic_1d` の density/VDF closure。`absorbing_maxwellian` / `zhao_charge_driven` |
| `zhao_branch` | `auto` | Zhao closure の branch。`auto` / `a` / `b` / `c` |
| `photoelectron_source_scale` | `1` | Zhao解析光電子sourceの独立倍率。UVなしは`0`。queue occupancyの$\eta$とは別物 |
| `photoelectron_density_model` | `none` | `none` / `kinetic_mean`。後者は `kinetic_1d` へ平均光電子密度を追加 |
| `photoelectron_histogram_enabled` | `false` | z-high を外向き通過する光電子の histogram と適用性検査を有効化 |
| `photoelectron_histogram_bins` | `32` | 法線運動エネルギー histogram の bin 数 |
| `photoelectron_histogram_energy_max` | histogram 有効時に必須 | histogram 上端 [J]。正値必須 |
| `photoelectron_ambient_charge_scale` | histogram 有効時に必須 | 線形モデルの適用性を比較する ambient signed-charge scale [C] |
| `max_photoelectron_charge_ratio` | `0.1` | `abs(Q_pe,batch)/Q_ambient,scale` 上限 |
| `return_model` | `none` | 1D 解析 return または unified 3D 明示軌道の ID |

#### 正規化後の `kinetic_1d` contract（標準・推奨）

公開 facade で `external_boundary.field.model="kinetic_1d"` を選ぶと、runtime の
`outer_plasma.model` も `kinetic_1d` に正規化されます。production では `cached_kneq0` と組み合わせます。
z-high に定義した負電荷と正電荷の `reservoir_face` species を、それぞれ無限遠の electron VDF と ion VDF
として使います。

| 項目 | 仕様 |
| --- | --- |
| gauge | `phi(infinity)=0`。`infinity_potential`の非ゼロ値を拒否 |
| far boundary | `absorbing_maxwellian`は`debye_length`のRobin tail。光電子ありZhao instantは$T_{pe}$と$n_{ref}$から導出した$\lambda_{D,pe}$、no-photo Zhaoはambient $T_e,n_\infty$から導出した$\lambda_{D,e}$、queueは$L=10\lambda_{D,pe}$の有限reservoir境界を使う |
| closure | 既定の `absorbing_maxwellian`、または蓄積電荷が決める interface 電場を保つ `zhao_charge_driven` |
| supported branch | `absorbing_maxwellian` は単調 branch。`zhao_charge_driven` は非単調 Type A を含む Zhao A/B/C |
| unsupported | `absorbing_maxwellian` では virtual cathode、trapped population、sub-Bohm inflow |
| nonlinear solve | 解析bordered-tridiagonal Jacobian + branch-preserving Newton |
| recovery | pseudo-transientとinterface-field continuation。元のPoisson residual合格後だけ受理 |
| fallback | 別sheath modelや前回解へfallbackしない |

公開設定の `external_boundary.particles.mode="same_batch"` は、
`return_model="kinetic_1d_profile_return"` と
`particle_transfer_mode="electrostatic_1d_instant_return"`、`outer_queue_enabled=false` へ正規化されます。
これらの raw ID は入力せず、解決後は次の instant 写像を表します。

1. 更新済みprofileの`phi_interface-phi_infinity`で無限遠VDFをinterfaceへ写像する。
2. 同じ離散profileとRobin tailでescapeまたはturning pointを判定する。
3. turning粒子について解析往復時間後に相当する復帰状態を作り、同じsimulation時刻・batchでinterfaceへ戻す。

instant写像の適用範囲:

- 対象は定常・準定常 sheath です。定常化後の平均電流と離脱力に使用できます。
- UV 照射開始時など、遅延した return current が効く過渡応答は表しません。
- 準定常条件は `tau_outer/field_evolution_timescale` で制限します。
- `tau_outer/batch_duration >= 1` では、batch 履歴を return current の物理時間履歴として解釈できません。

UVなしでは`photoelectron_source_scale=0`を指定します。この経路はenabledな`photo_raycast` species、
`sheath_photoelectron_ref_density_cm3`、$T_{pe}$を要求しません。z-highのambient electron/ionが与える準中性領域の
$n_\infty,T_e,u_e,u_i$から、入射electron reservoirの正規化、cutoff、速度写像を同じZhao VDFで解きます。
$E_I=0$は平坦なType B/C接続点、$E_I<0$はType Cです。電流ゼロは各batchのroot条件にせず診断として残すため、
表面電荷が定常化したときに旧no-photo Type Cのfloating rootへ近づくかを検証できます。

公開設定の `external_boundary.particles.mode="zhao_queue"` は runtime の `outer_queue_enabled=true` へ正規化され、
強 UV 立上がりなどの outer flight delay を batch 履歴へ反映します。この mode は
`external_boundary.field.kinetic_closure="zhao_charge_driven"`、`sim.batch_duration` または
`dt * batch_duration_step` から解決した正の `batch_duration` を要求します。

update stride は runtime で 1、histogram は無効に固定されるため、公開入力にはどちらの key も書きません。

1. batch開始時にdue eventをrank-local queueからpopする。
2. 残ったglobal photoelectron inventoryを水平面積で割り、$0\le z\le10\lambda_{D,pe}$の有限columnに一致する
   Zhao population scale $\eta$とprofileを更新する。
3. fresh sourceとdue returnを進め、外向き粒子は$L$より手前のreturnまたは$L$でのreservoir escapeを解き、
   batch中央を通過時刻として$t_{due}=t_{mid}+\tau_{outer}$でenqueueする。

出力名`outer_photoelectron_population_fraction`の$\eta$は確率ではなく定常reference populationに対するoccupancy scaleです。
$\eta=0$から連結するpathを$0\le\eta\le16$で探索し、$\eta>1$を許します。

clamp、targetを無視したfull-population解、
disconnected branchへのjump/fallbackは行いません。queue modeは`zhao_branch="auto"`を要求し、縮退条件を満たす連続branch
遷移とcolumnが$\eta$とともに単調増加するpathだけを許します。foldを含むpathやtargetに到達できないpathでは停止します。

eventはbatch開始時だけreleaseされ、enqueue後にouter fieldが変わってもterminal状態を再積分しません。queue modeは
$L$外のRobin tailを使わず、$L$到達をreservoirへの吸収/escapeとします。各eventでは`tau_outer`に次のbatch-start pollまでの
量子化遅延とmidpoint crossing時刻の誤差上限`batch_duration/2`を加えた時間が
`max_frozen_field_ratio * field_evolution_timescale`以下でなければ停止します。
`batch_duration`自体にも同じ上限を設定validationで課します。

`batch_duration`、tracked粒子数、水平面積、有効interface位置、profile gridの収束を確認します。

組合せ制約:

- `reservoir_potential_model`、Zhao 系 `sheath_injection_model`、`b0 != 0` との併用を拒否します。
- `kinetic_closure="zhao_charge_driven"` は `model="kinetic_1d"`、`infinity_potential=0`、
  `photoelectron_density_model="none"` を要求します。legacy `sheath_injection_model` や
  `reservoir_potential_model` と重ねません。
- `zhao_branch="a"` / `"b"` / `"c"` は `zhao_charge_driven` でのみ指定できます。`auto` は利用可能な branch を探索します。
- `zhao_charge_driven` は準中性なambient electron/ion density、正の内向きelectron/ion driftを要求し、
  `sheath_reference_coordinate`を拒否します。
- `photoelectron_source_scale>0`では正の`sheath_photoelectron_ref_density_cm3`と負電荷`photo_raycast` speciesを
  ちょうど1つ要求します。`photoelectron_source_scale=0`ではenabledな`photo_raycast` speciesを拒否し、queueも使えません。
- `sheath_electron_drift_mode="normal"`と`sheath_ion_drift_mode="normal"`だけを受理します。
- Zhaoに使う負電荷`photo_raycast` speciesは`normal_drift_speed=0`、ion温度はcold-ion近似$T_i\le0.1T_e$を要求します。
- 光電子ありZhao profileの電位/長さscaleは$T_{pe}$と基準密度から導出した$\lambda_{D,pe}$です。
  no-photoでは同じ式の正規化をambient $T_e,n_\infty$へ切り替え、$\lambda_{D,e}$をsummaryへ出力します。
- `debye_length`と`thermal_voltage`はZhao root/profileを変えません。ただしsplit-interfaceの`interface_eta_gap`、
  lateral phi/field、local-charge診断のreference inputとして現時点でも必要です。
- tracked `photo_raycast.emit_current_density_a_m2`は、$T_{pe}$をJへ換算し、$v_{th,pe}=\sqrt{2T_{pe}/m_{pe}}$として
  `photoelectron_source_scale * |q_{pe}|n_{ref}\sin(\alpha)v_{th,pe}/(2\sqrt{\pi})`と1%以内で一致する必要があります。
- analytic raw currentはtracked sourceの整合性検査とcurrent-density診断に使いますが、root、surface charge、ledgerへ
  別途加えません。tracked放出と再吸収だけが後二者を更新し、$\eta$はcurrent診断のraw photoelectron emission-current項を
  scaleしません。
- 初版は有効平面近似です。`ray_direction`やrough surfaceから到達するVDFをZhao outer populationへ自己無撞着に
  接続せず、`ray_direction`と`sheath_alpha_deg`はそれぞれ照射rayによる放出面samplingと解析sourceを独立に指定します。
- Zhaoの収束はprofile grid、有効interface位置、tracked粒子数、`dt`/batch解像度で調べます。
  `debye_length`や`thermal_voltage`の変更をprofile収束試験として解釈しません。
- `photoelectron_density_model="kinetic_mean"` は既定の `absorbing_maxwellian` で、先頭の負電荷
  `photo_raycast` species から half-Maxwellian flux を作ります。
- 平均密度モデルが供給するのは outer profile だけです。表面電荷は追跡粒子が更新するため、統計的な return current を重ねません。
- tracked return を使う全 species で `deposit_opposite_charge_on_emit=true` が必要です。

詳細は[粒子の escape と return](ParticleEscapeReturn.html)を参照してください。

既定 closure の例は
[`periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml)です。no-photo Zhao例は
[`periodic2_zhao_no_photo_outer.toml`](../examples/periodic2_zhao_no_photo_outer.toml)、UVありの定常warm-start例は
[`periodic2_zhao_charge_driven_outer.toml`](../examples/periodic2_zhao_charge_driven_outer.toml)です。

過渡queueの例は[`periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml)です。
これは記載した物理timescaleで
長いflightを拒否するexpected-fail guard fixtureであり、成功した物理検証例ではありません。

物理モデルの前提は
[ADR 0001](adr/0001-kinetic-outer-plasma.md)と[ADR 0003](adr/0003-zhao-charge-driven-outer-closure.md)にあります。

#### 正規化後の `unified_linear_response` contract（高度・限定用途）

`unified_linear_response`は`kinetic_1d`の上位互換ではありません。roughnessとplasma responseが同じ領域に重なり、
split windowを置けず、線形性gateを満たす場合だけ、rough surfaceの線形screeningとして選びます。

| 項目 | 仕様 |
| --- | --- |
| zero mode | 表面投影から遠方まで一つのPoisson grid |
| rough surface | plasma-accessible areaと線形mean-plasma電荷を含む |
| nonzero mode | 表面最高点直上から$\sqrt{k^2+\kappa^2}$ tailへ接続 |
| `interface_z` | field切断面ではなく粒子ownership面 |
| geometry/kernel | single-valued height fieldと`triangle_p0`が必須 |
| photoelectron mean | `photoelectron_density_model="none"`のみ |
| particle transfer | `none`または`electrostatic_3d_explicit_orbit` |
| 3D orbit | `b0=0`、固定刻み、energy/frozen-field errorの上限 |
| applicability | 線形性上限を超えたらfallbackせず停止 |

数値・適用性の規則:

- `unified_grid_points >= 17` が必要で、既定値は `129` です。
- `accessible_fraction_tolerance` は、両周期軸の高さ標本を 2 倍にしたときの accessible fraction の最大差を制限します。
- refinement 後の標本を solve に使い、許容差違反時は初期化で停止します。
- production study では、報告する物理量について grid refinement を確認します。

検証例は `examples/periodic2_unified_linear_response.toml`、明示粒子移送の例は
`examples/periodic2_unified_explicit_orbit.toml`、詳細は `docs/adr/0002-unified-periodic-outer-domain.md` です。

### `[coupling]`: 互換・正規化後 runtime リファレンス

> **新規設定用の API ではありません。** 新規ファイルでは `[external_boundary.particles]` を使ってください。
> 以下の raw 名は旧ファイルと正規化後の runtime contract を読むためのものです。

facade は `update_mode`、`particle_transfer_mode`、`outer_queue_enabled` を model と mode から導出します。

| キー | 型 | 既定値 | 説明 |
| --- | --- | ---: | --- |
| `update_mode` | string | `"explicit"` | 現在は`explicit`のみ。outer profileを明示的な更新点で再計算 |
| `particle_transfer_mode` | string | `"none"` | facade が導出する `none` / `electrostatic_1d_instant_return` / `electrostatic_3d_explicit_orbit` |
| `steady_start_mode` | string | `"none"` | `none` / `zhao_floating`。新規実行を Zhao 零電流定常根と対応する平面電荷から開始 |
| `steady_start_mesh_id` | int | `1` | `zhao_floating`で初期電荷を面積比で配る水平平面の`mesh_id` |
| `outer_update_stride` | int | `1` | outer profile更新batch間隔 |
| `field_evolution_timescale` | float | `0` | frozen-field比較時間 [s]。1D returnでは正値必須 |
| `max_frozen_field_ratio` | float | `0.1` | instantでは`tau_outer`、queueでは`tau_outer`、次のpollまでの遅延、midpoint時刻誤差上限の合計を`field_evolution_timescale`で割った上限。queueの`batch_duration`にも適用 |
| `outer_orbit_dt` | float | `0` | 3D outer orbit固定刻み [s]。3D modeでは正値必須 |
| `outer_orbit_max_steps` | int | `100000` | 3D outer orbit step上限。到達時はdiscardせず停止 |
| `outer_orbit_energy_tolerance` | float | `1e-4` | 3D outer orbit全エネルギー相対誤差上限 |
| `outer_queue_enabled` | bool | `false` | 対応するZhao構成でouter flightをbatch間queueへ保存し、queued photoelectron columnで過渡closureを解く |

facade が解決する return / transfer の runtime 対:

| `outer_plasma.model` | `outer_plasma.return_model` | `coupling.particle_transfer_mode` |
| --- | --- | --- |
| `linear_debye` | `electrostatic_1d_instant_return` | `electrostatic_1d_instant_return` |
| `kinetic_1d` | `kinetic_1d_profile_return` | `electrostatic_1d_instant_return` |
| `unified_linear_response` | `electrostatic_3d_explicit_orbit` | `electrostatic_3d_explicit_orbit` |

`kinetic_1d_profile_return`ではreturnとtransferの文字列は同一ではありません。1D transferを有効にした
`linear_debye`と`kinetic_1d`は同じprofileで流入も所有するため、正規化後の
`reservoir_potential_model`と`sheath_injection_model`はともに`none`です。
`kinetic_1d`と`unified_linear_response`の`infinity_potential`は0固定です。

定常 warm start:

公開 facade では `external_boundary.particles.mode="same_batch"`、
`steady_start_mode="zhao_floating"`、必要なら `steady_start_mesh_id` を指定します。
raw の transfer ID と `outer_queue_enabled` は facade が導出します。

`zhao_floating`は、新規実行の最初batch前に設定済み無限遠reservoirとUV sourceから Zhao 零電流定常根を解き、
その根から`phi(infinity)=0`のkinetic profileを構築します。水平断面積を$A$、定常根のinterface電場を$E_I$
とすると、選択平面へ与える総電荷は次です。

$$
Q_{seed}=
\begin{cases}
2\epsilon_0 A E_I, & \texttt{symmetric_vacuum},\\
\epsilon_0 A E_I, & \texttt{e_bottom_zero}.
\end{cases}
$$

電荷は`steady_start_mesh_id`の三角形へ面積比で配り、他のmeshは電荷0のまま開始します。従ってplane + sphereで
`steady_start_mesh_id=1`がplaneなら、planeだけをseedしsphereは中性です。最初のouter refresh、無限遠reservoirからの流入補正、
interfaceを出た粒子のinstant return / escapeは、この同じprofileを使います。analytic currentを表面電荷へ追加せず、
後続の電荷更新は通常どおりtracked粒子だけが行います。

このmodeは未帯電状態からの物理過渡を時間積分せず、定常branch上の初期条件を明示的に与えます。定常・準定常量の
評価用であり、UV立上がりやreturn-current遅延の主張には使いません。後者のためのqueue過渡closureは別modeとして残ります。

公開 facade の `zhao_floating` は次を要求します。

- `external_boundary.field.model="kinetic_1d"` +
  `kinetic_closure="zhao_charge_driven"` + `external_boundary.particles.mode="same_batch"`
- `zero_mode_policy="exclude_k0"` と対応する lower boundary model
- 新規実行では全meshの初期電荷0。`output.resume=true`ではcheckpointのmesh電荷とouter stateを復元し、再seedしない
- `mesh.mode="template"`で、選択meshが水平・同一高さの1平面としてperiodic cellの$A$を覆い、outer interfaceより下にあること

物理的な一意性や安定性はwarm startだけでは確立しません。publication用には、独立に緩和させた状態または摂動した
seedから同じ定常観測量へ返るかを確認します。

粒子移送の規則:

- 正規化後の `outer_plasma.return_model` と `coupling.particle_transfer_mode` は上表の対応する対になります。
- 1D transferはopenなz-high interface、x/y周期wrap、`b0=0`だけに対応します。
- `kinetic_1d`はenabledな負・正z-high `reservoir_face` speciesをそれぞれちょうど1つ要求します。
- instant modeは正の`field_evolution_timescale`を要求し、`max_frozen_field_ratio`を適用性上限に使います。
- 公開 facade の `zhao_queue` は `kinetic_1d` + `zhao_charge_driven` + `zhao_branch="auto"`、
  直接指定または`dt * batch_duration_step`から解決した正の`batch_duration`を要求し、
  runtime の `outer_update_stride` を 1 に固定します。
  各eventでは`tau_outer`、次のbatch-start pollまでの遅延、midpoint crossing時刻誤差上限の合計に
  `max_frozen_field_ratio * field_evolution_timescale`の上限を課し、`batch_duration`にも同じ上限を要求します。
- queue modeではphotoelectron histogramを公開設定できません。3D explicit orbitのpersistent queueは未実装です。

### periodic2 と外部境界の組合せ制約

関連exampleは目的別に選びます。

| 目的 | example | 主要制約 |
| --- | --- | --- |
| split linear reference | `periodic2_linear_outer_reference.toml` | 閾値違反時にfallbackしない |
| 1D instant return | `periodic2_outer_particle_transfer.toml` | z-high、`b0=0`、x/y periodic |
| 強UV Zhao過渡queue | `periodic2_zhao_transient_outer.toml` | expected-failのflight/batch time-scale guard、stride 1、有限$10\lambda_{D,pe}$領域、rank-local checkpoint |
| unified 3D orbit | `periodic2_unified_explicit_orbit.toml` | 全3D field、固定刻みouter orbit |
| tracked photoelectron return | `periodic2_photoelectron_return.toml` | instant return、放出元逆符号電荷、outgoing histogram |

光電子 histogram の規則:

- 全 MPI rank を集約し、前 batch と累積の値を checkpoint に保存します。
- histogram は診断と適用性検査だけを担当します。return / escape は `field.model` と `particles.mode` から解決されます。
- 公開 facade では `linear_debye + same_batch` の場合だけ有効化できます。
- `photoelectron_density_model="kinetic_mean"` とは、必要な return model が異なるため併用できません。
- z-high outward crossing の signed charge が適用性上限を超えると停止します。
- tracked outer transfer を使う全 `photo_raycast` species で `deposit_opposite_charge_on_emit=true` が必要です。

periodic2では、`sim.use_box=true`とちょうど2つのperiodic軸が必要です。outer transferを使う構成では、
その2軸をx/yとし、z-highをopen interfaceにします。同じ周期条件をfield、collision、`photo_raycast`に適用します。

| far correction | 意味 |
| --- | --- |
| `none` | 有限画像近似 |
| `auto` | 互換用。現在は`none` |
| `cached_kneq0` | versioned operatorを再利用するproduction非零モード |

`field_periodic_far_correction="auto"` は互換用に受理され、現在は `none` と同じ扱いです。
無限周期のproduction計算では、`cached_kneq0`を明示的に選択してください。
削除済みの`m2l_root_oracle`を指定すると、`cached_kneq0`を案内して起動時にrejectします。

`cached_kneq0`では`exclude_k0` providerが物理的$k=0$を1回加えます。
`symmetric_vacuum`は上下を$\pm Q/(2\epsilon_0A)$、`e_bottom_zero`は下側0・上側$Q/(\epsilon_0A)$とします。

### `[sim]`: 外部場と外部境界が参照する物理値

`phi_infty` と `sheath_*` はモデルが参照する物理値なので、新しい
`[external_boundary]` と併用します。一方、`reservoir_potential_model`、
`open_boundary_model`、`sheath_injection_model` は旧互換 selector です。
facade 使用時はこれら3つを書かず、`external_boundary.particles` と
`external_boundary.ordinary_open` から導出させます。

以下のraw selector説明は、旧設定と正規化後のruntime状態を読むために残しており、
新規設定向けの別APIではありません。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | 一様外部電場 [V/m] |
| `e0_abs` | float | 未指定 | 一様外部電場の大きさ [V/m] |
| `e0_phi_xy_deg` | float | `0.0` | `e0_abs` 指定時の xy 面内方位角 [deg] |
| `e0_phi_z_deg` | float | `0.0` | `e0_abs` 指定時の xy 面からの仰角 [deg] |
| `b0` | float[3] | `[0,0,0]` | 一様磁場 [T] |

一様外部電場は、`e0 = [Ex, Ey, Ez]` で直接指定するか、`e0_abs` と角度で指定します。
両形式の混在はエラーです。

#### 物理値と旧selectorの互換リファレンス

表中で「旧selector」と示した3行は互換入力専用です。それ以外はfacadeからも参照される
物理値または共通の境界処理値です。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `reservoir_potential_model` | string | `"none"` | 旧 selector。`none` / `infinity_barrier` |
| `phi_infty` | float | `0.0` | 無限遠基準電位 [V] |
| `open_boundary_model` | string | `"escape"` | 旧 selector。`escape` / `potential_barrier` |
| `multiple_box_events_policy` | string | `"abort"` | `abort` / `soft_discard`。1 step の box event 上限超過時の処理 |
| `multiple_box_events_soft_discard_count_limit` | int | `1000` | 累積 soft discard 件数がこの値を超えると停止 |
| `multiple_box_events_soft_discard_abs_charge_limit` | float | `1e-12` | 累積 soft discard 絶対 macro charge [C] がこの値を超えると停止 |
| `injection_face_phi_grid_n` | int | `3` | 注入面平均電位の `N x N` 評価格子 |
| `raycast_max_bounce` | int | `16` | `photo_raycast` の最大反射回数 |
| `sheath_injection_model` | string | `"none"` | 旧 selector。`none` / `zhao_auto` / `zhao_a` / `zhao_b` / `zhao_c` / `floating_no_photo` |
| `sheath_alpha_deg` | float | `60.0` | Zhao シースの太陽高度角 [deg] |
| `sheath_photoelectron_ref_density_cm3` | float | `64.0` | Zhao シースの基準光電子密度 [cm^-3] |
| `sheath_reference_coordinate` | float | 未指定 | シース 1D 座標の基準平面位置 [m] |
| `sheath_electron_drift_mode` | string | `"normal"` | `normal` / `full` |
| `sheath_ion_drift_mode` | string | `"normal"` | `normal` / `full` |

互換runtimeでの組合せと詳細:

- 公開facadeでは`particles.inflow_model="legacy_sheath"`と`"infinity_barrier"`を同時に選べません。
  正規化後のruntimeでは、これは`sheath_injection_model != "none"`と
  `reservoir_potential_model="none"`の組として表れます。
- legacy sheathの値と挙動は[シース流入補正](SheathInjectionClosures.html)にあります。
- 流束と速度は[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)を確認してください。
- シース補正と反射・return は[シース流入補正](SheathInjectionClosures.html)と
  [粒子の escape と return](ParticleEscapeReturn.html)に分けています。

`external_boundary.particles.inflow_model="infinity_barrier"`の評価
（互換runtimeでは`reservoir_potential_model="infinity_barrier"`）:

- 各 batch 冒頭で更新した電場・電位から、注入面平均電位を評価します。
- point / `triangle_p0` kernel、periodic2、zero mode、outer profile、`e0` は粒子運動と同じ規約で含めます。
- 同じ `N x N` 格子で母標準偏差・最小・最大も集計します。この診断による追加の電位評価はありません。
- Maxwellian reservoir で `abs(q_particle) * phi_std > 0.1 * (k_B*T + 0.5*m_particle*u_normal^2)` なら、
  MPI root が初回と最終 batch に面平均近似の警告を出します。

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

`external_boundary.ordinary_open.model="potential_barrier"`の判定
（互換runtimeでは`open_boundary_model="potential_barrier"`）:

1. 境界通過点の BEM 電位 `phi_boundary` を、粒子運動と同じ snapshot 規約で評価します。`e0` の局所電位も含みます。
2. 電位障壁 `q_particle * (phi_infty - phi_boundary)` を計算します。
3. 障壁が正で `0.5 * m_particle * v_normal^2` より大きければ、法線速度を反転します。それ以外は脱出です。

一様電場には有限な無限遠電位がありません。`e0` と併用する場合、`phi_infty` を有効な reservoir 基準として整合させます。

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
現状、`velocity_distribution="grid"`は
`external_boundary.particles.inflow_model!="legacy_sheath"`のときのみ有効です。
正規化後のruntimeでは`sim.sheath_injection_model="none"`に対応します。

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

生成した光電子は常に`w_hit`を重みに使い、通常粒子として追跡します。表面へ戻れば通常の衝突として吸収し、
open面では`external_boundary.ordinary_open.model`、z-high interfaceでは
`external_boundary.particles.mode`から解決した輸送がreturn / escapeを決めます。

---

### 旧互換runtime selector `sim.sheath_injection_model`

新規設定では`external_boundary.particles.inflow_model="legacy_sheath"`と
`legacy_sheath_model`を使います。以下は、旧設定と正規化後のruntime値を解釈するための説明です。
`sim.sheath_injection_model`は既存の`reservoir_face` / `photo_raycast` speciesを束ね、
シースに対応する流束や法線速度cutoffを上書きします。

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

#### `[[mesh.templates]]`: 組み込み形状

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

##### `kind = "plane"`

XY 平面上の長方形を `nx * ny` 個の矩形セルに分け、各セルを 2 三角形へ分割します。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | 平面中心 `[x, y, z]` [m] |
| `size_x` | float | `1.0` | x 方向サイズ [m]。`> 0` |
| `size_y` | float | `1.0` | y 方向サイズ [m]。`> 0` |
| `nx` | int | `1` | x 方向分割数。`>= 1` |
| `ny` | int | `1` | y 方向分割数。`>= 1` |

要素数は `2 * nx * ny` です。

##### `kind = "plate_hole"` / `"plane_hole"`

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

##### `kind = "disk"`

XY 平面上の円板です。
内部は極座標で分割され、中心から外周へ向かって三角形化されます。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | 円板中心 `[x, y, z]` [m] |
| `radius` | float | `0.5` | 円板半径 [m]。`> 0` |
| `n_theta` | int | `24` | 周方向分割数。`>= 3` |
| `n_r` | int | `4` | 半径方向分割数。`>= 1` |

内部的には `inner_radius=0` の `annulus` と同じ生成経路を使います。

##### `kind = "annulus"`

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

##### `kind = "box"`

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

##### `kind = "cylinder"`

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

##### `kind = "sphere"`

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

`conductor` の制約:

- `field_bc_mode="free"` の直接 Coulomb 係数で電荷を再配分します。
- `field_bc_mode="periodic2"` とは併用できません。
- 導体要素数が大きい場合、dense solve のバッチごとのコストが増えます。

---

### `[field]`: 要素核

`element_kernel="point"` が互換既定値です。`sim.softening` はこの point kernel に適用します。

`element_kernel="triangle_p0"` の規則:

| 項目 | 規則 |
| --- | --- |
| source | 各 `q_elem` を三角形上の一様な面電荷密度として扱う |
| solver | `direct` / `treecode` / `fmm` / `auto`。`auto` は `tree_min_nelem` で direct / FMM を選ぶ |
| 必須条件 | `sim.softening=0`、全表面が `insulator` |
| Treecode | 厳密 panel near + monopole far |
| FMM | 厳密 panel near + 厳密 triangle P2M |
| 面の向き | OBJ は `[mesh].surface_side`、template は各 `[[mesh.templates]].surface_side` を指定 |

`outward_closed` は、法線の向きが整合した閉じた two-manifold だけに使用できます。

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
| `outer_plasma_profile.csv` | outer stateが有効な`kinetic_1d` / `unified_linear_response`のprofile。条件付きcheckpoint |
| `photoelectron_histogram.csv` | `photoelectron_histogram_enabled=true`の前batch・累積histogram。条件付きcheckpoint |
| `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | `outer_queue_enabled=true`のactive event。serialでは前者、MPIではrankごとに後者を使う条件付きcheckpoint |
| `mesh_potential.csv` | `write_mesh_potential=true` のとき |
| `charge_history.csv` | `history_stride > 0` のとき |
| `potential_history.csv` | `write_potential_history=true` かつ `history_stride > 0` のとき |
| `performance_profile.csv` | `BEACH_PROFILE=1`のとき |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | serialまたはMPI rank別の乱数状態 |
| `macro_residuals.csv` | MPIでも単一のglobalマクロ粒子数残差 |
| `charge_ledger.csv` | 粒子種別の電荷収支、粒子数、再開用累積値 |

histogram state が ready な場合、`summary.txt` へ次を追加します。

| 種類 | キー |
| --- | --- |
| histogram 定義 | `photoelectron_histogram_bins`, `photoelectron_histogram_energy_max_J` |
| 進行 | `photoelectron_last_completed_batch` |
| 累積値 | `photoelectron_cumulative_signed_charge_C`, `photoelectron_cumulative_kinetic_energy_J`, `photoelectron_cumulative_count` |
| 前 batch | `photoelectron_previous_signed_current_A`, `photoelectron_previous_charge_ratio` |
| 適用性 | `photoelectron_max_charge_ratio`, `photoelectron_linear_applicability_status` |

`coupling_steady_start_mode`、`coupling_steady_start_mesh_id`、`coupling_outer_queue_enabled`は常に
`summary.txt`へ出力します。`zhao_floating`の新規開始時は、resolved branch、$E_I$、$Q_{seed}$、mesh IDを
`zhao_steady_start_branch=...` で始まる1行として標準出力に記録します。resume時は再seedせず、
`zhao_steady_start_restored_after_batches=...`に復元したbatch数を記録します。`coupling_outer_queue_enabled=T`のときだけ、
次のqueue stateを追加します。

| 種類 | キー |
| --- | --- |
| closure | `outer_photoelectron_population_fraction` |
| column | `outer_photoelectron_column_per_area_m2`, `outer_photoelectron_column_target_per_area_m2`, `outer_photoelectron_column_residual_per_area_m2` |
| queue stock | `outer_queue_event_count`, `outer_queue_signed_charge_C`, `outer_queue_fingerprint` |

各値の場所は[構成固有の出力](OutputGuide.html#構成固有の値を探す)にまとめています。

[再開に使う出力ファイル](OutputGuide.html#再開に使うファイル)に、列定義と条件別 checkpoint 要件を集約しています。

`mesh_potential.csv` の評価規約:

- 要素重心の電位 [V] を記録します。
- 自己項は `softening > 0` なら `1/softening`、それ以外は面積等価半径近似です。
- `periodic2` は explicit image shell を加えます。
- `cached_kneq0` は cached 非零モードと境界条件付き `k=0` を加えます。

`potential_history.csv`:

- `charge_history.csv` と同じ `history_stride` で記録します。
- 形式は `batch, elem_idx, potential_V` です。
- 履歴ごとに `field_solver%refresh` と `compute_mesh_potential` を実行するため、計算コストが増えます。

`resume=true` の要件:

| 条件 | 内容 |
|---|---|
| 出力 | `write_files=true` が必須 |
| 読み込み元 | `restart_from` 未指定なら `output.dir`、指定時は `restart_from` |
| 必須ファイル | `summary.txt`, `charges.csv`, serialの`rng_state.txt`またはMPI全rankの`rng_state_rankNNNNN.txt` |
| 条件付きファイル | ledger metadataがある場合の`charge_ledger.csv`、readyなouter stateの`outer_plasma_profile.csv`、histogram有効時の`photoelectron_histogram.csv`、queue有効時の`outer_event_queue.csv`またはMPI全rankの`outer_event_queue_rankNNNNN.csv` |
| 任意state | `macro_residuals.csv`が存在すればglobal残差を復元 |
| 挙動 | 必須 checkpoint がなければ新規実行にフォールバックせず停止 |

`restart_from` は checkpoint の読み込み元だけを変更します。新しい出力は常に `output.dir` に書きます。

MPI 実行時:

| ファイル | 内容 |
|---|---|
| `rng_state_rankNNNNN.txt` | rank 別乱数状態 |
| `outer_event_queue_rankNNNNN.csv` | Zhao過渡queueのrank-local active event。全rank分を保存 |
| `macro_residuals.csv` | 全rankで共有するglobal残差。rootが1個だけ書く |

再開時の整合条件:

- 旧形式の `macro_residuals_rankNNNNN.csv` がある checkpoint は、暗黙変換せず拒否します。
- `summary.txt` の `mpi_world_size` は現在の rank 数と一致させます。
- schema v2/v3/v4はmodel、ordered mesh、ordered speciesのfingerprint一致が必要です。
- schema v3 の outer profile は `field_V_m` と `charge_density_C_m3` が必須です。
- schema v4のqueue再開はtransient Zhao state、queue file schema 2、rank、world size、完了batch、global count、signed
  charge、全rank queue fingerprintが一致しなければ停止します。
- `[[particles.species]].species_key` は安定 ID です。省略時は `species_<1-based index>`、明示時は粒子種間で一意にします。

---

## 座標・配置の補助パラメータ

次のkeyも通常のTOML parameterです。ただし、読み込み時に右列の実座標・実寸を計算します。
同じ値を直接指定したときにエラーになるものと、計算値で上書きするものがあるため区別してください。

| key・指定値 | 型・条件 | 計算する値 | 直接指定したkeyとの関係 |
|---|---|---|---|
| `sim.box_origin`, `sim.box_size` | float[3]。2つを同時指定し、`box_size > 0` | `box_min = box_origin`、`box_max = box_origin + box_size` | `box_min`または`box_max`との併用はエラー |
| `inject_region_mode="face_fraction"`, `uv_low`, `uv_high` | `uv_*`はfloat[2]かつ`[0,1]`内。`reservoir_face` / `photo_raycast`のみ | `inject_face`上の`pos_low`, `pos_high` | `pos_low`または`pos_high`との併用はエラー |
| templateの`placement_mode="box_anchor"`, `anchor`, `offset`または`offset_frac` | `anchor`はbox centerまたは各face center。`offset`は[m]、`offset_frac`はbox幅に対する割合 | templateの`center` | `center`との併用はエラー。`offset`と`offset_frac`の併用もエラー |
| templateの`size_mode="box_fraction"`, `size_frac` | 対応kindは`plane` / `plane_hole` / `plate_hole` / `box` / `sphere` / `cylinder` | 下表の寸法key | 対応する寸法keyを明示していてもエラーにせず、計算値で上書き |
| `[mesh.groups.<name>]`の`placement_mode`, `anchor`, `offset`または`offset_frac` | `absolute`または`box_anchor` | grouped templateが共有する`group_origin` | `offset`と`offset_frac`の併用はエラー。`absolute`で`anchor`を指定するのもエラー |
| templateの`group`, `center_local` | `[mesh.groups.<name>]`の定義が必要 | `center = group_origin + group_scale * center_local` | `center`, `placement_mode`, `anchor`, `offset`, `offset_frac`, `size_mode`, `size_frac`との併用はエラー |
| `[mesh.groups.<name>]`の`scale`または`scale_from` + `scale_factor` | `scale > 0`。`scale_from`はbox幅の参照名 | grouped templateで明示した長さkeyをscale倍 | 下表の明示寸法を`scale * input`で上書き。未指定でdefaultになった寸法はscaleしない |

`size_mode="box_fraction"`が上書きする寸法は形状ごとに異なります。

| `kind` | `size_frac` | 上書きするkey |
|---|---|---|
| `plane`, `plane_hole`, `plate_hole` | float[2] | `size_x`, `size_y` |
| `box` | float[3] | `size` |
| `sphere` | float | `radius`。box 3辺の最小値を基準にする |
| `cylinder` | float[2] | `radius`, `height`。radiusはx/y幅の最小値、heightはz幅を基準にする |

補助パラメータの選択肢:

- group scale は、template で明示した `size_x`, `size_y`, `size`, `radius`, `inner_radius`, `height` にだけ掛けます。
- `anchor` は `box_center` または各軸の `*_low_face_center` / `*_high_face_center` です。
- `scale_from` は `box_x`, `box_y`, `box_z`, `box_min_xy`, `box_max_xy`, `box_min_xyz`, `box_max_xyz` から選びます。
- `placement_mode="absolute"`、`size_mode="absolute"`、`inject_region_mode="absolute"` は直接指定値をそのまま使います。

入力前に[`beach.toml`を作成・検証する](Configuration.html)の `beachx lint` で組合せを確認してください。

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
