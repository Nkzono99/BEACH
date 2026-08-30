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
| 開発実行 | `fpm run --target beach -- path/to/beach.toml` でも同じ |
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
| 長さ | `domain.box_min`, `domain.box_max`, `pos_low`, `pos_high` | m |
| 電荷 | `q_particle`, 要素電荷出力 | C |
| 質量 | `m_particle` | kg |
| 速度 | `drift_velocity`, `ray_direction` | m/s。ただし `ray_direction` は方向ベクトル |
| 電場 | `e0`, `e0_abs` | V/m |
| 磁場 | `b0` | T |
| 密度 | `number_density_cm3`, `number_density_m3` | cm^-3 または m^-3 |
| 温度 | `temperature_k`, `temperature_ev` | K または eV。両方の同時指定は不可 |
| 角度 | `e0_phi_xy_deg`, `e0_phi_z_deg` | degree |

`*_low` / `*_high` は各軸の下限・上限です。
`inject_face` は `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, `z_high` のいずれかを指定します。

---

## 公式入門ケース

最初の実行には [10 分チュートリアル](Tutorial.html) と
`examples/tutorial_insulator.toml` を使ってください。`beachx config init`
も同一の設定を生成します。この入門ケースは`domain.periodic_axes=["x","y"]`と
`field_boundary.mode="periodic2"`を使います。

境界reservoir流入、`plane_source`、閉じた光電子、無限周期補正は、そのケースが実行できた後に
[シミュレーションケースを設計する](ConfigurationRecipes.html)から追加する応用設定です。

---

## TOML の階層とセクション一覧

`[sim]`、`[domain]`、`[field_boundary]`、`[particle_boundary]`、`[reservoir]`、
`[surface_current_model]`、`[particles]`、`[mesh]`、`[periodic2]`、`[output]`が公開構成です。

```text
beach.toml
├── [sim]
├── [domain]
├── [field_boundary]
├── [particle_boundary]
├── [reservoir]
├── [surface_current_model]
├── [particles]
│   └── [[particles.species]]       # 1 件以上の array-of-tables
│       ├── [particles.species.boundary_inflow]
│       └── [particles.species.boundary]
├── [mesh]
│   ├── [mesh.groups.<name>]        # 名前付きの子 table
│   └── [[mesh.templates]]          # 0 件以上の array-of-tables
├── [periodic2]
└── [output]
```

本文中の`sim.dt`や`domain.periodic_axes`は「table名.key」の参照表記です。

| TOML table | 親 | 件数・必須条件 | 内容 |
|---|---|---|---|
| `[sim]` | root | 条件付き | 時間刻み、バッチ数、場ソルバ、外部場 |
| `[domain]` | root | box使用時 | box geometryと周期軸 |
| `[field_boundary]` | root | 任意 | 場の`free` / `periodic2` closure |
| `[particle_boundary]` | root | 任意 | 非周期面のglobal粒子作用 |
| `[reservoir]` | root | 任意 | 外部reservoirの流入障壁と基準電位 |
| `[surface_current_model]` | root | 任意 | species別`fixed_current` targetまたはmatching-plane応答を解く外部シースclosure |
| `[particles]` | root | 必須 | `[[particles.species]]` のコンテナ。直下に通常 key は置かない |
| `[[particles.species]]` | `[particles]` | 1 件以上 | 粒子種、注入方式、速度分布、マクロ粒子重み |
| `[particles.species.boundary]` | 最新の`[[particles.species]]` | 任意 | その粒子種だけの非周期面override |
| `[particles.species.boundary_inflow]` | 最新の`[[particles.species]]` | 任意 | 外部reservoirから流入する非周期面 |
| `[mesh]` | root | 任意 | OBJ または組み込み template の選択 |
| `[mesh.groups.<name>]` | `[mesh]` | 0 件以上 | 複数 template で共有する配置と scale |
| `[[mesh.templates]]` | `[mesh]` | 0 件以上 | `mode="template"` で使う組み込み形状 |
| `[periodic2]` | root | 条件付き | split periodic2 の非零モード・零モード・下側境界モデル |
| `[output]` | root | 任意 | 出力先、履歴、checkpoint 再開 |

`boundary_inflow`、`plane_source`、`reservoir_face`、`photo_raycast`を使う場合、`[sim]`と`[domain]`が必要です。
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
| `batch_duration` | float | `0.0` | 1 バッチの物理時間 [s]。適応的な非零モード進行では、各 accepted batch の最大試行幅 |
| `batch_duration_step` | float | `0.0` | `batch_duration = dt * batch_duration_step` として解決 |
| `max_step` | int | `400` | 粒子 1 個あたりの最大 push 回数 |
| `tol_rel` | float | `1.0e-8` | 相対変化量の監視値。停止条件には未使用 |
| `q_floor` | float | `1.0e-30` | `rel_change` 計算時の分母下限 |
| `multiple_box_events_policy` | string | `"abort"` | 1 step の境界event上限超過時に `abort` / `soft_discard` |
| `multiple_box_events_retry_backend` | string | `"none"` | `multiple_box_events` 後の再試行。`none` / `upper_panel_fourier` |
| `multiple_box_events_soft_discard_count_grace` | int | `1000` | 累積 soft discard 率の判定を開始する件数猶予。`>= 0` |
| `multiple_box_events_soft_discard_fraction_limit` | float | `1.0e-6` | 累積 soft discard 率の停止上限。`0 < value <= 1` |
| `multiple_box_events_soft_discard_abs_charge_limit` | float | `1.0e-12` | 累積 soft discard 絶対電荷の停止上限 [C] |
| `raycast_max_bounce` | int | `16` | `photo_raycast` の最大bounce数 |

`batch_duration` と `batch_duration_step` の同時指定はエラーです。
`boundary_inflow` / `plane_source` / `reservoir_face` / `photo_raycast` では、解決後の `batch_duration > 0` が必須です。

`upper_panel_fourier` は `cached_kneq0` の `periodic2` 構成だけで使えます。通常の場計算は変更せず、
境界event上限を超えた 1 step だけを元の状態から再計算します。詳細と成立域は
[粒子イベント](ParticleEvents.md)にまとめています。

`soft_discard` では、累積件数を $D$、accepted batch で処理した累積 macro particle 数を $P$、
`multiple_box_events_soft_discard_count_grace` を $G$、
`multiple_box_events_soft_discard_fraction_limit` を $f_{\mathrm{limit}}$ とすると、
$D>G$ かつ $D/P>f_{\mathrm{limit}}$ で停止します。いずれの閾値も等値なら許容します。

累積絶対 macro charge の上限はこの率判定とは独立で、超過すれば停止します。

`summary.txt` と checkpoint では
`multiple_box_events_soft_discarded`、`multiple_box_events_soft_discard_fraction`、
`multiple_box_events_soft_discarded_abs_charge_C` を監査してください。絶対電荷上限は物理的な誤差 budget です。
累積率は後半の burst を希釈しうるため、batch ごとの集約 log も併せて確認します。

#### 場ソルバ

`field_solver` は、境界要素電荷から評価点の Coulomb 電場を計算する方式です。
選択肢ごとの対応パラメータは下表のとおりです。

| `field_solver` | 用途 | 対応する場境界 |
|---|---|---|
| `direct` | 要素数が小さい場合の厳密な全対全評価、split reference | `free`、または条件を満たす`periodic2` split reference |
| `treecode` | 中規模以上の近似評価 | `field_boundary.mode="free"` |
| `fmm` | 大規模評価、`periodic2`、FMM コア検証 | `field_boundary.mode="free"` / `"periodic2"` |
| `auto` | 要素数に応じて direct / FMM を自動選択 | `field_boundary.mode="free"` |

solverと場境界の組合せは[場の評価の互換表](FieldSolvers.html#solverと場境界の互換表)で確認できます。
要素電荷は常に三角形上の一定面密度（P0 panel）として評価され、設定でkernelを選択しません。

共通キー:

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `direct` / `treecode` / `fmm` / `auto` |
| `field_normalization` | string | `"si"` | `si` / `box` / `mesh` / `length` |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` で使う長さスケール [m] |
| `field_periodic_image_layers` | int | `1` | `periodic2` の近傍 image shell |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | cache生成時のEwald分割値。`0`は自動 |
| `field_periodic_ewald_layers` | int | `4` | cache生成時の実空間・逆空間shell深さ |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | operator cacheの保存先 |
| `field_periodic_generation_tolerance` | float | `1.0e-8` | cache identityに含める生成許容誤差 |

`field_periodic_far_correction="auto"` は互換用に受理され、現在は `none` と同じ扱いです。

`field_normalization` は場計算内部の座標と周期 cell の正規化だけを変えます。
出力される電場・電位は SI に戻されます。

| `field_normalization` | 長さスケール |
|---|---|
| `si` | 入力 SI 座標をそのまま使う |
| `box` | `domain.box_max - domain.box_min`の最大幅。`[domain]`が必須 |
| `mesh` | mesh bounding box の最大幅。mesh が空なら `field_length_scale` |
| `length` | `field_length_scale` |

##### `field_solver = "direct"`

全 source 要素を直接足し上げます。
近似誤差はなく、計算量は評価点数を `M`、要素数を `N` として `O(MN)` です。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"direct"` を指定 |
| `field_normalization` | string | `"si"` | direct 評価前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |

`tree_theta`、`tree_leaf_max`、`tree_min_nelem` は `direct` では使いません。

##### `field_solver = "treecode"`

source octree を作り、遠方 node は monopole 近似、近傍 node は解析的 panel kernel で評価します。
FMM のような local expansion は使わず、評価点ごとに木を走査します。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"treecode"` を指定 |
| `field_normalization` | string | `"si"` | tree 構築前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_theta` | float | `0.5` | MAC パラメータ。`0 < theta <= 1`。大きいほど速く粗い |
| `tree_leaf_max` | int | `16` | 葉 node あたり最大 source 数。`>= 1` |

`tree_min_nelem` は `field_solver="auto"` 用のしきい値なので、明示 `treecode` では切替には使いません。

##### `field_solver = "fmm"`

simulator 非依存の Coulomb FMM コアを使います。
source 幾何の plan と、電荷更新ごとの state を分け、P2M/M2M/M2L/L2L/L2P と近傍 direct 和で評価します。
選択と精度確認は[FMM](FMM.html)、内部実装は
[Coulomb FMMコア詳細](FMMCore.html)にまとめています。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"fmm"` を指定 |
| `field_normalization` | string | `"si"` | FMM plan 構築前に座標を正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_theta` | float | `0.5` | near/far 判定の MAC パラメータ。`0 < theta <= 1` |
| `tree_leaf_max` | int | `16` | source tree の葉 node あたり最大 source 数。`>= 1` |
`tree_min_nelem` は明示 `fmm` では使いません。

##### `field_solver = "auto"`

要素数が `tree_min_nelem` 未満なら direct、以上なら FMM を使います。

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `"auto"` を指定 |
| `field_normalization` | string | `"si"` | 自動選択前に共通で使う正規化 |
| `field_length_scale` | float | `1.0` | `field_normalization="length"` または mesh fallback で使用 |
| `tree_min_nelem` | int | `256` | FMM へ切り替える要素数しきい値。`>= 1` |
| `tree_theta` | float | `0.5` | FMM 選択時の near/far MAC パラメータ |
| `tree_leaf_max` | int | `16` | FMM 選択時の葉 node あたり最大 source 数 |

`tree_theta` と `tree_leaf_max` は、明示指定がなければ要素数から次の値を使います。

| 要素数 `nelem` | `tree_theta` | `tree_leaf_max` |
|---:|---:|---:|
| `< 1500` | `0.40` | `12` |
| `1500 <= nelem < 10000` | `0.50` | `16` |
| `10000 <= nelem < 50000` | `0.58` | `20` |
| `50000 <= nelem` | `0.65` | `24` |

### `[domain]`: box geometryと周期topology

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]
periodic_axes = ["x", "y"]
```

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `box_min`, `box_max` | float[3] | なし | box下限・上限 [m]。両方を指定 |
| `box_origin`, `box_size` | float[3] | なし | originと正のsize [m]。両方を指定 |
| `periodic_axes` | string[] | `[]` | 重複しない`"x"` / `"y"` / `"z"` |

2組のgeometry表現は同時指定できません。`[domain]`を使う場合はどちらか1組が必須です。
周期性を指定できる公開キーは`domain.periodic_axes`だけです。
`field_boundary.mode="periodic2"`では、`periodic_axes=["x","y"]`が必要です。

### `[field_boundary]`: 場のclosure

```toml
[field_boundary]
mode = "periodic2"
```

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `mode` | string | `"free"` | `free` / `periodic2` |

`periodic2`では`[domain]`が必要です。非零・零モードoperatorを明示する場合は`[periodic2]`も指定します。
solverとの組合せは
[場の評価の互換表](FieldSolvers.html#solverと場境界の互換表)に従います。

### `[particle_boundary]`: global粒子境界

```toml
[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"
```

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `x_low`, `x_high` | string | domainに従う | 非周期x面の`open` / `reflect` / `redistributed_reflect` |
| `y_low`, `y_high` | string | domainに従う | 非周期y面の`open` / `reflect` / `redistributed_reflect` |
| `z_low`, `z_high` | string | domainに従う | 非周期z面の`open` / `reflect` / `redistributed_reflect` |
| `ordinary_open_model` | string | `"escape"` | effectiveなopen面で`escape` / `potential_barrier` |

省略した面はdomain topologyを継承し、非周期面ではopenになります。
周期面をこのtableで上書きできません。
`reflect`は法線速度だけを反転し、接線速度とevent位置の接線成分を維持します。

`redistributed_reflect`は同じ速度作用に加え、単一面では面内2軸をbox spanの両端guardを除く範囲から一様再標本化します。
edge / cornerの同時eventではevent maskに含まれない軸だけを再標本化します。

`potential_barrier`は境界通過点の電位`phi_boundary`と`reservoir.phi_infty`から
`q_particle * (phi_infty - phi_boundary)`を評価し、法線運動エネルギーを超える正の障壁だけを反射します。

### `[reservoir]`: 外部reservoir条件

```toml
[reservoir]
inflow_model = "source_vdf"
phi_infty = 0.0
face_potential_grid_n = 3
```

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `inflow_model` | string | `"source_vdf"` | `source_vdf` / `infinity_barrier` |
| `phi_infty` | float | `0.0` | barrierが参照する無限遠基準電位 [V] |
| `face_potential_grid_n` | int | `3` | 注入面平均電位の`N x N`評価格子。`>=1` |

`infinity_barrier`は`boundary_inflow`の面平均電位と`phi_infty`からreservoirの法線VDFを補正します。
内部の`plane_source`とdeprecatedな`reservoir_face`には適用しません。
一様電場には有限な無限遠電位がないため、併用時は`phi_infty`を有効なreservoir基準として整合させます。

### `[surface_current_model]`: 外部シースclosure

このトップレベルtableを省略すると、speciesごとの`target_*_current_a`を使う手動設定になります。
`model="zhao_stationary"`は、平面・無衝突・非磁化の外部シースについてZhaoのA/B/C零電流定常根を解き、
ambient electron、ion、PE emission、PE escape、PE returnの電流とz-high kinetic barrierを一度だけ決定します。

`model="matching_plane_quasistatic"`は、box上端をmatching planeとして、選択したresponse backendの外部シースと
準定常に連結します。`response_backend="table"`が既定で、`"zhao_online"`は有限$H$のcharge-driven Zhao
A/B/C responseをBEACH内で解きます。

`zhao_stationary`で`photoelectron_source_scale=0.0`にすると、PE channelを作らず、ambient electronとionだけの
零電流根を使えます。

```toml
[surface_current_model]
model = "zhao_stationary"
zhao_branch = "auto"
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_species = "photoelectron"
solar_elevation_deg = 60.0
photoelectron_ref_density_m3 = 6.4e7
photoelectron_source_scale = 1.0
# reference_area_m2 = 1.0e-8
```

光電子なしではPE固有キーをすべて省略し、source scaleを明示的にゼロにします。

```toml
[surface_current_model]
model = "zhao_stationary"
zhao_branch = "auto"
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_source_scale = 0.0
```

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `model` | string | `"none"` | `none` / `zhao_stationary` / `matching_plane_quasistatic` |
| `response_backend` | string | `"table"` | matchingの応答源。`table` / `zhao_online` |
| `zhao_branch` | string | `"auto"` | `auto` / `a` / `b` / `c`。stationaryまたはonline Zhaoのbranch。PEなしstationaryは`auto` / `c`のみ |
| `electron_species` | string | 必須 | ambient electronの`species_key` |
| `ion_species` | string | 必須 | cold ionの`species_key` |
| `photoelectron_species` | string | PE有効時に必須 | PE emission/returnを追跡する`photo_raycast`の`species_key`。matchingで省略するとPEなし |
| `solar_elevation_deg` | float | stationaryのPE有効時に必須 | Zhao sourceに使う太陽高度角 $\alpha$。$0<\alpha\le90$ degree |
| `photoelectron_ref_density_m3` | float | stationaryのPE有効時に必須 | PE基準密度 $n_{pe,ref}$ [m^-3] |
| `photoelectron_source_scale` | float | `1.0` | stationary Zhaoの$n_{pe,0}=s_{UV}n_{pe,ref}\sin\alpha$に使う$s_{UV}$。`0.0`はPEなし |
| `reference_area_m2` | float | domainのx-y面積 | Zhao電流密度を総電流へ変換する面積 [m^2]。matchingでは指定不可 |
| `response_table_path` | string | table matchingで必須 | 外部シース応答CSV v1。相対パスは`beach.toml`のディレクトリ基準、絶対パスはそのまま。onlineでは指定不可 |
| `implicit_zero_mode` | bool | `false` | matching tableの面平均$D_H$だけを後退Eulerで更新。`e_bottom_zero`、2 node以上の$D_H$軸、singleton feedback軸が必須 |
| `coupling_rtol` | float | `1.0e-4` | matching固定点反復の相対収束許容値。有限な$0<r\le1$ |
| `coupling_atol` | float[4] | `[0.0, 0.0, 0.0, 0.0]` | feedback成分ごとの絶対許容値。順にPE外向きflux [m^-2 s^-1]、PE平均法線energy [eV]、electron外向きflux [m^-2 s^-1]、ion外向きflux [m^-2 s^-1]。各値は有限かつ非負、inactive成分は0 |
| `coupling_max_iterations` | int | `20` | matching固定点反復の最大回数。`>=1` |
| `coupling_relaxation` | float | `0.5` | matching更新の緩和係数。有限な$0<\omega\le1$ |

#### Zhao stationary closure

参照するspeciesはenabledかつ相異なり、`surface_charge_closure="fixed_current"`を指定します。
手動の`target_absorbed_current_a` / `target_emission_current_a`は同時指定できません。ambient electronとionは
z-high reservoirから内向きに流入します。

`photoelectron_source_scale=0.0`では`photoelectron_species`、`solar_elevation_deg`、
`photoelectron_ref_density_m3`を指定せず、`zhao_branch`は`"auto"`または`"c"`にします。
Zhao Type Cの$J_e+J_i=0$を解き、electron/ionの吸収targetとz-high kinetic mapだけを生成します。
PE emission、return、escape targetは生成しません。

PEを有効にする場合、PEは負電荷`photo_raycast`、`inject_face="z_high"`、
`deposit_opposite_charge_on_emit=true`、有効なz-high particle boundaryの`open`を要求します。3 speciesは単価電荷で、ambient electronとPEの質量は一致、
$T_e>0$、$T_{pe}>0$、$T_i\le0.1T_e$でなければなりません。

外部closureは非磁化なので`sim.b0=[0,0,0]`も必須です。Zhao固有の0 V reservoirを使うため、
`reservoir.inflow_model="infinity_barrier"`とは併用できません。

`ion_species`の`number_density_*`を無限遠ion密度として使い、ambient electron密度は定常根が解きます。
electron側の設定密度とPEの`emit_current_density_a_m2`はraw Monte Carlo mapを生成するためのsampling入力であり、
fixed-current targetにはなりません。解いたelectron密度は`summary.txt`へ出力します。

解いたsigned電流密度を$J_e<0$、$J_i>0$、$J_{emit}>0$、$J_{escape}>0$とすると、PE returnは
$J_{return}=J_{escape}-J_{emit}\le0$です。BEACHは面積$A$を掛け、electron/ion/returnを各吸収channelへ、
emissionをPE放出反作用channelへ渡し、外向きPE粒子電流$-AJ_{escape}$を外部escape targetとして記録します。
escape targetは表面要素へdepositされません。rawな境界escape電荷との差は外部境界closureの補正としてledgerへ
記録されます。これによりPE連続条件と表面定常電流を別々に検証できます。

BEACH内のraycast・軌道追跡から得た要素別分布とraw escape統計は上書きされません。吸収・放出には
raw / target / applied、escapeにはraw / target / applied / correctionを出力します。

Zhaoを有効にすると、選択したambient speciesのz-high流入VDFは無限遠0 Vから現在の面平均電位まで
エネルギー保存で写像されます。Type A electronは$\phi_m$もaccess bottleneckとして通過するtailだけを生成し、
Type B/C electronとionの外部bottleneckは0 Vです。

PE放出VDFは`temperature_ev`と`normal_drift_speed`で指定した表面half-Maxwellianのままです。外向きPEは
Type Aでは$\phi_m$、Type B/Cでは0 Vまでの残りのbarrierをz-highで
判定し、不足分をreturn、高エネルギーtailをescapeへ分けます。fixed-current補正はこのraw分布を保って総量を
Zhao targetへ合わせます。このmodel固有写像は参照speciesのz-highにだけ適用され、一般の`reservoir.inflow_model`
とは独立です。

このmodelはbox外の電場、空間電荷、Debye shielding、return軌道・遅延を解かず、batch中の表面電位から
電流を更新しません。z-highでのreturnは外部turning pointまでの距離と時間を省略します。外部シースの
自己無撞着な過渡解ではなく、固定された定常外部電流closureです。
計算例は`examples/periodic2_zhao_fixed_current.toml`です。

#### Matching-plane quasistatic closure

table backendは次のように指定します。`response_backend`を省略した場合も`"table"`です。

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "table"
response_table_path = "responses/outer_sheath.csv"
implicit_zero_mode = false
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.0, 0.0, 0.0] # 既定値: 相対許容値だけを使う
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

online Zhao backendではresponse pathを指定しません。

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.05, 0.0, 0.0] # 有限ray samplingに対するPE energy許容値 [eV]
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

matchingでPEを使わない場合は`photoelectron_species`とPE speciesを省略します。online ZhaoはPE populationを
0として外部reservoirへ接続します。table backendではPE flux/energy入力軸を両方0にしてください。

active feedback成分$j$は、backend scaleを$s_j$、`coupling_rtol`を$r$、`coupling_atol[j]`を$a_j$として
$|X_{raw,j}-X_j|\le\max(r s_j,a_j)$で判定します。上の`0.05` eVは有限ray sampling向けの設定例であり、
macro-particle数を変えた収束試験から決めます。online Zhaoではambient electron / ion外向きfluxに対応する
第3・第4成分がinactiveです。tableのsingleton feedback軸もinactiveで、これらに非零の絶対許容値を指定すると
設定を拒否します。絶対許容値が支配する成分は残差を換算するため、accepted trialの
`matching_plane_residual`は引き続き`coupling_rtol`以下です。

matching planeは`domain.box_max`のz成分$H$、連結面積はdomainのx-y面積です。全mesh頂点は$H$より厳密に
下へ置きます。
`reference_area_m2`で別の面積へ置き換えることはできません。

tableの`response_table_path`は外部Zhao計算または1D PICから作った非線形応答CSV v1を指します。
テスト用のsynthetic tableは配線・補間の検証専用であり、
productionの物理入力として有効ではありません。相対pathは設定fileのdirectoryから解決し、解決後のpathを
256文字以内にします。`beachx lint`もこの解決後pathを検査します。

onlineは各queryで$E_H=D_H/\epsilon_0$を境界条件とし、$H$から上流0 V・零電場へ接続する有限$H$の
Sagdeev A/B/C rootを解きます。これは`zhao_stationary`の壁面零電流rootではなく、帯電途中にも$J=0$を
課しません。

`zhao_branch="auto"`は適用可能なbranchを探索し、`"a"` / `"b"` / `"c"`はbranchを固定します。
`auto`で複数の物理解を検出した場合、またはcompatibleなbranchの数値失敗で一意性を確認できない場合は停止します。
v1の複数根検出は有限個のmultistart結果をcluster化するもので、数学的なroot isolationではありません。

$H$は外部半無限領域のinterface原点、zero-mode gauge、PE moment測定面を固定します。平面・並進対称なので
$H$の絶対座標はSagdeev方程式の数値parameterには入らず、壁面から$H$までの距離拘束は解きません。

外向きPEのnumber fluxと平均法線energyは、その2 momentを再現するhalf-Maxwellianへ縮約します。PE fluxが0なら
PE populationは0のまま、PE roleがあればその設定温度、なければambient electron温度を数値scaleに使います。

ambient electron / ionの外向き軸はtransparentかつ固定点残差でinactiveです。online solverはstatelessで、
outer inventory、前rootのcontinuation、外部flight time、遅延return queueを持ちません。

設定したbranch policyで物理解がない場合やsolveが収束しない場合はfail closedとし、明示したbranchやbackendを
暗黙に切り替えません。

このmodelは次の設定契約をすべて要求します。

- `field_boundary.mode="periodic2"`、`domain.periodic_axes=["x","y"]`、z非周期open box
- 明示的な`[periodic2]` split設定。`nonzero_mode_backend`は`cached_kneq0`または
  `panel_spectral_reference`、`zero_mode_policy="exclude_k0"`、下側closureは
  `e_bottom_zero`または`symmetric_vacuum`
- `sim.e0=[0,0,0]`と`sim.b0=[0,0,0]`、genericなreservoir potential modelを使わず、
  `particle_boundary.ordinary_open_model="escape"`。`sim.multiple_box_events_policy`は`"abort"`または
  件数猶予・累積率・絶対電荷上限を指定した`"soft_discard"`
- electron、ion、および任意のphotoelectron roleはenabledかつ相異なり、それぞれ
  `surface_charge_closure="explicit"`。これら以外のenabled speciesは置かない
- electronは負電荷、ionは正電荷。electron/ionは`source_mode="volume_seed"`、
  `npcls_per_step=0`とし、`boundary_inflow`はz-highの`reservoir`だけを指定
- 指定する場合、photoelectronは負電荷で`source_mode="photo_raycast"`、`inject_face="z_high"`、
  `deposit_opposite_charge_on_emit=true`
- 全roleで有効なx/y particle boundaryが`periodic`、z-low/z-highが`open`
- 手動`fixed_current` targetを指定しない

backend固有の制約は次のとおりです。

- tableは`response_table_path`を必須とし、`zhao_branch`を指定しない
- onlineは`response_table_path`を禁止し、`zhao_branch="auto" / "a" / "b" / "c"`を受理する
- onlineは全roleの単価電荷、$T_e>0$、$0\le T_i\le0.1T_e$、正のion number density、
  ambient electron / ionの正の内向きdrift（`drift_velocity`のz成分は負）を要求する。
  PE指定時はambient electronとPEの同一質量および$T_{pe}>0$も要求する
- stationary Zhaoの`solar_elevation_deg`、`photoelectron_ref_density_m3`、
  `photoelectron_source_scale`はmatchingでは指定しない

online Zhaoは平面・無衝突・非磁化の低次元準定常closureです。full VDF、1D PIC、time-dependent outer sheath、
またはambient outward populationの外部returnを解くmodelではありません。

`model="none"`では`model`以外のキーを指定しません。
廃止済みトップレベル`[outer_plasma]` / `[coupling]`
は復活していません。連結設定はすべて`[surface_current_model]`に置きます。table CSVの列契約、online Zhaoの縮約、
固定点反復、検証手順の詳細は[matching-plane準定常連成](MatchingPlaneCoupling.html)に示します。

### `[periodic2]`: 非零モード・零モード・下側境界

`[periodic2]`はトップレベルtableです。`domain.periodic_axes=["x","y"]`と
`field_boundary.mode="periodic2"`を指定します。productionでは`field_solver="fmm"`を使い、
小規模検証用のsplit referenceに限り`field_solver="direct"`を使います。

| キー | 既定値 | 意味 |
|---|---:|---|
| `nonzero_mode_backend` | 必須 | `panel_spectral_reference` / `cached_kneq0` |
| `zero_mode_policy` | 必須 | `exclude_k0` |
| `lower_boundary_model` | 必須 | `symmetric_vacuum` / `e_bottom_zero` |
| `max_nonzero_mode_potential_step` | `0` | accepted trial で許す $k\ne0$ 電位変化の上限 [V]。省略または `0` で無効 |
| `reference_mode_layers` | `4` | Fourier mode cutoff |
| `panel_quadrature_order` | `12` | panel 面積積分次数 |

`max_nonzero_mode_potential_step > 0` は `nonzero_mode_backend="cached_kneq0"` でだけ使えます。
解決後の `sim.batch_duration` を $h_0$ とし、各 accepted batch で
$h_0,h_0/2,h_0/4,\ldots$ の固定 ladder を順に試します。各 trial の候補電荷
$\mathbf q_{\mathrm{candidate}}$ について、全 panel 重心で

$$
\max_j\left|
\left[P_{k\ne0}
  \left(\mathbf q_{\mathrm{candidate}}-\mathbf q_{\mathrm{current}}\right)
\right]_j
\right|
\le
\texttt{max\_nonzero\_mode\_potential\_step}
$$

を最初に満たした幅を受理します。`P_{k\ne0}` は cached 非零モード電位 operator です。
棄却 trial は RNG とmacro粒子数残差を trial 前へ戻し、
統計、履歴、charge ledger へ加えません。

この機能はtime-scaledな`boundary_inflow` / `plane_source` / `reservoir_face` / `photo_raycast`に対応します。
流束駆動のspeciesには`target_macro_particles_per_batch`が必要で、固定`w_particle`は使えません。
`volume_seed`とは併用できません。

`sim.batch_count` は accepted batch数であり、`simulated_time_s` は受理した幅の総和です。
適応区間ではdynamic teamを無効化し、同じ実行内では全MPI rankの実OpenMP team sizeを揃えます。
restart前後のteam size一致は要求せず、再開後の実team sizeを新しい診断値として記録します。

この判定は frozen-field近似中の局所電位変化を制限するtrust boundであり、局所打切り誤差を
保証しません。結果の時間幅収束は、この上限を半分にした計算と比較してください。

### periodic2の組合せ制約

periodic2では`[domain]`、`periodic_axes=["x","y"]`、`field_boundary.mode="periodic2"`が必要です。
`examples/periodic2_closed_photoelectron.toml`は、x/y周期、境界reservoir、閉じた光電子を組み合わせた基準例です。
同じ周期条件をfield、collision、`photo_raycast`に適用します。

| `nonzero_mode_backend` | 意味 |
| --- | --- |
| `panel_spectral_reference` | 小規模なsplit reference |
| `cached_kneq0` | versioned operatorを再利用するproduction非零モード |

`cached_kneq0`では`exclude_k0` providerが物理的$k=0$を1回加えます。
`symmetric_vacuum`は上下を$\pm Q/(2\epsilon_0A)$、`e_bottom_zero`は下側0・上側$Q/(\epsilon_0A)$とします。

### `[sim]`: 外部場

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | 一様外部電場 [V/m] |
| `e0_abs` | float | 未指定 | 一様外部電場の大きさ [V/m] |
| `e0_phi_xy_deg` | float | `0.0` | `e0_abs` 指定時の xy 面内方位角 [deg] |
| `e0_phi_z_deg` | float | `0.0` | `e0_abs` 指定時の xy 面からの仰角 [deg] |
| `b0` | float[3] | `[0,0,0]` | 一様磁場 [T] |

一様外部電場は、`e0 = [Ex, Ey, Ez]` で直接指定するか、`e0_abs` と角度で指定します。
両形式の混在はエラーです。

---

### `[[particles.species]]`: 粒子種

`[[particles.species]]` は 1 件以上必須です。
`source_mode` によって、使うキーと制約が変わります。

#### 共通キー

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `species_key` | string | `"species_<1-based index>"` | 再開fingerprintとsurface-current role参照に使う安定ID。省略時は定義順から生成し、明示値は粒子種間で一意 |
| `enabled` | bool | `true` | 種を有効化 |
| `source_mode` | string | `"volume_seed"` | `volume_seed` / `plane_source` / `photo_raycast` / deprecated `reservoir_face` |
| `q_particle` | float | `-1.602176634e-19` | 粒子電荷 [C] |
| `m_particle` | float | `9.10938356e-31` | 粒子質量 [kg] |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | 位置下限 [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | 位置上限 [m] |
| `drift_velocity` | float[3] | `[0,0,-8e5]` | ドリフト速度 [m/s] |
| `temperature_k` | float | `2.0e4` | 温度 [K] |
| `temperature_ev` | float | 未指定 | 温度 [eV]。`temperature_k` と排他 |
| `velocity_distribution` | string | `"maxwellian"` | `maxwellian` / `grid` |
| `inject_face` | string | 未指定 | `photo_raycast`の照射開口面。deprecatedな`reservoir_face`でも必須 |
| `source_normal` | float[3] | 未指定 | `plane_source`の一方向法線。axis-alignedな非ゼロベクトル |
| `boundary` | table | 未指定 | `[particles.species.boundary]` のspecies別6面override |
| `boundary_inflow` | table | 未指定 | `[particles.species.boundary_inflow]`のspecies別reservoir流入面 |
| `surface_charge_closure` | string | `"explicit"` | 表面 source 電荷 closure。`explicit` / `fixed_current` / `neutral_return` |
| `target_absorbed_current_a` | float | 未指定 | `fixed_current` の signed 吸収電流 [A]。符号は `q_particle` と一致 |
| `target_emission_current_a` | float | 未指定 | `fixed_current` の signed 放出反作用電流 [A]。符号は `q_particle` と逆 |

#### `[particles.species.boundary]`: species別override

このtableは直前の`[[particles.species]]`に属します。

```toml
[particles.species.boundary]
z_high = "reflect"
```

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `x_low`, `x_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `y_low`, `y_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `z_low`, `z_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |

`inherit`は`[particle_boundary]`のglobal作用を使います。周期面はoverrideできません。
`inject_face`は`photo_raycast`と旧`reservoir_face`の生成面、species境界は生成後の軌道作用です。
`boundary_inflow`は外部からの生成であり、外向き作用を上書きしません。

`surface_charge_closure="neutral_return"`は、負電荷`photo_raycast`、
`deposit_opposite_charge_on_emit=true`、effectiveな`inject_face`境界が`reflect`または
`redistributed_reflect`の組合せだけで使用できます。
解決済み帰還先depositをspecies別のglobal係数で補正し、光電子による表面総電荷増分を0にします。

実escapeまたは`soft_discard`が生じる条件とは併用できません。通常のtracked chargeを
そのままcommitする`"explicit"`が既定です。未帰還率が固定上限5%を超える場合は補正せず停止します。

`surface_charge_closure="fixed_current"` は、軌道追跡で得た要素別分布を保ったまま、species ごとの総電流を
外部 closure の値へ合わせます。targetはspeciesで手動指定するか、トップレベルの
`[surface_current_model]`で自動計算します。手動の吸収 channel は

```toml
surface_charge_closure = "fixed_current"
target_absorbed_current_a = -2.0e-6
```

と指定し、raw 吸収電荷 $R_s$ を $I_{s,\mathrm{abs}}^{\mathrm{target}}\Delta t/R_s$ 倍します。
`photo_raycast` では `deposit_opposite_charge_on_emit=true` として
`target_emission_current_a` も指定でき、放出反作用と帰還吸収を独立に補正します。二つの差である net PE 電流を
直接倍率化しません。target が非ゼロなのに対応する raw channel が空なら fail closed します。

`fixed_current`はraw標本の要素別経験分布を倍率化します。raw hitが1件ならtarget全量がその要素へ入ります。
必要なhit数や許容倍率はmeshと評価量の誤差要件からは自動決定できないため、固定の最小countや最大scaleはありません。

`charge_ledger.csv`の`absorbed_count` / `emitted_count`、raw charge、`fixed_*_weight_scale`を確認し、
マクロ粒子数・`rays_per_batch`・batch幅・乱数seedを変えた要素別電荷分布の収束を確認してください。
ledger residualは収支診断であり、空間分布の統計精度を保証しません。

`fixed_current` と `neutral_return` は同じ species では排他です。外部モデル由来の PE return VDF を別 species として
top 面から注入する場合、その return species に `target_absorbed_current_a` を設定し、top 面の full reflection と
`neutral_return` は使わないでください。同じ return current の二重計上を避けるためです。

#### `source_mode = "volume_seed"`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `npcls_per_step` | int | `0` | 1 バッチに生成するマクロ粒子数 |
| `w_particle` | float | `1.0` | マクロ粒子重み |

制約:

| 条件 | 内容 |
|---|---|
| 粒子数 | 境界流入を使わない場合、有効 species 全体で `npcls_per_step` 合計が 1 以上必要 |
| 重み自動解決 | `boundary_inflow`を持たない`volume_seed`では`target_macro_particles_per_batch`は使用不可 |

`boundary_inflow`を持つspeciesでは`npcls_per_step=0`を許容します。Maxwell分布で正の値なら、
体積seedと境界流入を同じspeciesに加えます。速度gridの境界流入は正の値と併用できません。

#### `[particles.species.boundary_inflow]`

このtableは直前の`[[particles.species]]`に属します。

```toml
[particles.species.boundary_inflow]
z_high = "reservoir"
```

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `x_low`, `x_high` | string | 省略 | `reservoir`。省略時は流入なし |
| `y_low`, `y_high` | string | 省略 | `reservoir`。省略時は流入なし |
| `z_low`, `z_high` | string | 省略 | `reservoir`。省略時は流入なし |

`reservoir`は選択したbox面全体から外部VDFを流入させます。周期面には指定できず、
外向き作用は`[particle_boundary]`と`[particles.species.boundary]`で独立に決まります。
流入面の有効作用は`open`でなければなりません。

Maxwell分布では`number_density_*`と温度、速度gridでは`particle_flux_m2_s`または
`current_density_a_m2`を使います。`w_particle`と`target_macro_particles_per_batch`は排他です。
`reservoir.inflow_model="infinity_barrier"`の場合だけ、面平均電位と`phi_infty`から流入を補正します。
複数面では`target_macro_particles_per_batch`を全流入面の合計targetとして重みを解決し、
speciesとfaceの組ごとにマクロ粒子端数を保持します。

初版では`source_mode="volume_seed"`とだけ併用できます。`plane_source`、`photo_raycast`、deprecatedな
`reservoir_face`と同じspeciesには指定できません。将来、sourceの複数化を拡張できるよう責務は分離しています。
詳しい流束と補正は[reservoir流入](ReservoirInjection.html)を参照してください。

#### `source_mode = "plane_source"`

`plane_source`はbox内部のaxis-aligned矩形面から`source_normal`方向へ一方向の流束を生成します。
粒子数と速度分布に使うキーはdeprecatedな`reservoir_face`と同じです。

| 条件 | 内容 |
|---|---|
| 領域 | `[domain]`が必須 |
| 時間 | `sim.batch_duration > 0`が必須 |
| 矩形面 | `pos_low` / `pos_high`はちょうど1軸で一致し、残る2軸は正の長さを持つ |
| 配置 | 法線座標はbox境界より厳密に内側。接線範囲はbox内で、境界との一致を許容 |
| 方向 | `source_normal`はzero-thickness軸に沿う非ゼロベクトル。内部で正規化し、入力は正負単位ベクトルを推奨 |
| 外部補正 | `[reservoir]`の`infinity_barrier`、`phi_infty`、`face_potential_grid_n`は適用しない |

Maxwell分布では片側流束をdensityとtemperatureから計算します。速度gridでは
`particle_flux_m2_s`または`current_density_a_m2`を使います。面積、`batch_duration`、`w_particle`
からマクロ粒子数を求め、端数を次のbatchへ持ち越します。

#### `source_mode = "reservoir_face"`（deprecated）

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
| 領域 | `[domain]`が必須 |
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

CSVの内容はordered species fingerprintに含まれます。同じパスのCSVを変更すると、既存checkpointからの再開は
分布変更として拒否されます。

CSVはrun内の初回利用時に1回だけ読み込んでsnapshot化し、以後のsamplingと
checkpoint fingerprintは同じsnapshotを使います。実行中にディスク上のCSVを置換しても、途中のbatchから
分布が切り替わることはありません。次のprocessは置換後の内容を新しいmodelとして読み込みます。

`velocity_distribution="grid"`にも、他の `reservoir_face` と同じ面・時間・流束の制約を適用します。

粒子数は次のように決まります。

```text
n_macro_expected = gamma_in * A * batch_duration / w_particle
n_injected = floor(residual + n_macro_expected)
```

残差は次バッチへ繰り越されます。
`target_macro_particles_per_batch > 0` のときは、その値に近づくよう `w_particle` を自動計算します。

このmodeは既存caseの挙動を維持する互換入力です。BEACHは`boundary_inflow`または`plane_source`へ
暗黙変換しません。新しい外部plasma条件には`boundary_inflow`、内部矩形面には`plane_source`を使います。

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
| 領域 | `[domain]`が必須 |
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
`field_boundary.mode="periodic2"`では、periodic imageに命中してもprimary cellへwrapしたhit座標から放出します。

生成した光電子は常に`w_hit`を重みに使い、通常粒子として追跡します。表面へ戻れば通常の衝突として吸収します。

species境界が`inherit`なら`[particle_boundary]`の作用を使います。closed PEでは
`[particles.species.boundary]`で`inject_face`を`reflect`または`redistributed_reflect`にします。

---

### `[mesh]`: メッシュ入力

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `mode` | string | `"template"` | `auto` / `obj` / `template` |
| `obj_path` | string | `"examples/simple_plate.obj"` | OBJ ファイルパス |
| `surface_model` | string | `"insulator"` | OBJ 全体の表面モデル |
| `surface_side` | string | `mode="obj"`または`"auto"`で必須 | OBJ panelの真空側: `normal_plus` / `normal_minus` / `outward_closed` |
| `obj_scale` | float | `1.0` | OBJ 読み込み後の一様スケール |
| `obj_rotation` | float[3] | `[0,0,0]` | OBJ 読み込み後の回転角 [deg] |
| `obj_offset` | float[3] | `[0,0,0]` | OBJ 読み込み後の平行移動 [m] |

`mode="auto"` では、`obj_path` が存在すれば OBJ、なければ template を使います。どちらを選んでも
設定が妥当になるよう、`auto`でもOBJ用の`surface_side`を必ず指定します。
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
| `surface_model` | string | `"insulator"` | `insulator` / `conductor` |
| `surface_side` | string | `enabled=true` で必須 | panel の真空側: `normal_plus` / `normal_minus` / `outward_closed` |
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

`dielectric`と`epsilon_r`は、分極を解かないmetadata aliasになっていたため入力から削除しました。
誘電体分極、誘電率interface条件、内部fieldは未実装です。

`conductor` の制約:

- `field_boundary.mode="free"`の直接Coulomb係数で電荷を再配分します。
- `field_boundary.mode="periodic2"`とは併用できません。
- 導体要素数が大きい場合、dense solve のバッチごとのコストが増えます。

---

### 要素 source の固定規則

要素 source は暗黙の P0 triangle panel に固定されています。旧 `[field]` table と
`sim.softening` は削除されており、入力に残すと unknown table / key として停止します。

| 項目 | 規則 |
| --- | --- |
| source | 各 `q_elem` を三角形上の一様な面電荷密度として扱う |
| solver | `direct` / `treecode` / `fmm` / `auto`。`auto` は `tree_min_nelem` で direct / FMM を選ぶ |
| 対象表面 | `insulator` / `conductor` の共通 source 離散化 |
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
| `write_potential_history` | bool | `false` | `potential_history.csv`を出力。`[domain]`があれば同じbatchの`top_reference_history.csv`も出力 |
| `dir` | string | `"outputs/latest"` | 出力先ディレクトリ |
| `history_stride` | int | `1` | 履歴 CSV の出力間隔 [batch] |
| `checkpoint_stride` | int | `0` | 再開用 checkpoint の出力間隔 [accepted batch]。`0` は定期出力なし |
| `resume` | bool | `false` | 既存 checkpoint から再開 |
| `restart_from` | string | なし | `resume=true` 時の checkpoint 読み込み元 |

出力ファイル:

| ファイル | 条件・内容 |
|---|---|
| `summary.txt` | 実行統計と設定概要 |
| `charges.csv` | 最終要素電荷 |
| `mesh_triangles.csv` | 要素 geometry。`mesh_id` 列を含む |
| `mesh_sources.csv` | `mesh_id` ごとの元メッシュ種別、表面モデル、互換用`epsilon_r`列、要素数。現行入力では`epsilon_r=1` |
| `mesh_potential.csv` | `write_mesh_potential=true` のとき |
| `charge_history.csv` | `history_stride > 0` のとき |
| `potential_history.csv` | `write_potential_history=true` かつ `history_stride > 0` のとき |
| `top_reference_history.csv` | 上記に加えて`[domain]`があるとき。z-high全断面の平均・標準偏差・最小・最大電位 |
| `matching_plane_history.csv` | `surface_current_model.model="matching_plane_quasistatic"`かつ`history_stride > 0`のとき。該当accepted batchの連成stateと残差 |
| `performance_profile.csv` | `BEACH_PROFILE=1`のとき |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | serialまたはMPI rank別の乱数状態 |
| `macro_residuals.csv` | MPIでも単一のglobalマクロ粒子数残差。species×faceを区別 |
| `charge_ledger.csv` | 粒子種別の電荷収支、粒子数、再開用累積値 |
| `checkpoint_complete.txt` | schema v8 以降の最終出力と各定期 slot の完了 manifest |
| `checkpoint_latest.txt` | `checkpoint_stride > 0` のとき。正常公開された最新 slot を示す advisory index |

各値の場所は[構成固有の出力](OutputGuide.html#構成固有の値を探す)にまとめています。

[再開に使う出力ファイル](OutputGuide.html#再開に使うファイル)に、列定義と条件別 checkpoint 要件を集約しています。

`mesh_potential.csv` の評価規約:

- 要素重心の電位 [V] を記録します。
- 自己項は解析的 P0 panel kernel で評価します。
- `periodic2` は explicit image shell を加えます。
- `cached_kneq0` は cached 非零モードと境界条件付き `k=0` を加えます。

`potential_history.csv`:

- `charge_history.csv` と同じ `history_stride` で記録します。
- 形式は `batch, elem_idx, potential_V` です。
- 履歴ごとに `field_solver%refresh` と `compute_mesh_potential` を実行するため、計算コストが増えます。

`top_reference_history.csv`:

- `potential_history.csv`と同じpost-commit snapshot、同じbatchで1行を記録します。
- 形式は`batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,potential_min_V,potential_max_V`です。
- `sample_n=reservoir.face_potential_grid_n`とし、全box z-high面のcell-centered格子を使います。
- z-high面平均は無限遠電位やプラズマ電位ではなく、相対電位を読むための診断値です。

`matching_plane_history.csv`:

- `charge_history.csv`と同じ位相で、batch 1から`history_stride`間隔のaccepted batchを1行ずつ記録します。
- 列は`batch,simulated_time_s,D_H_C_m2,phi_H_V`、electron / ionのinward fluxとaccess potential、
  PE barrier potential、4つのoutward feedback、PE return / escape flux、`iterations,residual`です。
- $D_H$と$\Phi_H$はそのbatchの粒子追跡に使ったcommit前state、`simulated_time_s`はtrial受理後の時刻です。

`resume=true` の要件:

| 条件 | 内容 |
|---|---|
| 出力 | `write_files=true` が必須 |
| 読み込み元 | `restart_from` 未指定なら `output.dir`、指定時は `restart_from` |
| 必須ファイル | `summary.txt`, `charges.csv`, serialの`rng_state.txt`またはMPI全rankの`rng_state_rankNNNNN.txt`。schema v8 以降は`checkpoint_complete.txt`も必須 |
| 条件付きファイル | ledger metadataがある場合の`charge_ledger.csv` |
| 条件付き state | schema v8 以降は`checkpoint_complete.txt`が宣言した`macro_residuals.csv`を必須とし、旧schemaでは存在時に復元 |
| 挙動 | 必須 checkpoint がなければ新規実行にフォールバックせず停止 |

`restart_from` は checkpoint の読み込み元だけを変更します。新しい出力は常に `output.dir` に書きます。

`checkpoint_stride > 0` では、accepted batch の commit 後に `checkpoints/slot0` と `slot1` を交互に更新します。
各 directory の `checkpoint_complete.txt` を全ファイルの書き込み後に原子的に完了状態へ切り替え、定期出力では続けて
`checkpoint_latest.txt` を切り替えるため、書き込み中に停止しても一つ前の slot を保持します。

再開時は `output.dir` または `restart_from` 直下の最終出力と両定期 slot を比較し、load 可能で完全な
checkpoint のうち `batches` が最大のものを選びます。`checkpoint_latest.txt` が欠落、破損、または古い場合も、
完了 manifest を持つ slot を直接検査します。最終出力は `checkpoint_stride` に関係なく従来どおり作成します。

MPI 実行時:

| ファイル | 内容 |
|---|---|
| `rng_state_rankNNNNN.txt` | rank 別乱数状態 |
| `macro_residuals.csv` | 全rankで共有するglobal残差。rootが1個だけ書く |

再開時の整合条件:

- 旧形式の `macro_residuals_rankNNNNN.csv` がある checkpoint は、暗黙変換せず拒否します。
- `summary.txt` の `mpi_world_size` は現在の rank 数と一致させます。
- 両 slot の回収時は、現行 loader が対応しない checkpoint schema を候補から除外します。
- schema v2 以降は model、ordered mesh、ordered species の fingerprint 一致が必要です。
- model fingerprint は境界eventのchord方向・離散電場work整合velocityと、表面注入位置のtrajectory契約versionも
  含みます。旧契約のcheckpointは、再開途中で運動則を切り替えないため意図的に拒否します。
- `tree_theta` / `tree_leaf_max` は値に加えて明示指定の有無も含みます。同じraw値でも、
  自動推定を使う設定と明示overrideを使う設定は異なるsolver契約です。
- schema v8 以降は `checkpoint_complete.txt` と summary の batch、MPI world size、条件付きファイル宣言が一致する必要があります。
- schema v5はneutral-return補正量、係数、未解決率を`charge_ledger.csv`から復元します。
- schema v6の`macro_residuals.csv`は`species_idx,face,residual`です。`face=0`は従来source、
  `1..6`はboundary faceを表します。旧`species_idx,residual`の2列形式も読み込めます。
- `[[particles.species]].species_key` は安定 ID です。省略時は `species_<1-based index>`、明示時は粒子種間で一意にします。

---

## 座標・配置の補助パラメータ

次のkeyも通常のTOML parameterです。ただし、読み込み時に右列の実座標・実寸を計算します。
同じ値を直接指定したときにエラーになるものと、計算値で上書きするものがあるため区別してください。

| key・指定値 | 型・条件 | 計算する値 | 直接指定したkeyとの関係 |
|---|---|---|---|
| `domain.box_origin`, `domain.box_size` | float[3]。2つを同時指定し、`box_size > 0` | `box_min = box_origin`、`box_max = box_origin + box_size` | `domain.box_min`または`domain.box_max`との併用はエラー |
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
