title: BEACH 入力パラメータリファレンス

Lang: [日本語](Parameters.md) | [English](Parameters.en.md)

# 入力パラメータリファレンス

本文書は、Fortran実行系が読む`beach.toml`のパラメータリファレンスです。
単位は、特に断りがない限り SI 単位です。

このページはすべての入力キーを掲載し、型、既定値、単位、値域、必須・排他条件、一文での効果を保持します。
モデルの導出、アルゴリズム、設計理由、出力の解釈、運用例は併記せず、専用の解説・ガイドへリンクします。

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

数値と配列の各成分は、明示的に許可されたキーがない限り有限値でなければなりません。
`*_low` / `*_high` は各軸の下限・上限です。
`inject_face` は `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, `z_high` のいずれかを指定します。

---

## 公式入門ケース

最初の実行には [10 分チュートリアル](Tutorial.html) と
`examples/tutorial_insulator.toml` を使ってください。`beachx config init`
も同一の設定を生成します。この入門ケースは、毎 batch 200 個のマクロ電子による表面電荷更新を
20 batch 追跡し、`field_solver="direct"`と`field_boundary.mode="free"`を使います。

境界reservoir流入、`plane_source`、周期境界、閉じた光電子、無限周期補正は、そのケースが実行できた後に
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
| `dt` | float | `1.0e-9` | 時間刻み [s]。有限かつ `>0` |
| `rng_seed` | int | `12345` | 乱数シード |
| `batch_count` | int | `1` | 通常実行のバッチ数。resume 時は累積到達数。`>=1` |
| `batch_duration` | float | `0.0` | 1 バッチの物理時間 [s]。有限値、流束 source では `>0` |
| `batch_duration_step` | float | `0.0` | `batch_duration=dt*batch_duration_step`。`>0`、`batch_duration` と排他 |
| `max_step` | int | `400` | 粒子 1 個あたりの最大 push 回数。`>=1` |
| `tol_rel` | float | `1.0e-8` | 相対変化量の監視値。有限かつ `>=0`、停止条件には未使用 |
| `q_floor` | float | `1.0e-30` | `rel_change` 計算時の分母下限。有限かつ `>0` |
| `multiple_box_events_policy` | string | `"abort"` | 1 step の境界event上限超過時に `abort` / `soft_discard` |
| `multiple_box_events_retry_backend` | string | `"none"` | `multiple_box_events` 後の再試行。`none` / `upper_panel_fourier` |
| `multiple_box_events_soft_discard_count_grace` | int | `1000` | 累積 soft discard 率の判定を開始する件数猶予。`>= 0` |
| `multiple_box_events_soft_discard_fraction_limit` | float | `1.0e-6` | 累積 soft discard 率の停止上限。`0 < value <= 1` |
| `multiple_box_events_soft_discard_abs_charge_limit` | float | `1.0e-12` | 累積 soft discard 絶対電荷の停止上限 [C]。有限かつ `>0` |
| `raycast_max_bounce` | int | `16` | `photo_raycast` の最大 bounce 数。有効時 `>=1` |

`batch_duration` と `batch_duration_step` の同時指定はエラーです。
`boundary_inflow` / `plane_source` / `reservoir_face` / `photo_raycast` では、解決後の `batch_duration > 0` が必須です。

`upper_panel_fourier` は `cached_kneq0` の `periodic2` 構成でのみ有効です。
`soft_discard` は、件数猶予の超過後に累積破棄率が率上限を超えた場合、または累積絶対電荷が電荷上限を超えた場合に停止します。
再試行の成立条件と診断値は[粒子の衝突・境界イベント](ParticleEvents.html#判定を完了できなければ停止する)を参照してください。

#### 外部場

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | 一様外部電場 [V/m] |
| `e0_abs` | float | 未指定 | 一様外部電場の大きさ [V/m]。`>=0` |
| `e0_phi_xy_deg` | float | `0.0` | `e0_abs` 指定時の xy 面内方位角 [deg]。明示指定には `e0_abs` が必須 |
| `e0_phi_z_deg` | float | `0.0` | `e0_abs` 指定時の xy 面からの仰角 [deg]。明示指定には `e0_abs` が必須 |
| `b0` | float[3] | `[0,0,0]` | 一様磁場 [T] |

`e0` での直接指定と、`e0_abs` と角度による指定は排他です。

#### 場ソルバ

`field_solver` は、要素電荷から Coulomb 場を評価する方式を選びます。

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
| `field_length_scale` | float | `1.0` | `field_normalization="length"` の長さ [m]。`>0` |
| `field_periodic_image_layers` | int | `1` | `periodic2` の近傍 image shell。`>=0`、`cached_kneq0` では `>=1` |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | cache 生成の Ewald 分割値。`>=0`、0 は自動 |
| `field_periodic_ewald_layers` | int | `4` | cache 生成の実・逆空間 shell 深さ。far correction 有効時 `>=1` |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | operator cacheの保存先。`cached_kneq0` では空文字不可 |
| `field_periodic_generation_tolerance` | float | `1.0e-8` | cache identity に含める生成許容誤差。`cached_kneq0` では有限かつ `>0` |
| `tree_theta` | float | 未指定時は下表の要素数 heuristic | treecode / FMM の MAC パラメータ。`0<theta<=1`、大きいほど速く粗い |
| `tree_leaf_max` | int | 未指定時は下表の要素数 heuristic | treecode / FMM の葉あたり最大 source 数。`>=1` |
| `tree_min_nelem` | int | `256` | `auto` が direct から FMM へ切り替える要素数。`>=1` |

`field_periodic_far_correction="auto"` は互換用に受理され、現在は `none` と同じ扱いです。

`field_normalization` は内部座標のみを変え、電場・電位の出力単位は SI のままです。

| `field_normalization` | 長さスケール |
|---|---|
| `si` | 入力 SI 座標をそのまま使う |
| `box` | `domain.box_max - domain.box_min`の最大幅。`[domain]`が必須 |
| `mesh` | mesh bounding box の最大幅。mesh が空なら `field_length_scale` |
| `length` | `field_length_scale` |

mode ごとの追加キー:

| `field_solver` | 評価 | 使用する追加キー |
|---|---|---|
| `direct` | 全 source を直接加算する `O(MN)` 参照計算 | なし |
| `treecode` | 遠方を monopole、近傍を panel kernel で評価 | `tree_theta`、`tree_leaf_max` |
| `fmm` | Coulomb FMM で大規模近似評価 | `tree_theta`、`tree_leaf_max`、periodic2 関連キー |
| `auto` | `tree_min_nelem` 未満は direct、以上は FMM | `tree_min_nelem`、FMM 選択時の `tree_theta`、`tree_leaf_max` |

FMM の選択と精度確認は [FMM を使う](FMM.html)、実装は
[Coulomb FMM コア詳細](FMMCore.html)を参照してください。

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
| `box_min`, `box_max` | float[3] | なし | box下限・上限 [m]。両方を指定し、全軸で `box_max > box_min` |
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

省略面は domain topology を継承し、非周期面は `open` になります。周期面は上書きできません。
反射・再分布反射・電位障壁の定義は[粒子の衝突・境界イベント](ParticleEvents.html#openreflectredistributed_reflectperiodic)を参照してください。

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

`infinity_barrier` は `boundary_inflow` だけに適用します。一様電場と併用する場合は、`phi_infty` を有効な reservoir 基準に設定します。
補正式と適用範囲は[境界から粒子を流入させる](ReservoirInjection.html#3-source_vdf-または-infinity_barrier-を選ぶ)を参照してください。

### `[surface_current_model]`: 外部シースclosure

省略時は species ごとの手動電流 target を使います。外部シースを使う場合は次の model を選びます。

| `model` | 効果 | 詳細 |
|---|---|---|
| `none` | 外部シース closure を使わない | species 側の `surface_charge_closure` を設定 |
| `zhao_stationary` | Zhao 定常根から固定電流と z-high のエネルギー障壁を決定 | [Zhao stationary closure](ZhaoStationaryClosure.html) |
| `matching_plane_quasistatic` | box 上端を外部シース応答と準定常連成 | [matching-plane 準定常連成を使う](MatchingPlaneCoupling.html) |

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `model` | string | `"none"` | `none` / `zhao_stationary` / `matching_plane_quasistatic` |
| `response_backend` | string | `"table"` | matchingの応答源。`table` / `zhao_online` |
| `zhao_branch` | string | `"auto"` | `auto` / `a` / `b` / `c`。stationaryまたはonline Zhaoのbranch。PEなしstationaryは`auto` / `c`のみ |
| `electron_species` | string | 未指定 | ambient electron の `species_key`。Zhao / matching で必須、1–64 文字 |
| `ion_species` | string | 未指定 | cold ion の `species_key`。Zhao / matching で必須、1–64 文字 |
| `photoelectron_species` | string | PE有効時に必須 | PE の `species_key`。1–64 文字、matching で省略すると PE なし |
| `solar_elevation_deg` | float | stationaryのPE有効時に必須 | Zhao sourceに使う太陽高度角 $\alpha$。$0<\alpha\le90$ degree |
| `photoelectron_ref_density_m3` | float | stationaryのPE有効時に必須 | PE基準密度 $n_{pe,ref}$ [m^-3]。`>0` |
| `photoelectron_source_scale` | float | `1.0` | stationary Zhao の $s_{UV}$。`>=0`、0 は PE なし |
| `reference_area_m2` | float | domainのx-y面積 | Zhao電流密度を総電流へ変換する面積 [m^2]。`>0`、matching では指定不可 |
| `response_table_path` | string | table matchingで必須 | 外部シース応答 CSV v1。解決後 1–256 文字、online では指定不可 |
| `implicit_zero_mode` | bool | `false` | matching の面平均 $D_H$ を後退 Euler 更新。`e_bottom_zero` が必須。table は有限 $D_H$ 軸と singleton feedback、online は CSV なしで選択 branch の終点を探索 |
| `coupling_rtol` | float | `1.0e-4` | matching固定点反復の相対収束許容値。有限な$0<r\le1$ |
| `coupling_atol` | float[4] | `[0.0, 0.0, 0.0, 0.0]` | feedback成分ごとの絶対許容値。順にPE外向きflux [m^-2 s^-1]、PE平均法線energy [eV]、electron外向きflux [m^-2 s^-1]、ion外向きflux [m^-2 s^-1]。各値は有限かつ非負、inactive成分は0 |
| `coupling_max_iterations` | int | `20` | matching固定点反復の最大回数。`>=1` |
| `coupling_relaxation` | float | `0.5` | matching更新の緩和係数。有限な$0<\omega\le1$ |

#### Zhao stationary closure

入力制約:

| 項目 | 必要条件 |
|---|---|
| role species | enabled、相互に異なる、`surface_charge_closure="fixed_current"` |
| ambient electron / ion | z-high reservoir から内向きに流入、手動 `target_*_current_a` は指定不可 |
| PE なし | `photoelectron_source_scale=0.0`、PE 固有キーは省略、`zhao_branch="auto"` または `"c"` |
| PE あり | 負電荷の `photo_raycast`、`inject_face="z_high"`、`deposit_opposite_charge_on_emit=true`、有効な z-high 境界は `open` |
| species 物性 | 単価電荷、ambient electron と PE の質量は同一、$T_e>0$、$T_{pe}>0$、$T_i\le0.1T_e$ |
| 外部場 | `sim.b0=[0,0,0]`、`reservoir.inflow_model="infinity_barrier"` は使用不可 |

PEなしType Cは$J_e+J_i=0$を満たすelectron/ion吸収targetとz-high kinetic barrier mapだけを生成し、
PE emission / return / escape targetは生成しません。
`ion_species.number_density_*` は無限遠 ion 密度です。electron 密度と PE 放出電流密度の入力は粒子分布の標本化に使い、電流 target は closure が決めます。

この model は定常電流 closure であり、box 外の過渡シースは解きません。電流・障壁・PE return の定義と出力は
[Zhao stationary closure](ZhaoStationaryClosure.html)を参照してください。

計算例は `examples/periodic2_zhao_fixed_current.toml` です。

#### Matching-plane quasistatic closure

`response_backend="table"` は外部応答 CSV、`"zhao_online"` は BEACH 内の有限 $H$ Zhao 応答を使います。
PE なしでは `photoelectron_species` を省略します。table の PE flux / energy 入力軸も 0 にします。
matching plane の $H$ は `domain.box_max` の z 成分で、面積は domain の x-y 面積です。すべての mesh 頂点を $H$ より下に置きます。

次の入力契約をすべて満たす必要があります。

| 項目 | 必要条件 |
|---|---|
| box / 場 | x/y 周期・z 非周期 open、`field_boundary.mode="periodic2"`、明示的な `[periodic2]` split 設定 |
| split | `nonzero_mode_backend` は `cached_kneq0` または `panel_spectral_reference`、`zero_mode_policy="exclude_k0"`、下側は `e_bottom_zero` または `symmetric_vacuum` |
| 外部場・開放面 | `sim.e0=sim.b0=[0,0,0]`、`ordinary_open_model="escape"`、generic reservoir potential model は使用不可 |
| event policy | `abort`、または率・件数猶予・絶対電荷上限を指定した `soft_discard` |
| role | enabled かつ相互に異なる electron / ion / 任意の PE だけを置き、`surface_charge_closure="explicit"` |
| electron / ion source | それぞれ負 / 正電荷、`source_mode="volume_seed"`、`npcls_per_step=0`、z-high の `boundary_inflow="reservoir"` のみ |
| PE source | 負電荷の `photo_raycast`、`inject_face="z_high"`、`deposit_opposite_charge_on_emit=true` |
| 粒子境界 | 全 role で x/y は `periodic`、z-low/z-high は `open` |
| 電流 target | 手動 `fixed_current` target は指定不可 |

| backend | 必須・禁止・物性制約 |
|---|---|
| `table` | `response_table_path` が必須、`zhao_branch` は指定不可 |
| `zhao_online` | `response_table_path` は指定不可、`zhao_branch` は `auto` / `a` / `b` / `c`。`implicit_zero_mode=true` では response/query CSV なしで選択 branch の終点を探索。`auto` の多重性は解消しない |
| `zhao_online` の species | 全 role は単価電荷、$T_e>0$、$0\le T_i\le0.1T_e$、ion 密度は正、electron / ion の `drift_velocity` の z 成分は負。PE 指定時は electron と同一質量かつ $T_{pe}>0$ |
| matching 共通 | stationary 専用の `solar_elevation_deg`、`photoelectron_ref_density_m3`、`photoelectron_source_scale` は指定不可 |

`model="none"` では `model` 以外を指定しません。廃止済みの `[outer_plasma]` / `[coupling]` は未対応です。
model の選択、物理的意味、適用限界は[matching-plane 準定常連成を使う](MatchingPlaneCoupling.html)、
CSV 契約、陰的更新、固定点反復は[matching-plane 数値・応答表リファレンス](MatchingPlaneReference.html)を参照してください。

### `[periodic2]`: 非零モード・零モード・下側境界

`[periodic2]`はトップレベルtableです。`domain.periodic_axes=["x","y"]`と
`field_boundary.mode="periodic2"`を指定します。productionでは`field_solver="fmm"`を使い、
小規模検証用のsplit referenceに限り`field_solver="direct"`を使います。

| キー | 型 | 既定値 | 意味・制約 |
|---|---|---:|---|
| `nonzero_mode_backend` | string | 必須 | `panel_spectral_reference` / `cached_kneq0` |
| `zero_mode_policy` | string | 必須 | `exclude_k0` |
| `lower_boundary_model` | string | 必須 | `symmetric_vacuum` / `e_bottom_zero` |
| `max_nonzero_mode_potential_step` | float | `0` | accepted trial の $k\ne0$ 電位変化上限 [V]。`>=0`、0 は無効、`cached_kneq0` 専用 |
| `reference_mode_layers` | int | `4` | Fourier mode cutoff。`>=1` |
| `panel_quadrature_order` | int | `12` | panel 面積積分次数。`>=2` |

適応的な電位変化上限は `boundary_inflow` / `plane_source` / `reservoir_face` / `photo_raycast` に対応します。
最初の 3 source は `target_macro_particles_per_batch` が必須で固定 `w_particle` は使用不可、正の `volume_seed` も併用不可です。
`sim.batch_count` は accepted batch 数、`simulated_time_s` は受理幅の総和です。受理条件、rollback、収束確認は
[`batch_duration` をどう決めるか](BatchDurationStability.html)を参照してください。

### periodic2の組合せ制約

periodic2では`[domain]`、`periodic_axes=["x","y"]`、`field_boundary.mode="periodic2"`が必要です。
`examples/periodic2_closed_photoelectron.toml`は、x/y周期、境界reservoir、閉じた光電子を組み合わせた基準例です。
同じ周期条件をfield、collision、`photo_raycast`に適用します。

| `nonzero_mode_backend` | 意味 |
| --- | --- |
| `panel_spectral_reference` | 小規模なsplit reference |
| `cached_kneq0` | versioned operatorを再利用するproduction非零モード |

非零・零modeの分割、`exclude_k0`の役割、下側境界から決まる平均場は
[periodic2静電場](PeriodicElectrostatics.html)を参照してください。

### `[[particles.species]]`: 粒子種

`[[particles.species]]` は 1 件以上必須です。
`source_mode` によって、使うキーと制約が変わります。

#### 共通キー

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `species_key` | string | `"species_<1-based index>"` | 安定 ID。1–64 文字、粒子種間で一意 |
| `enabled` | bool | `true` | 種を有効化 |
| `source_mode` | string | `"volume_seed"` | `volume_seed` / `plane_source` / `photo_raycast` / deprecated `reservoir_face` |
| `q_particle` | float | `-1.602176634e-19` | 粒子電荷 [C]。enabled speciesでは非0 |
| `m_particle` | float | `9.10938356e-31` | 粒子質量 [kg]。`>0` |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | 位置下限 [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | 位置上限 [m] |
| `drift_velocity` | float[3] | `[0,0,-8e5]` | ドリフト速度 [m/s] |
| `temperature_k` | float | `2.0e4` | 温度 [K]。`>=0` |
| `temperature_ev` | float | 未指定 | 温度 [eV]。`>=0`、`temperature_k` と排他 |
| `velocity_distribution` | string | `"maxwellian"` | `maxwellian` / `grid` |
| `inject_face` | string | 未指定 | `photo_raycast`の照射開口面。deprecatedな`reservoir_face`でも必須 |
| `source_normal` | float[3] | 未指定 | `plane_source`の一方向法線。axis-alignedな非ゼロベクトル |
| `boundary` | table | 未指定 | `[particles.species.boundary]` のspecies別6面override |
| `boundary_inflow` | table | 未指定 | `[particles.species.boundary_inflow]`のspecies別reservoir流入面 |
| `surface_charge_closure` | string | `"explicit"` | 表面 source 電荷 closure。`explicit` / `fixed_current` / `neutral_return` |
| `target_absorbed_current_a` | float | 未指定 | `fixed_current` の signed 吸収電流 [A]。0または`q_particle`と同符号 |
| `target_emission_current_a` | float | 未指定 | `fixed_current` の signed 放出反作用電流 [A]。0または`q_particle`と逆符号 |

流束駆動 source（`boundary_inflow`、`plane_source`、旧 `reservoir_face`）の共通キー:

| キー | 型 | 既定値 | 説明・制約 |
|---|---|---:|---|
| `number_density_cm3` | float | 未指定 | Maxwell 分布の上流密度 [cm^-3]。`>0`、`number_density_m3` と排他 |
| `number_density_m3` | float | 未指定 | Maxwell 分布の上流密度 [m^-3]。`>0`、`number_density_cm3` と排他 |
| `w_particle` | float | `target_macro_particles_per_batch` といずれか一方必須 | マクロ粒子重み。`>0` |
| `target_macro_particles_per_batch` | int | `w_particle` といずれか一方必須 | マクロ粒子数 target。`>0`、または species 2 以降で species 1 の重みを共有する `-1` |
| `velocity_grid_path` | string | 未指定 | `velocity_distribution="grid"` の非空 CSV path。列は `vx_m_s,vy_m_s,vz_m_s,f` |
| `velocity_grid_pdf_kind` | string | `"phase_space"` | `phase_space` / `flux_weighted` |
| `velocity_grid_sampling` | string | `"auto"` | `auto` / `rectilinear` / `discrete` |
| `particle_flux_m2_s` | float | 未指定 | grid 分布の入射粒子 flux [m^-2 s^-1]。`>0`、`current_density_a_m2` と排他 |
| `current_density_a_m2` | float | 未指定 | grid 分布の入射電流密度 [A/m^2]。非0、`particle_flux_m2_s` と排他 |

`w_particle` と `target_macro_particles_per_batch` はちょうど一方を明示します。
Maxwell 分布では密度2形式のちょうど一方と温度を、grid 分布では CSV と flux 2形式のちょうど一方を指定します。
分布の意味と CSV の sampling は[境界から粒子を流入させる](ReservoirInjection.html#2-maxwell-分布または速度-grid-を選ぶ)を参照してください。

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

| `surface_charge_closure` | 効果 | 必要条件 |
|---|---|---|
| `explicit` | 追跡した電荷をそのまま反映 | 既定 |
| `neutral_return` | PE による表面総電荷増分を 0 に補正 | 負電荷 `photo_raycast`、放出反作用あり、注入面は反射、escape / `soft_discard` なし、未帰還率 $\le5\%$ |
| `fixed_current` | 要素別の標本分布を保ち、総電流を target へ倍率化 | 正の `batch_duration`、手動 `target_*_current_a` または `[surface_current_model]`、非ゼロ target に対応する raw channel |

`fixed_current` と `neutral_return` は同じ species では排他です。
補正式、PE return の二重計上を避ける条件、標本数の収束確認は
[表面電荷更新の数値仕様](SurfaceChargeNumerics.html)、[光電子の放出とライフサイクル](PhotoelectronEmission.html)、
[出力形式リファレンス](OutputReference.html#charge-ledger)を参照してください。

#### `source_mode = "volume_seed"`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `npcls_per_step` | int | `0` | 1 バッチに生成するマクロ粒子数。`>=0` |
| `w_particle` | float | `1.0` | マクロ粒子重み。`>0` |

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

`reservoir` は選択した非周期 box 面全体から流入させます。対象面の有効な粒子境界は `open` でなければなりません。
`source_mode="volume_seed"` とのみ併用でき、複数面の `target_macro_particles_per_batch` は全流入面の合計です。
流束と電位障壁の詳細は[境界から粒子を流入させる](ReservoirInjection.html)を参照してください。

#### `source_mode = "plane_source"`

`plane_source` は box 内部の axis-aligned 矩形面から `source_normal` 方向へ流束を生成します。

| 条件 | 内容 |
|---|---|
| 領域 | `[domain]`が必須 |
| 時間 | `sim.batch_duration > 0`が必須 |
| 矩形面 | `pos_low` / `pos_high`はちょうど1軸で一致し、残る2軸は正の長さを持つ |
| 配置 | 法線座標はbox境界より厳密に内側。接線範囲はbox内で、境界との一致を許容 |
| 方向 | `source_normal`はzero-thickness軸に沿う非ゼロベクトル。内部で正規化し、入力は正負単位ベクトルを推奨 |
| 外部補正 | `[reservoir]`の`infinity_barrier`、`phi_infty`、`face_potential_grid_n`は適用しない |

流束と速度分布は上記の共通キーで指定します。

#### `source_mode = "reservoir_face"`（deprecated）

キーは上記の流束駆動 source 共通表を使います。互換入力の追加制約:

| 条件 | 内容 |
|---|---|
| 領域 | `[domain]`が必須 |
| 時間 | `sim.batch_duration > 0` が必須 |
| 注入面 | `inject_face` が必須 |
| 注入範囲 | `pos_low` / `pos_high` は指定 face 上にある必要がある |
| 重み | `w_particle` と `target_macro_particles_per_batch` は同時指定不可 |
| 重み共有 | `target_macro_particles_per_batch=-1` は species 2 以降だけ可。species 1 の `w_particle` を共有 |

新しいケースでは、外部 plasma に `boundary_inflow`、内部矩形面に `plane_source` を使います。
BEACH は旧 mode を暗黙変換しません。CSV と重みの契約は[境界から粒子を流入させる](ReservoirInjection.html#6-旧-reservoir_face-から移行する)を参照してください。

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

レイ重み、周期 image、再吸収、closed PE の設定は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)を参照してください。

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
| `history_stride` | int | `1` | 履歴 CSV の出力間隔 [batch]。`>=0`、0 は無効 |
| `checkpoint_stride` | int | `0` | 再開用 checkpoint の間隔 [accepted batch]。`>=0`、0 は定期出力なし、正値では `write_files=true` |
| `resume` | bool | `false` | 既存 checkpoint から再開 |
| `restart_from` | string | なし | `resume=true` 時の checkpoint 読み込み元。指定時は `write_files=true` |

生成条件、列定義、電位の評価規約、matching-plane の状態、ledger の解釈は
[出力形式リファレンス](OutputReference.html)を参照してください。

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

定期 slot の選択、MPI の必須ファイル、fingerprint と schema の互換条件は
[再開に使うファイル](OutputReference.html#再開に使うファイル)を参照してください。

---

## 座標・配置の補助パラメータ

次の入力キーは box 基準の値を実座標・実寸へ変換します。

| キー | 型 | 既定値 | 効果・制約 |
|---|---|---:|---|
| `domain.box_origin` | float[3] | 未指定 | `box_min` を設定。`box_size` と組にし、`box_min` / `box_max` と排他 |
| `domain.box_size` | float[3] | 未指定 | `box_max=box_origin+box_size`。全成分 `>0` |
| `inject_region_mode` | string | `"absolute"` | `absolute` / `face_fraction`。`face_fraction` は `reservoir_face` / `photo_raycast` 専用 |
| `uv_low`, `uv_high` | float[2] | 未指定 | `face_fraction` で両方必須。各成分は `[0,1]`、`pos_low` / `pos_high` と排他 |
| template `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| template `anchor` | string | 未指定 | `box_anchor` で使う box 中心または 6 面中心 |
| template `offset` | float[3] | 未指定 | anchor からの差 [m]。`offset_frac` と排他 |
| template `offset_frac` | float[3] | 未指定 | box 幅に対する差。`offset` と排他 |
| template `size_mode` | string | `"absolute"` | `absolute` / `box_fraction` |
| template `size_frac` | float / float[2] / float[3] | 未指定 | `box_fraction` で必須。全成分 `>0`、kind により下表の寸法を上書き |
| template `group` | string | 未指定 | `[mesh.groups.<name>]` の name。空文字不可 |
| template `center_local` | float[3] | 未指定 | group 使用時に必須。`center=group_origin+group_scale*center_local` |
| group `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| group `anchor` | string | 未指定 | `box_anchor` で使う box 中心または 6 面中心 |
| group `offset` | float[3] | 未指定 | group 原点へ加える差 [m]。`offset_frac` と排他 |
| group `offset_frac` | float[3] | 未指定 | box 幅に対する差。`offset` と排他 |
| group `scale` | float | `1.0` | group 座標と明示寸法の倍率。`>0`、`scale_from` / `scale_factor` と排他 |
| group `scale_from` | string | 未指定 | box 幅の参照名。`scale_factor` と組で指定 |
| group `scale_factor` | float | 未指定 | `scale_from` に掛ける正値。`>0` |

group 使用時は、template の `center`、直接 placement キー、`size_mode`、`size_frac` を併用できません。

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
