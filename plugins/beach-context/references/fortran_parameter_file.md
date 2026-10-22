title: BEACH 入力パラメータリファレンス

# 入力パラメータリファレンス

設定ファイルは `run / domain / mesh / particles / fields / sheath / output` の7グループ形式だけを受け付けます。完全な入力例は `examples/` を参照してください。

`beach.toml` のキーを調べるための一覧です。型、既定値、単位、必須条件と排他条件を記載します。
単位は特記しない限り SI です。最初の設定は[ケース設計](../guide/ConfigurationRecipes.md)、
編集後の確認は[設定の検証](../guide/Configuration.md)から始めてください。

| 調べる対象 | このページの項目 |
|---|---|
| 実行制御と再開 | [`run`](#run) |
| 領域と周期軸 | [`domain`](#domain) |
| 表面メッシュ | [`mesh`](#mesh) |
| 粒子追跡・粒子源・境界 | [`particles`](#particles) |
| 場ソルバ・周期場・外部場 | [`fields`](#fields) |
| 外部シース | [`sheath`](#sheath) |
| 出力と診断 | [`output`](#output) |

設定ファイルへ全項目をコピーする必要はありません。使う粒子源・境界・solver に必要な項目を選び、
省略した項目には表の既定値を使います。「必須」「未指定」と書かれた項目は、各表の適用条件も確認してください。

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
| 検証 | `beachx lint beach.toml` は同梱 schema でキー・型・値域を、`beach --check-config` で組合せと派生値を検証する。組合せ制約の正本は Fortran の preflight |
| Fortran の検証 | 開発・診断用の `beach --check-config beach.toml` は読込・正規化・モデルの最小条件を検査し、パスが必須。実行前の lint の代わりにはしない |
| 基本値 | 実数と配列成分は有限値。整数項目に実数・真偽値は指定不可 |
| 整数の格納範囲 | 符号付き 32 bit の `-2147483648..2147483647`。範囲外は変換前に拒否。各項目固有の値域は lint で検証 |
| 文字列 | Fortran の格納先を超える長さはエラー。暗黙の切り詰めはしない |

Editor schema を使う場合は、`beach.toml` の先頭に GitHub Raw URL のコメント directive を置きます。
Fortran パーサは最初のセクションより前の通常キーを受け付けないため、`"$schema" = "..."` は使いません。

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

---

## 単位と座標

| 種類 | 代表キー | 単位・向き |
|---|---|---|
| 時間 | `particles.tracking.dt_s`, `run.batch.duration_s` | 秒 |
| 長さ | `domain.box_min`, `domain.box_max`, `source.pos_low`, `source.pos_high` | m |
| 電荷 | `charge_c`, 要素電荷出力 | C |
| 質量 | `mass_kg` | kg |
| 速度 | `distribution.drift_velocity_m_s`, `source.ray_direction` | m/s。ただし `source.ray_direction` は方向ベクトル |
| 電場 | `fields.external.electric_v_m`, `fields.external.electric_magnitude_v_m` | V/m |
| 磁場 | `fields.external.magnetic_t` | T |
| 密度 | `distribution.number_density_cm3`, `distribution.number_density_m3` | cm^-3 または m^-3 |
| 温度 | `distribution.temperature_k`, `distribution.temperature_ev` | K または eV。両方の同時指定は不可 |
| 角度 | `fields.external.electric_azimuth_deg`, `fields.external.electric_elevation_deg` | degree |

数値と配列の各成分は有限値でなければなりません。無効化した species も、入力された値の型・有限性・
文字列長の検査を受けます。整数の格納範囲も同じです。
列挙値、宣言上の値域、無効な機能への指定などの詳細な入力診断は Schema / Python が担当します。
Fortran はその全診断を繰り返さず、派生値・参照解決と実行に必要な最小条件を確認します。
不正入力の採否は一致を保証しません。
物理的な組合せ制約は各項目の適用条件に従います。
`*_low` / `*_high` は各軸の下限・上限です。
`source.inject_face` は `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, `z_high` のいずれかを指定します。

---

## 公式入門ケース

最初の実行には [入門チュートリアル](../start/Tutorial.md) と
`examples/tutorial_electron_beam.toml` を使ってください。`beachx config init`
も同一の設定を生成します。この入門ケースは、5 eV の電子ビームによる浮遊導体板の帯電を
150 batch（1.5 μs）追跡し、`fields.solver.method="direct"`と`fields.boundary="free"`を使います。

境界reservoir流入、周期境界、閉じた光電子、無限周期補正は、そのケースが実行できた後に
[シミュレーションケースを設計する](../guide/ConfigurationRecipes.md)から追加する応用設定です。

---

## TOML の階層とセクション一覧

公開設定は 7 グループで構成します。`run.batch.duration_s` は「table.key」の参照表記です。

```text
beach.toml
├── [run]
│   ├── [run.batch]
│   │   └── [run.batch.adaptive]
│   └── [run.restart]
├── [domain]
├── [mesh]
│   ├── [mesh.groups.<name>]
│   └── [[mesh.parts]]
├── [particles]
│   ├── [particles.tracking]
│   │   └── [particles.tracking.events]
│   ├── [particles.raycast]
│   ├── [particles.boundary]
│   ├── [particles.reservoir]
│   ├── [particles.secondary_emission]
│   └── [[particles.species]]
│       ├── [particles.species.distribution]
│       ├── [particles.species.source]
│       ├── [particles.species.sampling]
│       ├── [particles.species.charging]
│       ├── [particles.species.boundary]
│       ├── [particles.species.inflow]
│       └── [particles.species.surface_interaction]
│           └── [particles.species.surface_interaction.secondary_emission]
├── [fields]
│   ├── [fields.solver]
│   │   └── [fields.solver.tree]
│   ├── [fields.periodic]
│   │   ├── [fields.periodic.upper_fourier]
│   │   └── [fields.periodic.reference]
│   └── [fields.external]
├── [sheath]
│   ├── [sheath.zhao]
│   ├── [sheath.species]
│   ├── [sheath.photoelectrons]
│   └── [sheath.coupling]
└── [output]
    ├── [output.diagnostics]
    ├── [output.final]
    ├── [output.history]
    └── [output.checkpoint]
```

| TOML table | 親 | 件数・必須条件 | 内容 |
|---|---|---|---|
| `[run]` | root | 条件付き | 乱数シードとバッチ数 |
| `[run.batch]` | `[run]` | 流束源で正の時間幅が必須 | バッチの物理時間 |
| `[run.batch.adaptive]` | `[run.batch]` | 任意、`cached_kneq0` 専用 | 非零モード電位変化の上限 |
| `[run.restart]` | `[run]` | 再開時に指定 | checkpoint の読み込み元 |
| `[domain]` | root | box 使用時 | box geometry と周期軸 |
| `[mesh]` | root | 任意 | OBJ と template の要素リスト |
| `[mesh.groups.<name>]` | `[mesh]` | 0 件以上 | 複数 template で共有する配置と scale |
| `[[mesh.parts]]` | `[mesh]` | 0 件以上 | OBJ または組み込み形状 |
| `[particles]` | root | 必須 | 粒子種・追跡・境界・SEE 世代制御 |
| `[particles.tracking]` | `[particles]` | 任意 | 粒子の時間刻みと追跡長 |
| `[particles.tracking.events]` | `[particles.tracking]` | 任意 | 境界イベントの再試行と破棄 |
| `[particles.raycast]` | `[particles]` | 任意 | 光電子 ray の bounce 上限 |
| `[particles.boundary]` | `[particles]` | 任意 | 非周期面の global 粒子作用 |
| `[particles.reservoir]` | `[particles]` | 任意 | 外部 reservoir の流入障壁と基準電位 |
| `[particles.secondary_emission]` | `[particles]` | 任意 | 二次電子放出の世代上限 |
| `[[particles.species]]` | `[particles]` | 1 件以上 | 粒子種の識別子・電荷・質量 |
| `[particles.species.distribution]` | `[[particles.species]]` | 粒子分布の指定時 | 密度・温度・ドリフト |
| `[particles.species.source]` | `[[particles.species]]` | 独立した粒子源の使用時 | 生成・表面放出方式と位置 |
| `[particles.species.sampling]` | `[[particles.species]]` | 独立した粒子源・境界流入の使用時 | マクロ粒子数・重み・光線数 |
| `[particles.species.charging]` | `[[particles.species]]` | 任意 | 表面電流の補正方式と目標電流 |
| `[particles.species.boundary]` | `[[particles.species]]` | 任意 | その粒子種だけの非周期面 override |
| `[particles.species.inflow]` | `[[particles.species]]` | 任意 | 外部 reservoir から流入する非周期面 |
| `[particles.species.surface_interaction]` | `[[particles.species]]` | 任意 | 吸収時の表面作用 |
| `[particles.species.surface_interaction.secondary_emission]` | `[particles.species.surface_interaction]` | 任意 | 入射種の二次電子放出条件 |
| `[fields]` | root | 任意 | `free` / `periodic2` 場境界 |
| `[fields.solver]` | `[fields]` | 任意 | 場ソルバと長さの正規化 |
| `[fields.solver.tree]` | `[fields.solver]` | 任意 | Treecode / FMM の精度・木の分割 |
| `[fields.periodic]` | `[fields]` | periodic2 使用時 | 周期 backend と下側境界 |
| `[fields.periodic.upper_fourier]` | `[fields.periodic]` | 任意 | 上部真空の Fourier 切替高さ・次数 |
| `[fields.periodic.reference]` | `[fields.periodic]` | split reference 使用時 | Fourier cutoff と面積分 |
| `[fields.external]` | `[fields]` | 任意 | 一様外部電場・磁場 |
| `[sheath]` | root | 任意 | 外部シースの零電流根による closure |
| `[sheath.zhao]` | `[sheath]` | 零電流 closure 使用時 | Zhao 根の branch と受理条件 |
| `[sheath.species]` | `[sheath]` | 零電流 closure 使用時 | 電子・イオン・PE・SEE の役割 |
| `[sheath.photoelectrons]` | `[sheath]` | 零電流 closure 使用時 | 外部シースの光電子源 |
| `[sheath.coupling]` | `[sheath]` | 任意 | 観測流出による外部根の更新 |
| `[output]` | root | 任意 | 出力の有効化と出力先 |
| `[output.diagnostics]` | `[output]` | 任意 | 電荷変化の監視値 |
| `[output.final]` | `[output]` | 任意 | 最終電位の出力 |
| `[output.history]` | `[output]` | 任意 | 履歴の種類と間隔 |
| `[output.checkpoint]` | `[output]` | 任意 | 定期 checkpoint の間隔 |

---

## パラメータ詳細リファレンス

### `[run]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `rng_seed` | int | `12345` | 乱数シード |
| `batch_count` | int | `1` | 通常実行のバッチ数。resume 時は累積到達数。`>=1` |

#### `[run.batch]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `duration_s` | float | `0.0` | 1 バッチの物理時間 [s]。有限値、流束 source では `>0` |
| `duration_steps` | float | `0.0` | `run.batch.duration_s=particles.tracking.dt_s*run.batch.duration_steps`。`>0`、`run.batch.duration_s` と排他 |
| `particle_distribution` | string | `"static"` | バッチ粒子の MPI rank 配分。バッチの粒子は種・source・面ごとの global slot に並び、各 slot は (seed, batch, species, source, face, index) の PCG64 counter-based 乱数で生成されるので、どちらでもバッチの粒子集合は rank 数に依らない。`static` は slot を rank に連続区間で割り振る（決定的）。`shared` は全 rank 共有の slot pool を rank 0 のカウンタ（MPI RMA の fetch-and-add）から `particle_pool_chunk` 個ずつ取り合い、rank 間の粒子追跡負荷を動的に均す（`MPI_THREAD_SERIALIZED` を使う）。粒子ごとの境界・二次電子の乱数も global slot で決まる。`summary.txt` に `particle_distribution=`、`particle_pool_claims=`、`particle_pool_claim_seconds=` を記録する |
| `particle_pool_chunk` | int | `64` | `shared` で 1 回の claim が取る slot 数。`>=1` |

`run.batch.duration_s` と `run.batch.duration_steps` の同時指定はエラーです。
`inflow` / `photo_raycast` では、解決後の `run.batch.duration_s > 0` が必須です。
刻み・追跡長・粒子数を一緒に設計する手順は[ケース設計](../guide/ConfigurationRecipes.md#7-時間と粒子数を決める)を参照してください。

##### `[run.batch.adaptive]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `max_nonzero_mode_potential_step_v` | float | `0` | accepted trial の $k\ne0$ 電位変化上限 [V]。`>=0`、0 は無効、`cached_kneq0` 専用 |

適応的な電位変化上限は `inflow` / `photo_raycast` に対応します。
境界流入は `sampling.target_macro_particles_per_batch` が必須で固定 `sampling.weight` は使用不可、正の `volume_seed` も併用不可です。
`run.batch_count` は accepted batch 数、`simulated_time_s` は受理幅の総和です。受理条件、rollback、収束確認は
[`run.batch.duration_s` をどう決めるか](../guide/BatchDurationStability.md)を参照してください。

#### `[run.restart]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `from` | string | `output.dir` | checkpoint 読み込み元。`output.enabled=true` が必須。新しい出力は `output.dir` へ書く |

`[run.restart]` を書くと再開を有効にします。空の table では `output.dir` から読み込みます。

`[run.restart]` の要件:

| 条件 | 内容 |
|---|---|
| 出力 | `output.enabled=true` が必須 |
| 読み込み元 | `run.restart.from` 未指定なら `output.dir`、指定時は `run.restart.from` |
| 必須ファイル | `summary.txt`, `charges.csv`。schema v8 以降は`checkpoint_complete.txt`も必須。schema v10 以前は serialの`rng_state.txt`またはMPI全rankの`rng_state_rankNNNNN.txt`の存在も確認する（内容は使わない） |
| MPI rank 数 | 保存時と違ってよい（乱数は状態を持たない counter-based の PCG64） |
| 条件付きファイル | ledger metadataがある場合の`charge_ledger.csv` |
| 条件付き state | schema v8 以降は`checkpoint_complete.txt`が宣言した`macro_residuals.csv`を必須とし、旧schemaでは存在時に復元 |
| 挙動 | 必須 checkpoint がなければ新規実行にフォールバックせず停止 |

`run.restart.from` は checkpoint の読み込み元だけを変更します。新しい出力は常に `output.dir` に書きます。

定期 slot の選択、MPI の必須ファイル、mesh 識別子と schema の互換条件は
[再開に使うファイル](OutputReference.md#再開に使うファイル)を参照してください。

### `[domain]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `box_min`, `box_max` | float[3] | なし | box下限・上限 [m]。両方を指定し、全軸で `box_max > box_min` |
| `box_origin`, `box_size` | float[3] | なし | originと正のsize [m]。両方を指定し、`box_max=box_origin+box_size`、全成分 `box_size>0` |
| `periodic_axes` | string[] | `[]` | 重複しない`"x"` / `"y"` / `"z"` |

2組のgeometry表現は同時指定できません。`[domain]`を使う場合はどちらか1組が必須です。
周期性を指定できる公開キーは`domain.periodic_axes`だけです。
`fields.boundary="periodic2"`では、`periodic_axes=["x","y"]`が必要です。

### `[mesh]`

`parts` は 0 件以上の table 配列で、OBJ と組み込み template をリスト順に混在させられます。
`groups` は template kind 専用の共通配置と scale を定義します。
有効な要素をリスト順に 1 から連番の `mesh_id` へ割り当て、disabled 要素は番号を消費しません。
1 つの OBJ 要素の全三角形は、ファイル内で形状が分離していても同じ番号を持ちます。
`surface_model="conductor"` の要素は `mesh_id` ごとに 1 つの浮遊導体 group になります。独立導体は別々の parts に分けます。
有効な要素が 0 件の設定からは mesh を構築できません。

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `parts` | table 配列 | `[]` | OBJ と template の順序付きリスト |
| `groups` | 名前付き table | なし | template 共通の原点と scale |

#### `[mesh.groups.<name>]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| `anchor` | string | 未指定 | `box_anchor` で使う box 中心または 6 面中心 |
| `offset` | float[3] | 未指定 | group 原点へ加える差 [m]。`offset_frac` と排他 |
| `offset_frac` | float[3] | 未指定 | box 幅に対する差。`offset` と排他 |
| `scale` | float | `1.0` | group 座標と明示寸法の倍率。`>0`、`scale_from` / `scale_factor` と排他 |
| `scale_from` | string | 未指定 | box 幅の参照名。`scale_factor` と組で指定 |
| `scale_factor` | float | 未指定 | `scale_from` に掛ける正値。`>0` |

#### `[[mesh.parts]]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `enabled` | bool | `true` | この要素を有効化 |
| `kind` | string | `"plane"` | `obj` / `plane` / `plate_hole` / `plane_hole` / `disk` / `annulus` / `box` / `cylinder` / `sphere` |
| `surface_model` | string | `"insulator"` | `insulator` / `conductor` |
| `surface_side` | string | `enabled=true` で必須 | panel の真空側: `normal_plus` / `normal_minus` / `outward_closed` |
| `center` | float[3] | `[0,0,0]` | 形状中心 [m] |
| `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| `anchor` | string | 未指定 | `box_anchor` で使う box 中心または 6 面中心 |
| `offset` | float[3] | 未指定 | anchor からの差 [m]。`offset_frac` と排他 |
| `offset_frac` | float[3] | 未指定 | box 幅に対する差。`offset` と排他 |
| `size_mode` | string | `"absolute"` | `absolute` / `box_fraction` |
| `size_frac` | float / float[2] / float[3] | 未指定 | `box_fraction` で必須。全成分 `>0`、kind により下表の寸法を上書き |
| `group` | string | 未指定 | `[mesh.groups.<name>]` の name。空文字不可 |
| `center_local` | float[3] | 未指定 | group 使用時に必須。`center=group_origin+group_scale*center_local` |

共通キーは `kind`、`enabled`、`surface_model`、`surface_side` です。表の `center` から `center_local` までは、`offset` を除いて template kind 専用です。
`obj` に形状・group・placement キーを指定した場合、および template kind に `path` / `scale` / `rotation` を指定した場合は preflight エラーです。
`offset` は OBJ では平行移動、template では従来どおり anchor からの差を表します。

`kind` の概要:

| `kind` | 生成形状 | 基準面・軸 |
|---|---|---|
| `obj` | OBJ の三角形 | 要素ごとの scale・rotation・offset |
| `plane` | 長方形平面 | XY 平面、`z=center[3]` |
| `plate_hole`, `plane_hole` | 中央に円形穴を持つ長方形平面 | XY 平面、穴中心は `center` |
| `disk` | 円板 | XY 平面、中心は `center` |
| `annulus` | 同心リング | XY 平面、中心は `center` |
| `box` | 閉じた直方体表面 | 各軸に平行な 6 面 |
| `cylinder` | 円柱側面と任意の上下キャップ | z 軸方向 |
| `sphere` | 球面 | 中心は `center` |

##### `kind = "obj"`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `path` | string | 必須 | OBJ パス。相対パスは実行時のカレントディレクトリから解決 |
| `scale` | float | `1.0` | 一様スケール |
| `rotation` | float[3] | `[0,0,0]` | x→y→z の外因性回転 [deg] |
| `offset` | float[3] | `[0,0,0]` | 平行移動 [m] |

変換は `v_new = R(rotation) * (v_old * scale) + offset` で、この要素だけに適用します。
LF / CRLF、面行 `f v` / `f v/vt` / `f v/vt/vn` / `f v//vn` に対応し、多角形はファン三角形分割します。
平面と導体の四面体 OBJ の混在例は [`examples/mixed_parts.toml`](../../examples/mixed_parts.toml) です。

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

誘電体分極、誘電率interface条件、内部fieldは未実装です。

`conductor` の制約:

- その batch の場 snapshot の電位で等電位化します。free-space と `periodic2` の両方で使えます。
- 影響行列は初期化時に導体要素数回の場 refresh で組みます。導体要素数が大きい場合、初期化と dense solve のコストが増えます。

### `[particles]`

`[[particles.species]]` は 1 件以上必要です。

#### `[particles.tracking]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `dt_s` | float | `1.0e-9` | 時間刻み [s]。有限かつ `>0` |
| `max_steps_per_particle` | int | `400` | 粒子 1 個あたりの最大 push 回数。`>=1` |
| `substep_field_change_max` | float | `0` | 粒子ごとの適応サブステップ。step 候補の中点電場が前回受理した中点電場からこの割合より大きく変わり、その差 × 変位が 1e-4 eV を超えるときは、その step を半分にして取り直す。帯電した要素の近傍で 1 step の変位が場の変化スケールに近い case で、エネルギー保存を保つ。`0` で無効、有限かつ `>=0`。`0.1`〜`0.2` から試し、`summary.txt` の `adaptive_substep_*` で refinement の量を確認する |
| `substep_max_halvings` | int | `6` | 1 つの `dt_s` を半分にする回数の上限（最小 sub-step は `dt_s/2^n`）。`0..16` |
| `step_growth` | bool | `false` | 場の変化が小さい領域で step を `dt_s` より長くする。受理した step の中点電場の変化が `substep_field_change_max` の 4 分の 1 未満（または変化 × 変位が床値以下）なら次の step を 2 倍にし、変化が大きければ既存の半分化で戻す。伸長時の床値は 1e-4 eV と粒子の運動エネルギーの `step_growth_energy_fraction` 倍の大きい方（遠方の弱い場では FMM 打ち切り雑音が場と同程度で、絶対床値だけでは雑音が半分化を起こす）。半分化の後 8 chunk は伸長を試みない。`substep_field_change_max` が `0` のときは `0.05` を自動で使う。衝突判定は線分なので長い step でも面をすり抜けない。粒子層から離れた弱い場の領域（本番相当の入力では step の 4 分の 3）で効き、長寿命の捕捉電子の追跡時間を大きく縮める。`summary.txt` の `adaptive_grown_steps` が伸ばした push の数 |
| `step_growth_max` | int | `5` | `step_growth` で 1 つの step を 2 倍にする回数の上限（最大 step は `dt_s*2^n`）。`1..16` |
| `step_growth_energy_fraction` | float | `0.01` | `step_growth` の仕事の床値を粒子の運動エネルギーのこの割合にする（1e-4 eV との大きい方）。本番相当の入力では `0.001` で粒子追跡が 1.9 倍、`0.01` で 4.4 倍速くなり、吸収・脱出数と表面電荷は実行間の揺らぎの範囲で一致した。`(0, 1]` |

##### `[particles.tracking.events]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `policy` | string | `"abort"` | 1 step の境界event上限超過時に `abort` / `soft_discard` |
| `retry_backend` | string | `"none"` | `multiple_box_events` 後の再試行。`none` / `upper_panel_fourier` |
| `discard_count_grace` | int | `1000` | 累積 soft discard 率の判定を開始する件数猶予。`>= 0` |
| `discard_fraction_limit` | float | `1.0e-6` | 累積 soft discard 率の停止上限。`0 < value <= 1` |
| `discard_charge_warning_c` | float | `1.0e-12` | 累積 soft discard 絶対電荷の警告閾値 [C]。有限かつ `>0` |

境界通過の再試行・破棄、または光電子の ray tracing を調整する場合に使います。

`upper_panel_fourier` は `cached_kneq0` の `periodic2` 構成でのみ有効です。
`soft_discard` は、件数猶予の超過後に累積破棄率が率上限を超えた場合だけ停止します。
累積絶対電荷は実行時間とともに増えるため、電荷閾値の初回超過は警告として記録します。
再試行の成立条件と診断値は[粒子の衝突・境界イベント](../developer/ParticleEvents.md#判定を完了できなければ停止する)を参照してください。

#### `[particles.raycast]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `max_bounces` | int | `16` | `photo_raycast` の最大 bounce 数。有効時 `>=1` |

#### `[particles.boundary]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `x_low`, `x_high` | string | domainに従う | 非周期x面の`open` / `reflect` / `redistributed_reflect` |
| `y_low`, `y_high` | string | domainに従う | 非周期y面の`open` / `reflect` / `redistributed_reflect` |
| `z_low`, `z_high` | string | domainに従う | 非周期z面の`open` / `reflect` / `redistributed_reflect` |

省略面は domain topology を継承し、非周期面は `open` になります。周期面は上書きできません。
反射・再分布反射の定義は[粒子の衝突・境界イベント](../developer/ParticleEvents.md#openreflectredistributed_reflectperiodic)を参照してください。

#### `[particles.reservoir]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `inflow_model` | string | `"source_vdf"` | `source_vdf` / `infinity_barrier` |
| `phi_infty_v` | float | `0.0` | barrierが参照する無限遠基準電位 [V] |
| `face_potential_grid_n` | int | `3` | 注入面平均電位の`N x N`評価格子。`>=1` |

`infinity_barrier` は `inflow` だけに適用します。一様電場と併用する場合は、`particles.reservoir.phi_infty_v` を有効な reservoir 基準に設定します。
補正式と適用範囲は[境界から粒子を流入させる](../guide/ReservoirInjection.md#3-流入の写像を選ぶ)を参照してください。

#### `[particles.secondary_emission]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `max_generations` | int | `3` | `>=0`、一次粒子は世代 0。上限でさらに放出せず、放出するはずだった個数を記録 |

### `[[particles.species]]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `species_key` | string | `"species_<1-based index>"` | 安定 ID。1–64 文字、粒子種間で一意 |
| `enabled` | bool | `true` | 種を有効化 |
| `charge_c` | float | `-1.602176634e-19` | 粒子電荷 [C]。enabled speciesでは非0 |
| `mass_kg` | float | `9.10938356e-31` | 粒子質量 [kg]。`>0` |

#### `[particles.species.distribution]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `drift_velocity_m_s` | float[3] | `[0,0,-8e5]` | ドリフト速度 [m/s] |
| `temperature_k` | float | `2.0e4` | 温度 [K]。`>=0` |
| `temperature_ev` | float | 未指定 | 温度 [eV]。`>=0`、`distribution.temperature_k` と排他 |
| `number_density_cm3` | float | 未指定 | Maxwell 分布の上流密度 [cm^-3]。`>0`、`distribution.number_density_m3` と排他 |
| `number_density_m3` | float | 未指定 | Maxwell 分布の上流密度 [m^-3]。`>0`、`distribution.number_density_cm3` と排他 |

境界流入では、密度 2 形式のちょうど一方と温度を指定します。流入速度は drift Maxwell 分布です。
分布の意味は[境界から粒子を流入させる](../guide/ReservoirInjection.md#2-分布を選ぶ)を参照してください。

#### `[particles.species.source]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `mode` | string | `[particles.species.source]` 指定時に必須 | `volume_seed` / `photo_raycast` / `secondary_emission` |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | 位置下限 [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | 位置上限 [m] |
| `inject_face` | string | 未指定 | `photo_raycast`の照射開口面 |
| `emit_current_density_a_m2` | float | `0.0` | レイ垂直面基準の放出電流密度 [A/m^2] |
| `deposit_opposite_charge_on_emit` | bool | `false` | 放出元要素に逆符号電荷を堆積 |
| `normal_drift_speed_m_s` | float | `0.0` | 放出法線方向ドリフト [m/s] |
| `ray_direction` | float[3] | 注入面内向き法線 | レイ方向 |
| `region_mode` | string | `"absolute"` | `absolute` / `face_fraction`。`face_fraction` は `photo_raycast` 専用 |
| `uv_low`, `uv_high` | float[2] | 未指定 | `face_fraction` で両方必須。各成分は `[0,1]`、`source.pos_low` / `source.pos_high` と排他 |

##### `mode = "photo_raycast"`

制約:

| 条件 | 内容 |
|---|---|
| 領域 | `[domain]`が必須 |
| 時間 | `run.batch.duration_s > 0` が必須 |
| 放出量 | `source.emit_current_density_a_m2 > 0`, `sampling.rays_per_batch > 0` が必須 |
| 注入面 | `source.inject_face` が必須 |
| 粒子属性 | `charge_c` は非ゼロ、`mass_kg > 0` |
| レイ方向 | 正規化可能で、注入面内向き法線との内積が正 |
| 使用不可 | `sampling.volume_macro_particles_per_batch`, `distribution.number_density_cm3` / `distribution.number_density_m3`, `sampling.weight`, `sampling.target_macro_particles_per_batch` |

レイ重み、周期 image、再吸収、closed PE の設定は
[光電子の放出とライフサイクル](../models/PhotoelectronEmission.md)を参照してください。

#### `[particles.species.sampling]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `weight` | float | `1.0` | マクロ粒子重み。`>0`。境界流入では `target_macro_particles_per_batch` といずれか一方を明示 |
| `target_macro_particles_per_batch` | int | `sampling.weight` といずれか一方必須 | マクロ粒子数 target。`>0`、または species 2 以降で species 1 の重みを共有する `-1` |
| `volume_macro_particles_per_batch` | int | `0` | この粒子種が 1 バッチに生成するマクロ粒子数。`>=0` |
| `rays_per_batch` | int | `0` | 1 バッチの発射レイ数 |

境界流入では `sampling.weight` と `sampling.target_macro_particles_per_batch` のちょうど一方を明示します。

制約:

| 条件 | 内容 |
|---|---|
| 粒子数 | 境界流入を使わない場合、有効な粒子種全体で `sampling.volume_macro_particles_per_batch` の合計が 1 以上必要 |
| 重み自動解決 | `inflow` を持たない `volume_seed` では `sampling.target_macro_particles_per_batch` は使用不可 |

境界流入だけを使う粒子種では `source` と `sampling.volume_macro_particles_per_batch` を省略し、体積生成は行いません。
Maxwell 分布で体積生成を加える場合は `source.mode="volume_seed"` と正の
`sampling.volume_macro_particles_per_batch` を指定します。

#### `[particles.species.charging]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `closure` | string | `"explicit"` | 表面 source 電荷 closure。`explicit` / `fixed_current` / `neutral_return` |
| `target_absorbed_current_a` | float | 未指定 | `fixed_current` の signed 吸収電流 [A]。0または`charge_c`と同符号 |
| `target_emission_current_a` | float | 未指定 | `fixed_current` の signed 放出反作用電流 [A]。0または`charge_c`と逆符号 |

| `charging.closure` | 効果 | 必要条件 |
|---|---|---|
| `explicit` | 追跡した電荷をそのまま反映 | 既定 |
| `neutral_return` | PE による表面総電荷増分を 0 に補正 | 負電荷 `photo_raycast`、放出反作用あり、注入面は反射、escape / `soft_discard` なし、未帰還率 $\le5\%$ |
| `fixed_current` | 要素別の標本分布を保ち、総電流を target へ倍率化 | 正の `run.batch.duration_s`、手動 `target_*_current_a` または `[sheath]`、非ゼロ target に対応する raw channel |

`fixed_current` と `neutral_return` は同じ species では排他です。
補正式、PE return の二重計上を避ける条件、標本数の収束確認は
[表面電荷更新の数値仕様](../developer/SurfaceChargeNumerics.md)、[光電子の放出とライフサイクル](../models/PhotoelectronEmission.md)、
[出力形式リファレンス](OutputReference.md#charge-ledger)を参照してください。

#### `[particles.species.boundary]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `x_low`, `x_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `y_low`, `y_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `z_low`, `z_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |

このtableは直前の`[[particles.species]]`に属します。

`inherit`は`[particles.boundary]`のglobal作用を使います。周期面はoverrideできません。
`source.inject_face`は`photo_raycast`の生成面、species境界は生成後の軌道作用です。
`inflow`は外部からの生成であり、外向き作用を上書きしません。

#### `[particles.species.inflow]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `x_low`, `x_high` | string | 省略 | `reservoir`。省略時は流入なし |
| `y_low`, `y_high` | string | 省略 | `reservoir`。省略時は流入なし |
| `z_low`, `z_high` | string | 省略 | `reservoir`。省略時は流入なし |

`reservoir` は選択した非周期 box 面全体から流入させます。対象面の有効な粒子境界は `open` でなければなりません。
境界流入だけなら `source` と `sampling.volume_macro_particles_per_batch` は省略します。
Maxwell 分布で体積生成を追加する場合は `source.mode="volume_seed"` と正の
`sampling.volume_macro_particles_per_batch` を指定します。
粒子源が `volume_seed` の場合だけ併用でき、複数面の
`sampling.target_macro_particles_per_batch` は全流入面の合計です。
流束と電位障壁の詳細は[境界から粒子を流入させる](../guide/ReservoirInjection.md)を参照してください。

#### `[particles.species.surface_interaction]`

##### `[particles.species.surface_interaction.secondary_emission]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `enabled` | bool | `false` | 吸収後に放出する |
| `model` | string | `"whipple"` | `whipple` / `constant` |
| `yield_max` | float | `1.5` | 無次元。有効な Whipple で `>0` |
| `energy_max_ev` | float | `500.0` | [eV]。有効な Whipple で `>0` |
| `const_yield` | float | `0.0` | 無次元。有効な constant で `>=0` |
| `target_species` | string | 未指定 | 有効時必須、有効な `secondary_emission` 種を参照 |

`source.mode="secondary_emission"` では `sampling.*`、`inflow.*`、`source.emit_current_density_a_m2`、光線数、`charging.closure` を省略します。
電荷は負、質量は正、温度は 0 以上とし、重みは入射粒子から受け継ぎます。
外部シースの役割へ指定する場合は電荷 -e、ambient electron と同じ質量、open な z-high 面も要求します。
モデルと Phase 1 の暫定電流配分は [二次電子放出](../models/SecondaryElectronEmission.md)を参照してください。

### `[fields]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `boundary` | string | `"free"` | `free` / `periodic2` |

`periodic2`では`[domain]`が必要です。非零・零モードoperatorを明示する場合は`fields.periodic.backend` も指定します。
solverとの組合せは
[場の評価の互換表](../guide/FieldSolvers.md#solverと場境界の互換表)に従います。

#### `[fields.solver]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `method` | string | `"auto"` | `direct` / `treecode` / `fmm` / `auto` |
| `normalization` | string | `"si"` | `si` / `box` / `mesh` / `length` |
| `length_scale_m` | float | `1.0` | `fields.solver.normalization="length"` の長さ [m]。`>0` |
| `auto_fmm_min_elements` | int | `256` | `auto` が direct から FMM へ切り替える要素数。`>=1` |
| `near_field_precision` | string | `"double"` | FMM 近接場の panel 和の精度。`double` は従来と bit 単位で同じ。`single` は葉ごとの panel 表を単精度で持ち 8 レーン SIMD で評価する（その和が約 1.5 倍速く、場の相対誤差は 10⁻⁶ 程度。展開次数 4 の FMM 打ち切り誤差は 2×10⁻⁴）。出力 `summary.txt` の `field_near_field_precision` に記録される |

`fields.solver.method` は、要素電荷から Coulomb 場を評価する方式を選びます。

| `fields.solver.method` | 用途 | 対応する場境界 |
|---|---|---|
| `direct` | 要素数が小さい場合の厳密な全対全評価、split reference | `free`、または条件を満たす`periodic2` split reference |
| `treecode` | 中規模以上の近似評価 | `fields.boundary="free"` |
| `fmm` | 大規模評価、`periodic2`、FMM コア検証 | `fields.boundary="free"` / `"periodic2"` |
| `auto` | 要素数に応じて direct / FMM を自動選択 | `fields.boundary="free"` |

solverと場境界の組合せは[場の評価の互換表](../guide/FieldSolvers.md#solverと場境界の互換表)で確認できます。
要素電荷は常に三角形上の一定面密度（P0 panel）として評価され、設定でkernelを選択しません。

`fields.solver.normalization` は内部座標のみを変え、電場・電位の出力単位は SI のままです。

| `fields.solver.normalization` | 長さスケール |
|---|---|
| `si` | 入力 SI 座標をそのまま使う |
| `box` | `domain.box_max - domain.box_min`の最大幅。`[domain]`が必須 |
| `mesh` | mesh bounding box の最大幅。mesh が空なら `fields.solver.length_scale_m` |
| `length` | `fields.solver.length_scale_m` |

mode ごとの追加キー:

| `fields.solver.method` | 評価 | 使用する追加キー |
|---|---|---|
| `direct` | 全 source を直接加算する `O(MN)` 参照計算 | なし |
| `treecode` | 遠方を monopole、近傍を panel kernel で評価 | `fields.solver.tree.theta`、`fields.solver.tree.leaf_max` |
| `fmm` | Coulomb FMM で大規模近似評価 | `fields.solver.tree.theta`、`fields.solver.tree.leaf_max`、periodic2 関連キー |
| `auto` | `fields.solver.auto_fmm_min_elements` 未満は direct、以上は FMM | `fields.solver.auto_fmm_min_elements`、FMM 選択時の `fields.solver.tree.theta`、`fields.solver.tree.leaf_max` |

FMM の選択と精度確認は [FMM を使う](../guide/FMM.md)、実装は
[Coulomb FMM コア詳細](../developer/FMMCore.md)を参照してください。

`fields.solver.tree.theta` と `fields.solver.tree.leaf_max` は、明示指定がなければ要素数から次の値を使います。

| 要素数 `nelem` | `fields.solver.tree.theta` | `fields.solver.tree.leaf_max` |
|---:|---:|---:|
| `< 1500` | `0.40` | `12` |
| `1500 <= nelem < 10000` | `0.50` | `16` |
| `10000 <= nelem < 50000` | `0.58` | `20` |
| `50000 <= nelem` | `0.65` | `24` |

##### `[fields.solver.tree]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `theta` | float | 未指定時は上記の要素数 heuristic | treecode / FMM の MAC パラメータ。`0<theta<=1`、大きいほど速く粗い |
| `leaf_max` | int | 未指定時は上記の要素数 heuristic | treecode / FMM の葉あたり最大 source 数。`>=1` |

#### `[fields.periodic]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `image_layers` | int | `1` | `periodic2` の近傍 image shell。`>=0`、`cached_kneq0` では `>=1` |
| `backend` | string | `"finite_images"` | `finite_images` / `cached_kneq0` / `panel_spectral_reference` |
| `ewald_alpha` | float | `0.0` | cache 生成の Ewald 分割値。`>=0`、0 は自動 |
| `ewald_layers` | int | `4` | cache 生成の実・逆空間 shell 深さ。far correction 有効時 `>=1` |
| `cache_dir` | string | 環境変数から解決 | 演算子の保存先。明示指定を優先し、省略時は `$BEACH_CACHE_DIR/periodic2` → `$XDG_CACHE_HOME/beach/periodic2` → `$HOME/.cache/beach/periodic2`。空の環境変数は未設定と扱う。`cached_kneq0` で明示する場合は空文字不可 |
| `generation_tolerance` | float | `1.0e-8` | cache identity に含める生成許容誤差。`cached_kneq0` では有限かつ `>0` |
| `lower_boundary_model` | string | 必須 | `symmetric_vacuum` / `e_bottom_zero` / `e_top_from_sheath`。最後は `sheath.closure="zero_current"` 専用 |

`fields.boundary="periodic2"` では `domain.periodic_axes=["x","y"]` が必要です。通常は `fields.solver.method="fmm"` を使い、小規模な split reference に限り `"direct"` を使います。

##### `[fields.periodic.upper_fourier]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `switch_height` | string / float | `"off"` | `"off"` / `"auto"` / 高さ [m]。有効時は periodic2 + cached 非零モード + 零モード除外が必要。$z_\mathrm{top}+0.1\min(L_x,L_y)\le z<H$ |
| `mode_layers` | int | `0` | 0 は減衰から自動、明示値は1〜32。自動次数は reference の次数以上、最大32 |

auto の候補面と代理電荷、起動・再開時の選択は[periodic2静電場](../models/PeriodicElectrostatics.md#上部真空の-fourier-評価への切替)を参照してください。
外部シース接続で switch が off の場合は警告し、計算を続けます。

##### `[fields.periodic.reference]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `mode_layers` | int | `4` | Fourier mode cutoff。`>=1` |
| `panel_quadrature_order` | int | `12` | panel 面積積分次数。`>=2` |

#### periodic2 の組合せ制約

periodic2では`[domain]`、`periodic_axes=["x","y"]`、`fields.boundary="periodic2"`が必要です。
`examples/closed_photoelectron.toml`は、x/y周期、境界reservoir、閉じた光電子を組み合わせた基準例です。
同じ周期条件をfield、collision、`photo_raycast`に適用します。

| `fields.periodic.backend` | 意味 |
| --- | --- |
| `panel_spectral_reference` | 小規模なsplit reference |
| `cached_kneq0` | versioned operatorを再利用するproduction非零モード |

非零・零modeの分割、`exclude_k0`の役割、下側境界から決まる平均場は
[periodic2静電場](../models/PeriodicElectrostatics.md)を参照してください。

#### `[fields.external]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `electric_v_m` | float[3] | `[0,0,0]` | 一様外部電場 [V/m] |
| `electric_magnitude_v_m` | float | 未指定 | 一様外部電場の大きさ [V/m]。`>=0` |
| `electric_azimuth_deg` | float | `0.0` | `fields.external.electric_magnitude_v_m` 指定時の xy 面内方位角 [deg]。明示指定には `fields.external.electric_magnitude_v_m` が必須 |
| `electric_elevation_deg` | float | `0.0` | `fields.external.electric_magnitude_v_m` 指定時の xy 面からの仰角 [deg]。明示指定には `fields.external.electric_magnitude_v_m` が必須 |
| `magnetic_t` | float[3] | `[0,0,0]` | 一様磁場 [T] |

`fields.external.electric_v_m` での直接指定と、`fields.external.electric_magnitude_v_m` と角度による指定は排他です。

#### 要素 source の固定規則

| 項目 | 規則 |
| --- | --- |
| source | 各 `q_elem` を三角形上の一様な面電荷密度として扱う |
| solver | `direct` / `treecode` / `fmm` / `auto`。`auto` は `fields.solver.auto_fmm_min_elements` で direct / FMM を選ぶ |
| 対象表面 | `insulator` / `conductor` の共通 source 離散化 |
| Treecode | 厳密 panel near + monopole far |
| FMM | 厳密 panel near + 厳密 triangle P2M |
| 面の向き | 各要素の `mesh.parts[i].surface_side` を指定 |

`outward_closed` は、法線の向きが整合した閉じた two-manifold だけに使用できます。

### `[sheath]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `closure` | string | `"none"` | `none` / `zero_current` |
| `reference_area_m2` | float | domainのx-y面積 | 電流密度を総電流へ変換する面積 [m^2]。`>0` |

省略時は species ごとの手動電流 target を使います。外部シースを使う場合は次の model を選びます。

| `sheath.closure` | 効果 | 詳細 |
|---|---|---|
| `none` | 外部シース closure を使わない | species 側の `charging.closure` を設定 |
| `zero_current` | 外部 1-D シースの Zhao 零電流根から species 別電流と z-high のエネルギー障壁を決定 | [Zhao closure](../models/ZhaoStationaryClosure.md) |

入力制約:

| 項目 | 必要条件 |
|---|---|
| role species | enabled、相互に異なる。電子・イオン・PE は `charging.closure="fixed_current"`、SEE は `charging.closure` を省略 |
| ambient electron / ion | z-high reservoir から流入。electron の `distribution.drift_velocity_m_s` の z 成分は 0 以下、ion は負。electron にもイオンと同じ drift を与え、`sheath.zhao.upstream_field_floor_v_m=2.0e-3` とする（drift 0 なら不要。[理由](../models/ZhaoUpstreamBand.md)）。手動 `target_*_current_a` は指定不可 |
| PE なし | `sheath.photoelectrons.source_scale=0.0`、PE 固有キーは省略、`sheath.zhao.branch="auto"` または `"c"` |
| PE あり | 負電荷の `photo_raycast`、`source.inject_face="z_high"`、`source.deposit_opposite_charge_on_emit=true`、有効な z-high 境界は `open`。`source.emit_current_density_a_m2`（レイ垂直面基準）× ray 仰角の sin が外部源 $e\,s_{UV}\,n_{pe,ref}\sin\alpha\,v_{th}/(2\sqrt\pi)$ から 5% 以上、または ray の仰角が $\alpha$ から 1 度以上ずれると起動時に警告（停止しない） |
| species 物性 | 単価電荷、ambient electron と PE の質量は同一、$T_e>0$、$T_{pe}>0$、$T_i\ge0$ かつイオンの入口が超音速 $m_i u_z^2/e >$ `sheath.zhao.ion_pressure_factor` $\cdot T_i$ |
| 外部場 | `fields.external.magnetic_t=[0,0,0]`、`particles.reservoir.inflow_model="infinity_barrier"` は使用不可 |
| 場 | `fields.boundary="periodic2"`、x/y 周期、明示的な split `[fields.periodic]`（`fields.periodic.backend` が `cached_kneq0` または `panel_spectral_reference`）。上端面の平均電位を壁電位 $\phi_0$ に固定するため |
| 外部根の更新 | `outflow_refresh_batches>0` では PE、SEE、反射電子（既定）のいずれかが外部源になる。`include_reflected_electrons=false` なら PE または SEE が必要 |

PEなしType Cは$J_e+J_i=0$を満たすelectron/ion吸収targetとz-high kinetic barrier mapだけを生成し、
PE emission / return / escape targetは生成しません。SEE があれば PE なしでも更新でき、初期根だけは源 0 の Type C です。
`sheath.species.ion` が選ぶ粒子種の `distribution.number_density_m3` は無限遠の太陽風密度です。`sheath.species.electron` にも同じ値を書きます（違えば設定エラー）。電子の流入は、Zhao 根が無限遠の準中性から解いた上流電子 Maxwellian 密度で注入します（[reservoir の密度](../models/ZhaoStationaryClosure.md#上流電子密度で電子を入れる)）。PE 放出電流密度の入力は粒子分布の標本化に使い、電流 target は closure が決めます。

split `[fields.periodic]` では z-high 面の水平平均電位を外部シースの壁電位 $\phi_0$ に固定し、x/y 周期では外部障壁で
戻る粒子をセル内の一様な位置へ戻します。`outflow_refresh_batches>0` は外部根を観測 PE + SEE 流出から準定常に解き直す
弱連成です。box 外の過渡シースは解きません。電流・障壁・PE return・更新の定義と出力は
[Zhao closure](../models/ZhaoStationaryClosure.md)を参照してください。`sheath.closure="none"` では `sheath.closure` 以外を指定しません。

計算例は `examples/zero_current.toml`、外部根を更新する例は
`examples/zero_current_refresh.toml` です。

#### `[sheath.zhao]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `branch` | string | `"auto"` | `auto` / `a` / `b` / `c`。PE なしは `auto` / `c` のみ |
| `upstream_band_tolerance` | float | `0.0` | 電子の drift がある A/C の根で、上流に接する $E^2<0$ の幅を A は $\lvert\phi_m\rvert$、C は $\lvert\phi_0\rvert$ の何割まで許すか。`0<=x<1`、0 はこの判定を使わない。`sheath.zhao.upstream_field_floor_v_m` と両方 0 なら drift のある A/C を棄却。設定例は幅ではなく場の下限で受理する。[電子の drift と上流の帯](../models/ZhaoUpstreamBand.md) |
| `upstream_field_floor_v_m` | float | `0.0` | 上流境界を「場がこの外部レベル [V/m] まで落ちた点」に置く。電子の drift がある A/C の根で、上流の帯の内側の場の不足 $\sqrt{-\min E^2}$ がこれ以下なら受理する。太陽風の運動電場 $\lvert\mathbf{u}\times\mathbf{B}\rvert$（400 km/s、5 nT で約 2 mV/m）が目安。`>=0`、0 はこの判定を使わない。`sheath.zhao.upstream_band_tolerance` と両方 0 なら drift のある A/C を棄却。[電子の drift と上流の帯](../models/ZhaoUpstreamBand.md) |
| `ion_pressure_factor` | float | `1.0` | 外部根のイオン流体（Sana & Mishra 2026）の圧力係数。1 は等温、3 は 1 次元断熱。イオン温度は `sheath.species.ion` の温度。`>0`。太陽風（入口エネルギー約 1 keV）では根への影響は 0.1 mV 以下。[外部シースとの接続](../models/ZhaoStationaryClosure.md#仮定) |

#### `[sheath.species]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `electron` | string | 未指定 | ambient electron の `species_key`。必須、1–64 文字 |
| `ion` | string | 未指定 | 太陽風イオンの `species_key`。必須、1–64 文字。温度は外部根のイオン流体にも使う |
| `photoelectron` | string | PE有効時に必須 | PE の `species_key`。1–64 文字 |
| `secondary_electron` | string | 未指定 | 有効な SEE 種の識別子。省略可、1–64 文字、完全に追跡。省略時、`secondary_emission` 種が 1 つだけならその種を当てる（複数なら設定エラー） |

#### `[sheath.photoelectrons]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `solar_elevation_deg` | float | PE有効時に必須 | 表面放出に使う太陽高度角 $\alpha$。$0<\alpha\le90$ degree |
| `ref_density_m3` | float | PE有効時に必須 | PE基準密度 $n_{pe,ref}$ [m^-3]。`>0` |
| `source_scale` | float | `1.0` | $s_{UV}$。`>=0`、0 は PE なし |

#### `[sheath.coupling]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `outflow_refresh_batches` | int | `0` | 外部根を、accepted batch N 個ごとに z-high で観測した PE + SEE 流出スペクトルから解き直す。`>=0`、0 は初期根を固定 |
| `include_reflected_electrons` | bool | `true` | 外部根を更新するとき（`outflow_refresh_batches>0`）だけ効く。外部根の源に、セル内で反射して z-high を外向きに抜けた上流電子も含める（外部モデルから見れば由来に関わらず境界からの放出電子）。明示的に `true` と書くなら `outflow_refresh_batches>0` が必要。各バッチで追跡したその脱出を電子の吸収目標に足し戻し、PE（なければ SEE）の脱出・帰還目標から差し引く。PE も SEE もなければ反射電子だけが外部源で、電子の吸収を根の流入 − 脱出、脱出を根の脱出に固定する（この場合は放出電子の種なしで `outflow_refresh_batches>0` を使える）。窓の反射束は `electron_outward_flux_m2_s` に出る |

### `[output]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `enabled` | bool | `true` | ファイル出力の有効/無効 |
| `dir` | string | `"outputs/latest"` | 出力先ディレクトリ |

生成条件、列定義、電位の評価規約、外部シース状態、ledger の解釈は
[出力形式リファレンス](OutputReference.md)を参照してください。

#### `[output.diagnostics]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `charge_rel_change_threshold` | float | `1.0e-8` | 相対変化量の監視値。有限かつ `>=0`、停止条件には未使用 |
| `charge_floor_c` | float | `1.0e-30` | `rel_change` 計算時の分母下限。有限かつ `>0` |

以下は電荷変化の診断用です。計算終了は `run.batch_count` で決まり、これらを変えても早期停止しません。

#### `[output.final]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `mesh_potential` | bool | `false` | `mesh_potential.csv` を出力 |

#### `[output.history]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `potential` | bool | `false` | `potential_history.csv`を出力。`[domain]`があれば同じbatchの`top_reference_history.csv`も出力 |
| `stride_batches` | int | `1` | 履歴 CSV の出力間隔 [batch]。`>=0`、0 は無効 |

#### `[output.checkpoint]`

| キー | 型 | 既定値 | 説明・条件 |
|---|---|---:|---|
| `stride_batches` | int | `0` | 再開用 checkpoint の間隔 [accepted batch]。`>=0`、0 は定期出力なし、正値では `output.enabled=true` |

## 座標・配置の補助パラメータ

`domain.box_origin` / `domain.box_size`、`particles.species.source.region_mode`、`[mesh.groups.<name>]`、`[[mesh.parts]]` の配置キーは box 基準の値を実座標・実寸へ変換します。各 table の型・既定値は上の表を参照してください。

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
- `placement_mode="absolute"`、`size_mode="absolute"`、`source.region_mode="absolute"` は直接指定値をそのまま使います。

入力前に[`beach.toml`を作成・検証する](../guide/Configuration.md)の `beachx lint` で組合せを確認してください。

---

## 検証ルール

| 項目 | ルール |
|---|---|
| 未知キー | すべてエラー |
| `[particles]` | 粒子種と `particles.secondary_emission.max_generations` を保持。SEE の独立注入は不可 |
| 型 | schema と Fortran パーサの両方で検証 |
| 列挙値・値域・無効な機能への指定 | `beachx lint` が詳細を診断。Fortran は実行に必要な最小条件を確認 |

実行前には次を推奨します。

```bash
beachx lint beach.toml
```
