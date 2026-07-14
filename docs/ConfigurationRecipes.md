title: 設定レシピ

Lang: [日本語](ConfigurationRecipes.md) | [English](ConfigurationRecipes.en.md)

# 設定レシピ

典型的なケースは、公式入門ケースの`beach.toml`を基に、必要なtableやsectionを置き換えて作ります。
全キーの定義は [入力パラメータリファレンス](Parameters.html)、高水準記法の詳細は [beachx config / 高水準記法ガイド](Configuration.html) にまとめています。

以下は[公式入門ケース](Tutorial.html)を基準にした差分例です。
断片だけでは実行できない場合があります。

## 公式の実行手順

```bash
beachx config init
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

`beach.toml` はそのまま Fortran 実行系に渡せます。高水準記法も Fortran parser が読み込み時に正規化します。

## レシピ一覧

| レシピ | 変更する主な場所 | 使う場面 |
| --- | --- | --- |
| 平面メッシュで最小実行 | `[mesh]`, `[[mesh.templates]]` | まず動作確認する |
| 粒子注入を選ぶ | `[[particles.species]]` | 流入、初期粒子、光電子放出を切り替える |
| 2軸周期境界 | `[sim]` | 周期セルを使う |
| 無限周期 + 外部kinetic sheath | `[sim]`, `[periodic2]`, `[outer_plasma]`, `[coupling]` | 月面レゴリスのproduction設定 |
| UV光電子を外部sheathへ含める | `[outer_plasma]`, `[[particles.species]]` | 放出・帰還とouter空間電荷を整合させる |
| OBJ メッシュ | `[mesh]` | 外部形状を使う |
| 履歴出力 | `[output]` | 時間発展を可視化する |
| 再開実行 | `[output]` | checkpoint から続ける |

## 平面メッシュの解像度を変える

```toml
[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_model = "insulator"
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]
```

`nx` / `ny` を増やすと要素数が増えます。初回は小さく動かし、出力と粒子数を確認してから解像度を上げてください。

## 粒子注入を選ぶ

| `source_mode` | 使う場面 | 主なキー |
| --- | --- | --- |
| `volume_seed` | 箱内に初期粒子を置く小テスト | `npcls_per_step`, `pos_low`, `pos_high` |
| `reservoir_face` | 面から Maxwellian などの流入を与える通常ケース | `number_density_cm3`, `temperature_ev`, `inject_face`, `target_macro_particles_per_batch` |
| `photo_raycast` | 表面からの光電子放出を raycast で扱う | `rays_per_batch`, `emit_current_density_a_m2`, `ray_direction` |

`reservoir_face` の最小例:

```toml
[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
inject_face = "z_high"
pos_low = [0.0, 0.0, 4.0]
pos_high = [1.0, 1.0, 4.0]
drift_velocity = [0.0, 0.0, -4.0e5]
```

`target_macro_particles_per_batch` は 1 batch あたりの計算粒子数を固定し、粒子重みを流入量から解きます。
物理粒子数を重みで直接指定したい場合は `w_particle` を使います。

## 2軸周期境界

`periodic2`は、3軸のうち2軸を周期軸とする場の境界条件です。productionの無限周期operatorには、
`field_solver = "fmm"`と`cached_kneq0`を使います。小規模なcorrectness fixtureに限り、
`field_solver = "direct"`と`panel_spectral_reference`も選択できます。

```toml
[sim]
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]

bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"

field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "cached_kneq0"
field_periodic_cache_dir = ".beach_cache/periodic2"

[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "symmetric_vacuum"
```

要件:

- `sim.use_box = true`
- ちょうど 2 軸が periodic
- 各周期軸の `box_max - box_min > 0`
- production cacheでは `sim.field_solver = "fmm"`

cache miss時に無限周期operatorを生成し、次回以降はfingerprintを検証して再利用します。
物理設定とmesh設定が同じcaseでは、cache directoryを共用できます。productionで
`field_periodic_far_correction = "none"`を使うのは、有限個の周期像だけを使う旧設定を再現する場合に限ります。

`symmetric_vacuum` は追加パラメータを必要とせず、非中性セルの面平均電場を上下へ等分します。
総表面電荷を `Q`、周期面積を `A` とすると、上下遠方場はそれぞれ
`E_z = +/- Q/(2 epsilon0 A)` です。旧計算の再現には `e_bottom_zero` を指定できますが、これは
下側電場を0に固定し、上側場を `Q/(epsilon0 A)` とする強い遮蔽近似です。

## 無限周期 + 外部kinetic sheath

月面レゴリスのproduction設定では、表面電荷の`k != 0`成分を無限周期cacheで解きます。
一方、面平均の`k = 0`成分は外部の1D kinetic sheathに接続します。次は、設定済みのbox、mesh、
ambient electron/ion speciesに追加する主要なsectionです。

```toml
[sim]
b0 = [0.0, 0.0, 0.0]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_far_correction = "cached_kneq0"
field_periodic_cache_dir = ".beach_cache/periodic2"
reservoir_potential_model = "none"
sheath_injection_model = "none"

[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "symmetric_vacuum"

[outer_plasma]
model = "kinetic_1d"
photoelectron_closure = "none"
return_model = "kinetic_1d_profile_return"
interface_z = 9.899494936611664e-4
infinity_potential = 0.0
debye_length = 10.5132
thermal_voltage = 10.0

[coupling]
update_mode = "explicit"
particle_transfer_mode = "electrostatic_1d_instant_return"
outer_update_stride = 1
field_evolution_timescale = 6.060915267313266e-8
max_frozen_field_ratio = 0.1
outer_queue_enabled = false
```

`interface_z`はz-high box面と一致させます。負電荷と正電荷の`reservoir_face` speciesを、それぞれ無限遠の
ambient electron VDFとion VDFとして使います。両speciesをz-highに定義し、準中性条件とion Bohm流入条件を
満たすように設定してください。

`infinity_potential = 0`は、無限遠電位のゲージを固定する設定です。流入粒子の加減速と、外向き粒子の
escape/returnは、共通のkinetic profileから得られる`phi_interface - phi_infinity`を使って計算します。

`debye_length`、`thermal_voltage`、`field_evolution_timescale`は、例の数値をそのままコピーせず、
対象とするplasmaと時間scaleから決めてください。まず`outer_update_stride = 1`で検証します。実運用で
`100`などに増やす場合は、表面電位、吸収・escape flux、電荷収支、離脱力に主要な変化がないことを確認してください。

非単調分枝、sub-Bohm ion、frozen-field上限超過を検出した場合は、別のmodelにfallbackせず停止します。
完全な小規模例は`examples/periodic2_kinetic_outer.toml`にあります。全パラメータの定義は
[入力パラメータリファレンス](Parameters.html)から確認できます。

## UV光電子を外部sheathへ含める

UVを有効にする場合、外部空間の平均光電子密度は`kinetic_mean`で解きます。
表面電荷の更新には、明示的に追跡する`photo_raycast`粒子を使います。上の設定に次の変更を加えてください。

```toml
[outer_plasma]
model = "kinetic_1d"
photoelectron_closure = "kinetic_mean"
return_model = "kinetic_1d_profile_return"
interface_z = 9.899494936611664e-4
infinity_potential = 0.0
debye_length = 10.5132
thermal_voltage = 10.0

[[particles.species]]
enabled = true
source_mode = "photo_raycast"
emit_current_density_a_m2 = 4.5e-6
rays_per_batch = 5000
deposit_opposite_charge_on_emit = true
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
temperature_ev = 2.2
normal_drift_speed = 0.0
inject_face = "z_high"
ray_direction = [0.0, 0.0, -1.0]
```

先頭の負電荷`photo_raycast` speciesに設定した放出電流と温度から、mean closureが決まります。
`deposit_opposite_charge_on_emit = true`は必須です。legacy `photo_escape_model`とは併用できません。
`kinetic_mean`はouter profileだけを供給し、帰還電流を表面電荷に二重加算しません。

まず、mesh、batch duration、ambient流入を同じにしたUVなしとUVありのcaseを比較します。
`outer_plasma_profile.csv`、`summary.txt`のsolver residual、species別電流、電荷収支を確認してください。

`sim.sheath_injection_model = "zhao_*"`は流入分布だけを補正する旧モデルで、ここで使う外部
`kinetic_1d` Poisson profileとは別物です。両者、および`reservoir_potential_model`との同時利用は拒否されるため、
このレシピではどちらも`"none"`にします。

## OBJ メッシュを使う

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj"
surface_model = "insulator"
```

OBJ では全要素に同じ `surface_model` が割り当てられます。
template mesh のように object ごとに surface model を分けたい場合は、複数 template を使うか、現行の制約を確認してください。

## 履歴を出す

```toml
[output]
dir = "outputs/latest"
history_stride = 1
write_mesh_potential = true
write_potential_history = true
```

`write_potential_history` は履歴ごとに電位を再評価するため、要素数が多いケースでは重くなります。
まず `charge_history.csv` だけで収束傾向を確認し、必要なときに電位履歴を有効にしてください。

## 再開する

同じ出力ディレクトリから続ける場合:

```toml
[output]
dir = "outputs/latest"
resume = true
```

checkpoint 読み込み元と新しい出力先を分ける場合:

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../previous/outputs/latest"
```

`sim.batch_count` は累積到達 batch 数です。checkpoint が `batches=100` で、新しい設定が `batch_count=150` なら追加で 50 batch だけ進みます。
