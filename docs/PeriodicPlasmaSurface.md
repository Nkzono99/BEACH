title: プラズマ中の周期表面を設定する

Lang: [日本語](PeriodicPlasmaSurface.md) | [English](PeriodicPlasmaSurface.en.md)

# プラズマ中の周期表面を設定する

月面のレゴリスのように、太陽風と紫外線にさらされる表面の一部を、x/y 周期のセル 1 つで表すケースの作り方です。
セルの高さは Debye 長よりはるかに小さいので、上端面（z-high 面）の外にあるプラズマとシースをどう扱うかを
決める必要があります。このページでは上端の閉じ方を 4 つの方式から選び、方式ごとの設定の差分と、
実行後に確認する出力を示します。

外部シースと整合した帯電を求めるなら、外部シースとの接続（Zhao 定常シース）を使います。
表面内での電荷の再分配だけを比べるなら、閉じた光電子が最も単純です。

## 共通の構成

どの方式も、次の構成から始めます。完全な設定は [`examples/grouped/closed_photoelectron.toml`](../examples/grouped/closed_photoelectron.toml) です。

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0e-4, 1.0e-4, 1.0e-3]
periodic_axes = ["x", "y"]

[fields]
boundary = "periodic2"

[fields.solver]
method = "fmm"

[particles.boundary]
z_low = "open"
z_high = "open"

[[particles.species]]          # 太陽風電子（イオンも同じ形）
species_key = "solar_wind_electron"

[particles.species.inflow]
z_high = "reservoir"

[[particles.species]]          # 光電子
species_key = "photoelectron"

[particles.species.source]
mode = "photo_raycast"
inject_face = "z_high"
deposit_opposite_charge_on_emit = true
```

- **太陽風:** 電子とイオンを、上端面から境界流入で入れます（[境界から粒子を流入させる](ReservoirInjection.html)）。
- **光電子:** 上端面から照射して、当たった表面から放出します（[光電子の放出と電荷の閉じ方](PhotoelectronEmission.html)）。
- **場:** x/y 周期の場を FMM で計算します（[periodic2 の静電場](PeriodicElectrostatics.html)）。

## 上端の閉じ方を選ぶ

| 方式 | 上端面で起きること | 表面の総電荷 | 電位の基準 | 使う場面 |
|---|---|---|---|---|
| 1. 開放 | 上端に達した粒子はすべて脱出する | 追跡のまま | なし（上端面の平均との差で読む） | 比較の基準 |
| 2. 閉じた光電子 | 光電子だけ上端で反射し、光電子の正味の電流を 0 に合わせる | 太陽風の分だけ変わる | なし（同上） | 表面内の再分配を調べる |
| 3. スカラー障壁 | 上流電位と上端面の電位の差で、流入する粒子と出ていく粒子をふるい分ける | 追跡のまま | 上流電位 `phi_infty_v` | 単一の障壁による比較 |
| 4. 外部シース | 外部シースの零電流根から、電流・障壁・電位基準を与える | 浮遊条件で 0 に拘束 | 上流のプラズマを 0 V | 外部シースと整合した帯電 |

方式は 1 つだけ選びます。閉じた光電子とスカラー障壁を同じ計算で重ねず、外部シースでは
`particles.reservoir.inflow_model="infinity_barrier"` を使えません。

### 1. 開放

共通の構成のままです。上端面の開放面に達した粒子は脱出します（`particles.boundary.open_model="escape"`、既定値）。
セルの高さは Debye 長より小さいので、上端に達した光電子の多くは、実際には外部のシースで表面へ押し戻されます。
開放は、他の方式と比べる基準としてだけ使います。

### 2. 閉じた光電子

光電子の粒子種で、上端面を反射にし、電荷の閉じ方を `neutral_return` にします。

```toml
[[particles.species]]
species_key = "photoelectron"

[particles.species.charging]
closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"          # 戻る位置を面内で一様にするなら "redistributed_reflect"
```

光電子の正味の電流が 0 になるよう、帰還した光電子の電荷を同じバッチの放出量に合わせます。太陽風による
総電荷の変化は拘束しません。仕組みと停止条件は
[光電子の放出と電荷の閉じ方](PhotoelectronEmission.html#閉じた光電子neutral_return)にあります。

### 3. スカラー障壁

外部プラズマの上流電位を決め、流入と流出の両方を、その電位との差で判定します。

```toml
[particles.reservoir]
inflow_model = "infinity_barrier"
phi_infty_v = 0.0
face_potential_grid_n = 5

[particles.boundary]
open_model = "potential_barrier"
```

流入は上端面の平均電位で、流出は粒子が横切った点の電位で判定します
（[流入の写像](ReservoirInjection.html#3-流入の写像を選ぶ)、[電位障壁](ParticleEscapeReturn.html#電位障壁potential_barrier)）。
上端面の電位は BEACH 内の電荷だけで決まり、外部シースの電位降下を含みません。

### 4. 外部シース（Zhao 定常シース）

外部シースの零電流根を解き、太陽風と光電子の電流を根の値に固定します。役割を持つ 3 つの粒子種は
`charging.closure="fixed_current"` にし、場は面平均成分を分けて扱う `cached_kneq0` にします。

```toml
[sheath]
closure = "zero_current"

[sheath.species]
electron = "solar_wind_electron"
ion = "solar_wind_ion"
photoelectron = "photoelectron"

[sheath.photoelectrons]
solar_elevation_deg = 60.0
ref_density_m3 = 6.4e7

[sheath.coupling]
outflow_refresh_batches = 50   # 外部根を更新しないなら省略

[fields.periodic]
backend = "cached_kneq0"
lower_boundary_model = "e_bottom_zero"
```

太陽風電子の drift は 0、イオンの drift は太陽風の法線速度（内向き、負の z）にします。電子とイオンの数密度には、
同じ太陽風密度を書きます。完全な設定は [`examples/grouped/zero_current.toml`](../examples/grouped/zero_current.toml)、
外部根を更新する例は [`zero_current_refresh.toml`](../examples/grouped/zero_current_refresh.toml)、
光電子なし（Type C）の例は [`zero_current_no_photo.toml`](../examples/grouped/zero_current_no_photo.toml) です。
モデルの仮定と既知の制約は[外部シースとの接続](ZhaoStationaryClosure.html)にあります。

## 結果を確認する

| 方式 | 見る出力 | 判断 |
|---|---|---|
| すべて | `charge_ledger.csv` の `charge_ledger_residual_C`、未解決粒子の電荷 | 残差が丸め誤差程度で、未解決粒子の電荷が結論に効かない |
| すべて | `top_reference_history.csv` の `potential_mean_V` と `potential_std_V` | 方式 1〜3 では各要素の電位からこの平均を引いて読む。ばらつきが大きいと、上端面を 1 つの面として扱う近似が弱い |
| 2 | `charge_ledger.csv` の `neutral_return_weight_scale`、`neutral_return_unresolved_fraction` | 倍率が 1 に近く、未帰還の割合が小さい（5% を超えると停止） |
| 3 | `summary.txt` の `reservoir_inflow_map`、`particle_ordinary_open_model` | 選んだ写像と判定が使われている |
| 4 | `summary.txt` の `surface_current_model_zhao_branch`、`surface_current_model_phi0_V` | 想定した branch と壁電位になっている |
| 4 | `fixed_current_history.csv` の `target_over_tracked` | 各粒子種の倍率が 1 の近くにある |
| 4（外部根の更新） | `matching_plane_history.csv` の `phi_H_V` と光電子の流出 | 壁電位と流出が落ち着いている |

`top_reference_history.csv` は `output.history.potential=true` と `output.history.stride_batches>0` で書かれます。
列の定義は[出力形式](OutputReference.html#履歴)にあります。

## 数値条件への依存を確かめる

方式によらず、少なくとも次を変えて、注目する量が変わらないことを確かめます。

1. 周期の場: 有限画像なら `fields.periodic.image_layers` を $N, N+1, N+2$ と増やす。`cached_kneq0` なら
   [遠方補正の精度](PeriodicFarCorrection.html)を確認する。
2. セルの高さ: 上端面を上下に動かす。
3. バッチ幅: $T, T/2, T/4$ と変え、同じ物理時刻で比べる（[バッチ幅を決める](BatchDurationStability.html)）。
4. 粒子の追跡: `particles.tracking.dt_s` を半分にし、`max_steps_per_particle` を 2 倍にする。
5. 標本数: マクロ粒子数、光線数（`sampling.rays_per_batch`）、乱数 seed を変える。

判断の基準の決め方は[計算結果の妥当性確認](ValidationGuide.html)にあります。

## 適用範囲

- セルは平らな表面の一部を表します。方式 1〜3 はセルの外のシースを解かず、方式 4 はそれを平面の定常解で表します。
- どの方式も、box 外での粒子の軌道、飛行時間、空間電荷を解きません。
- 方式 1〜3 の電位は、上端面の平均との差としてだけ意味を持ちます。上流のプラズマを基準にした電位が必要なら、方式 4 を使います。
