title: 境界から粒子を流入させる

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# 境界から粒子を流入させる

> **質問:** simulation box の外にある plasma から、密度・温度または速度分布に対応する粒子を
> box 境界へ流入させるには、何を設定すればよいですか。
>
> **一文での回答:** 非周期の box 面を `open` にし、species の
> `[particles.species.boundary_inflow]` に `"reservoir"` を指定します。

読了後には、最小設定を作り、Maxwell 分布と速度 grid、境界上の分布と無限遠からの電位補正を
選び、出力から実際の注入量を確認できます。粒子を box 内部の面から供給したい場合は、先に
[粒子をどこから入れるか](ParticleSourcesBoundaries.html)で `plane_source` との違いを確認してください。

## 1. 最小設定を作る

次は、既存の有効なケースへ z-high 面から電子を流入させる最小差分です。値は設定方法を示す例であり、
特定の plasma 環境を表す参照値ではありません。

```toml
[sim]
batch_duration = 1.0e-6

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]
periodic_axes = []

[particle_boundary]
z_high = "open"
ordinary_open_model = "escape"

[[particles.species]]
species_key = "electron"
source_mode = "volume_seed"
npcls_per_step = 0
velocity_distribution = "maxwellian"
number_density_m3 = 5.0e6
temperature_ev = 10.0
drift_velocity = [0.0, 0.0, -4.0e5]
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
w_particle = 1.0e5

[particles.species.boundary_inflow]
z_high = "reservoir"
```

z-high 面の内向き法線は $-z$ 方向なので、この例の負の z ドリフトは box 内向きです。
`x_low`、`x_high`、`y_low`、`y_high`、`z_low`、`z_high` の複数面を同時に選べます。

設定時には次の条件を守ります。

- 解決後の `sim.batch_duration` を正にする。
- 流入面を `domain.periodic_axes` に含めず、species override を含む有効な外向き作用を `open` にする。
- `boundary_inflow` は box 面全体を使う。`pos_low`、`pos_high`、`inject_face` では開口を切り出さない。
- 現行 schema では `source_mode="volume_seed"` が必要である。境界流入だけなら
  `npcls_per_step=0` とし、これは体積中に粒子を生成する指定ではない。

設定を保存したら、通常の実行前に構造と組合せを検査します。

```bash
beachx lint beach.toml
beach beach.toml
```

`lint` が成功しても、流束、マクロ粒子重み、電位基準が物理的に適切であることまでは保証しません。

## 2. Maxwell 分布または速度 grid を選ぶ

| 手元にある外部 plasma 情報 | 選ぶ設定 |
| --- | --- |
| 密度、温度、ドリフト速度 | `velocity_distribution="maxwellian"` |
| 計測・別計算による速度点と分布値 | `velocity_distribution="grid"` |

### Maxwell 分布

Maxwell 分布では、`number_density_m3` または `number_density_cm3`、温度、`drift_velocity` から
各面の片側流束 $\Gamma_\mathrm{in}$ を求めます。境界を通過する確率は内向き法線速度に比例するため、
法線速度は体積中の Maxwell 分布をそのまま置くのではなく、流束で重み付けした分布から標本化します。

1 batch の期待マクロ粒子数は、利用時の判断に必要な次の関係で決まります。

$$
N_\mathrm{macro,expected}
=\frac{\Gamma_\mathrm{in} A\,\Delta t_\mathrm{batch}}{w},
$$

ここで $A$ は選択面の面積、$\Delta t_\mathrm{batch}$ は `batch_duration`、$w$ は
`w_particle` です。整数化の端数は次の batch へ持ち越されるため、各 batch の注入数は同じとは限りません。
完全な流束式と制約は[入力パラメータ](Parameters.html#particlesspeciesboundary_inflow)と
[`SPEC.md`](../SPEC.md)を参照してください。

物理的なマクロ粒子重みが決まっている場合は `w_particle` を指定します。計算量から 1 batch の標本数を
決めたい場合は、代わりに `target_macro_particles_per_batch` を指定します。後者は重みを解決するための
target であり、物理流束を変更しません。二つは同時に指定できません。

### 速度 grid

速度 grid では、Maxwell 用の密度・温度を削除し、CSV と物理流束を指定します。

```toml
velocity_distribution = "grid"
velocity_grid_path = "inflow_vdf.csv"
velocity_grid_pdf_kind = "phase_space"
velocity_grid_sampling = "auto"
particle_flux_m2_s = 1.0e12
```

CSV の必要列は `vx_m_s,vy_m_s,vz_m_s,f` です。`f` は非負とし、次のどちらを表すかを明示します。

| `velocity_grid_pdf_kind` | CSV の `f` | BEACH の重み |
| --- | --- | --- |
| `phase_space` | 位相空間分布 | 内向き法線速度を掛けた $\max(v_n,0)f$ |
| `flux_weighted` | 境界を通過する粒子について既に重み付けした分布 | $f$ |

流入量には `particle_flux_m2_s` または `current_density_a_m2` の片方だけを指定します。電流密度は
$|J/q|$ で粒子流束へ変換されるため、その符号で速度方向を指定することはできません。方向は CSV の速度と
選択面の内向き法線で決まります。相対 CSV path は実行時の working directory を基準にします。

## 3. `source_vdf` または `infinity_barrier` を選ぶ

| 外部 VDF をどこで定義したか | `reservoir.inflow_model` | 挙動 |
| --- | --- | --- |
| simulation 境界上 | `"source_vdf"` | 設定した VDF をそのまま境界上の分布として使う。既定値 |
| 無限遠または上流 | `"infinity_barrier"` | 上流基準電位から境界面までの scalar な電位差で到達流束と法線速度を補正する |

`infinity_barrier` を使う場合は、外部基準電位を明示します。

```toml
[reservoir]
inflow_model = "infinity_barrier"
phi_infty = 0.0
face_potential_grid_n = 5
```

上流電位を $\phi_\infty$、batch 開始時の流入面平均電位を $\phi_f$ とすると、境界での法線速度は

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}
$$

で決まります。$q$ は符号を含む粒子電荷なので、electron と正 ion では同じ電位差に対する $B$ の符号が
逆になります。

| $B$ | 到達条件と境界での変化 |
| ---: | --- |
| $B>0$ | $v_{n,\infty}\ge\sqrt B$ の粒子だけが到達し、境界までに減速する |
| $B=0$ | 法線速度を変更しない |
| $B<0$ | すべての内向き粒子が到達でき、境界までに加速する |

接線速度は変更しません。$\phi_f$ は `face_potential_grid_n` による面内標本の平均であり、粒子ごとの
局所電位ではありません。

流入と流出は別の判断です。外向き粒子にも電位障壁を適用する場合は、`[particle_boundary]` の
`ordinary_open_model="potential_barrier"` を選びます。こちらは粒子が実際に境界を通過した位置の
局所電位を使って return または escape を決めます。

## 4. 出力で注入を確認する

`output.dir` を `outputs/latest` とした場合は、実行後に次を確認します。

```bash
beachx inspect outputs/latest
grep -E '^(reservoir_inflow_map|particle_ordinary_open_model|charge_ledger_residual_C)=' \
  outputs/latest/summary.txt
head -n 2 outputs/latest/charge_ledger.csv
```

`summary.txt` の `reservoir_inflow_map` は、選択に応じて `source_vdf` または `infinity_barrier` になります。
`particle_ordinary_open_model` は独立に `escape` または `potential_barrier` を示します。

`charge_ledger.csv` では、species ごとに少なくとも次を確認します。

- `injected_count`: この例のように境界流入だけを使う species では、box 外から生成されたマクロ粒子数
- `injected_from_remote_C`: その符号付き注入電荷。符号は `q_particle` と一致する
- `absorbed_count`、`escaped_count`、`discarded_unresolved_count`: 注入後の主な行き先

期待数が 1 より十分小さい場合、端数が蓄積するまで初期 batch の `injected_count` は 0 になり得ます。
必要な標本数が得られないときは、物理流束を変える前に `w_particle` または
`target_macro_particles_per_batch` を見直します。`batch_duration` の変更は粒子数だけでなく場の更新間隔も
変えるため、[batch 幅と安定性](BatchDurationStability.html)に従って比較してください。

出力が生成されたことは実行完了を示します。流束、電荷収支、統計誤差、結果の数値・物理的妥当性は
[出力ファイルを調べる](OutputGuide.html)と[計算結果の妥当性確認](ValidationGuide.html)で別に評価します。

## 5. 適用範囲を確認する

- `boundary_inflow` が表すのは、box 外の条件を境界上の粒子へ写す局所モデルです。box 外の軌道、
  途中の $E(z)$、折り返し位置、飛行時間、空間電荷、自己無撞着な外部 sheath は解きません。
- 一様外部電場には有限な無限遠電位がありません。`infinity_barrier` と併用する場合は、`phi_infty` を
  有効な reservoir 基準として定義し、その物理的意味を別途検証します。
- 流入位置は選択した box 面全体です。box 内側に置いた矩形面から供給する場合は
  `source_mode="plane_source"` を使います。
- `boundary_inflow` は外向き境界作用を上書きせず、`[outer_plasma]` や `[coupling]` を有効にする設定でも
  ありません。これらの削除済み table は入力として拒否されます。

全 key、既定値、排他条件、端数と checkpoint の形式は[入力パラメータ](Parameters.html)と
[`SPEC.md`](../SPEC.md)が正本です。MPI 分配、乱数状態、端数状態、実装の責務は
[開発者向けアーキテクチャ](Architecture.html)へ分離しています。

## 6. 旧 `reservoir_face` から移行する

`source_mode="reservoir_face"` は既存 case のための deprecated 入力です。新しい外部 plasma 条件は
`[particles.species.boundary_inflow]`、box 内部の矩形放出面は `source_mode="plane_source"` へ移します。

旧形式は `inject_face` と `pos_low` / `pos_high` で開口を持てますが、`boundary_inflow` は box 面全体を
使います。移行前後で面積が変わる場合は、設定密度または流束だけでなく、`injected_count` と
`injected_from_remote_C` を比較して総流入量が意図どおりか確認してください。BEACH は旧設定を
どちらの新形式にも暗黙変換しません。互換動作と完全な移行制約は
[入力パラメータ](Parameters.html#source_mode--reservoir_facedeprecated)と [`SPEC.md`](../SPEC.md)を参照してください。
