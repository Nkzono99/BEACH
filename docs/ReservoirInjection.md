title: 境界から粒子を流入させる

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# 境界から粒子を流入させる

box の外にあるプラズマ（外部プラズマ、reservoir）から、密度・温度または速度分布に対応する粒子を box の面から入れる手順です。
非周期の面を開放にし、粒子種の `[particles.species.inflow]` でその面を `"reservoir"` にします。
このページを読むと、最小の設定を作り、分布と流入の写像を選び、実際に入った量を出力で確かめられます。
box の内側の面から粒子を出したい場合は、[粒子源を選ぶ](ParticleSourcesBoundaries.html)の `plane_source` を使います。

## 1. 最小の設定を作る

既存のケースに、上端面から電子を入れる差分です。値は例で、特定のプラズマ環境を表すものではありません。

```toml
[run.batch]
duration_s = 1.0e-6

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]
periodic_axes = []

[particles.boundary]
z_high = "open"

[[particles.species]]
species_key = "electron"
charge_c = -1.602176634e-19
mass_kg = 9.1093837139e-31

[particles.species.distribution]
model = "maxwellian"
number_density_m3 = 5.0e6
temperature_ev = 10.0
drift_velocity_m_s = [0.0, 0.0, -4.0e5]

[particles.species.sampling]
weight = 1.0e5

[particles.species.inflow]
z_high = "reservoir"
```

上端面の内向き法線は $-z$ なので、負の z の drift は box の内向きです。6 つの面のうち複数を同時に選べます。

- バッチ幅（`run.batch.duration_s`）を正にする。入る粒子数はバッチ幅に比例する。
- 流入面は非周期で、その粒子種にとって開放面にする。
- 流入は選んだ面の全体から入る。面の一部だけから入れることはできない。
- 境界流入だけを使う粒子種は、`[particles.species.source]` を書かない。

設定を保存したら、検査してから実行します。

```bash
beachx lint beach.toml
beach beach.toml
```

`lint` が通っても、流束、重み、電位の基準が物理的に適切だとは限りません。

## 2. 分布を選ぶ

| 手元にある外部プラズマの情報 | 選ぶ設定 |
|---|---|
| 密度、温度、drift 速度 | `distribution.model = "maxwellian"` |
| 計測や別の計算による速度点と分布の値 | `distribution.model = "grid"` |

### Maxwell 分布

密度（`number_density_m3` または `number_density_cm3`）、温度、drift 速度から、各面を内向きに横切る束
$\Gamma_\mathrm{in}$ を求めます。面を横切る確率は内向きの法線速度に比例するので、法線速度は束で重み付けした分布から選びます。
1 バッチで入るマクロ粒子数の期待値は

$$
N_\mathrm{macro}
=\frac{\Gamma_\mathrm{in} A\,\Delta t_\mathrm{batch}}{w}
$$

です。$A$ は面の面積、$\Delta t_\mathrm{batch}$ はバッチ幅、$w$ は重みです。端数は次のバッチへ持ち越すので、
バッチごとの粒子数は一定とは限りません。

重みを物理的に決めたい場合は `sampling.weight`、1 バッチの標本数から決めたい場合は
`sampling.target_macro_particles_per_batch` を指定します。後者は重みを決めるための目標で、物理的な束は変えません。
両方は指定できません。

### 速度 grid

Maxwell 分布用の密度と温度の代わりに、CSV と物理的な束を指定します。

```toml
[particles.species.distribution]
model = "grid"
grid_path = "inflow_vdf.csv"
grid_pdf_kind = "phase_space"
particle_flux_m2_s = 1.0e12

[particles.species.sampling]
velocity_grid_sampling = "auto"
```

CSV の列は `vx_m_s,vy_m_s,vz_m_s,f` です。`f` は非負で、何を表すかを `grid_pdf_kind` で指定します。

| `grid_pdf_kind` | CSV の `f` | BEACH の重み |
|---|---|---|
| `phase_space` | 位相空間の分布 | 内向き法線速度を掛けた $\max(v_n,0)f$ |
| `flux_weighted` | 面を横切る粒子について重み付け済みの分布 | $f$ |

流入量は `particle_flux_m2_s` か `current_density_a_m2` の一方だけで指定します。電流密度は $|J/q|$ で束に直すので、
符号で向きは決まりません。向きは CSV の速度と面の内向き法線で決まります。CSV の相対パスは、実行時の作業ディレクトリが基準です。

## 3. 流入の写像を選ぶ

外部プラズマの分布をどこで定義したかで、`particles.reservoir.inflow_model` を選びます。

| 分布を定義した場所 | `inflow_model` | 動作 |
|---|---|---|
| 流入面の上 | `"source_vdf"`（既定値） | 設定した分布をそのまま流入面の分布として使う |
| 無限遠（上流） | `"infinity_barrier"` | 上流電位と流入面の平均電位の差で、届く粒子を選び、法線速度を変える |

```toml
[particles.reservoir]
inflow_model = "infinity_barrier"
phi_infty_v = 0.0
face_potential_grid_n = 5
```

上流電位を $\phi_\infty$、バッチの初めの流入面の平均電位を $\phi_f$ とすると、流入面での法線速度は

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}
$$

です。$q$ は符号付きの電荷なので、電子と正イオンでは同じ電位差に対して $B$ の符号が逆になります。

| $B$ | 届く粒子と速度の変化 |
|---:|---|
| $B>0$ | $v_{n,\infty}\ge\sqrt B$ の粒子だけが届き、流入面までに減速する |
| $B=0$ | 法線速度を変えない |
| $B<0$ | すべての粒子が届き、流入面までに加速する |

接線速度は変えません。$\phi_f$ は流入面上の `face_potential_grid_n` × `face_potential_grid_n` 点の平均で、粒子ごとの局所電位ではありません。
出ていく粒子も同じ上流電位で判定するには、開放面の扱いを電位障壁にします（[box 境界での粒子の扱い](ParticleEscapeReturn.html#電位障壁potential_barrier)）。
外部シースとの接続を使う場合は、流入の写像を外部シースが決めるので、`infinity_barrier` は使えません。

## 4. 入った量を出力で確かめる

```bash
beachx inspect outputs/latest
grep -E '^(reservoir_inflow_map|particle_ordinary_open_model|charge_ledger_residual_C)=' \
  outputs/latest/summary.txt
head -n 2 outputs/latest/charge_ledger.csv
```

`summary.txt` の `reservoir_inflow_map` は選んだ写像（`source_vdf` または `infinity_barrier`）、
`particle_ordinary_open_model` は開放面の扱いです。`charge_ledger.csv` では、粒子種ごとに次を確かめます。

- `injected_count`: box の外から入ったマクロ粒子数
- `injected_from_remote_C`: その電荷。符号は粒子の電荷と同じ
- `absorbed_count`、`escaped_count`、`discarded_unresolved_count`: 入った粒子の行き先

期待値が 1 より十分小さいと、端数がたまるまで最初のバッチの `injected_count` は 0 になり得ます。標本が足りないときは、
物理的な束を変える前に重みか目標の標本数を見直します。バッチ幅を変えると場を更新する間隔も変わるので、
[バッチ幅を決める](BatchDurationStability.html)に従って比べます。

## 5. 適用範囲

- 境界流入は、box の外の条件を境界上の粒子に置き換える局所的なモデルです。box の外の軌道、途中の電場、
  折り返し位置、飛行時間、空間電荷、外部シースは解きません。
- 一様な外部電場には無限遠の電位が決まりません。`infinity_barrier` と併用する場合は、`phi_infty_v` を外部プラズマの
  基準として決め、その意味を別に確かめます。

全キーと制約は[入力パラメータ](Parameters.html#particlesspeciesboundary_inflow)にあります。
