title: シミュレーション境界からの reservoir 流入

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# シミュレーション境界からの reservoir 流入

`[particles.species.boundary_inflow]` は、ボックス外側に与えた plasma reservoir からシミュレーション境界を
横切る粒子流入を定義します。上流 VDF、ボックス面積、`batch_duration`、マクロ粒子重みから、面ごとの注入数、
初期位置、面到達速度、次のバッチへ持ち越す端数を求めます。

この流入は `source_mode` と責務を分離しています。初版で同じ species に併用できるのは
`source_mode="volume_seed"` だけです。内部の明示的な放出面には `source_mode="plane_source"` を使いますが、
同じ species の境界流入との併用は fail closed です。
粒子源の選択は[粒子源と境界流入](ParticleSourcesBoundaries.html)、注入後の open 面処理は
[粒子の escape と局所 return](ParticleEscapeReturn.html)を参照してください。

## 流入面を選ぶ

各面の値に `"reservoir"` を指定すると、そのボックス面全体から流入します。

```toml
[[particles.species]]
species_key = "electron"
source_mode = "volume_seed"
npcls_per_step = 100
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
w_particle = 1.0e5
drift_velocity = [0.0, 0.0, -4.0e5]

[particles.species.boundary_inflow]
z_high = "reservoir"
```

`x_low`、`x_high`、`y_low`、`y_high`、`z_low`、`z_high` を指定でき、複数面を同時に有効化できます。
`domain.periodic_axes` に属する面には reservoir 流入を設定できません。流入位置は選択したボックス面内で一様に
生成されるため、`pos_low`、`pos_high`、`inject_face` で開口を切り出しません。

外向き粒子の処理は別の設定です。global 作用は `[particle_boundary]`、species override は
`[particles.species.boundary]` で指定します。reservoir 流入面の有効な外向き作用は `open` でなければ
ならず、`boundary_inflow` 自体は外向き作用を上書きしません。

## Maxwell 分布を流入流束で重み付けする

選択面の内向き単位法線を $\mathbf n$、面積を $A$ とします。上流のドリフト付き Maxwell 分布について、
法線ドリフトを $u_n=\mathbf u\cdot\mathbf n$、熱速度の標準偏差を
$\sigma=\sqrt{k_\mathrm{B}T/m}$ とします。上流法線速度の下限が $v_{\min}$ なら、片側流入流束は

$$
\Gamma_\mathrm{in}
=n\int_{v_n\ge v_{\min}}v_n f_n(v_n)\,dv_n
=n\left[u_n\{1-\Phi(a)\}+\sigma\varphi(a)\right],
\qquad
a=\frac{v_{\min}-u_n}{\sigma}
$$

です。$\Phi$ と $\varphi$ は標準正規分布の累積分布関数と確率密度です。温度 0 では、
$u_n\ge v_{\min}$ のとき $\Gamma_\mathrm{in}=nu_n$、それ以外は 0 です。

接線速度は Gaussian 分布、法線速度は $v_n f_n(v_n)$ に比例する片側分布から生成します。これは面通過確率が
法線速度に比例するためで、体積中の Maxwell 分布をそのまま境界へ置く方法とは異なります。

## 速度グリッドの値を流束へ変換する

`velocity_distribution="grid"` では、CSV の速度点と非負の値 $f$ から速度を生成します。

| `velocity_grid_pdf_kind` | CSV の $f$ の意味 | 生成時の重み |
| --- | --- | --- |
| `phase_space` | 位相空間分布 | $\max(v_n,0)f$ |
| `flux_weighted` | 流入粒子に対して重み付け済み | $f$ |

どちらも $v_n\ge v_{\min}$ かつ $v_n>0$ の点だけを使います。`velocity_grid_sampling="rectilinear"` は
完全な直交グリッドを補間し、`discrete` は入力点から離散的に選びます。`auto` は可能な場合に
`rectilinear` を使います。

速度グリッドでは、粒子数を `particle_flux_m2_s` または `current_density_a_m2` から決めます。
密度と温度は粒子数の計算に使いません。

## 面ごとの端数を持ち越して長時間の流束を保つ

Maxwell 分布から求めた流束、または速度グリッドに指定した流束を $\Gamma_\mathrm{in}$ とすると、
実粒子の期待流入数とマクロ粒子の期待数は

$$
N_\mathrm{phys}=\Gamma_\mathrm{in}A\,\Delta t_\mathrm{batch},
\qquad
N_\mathrm{macro,expected}=\frac{N_\mathrm{phys}}{w}
$$

です。整数化で残った端数 $r$ は species と face の組ごとに次のバッチへ持ち越します。

$$
N_\mathrm{macro}=\left\lfloor r+N_\mathrm{macro,expected}\right\rfloor,
\qquad
r\leftarrow r+N_\mathrm{macro,expected}-N_\mathrm{macro}
$$

各バッチの個数は変動しても、長時間の総流入量は期待流束へ収束します。正の `w_particle` を直接与えるか、
`target_macro_particles_per_batch` から逆算します。後者は有効な全流入面の合計 sample 数を調整し、
物理流束は変えません。

## 境界条件として流入補正を選ぶ

`[reservoir]` は境界流入と、open 面から出る粒子の scalar barrier に共通する外部 plasma 条件です。
任意の内部平面では、面のどちら側を外部 plasma とみなし、どの電位を上流基準へ対応させるかが一意に決まりません。
simulation boundary なら box 外側を reservoir と定義できるため、`phi_infty` と面電位による
`infinity_barrier` を外部から境界までの補正として自然に適用できます。この理由から、
新しい設定の流入側 potential/barrier 補正は `boundary_inflow` に適用し、`plane_source` には適用しません。
非推奨の `reservoir_face` は既存 case の数値挙動を保つため、従来の補正を互換動作として維持します。

```toml
[particle_boundary]
z_high = "open"
ordinary_open_model = "potential_barrier"

[reservoir]
inflow_model = "infinity_barrier"
phi_infty = 0.0
face_potential_grid_n = 5
```

| 設定 | 境界流入での作用 |
| --- | --- |
| `reservoir.inflow_model="source_vdf"` | 設定した VDF を境界上の分布として使い、電位補正を加えない |
| `reservoir.inflow_model="infinity_barrier"` | 面平均電位と `reservoir.phi_infty` の電位差で、到達流束と法線速度を補正する |

無限遠または上流の電位を $\phi_\infty$、境界面の電位を $\phi_f$ とすると、1 次元の静電写像は
接線速度を保ち、法線エネルギーを保存します。

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}
$$

| $B$ | 上流粒子の到達条件と境界面速度 |
| ---: | --- |
| $B>0$ | $v_{n,\infty}\ge\sqrt B$ の粒子だけが到達し、境界までに減速する |
| $B=0$ | 法線速度を変更しない |
| $B<0$ | すべての流入粒子が到達でき、境界までに加速する |

`infinity_barrier` は、各流入面を `reservoir.face_potential_grid_n` による $N\times N$ cell-center で評価し、
batch 開始時の平均値 $\bar\phi_f$ を使います。固定 P0 三角形 kernel、周期場、zero mode、`sim.e0` を含む
同じ電位規約を使います。

`particle_boundary.ordinary_open_model="potential_barrier"` は、同じ `phi_infty` と粒子が実際に通過した位置の
局所電位を使い、外向き粒子の return または escape を判定します。流入側の `infinity_barrier` は面平均、
流出側の `potential_barrier` は event 位置という違いがあります。

`[surface_current_model] model="zhao_stationary"`は、参照したambient electron/ionのz-highにmodel固有の写像を
追加します。無限遠0 Vから現在の面平均電位まで速度を写像し、Type A electronでは$\phi_m$もaccess bottleneckに
含めます。この写像は選択speciesのz-highだけに作用し、`reservoir.inflow_model`の一般設定を他面について変更しません。

## このモデルが表す範囲

境界流入として定義することで、密度、温度、VDF、`phi_infty` を「ボックス外側の条件」として一貫して解釈できます。
ただし、BEACH が解くのは境界までの scalar な速度写像と、境界通過後の粒子運動です。ボックス外側の軌道、
途中の $E(z)$、折り返し位置、飛行時間、空間電荷、自己無撞着な外部 sheath は解きません。
`[outer_plasma]` や `[coupling]` を復活させる設定でもありません。

一様外部電場には有限な無限遠電位がないため、`infinity_barrier` と併用する場合は `phi_infty` を有効な
reservoir 基準として定義し、物理的な妥当性を別途確認します。

## 初期位置、MPI、再開

同じ時刻に全粒子を境界へ置く人工的な整列を避けるため、各粒子へ
$\tau\in[0,\texttt{sim.dt})$ の仮想飛行時間を一様に与え、追跡前に
$\mathbf x\leftarrow\mathbf x+\mathbf v\tau$ として初期位置だけを速度方向へずらします。
全体のシミュレーション時刻と各粒子の追跡時間は変更しません。

MPI では全 rank 合計の粒子数と端数を root rank で決めてから分配します。端数は checkpoint に保存し、
再開時に復元します。結果を確認するときは、species・face 別の注入個数と注入電荷、解決後のマクロ粒子重み、
端数、適用した $v_{\min}$ と電位差、吸収電流と escape 電流を区別します。

checkpoint schema v6 の `macro_residuals.csv` は `species_idx,face,residual` で、`face=0` は従来 source、
`1..6` は boundary face です。旧 2 列形式は読み込み互換です。

`batch_duration` を変えると、1 バッチの期待粒子数と場の更新間隔が同時に変わります。
[バッチ幅と安定性](BatchDurationStability.html)に、定常値の確認方法を示しています。

## 旧 `reservoir_face` から移行する

`source_mode="reservoir_face"` は互換入力として残りますが、新しい設定では
`[particles.species.boundary_inflow]` を使ってください。旧設定は `inject_face` と
`pos_low` / `pos_high` による開口を持つため、境界全面を使う新形式への移行時は粒子流入量が変わらないか確認します。
内部に矩形放出面を置く用途は `source_mode="plane_source"` へ移します。BEACH は旧設定をどちらにも暗黙変換しません。

## Code reference

- 流束積分、マクロ粒子数、Maxwell／速度グリッドからの生成: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- 境界流入と電位障壁による速度補正: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- 入力組合せの検証: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- MPI 全体のマクロ粒子端数: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
