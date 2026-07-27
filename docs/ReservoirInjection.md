title: reservoir_face の流入量と速度サンプリング

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# `reservoir_face` の流入量と速度サンプリング

このページは、[粒子源の全体像](ParticleSourcesBoundaries.html)で `source_mode="reservoir_face"` を選んだ後の
注入アルゴリズムを説明します。このページで扱うのは、次の 2 点です。

1. 上流分布関数（VDF）から、1 バッチに流入する物理粒子数とマクロ粒子数をどう決めるか。
2. 選ばれた粒子の初期位置を注入開口から作り、面到達時の速度をどう生成するか。

```text
上流 VDF と開口
       ↓ 到達可能な流束を求める
1 バッチの流入量
       ↓ マクロ粒子数へ変換する
開口から作る初期位置・面到達速度
       ↓
共通の粒子追跡へ渡す
```

`reservoir_face` 自体は、粒子源の選択、シースモデルの選択、ボックス外へ出た粒子の脱出・帰還を担当しません。
これらは、このページの末尾にある関連ページで分けて説明します。

## このページの入力と出力

| 種類 | 内容 |
| --- | --- |
| 入力 | 上流 VDF、注入開口、`batch_duration`、マクロ粒子の重み、必要に応じて上流と注入面の電位差 |
| 出力 | 全ランク合計のマクロ粒子数、注入開口から作る初期位置、面到達時の速度、次のバッチへ持ち越す端数 |

上流 VDF は Maxwell 分布または速度グリッドで与えます。どちらも面通過流束として扱い、電位差による補正は必要な場合だけ加えます。

## 開口の位置と流入面積を決める

`inject_face` は `x_low`、`x_high`、`y_low`、`y_high`、`z_low`、`z_high` のいずれかです。
`pos_low` と `pos_high` は、そのボックス面上の矩形開口を定めます。面の内向き単位法線を $\mathbf n$、
開口面積を $A$ とします。

`reservoir_face` は `sim.use_box=true` と `sim.batch_duration>0` を要求します。開口の法線座標はボックス面と一致し、
2 つの接線方向は正の長さを持つ必要があります。

## Maxwell 分布を流入流束で重み付けする

上流のドリフト付き Maxwell 分布で、法線ドリフトを $u_n=\mathbf u\cdot\mathbf n$、
熱速度の標準偏差を $\sigma=\sqrt{k_\mathrm{B}T/m}$ とします。上流法線速度の下限が $v_{\min}$ なら、
片側流入流束は

$$
\Gamma_\mathrm{in}
=n\int_{v_n\ge v_{\min}}v_n f_n(v_n)\,dv_n
=n\left[u_n\{1-\Phi(a)\}+\sigma\varphi(a)\right],
\qquad
a=\frac{v_{\min}-u_n}{\sigma}
$$

です。$\Phi$ と $\varphi$ は標準正規分布の累積分布関数と確率密度です。温度 0 では、
$u_n\ge v_{\min}$ のとき $\Gamma_\mathrm{in}=nu_n$、それ以外は 0 です。

接線速度は通常の Gaussian 分布、法線速度は $v_n f_n(v_n)$ に比例する片側の流束重み付き分布から生成します。
体積中の Maxwell 分布をそのまま面へ置くのではなく、面を横切る確率が法線速度に比例することを反映しています。

## 速度グリッドの値を流束へ変換する

`velocity_distribution="grid"` では、CSV の速度点と非負の値 $f$ から速度を生成します。

| `velocity_grid_pdf_kind` | CSV の $f$ の意味 | 生成時の重み |
| --- | --- | --- |
| `phase_space` | 位相空間分布 | $\max(v_n,0)f$ |
| `flux_weighted` | 流入粒子に対して重み付け済み | $f$ |

どちらも $v_n\ge v_{\min}$ かつ $v_n>0$ の点だけを使います。`velocity_grid_sampling="rectilinear"` は
完全な直交グリッドを補間し、`discrete` は入力点から離散的に選びます。`auto` は、可能な場合に `rectilinear` を使います。

速度グリッドでは、粒子数を `particle_flux_m2_s` または `current_density_a_m2` から決めます。
密度と温度は粒子数の計算に使いません。現行実装では、Zhao 注入補正と併用できません。

## 端数を持ち越して長時間の流束を保つ

Maxwell 分布から求めた流束、または速度グリッドに指定した流束を $\Gamma_\mathrm{in}$ とすると、
実粒子の期待流入数とマクロ粒子の期待数は

$$
N_\mathrm{phys}=\Gamma_\mathrm{in}A\,\Delta t_\mathrm{batch},
\qquad
N_\mathrm{macro,expected}=\frac{N_\mathrm{phys}}{w}
$$

です。整数化で残った端数を $r$ として次のバッチへ持ち越し、

$$
N_\mathrm{macro}=\left\lfloor r+N_\mathrm{macro,expected}\right\rfloor,
\qquad
r\leftarrow r+N_\mathrm{macro,expected}-N_\mathrm{macro}
$$

とします。各バッチの個数は変動しても、多数のバッチを通じた総流入量は期待流束に一致します。

重みは、正の `w_particle` を直接与えるか、`target_macro_particles_per_batch` から逆算します。
後者は計算で追跡する粒子数を調整する設定であり、物理流束を変更しません。`-1` は粒子種 1 で解決した重みを共有します。
2 つの方法を同時に指定することはできません。

## 1 つの電位差で到達条件と注入面速度を決める

無限遠または上流の電位を $\phi_\infty$、注入面または界面の電位を $\phi_f$ とします。
1 次元の静電写像は接線速度を保ち、法線エネルギーを保存します。

$$
\frac12m v_{n,f}^2+q\phi_f
=\frac12m v_{n,\infty}^2+q\phi_\infty,
$$

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}.
$$

| $B$ | 上流粒子の到達条件と注入面速度 |
| ---: | --- |
| $B>0$ | $v_{n,\infty}\ge\sqrt B$ の粒子だけが到達し、注入面までに減速する |
| $B=0$ | 法線速度を変更しない |
| $B<0$ | すべての流入粒子が到達でき、注入面までに加速する |

上流 VDF は $v_n\ge v_{\min}=\sqrt{\max(B,0)}$ を満たす粒子に制限し、採用した速度を同じ $B$ で注入面速度へ
写します。これにより、生成する粒子数と面到達時の速度を 1 つの電位差から矛盾なく決めます。

## 電位差や上流 VDF の補正は別のモデルから受け取る

このページが扱うのは、与えられた上流 VDF と電位差を粒子へ変換する処理です。どの外部モデルを使うかは
[境界・外部領域の構成を選ぶ](OuterPlasmaModels.html)を正本とします。

| 構成 | `reservoir_face` が受け取るもの |
| --- | --- |
| 補正なし | $B=0$ として、設定した VDF を注入面上の分布として使う |
| `external_boundary.particles.inflow_model="infinity_barrier"` | 開口の平均電位と `sim.phi_infty` から求めた 1 つの電位差 |
| split `kinetic_1d` | 収束した外部プロファイルの $\phi_I-\phi_\infty$ |

`infinity_barrier` は、`injection_face_phi_grid_n` を $N$ とする開口内の $N\times N$ セル中心格子で
バッチ開始時の電位を評価し、平均値 $\bar\phi_f$ を使います。固定 P0 三角形カーネル、周期場、
ゼロモード、外部状態、`sim.e0` を含む
同じ電位規約を使いますが、注入面までの途中の $E(z)$、折り返し位置、飛行時間、空間電荷は解きません。
同じ評価から電位の母標準偏差、最小値、最大値も集計し、Maxwell 分布のリザーバーでは面内のばらつきが特徴エネルギーに対して
大きい場合に、初回と最終バッチで警告します。

split 外部モデルでは、z-high の `reservoir_face` を無限遠 VDF として解釈します。`kinetic_1d` では界面電場を
Poisson 境界条件に使いますが、粒子速度を変える量は $\phi_I-\phi_\infty$ です。

## 仮想飛行時間で初期位置をずらす

すべての粒子を注入面上の同じ時刻位置から生成すると、粒子配置に人工的な揃い方が生じます。BEACH は各粒子に
$\tau\in[0,\texttt{sim.dt})$ の仮想飛行時間を一様に与え、追跡開始前に
$\mathbf x\leftarrow\mathbf x+\mathbf v\tau$ として初期位置だけを速度方向へずらします。

実装内部では、この時間幅を `position_jitter_dt=sim.dt` としてサンプラーへ渡します。この内部識別子の `jitter` は
初期位置のずらしを指し、利用者が設定する入力項目ではありません。速度へ乱数ノイズを加える処理でもありません。
全体のシミュレーション時刻と、各粒子を追跡できる時間は変更しません。
両面が周期境界の軸は基本ボックスへ戻し、その他の軸はボックス内へ収めます。

## 全体の端数を MPI と再開で引き継ぐ

MPI では、全ランク合計の粒子数と端数をルートランクで決めてから各ランクへ分配します。
端数は `macro_residuals.csv` に保存し、再開時に全ランクへ配布します。このため、期待流入量と端数列は MPI ランク数に依存しません。

結果を確認するときは、粒子種別の注入個数と注入電荷、解決後のマクロ粒子重み、端数、適用した $v_{\min}$ と電位差、
バッチ平均の吸収電流と脱出電流を区別して確認します。

`batch_duration` を変えると、1 バッチの期待粒子数と場の更新間隔が同時に変わります。
[バッチ幅と安定性](BatchDurationStability.html)に、物理的な定常値の確認方法を示しています。

## 関連ページ

- 粒子源を選ぶ: [粒子源の全体像](ParticleSourcesBoundaries.html)
- 流入補正、流出境界、外部場を組み合わせる: [境界・外部領域の構成を選ぶ](OuterPlasmaModels.html)
- 注入後の粒子が開放面から出たときの処理: [粒子の脱出と帰還](ParticleEscapeReturn.html)

## Code reference

- 流束積分、マクロ粒子数、Maxwell／速度グリッドからの生成: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- 電位障壁と外部プロファイルによる速度補正: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- 開口形状と入力組み合わせの検証: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- MPI 全体のマクロ粒子端数: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
