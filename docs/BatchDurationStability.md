title: `batch_duration` の安定性と定常値

Lang: [日本語](BatchDurationStability.md) | [English](BatchDurationStability.en.md)

# `batch_duration` の安定性と定常値

通常の固定幅実行では、`sim.batch_duration`は1 batchで進める物理時間であり、表面電荷を更新する時間幅です。
`sim.batch_duration_step`を指定した場合は、
`sim.batch_duration = sim.dt * sim.batch_duration_step`として計算されます。
`reservoir_face`注入では、この時間内の物理流入量からマクロ粒子数または粒子重みを決めます。

`[periodic2].max_nonzero_mode_potential_step > 0`を指定した場合、`sim.batch_duration`は
各accepted batchの**最大**試行幅です。実際の進行幅は、`batch_duration`から2分の1ずつ短くする
固定ladder上で選ばれます。

`batch_duration`を小さくすると表面電荷の更新は安定しやすくなりますが、定常状態までに必要なbatch数は
増えます。値を決めるには、更新の安定性とMonte Carloノイズを分けて確認します。

## 実用ガイド

まず、次の手順で実用的な値を絞り込みます。

1. `sim.batch_duration_step` を小さめにして短い計算を走らせる。
2. `charge_history.csv`、`summary.txt` の `last_rel_change`、吸収数・脱出数を見る。
3. `batch_duration_step` を 2 倍と 1/2 倍に振り、最終電荷分布と履歴の形を比較する。
4. 最終電荷と履歴が大きく変わらず、振動・発散傾向がなければ採用する。
5. 揺らぎが大きい場合は、`batch_duration` だけでなく `target_macro_particles_per_batch`、`w_particle`、`batch_count`、`history_stride` も調整する。

| 症状 | 主な確認先 | 典型的な対応 |
| --- | --- | --- |
| 電荷履歴が batch ごとに大きく振動する | `charge_history.csv` | `batch_duration_step` を下げる |
| 最終電荷が step size で大きく変わる | 1/2 倍・2 倍比較 | より小さい `batch_duration_step` で再計算 |
| 履歴がノイズで読みにくい | `target_macro_particles_per_batch`, `w_particle` | マクロ粒子数や重みを調整 |
| 収束前に打ち切られる | `summary.txt` の `batches` | `batch_count` を増やす |

固定幅実行では、`batch_duration`はdeterministicなexplicit更新の時間刻みです。一方、マクロ粒子数と
粒子重みはMonte Carloノイズを主に左右します。

## periodic2非零モードの適応的な進行

productionの`cached_kneq0`計算では、局所電位の1 batch変化へ上限を設定できます。

```toml
[periodic2]
nonzero_mode_backend = "cached_kneq0"
max_nonzero_mode_potential_step = 1.0e-2 # V
```

解決後の`sim.batch_duration`を$h_0$とし、各accepted batchで
$h_0,h_0/2,h_0/4,\ldots$を順に試します。trial幅$h$から粒子数または重みを再計算し、
候補電荷とbatch開始電荷の差が作る$k\ne0$電位を全panel重心で評価します。最大絶対値が
`max_nonzero_mode_potential_step`以下となる最初のtrialを受理します。

棄却trialはRNG、macro粒子数残差、outer state、`implicit_mean` transactionを含めて完全に
rollbackされ、統計、履歴、charge ledgerへ現れません。したがって`batch_count`はaccepted batch数であり、
物理的な終了時刻は`batch_count * batch_duration`ではなく`summary.txt`の`simulated_time_s`です。
再開時は、加算順と受理ladderを再現するため、checkpoint作成時と同じ実OpenMP team sizeを使います。

この上限は、batch内で$k\ne0$場を凍結する近似に対する**局所電位trust bound**です。
局所打切り誤差や大域誤差の次数を保証しません。値を採用する前に
`max_nonzero_mode_potential_step`を半分にし、同じ`simulated_time_s`付近で表面電荷、局所電位幅、
総電荷、interface電位を比較します。固定幅計算との比較では、上限を省略または`0`にします。
この判定は$k=0$更新の安定性を制御しません。`sim.batch_duration`の最大値は、
explicitな平均電荷更新または選択したimplicit closureについて、別途安定である必要があります。

この経路は`cached_kneq0`、time-scaledな`reservoir_face` / `photo_raycast`を要求します。
`reservoir_face`は`target_macro_particles_per_batch`を指定する必要があり、固定`w_particle`は
macro粒子1個の電荷量がhalvingで縮まないため拒否します。explicit SW/UV更新と`implicit_mean`
PE更新に対応し、`volume_seed`とouter event queueは拒否します。

`implicit_mean` Zhao closureでは、大きいtrial幅によってfrozen ambient cohortまたは
frozen interfaceの電場・障壁座標がtrust regionを外れた場合に限り、同じladderで再試行します。
時間幅で縮まらない計測source規格化の変化は再試行しません。
nonmonotone barrier、物理root不在、数値failure、入力不正はhalvingで回復すると仮定せず、その場で停止します。
`BEACH adaptive-kneq0 reject`行の`max_delta_phi_V`は$k\ne0$上限による棄却、
`implicit_status`と`reason`はfrozen-cohort trust-region recoveryを表すため、検証時は別々に数えます。

target方式のreservoirと固定`rays_per_batch`のphoto sourceでは、幅を半分にするとmacro粒子電荷が
半分になり、同じ物理時間あたりの標本数は増えます。したがって上限半減比較には時間離散化差だけでなく
Monte Carlo分散の低下も含まれます。同じ乱数seedを使い、電荷分布のnormと粒子統計を併記してください。

## 連続時間モデルとの対応

絶縁体壁面の要素$j$に蓄積した電荷を$q_j(t)$、単位面積あたりの入射電荷fluxを
$J_j(\mathbf q)$とします。吸収だけを考える基本modelは、次のODEで表せます。

$$
\frac{dq_j}{dt} \;=\; J_j(\mathbf q)\, A_j
$$

ここで$A_j$は要素面積です。$J$は壁面電荷が作る電場に依存するため、一般には
**非線形**です（$J = J(\mathbf q)$）。

BEACHの1 batchは、batch先頭の場を固定したexplicit更新に対応します。平均更新は

$$
\mathbf q^{n+1} \;=\; \mathbf q^n \;+\; \Delta t_b \cdot \mathbf J(\mathbf q^n)\,\mathbf A \;+\; \boldsymbol\eta^n
$$

と表せます。ここで

- $\Delta t_b = $ 固定幅実行では`sim.batch_duration`、適応進行ではaccepted trial幅
- $\mathbf A$ は要素面積ベクトル
- $\boldsymbol\eta^n$ はバッチ内 Monte Carlo サンプリング誤差

です。batch内の粒子は同じ場$E(\mathbf q^n)$を使い、batch末尾で電荷差分をまとめて壁面へ反映します。
したがって、固定幅実行では`batch_duration`、適応進行では受理されたtrial幅がこのexplicit更新の
時間刻みに相当します。

## 定常値が`batch_duration`に依存しない条件

平均更新写像を

$$
F_{\Delta t_b}(\mathbf q) \;=\; \mathbf q \;+\; \Delta t_b\, \mathbf J(\mathbf q)\,\mathbf A
$$

と書くと、その不動点 $\mathbf q^{\ast}$ は

$$
F_{\Delta t_b}(\mathbf q^{\ast}) = \mathbf q^{\ast}
\quad\Longleftrightarrow\quad
\mathbf J(\mathbf q^{\ast}) = 0
$$

で与えられます。したがって、**平均モデルの不動点自体は $\Delta t_b$ に依存しません**。

したがって、反復が安定に収束し、Monte Carlo誤差が十分に平均化されていれば、
`batch_duration`を変えても連続時間モデルの定常解は変わりません。

この結論は、平均モデルの固定点に対するものです。実際の計算には

- バッチごとの有限サンプル誤差
- 収束判定に用いる監視量の揺らぎ
- 有限バッチ数で打ち切ることによる残差

が加わるため、実際に観測される収束値には弱い`batch_duration`依存性が残ることがあります。
平均モデルの固定点と、有限サンプル・有限時間で得た計算値を区別して解釈する必要があります。

## 線形安定性を決める条件

不動点 $\mathbf q^{\ast}$ 近傍で摂動 $\delta\mathbf q^n = \mathbf q^n - \mathbf q^{\ast}$ を考えると、平均更新の線形化は

$$
\delta \mathbf q^{n+1} \;=\; \bigl(I + \Delta t_b\, M\bigr)\,\delta \mathbf q^n,
\qquad
M_{ij} \;\equiv\; \frac{\partial (J_i A_i)}{\partial q_j}\bigg|_{\mathbf q^{\ast}}
$$

となります。一般の多自由度系での安定条件は、この更新行列のスペクトル半径に対する条件

$$
\rho\!\left(I + \Delta t_b\, M\right) < 1
$$

であり、各固有値 $\lambda_k$ に対して

$$
|1 + \Delta t_b\, \lambda_k| < 1
$$

が必要です。これが多自由度系におけるbatch更新の安定条件です。

絶縁体壁に蓄積した電荷が同種粒子の流入を抑える場合、負のフィードバックが働きます。このとき、$M$の
主要な固有値$\lambda_k$は実部が負（$\mathrm{Re}(\lambda_k) < 0$）になると期待されます。
この実負優勢モードの仮定が成り立つ場合に限り、応答時間スケール
$\tau_k \equiv 1/|\lambda_k|$を用いて、最速モードに対する条件を次のように表せます。

$$
0 \;<\; \Delta t_b \;<\; \frac{2}{|\lambda_{\max}|} \;=\; 2\,\tau_{\min}
$$

で発散を避けられ、

$$
0 \;<\; \Delta t_b \;<\; \frac{1}{|\lambda_{\max}|} \;=\; \tau_{\min}
$$

で単調収束（過減衰）になります。

実際の設定では、次のように使い分けます。

- $\Delta t_b < 2\,\tau_{\min}$ : 実負優勢モード仮定のもとでの非発散条件
- $\Delta t_b < \tau_{\min}$ : 同仮定のもとでの単調収束条件
- 一般の結合系では、厳密には $\rho(I + \Delta t_b\, M) < 1$ が本体

$2\tau$ / $\tau$条件は、実負優勢モードを仮定したexplicit Euler型の安定目安です。
一般の結合系では、更新行列のスペクトル半径を使う条件が優先されます。

## Monte Carloノイズとの関係

ノイズ込みの 1 モード近似として

$$
\delta q^{n+1} \;=\; \left(1 - \frac{\Delta t_b}{\tau}\right)\,\delta q^n \;+\; \xi^n
$$

を考えると、定常分散は$\xi^n$の分散に依存します。$\mathrm{Var}(\xi^n)$が$\Delta t_b$に
どう依存するかは、注入の正規化方式によって変わります。BEACHの`reservoir_face`には2つの方式があります。

### `w_particle`を固定する場合

`w_particle` を直接指定し、物理流入数だけが $\Delta t_b$ に比例して増減する方式では、1 バッチあたりの期待マクロ粒子数は

$$
N_\text{macro} \;\propto\; \Delta t_b
$$

となります。バッチ電荷増分のショットノイズ分散も概ね $\propto \Delta t_b$ とみなせ、

$$
\mathrm{Var}(\xi^n) \;\approx\; \alpha\, \Delta t_b
$$

と置けます。$\Delta t_b \ll \tau$ の極限では定常分散は `batch_duration` に強くは依存しません。

### `target_macro_particles_per_batch`を固定する場合

一方、`target_macro_particles_per_batch`から`w_particle`を求める方式では、
`src/config/bem_app_config_runtime.f90:644`に示すように

$$
w_\text{particle} \;\propto\; \frac{\Gamma\, A\, \Delta t_b}{N_\text{target}}
$$

の形で重みが決まります。マクロ粒子数が固定され、1粒子あたりの寄与が$\Delta t_b$に比例するため、
ノイズの$\Delta t_b$依存性は`w_particle`を固定する場合とは異なります。

### 安定性とノイズを別々に調整する

設定では、次の2つを分けて調整します。

- `batch_duration` は主に **deterministic な安定性** のつまみ
- 統計ノイズの主な制御つまみは **`w_particle` または `target_macro_particles_per_batch`**

`batch_duration`を変えたときにノイズが増減する方向は、注入の正規化方式によって
異なります。

## $\tau_{\min}$を物理スケールから見積もる

$\tau_{\min}$は、数値安定性を支配する最速の有効応答時間です。その値は幾何、電位分布、上流分布関数、
注入モデルに依存します。充電・シース緩和時間と電子プラズマ周波数の逆数を、異なる物理スケールとして
見積もります。

### 充電・シース緩和時間

代表的には、ある有効容量 $C_\text{eff}$ と有効コンダクタンス $G_\text{eff}$ を用いて

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}}{G_\text{eff}}
$$

あるいは典型電位変化 $\Delta\phi$ と有効電流 $I_\text{eff}$ を用いて

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}\,\Delta\phi}{I_\text{eff}}
$$

と見積もるのが自然です。これは幾何や遮蔽の影響を受ける、比較的遅い charging timescale です。

### プラズマ周波数の逆数

もうひとつの速い基準は

$$
\tau_{pe} \;=\; \omega_{pe}^{-1} \;=\; \sqrt{\frac{\varepsilon_0\, m_e}{n_e\, e^2}}
$$

です。これは電子プラズマの微視的な時間スケールであり、系が応答し得る速さの基準になります。

$\omega_{pe}^{-1}$は速い側の物理基準であり、そのまま$\tau_{\min}$の上界にはなりません。
実際に`batch_duration`を制限する有効時定数は、幾何や入射律速を含む$\tau_\text{charge}$側で
決まることがあります。

### 2つの時間スケールを使い分ける

したがって $\tau_{\min}$ については

- $\omega_{pe}^{-1}$ : 微視的な速い基準
- $\tau_\text{charge}$ : 系固有の charging / sheath relaxation timescale

を別々に見積もり、最終的な`batch_duration`は数値実験で決めます。安定性の確認には、
$\tau_\text{charge}$を含む系の有効応答時間が必要です。

## 設定値を決める

1. まず物理スケールとして $\omega_{pe}^{-1}$ と $\tau_\text{charge}$ の両方を概算する。
2. 初期の `batch_duration` は、振動を避けたいなら保守的に小さめから始める。
3. `batch_duration`を1/2倍と2倍に振り、電荷履歴、`last_rel_change`、要素電荷時系列のジッタを比較する。
4. 収束値がほぼ一致し、かつ振動や発散傾向が見えなければ、その `batch_duration` は実用上十分と判断できる。
5. ノイズが大きい場合は、まず `w_particle` または `target_macro_particles_per_batch` を調整する。`batch_duration` の変更だけでノイズを解決しようとしない。

適応進行を使う場合は、手順3を
`max_nonzero_mode_potential_step`とその1/2の比較に置き換えます。`batch_count`を揃えるだけでは終了時刻が
一致しないため、`simulated_time_s`が共通になる範囲を比較します。

この比較は、誤差の冪乗則を仮定するRichardson外挿ではなく、**step-size sensitivity check**です。

## 判断基準

| 項目 | 結論 |
|---|---|
| 定常値 | 平均更新の固定点は `batch_duration` に依存しない |
| 厳密な安定条件 | $\rho(I + \Delta t_b\, M) < 1$ |
| $2\tau$, $\tau$ 条件 | 実負優勢モードを仮定した explicit Euler の近似目安 |
| $\omega_{pe}^{-1}$ の位置づけ | 微視的な速い基準であり、一般にそのまま安定上界ではない |
| ノイズと `batch_duration` | 依存性は注入の正規化方式に依存する |
| ノイズ低減の主手段 | `w_particle` または `target_macro_particles_per_batch` の調整 |
| 確認方法 | `batch_duration` を振った step-size sensitivity check |
| 適応的$k\ne0$進行 | `batch_duration`は最大幅。電位上限を1/2にして収束を確認 |

平均モデルの定常値は`batch_duration`に依存しませんが、安定に到達できる時間幅は
$\rho(I + \Delta t_b M) < 1$で制限されます。$\tau_{\min}$はケースごとに異なるため、物理スケールの
見積もりとstep-size sensitivity checkの両方で決めます。

## 関連文書

- [Fortran パラメータファイル仕様](Parameters.html) — `sim.batch_duration` / `sim.batch_duration_step` の指定方法
- [Fortran 中心ワークフロー](Workflow.html) — バッチループの実行制御
- [計算モデルの全体像](Algorithms.html) — 物理モデルと数値手法の関係

## Code reference

- 注入での使用: `src/particles/bem_injection.f90`（`reservoir_face` / `photo_raycast`）
- batch生成と粒子重みの決定: `src/config/bem_app_config_runtime.f90`
