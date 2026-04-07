# `batch_duration` の安定性と定常値の理論

この文書は、BEACH のバッチループにおける `sim.batch_duration`（または `sim.batch_duration_step` から決まる 1 バッチの物理時間）と、収束した壁面電荷分布の正当性・安定性の関係を理論的に整理したものです。現行実装では、`batch_duration_step` を使う場合は `sim.batch_duration = sim.dt * sim.batch_duration_step` と解釈され、`reservoir_face` 注入では 1 バッチの物理流入量からマクロ粒子数または重みが決まります。  

実装上の入口は次のとおりです。

* バッチ計算手順: `SPEC.md`
* パラメータ定義: `sim.batch_duration` / `sim.batch_duration_step`
* 注入での使用: `src/particles/bem_injection.f90`
* バッチ生成・重み解決: `src/config/bem_app_config_runtime.f90`

## 1. 連続時間モデルへの還元

絶縁体壁面の各要素 `j` の蓄積電荷を (q_j(t))、そのときの入射電荷フラックスを (J_j(\mathbf q)) とすると、吸収のみを考える基本モデルは

$$
\frac{dq_j}{dt} = J_j(\mathbf q),A_j
$$

と書けます。ここで (A_j) は要素面積です。(J) は壁面電荷が作る電場に依存するため、一般には非線形です。

BEACH の 1 バッチ更新は、この連続時間モデルを **場をバッチ先頭で凍結した explicit な更新**として評価するものとみなせます。平均的には

$$
\mathbf q^{n+1}
===============

\mathbf q^n
+
\Delta t_b,\mathbf J(\mathbf q^n)\mathbf A
+
\boldsymbol\eta^n
$$

と表せます。ここで

* (\Delta t_b =) `sim.batch_duration`
* (\mathbf A) は要素面積ベクトル
* (\boldsymbol\eta^n) はバッチ内の Monte Carlo サンプリング誤差

です。実装上も、バッチ内で粒子は同じ場を見て進み、バッチ末尾で電荷差分がまとめて壁面へ反映されます。したがって、`batch_duration` はこの explicit 更新の時間刻みに相当します。 

## 2. 定常値の正当性

平均更新写像を

$$
F_{\Delta t_b}(\mathbf q)
=========================

\mathbf q + \Delta t_b,\mathbf J(\mathbf q)\mathbf A
$$

と書くと、その不動点 (\mathbf q^*) は

$$
F_{\Delta t_b}(\mathbf q^*)=\mathbf q^*
\quad\Longleftrightarrow\quad
\mathbf J(\mathbf q^*)=0
$$

で与えられます。したがって、**平均モデルの不動点自体は (\Delta t_b) に依存しません**。

この意味で、反復が安定に収束しており、かつ Monte Carlo 誤差が十分に平均化されているなら、`batch_duration` を変えても目指している連続時間の定常解は同じです。

ただし、ここで言えるのはあくまで **平均モデルの固定点** についてです。実際の計算では

* バッチごとの有限サンプル誤差
* 収束判定に用いる監視量の揺らぎ
* 有限バッチ数で打ち切ることによる残差

が入るため、観測される収束値には弱い `batch_duration` 依存が残りえます。したがって、厳密には

> 反復の平均的な固定点は `batch_duration` に依存しないが、有限サンプル・有限時間の実計算では小さな step-size dependence が観測されうる

と書くのが安全です。

## 3. 線形安定性

不動点 (\mathbf q^*) 近傍で摂動 (\delta\mathbf q^n=\mathbf q^n-\mathbf q^*) を考えると、平均更新の線形化は

$$
\delta \mathbf q^{n+1}
======================

\left(I+\Delta t_b M\right)\delta \mathbf q^n,
\qquad
M_{ij}
======

\frac{\partial (J_iA_i)}{\partial q_j}\bigg|_{\mathbf q^*}
$$

となります。

このとき一般の多自由度系での安定条件は

$$
\rho!\left(I+\Delta t_b M\right)<1
$$

です。すなわち、各固有値 (\lambda_k) に対して

$$
|1+\Delta t_b \lambda_k|<1
$$

が必要です。

したがって、安定条件として本質的なのは **explicit Euler の安定条件** であり、これを単純な scalar 問題のように

$$
\Delta t_b < 2\tau_{\min}
$$

と書けるのは、**支配的な固有値が実負で近似できる場合**です。もし主要モードが実負であり、

$$
\tau_k \equiv \frac{1}{|\lambda_k|}
$$

と置けるなら、最速モードに対して

$$
0 < \Delta t_b < \frac{2}{|\lambda_{\max}|}=2\tau_{\min}
$$

で発散を避けられ、

$$
0 < \Delta t_b < \frac{1}{|\lambda_{\max}|}=\tau_{\min}
$$

で単調収束になります。

したがって、実務上は次のように読むのが適切です。

* (\Delta t_b < 2\tau_{\min}): 実負優勢モード仮定のもとでの非発散条件
* (\Delta t_b < \tau_{\min}): 同仮定のもとでの単調収束条件
* 一般の結合系では、厳密には (\rho(I+\Delta t_b M)<1) が本体

この意味で、元の (2\tau)・(\tau) 条件は「BEACH 版の CFL 条件」と断言するより、**実負優勢モードを仮定した explicit Euler 型の安定目安**と呼ぶ方が正確です。

## 4. Monte Carlo ノイズとの関係

ノイズ込みの 1 モード近似として

$$
\delta q^{n+1}
==============

\left(1-\frac{\Delta t_b}{\tau}\right)\delta q^n + \xi^n
$$

を考えると、定常分散は (\xi^n) の分散に依存します。ただし、ここで重要なのは **(\mathrm{Var}(\xi^n)) の (\Delta t_b) 依存が実装の正規化方式に依存する** ことです。

### 4.1 `w_particle` 固定の場合

`reservoir_face` で `w_particle` を固定し、物理流入数だけが (\Delta t_b) に比例して増減するなら、1 バッチあたりの期待マクロ粒子数は概ね

$$
N_{\rm macro}\propto \Delta t_b
$$

です。このときバッチ電荷増分のショットノイズ分散も概ね (\propto \Delta t_b) とみなせます。

この仮定のもとでは

$$
\mathrm{Var}(\xi^n)\approx \alpha,\Delta t_b
$$

と置け、さらに (\Delta t_b\ll\tau) の極限では定常分散は `batch_duration` に強くは依存しません。

### 4.2 `target_macro_particles_per_batch` 固定の場合

一方、現行の BEACH では `reservoir_face` に対して `w_particle` を直接与える方法と、`target_macro_particles_per_batch` から `w_particle` を解く方法の両方があります。後者では runtime 側で

$$
w_{\rm particle}
\propto
\frac{\Gamma A \Delta t_b}{N_{\rm target}}
$$

という形で重みが決まるため、ノイズの (\Delta t_b) 依存は上の単純な (\mathrm{Var}(\xi^n)\propto\Delta t_b) とは変わります。  

### 4.3 実務上の解釈

したがって、理論整理としては

* `batch_duration` は主に deterministic な安定性のつまみ
* 統計ノイズの主な制御つまみは `w_particle` または `target_macro_particles_per_batch`

と考えるのがよいです。

特に、

> `batch_duration` を小さくすれば必ずノイズが下がる
> `batch_duration` を大きくしてもノイズはほとんど変わらない

のいずれも一般論としては言い切れません。そこは注入の正規化方式に依存します。

## 5. (\tau_{\min}) の物理的見積もり

(\tau_{\min}) は、数値安定性を支配する最速の有効応答時間です。ただし、これを 1 つの物理式で一般に与えるのは難しく、幾何・電位分布・上流分布関数・注入モデルに依存します。実務上は、次の 2 つを分けて考えるのが安全です。

### 5.1 充電・シース緩和時間

代表的には、ある有効容量 (C_{\rm eff}) と有効コンダクタンス (G_{\rm eff}) を用いて

$$
\tau_{\rm charge}\sim \frac{C_{\rm eff}}{G_{\rm eff}}
$$

あるいは、典型電位変化 (\Delta\phi) と有効電流 (I_{\rm eff}) を用いて

$$
\tau_{\rm charge}\sim \frac{C_{\rm eff},\Delta\phi}{I_{\rm eff}}
$$

と見積もるのが自然です。これは幾何や遮蔽の影響を受ける、比較的遅い charging timescale です。

### 5.2 プラズマ周波数の逆数

もうひとつの速い基準は

$$
\tau_{pe}=\omega_{pe}^{-1}
==========================

\sqrt{\frac{\varepsilon_0 m_e}{n_e e^2}}
$$

です。これは電子プラズマの微視的な速い時間スケールであり、系がどこまで急峻に応答しうるかを見る基準にはなります。

ただし、(\omega_{pe}^{-1}) をそのまま (\tau_{\min}) の上界とみなすのは強すぎます。むしろこれは **速い側の物理基準** であり、実際に `batch_duration` を制限する有効時定数は、幾何や入射律速を含んだ (\tau_{\rm charge}) 側で決まることが少なくありません。

### 5.3 実用上の選び方

したがって、(\tau_{\min}) については

* (\omega_{pe}^{-1}): 微視的な速い基準
* (\tau_{\rm charge}): 系固有の charging / sheath relaxation timescale

を別物として見積もり、最終的には数値実験で詰めるのが安全です。

元の文書のように

$$
\tau_{\min}\approx \min(\tau_{\rm sheath},\tau_{pe})
$$

と断定するより、

> (\omega_{pe}^{-1}) は速い基準、実際の安定制約は (\tau_{\rm charge}) を含む有効応答時間で決まる

と書く方が理論的に無理がありません。

## 6. 実用的な使い方

1. まず物理スケールとして (\omega_{pe}^{-1}) と (\tau_{\rm charge}) の両方を概算する。
2. 初期の `batch_duration` は、振動を避けたいなら保守的に小さめから始める。
3. `batch_duration` を 1/2 倍と 2 倍に振って、電荷履歴や監視量の振る舞いを比較する。
4. 収束値がほぼ一致し、かつ振動や発散傾向が見えなければ、その `batch_duration` は実用上十分と判断できる。
5. ノイズが大きい場合は、まず `w_particle` または `target_macro_particles_per_batch` を調整する。`batch_duration` の変更だけでノイズを解決しようとしない。
6. `last_rel_change` の振動や要素電荷時系列のジッタは有用な診断量だが、これは厳密な Richardson 外挿というより **step-size sensitivity check** と呼ぶのが適切である。

## 7. まとめ

| 項目                       | 修正版の結論                                                  |
| ------------------------ | ------------------------------------------------------- |
| 定常値の正当性                  | 平均更新の固定点は `batch_duration` に依存しない                       |
| 厳密な安定条件                  | (\rho(I+\Delta t_b M)<1)                                |
| (2\tau), (\tau) 条件       | 実負優勢モードを仮定した explicit Euler の近似目安                       |
| (\omega_{pe}^{-1}) の位置づけ | 微視的な速い基準であり、一般にそのまま安定上界ではない                             |
| ノイズと `batch_duration`    | 依存性は注入の正規化方式に依存する                                       |
| ノイズ低減の主手段                | `w_particle` または `target_macro_particles_per_batch` の調整 |
| 実務上の確認方法                 | `batch_duration` を振った step-size sensitivity check       |
