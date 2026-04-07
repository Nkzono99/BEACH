title: batch_duration の安定性と定常値の理論

# `batch_duration` の安定性と定常値の理論

この文書は、BEACH のバッチループにおける `sim.batch_duration`（= 1 バッチで反映する物理時間）と、
収束した壁面電荷分布の正当性・安定性の関係を理論的にまとめたものです。

実装上の入口は次のとおりです。

- バッチ計算手順: [`SPEC.md` §4](https://github.com/Nkzono99/BEACH/blob/main/SPEC.md)
- パラメータ定義: [`docs/fortran_parameter_file.md`](fortran_parameter_file.html) の `sim.batch_duration` / `sim.batch_duration_step`
- 注入での使用: `src/particles/bem_injection.f90`（`reservoir_face` / `photo_raycast`）
- 1 バッチの場凍結: `src/core/bem_types.f90` と `src/config/bem_app_config_runtime.f90`

## 1. 連続時間モデルへの還元

絶縁体壁面の各要素 `j` の蓄積電荷 $q_j(t)$ と入射電荷フラックス（壁単位面積あたり）$J_j(\mathbf q)$ を用いると、
吸収のみ（v1.0 の正規仕様）の動力学は次の ODE に書けます。

$$
\frac{dq_j}{dt} \;=\; J_j(\mathbf q)\, A_j
$$

$J$ は壁電荷を介した電場分布に依存するので **非線形** です（$J = J(\mathbf q)$）。

BEACH の 1 バッチ更新は、この ODE に対する **前進 Euler 法** を **Monte Carlo サンプリング**で評価したものに一致します。

$$
\mathbf q^{n+1} \;=\; \mathbf q^n \;+\; \Delta t_b \cdot \widehat{\mathbf J}(\mathbf q^n)\,\mathbf A
$$

- $\Delta t_b = $ `sim.batch_duration`
- $\widehat{\mathbf J}$ はバッチ内 Monte Carlo サンプリングによる $J$ の不偏推定
- 場 $E$ をバッチ先頭で凍結 → バッチ内の粒子はすべて同じ $E(\mathbf q^n)$ を見る

これが「離散的に反映」の正体です。

## 2. 定常値の正当性（バイアス）

不動点 $\mathbf q^\*$ は

$$
\widehat{\mathbf J}(\mathbf q^\*) \;=\; 0
$$

で定義されます。期待値レベルでは $\mathbf J(\mathbf q^\*) = 0$ は $\Delta t_b$ に依存しません。すなわち

> **離散反復が収束した値は、$\Delta t_b$ の選び方によらず連続時間の定常解に一致する**（バイアスは生じない）

これが「正当性」の理論的な根拠です。反復が安定である限り常に成り立ちます。

注意点:

- $\Delta t_b$ が大きすぎて発散すると当然この性質は無意味になります（§3）。
- Monte Carlo ノイズに起因する確率的バイアスは通常 $O(\sigma^2)$ で、入射粒子が線形的にサンプリングされる絶縁体モデルでは 0 次バイアスは無いと見なせます（§4）。

## 3. 線形安定性（CFL 型条件）

不動点近傍の摂動 $\delta \mathbf q^n$ について、

$$
\delta \mathbf q^{n+1} \;=\; \bigl(I + \Delta t_b\, M\bigr)\, \delta \mathbf q^n,
\qquad
M_{ij} \;\equiv\; \frac{\partial (J_i A_i)}{\partial q_j}\bigg|_{\mathbf q^\*}
$$

絶縁体壁が電荷を貯めるとそれ以上同種粒子を引き寄せにくくなる（負のフィードバック）ため、
$M$ の固有値 $\lambda_k$ はおおむね $\mathrm{Re}(\lambda_k) < 0$ です。
応答時間スケールを $\tau_k \equiv 1/|\lambda_k|$ と書くと、最も速いモードに対する **安定条件** は次のとおりです。

$$
\boxed{\;\Delta t_b \;<\; \frac{2}{|\lambda_{\max}|} \;=\; 2\,\tau_{\min}\;}
\quad\text{(発散しない条件)}
$$

$$
\boxed{\;\Delta t_b \;<\; \tau_{\min}\;}
\quad\text{(振動せず単調収束する条件)}
$$

これが BEACH 版の **CFL 条件** です。前進 Euler なので

- $\Delta t_b > 2\,\tau_{\min}$ : 発散
- $\tau_{\min} < \Delta t_b < 2\,\tau_{\min}$ : 振動収束
- $\Delta t_b < \tau_{\min}$ : 過減衰（単調収束）

となります。

## 4. Monte Carlo ノイズとのトレードオフ

ノイズも入れた線形 SDE 描像で、

$$
\delta q^{n+1} \;=\; \bigl(1 - \Delta t_b/\tau\bigr)\,\delta q^n + \xi^n,
\qquad
\mathrm{Var}(\xi^n) \;\approx\; \alpha\, \Delta t_b
$$

（バッチあたりの入射数 $\propto \Delta t_b$ なのでショットノイズ分散も $\propto \Delta t_b$）と書けます。
定常分散は

$$
\sigma_q^2 \;\approx\; \frac{\alpha\,\Delta t_b}{1 - (1-\Delta t_b/\tau)^2}
\;\xrightarrow{\Delta t_b \ll \tau}\;
\frac{\alpha\,\tau}{2}
$$

つまり **小さい $\Delta t_b$ ではノイズはほぼ $\Delta t_b$ に依らず一定** で、$\Delta t_b$ を大きくしても短期的なゆらぎは減りません。
逆に大きすぎると安定限界に当たって発散します。

ノイズを下げる正攻法は `target_macro_particles_per_batch` を増やす（= $w_\text{particle}$ を下げる）ことで、
$\alpha \propto w^2$ なので $\sigma_q^2 \propto w^2$ がそのまま効きます。
**`batch_duration` は安定性側、`w_particle` はノイズ側** と分離して考えるのが理論的に整理しやすいです。

## 5. $\tau_{\min}$ の物理的見積もり

問題依存ですが、プラズマ吸収壁では次のいずれかで上から押さえられます。

### 5.1 シース充電時間

壁電位 $\phi$ が熱電位 $T_e/e$ 程度になるまでに必要な蓄積電荷を、入射電子フラックスで割ったもの:

$$
\tau_\text{sheath} \;\sim\; \frac{C\,T_e/e}{e\, \Gamma_e\, A}
\;\sim\; \frac{\varepsilon_0\, T_e}{e^2\, n_e\, v_{th,e}\, L}
$$

ここで $C$ は注目要素群の自己キャパシタンス、$L$ は系の特徴長、$\Gamma_e$ は電子フラックス。

### 5.2 プラズマ周波数の逆数

$$
\tau_\text{pe} \;=\; \omega_{pe}^{-1}
\;=\; \sqrt{\frac{\varepsilon_0\, m_e}{n_e\, e^2}}
$$

シース応答はこれより速くはなり得ない物理下限です。

### 5.3 実用上の選び方

実用上は $\tau_\text{sheath}$ の方が長いことが多く、$\tau_{\min} \approx \min(\tau_\text{sheath},\, \tau_\text{pe})$ を上限の見積もりに使います。

具体例: `examples/beach.toml` 系列（$n_e \sim 5\,\mathrm{cm^{-3}}$, $T_e = 10\,\mathrm{eV}$, $v_{th,e} \sim 1.3\times 10^6\,\mathrm{m/s}$）では
$\omega_{pe}^{-1} \sim 8\,\mu\mathrm{s}$、`batch_duration_step = 60000` × `dt = 2\times10^{-8}` $= 1.2\,\mathrm{ms}$ なので、

$$
\Delta t_b \;\approx\; 1.2\,\mathrm{ms} \;\gg\; \omega_{pe}^{-1}
$$

つまりプラズマ振動の時間スケールはとっくに超えており、決め手は $\tau_\text{sheath}$ になります。
このオーダーは ms 程度であることが多く、現状の `examples/beach.toml` は安定限界の少し下を狙う設定になっています。

## 6. 実用的な使い方（理論からの処方箋）

1. **$\tau_{\min}$ を 2 通りで見積もる**: シース充電時間と $\omega_{pe}^{-1}$。安全側に短い方を使う。
2. **初期 $\Delta t_b$ は $\tau_{\min}/3$ ～ $\tau_{\min}/2$ 程度**: 過減衰側で振動を避ける。
3. **`batch_duration` を 2 倍 / 1/2 倍した実験を回す**:
   `charge_history.csv` の `last_rel_change` または要素電荷時系列に振動が出始めたら不安定の兆候。
   両条件の収束値を比較し、一致すれば $\Delta t_b$ に関する系統誤差は無視できる（**Richardson 比較**）。
4. **ノイズが大きいときは `target_macro_particles_per_batch` を増やす**: `batch_duration` で対処しない（§4 の通り効果が薄い）。
5. **`tol_rel`（= `last_rel_change`）の振動振幅** は $\sigma_q/q$ の指標。
   $\Delta t_b$ を変えたときに振動が「収まる方向」か「増える方向」かで $\Delta t_b < \tau$ か $\Delta t_b > \tau$ かを大まかに判定できる。

## 7. まとめ

| 項目 | 理論的な答え |
|---|---|
| 定常値の正当性 | $\Delta t_b$ に依らず連続時間定常解と一致（不動点の定義から） |
| 発散しない上限 | $\Delta t_b < 2\,\tau_{\min}$ |
| 単調収束の上限 | $\Delta t_b < \tau_{\min}$ |
| $\tau_{\min}$ の物理 | $\min\!\bigl(\tau_\text{sheath},\;\omega_{pe}^{-1}\bigr)$ で上から見積もる |
| ノイズ低減への効き | 弱い（$\Delta t_b$ では下がらない）。`target_macro_particles_per_batch` で対処 |
| Richardson 検証 | $\Delta t_b$ を 1/2 と 2 倍で回し、定常電荷分布が一致すれば OK |

「定常値が $\Delta t_b$ の取り方に依らない」ことは理論的にきれいに言えます。
「発散しない上限」も古典的な安定性解析でほぼ言えます。
残る不確かさは $\tau_{\min}$ の値そのもので、ここだけはケース依存の物理見積もりか Richardson 数値比較で詰める必要があります。

## 関連文書

- [Fortran パラメータファイル仕様](fortran_parameter_file.html) — `sim.batch_duration` / `sim.batch_duration_step` の指定方法
- [Fortran 中心ワークフロー](fortran_workflow.html) — バッチループの実行制御
- `SPEC.md` — 1 バッチの計算手順と停止条件
