title: Zhao stationary closure を理解する

Lang: [日本語](ZhaoStationaryClosure.md) | [English](ZhaoStationaryClosure.en.md)

# Zhao stationary closure を理解する

`surface_current_model.model="zhao_stationary"` は、BEACH の計算領域外に定常シースを仮定し、
その電流と速度空間の障壁を box 上端へ接続します。外部シースを時間発展させる model ではありません。
run 開始時に一度だけ零電流根を解き、固定した総電流を BEACH 内で追跡した hit・return 分布へ割り当てます。

このページでは、零電流根、0 V reservoir からの流入、光電子（PE）の return / escape がどうつながるかを説明します。
全入力キーと制約は[入力パラメータ](Parameters.html)、出力名は
[出力形式リファレンス](OutputReference.html#zhao_stationary)を正本とします。

## 起動時に解く問題

Zhao stationary closure は平面、無衝突、非磁化の外部シースを仮定します。ambient electron、cold ion、
PE を有効にした場合は photoelectron を用いて、Zhao Type A / B / C の零電流定常根を解きます。
解から branch、表面電位 $\phi_0$、電位極小 $\phi_m$、ambient electron 密度、species 別電流密度を得ます。

零電流根は run 中に解き直しません。各 batch の表面電荷が変わっても、branch、$\phi_0$、$\phi_m$、
電流 target は同じです。各 batch の外向き flux と外部応答を反復する
[matching-plane 準定常連成](MatchingPlaneCoupling.html)とは、この点が異なります。

### PE なしは Type C

`photoelectron_source_scale=0.0` では PE species と PE 固有入力を使いません。`zhao_branch="auto"` または
`"c"` は Type C として

$$
J_e + J_i = 0
$$

を解きます。生成されるのは electron / ion の吸収 target と z-high の kinetic map だけです。
PE emission、return、escape target と PE 粒子は生成されません。

### PE ありは大きな channel を別々に閉じる

PE source density は太陽高度 $\alpha$、基準密度 $n_{pe,ref}$、scale $s_{UV}$ から

$$
n_{pe,0}=s_{UV}n_{pe,ref}\sin\alpha
$$

と定めます。ion species の数密度は無限遠 ion 密度です。electron の設定密度と PE species の
`emit_current_density_a_m2` は raw 粒子分布の標本化に使い、Zhao 根が固定電流 target を決めます。
電流の符号は表面帯電への寄与で表します。

| channel | 符号 | BEACH での扱い |
|---|---:|---|
| ambient electron 吸収 $J_e$ | 負 | 表面への吸収 target |
| ion 吸収 $J_i$ | 正 | 表面への吸収 target |
| PE emission $J_{emit}$ | 正 | 放出反作用 target |
| PE return $J_{return}$ | 0 以下 | 表面への再吸収 target |
| PE escape $J_{escape}$ | 正 | 表面に deposit しない外部境界 target |

Zhao 根は次の二つの収支を与えます。

$$
J_{return}=J_{escape}-J_{emit}\le0,
\qquad
J_e+J_i+J_{escape}=0.
$$

したがって、PE channel と表面 channel はそれぞれ

$$
J_{return}+J_{emit}-J_{escape}=0,
\qquad
J_e+J_i+J_{return}+J_{emit}=0
$$

で閉じます。BEACH は emission と return を別々に倍率補正し、小さな net PE 電流を倍率の分母にしません。
PE 粒子が外へ運ぶ signed current は $-A J_{escape}$ で、表面には加算しません。

## 0 V reservoir を現在の box 上端へ写像する

ambient の設定 Maxwell VDF は、無限遠のプラズマ電位を 0 V とした reservoir 分布です。
各 batch の開始時に、BEACH は z-high 面の平均電位 $\phi_f$ を現在の表面電荷から評価します。
流入粒子は外部 access bottleneck と $\phi_f$ の両方へ到達できる reservoir tail から選ばれ、法線速度は

$$
\frac12 m v_{n,f}^{2}=\frac12 m v_{n,\infty}^{2}-q(\phi_f-0)
$$

となるようエネルギー保存で注入面へ写像されます。この面平均近似が不十分でないかは、z-high 面内の
電位ばらつきも併せて確認します。

| branch | ambient electron の access | electron / PE の外向き barrier | ion の access / barrier |
|---|---:|---:|---:|
| Type A | $\phi_m$ | $\phi_m$ | 0 V |
| Type B / C | 0 V | 0 V | 0 V |

粒子が z-high を外向きに横切るときは、面平均ではなく横切り位置の局所電位と法線運動エネルギーを使います。
固定 barrier へ到達できない electron / PE は z-high で鏡面反射し、到達できる粒子だけを escape と分類します。
PE の放出速度そのものは、PE species に設定した表面 half-Maxwellian のままです。

## 固定される量と batch ごとに変わる量

| run 開始時に固定 | 各 batch で再評価・再標本化 |
|---|---|
| branch、$\phi_0$、$\phi_m$、ambient electron 密度 | 表面電荷と BEACH 領域内の場 |
| species 別の signed 電流密度と reference area | z-high 面平均電位 $\phi_f$ と流入 tail |
| 吸収・放出・escape の電流 target | 粒子軌道、hit 位置、局所 return / escape 分類 |
| branch 別の access / barrier 電位 | $I\Delta t$ の target 電荷と raw 分布への倍率 |

この分離により総電流は定常根へ合わせられますが、要素別の帯電分布は各 batch の Monte Carlo 軌道に依存します。

## 適用範囲

- box 外の電場、空間電荷、Debye shielding、turning point までの距離や飛行時間は解きません。
- z-high 反射は外部 return 軌道を境界へ縮約した断熱的な近似です。
- 非磁化 closure なので一様磁場はゼロでなければなりません。
- 平面定常解なので、曲率、衝突、外部過渡、run 中に変化する照射や plasma 条件は自己整合には応答しません。
- 零電流 residual が小さくても、有限粒子で得た hit / return の空間分布の収束は保証されません。

species、境界、電荷、温度の全入力条件は[入力パラメータ](Parameters.html)で確認してください。
raycast の放出反作用と VDF は[光電子の放出とreturn](PhotoelectronEmission.html)で説明しています。

## 例を実行する

PE ありの完全な例は `examples/periodic2_zhao_fixed_current.toml` です。

```bash
beach examples/periodic2_zhao_fixed_current.toml
```

PE なしの Type C は `examples/periodic2_zhao_no_photo_fixed_current.toml` で確認できます。
どちらも `[surface_current_model]` の設定だけでなく、必要な reservoir、open z-high、fixed-current species を含みます。

完了後、`summary.txt` で model と選択 branch、二つの budget residual を確認します。次に
`charge_ledger.csv` で raw / target / applied 電荷、hit 数、`fixed_*_weight_scale` を比較します。
PE ありでは ray 数、batch 幅、乱数 seed を変え、要素別 return 分布が収束することを確認してください。
実行完了と零電流収束だけでは、外部シース近似や空間分布の物理妥当性は確定しません。

実装上の定義は [SPEC 7.7](../SPEC.md#77-自動表面電流model) と
[`bem_surface_current_model.f90`](../src/physics/sheath/bem_surface_current_model.f90) にあります。
