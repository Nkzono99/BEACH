title: `batch_duration` の理論的背景

Lang: [日本語](BatchDurationTheory.md) | [English](BatchDurationTheory.en.md)

# `batch_duration` の理論的背景

このページは、`batch_duration` が帯電更新の安定性と観測される定常値にどう影響するかを説明します。
値を決める実行手順は [`batch_duration` をどう決めるか](BatchDurationStability.html)を参照してください。

一文で言えば、`batch_duration` は平均帯電モデルの陽的な更新幅です。安定に収束する平均モデルの
固定点は幅に依存しませんが、離散更新の安定性、有限時間誤差、Monte Carlo 誤差は幅に依存します。

## 平均帯電モデル

insulator 要素 $j$ の蓄積電荷を $q_j$、面積を $A_j$、現在の全要素電荷
$\mathbf q$ に対する正味の表面電荷 flux を $J_j(\mathbf q)$ とします。$J_j$ は入射と放出を
符号付きで含みます。平均の帯電速度を

$$
F_j(\mathbf q)=J_j(\mathbf q)A_j
$$

と書くと、連続時間の平均モデルは

$$
\frac{dq_j}{dt}=F_j(\mathbf q)
$$

です。この式は、同じ表面電荷から得られる場に対して粒子統計を十分に平均した応答を表します。

BEACH は 1 batch 中の場を batch 開始時の状態に固定します。そのため、1 batch の更新は概念的に

$$
\mathbf q^{n+1}
=\mathbf q^n+\Delta t_b\,\mathbf F(\mathbf q^n)+\boldsymbol\eta^n
$$

とみなせます。$\Delta t_b$ は `batch_duration`、$\boldsymbol\eta^n$ は有限個の macro 粒子から生じる
Monte Carlo 誤差です。

## 固定点が幅に依存しない条件

平均更新の固定点 $\mathbf q^\ast$ は

$$
\mathbf F(\mathbf q^\ast)=\mathbf 0
$$

を満たします。全要素の面積が正なら、これは $\mathbf J(\mathbf q^\ast)=\mathbf 0$ と同じ条件です。
$\Delta t_b>0$ は固定点の式に現れないため、平均モデルが同じ固定点へ安定に収束する限り、
固定点自体は `batch_duration` に依存しません。

実際の run では、有限標本の $\boldsymbol\eta^n$、有限の物理終了時刻、非線形応答が残ります。
したがって「理論上の固定点が同じ」であっても、観測した最終電荷が時間幅に依存しないとは限りません。

## 固定点近傍の線形安定性

$\mathbf F$ が固定点近傍で微分可能とし、Jacobian を

$$
M=\left.\frac{\partial\mathbf F}{\partial\mathbf q}\right|_{\mathbf q^\ast}
$$

とします。Monte Carlo 誤差を一度無視し、摂動を
$\delta\mathbf q^n=\mathbf q^n-\mathbf q^\ast$ と置くと、線形化した更新は

$$
\delta\mathbf q^{n+1}
\simeq (I+\Delta t_b M)\delta\mathbf q^n
$$

です。したがって一般的な線形安定条件は

$$
\rho(I+\Delta t_b M)<1
$$

です。$\rho$ は spectral radius を表します。

すべての支配固有値が実負で、最速の応答時間を $\tau_{\min}$ で近似できる場合に限り、

$$
\Delta t_b<2\tau_{\min}
$$

が非発散の目安、

$$
\Delta t_b<\tau_{\min}
$$

が単調収束の目安になります。複素固有値、非正規な結合、強い非線形性がある場合、この 1 時間 scale の
規則だけでは安定性を判断できません。これは一般的な BEACH の CFL 条件でもありません。

## 物理的な時間 scale

初期値の候補を考えるとき、電子プラズマ周波数の逆数 $\omega_{pe}^{-1}$ と、等価容量・コンダクタンスから得る
charging time

$$
\tau_\text{charge}\sim\frac{C_\text{eff}}{G_\text{eff}}
$$

は別の時間 scale です。前者は速いプラズマ応答、後者は表面電荷が変わる遅い応答の目安になります。
実際の上限は geometry、電位、流入分布、表面応答にも依存するため、どちらか一方をそのまま
`batch_duration` の上限にはできません。

## Monte Carlo 誤差との結合

$\boldsymbol\eta^n$ の分布は、macro 粒子数と重みに依存します。時間幅を変えたときに粒子重みも変わる
source では、比較結果に時間離散化差と標本分散差の両方が含まれます。

このため、固定幅比較では乱数 seed と並列構成を揃え、表面電荷の norm だけでなく absorbed / escaped 数、
macro 粒子数、重みも記録します。ノイズが支配的なら、線形安定性の結論を出す前に粒子統計を改善します。

## adaptive $k\ne0$ 条件が制御するもの

適応進行は候補電荷と batch 開始電荷の差 $\Delta\mathbf q$ が作る $k\ne0$ 電位について、
全 panel 重心 $\mathbf x_i$ で

$$
\max_i\left|\Delta\phi_{k\ne0}(\mathbf x_i;\Delta\mathbf q)\right|
\le \Phi_{\max}
$$

を満たす幅を受理します。$\Phi_{\max}$ が `max_nonzero_mode_potential_step` です。

この条件は、batch 中に固定した $k\ne0$ 場が 1 回の更新で大きく変わりすぎないための trust bound です。
局所打切り誤差の推定器ではなく、大域精度の次数も与えません。$k=0$ 更新や、応答表・粒子 sampling の
適用範囲も制御しません。

## 前提と適用限界

- 平均モデルは、同じ $\mathbf q$ に対する粒子応答を統計平均できると仮定します。
- 線形条件は、調べる固定点の近傍で $\mathbf F$ が微分可能である場合に限ります。
- 複数の固定点、hysteresis、非定常な attractor がある場合、固定点解析だけでは最終状態を予測できません。
- $\boldsymbol\eta^n$ を無視した安定条件は、有限 macro 粒子 run の揺らぎ幅を保証しません。
- adaptive $k\ne0$ 条件は、matching-plane の内側反復や implicit zero-mode 更新を含む全系の安定条件ではありません。
- この説明は v1.0 の insulator accumulation を対象とし、未実装の resistive / dielectric 応答へ拡張しません。

したがって最終的な値は、理論 scale だけでなく、同じ物理時刻での step-size sensitivity check で決めます。
実行完了、数値収束、物理的妥当性は別々に判定してください。

## 関連文書

- [`batch_duration` をどう決めるか](BatchDurationStability.html) — 固定幅比較と adaptive 設定
- [BEACH の計算サイクル](Algorithms.html) — batch 中に場を固定する更新順序
- [入力パラメータリファレンス](Parameters.html) — 設定契約
- [結果を検証する](ValidationGuide.html) — 数値収束と物理妥当性
