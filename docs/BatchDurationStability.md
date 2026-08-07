title: `batch_duration`の安定性と定常値

Lang: [日本語](BatchDurationStability.md) | [English](BatchDurationStability.en.md)

# `batch_duration`の安定性と定常値

`sim.batch_duration`を選び、帯電結果の時間幅依存性を確認するためのタスクガイドです。
固定幅では1 batchの物理時間と表面電荷更新幅を表します。
`sim.batch_duration_step`を指定した場合は、`sim.batch_duration = sim.dt * sim.batch_duration_step`です。

> `sim.tol_rel`は監視・出力値です。現行実装は`tol_rel`を早期停止条件に使わず、
> `sim.batch_count`で指定したaccepted batch数まで実行します。

## 固定幅を選ぶ

### 前提

- `boundary_inflow`、`plane_source`、`reservoir_face`、`photo_raycast`には正の`sim.batch_duration`が必要です。
- 比較用に`history_stride > 0`を設定し、`charge_history.csv`を保存します。
- 同じメッシュ、粒子分布、乱数seed、OpenMP/MPI構成を使います。

### 操作

1. 保守的に小さい`batch_duration`または`batch_duration_step`から始めます。
2. 基準値の1/2倍、1倍、2倍で3 runを実行します。
3. 3 runを同じ物理時刻付近で比較します。
4. `summary.txt`の`last_rel_change`、総電荷、absorbed/escaped数と、
   `charge_history.csv`の要素電荷を確認します。

### 期待する出力

各runに`summary.txt`と`charge_history.csv`が作成され、次を比較できます。

- 最終表面電荷分布
- 総電荷と局所電位幅
- 電荷履歴の振動、発散、Monte Carloジッタ
- `simulated_time_s`とaccepted batch数

### 解釈

| 観測 | 判断 |
| --- | --- |
| 1/2倍と1倍で最終電荷と履歴がほぼ一致 | 1倍は実用上十分な候補 |
| 時間幅を増やすと振動または発散 | `batch_duration`を下げる |
| 時間幅で最終電荷が大きく変化 | より小さい時間幅で再計算 |
| 履歴がノイズに埋もれる | `w_particle`または`target_macro_particles_per_batch`を調整 |
| 終了時にも変化が続く | `batch_count`を増やす |

正常終了だけでは安定性や定常到達を示しません。この比較は誤差の冪乗則を仮定する
Richardson外挿ではなく、step-size sensitivity checkです。

### 次の選択

- 1倍と1/2倍が一致する: 小さい方を検証基準とし、コストを優先する場合は1倍を採用します。
- deterministicな振動が残る: `batch_duration`を下げます。
- ノイズが支配的: 時間幅ではなくmacro粒子数または重みを先に調整します。
- `cached_kneq0`で1 batchの局所電位変化を制限する: 次の適応進行を使います。

## 適応的な$k\ne0$進行を使う

### 前提

この経路は次を要求します。

- `[periodic2].nonzero_mode_backend = "cached_kneq0"`
- time-scaledな`boundary_inflow`、`plane_source`、`reservoir_face`、`photo_raycast`
- reservoir流入と面sourceでは固定`w_particle`ではなく`target_macro_particles_per_batch`
- 正の`sim.batch_duration`

`volume_seed`はこの経路では使えません。

### 操作

```toml
[periodic2]
nonzero_mode_backend = "cached_kneq0"
max_nonzero_mode_potential_step = 1.0e-2 # V
```

解決済み`sim.batch_duration`を$h_0$とすると、各accepted batchで
$h_0,h_0/2,h_0/4,\ldots$を順に試します。BEACHはcandidate電荷とbatch開始電荷の差が作る
$k\ne0$電位を全panel重心で評価し、最大絶対値が上限以下となる最初のtrialを受理します。

### 期待する出力

- 棄却trialはRNGとmacro粒子残差をrollbackし、統計、履歴、charge ledgerへ記録しません。
- `summary.txt`の`simulated_time_s`が実際の物理終了時刻です。
- `batch_count`はaccepted batch数であり、物理時間は一般に
  `batch_count * batch_duration`ではありません。

同じ実行内の棄却trialは固定OpenMP team sizeで再現します。restart後は別のteam sizeを使用できるため、
restart前後のbitwise一致ではなく、電荷分布とaccepted widthの数値的一致を確認します。

### 解釈

`max_nonzero_mode_potential_step`は、batch内で$k\ne0$場を凍結するための局所電位trust boundです。
局所打切り誤差や大域精度の次数は保証しません。また、$k=0$更新の安定性も制御しません。

target-count reservoirと固定`rays_per_batch`のphoto sourceでは、trial幅を半分にするとmacro粒子電荷も半分になります。
上限半減比較には時間離散化差とMonte Carlo分散差の両方が含まれるため、同じ乱数seedを使い、
電荷分布normと粒子統計を併記します。

### 次の選択

1. `max_nonzero_mode_potential_step`を1/2にします。
2. 同じ`simulated_time_s`付近で表面電荷、局所電位幅、総電荷、粒子統計を比較します。
3. 固定幅controlにはこのキーを省略するか`0`を指定します。

## 判断に必要な理論

insulator要素$j$の蓄積電荷を$q_j$、面積を$A_j$、入射電荷fluxを$J_j(\mathbf q)$とすると、
平均モデルは

$$
\frac{dq_j}{dt}=J_j(\mathbf q)A_j
$$

です。1 batchはbatch開始時の場を凍結するexplicit updateとして

$$
\mathbf q^{n+1}
=\mathbf q^n+\Delta t_b\,\mathbf J(\mathbf q^n)\mathbf A+\boldsymbol\eta^n
$$

とみなせます。$\boldsymbol\eta^n$はMonte Carlo誤差です。

平均更新の固定点$\mathbf q^\ast$は$\mathbf J(\mathbf q^\ast)=0$を満たすため、
平均モデルが安定に収束する限り固定点自体は$\Delta t_b$に依存しません。ただし実runには有限標本誤差と
有限時間での停止誤差があるため、観測結果の時間幅依存性は残ります。

固定点近傍の一般的な線形安定条件は

$$
\rho(I+\Delta t_b M)<1
$$

です。支配固有値が実負で、最速応答時間を$\tau_{\min}$と近似できる場合に限り、
$\Delta t_b<2\tau_{\min}$が非発散、$\Delta t_b<\tau_{\min}$が単調収束の目安です。
これは一般的なBEACHのCFL条件ではありません。

物理scaleとして$\omega_{pe}^{-1}$とcharging time
$\tau_\text{charge}\sim C_\text{eff}/G_\text{eff}$を別々に見積もれますが、
実際の上限は幾何、電位、流入分布にも依存します。最終値はstep-size sensitivity checkで決めてください。

## 関連文書

- [入力パラメータリファレンス](Parameters.html) — `sim.batch_duration`と`sim.batch_duration_step`
- [境界reservoirの流入量と速度サンプリング](ReservoirInjection.html) — 粒子数と重み
- [結果を検証する](ValidationGuide.html) — 数値収束と物理妥当性
- [計算モデルの全体像](Algorithms.html) — batch loop
