# ADR 0001: 単調 1D kinetic outer-plasmaモデル

- Status: accepted for Phase 7
- Date: 2026-07-11

> `kinetic_1d` に関する決定は引き続き有効である。legacy Zhao injection を独立維持する記述は
> [ADR 0009](0009-remove-linear-debye-and-legacy-sheath.md)で置き換えた。

## Context

BEACH の `linear_debye` outer model は小振幅応答だけを表す。強い表面電位と
photoelectron 空間電荷を扱うには、interface から無限遠までの Poisson 問題と粒子軌道を
同じ電位 profile で閉じる必要がある。一方、任意 VDF、衝突、磁化、非単調 potential、
trapped population を同時に導入すると、境界データだけでは解が一意に定まらない。

## Decision

Phase 7 の `kinetic_1d` は、静電・無衝突・非磁化・平面平均された単調 profile に限定する。
座標は interface を `z=0`、outer plasma を `z>0` とし、`phi(infinity)=0` を gauge とする。
surface charge が与える `E(0)=-phi'(0)` と far Robin 条件

```text
phi'(L) + phi(L) / lambda_tail = 0
```

で Poisson 方程式を解く。`phi(0)` は未知量であり、独立に指定しない。

### Ambient electrons

無限遠から interface へ向かう 1D half-Maxwellian を入射 VDF とする。interface では粒子を完全吸収し、
specular reflected population には electrostatic turning point を持つ軌道だけを含める。
単調な electron-repelling branch `q_e phi(0) >= 0` では

```text
n_e(phi; phi_0) / n_e_inf =
  exp(-q_e phi / T_e)
  [1 + erf(sqrt(q_e (phi_0 - phi) / T_e))]
  / [1 + erf(sqrt(q_e phi_0 / T_e))].
```

分母により、infinity での loss-cone を含む全密度が `n_e_inf` になるよう正規化する。これにより
passing/reflected population を二重計上しない。`q_e phi(0) < 0` の electron-attracting branch、
interface から供給される ambient-electron outgoing VDF、任意 VDF は Phase 7 では未対応とする。

### Ambient ions

無限遠で密度 `n_i_inf`、interface 向き速度 `u_i_inf > 0` の cold drifting beam とする。
粒子数 flux とエネルギー保存から

```text
u_i(phi)^2 = u_i_inf^2 - 2 q_i phi / m_i,
n_i(phi)   = n_i_inf u_i_inf / u_i(phi).
```

平方根の中が非正なら、ion-inaccessible として `no_physical_solution` を返す。入口では
`u_i_inf >= sqrt((T_e + gamma_i T_i) / m_i)` を要求する。これは presheath を解く代わりに
課すkinetic Bohm入口条件であり、既定を`gamma_i=1`とする。

### Photoelectrons

photoelectron は interface から外向きに放出される half-Maxwellian flux とする。入力は
emitted number flux `Gamma_pe` と温度 `T_pe` であり、surface から outer へ向かう電流は
`J_pe,out = q_e Gamma_pe` である。軌道エネルギーから infinity へ到達できる割合を

```text
f_escape = exp(-max(0, q_e (phi_inf - phi_0)) / T_pe),
f_return = 1 - f_escape
```

とする。空間密度は同じ flux VDF の Liouville 写像から評価し、returning branch を一度だけ
加える。tracked return と併用する場合、`photoelectron_density_model="kinetic_mean"`は空間電荷だけを供給し、return
current を表面の電荷収支へ再加算しない。粒子のreturnは`return_model`と`particle_transfer_mode`が決める。

### 電荷と電流の整合条件

Poisson source は `rho=sum(q_s n_s)` とする。解は離散 Poisson residual、far Robin residual、
および Gauss則の整合条件

```text
integral_0^L rho dz - epsilon_0 [E(L)-E(0)] = 0
```

を満たす必要がある。surface の時間発展は current-driven であり、kinetic solve 自体は浮遊電位を強制しない。
診断 current は「outer から surface へ入る正電荷」を正とする。ambient electron/ion の吸収 flux、
photoelectron escape/return、設定された external current は別々に出力する。`current_mode="floating"` は許容差内の total current を要求し、
`current_mode="driven"` は不均衡を物理的な surface-charge evolution として許す。

### Branch and failure policy

初版は `phi(z)` が interface から infinity へ単調に増加する electron-repelling branch のみを
選ぶ。Newton step は line search で positivity、ion accessibility、単調性を維持する。

- 設定またはVDFが適用条件外: `not_applicable`
- Bohm 条件違反、inaccessible ion、要求された floating root が存在しない: `no_physical_solution`
- residual/Jacobian/反復の数値的失敗: `numerical_failure`
- 収束: `converged`

別モデルへの silent fallback は行わない。非単調 virtual cathode、二重層、明示的 trapped
population、衝突 presheath、磁化 orbit は将来の別 ADR と solver を必要とする。

## 数値手法に求める条件

- 伸長格子 `z_j=L[(exp(a*j/(N-1))-1)/(exp(a)-1)]`、`a=0` は一様格子
- conservative finite-volume Poisson residual
- 密度モデルの解析微分から組み立てるtridiagonal + interface-potential border Jacobian
- bordered-tridiagonal solveによる格子点数`N`に対して`O(N)`のNewton step
- 単調分枝を維持するbacktrackingと、Newton停滞時のpseudo-transient continuation
- 前batch fieldとの差が大きい場合の適応interface-field continuation
- pseudo-transient項は反復経路だけを変え、収束判定には元の未正則化Poisson residualを使う
- UV放出が有効で初回のinterface fieldがゼロの場合は、far Robin条件を満たす微小な単調負電位profileをNewton初期値に使う。境界電位は固定せず、収束解は同じPoisson residualで決める
- 前 batch profile は同一 grid/config fingerprint の場合だけ初期値に使う
- MPI は root solve 後に status、scalar diagnostics、grid/profile を broadcast する
- restart は grid/profile、status、反復数、前 batch identity を保存する

## Consequences

強い単調 sheath と photoelectron return を線形 Debye 仮定なしで扱える。対象は決定した単調分枝に限り、
一般的な "kinetic sheath" 全体は扱わない。適用域外では明示的に停止する。legacy Zhao model は
独立の injection reference として維持し、この state ownership へ混在させない。
