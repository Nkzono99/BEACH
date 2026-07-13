# ADR 0002: unified periodic outer domain と handoff-only interface

- Status: accepted for Phase 8 field domain and Phase 9 explicit orbit
- Date: 2026-07-11

## Context

従来の split model は、surface と outer interface の間を vacuum、interface 上の `k!=0` を
無視できると仮定した。rough surface では、local plasma charge を無視するには interface を
低く、横方向 mode を減衰させるには高く置く必要があり、両条件を満たす window が存在しない。

## Decision

field domain は surface projection の下端から far Robin boundary まで連続させる。particle
interface は field の境界ではなく ownership/handoff 面だけとし、位置を変えても同じ field
profile と orbit energy を参照する。

### Accessible area

Phase 8 の 1D mean closure は single-valued topography に限定する。周期セル内の各 `(x,y)` に
対し、plasma-facing surface の最上交点を `h(x,y)` とし、

```text
A_access(z) / A_xy = mean_xy[ I(z > h(x,y)) ]
```

と定義する。数値実装は cell-centered uniform `(x,y)` sample と vertical ray intersection を
使い、sample 数を倍増した収束を要求する。複数の上向き交点、closed cavity、overhang、横方向
迂回だけで reservoir へ接続する領域はこの定義では表せないため `not_applicable` とする。

plasma density closure `n_closure(z)` は accessible volume 内の条件付き密度であり、full-cell
zero-mode source は

```text
rho_mean(z) = f_access(z) rho_closure(z)
```

とする。tracked-particle residence histogram は
`n_hist = sum(weight)/(A_access dz observation_time)` という独立診断であり、closure source と
加算しない。configured tolerance を超える不一致は production validity を失う。

### Unified zero mode

surface panel source は既存の exact height projection `F_i(z)` を用いる。plasma source は同じ
1D grid の finite volume へ `f_access rho_closure` として入れる。bottom field と far Robin
condition を一つの Poisson solve に課し、interface Neumann condition は廃止する。Gauss residual
は surface、local mean plasma、far tail を含む domain 全体で評価する。

### Nonzero-mode plasma tail

表面最高点から幾何的に決める応答開始面 `z_r` より下を vacuum、上を線形応答 plasma とする
単一 mode を考える。`z_r` は particle ownership/handoff 面から独立であり、全 source panel より
厳密に上へ置く。
`k=sqrt(kx^2+ky^2)`、`alpha=sqrt(k^2+kappa^2)` とする。surface source が vacuum 解として
応答開始面へ運ぶ incident amplitude を `I_k` とすると、下側反射と上側透過は

```text
R_k = (k - alpha) / (k + alpha) I_k
T_k = 2 k / (k + alpha) I_k.
```

`z<=z_r` では free-space periodic field に `R_k exp(k(z-z_r))` を加え、`z>=z_r` では
free-space continuation を `T_k exp(-alpha(z-z_r))` へ置換する。この構成は potential、normal
field、tangential field を応答開始面で連続にする。`kappa=0` では補正は厳密に 0 になる。

finite mode truncation は configured `mode_layers` で制御する。neglected-amplitude bound は
将来追加する診断であり、現行実装は出力しない。このため現状は`mode_layers`とpanel quadratureを増やした
目的量の収束確認をproduction受理条件とする。各retained modeについて`max |q_s phi_k|/T_s`が
linearity toleranceを超える場合、またはmode間 nonlinear couplingが必要な場合は
`not_applicable`とし、1D closureで代用しない。

### Outer orbit

`electrostatic_3d_explicit_orbit` は、zero mode と継続した nonzero tail の全 3D electric field を
固定刻み velocity-Verlet で追跡する。ownership 面へ戻れば periodic x/y を wrap して local stepper
へ返し、unified grid 上端の far plane を外向きに横切れば infinity escape とする。全エネルギーの
相対誤差、outer flight time、frozen-field ratio を実測し、設定上限を超えた場合は停止する。
`max_steps` 到達粒子は discard せず、persistent queue が必要な未解決状態として停止する。
外部領域の `B0` を無視できる定量条件は未決のため、この mode は `B0=0` だけを受理する。

Phase-mixed statistical transfer を使う場合は、handoff 面で横位相を一様化し、mode ごとの
ponderomotive correction を無視できる linearity bound を要求する。両方式を同じ粒子へ重複適用
しない。statistical return は着地点分布と遅延が未仕様なので、擬似乱数だけで代用せず
fail-closed とする。

## Failure policy

- height-field geometry でない、accessible-grid 未収束: `not_applicable`
- nonzero linearity または truncation tolerance 超過: `not_applicable`
- ion accessibility/Bohm/kinetic branch 不成立: `no_physical_solution`
- Poisson、mode reconstruction、orbit integration の数値失敗: `numerical_failure`
- 3D nonlinear particle space charge が必要: 明示停止し、将来の 3D PIC/Poisson model へ送る

interface を動かして解が変化する場合は tolerance を緩めて受理せず、実装または model contract の
不整合として扱う。
