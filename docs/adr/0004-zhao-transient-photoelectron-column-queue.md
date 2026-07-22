# ADR 0004: Zhao transient photoelectron column queue

- Status: experimental
- Date: 2026-07-21

## Context

ADR 0003のcharge-driven Zhao closureは、現在のinterface電場を境界条件として準定常A/B/C populationを解く。
強いUV条件では、未帯電のambient-only stateとfull-photoelectron定常branchの間に連続解がない電場域がある。
full populationを即座に適用するとこのgapで停止し、定常branchを数値的に補間するとPoisson解ではない状態を作る。

tracked `photo_raycast`粒子は放出と再吸収による表面電荷を既に更新しているが、従来はinterface通過後の滞在粒子数を
outer closureの状態として保持しなかった。過渡状態には、放出済みでまだreturnまたはescapeしていない光電子cloudの
有限inventoryとflight timeが必要である。

## Decision

`kinetic_closure="zhao_charge_driven"`で`outer_queue_enabled=true`を選んだ場合、rank-localで永続化するevent queueを使う。
queue内の全rank合計macro粒子数を水平面積$A_{xy}$で割り、photoelectron column target $N_{pe,q}$とする。
Zhao profileを$0\le z\le L$、$L=10\lambda_{D,pe}$へ128点で再標本化し、

```text
integral_0^L [n_pe,f(z; eta) + n_pe,c(z; eta)] dz = N_pe,q
```

を満たすoccupancy scale $\eta$を解く。$\eta$は確率ではなく、定常reference populationに対するscaleである。
$\eta=0$から連結するpathを$0\le\eta\le16$で探索し、1を超えるovershootを許す。queue modeは`zhao_branch="auto"`を
要求する。現行bisectionはcolumnが$\eta$とともに単調増加するpathだけを受理し、foldを含むpathは未対応として停止する。
bracket拡張中にtarget columnを初めて横切ったaccepted substepで探索を戻し、その実際の$\eta$をbracket端点に使う。
したがってtargetより先にあるfoldや別rootを探索する必要はない。
targetを無視したfull-population解、`[0,1]` clamp、disconnected branchへのjump、別closureへのsilent fallbackは行わない。

$\eta$はphotoelectron density、無限遠準中性、Sagdeev項をscaleする。旧Zhaoのzero-current式はcharge-driven rootには
含まれず、current-density診断として評価する。そのraw photoelectron emission-current項は$\eta$でscaleせず、
tracked sourceのfull currentを表す。analytic currentをsurface chargeやledgerへdepositせず、tracked粒子だけが更新する。

同じ$[0,L]$をqueueのownership領域とする。離散profile上で$L$より手前にturning pointがあればreturn event、$L$へ
到達すれば外部reservoirへ吸収されるescape eventとする。event timeはbatch内の放出時刻をmidpointで近似する。
各batchは、開始時にdue eventをpopしてclosureを予測し、粒子更新と新規event enqueue後に同じbatch indexでclosureを
再計算するpredictor/corrector順序を取る。これにより連続実行とsplit/restart実行が同じstate transitionを持つ。

flight time中にfieldを凍結する近似は、代表flight time $\tau_{outer}$だけでなく、midpoint due timeから次のbatch-start pollまでの
量子化遅延$\delta_{poll}$と、batch内crossing時刻のmidpoint近似誤差上限$\Delta t_b/2$も含め、

```text
(tau_outer + delta_poll + batch_duration/2) / field_evolution_timescale <= max_frozen_field_ratio
```

をruntimeで検査する。queue countとsigned chargeをcharge ledgerへ含める。restartはrankごとのqueue shardをatomicに書き、
schema、world size、rank、batch、global count、signed chargeを検証してから復元する。
queue-file schema 2はrank-local payloadのfingerprintを持ち、summaryの全rank fingerprintとも照合して、同じcountと
signed chargeを持つ異なるqueue shardの混在もfail closedで拒否する。

## Consequences

UV開始直後の未形成cloudから、部分populationを経て準定常branchへ連続的に進められる。return/escapeの遅延とqueue stockを
明示的に記録でき、surface chargeの二重加算を避けられる。branch連続性やcolumn targetを満たせない場合は
`no_physical_solution`で停止し、物理的に未定義な補間を隠さない。

このclosureは動的Vlasov/Poisson解ではない。$L=10\lambda_{D,pe}$、128点、batch midpoint時刻、flight中の凍結field、
平面・無衝突・非磁化という近似を持つ。production利用前に$L$、profile grid、tracked粒子数、batch duration、
`field_evolution_timescale`、interface位置の収束と、queue inventory・current・charge ledgerの保存則を確認する必要がある。
強UV runで判明したcontinuationの診断と、連成rootおよび動的outer stateへの段階移行は
[ADR 0005](0005-zhao-continuation-and-dynamic-outer.md)で定義する。
