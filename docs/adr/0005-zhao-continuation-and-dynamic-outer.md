# ADR 0005: Zhao continuationの診断と動的outerへの段階移行

- Status: accepted for staged implementation
- Date: 2026-07-22

## Context

ADR 0004のqueue closureは、interface外に滞在する光電子の総column
$N_{pe,q}$を永続stateとして持ち、各batchで同じcolumnを持つZhao準定常profileを解く。
これは未帯電状態へfull photoelectron populationを即時適用する問題を避けるが、profileの
free、captured、outgoing、returning population、normal energy、飛行phaseを1個のoccupancy
$\eta$へ縮約している。

強UVのsingle-sphere plus plane診断runは、15 batchまでType Bを追跡した後、batch 16の
interface field $E_I=0.9072962759\,\mathrm{V/m}$で停止した。最後に受理した状態は
$N_{pe,q}=9.32021\times10^7\,\mathrm{m^{-2}}$、$\eta=0.232040$で、column相対残差は
約$3.6\times10^{-9}$だった。停止時のmessageは次である。

```text
photoelectron-column Zhao path ended while increasing eta:
bounded Zhao eta continuation exhausted step-halving:
automatic Zhao continuation step moved too far on one root branch
```

この時点までreturn eventはまだ1件もpopされていない。したがって今回の停止は、returned
photoelectron feedbackやcolumn residualの収束不良ではない。固定$\eta$のNewton solveは同じ
Type Bラベルのcandidateを返したが、直前rootに対する正規化jumpがguard 0.25を超え、刻みを
下限まで半減しても近傍rootを回復できなかった。

現行algorithmは、前状態から$\eta$を固定して$E_I$を変更し、その後$E_I$を固定して
$\eta$を変更する。さらに$N_{pe}(\eta)$が単調なpathだけを受理する。A/B/Cはprofile形状の分類で
あってroot identityではないため、同じB内に複数rootがある場合や、$\eta$を座標としたroot
manifoldにfoldがある場合を区別できない。このため現在のstatusだけでは、次を判別できない。

- 近傍rootは存在するがNewton basinが遠いrootへ切り替わった
- $\eta$を座標にしたfoldであり、columnを座標にすれば解を追跡できる
- target columnとの交点を持たない物理的なbranch終端である
- Zhao準定常manifold自体が過渡stateを表せない

一方、UVなしの対応runは5000 batchを完走し、Type C、
$\phi_I=-7.25080\,\mathrm{V}$、$E_I=-0.688793\,\mathrm{V/m}$、total-current mismatch 0.93%へ
到達した。従ってno-photo Zhao経路を置き換える理由はなく、強UVのbranch追跡と過渡closureを
局所的に改修すべきである。

## Decision

改修を次の三段階に分ける。単純なiteration増加、jump guardの緩和、または別closureへの
silent fallbackは行わない。

### 1. Branch atlasと失敗分類

Zhao branchごとの残差$G_b$とJacobian評価をsolverから分離する。固定$E_I$では、encoded root
$y$と$\eta$をまとめた$z=(y,\eta)$に対し、pseudo-arclength条件

```text
G_b(z; E_I) = 0
t_prev dot (z - z_predictor) = 0
```

を使って連結root curveを追跡する。Type B/Cの$y$は$(\phi_0,\log n_{e,\infty})$、Type Aは
$(\phi_0,\phi_m,\log n_{e,\infty})$とする。branch atlasは最初にtestまたは診断utilityとして
実装し、runtime fallbackにはしない。

各accepted pointで$E_I$、$\eta$、column、current、root residual、jump metric、Jacobianの
conditioning指標、接線と$dN_{pe}/ds$を記録する。停止statusを少なくとも次に分ける。

- `numerical_failure`: residualまたはlinear solveが数値的に収束しない
- `fold_detected`: 選んだ座標のfoldだが、連結root curveは継続する
- `physical_endpoint`: density正値性、Sagdeev実数条件、branch topologyなどの物理制約で終端する
- `target_unreachable`: 追跡した連結curveにtarget columnとの交点がない

Newton不収束だけを`no_physical_solution`とは判定しない。pseudo-arclengthがfold後の数学的rootを
見つけても、そのrootの動的安定性が不明なまま自動採用しない。

### 2. 到達可能な準定常rootの連成solve

branch atlasが前batchのrootとtargetを結ぶ連結curveを確認した場合、runtime column closureを
「固定$\eta$ root＋外側bisection」から次の連立問題へ変更する。

```text
G_b(y; E_I, eta) = 0
N_pe(y, eta; E_I) - N_pe,q = 0
```

前batchから$E_I$と$N_{pe,q}$を同時に

```text
E(lambda) = E_old + lambda * (E_new - E_old)
N(lambda) = N_old + lambda * (N_new - N_old)
```

とhomotopyし、previous tangentをpredictorに使うsafeguarded bordered Newtonで解く。root identityは
A/B/Cの文字だけでなく、連結curve、予測点、接線方向で決める。A/B/C間の接続は既存の退化条件を
満たす点だけで許可する。continuation stateはrestartへ保存し、straight runとsplit runを一致させる。

### 3. 準定常manifold外の動的outer state

branch atlasにtarget交点がない場合、またはouter flight timeに対してfieldが速く変化してfrozen-field
guardを満たさない場合は、stationary Zhao profileを補間しない。この場合は現行scalar queueを、少なくとも
normal-energy bin、terminal class、flight phaseまたはage cohortを持つdynamic inventoryへ拡張する。
outgoing/free、captured/outgoing、captured/returningを別stateとし、同じ総columnでも異なるVDF履歴を区別する。

このreduced dynamic modelでは、Zhao解を毎batch強制するのではなく、定常fixed pointまたはfar-terminal
conditionとして使う。particle number、signed charge、species別face current、Gauss residual、outer kinetic
energyとfield energyを保存・出力する。enqueue時に固定したterminal stateだけでなく、interface entry state、
normal energy、profile generationもrestart fingerprintへ含める。

full 1D kinetic Vlasov/particle-Poisson modelは、まずflat-planeの基準解として構築する。reduced modelが
energy-bin、outer timestep、domain lengthへ収束しない、またはoscillatory/trapped topologyを再現できない場合に
限ってonline closureへ昇格する。

## Rejected alternatives

- `branch_same_root_step_limit`や最小刻みを緩和する。遠い同ラベルrootへのjumpを受理し得る。
- continuation試行回数だけを増やす。座標foldでは刻みを小さくしても解決しない。
- full $\eta=1$、branch 0、linear Debye、instant return、`absorbing_maxwellian`へfallbackする。保存則と
  provenanceを保ったまま物理modelが変わったことを識別できない。
- disconnected branch間を補間する。Poissonを満たさないprofileを作る。
- pseudo-arclengthで得たfold後のrootを無条件に採用する。存在と動的安定性は同じではない。
- 最初からfull 1D kineticをproduction pathへ統合する。検証範囲と計算費用を同時に拡大し過ぎる。

## Validation

実装は次の順で検証する。

1. 既知のA/B/C/no-photo rootと退化点を、direct scanとbranch atlasで一致させる。今回の直前rootと
   failure targetをfixture化し、`numerical_failure`、fold、endpointを再現可能に分類する。
2. 同じBラベルに複数rootがあるfixtureで、previous tangentに連結したrootだけを選ぶ。column自体にfoldがあり
   target交点がない場合は`target_unreachable`でfail closedする。
3. 同一seedの強UV runで少なくともbatch 32まで進み、$E_I=0.9073\,\mathrm{V/m}$および
   $N_{pe,q}=9.32021\times10^7\,\mathrm{m^{-2}}$を越える。受理rootのcolumn相対残差を$10^{-8}$以下とし、
   rootまたはbranchの不連続jumpを許さない。
4. 実際のreturn eventがpopされるまで、少なくとも100--200 batchを実行し、inventory、signed charge、
   current、Gauss residual、energy budgetを確認する。
5. UVなしrunを回帰させ、$\phi_I$と$E_I$を上記baselineの1%以内、total-current mismatchを1%以下、
   charge-ledger相対残差を$10^{-10}$以下に保つ。定常判定は単一step差ではなく、$Q$、$\phi$、$E$、$J$の
   window統計で行う。max-step unresolved particle率はbaselineの0.884%を悪化させず、publication用には
   別途0.1%未満を目標に収束確認する。
6. flat-planeでUV turn-on/off、batch duration、outer domain length、profile grid、energy bin、particle数、
   restart分割を収束確認した後にsingle-sphere plus planeへ戻す。

既存のZhao unit test 15件が計算ノードで通る状態を改修前baselineとする。現行queueとno-photo Type Cは、
上記の新statusまたはstateを使わない既存caseで挙動を変えない。

## Consequences

第一段階だけでは強UV runを完走させないが、数値的なroot追跡失敗を物理解なしと誤認しなくなる。
第二段階は、解が存在する場合に現在のscalar queueを保ったまま人工的なfield-then-$\eta$経路を除去する。
第三段階は、stationary Zhao manifoldに存在しない過渡状態を明示的なdynamic stateとして扱う。

この順序により、最小のsolver改修で到達可能な準定常解を回復しつつ、強UVで本当に必要な動的physicsを
silent approximationで隠さない。各段階のartifactとstatusを残すため、後続実装がどの仮定を追加したかを
run provenanceから判別できる。
