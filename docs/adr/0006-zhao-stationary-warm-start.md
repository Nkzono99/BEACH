# ADR 0006: Zhao 定常 warm start と instant outer boundary

- Status: experimental
- Date: 2026-07-22

## Context

[ADR 0003](0003-zhao-charge-driven-outer-closure.md)は、現在の表面電荷が決めるinterface電場$E_I$を
満たす Zhao A/B/C profileを準定常outer stateとし、zero-currentを帯電過程の診断に残した。強UVでは、
未帯電のambient-only stateとfull-photoelectron定常branchの間に Zhao 解のない$E_I$領域がある。

[ADR 0004](0004-zhao-transient-photoelectron-column-queue.md)は、このgapを物理過渡として扱うため、外部光電子columnと
flight delayをqueue stateにした。しかし[ADR 0005](0005-zhao-continuation-and-dynamic-outer.md)の強UV診断では、現行の
scalar column closureがbatch 16のtargetへ連続する定常 Zhao componentを確認できず、過渡応答にはより動的な
outer stateが必要と判断した。

一方、直近の研究目的には、小さなBEACH領域の定常帯電、平均電流、離脱力、仕事の評価がある。
この問いでは、定常 Zhao sheathをouter boundary operatorとして使い、無限遠reservoir補正、粒子の
escape/return、無限遠基準電位を同じprofileで接続できればよい。未帯電からのUV立上がりは別の物理問題である。

## Decision

`coupling.steady_start_mode="zhao_floating"`をopt-inの定常初期化として追加する。既定は`none`とし、
queue過渡closureを置換しない。

新規実行の最初batch前に次を行う。

1. 設定済みの無限遠ambient density、temperature、normal drift、UV sourceから Zhao 零電流定常根を解く。
2. その根から`phi(infinity)=0`のkinetic profileとinterface電場$E_I$を構成する。
3. 水平periodic cellの面積を$A$とし、zero-mode境界条件に対応する平面電荷を与える。

```text
symmetric_vacuum: Q_seed = 2 epsilon_0 A E_I
e_bottom_zero:    Q_seed =   epsilon_0 A E_I
```

4. `coupling.steady_start_mesh_id`で選んだ平面の三角形に$Q_{seed}$を面積比で配る。
5. 同じprofileを初回outer stateとし、無限遠reservoirからinterfaceへの流入VDF、
   `phi_interface-phi_infinity`、およびinstant escape/returnの全てに使う。

選択meshは、水平、同一高さ、outer interfaceより下で1つの平面をなし、面積が$A$と一致しなければならない。
実験的な初期実装は`mesh.mode="template"`だけを受け入れ、生成済みplane templateの非重複・無欠損tilingを
構築時の不変条件として使う。任意OBJ投影のunionを面積と外形だけから推測しない。
新規実行は全meshの初期電荷0を要求する。plane + sphereでplaneを選んだ場合、planeだけをseedし、sphereは
中性で開始する。

零電流条件は初期状態を構成するときにだけ課す。後続batchはADR 0003のcharge-driven contractに戻り、
現在の表面電荷が$E_I$を決め、zero-currentは診断である。analytic currentはsurface chargeやledgerへdepositせず、
tracked放出、流入、再吸収と二重計上しない。

このmodeは次の構成に限る。

- `outer_plasma.model="kinetic_1d"`
- `outer_plasma.kinetic_closure="zhao_charge_driven"`
- `outer_plasma.return_model="kinetic_1d_profile_return"`
- `coupling.particle_transfer_mode="electrostatic_1d_instant_return"`
- `coupling.outer_queue_enabled=false`
- `periodic2.zero_mode_policy="exclude_k0"`と対応するlower boundary model

同一configで`output.resume=true`とした場合は初期化を再実行せず、checkpointからmesh電荷と完全なouter stateを
復元する。resumeでouter stateが欠落している場合は再seedや新規実行へfallbackせず停止する。標準出力には、
新規実行でresolved branch、$E_I$、$Q_{seed}$、mesh ID、resumeで復元したときは処理済みbatch数を記録する。

## Rejected alternatives

- 1つのcoarse batchの粒子weightでType A可解域へ飛ぶ。batch durationとMonte Carlo sourceに依存し、初期条件の意味が不明確になる。
- branch `0`とfull-population branchの間を数値補間する。Zhao-Poisson解でないprofileを定常解として扱うことになる。
- queueが失敗したときにwarm startへfallbackする。過渡modelから定常初期条件へ暗黙に物理を変える。
- planeとsphereの両方へ一様にseedする。1D zero modeが決めるのはcell平均平面電荷であり、局所物体の電荷分配ではない。

## Consequences

強UVでも、未帯電からのbranch gapを積分せず、解析的に定義された Zhao 定常branchから小さな通常batchを
開始できる。BEACHの局所帯電とouter reservoir補正、電位基準、return/escapeが一つのprofileで整合する。

一方、このRUNは未帯電状態からの過渡応答、UV turn-on時間、return-current遅延、または定常branchの動的選択を
予測しない。queue過渡closureと将来のdynamic outer modelの動機は残る。warm startで得た解の存在は、物理的な一意性や
安定性を保証しない。

## Validation

1. flat planeで、seed後の$E_I$が零電流根の$E_I$と一致し、total currentが解析tolerance内で0になること。
2. `symmetric_vacuum`と`e_bottom_zero`の両方で$Q_{seed}$式、面積配分、Gauss residualを確認すること。
3. plane + sphereでplaneだけがseedされ、sphereが電荷0で開始すること。
4. reservoirの流入数・cutoff・速度写像とoutward粒子のreturn/escapeが同じprofile barrierを使うこと。
5. straight runとsplit-resumeでmesh電荷、outer state、観測量が一致し、resume時に再seedされないこと。
6. queue過渡case、no-photo経路、`steady_start_mode="none"`の既存caseが回帰すること。
7. publication用には、独立に緩和させた状態または摂動したseedから同じ定常電位、電流、帯電、力へ戻るかを確認すること。
