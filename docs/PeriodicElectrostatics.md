title: periodic2静電場

Lang: [日本語](PeriodicElectrostatics.md) | [English](PeriodicElectrostatics.en.md)

# periodic2静電場

`field_bc_mode="periodic2"`の静電場は、近傍の有限画像、無限周期の横方向mode、平面平均`k=0`、
外部plasma responseという4つの成分から構成されます。x/y周期・z開放のslabでは、各成分を担当する
計算経路を分け、最終的な電場・電位へ一度ずつ加算します。

## 場を4つの成分に分ける

| 成分 | 物理的な意味 | production経路 |
| --- | --- | --- |
| primary + near images | primary cell近傍の強い局所場 | Direct/FMMの有限画像和 |
| far `k\ne0` | x/y方向に変化する無限周期遠方場 | `cached_kneq0` operator |
| surface `k=0` | 各高さより下の総電荷が作る平面平均場 | triangle-height累積多項式 |
| plasma response | 外部plasmaによるzero/nonzero応答 | 選択したouter model |

`cached_kneq0`はsurface `k=0`を除いたnonzero modeを返します。境界条件を反映した物理的な`k=0`は、
場の合成処理が1回だけ加えます。この分担により、横方向の無限周期補正とz方向のboundary/sheath modelを
独立に選べます。

実際の組合せは、[periodic2有限画像構成](FinitePeriodicConfiguration.html)と
[periodic2無限周期＋outer plasma構成](InfinitePeriodicOuterConfiguration.html)で説明します。

## nonzero modeの計算経路を選ぶ

| 非零modeの構成 | 用途 | 制約 |
| --- | --- | --- |
| finite images | 小規模比較、有限画像model | 画像範囲外は含まない |
| `panel_spectral_reference` | triangle P0の小規模reference | Direct、mode/quadrature収束が必要 |
| `cached_kneq0` | FMM productionの無限周期nonzero mode | x/y periodic、z nonperiodic、`exclude_k0` |

`field_periodic_far_correction="auto"`は、現在は互換性のため`none`として動作します。無限周期のproduction計算では
`cached_kneq0`を明示します。起動時には、高水準設定とtyped `[periodic2]`の整合性を検証します。
zero-mode ownershipが矛盾している場合や、未対応のouter modelが選ばれている場合は停止します。
`m2l_root_oracle`は削除済みで、設定すると起動時にrejectされます。

Ewald2P teacher、root multipoleからlocal展開へのoperator、cacheとFMM stateの接続は
[periodic2遠方補正](PeriodicFarCorrection.html)に分けています。このページでは、それを場全体のnonzero成分として扱います。

## 有限画像和が含む範囲を定める

primary cellと設定した画像層$N$について

$$
(i,j)\in[-N,N]^2
$$

のsourceを陽に加えます。near fieldは元のkernelで評価できますが、範囲外に続く画像が作る滑らかなfar fieldは
含まれません。したがって、`field_periodic_far_correction="none"`は有限画像modelです。$N$を増やした結果が
収束するまでは、無限周期解として扱えません。

FMMのnear image層はfar operatorを作るときに差し引くshellと一致する必要があります。cache fingerprintが画像層を
identityに含むのはこのためです。

## Ewald2Pで無限周期の遠方場を分離する

`cached_kneq0`はEwald2P teacherと有限画像shellの差を、root multipoleからtarget localへのoperatorとして
適用します。cached結果からteacher由来の対称`k=0`を除き、場の合成処理が選択した物理`k=0`を一度だけ加えます。

$$
K_\mathrm{surface}
=\left(K_\mathrm{shell}+R_\mathrm{Ewald}^{\mathrm{full}}-K_0^\mathrm{sym}\right)
+K_0^\mathrm{physical}
$$

この式の括弧内がnonzero backendの責務です。`zero_mode_policy="exclude_k0"`は平均場を捨てる指定ではなく、
二重加算を防ぐownership規則です。Ewald分割、operator fit、FMMへの注入位置、cache lifecycleは
[periodic2遠方補正](PeriodicFarCorrection.html)に分離しています。

## 物理`k=0`を一度だけ加える

triangle $i$の総電荷を$q_i$、面積のうち高さ$z$以下にある割合を$F_i(z)$とすると、平面平均された累積電荷は

$$
C(z)=\sum_iq_iF_i(z)
$$

です。$F_i$は3頂点の高さの間で区分二次関数になります。geometry planは全頂点の高さをsortしてbreakpointを作り、
各triangleがそれぞれの区間へ加える二次係数を保存します。水平triangleはsheet chargeとして別に保持します。
面上評価ではminus trace、plus trace、principal valueを区別します。

batchで$q_i$が変わると、保存済みgeometry係数へ$q_i$を掛けて区間差分を作り、prefix sumで

$$
C(z)=a_0+a_1z+a_2z^2
$$

の区間係数とそのprimitiveを更新します。geometry planは再構築しません。

## Gauss則からzero-mode fieldとpotentialを求める

cell面積$A=L_xL_y$、下側far fieldを$E_\mathrm{bottom}$とするとGauss則から

$$
E_0(z)=E_\mathrm{bottom}+\frac{C(z)}{\epsilon_0A}
$$

です。gauge点$(z_g,\phi_g)$からのpotentialは

$$
\phi_0(z)=\phi_g-E_\mathrm{bottom}(z-z_g)
-\frac1{\epsilon_0A}\int_{z_g}^zC(\zeta)\,d\zeta
$$

です。breakpoint区間を二分探索し、区間内の二次式と三次primitiveを使うため、1点評価は$O(\log N_z)$です。

非中性cellではz遠方に一定fieldと線形potentialが残り得ます。zero modeは数値的に消してよい成分ではなく、Gauss則と
境界条件を満たす物理成分です。

## z方向の境界条件で平均場を閉じる

総表面電荷$Q=\sum_iq_i$について、現行選択は次のとおりです。

| `lower_boundary_model` | $E_\mathrm{bottom}$ | $E_\mathrm{top}$ | 意味 |
| --- | ---: | ---: | --- |
| `symmetric_vacuum` | $-Q/(2\epsilon_0A)$ | $+Q/(2\epsilon_0A)$ | 上下に同じvacuum半空間がある無外場境界条件 |
| `e_bottom_zero` | $0$ | $Q/(\epsilon_0A)$ | 下側電束を0に固定するlegacy境界条件 |

どちらのmodelも、誘電体内部のscreeningやpolarizationは解きません。`symmetric_vacuum`は、追加のinterfaceや
誘電率を持たない最小の対称境界条件です。`e_bottom_zero`は過去の計算を再現するための設定であり、一般的な
物理defaultではありません。

outer modelを接続すると、surface zero modeのfield/interface条件を使ってplasma profileを作ります。標準のsplit kinetic
構成は[kinetic 1D外部プラズマ](KineticOuterPlasma.html)で説明します。

## 粒子衝突では軌道が届く周期画像を調べる

field targetはprimary periodic cellへwrapして評価しますが、軌道上の衝突・境界位置は物理座標のまま保持します。
mesh collisionでは、軌道線分が到達し得るperiodic imageを幾何的に探索します。fieldのnear-image layerと
collisionのimage boundは、それぞれ独立に決まります。[<sup>1</sup>](ParticleEvents.html)

## 成分ごとに収束を確認する

- finite image modelではimage layerを増やして目的量を収束させる。
- cached modelではcache miss/hit、thread/MPI構成で同じoperator結果を確認する。
- Ewald $\alpha$、real/reciprocal layer、proxy/check設定に対するteacher/operator誤差を確認する。
- primary、near、far、symmetric `k=0` subtraction、physical `k=0`の二重加算がないことをoracleと比較する。
- Gauss residualと上下の境界条件を確認する。
- 非中性cellの有限高さpotential差を、そのまま無限遠escape energyと解釈しない。

FMM内部のEwald式とoperator APIは[FMM内部実装](FMMCore.html)にまとめています。

## Code reference

- periodic FMM plan/state/evaluation: [`bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90)
- Ewald teacherとcached root operator: [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
- cached symmetric `k=0` subtraction: [`bem_coulomb_fmm_eval_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90)
- surface zero-mode plan/state: [`bem_periodic_zero_mode_plan.f90`](../src/physics/periodic_zero_mode/bem_periodic_zero_mode_plan.f90)
- zero-mode evaluation: [`bem_periodic_zero_mode_eval.f90`](../src/physics/periodic_zero_mode/bem_periodic_zero_mode_eval.f90)
- component ownershipと場の合成: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
