title: periodic2静電場

Lang: [日本語](PeriodicElectrostatics.md) | [English](PeriodicElectrostatics.en.md)

# periodic2静電場

`field_bc_mode="periodic2"`の静電場は、近傍の有限画像、無限周期の横方向mode、平面平均`k=0`、
外部plasma responseという4つの成分から構成されます。x/y周期・z開放のslabでは、各成分を担当する
計算経路を分け、electrostatic snapshotへ一度ずつ加算します。

## 場を4つの成分に分ける

| 成分 | 物理的な意味 | production経路 |
| --- | --- | --- |
| primary + near images | primary cell近傍の強い局所場 | Direct/FMMの有限画像和 |
| far `k\ne0` | x/y方向に変化する無限周期遠方場 | `cached_kneq0` operator |
| surface `k=0` | 各高さより下の総電荷が作る平面平均場 | triangle-height累積多項式 |
| plasma response | 外部plasmaによるzero/nonzero応答 | 選択したouter model |

`cached_kneq0`はsurface `k=0`を除いたnonzero modeを返します。境界条件を反映した物理的な`k=0`は、
snapshotが1回だけ加えます。この分担により、横方向の無限周期補正とz方向のboundary/sheath modelを
独立に選べます。

実際の組合せは、[periodic2有限画像構成](FinitePeriodicConfiguration.html)と
[periodic2無限周期＋outer plasma構成](InfinitePeriodicOuterConfiguration.html)で説明します。

## nonzero modeの計算経路を選ぶ

| 非零modeの構成 | 用途 | 制約 |
| --- | --- | --- |
| legacy finite images | 小規模比較、有限画像model | 画像範囲外は含まない |
| `panel_spectral_reference` | triangle P0の小規模reference | Direct、softening 0、mode/quadrature収束が必要 |
| `cached_kneq0` | FMM productionの無限周期nonzero mode | x/y periodic、z nonperiodic、`exclude_k0` |
| `m2l_root_oracle` | far correctionの診断 | production hot pathには使わない |

`field_periodic_far_correction="auto"`は、現在は互換性のため`none`として動作します。無限周期のproduction計算では
`cached_kneq0`を明示します。起動時には、高水準設定とtyped `[periodic2]`の整合性を検証します。
zero-mode ownershipが矛盾している場合や、未対応のouter modelが選ばれている場合は停止します。

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

Coulomb kernelの周期画像和は収束が遅く、非中性slabの平均場はz方向のboundary closureを決めなければ
一意になりません。Ewald法は数値parameter $\alpha$を使って

$$
\frac1r=
\frac{\operatorname{erfc}(\alpha r)}r
+\frac{\operatorname{erf}(\alpha r)}r
$$

へ分けます。第1項はreal spaceで急速に減衰し、第2項は滑らかなreciprocal-space modeへ展開できます。
Ewald2Pはx/yだけを逆格子へ変換し、zを開放座標として残します。

| Ewald部分 | 評価 |
| --- | --- |
| real space | screened Coulombの有限画像和 |
| reciprocal `k\ne0` | x/y Fourier modeの有限層和 |
| `k=0` | 開放z方向の平均場として別項 |

$\alpha$は、計算量をreal spaceとreciprocal spaceへ配分するためのparameterです。Debye長や物理的なscreening率を
表すものではありません。cutoffを十分に収束させれば、結果は$\alpha$に依存しません。実装上のEwald referenceは、
設定したreal/reciprocal cutoffで評価する高精度teacherであり、無限和の解析的な厳密値ではありません。

## 遠方補正を再利用可能なoperatorにする

cold buildはEwald2P teacherから有限画像shellを引いたresidual

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N)
$$

をproxy sourceとcheck pointでsampleします。geometry、periodic length、FMM order、target topologyが固定されていれば、
root source multipoleからfar local expansionへの写像は線形になります。

$$
\mathbf L_t^\mathrm{far}=\mathbf A_t\mathbf M_\mathrm{root}
$$

$\mathbf A_t$を一度QR fitし、cacheへ保存します。電荷が変わったときはoperatorを再fitせず、現在の
$\mathbf M_\mathrm{root}$へ保存済みの行列を適用します。Ewald sumを実行するのはteacherを作るcold pathだけです。
warm particle evaluationでは、全sourceのEwald和を計算しません。

## 物理`k=0`を一度だけ加える

Ewald teacherのfull residualには、operatorを構成するための対称vacuum `k=0`が含まれます。一方、物理的な
zero modeは`lower_boundary_model`やouter modelによって決まります。FMM coreは同じsource stateから
symmetric `k=0`を再構成し、cached結果から差し引きます。

$$
K_{k\ne0}=K_\mathrm{shell}+R_\mathrm{Ewald}^{\mathrm{full}}-K_0^\mathrm{sym}
$$

その後、snapshotが選択したboundary closureの$K_0^\mathrm{physical}$を加えます。

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}
$$

`zero_mode_policy="exclude_k0"`は平均場を除外する指定ではありません。nonzero backendによる二重加算を防ぐ
ownership規則です。

## cold buildしたoperatorをbatch間で再利用する

| 時期 | 実行する処理 | 実行しない処理 |
| --- | --- | --- |
| cache miss | proxy/check配置、Ewald teacher評価、QR fit、checksum付きpublish | particle追跡 |
| cache hit | fingerprint、shape、checksum検証とoperator読込 | Ewald再評価、再fit |
| charge refresh | P2M/M2M、cached行列適用、zero-mode state更新 | operator再生成 |
| particle evaluation | near、local expansion、cached symmetric `k=0`減算、physical `k=0`加算 | all-source Ewald和 |

fingerprintにはgeometry、periodic length、FMM order、image/Ewald層、source kernel、target topology、
generator/build versionなどを含みます。fingerprintが一致しないcacheは再利用しません。cold/warmの結果は、
丸め誤差の範囲で一致する必要があります。

## 表面電荷から`k=0`を構築する

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
boundary closureを満たす物理成分です。

## z方向の境界条件で平均場を閉じる

総表面電荷$Q=\sum_iq_i$について、現行選択は次のとおりです。

| `lower_boundary_model` | $E_\mathrm{bottom}$ | $E_\mathrm{top}$ | 意味 |
| --- | ---: | ---: | --- |
| `symmetric_vacuum` | $-Q/(2\epsilon_0A)$ | $+Q/(2\epsilon_0A)$ | 上下に同じvacuum半空間がある無外場closure |
| `e_bottom_zero` | $0$ | $Q/(\epsilon_0A)$ | 下側電束を0に固定するlegacy closure |

どちらのmodelも、誘電体内部のscreeningやpolarizationは解きません。`symmetric_vacuum`は、追加のinterfaceや
誘電率を持たない最小の対称closureです。`e_bottom_zero`は過去の計算を再現するための設定であり、一般的な
物理defaultではありません。

outer modelを接続すると、surface zero modeのfield/interface条件を使ってplasma profileを作ります。split kinetic modelは
[kinetic 1D外部プラズマ](KineticOuterPlasma.html)、surfaceから統合する線形modelは
[unified linear response](UnifiedLinearResponse.html)を参照してください。

## 粒子衝突では軌道が届く周期画像を調べる

field targetはprimary periodic cellへwrapして評価しますが、軌道上の衝突・境界位置は物理座標のまま保持します。
mesh collisionでは、軌道線分が到達し得るperiodic imageを幾何的に探索します。fieldのnear-image layerと
collisionのimage boundは、それぞれ独立に決まります。詳細は
[粒子の衝突・境界イベント](ParticleEvents.html)を参照してください。

## 成分ごとに収束を確認する

- finite image modelではimage layerを増やして目的量を収束させる。
- cached modelではcache miss/hit、thread/MPI構成で同じoperator結果を確認する。
- Ewald $\alpha$、real/reciprocal layer、proxy/check設定に対するteacher/operator誤差を確認する。
- primary、near、far、symmetric `k=0` subtraction、physical `k=0`の二重加算がないことをoracleと比較する。
- Gauss residualとlower/upper boundary closureを確認する。
- 非中性cellの有限高さpotential差を、そのまま無限遠escape energyと解釈しない。

FMM内部のEwald式とoperator APIは[FMM内部実装](FMMCore.html)を参照してください。

## Code reference

- periodic FMM plan/state/evaluation: [`bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90)
- Ewald teacherとcached root operator: [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
- cached symmetric `k=0` subtraction: [`bem_coulomb_fmm_eval_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90)
- surface zero-mode plan/state: [`bem_periodic_zero_mode_plan.f90`](../src/physics/periodic_zero_mode/bem_periodic_zero_mode_plan.f90)
- zero-mode evaluation: [`bem_periodic_zero_mode_eval.f90`](../src/physics/periodic_zero_mode/bem_periodic_zero_mode_eval.f90)
- component ownershipとsnapshot合成: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
