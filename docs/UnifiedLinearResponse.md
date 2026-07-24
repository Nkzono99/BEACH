title: 高度な粗面線形screening

Lang: [日本語](UnifiedLinearResponse.md) | [English](UnifiedLinearResponse.en.md)

# 高度な粗面線形screening（`unified_linear_response`）

`external_boundary.field.model="unified_linear_response"`は、rough surfaceの高さ範囲、表面電荷の平面平均source、
plasmaが占有できる面積、線形Debye応答を一つのfield modelにまとめます。surfaceとownership interfaceの間に
横方向modeが十分に減衰するvacuum windowを仮定できない場合の線形modelです。

これは標準の外部シースモデルではなく、roughnessとplasma responseが同じ領域に重なる場合の高度な線形screening
です。species VDF、Bohm条件、電流balanceを含む自己整合な外部シースには`kinetic_1d`を使います。
unified modelはfieldだけを所有し、source VDFは所有しません。流入補正は
`external_boundary.particles.inflow_model`で独立に選びます。
`unified_linear_response`を`kinetic_1d`の高精度版や自動fallbackとして選んではいけません。

## split windowを一つのfield solveへ置き換える

`kinetic_1d`は、ownership interfaceより下のsurface fieldと上の1D outer profileを接続します。
unified modelはsplit interfaceをfield境界にせず、surface projectionからfar boundaryまで一つのzero-mode Poisson gridを解き、nonzero modeも
plasma tailへ連続接続します。

`sim.box_max[2]`から導出されるinterfaceはz-high box面上の粒子ownership境界であり、field solveの境界でもplasma responseの
開始面でもありません。mesh geometryとDebye長が同じなら、interface高さだけを動かしてlocal field profileが
変わらないことが、このmodelに求められる条件です。

plasma responseにはDebye–Hückel型の線形関係

$$
\rho_\mathrm{DH}(z)=-\epsilon_0\kappa^2\phi(z),
\qquad
\kappa=\lambda_D^{-1}
$$

を用います。ambient species別VDF、Bohm条件、浮遊電流balance、photoelectron mean densityを状態に含む構成は、
[kinetic 1D外部プラズマ](KineticOuterPlasma.html)です。

## rough surfaceからplasma占有率を求める

periodic cell内の各$(x,y)$で、plasma側から最初に見える最上面高さを$h(x,y)$とします。高さ$z$でplasmaが占有できる
水平面積率は

$$
f_\mathrm{access}(z)=\frac{1}{A_{xy}}
\int_{A_{xy}}I[z>h(x,y)]\,dx\,dy
$$

です。full-cell平均のplasma chargeは

$$
\rho_\mathrm{plasma}(z)
=f_\mathrm{access}(z)\rho_\mathrm{DH}(z)
=-\epsilon_0f_\mathrm{access}(z)\kappa^2\phi(z)
$$

となります。surfaceより下で$f_\mathrm{access}=0$、全surfaceより上で1、roughness範囲で0と1の間を取り、
solidが占める水平面積へplasma chargeを置くことを避けます。

$h(x,y)$は`interface_sample_n`四方のcell-centered vertical rayでsampleし、各軸2倍のresolutionでも再計算します。

$$
\max_z|\Delta f_\mathrm{access}(z)|
\le\texttt{accessible\_fraction\_tolerance}
$$

を満たさなければ停止します。overhang、closed cavity、同じvertical rayに複数のplasma-facing交点があるgeometry、
横方向へ迂回してだけreservoirへ接続する細孔は、single-valued height fieldで表せないため適用外です。

## 表面電荷とplasma応答からzero modeを解く

zero-mode gridは

$$
z_{\min}=z_{\mathrm{mesh,min}}-\lambda_D,
\qquad
z_{\max}=z_{\mathrm{mesh,max}}+10\lambda_D
$$

で定まる区間を`unified_grid_points`点で分割します。既定は129点、最小は17点です。triangleの高さ方向への厳密projectionから
surface zero-mode field $E_s(z)$を同じgridで評価し、離散微分

$$
\rho_s(z)=\epsilon_0\frac{dE_s}{dz}
$$

をsurface sourceとして

$$
\frac{d^2\phi}{dz^2}
-f_\mathrm{access}(z)\kappa^2\phi(z)
=-\frac{\rho_s(z)}{\epsilon_0}
$$

を解きます。下端には`lower_boundary_model`が選ぶ$E_\mathrm{bottom}$を与え、上端には

$$
\phi'(z_{\max})+\frac{\phi(z_{\max})}{\lambda_D}=0
$$

というfar Robin条件です。nonuniform spacingに対応したtridiagonal離散式を$O(N_z)$で解きます。上端より先は
同じDebye長で指数外挿するため、ownership面がgrid上端より高くてもzero-mode fieldを連続に評価できます。

flat surface、$f_\mathrm{access}=1$、surface sourceなしなら$\phi\propto e^{-z/\lambda_D}$へ戻ります。
回帰検証はこの解析limit、grid refinement、surfaceとplasma chargeを合わせたGauss則の整合性を使います。

## 真空nonzero modeをscreened tailへ接続する

periodic solverが作る真空の$k\ne0$ modeを$e^{-kz}$のまま無限遠まで延長せず、surface最高点の直上

$$
z_r=z_{\mathrm{mesh,max}}+
\max(\epsilon_\mathrm{roundoff},10^{-6}\lambda_D)
$$

で線形plasma responseへ接続します。response start $z_r$もownership interfaceとは独立です。

各Fourier modeについて

$$
k=\sqrt{k_x^2+k_y^2},
\qquad
\alpha=\sqrt{k^2+\kappa^2}
$$

とし、真空側incident amplitude $I_k$に対して

$$
R_k=\frac{k-\alpha}{k+\alpha},
\qquad
T_k=\frac{2k}{k+\alpha}
$$

を使います。$z\le z_r$ではbase真空場へ$R_kI_ke^{k(z-z_r)}$を加えます。$z>z_r$ではbaseに含まれる
$I_ke^{-k(z-z_r)}$を差し引き、$T_kI_ke^{-\alpha(z-z_r)}$へ置き換えます。これにより$z_r$で電位、法線電場、
接線電場が連続します。$\kappa\to0$では$R_k\to0$、$T_k\to1$となり補正は消えます。

mode amplitudeは`triangle_p0` panelをDuffy quadratureで面積積分し、`periodic2.reference_mode_layers`まで保持します。
省略modeの自動誤差上限は現行実装にないため、`reference_mode_layers`と`panel_quadrature_order`を増やし、電位、電場、
力など目的量の収束を確認します。

## 線形近似の範囲内かを判定する

zero modeは

$$
\eta_0=\frac{\max_z|\phi_0(z)|}{V_T}
$$

を、retained nonzero modeは透過振幅から$\eta_k=|q\phi_{k,\mathrm{transmitted}}|/T$を評価します。
$\max(\eta_0,\max_k\eta_k)$が`max_linearity_ratio`を越えた場合は、非線形な横方向結合が必要になるため停止します。
これはsmall-perturbationの運用gateであり、誤差そのものの上限ではありません。

## geometry・source・粒子transferの適用範囲

| 項目 | 適用条件 |
| --- | --- |
| geometry | x/y periodic、z open、single-valued plasma-facing height field |
| source kernel | `triangle_p0` |
| prescribed field | `sim.e0=0` |
| mean plasma | scalar linear Debye responseのみ |
| species model | species別VDF、Bohm条件、光電子の平均密度モデルなし |
| particle transfer | `external_boundary.particles.mode="local_source"`、または3D外部軌道を使う`"same_batch"` |
| magnetic field | `particles.mode="same_batch"`は`sim.b0=0` |
| failure | 適用条件を満たさなければ別モデルへfallbackせず停止 |

`max_gap_ratio`と`max_local_charge_ratio`はsplit scalar-interface model向けdiagnosticです。unified固有の受理判定は
accessible-area収束とzero/nonzero linearityを中心に行います。

## zero/nonzero responseを場へ一度ずつ加える

1回のfield評価は次を一度ずつ合成します。

1. production FMMまたはspectral referenceのperiodic $k\ne0$ base field。
2. unified Poisson profileの$k=0$ field。
3. retained $k\ne0$ modeのreflection/transmission correction。

surface chargeをcommitして場を更新するたび、surface zero mode、unified linear solve、nonzero-tail amplitudeを更新します。
現行のunified経路は、`outer_update_stride`によるsolve skipを行いません。MPI rootがtridiagonal solveを実行し、
statusと$z,\phi,E,\rho$を全rankへbroadcastします。

`external_boundary.particles.mode="local_source"`ならz-high粒子は通常のopen境界に従います。`"same_batch"`なら
同じzero/nonzero field中で外部軌道を追跡します。[<sup>1</sup>](ParticleEscapeReturn.html)

## geometry・grid・modeを独立に収束させる

`outer_plasma_profile.csv`にはzero-modeの$z,\phi,E,\rho_\mathrm{plasma}$を保存します。`summary.txt`では
`outer_linearity_ratio`、`outer_nonzero_tail_linearity`、`outer_accessible_fraction_min/max`、
`outer_accessible_fraction_refinement_error`、`outer_response_start_z_m`、Gauss residualを確認します。

productionでは少なくとも次を変えます。

- `unified_grid_points`によるzero-mode grid refinement。
- `interface_sample_n`と`accessible_fraction_tolerance`によるgeometry sampling。
- `reference_mode_layers`と`panel_quadrature_order`によるnonzero tail。
- explicit orbit使用時の`outer_orbit_dt`。

## Code reference

- accessible fractionとunified Poisson solve: [`bem_outer_plasma_unified.f90`](../src/physics/outer_plasma/bem_outer_plasma_unified.f90)
- 1D grid演算: [`bem_outer_plasma_grid.f90`](../src/physics/outer_plasma/bem_outer_plasma_grid.f90)
- zero/nonzero fieldの合成: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- explicit outer orbit: [`bem_outer_plasma_orbit.f90`](../src/physics/outer_plasma/bem_outer_plasma_orbit.f90)
- periodic nonzero operator: [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
