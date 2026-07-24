# BEACH 仕様書（現行 Fortran 実装）

## 1. 目的

BEACH は、三角形境界要素上の電荷蓄積とテスト粒子追跡を行うシミュレータです。

- 境界は三角形メッシュ
- 粒子運動は Boris push
- 境界衝突時は吸収して要素へ電荷を堆積
- `batch_count` 単位でバッチ処理し、統計と履歴を更新

計算の主系は Fortran（`src/`, `app/`）で、Python（`beach/`）は後処理・可視化を担います。

## 2. スコープ

### 2.1 実装済み（現行）

- 三角形メッシュ（template / OBJ）
- 静電場（既定は要素重心の point kernel + softening。`field.element_kernel="triangle_p0"` では要素総電荷を三角形上の一定面密度として扱い、direct/treecode/FMM の panel kernel で評価）
- 一様外部磁場 `b0`（任意）
- Boris 法による粒子更新
- 線分 vs 三角形の最初の交差判定
- 衝突時の吸収 + 要素電荷加算（insulator accumulation）
- `surface_model="conductor"` の mesh_id ごとの浮遊導体再配分（`field_bc_mode="free"` のみ）
- ボックス境界条件（open / reflect / periodic）
- 粒子注入モード
  - `volume_seed`
  - `reservoir_face`
  - `photo_raycast`
- チェックポイント再開（`resume`）

### 2.2 現行で未実装・予約

- 表面導電/拡散モデル
- periodic2 と conductor の併用
- 誘電体分極・誘電境界条件（`surface_model="dielectric"` は現状 metadata）
- 反射/二次電子放出を伴う一般化衝突モデル

## 3. データモデル

### 3.1 メッシュ（`mesh_type`）

- 要素数 `nelem`
- 三角形頂点 `v0`, `v1`, `v2`
- 要素重心 `centers`（`center_x/y/z`）
- 要素法線 `normals`
- 要素電荷 `q_elem`
- 要素メッシュ ID `elem_mesh_id`
- 衝突高速化用 AABB / grid キャッシュ

OBJ メッシュ読み込み時、`obj_scale` / `obj_rotation` / `obj_offset` が指定されている場合は scale → rotate → offset の順で頂点座標を変換してからメッシュを初期化する。CRLF 改行および `f v/vt/vn`, `f v//vn` 形式の面行に対応する。

### 3.2 粒子（`particles_soa`）

- 位置 `x(3,n)`、速度 `v(3,n)`
- 電荷 `q(n)`、質量 `m(n)`、重み `w(n)`
- 生存フラグ `alive(n)`

### 3.3 統計（`sim_stats`）

- `processed_particles`
- `absorbed`
- `escaped`
- `escaped_boundary`
- `survived_max_step`
- `batches`
- `last_rel_change`
- フェーズ別計測は `bem_performance_profile` モジュールに分離。`BEACH_PROFILE=1` 環境変数で `performance_profile.csv` を出力

## 4. 1バッチの計算手順

run開始時に、`sim`、`outer_plasma`、`coupling`へ分かれた外部境界設定を、流入写像、通常open面、
z-high輸送、queue ownershipからなる単一の内部契約へ正規化します。粒子loopではこの契約をread-onlyに共有し、
batch injectionもhot loop外で同じresolverを使います。外部場の構築は既存のsnapshot/couplerが担います。

1. delayed outer queueが有効なら、batch開始時刻までにdueとなったreturn/escape eventを取り出す
2. 現在の要素電荷とqueue inventoryに基づいて静電snapshotとouter stateをリフレッシュ
3. refresh済みouter stateを使ってそのバッチの粒子群を生成し、due return粒子を追加
4. 各粒子を `max_step` まで前進（OpenMP スレッド並列）
5. 各ステップで
   - 同一時刻の状態 `x0, v0` から予測中点 `x_mid = x0 + 0.5*v0*dt` を計算
   - 場評価点だけをsolverのprimitive target boxへ写像する（周期軸はwrap、非周期軸はbox面へclamp）。軌道候補座標は写像しない
   - 境界要素電場 `E(x_mid)` を1回評価し、一様外部電場 `sim.e0` を1回加える
   - 中点場を使う Boris push と台形位置更新で `x1, v1` を計算
   - `x1` がbox内部なら `x0 -> x1` のmesh collisionを1回探索
   - box faceへ到達する場合は、full chordのmesh hit parameterと最初のface event parameterを比較して最早順序を決める
   - periodic2のfull-chord queryがbox外区間でrange limitに達した場合は、最初のface eventまでに制限して再照会する
   - reflect/periodic後は残り時間を同じBoris規約で再積分し、そのchordのmesh hitを探索（各local continuationにつき最大8 box event）
   - 衝突時: 粒子を消滅し `q_particle * w_particle` をスレッド別バッファ `dq_thread(elem_idx, tid)` へ加算
   - 残り時間中に9回目のbox eventへ到達し、それ以前にmesh hitがなければ、状態をcommitせず `dt` 縮小を要求する明示的failureとする
6. バッチ終了時に要素電荷差分をコミット: 全スレッドの `dq_thread` を合算し、`photo_emission_dq` を加算した後、MPI allreduce を行い `mesh%q_elem` に反映
7. `rel_change = ||dq|| / max(||q||, q_floor)` を更新
8. 統計と履歴を更新

`photo_raycast` で `deposit_opposite_charge_on_emit=true` の場合は、放出元要素に `-q_particle * w_hit` も加算します。
生成する光電子の重みは常に `w_hit` です。生成後は通常粒子として追跡し、box 内の再吸収、open 面の
`potential_barrier`、または outer particle transfer によって return / escape を決めます。

## 5. 物理モデル

### 5.1 電場

互換既定の `field.element_kernel="point"` は、要素重心点電荷と softening による次式を使います:

- `E(r) = k * Σ_j q_j * (r - c_j) / (|r - c_j|^2 + softening^2)^(3/2)`

ここで `c_j` は要素 `j` の重心です。
`field_solver="treecode"` のときはこの核を遠方で monopole 近似し、近傍は direct 和を使います。  
`field_solver="fmm"` のときは simulator 非依存の Coulomb FMM コアを使い、source octree、optional target tree、Cartesian tensor による multipole/local 展開、近傍 direct 和で電場を評価します。現行 adapter の内部既定次数は 4 です。詳しくは `docs/Algorithms.md` の FMM コア詳細を参照してください。

`field.element_kernel="triangle_p0"` は `field_solver="direct" | "treecode" | "fmm" | "auto"`、`softening=0`、`surface_model="insulator"` に限定します。direct は厳密 free-space panel 和、treecode は全頂点を含む node 半径で near/far を判定して近傍を厳密 panel 核、遠方を monopole で電場・電位とも評価します。FMM は全頂点を含む topology、近傍の厳密 panel 核、三角形 monomial の厳密 P2M を使います。auto は `tree_min_nelem` 未満で direct、以上で FMM を選びます。各 OBJ または template に `surface_side="normal_plus" | "normal_minus" | "outward_closed"` が必要です。`q_elem` は要素総電荷 [C]、面密度は `q_elem/area` です。面上電位は連続、法線電場は `sigma/epsilon0` だけ跳び、重心電位と principal-value 電場を自己項として用います。非対応 solver へ点電荷 fallback はしません。

`sim.field_bc_mode="periodic2"` かつ `field_solver="fmm"` では、`bc_low/high` が `periodic` の2軸を周期軸として扱います（第三軸は開放）。  
近傍画像和は `sim.field_periodic_image_layers = N` に対して各周期軸 `[-N, N]` を評価します。`periodic2` の遠方補正の既定は `field_periodic_far_correction="none"`（`sim` table）です。`auto` は互換用に受理され、`none` に正規化されます。`none` は explicit image shell だけを評価する有限画像近似であり、完全な周期遠方場を与えるものではありません。`m2l_root_oracle` は削除済みで、指定時は reject します。無限周期の非零モードには `cached_kneq0` を使用します。

`field_periodic_far_correction="cached_kneq0"` は production 用の無限 periodic2 非零モード backend です。runtime が加算する有限画像 kernel を `K_shell(N)` とすると、cache は滑らかな full-periodic Ewald residual を root-local operator として保持します。charge refresh 時に source 高さ分布から対称 `k=0` state を一度構築し、各 eval で O(log n) で差し引くため、runtime total は代数的に `K_periodic,k!=0` になります。Ewald all-source 和は cache miss 時の operator generation にだけ使い、particle eval hot path では使いません。物理的な `k=0` は `exclude_k0` provider が場の合成時に一度だけ加えます。`lower_boundary_model="symmetric_vacuum"` は均質真空の無外場境界条件として `E_bottom=-Q/(2 epsilon0 A)`、`E_top=+Q/(2 epsilon0 A)` を選び、interface位置や誘電率を必要としません。`e_bottom_zero` は下側場を0に固定する旧計算再現用境界条件です。外部シースのGauss残差は上側へ入る電束 `Q + epsilon0 A E_bottom` とouter chargeの和で評価します。したがって non-neutral cell も暗黙の charged walls ではなく、この明示的なzero-mode boundary conditionで閉じます。cache fingerprintは周期長、FMM order、画像/Ewald層、source/target topology、softening、generator version、tolerance、real kind、build versionを含みます。MPIではrank 0だけがlock、検証、cache I/O、atomic publishを担当します。cache missのoperator生成はtarget sliceを全rankに分配し、各rank内でproxy RHSをOpenMP並列評価した後、`MPI_Allreduce(SUM)`で全rankに組み立てます。
`tree_theta`/`tree_leaf_max` を未指定の場合は、`periodic2` でも通常の自動推定値を使います。現行実装の推定値は `nelem < 1500` で `theta=0.40`, `leaf_max=12`、`1500 <= nelem < 10000` で `0.50` / `16`、`10000 <= nelem < 50000` で `0.58` / `20`、`50000 <= nelem` で `0.65` / `24` です。

`periodic2.nonzero_mode_backend="panel_spectral_reference"` は、P0 panelのFourier `k!=0`成分、triangle-heightの厳密`k=0`成分、選択したouter responseを合成する小規模correctness referenceです。この経路だけは`field_solver="direct"`を用い、`zero_mode_policy="exclude_k0"`、対応するlower boundary model、x/y periodic・z open、`e0=0`を必須とします。有限image shellや`charged_walls`とは混用しません。interface面の`k!=0`電位・電場減衰、gap、局所平均plasma電荷推定を実測し、設定閾値を超えた場合は`not_applicable`として停止します。外部状態は`outer_update_stride`とともにcheckpointされ、restart後も更新位相を保存します。

運用上、自己整合な外部シースの対応モデルは
`external_boundary.field.model="kinetic_1d"`とする。正規化後の実行時表現は
`outer_plasma.model="kinetic_1d"`である。
外部シースを暗黙に有効化しないため、設定上の既定は引き続き`none`とする。

公開TOMLでは、外部条件を`[external_boundary.field]`、`[external_boundary.particles]`、
`[external_boundary.ordinary_open]`の3責務で指定する。fieldは外部場、particlesはz-high粒子のlifecycleと
reservoir流入、ordinary_openはouterが所有しないopen面を表す。`particles.mode`は
`local_source | same_batch | zhao_queue`、`particles.inflow_model`は
`auto | source_vdf | infinity_barrier`とする。

`same_batch`は`kinetic_1d`の`kinetic_1d_profile_return`と
`electrostatic_1d_instant_return`へ正規化する。

`zhao_queue`は`kinetic_1d + zhao_charge_driven`へ同じkinetic return/transfer対とpersistent queueを加える。
`interface_z`は`sim.box_max[2]`、更新方式は`explicit`へ正規化する。kinetic 1D tracked構成では
`inflow_model="auto"`が同じprofileへ流入ownershipを渡し、local scalar補正との併用を拒否する。

旧`sim.reservoir_potential_model`、`sim.open_boundary_model`、
`[outer_plasma]`、`[coupling]`はcompatibility inputとして読み込むが、`[external_boundary]`との混在は拒否する。
以下の節で使うreturn/transfer/queue名は正規化後の実行時contractを表し、通常利用者が組み立てる公開設定ではない。

`coupling.particle_transfer_mode="electrostatic_1d_instant_return"`では、z-high面を唯一のouter particle interfaceとします。無限遠reservoirの法線VDFはLiouville/flux保存と法線エネルギー保存でinterfaceへ写像します。`kinetic_1d`は`return_model="kinetic_1d_profile_return"`を使います。流入障壁は各batchで先に更新したouter stateの`phi_interface-phi_infinity`から計算し、別の近似へfallbackしません。instant経路の外向き粒子は同じ離散kinetic profileとfar Robin tail上でescape/turning pointを判定し、区分線形電位と指数tailを解析積分して往復時間後に相当するlocal復帰状態を構成します。return位置のx/yには`v_t*tau_outer`を加えて周期wrapします。

`outer_queue_enabled=false`では、outer flightをglobal simulation timeへ加算せず、return粒子の残り`dt`を既存stepperで再積分し、turning粒子のoutward/returned chargeを同じbatchに計上します。これは定常・準定常sheathを消去した簡略化モデルで、定常化後の平均電流と離脱力には適用できますが、UV照射開始などの遅延return currentを含む過渡応答は表しません。`tau_outer/field_evolution_timescale`が上限を超える場合は停止し、`tau_outer/batch_duration >= 1`ではbatch履歴を物理的なreturn-current時間履歴として解釈しません。

`outer_queue_enabled=true`は、`model="kinetic_1d"`、`kinetic_closure="zhao_charge_driven"`、`zhao_branch="auto"`、`return_model="kinetic_1d_profile_return"`、`particle_transfer_mode="electrostatic_1d_instant_return"`の組合せだけで使えます。正の`sim.batch_duration`と`outer_update_stride=1`を要求し、`batch_duration <= max_frozen_field_ratio * field_evolution_timescale`を満たさなければ停止します。各eventでは、`t_due=t_mid+tau_outer`から最初のbatch-start pollまでの量子化遅延`delta_poll`と、batch内crossing時刻のmidpoint近似誤差上限`batch_duration/2`も含め、`tau_outer + delta_poll + batch_duration/2 <= max_frozen_field_ratio * field_evolution_timescale`を課します。超過時はenqueueせず停止します。queueが所有する外部領域はinterfaceから$L=10\lambda_{D,pe}$までの有限control volumeです。turning pointが$L$より手前ならreturn、$L$へ到達すればreservoirへの吸収/escapeとし、queue modeでは$L$外のRobin tailを使ってreturnを判定しません。各batch開始時にdue eventをrank-local queueから取り出し、残った光電子inventoryを面積で割ってZhao closureをpredictor更新します。そのbatchで外向き通過したeventはbatch中央を通過時刻とし、interfaceへのreturnまたは$L$へのescapeまでの`tau_outer`を使って`t_due=t_mid+tau_outer`でqueueへ追加します。surface chargeのcommit後、post-enqueue inventoryでもう一度Zhao closureをcorrector更新し、次batchとcheckpointのcontinuation seedにします。straight runとsplit-resumeは全batchで同じpredictor/corrector列を通ります。return/escapeはdueとなったbatchで計上するため過渡遅延を表しますが、eventはbatch開始時だけreleaseされ、enqueue時のterminal状態を後の場で再積分しません。両modeとも`b0=0`のみを許し、tracked kinetic構成では`inflow_model="auto"`を要求します。

`coupling.steady_start_mode="zhao_floating"`は、定常・準定常研究用の明示的な初期条件です。新規実行の最初batch前に、
設定済みの無限遠reservoirとUV sourceで Zhao 零電流定常根を解き、`phi(infinity)=0`のprofileを構築します。定常根の
interface電場を$E_I$、水平cell面積を$A$として、`symmetric_vacuum`では$Q_{seed}=2\epsilon_0AE_I$、
`e_bottom_zero`では$Q_{seed}=\epsilon_0AE_I$を`coupling.steady_start_mesh_id`の水平平面へ面積比で配ります。選択meshは
同一高さの水平平面でperiodic cell全体を覆い、outer interfaceより下でなければなりません。現在は非重複・無欠損tilingを
構築時に保証できる`mesh.mode="template"`だけを受け入れます。他のmeshは電荷0のままなので、
plane + sphereでplaneを選べばsphereは中性で開始します。この同一profileを初回outer state、reservoir流入補正、instant
return / escapeに使います。後続のbatchは通常のcharge-driven更新に戻り、analytic currentを表面電荷に二重加算しません。
`kinetic_1d` + `zhao_charge_driven` + `kinetic_1d_profile_return` + `electrostatic_1d_instant_return`、
`outer_queue_enabled=false`、`zero_mode_policy="exclude_k0"`、対応するlower boundaryを要求します。新規実行では既存初期電荷を拒否します。
同一configのresumeではcheckpointのmesh電荷と完全なouter stateを復元し、再seedしません。これは未帯電状態からのphysical transientを表せず、
queue過渡closureを置換しません。publication用の定常結果でも、独立な緩和状態または摂動seedに対する感度を確認します。

上記の`sim.batch_duration`は実行時に解決された値を指し、直接指定の代わりに正の
`dt * batch_duration_step`を使ってもよいものとします。

非queue 1D経路はlong flightのfrozen-field上限違反で停止し、Zhao 1D queueはflight、batch-start poll遅延、
midpoint crossing時刻誤差上限の合計へ上限を課し、さらにbatch幅も設定時に制限します。外部磁場を無視する条件は未決なので`b0=0`のみを許可します。

`outer_plasma.model="kinetic_1d"`は、enabledな負・正z-high `reservoir_face` speciesをそれぞれちょうど1つ要求し、無限遠の
electron half-Maxwellian / cold drifting ion VDFとして用い、伸長1D格子上のPoisson方程式を
interface Neumann条件と遠方Robin条件で解きます。初版は単調・無衝突・非磁化分枝に限定し、
ionにはkinetic Bohm入口条件を課します。無限遠電位は`phi(infinity)=0`をゲージとして内部で固定し、
公開・互換入力では変更できません。`photoelectron_density_model="kinetic_mean"`は負電荷
`photo_raycast` speciesの放出fluxからoutgoing/returning平均密度を構成します。解状態は
`converged`、`not_applicable`、`no_physical_solution`、`numerical_failure`を区別し、別のfield modelへ
silent fallbackしません。profileは`outer_plasma_profile.csv`へ保存し、restart時のNewton初期値に使います。
非線形solveは密度モデルの解析微分から構成したbordered-tridiagonal Jacobianを使い、1反復を
格子点数に対して線形時間で解きます。Newton backtrackingが停滞した場合は同じPoisson残差に対する
pseudo-transient stepへ切り替え、前回profileとのinterface field差が大きい場合は適応continuationで
目標fieldへ進みます。pseudo-transient項や中間fieldを収束解として受理せず、最終目標fieldにおける
元の未正則化残差が設定許容値以下の場合だけ`converged`とします。
`photoelectron_density_model="kinetic_mean"`とtracked `kinetic_1d_profile_return`を併用しても、平均密度モデルはouter空間電荷と
current診断だけを供給し、表面へreturn chargeを再加算しません。表面の電荷収支はtracked粒子の放出と
再吸収だけで更新します。

`outer_plasma.kinetic_closure="ambient_linear_debye"`は、surface zero modeが与えるinterface電場$E_I$から
$\phi_I=\lambda_D E_I$、$\phi(z)=\phi_I\exp(-z/\lambda_D)$、
$\rho_{\mathrm{amb}}=-\epsilon_0\phi/\lambda_D^2$を解析的に構成します。
`photoelectron_density_model="none"`を要求し、enabledな`photo_raycast`光電子はtracked放出・再吸収・escapeと
表面電荷更新には残しますが、1D平均密度、平均outer電流、outer space chargeへ含めません。同じprofileを
ambient reservoir流入とz-high return/escapeへ使用します。このclosureは単調なambient線形応答だけを表し、
光電子space charge、非線形sheath、virtual cathode、space-charge-limited / inverse sheath、trapped populationは
適用外です。

`outer_plasma.kinetic_closure="zhao_charge_driven"`は、同じsurface zero modeのinterface電場を境界条件として、
Zhao Type A/B/Cのfree/reflected ambient electron、free/captured photoelectron、cold ion populationを使う選択肢です。
無限遠準中性、Sagdeev積分のfield条件、Type Aではfar-field条件をrootとして解きます。定常Zhaoのzero-current式は
rootに課さず、species別とtotalのcurrent densityを診断として出力します。
`zhao_branch="auto"`または`a`、`b`、`c`を選べます。Type Aの非単調profileでは、流入とescape/returnの両方が離散profile
全体の最初のenergy barrierを走査します。`photoelectron_density_model="kinetic_mean"`との併用を拒否し、
tracked粒子だけが表面電荷を更新します。初版はz-high interfaceを
Zhaoの有効放出面とみなす平面・無衝突・非磁化の準定常modelであり、実放出面からinterfaceまでの一般VDF接続は対象外です。
ambient electron、ion、photoelectronにはそれぞれenabledな負電荷z-high `reservoir_face`、正電荷z-high
`reservoir_face`、負電荷`photo_raycast` speciesをちょうど1つ要求し、0個または複数を拒否します。
electron/ion drift modeは`normal`だけを受理し、photoelectronは`normal_drift_speed=0`、ionは$T_i\le0.1T_e$を要求します。

Zhao profileの物理scaleはphotoelectron温度$T_{pe}$と、$T_{pe}$および$n_{ref}$から導出した$\lambda_{D,pe}$です。
`outer_plasma.debye_length`と`thermal_voltage`はZhao root/profileへ入らず、現時点では`interface_eta_gap`、横方向phi/field比、
local-charge推定などsplit-interface適用性診断のreference inputとして残ります。tracked
`photo_raycast.emit_current_density_a_m2`は、$T_{pe}$をJへ換算して$v_{th,pe}=\sqrt{2T_{pe}/m_{pe}}$とし、有効平面のraw source
$|q_{pe}|n_{ref}\sin(\alpha)v_{th,pe}/(2\sqrt{\pi})$と1%以内で一致しなければなりません。analytic raw currentは
tracked sourceの整合性検査とcurrent-density診断に使いますが、root、surface charge、ledgerへ別途加えず、tracked放出と
再吸収だけが後二者を更新します。過渡population scale $\eta$はphotoelectron密度、無限遠準中性、Sagdeev項を
scaleしますが、current診断のraw photoelectron emission-current項はscaleせず、full tracked sourceのまま使います。
初版は`ray_direction`やrough surfaceからinterfaceへ到達するVDFをZhao outer populationへ
自己無撞着に接続せず、`ray_direction`と`sheath_alpha_deg`はそれぞれ照射rayによる放出面samplingと解析sourceを独立に指定します。
Zhaoの収束確認はprofile grid、有効interface位置、tracked粒子数、時間解像度で行い、汎用の`debye_length`や
`thermal_voltage`をprofile収束parameterとして扱いません。

過渡queueでは、pop後にqueueへ残るphotoelectron macro粒子数を水平面積$A_{xy}$で割ったcolumn targetを使います。
Zhao profileは有限領域$0\le z\le10\lambda_{D,pe}$へ再標本化し、free/captured photoelectron密度の積分がtargetに
一致するpopulation scale $\eta$を解きます。`outer_photoelectron_population_fraction`という出力名でも、$\eta$は確率ではなく
定常reference populationに対するoccupancy scaleです。$\eta=0$から連結するpathを$0\le\eta\le16$で探索し、1を超える
overshootを許します。`[0,1]`へclampせず、targetを無視したfull-population解やdisconnected branchへjump/fallbackしません。
queue modeは`zhao_branch="auto"`を要求し、縮退条件を満たす連続的なA/B等のbranch遷移だけを許します。現在の
bisectionはcolumnが$\eta$とともに単調増加するpathだけを受理し、foldを含む連続pathは未対応として停止します。
targetへ到達する連結・単調解がなければ`no_physical_solution`で停止します。Zhao continuationが失敗した場合は、
MPI rootが`BEACH zhao-continuation` prefixの5行へcall stage、batch、reason/status、target、直前root、拒否candidate、
root jumpをfull-range scientific formatで出力し、flushと全rank barrierの後にfail closedします。このtelemetryは
root選択を変更しません。固定branchのpseudo-arclength atlasは診断APIであり、runtime fallbackではありません。
`diagnose_zhao_ab_degeneracy`は$q=\sqrt{-\phi_m/T_{pe}}$を使い、Type B密度ゼロ端におけるType A準中性curve上の
far-field $q^3$係数、A/Bのinterface field積分差、有限$q$ probeを返します。
`regular_connection_conditions_met`は局所接続の必要条件であり、独立componentの存在や安定性を判定しません。
このA/B診断もruntime root選択には使いません。
`trace_zhao_field_column_homotopy`は、photoelectron sourceが非零のType B/C branchを1本固定し、前後のbatch状態を結ぶ
$E_I(\lambda)$と$N_{pe}(\lambda)$の直線homotopy上で、Zhao残差と有限長column残差を
pseudo-arclength追跡する診断APIです。`target_reached`は$\lambda=1$固定correctorでtargetへ着地した場合だけtrueとなり、
`homotopy_fold_detected`は接線の$\lambda$成分が反転したことを表します。有限density floor、$\eta$範囲、homotopy範囲、
point数への到達は`search_limit`であり、全Zhao manifoldでの不可達を意味しません。非単調Type Aの5座標系と、
columnが恒等的に0となるno-photo Type CはこのAPIの対象外です。このAPIもproduction continuation、branch選択、
fallbackを変更しません。
これはflight delayと有限column inventoryのclosureであり、時間依存Vlasov--Poisson、
outer collision、energy-resolved cloud evolutionではありません。`batch_duration`、tracked粒子数、水平面積、有効interface位置、
profile gridについて収束を確認します。

`sim.field_normalization` で場計算内部の長さを正規化できます。`"si"` が既定で従来どおり、`"box"` は最大 box 幅、`"mesh"` は mesh bbox 最大幅、`"length"` は `sim.field_length_scale` を長さ基準 `L0` とします。direct/treecode/FMM の Coulomb kernel は座標・softening・periodic cell を `L0` で割った無次元距離で評価し、電場で `k_coulomb/L0^2`、電位で `k_coulomb/L0` を掛けて SI に戻します。入力ファイルと出力 CSV は SI 単位のままです。

### 5.2 粒子前進

- Boris 法（`E`, `B`）
- `B` は `sim.b0` の一様場
- public な粒子入力 `x, v` と出力 `x_new, v_new` は同一時刻の状態であり、half-step staggered 状態ではない
- production の空間電場は予測中点 `x_mid = x + 0.5*v*dt` で1回評価する。box crossing候補では評価点だけを周期軸でwrap・非周期軸でbox面へclampし、`sim.e0` はその評価結果へ1回だけ加える
- public pure procedure `boris_update_velocity(v, q, m, dt, e, b, v_new)` が電場 half kick、磁場回転、電場 half kick による速度更新を行う
- 既存の public call `boris_push(x, v, q, m, dt, e, b, x_new, v_new)` は署名を変えず、速度計算を `boris_update_velocity` に委譲し、位置を `x_new = x + 0.5*(v + v_new)*dt` で更新する
- 予測中点の空間電場評価と台形位置更新により candidate kinematics は二次精度であり、一様電場の一定加速度変位は丸め誤差まで解析解と一致する

### 5.3 衝突判定

- 線分 `x0 -> x1` と三角形群の交差（Möller–Trumbore 法）
- AABB によるブロードフェーズ枝刈り。`use_collision_grid=true`（デフォルト）では一様グリッド + 3D-DDA による高速化
- collision query は非有限な線分端点を探索前に拒否する。3D-DDA は grid geometry と各増分の有限性を検証し、反復数を `nx+ny+nz+3` 以下に制限する。進行不能時は命中なしへfallbackせず明示 status でfail closedとする
- 複数候補がある場合は最小 `t`（最初の衝突）を採用
- `sim.field_bc_mode="periodic2"` では、mesh は primitive cell の base element のみ保持しつつ、衝突判定だけ periodic image を implicit に探索する
- periodic2 collision 用 mesh は runtime で canonical unwrapped 形へ平行移動して collision grid を再構築する。raw 頂点は periodic 軸で box 外を含んでよい
- periodic2 hit の `elem_idx` は常に base element index を返し、吸収・堆積先も primitive cell 側へ集約する

### 5.4 ボックス境界

- `open`: 粒子を消滅（`escaped_boundary`）
- `reflect`: 法線成分反転
- `periodic`: 反対側へラップ
- additive な `find_first_boundary_event` は、box 内の始点から候補終点までの最初の交差 fraction と、corner/edge で同時に交差する全 face を bit mask で返す
- additive な `apply_escape_reflect_periodic_event` は同時 face を軸順序に依存せず一括適用し、reflect/periodic 後の位置をbox座標とspanに応じたnormal値のguard幅だけ内側へ置く。これによりzero-valued faceのsubnormalな1 ULP offsetと、直後のevent fractionが0へunderflowする境界chatterを避ける。非有限値、不正な box/face、event/config 不一致は state を変更せず明示 status を返す
- production particle loop はcandidate生成とmesh queryを先行し、box crossing時だけevent resolverへ進む。候補終点がstrictなbox内部なら追加event geometryを行わず、場評価1回・collision query 1回のfast pathとなる
- production particle loop は candidate の位置または速度が非有限なら collision query を呼ばず `particle_step_invalid_boundary` としてfail closedとする
- reflect/periodic crossingだけ残り時間を最大8回再積分する。各eventはmeshとの最早順序を保って処理し、9回目のbox eventまでにmesh hitがなければ `particle_step_multiple_box_events` でfail closedとする。上限なしのevent loopやadaptive substepは行わない
- `open_boundary_model="potential_barrier"` は単一open面のevent位置で局所電位を評価し、補間法線速度の運動エネルギーと `q_particle * (phi_infty - phi_boundary)` を比較する。エネルギー不足では法線速度を反転して残りstepを再積分し、それ以外はescapeとする。複数open faceの同時crossingはfail closedとする
- outer transferが所有するz-highを含む同時eventはface maskのmembershipで判定する。z-highの二次軌道補正後も同時刻と判定できる周期・反射面は横方向作用を先に合成してからouterへ渡す。補正によってz-highと横面の先後・同時関係が元のface maskから変わる場合は、作用順序を推測せずfail closedとする。別open面も同時に含む場合もownerが一意でないためfail closedとする
- z-highの二次補正は、候補終点がz-high外側にあるためchordが検出したcrossingの法線時刻だけを再評価する。候補終点がbox内へ戻る途中の一時的越境は探索せず、x/y面とmesh hitは従来のchord判定を維持する
- instant outer return後の`dt_remaining`でz-highへ再到達した場合も同じ契約へ再dispatchする。1 local stepあたり最大8 external eventとし、9回目は`particle_step_multiple_external_events`でfail closedとする。通常box eventのsoft-discard policyとは混同しない
- legacy `apply_box_boundary` はphoto rayとsource compatibilityのため残す

## 6. 注入モード

### 6.1 `volume_seed`

- 各 species で `npcls_per_step` 個を毎バッチ生成
- 総生成数は `batch_count * Σ npcls_per_step`

### 6.2 `reservoir_face`

- 上流ドリフト Maxwell 分布から流入フラックス `gamma_in` を計算
- `n_macro_expected = gamma_in * A * batch_duration / w_particle`
- 残差繰越つきで `floor` して今バッチのマクロ粒子数を決定
- MPI 実行時も全 rank 合計の期待値と残差を一度だけ更新し、確定した整数個数を rank 間で分配する
- `target_macro_particles_per_batch` 指定時は `w_particle` を自動解決
- `position_jitter_dt=sim.dt` の速度方向ジッタ後、周期軸はprimitive cellへwrapし、非周期軸はbox面へclampして全粒子を有効box内から開始する
- `reservoir_potential_model="infinity_barrier"` 時は注入面平均電位を使って法線速度下限とface到達速度を同じenergy mapで補正する
- 1D outer transfer時は`infinity_barrier`を重ねず、refresh済みouter stateの
  `phi_interface - phi_infinity`（kineticでは全profile上の最大障壁）から流入cutoffを計算する
- 注入面平均の `N x N` 電位評価では、追加の電位評価を行わずに母標準偏差・最小・最大も集計する。Maxwellian reservoirで `abs(q_particle) * phi_std` が `k_B*T + 0.5*m*u_n^2` の10%を超える場合、MPI rootは初回と最終batchに面平均近似の警告を出す

### 6.3 `photo_raycast`

- 指定面から `rays_per_batch` 本を照射
- 各レイは最初の命中要素からのみ放出（命中しなければ放出なし）
- 1ヒット重み:
  - `w_hit = J_perp * A_perp * batch_duration / (|q| * rays_per_batch)`
  - MPI 実行時は `rays_per_batch` の代わりに `global_rays_per_batch`（全 rank 合計）を使用
- 生成粒子の重みには常に `w_hit` を使う
- box内では通常粒子として追跡し、open面では共通の`open_boundary_model`またはouter particle transferを適用する
- `sim.field_bc_mode="periodic2"` で periodic image に命中した場合も、放出位置は primary cell に wrap した hit 座標を使う

## 7. 実行制御と停止条件

現行実装の停止条件は次の通りです。

- 通常実行では `batch_count` バッチ実行し、再開実行では処理済みバッチ数が `batch_count` に達するまで実行
- 各粒子は `max_step` で打ち切り

補足:

- `tol_rel` は `rel_change` の監視値として出力されますが、早期終了には使いません。

## 8. 入出力仕様

### 8.1 出力

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv`（`history_stride > 0`）
- `potential_history.csv`（`output.write_potential_history = true` 時、`history_stride` に従う）
- `mesh_potential.csv`（`output.write_mesh_potential = true` 時）
- `rng_state.txt`
- `macro_residuals.csv`
- `charge_ledger.csv`（粒子種別の電荷収支と粒子数）
- `outer_plasma_profile.csv`（readyなouter stateの条件付きcheckpoint）
- `outer_event_queue.csv`（serialのZhao過渡queue。queue有効時の条件付きcheckpoint）
- `performance_profile.csv`（`BEACH_PROFILE=1` 環境変数設定時）

`mesh_triangles.csv` は要素ごとの `mesh_id` を含み、`mesh_sources.csv` で `mesh_id` と元メッシュ設定を対応付けます。

MPI 実行時はRNGを`rng_state_rankNNNNN.txt`、Zhao過渡queueを`outer_event_queue_rankNNNNN.csv`としてrank別に保存します。
マクロ粒子残差は全rankで共有する状態なので、rootが単一の`macro_residuals.csv`を保存します。

### 8.2 再開（`output.resume = true`）

再開時は既定で `output.dir` から以下を読み込みます。`output.restart_from` を指定した場合は、そのディレクトリから checkpoint を読み込み、新しい出力は引き続き `output.dir` に書きます。

- 必須: `summary.txt`, `charges.csv`, `rng_state.txt`（MPI 時は `rng_state_rankNNNNN.txt`）
- 任意: `macro_residuals.csv`（MPI 時も単一の global ファイル）
- schema v2/v3/v4では電荷収支出力時の`charge_ledger.csv`
- queue有効時はserialの`outer_event_queue.csv`またはMPI全rankの`outer_event_queue_rankNNNNN.csv`

`sim.batch_count` は累積の到達バッチ数です。例えば checkpoint が `batches=100` のとき `batch_count=150` で再開すると、追加で50バッチだけ実行します。`batch_count` が checkpoint の処理済みバッチ数より小さい場合は停止します。MPI 実行時の再開では、前回と同一の `mpi_world_size` が必要です。
`output.resume=true` で必須 checkpoint が存在しない場合は新規実行へフォールバックせず停止します。`summary.txt` の統計値、`charges.csv` の電荷、`macro_residuals.csv` の残差は resume 時に有限性と基本範囲を検証します。
新規出力の`summary.txt`は`checkpoint_schema_version=4`とmodel / ordered mesh / ordered species fingerprintを持ちます。schema v4はschema v3のreadyなouter held stateに加え、Zhao過渡closureのpopulation fraction、column target/residual、queue inventoryを復元します。queue本体はserialの`outer_event_queue.csv`またはMPI rank別ファイルにactive event、terminal outcome、due時刻、`next_event_id`を保存します。queue有効時はschema、rank、world size、完了batchの不一致や欠落ファイルをfail closedで拒否します。schema v3は`outer_plasma_profile.csv`に`z, phi, E, rho`を保存し、outer solverのstatus、反復数、residual、積分電荷、species別電流とともにheld stateを復元します。schema v2は3 fingerprintを照合した上でread-only migrationとして受理しますが、旧3列outer profileは初期値にだけ使い、次のouter refreshでroot solveを強制します。それより古いschemaはPhase 0で実装済みのlegacy point-source modelに限って読み込めます。
`charge_ledger_residual_C` は surface / flight / unresolved stock の差分と外部 flux から作る電荷保存残差です。species 間相殺を避けた `discarded_unresolved_abs_C` は別診断であり、残差が 0 でも max-step discard が物理的に許容されることを意味しません。
旧形式の `macro_residuals_rankNNNNN.csv` が残っている場合は、global 残差との対応が曖昧なため停止します。

## 9. 設計方針

- v1.0 の基本系は insulator accumulation とし、conductor は制限付き拡張として扱う
- 拡張点は維持しつつ、現行利用者向けには実装済み挙動を優先して文書化する
- 設定追加・削除時は `docs/Parameters.md` を同時更新する

## 10. 実行運用（推奨）

通常利用では、`pip install beach-bem` で導入した `beach` コマンドの利用を推奨します。

- 推奨導線: `pip install beach-bem` で導入し、`beach [config.toml]` で実行
- 開発版導線: `python -m pip install "git+https://github.com/Nkzono99/BEACH.git"` で導入
- 開発時の導線: `python -m pip install -e .` + `make` を利用し、必要に応じて `fpm run -- ...` で直接実行

## 11. Python 後処理レイヤ

Python パッケージ `beach` は Fortran 出力を読み込み、後処理・可視化を行います。
詳細な API リファレンスは `docs/PythonPostprocessAPI.md` を参照してください。

### 11.1 電位再構成・電場計算

- `compute_potential_mesh` / `compute_potential_points`: 重心点電荷近似による電位の Python 側再構成
- `compute_electric_field_points`: 任意 3D 点での電場ベクトル。`E(r) = K * sum_j q_j * (r - r_j) / |r - r_j|^3` の direct 和で計算

### 11.2 Coulomb 力と periodic2 対応

- `calc_coulomb`: メッシュグループ間の Coulomb 力/トルク計算
- `periodic2` パラメータ: Fortran 側の `sim.field_bc_mode="periodic2"` に対応し、Python 側でも 2 軸周期境界での画像シェル和をサポート。ソース電荷を `ix in [-nimg, nimg], iy in [-nimg, nimg]` でシフトして直接和を取る
- `periodic2=None` では出力ディレクトリ近傍の `beach.toml` から `field_bc_mode` を自動判定

### 11.3 電気力線追跡

- `trace_field_lines`: シード点から RK4 積分で電気力線を追跡。`direction` で順方向 / 逆方向 / 両方向を選択可能
- `plot_field_lines_3d`: 力線を 3D 描画し、三角形メッシュを電荷密度で着色してオーバーレイ

### 11.4 Python 側の制限

- Python 側の電場・電位計算は要素重心の点電荷による direct 和のみ（Fortran 側の treecode / fmm は使用しない）
- pure Python の電位再構成は `cached_kneq0` operator を再現せず、explicit image shell のみ。native `FieldKernel` は `cached_kneq0` を利用できる
