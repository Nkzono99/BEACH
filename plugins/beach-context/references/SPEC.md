# BEACH 仕様書（現行 Fortran 実装）

## 1. 目的

BEACH は、三角形境界要素上の電荷蓄積とテスト粒子追跡を行うシミュレータです。

- 境界は三角形メッシュ
- 粒子運動は Boris push
- 境界衝突時は吸収して要素へ電荷を堆積
- `batch_count` 個の accepted batch を処理し、統計と履歴を更新

計算の主系は Fortran（`src/`, `app/`）で、Python（`beach/`）は後処理・可視化を担います。

## 2. スコープ

### 2.1 実装済み（現行）

- 三角形メッシュ（template / OBJ）
- 静電場（要素総電荷を三角形上の一定面密度として扱う `triangle_p0` panel kernel。direct/treecode/FMM で評価）
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
- `simulated_time_s`
- `adaptive_nonzero_mode_rejected_trials`
- `adaptive_nonzero_mode_last_batch_duration_s`
- `adaptive_nonzero_mode_last_potential_step_V`
- `adaptive_nonzero_mode_omp_threads`
- `last_rel_change`
- フェーズ別計測は `bem_performance_profile` モジュールに分離。`BEACH_PROFILE=1` 環境変数で `performance_profile.csv` を出力

## 4. 1バッチの計算手順

run開始時に、`sim`、`outer_plasma`、`coupling`へ分かれた外部境界設定を、流入写像、通常open面、
z-high輸送、queue ownershipからなる単一の内部契約へ正規化します。粒子loopではこの契約をread-onlyに共有し、
batch injectionもhot loop外で同じresolverを使います。外部場の構築は既存のsnapshot/couplerが担います。

通常は解決後の`sim.batch_duration`を1 batchの幅とします。
`periodic2.max_nonzero_mode_potential_step > 0`では、この値を最大幅$h_0$とし、各accepted batchで
$h_0,h_0/2,h_0/4,\ldots$の固定ladderを先頭から試します。以下の手順3から電荷候補の作成までが
1 trialであり、受理されたtrialだけが手順6以降をcommitします。

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
6. バッチ終了時に要素電荷候補を作成: 全スレッドの `dq_thread` を合算し、`photo_emission_dq` と
   `implicit_mean` の候補更新を含めてMPI間で確定する。適応進行では、全panel重心で
   $\max_j|[P_{k\ne0}(\mathbf q_{\mathrm{candidate}}-\mathbf q_{\mathrm{current}})]_j|$を評価し、
   `max_nonzero_mode_potential_step`以下なら受理する。超えるtrialはcommitせず、短い次のladder幅で再試行する
7. 受理した候補だけを `mesh%q_elem` に反映し、`rel_change = ||dq|| / max(||q||, q_floor)` を更新
8. 受理した幅を`simulated_time_s`へ加え、accepted batch数、統計、履歴、charge ledgerを更新

適応trialを棄却すると、RNG、macro粒子数残差、outer state、`implicit_mean`のmean/outer transactionを
trial開始前へ完全にrollbackします。棄却trialの粒子outcome、統計、履歴、charge ledgerは外部へ出しません。
`sim.batch_count`と`sim_stats%batches`はaccepted batch数を表します。
適応trialの加算順を再現するため、再開時はcheckpointの
`adaptive_nonzero_mode_omp_threads`と同じ実OpenMP team sizeを使います。適応区間ではdynamic teamを
無効化し、全MPI rankでteam sizeが一致しない実行と、checkpointとの不一致を拒否します。
`implicit_mean`のZhao更新が数値的に有効な候補を作れない場合も、設定不正ではなく短い時間幅で
回復可能なstatusに限って同じladder上で棄却・再試行します。標準出力は$k\ne0$上限超過と
この$k=0$ closure recoveryを別のreasonとして記録します。

`photo_raycast` で `deposit_opposite_charge_on_emit=true` の場合は、放出元要素に `-q_particle * w_hit` も加算します。
生成する光電子の重みは常に `w_hit` です。生成後は通常粒子として追跡し、box 内の再吸収、open 面の
`potential_barrier`、または outer particle transfer によって return / escape を決めます。

## 5. 物理モデル

### 5.1 電場

要素 source model は `triangle_p0` に固定され、設定では選択しません。公開設定に
`[field]` table と `sim.softening` はありません。旧設定を残した入力は alias や暗黙変換を行わず、
unknown table / key として停止します。

`triangle_p0` は各要素の総電荷 `q_elem` [C] を面積で割り、三角形上の一定面密度として扱います。
`field_solver="direct" | "treecode" | "fmm" | "auto"` と
すべてのsurface modelで共通して使います。direct は厳密 free-space panel 和、treecode は全頂点を含む
node 半径で near/far を判定して近傍を厳密 panel 核、遠方を monopole で電場・電位とも評価します。
FMM は全頂点を含む topology、近傍の厳密 panel 核、三角形 monomial の厳密 P2M を使います。
auto は `tree_min_nelem` 未満で direct、以上で FMM を選びます。各 OBJ または有効な template に
`surface_side="normal_plus" | "normal_minus" | "outward_closed"` が必要です。面上電位は連続、
法線電場は `sigma/epsilon0` だけ跳び、重心電位と principal-value 電場を自己項として用います。

`sim.field_bc_mode="periodic2"` かつ `field_solver="fmm"` では、`bc_low/high` が `periodic` の2軸を周期軸として扱います（第三軸は開放）。  
近傍画像和は `sim.field_periodic_image_layers = N` に対して各周期軸 `[-N, N]` を評価します。`periodic2` の遠方補正の既定は `field_periodic_far_correction="none"`（`sim` table）です。`auto` は互換用に受理され、`none` に正規化されます。`none` は explicit image shell だけを評価する有限画像近似であり、完全な周期遠方場を与えるものではありません。`m2l_root_oracle` は削除済みで、指定時は reject します。無限周期の非零モードには `cached_kneq0` を使用します。

`field_periodic_far_correction="cached_kneq0"` は production 用の無限 periodic2 非零モード backend です。runtime が加算する有限画像 kernel を `K_shell(N)` とすると、cache は滑らかな full-periodic Ewald residual を root-local operator として保持します。charge refresh 時に source 高さ分布から対称 `k=0` state を一度構築し、各 eval で O(log n) で差し引くため、runtime total は代数的に `K_periodic,k!=0` になります。Ewald all-source 和は cache miss 時の operator generation にだけ使い、particle eval hot path では使いません。物理的な `k=0` は `exclude_k0` provider が場の合成時に一度だけ加えます。`lower_boundary_model="symmetric_vacuum"` は均質真空の無外場境界条件として `E_bottom=-Q/(2 epsilon0 A)`、`E_top=+Q/(2 epsilon0 A)` を選び、interface位置や誘電率を必要としません。`e_bottom_zero` は下側場を0に固定する旧計算再現用境界条件です。外部シースのGauss残差は上側へ入る電束 `Q + epsilon0 A E_bottom` とouter chargeの和で評価します。したがって non-neutral cell も暗黙の charged walls ではなく、この明示的なzero-mode boundary conditionで閉じます。cache fingerprintは周期長、FMM order、画像/Ewald層、source/target topology、generator version、tolerance、real kind、build versionを含みます。MPIではrank 0だけがlock、検証、cache I/O、atomic publishを担当します。cache missのoperator生成はtarget sliceを全rankに分配し、各rank内でproxy RHSをOpenMP並列評価した後、`MPI_Allreduce(SUM)`で全rankに組み立てます。
`tree_theta`/`tree_leaf_max` を未指定の場合は、`periodic2` でも通常の自動推定値を使います。現行実装の推定値は `nelem < 1500` で `theta=0.40`, `leaf_max=12`、`1500 <= nelem < 10000` で `0.50` / `16`、`10000 <= nelem < 50000` で `0.58` / `20`、`50000 <= nelem` で `0.65` / `24` です。

`periodic2.max_nonzero_mode_potential_step`は[V]単位の非負値で、省略または`0`では無効です。
正値は`nonzero_mode_backend="cached_kneq0"`、正の`sim.batch_duration`、time-scaledな
`reservoir_face` / `photo_raycast` sourceを要求し、`volume_seed`とouter event queueを拒否します。
`reservoir_face`はtrial幅とともにmacro粒子電荷も縮小できる
`target_macro_particles_per_batch`方式に限り、固定`w_particle`は拒否します。
explicit SW/UV更新と、`kinetic_1d + same_batch` PEが内部で使う`implicit_mean`更新に対応します。
判定量は候補電荷差のcached $k\ne0$電位を全panel重心で評価した最大絶対値です。これはbatch内の
frozen-field局所電位変化を制限するtrust boundであり、局所打切り誤差を推定または保証しません。
時間幅収束は、この上限を半分にした計算との比較で確認します。
適応時はOpenMPの粒子割当を同じthread数に対して決定的なstatic partitionとし、受理境界に
丸め誤差幅を設けます。thread数を変えたrun間のbitwise同一性は契約しません。

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
`interface_z`は`sim.box_max[2]`へ正規化する。通常の更新方式は`explicit`だが、
`kinetic_1d + same_batch`にenabledな負電荷`photo_raycast` speciesがあり、closureが
`ambient_linear_debye`または`zhao_charge_driven`の場合は、公開keyを増やさず内部の
`implicit_mean`へ正規化する。後者は`steady_start_mode="zhao_floating"`でcommit済みA/B/C branchを
初期化することも要求する。kinetic 1D tracked構成では
`inflow_model="auto"`が同じprofileへ流入ownershipを渡し、local scalar補正との併用を拒否する。

旧`sim.reservoir_potential_model`、`sim.open_boundary_model`、
`[outer_plasma]`、`[coupling]`はcompatibility inputとして読み込むが、`[external_boundary]`との混在は拒否する。
以下の節で使うreturn/transfer/queue名は正規化後の実行時contractを表し、通常利用者が組み立てる公開設定ではない。

`coupling.particle_transfer_mode="electrostatic_1d_instant_return"`では、z-high面を唯一のouter particle interfaceとします。無限遠reservoirの法線VDFはLiouville/flux保存と法線エネルギー保存でinterfaceへ写像します。`kinetic_1d`は`return_model="kinetic_1d_profile_return"`を使います。流入障壁は各batchで先に更新したouter stateの`phi_interface-phi_infinity`から計算し、別の近似へfallbackしません。instant経路の外向き粒子は同じ離散kinetic profileとfar Robin tail上でescape/turning pointを判定し、区分線形電位と指数tailを解析積分して往復時間後に相当するlocal復帰状態を構成します。return位置のx/yには`v_t*tau_outer`を加えて周期wrapします。

`outer_queue_enabled=false`では、outer flightをglobal simulation timeへ加算せず、return粒子の残り`dt`を既存stepperで再積分し、turning粒子のoutward/returned chargeを同じbatchに計上します。これは定常・準定常sheathを消去した簡略化モデルで、定常化後の平均電流と離脱力には適用できますが、UV照射開始などの遅延return currentを含む過渡応答は表しません。通常の`same_batch`粒子では`tau_outer/field_evolution_timescale`が上限を超える場合に停止し、`tau_outer/batch_duration >= 1`ではbatch履歴を物理的なreturn-current時間履歴として解釈しません。

`implicit_mean`には2つの兄弟経路がある。`ambient_linear_debye`は解析Maxwellianの平均escape率を使う
backward Euler解でcell総電荷とinterface電位を確定し、粒子標本から平均escape率を解き直さない。
`zhao_charge_driven`は、実測したinterface法線energy分布と非単調profile全体のbarrierを使う一般化電荷根を解く。
どちらも確定したprofile上でz-high crossingを追跡してreturn軌道と再吸収位置の準定常shadow標本を作るが、
平均escape / return総量の決め方とray weightの規約は共有しない。outer flightをglobal simulation timeへ加えず、
shadowの`tau_outer`とfrozen-field ratioは診断するがratio超過を停止条件にしない。他speciesのnonqueue instant
returnは上記のfail-closed contractを保つ。

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
`dt * batch_duration_step`を使ってもよいものとします。適応的な非零モード進行ではこの値は最大trial幅であり、
各accepted batchの物理時間は実際に受理したladder幅です。

非queue 1D経路の通常粒子はlong flightのfrozen-field上限違反で停止します。`implicit_mean`のdeferred
`photo_raycast` shadowだけはflight timeとratioを診断へ残したまま、この上限を停止条件にしません。
Zhao 1D queueはflight、batch-start poll遅延、midpoint crossing時刻誤差上限の合計へ上限を課し、
さらにbatch幅も設定時に制限します。外部磁場を無視する条件は未決なので`b0=0`のみを許可します。

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
公開TOMLでは`photoelectron_density_model` keyを受理せず、facadeが内部値を`none`へ解決します。
`"none"`の明記も拒否します。enabledな`photo_raycast`光電子は3D領域内のtracked放出・再吸収と
interfaceへの外向き通過に残しますが、1D平均密度とouter space chargeへ含めません。同じprofileを
ambient reservoir流入と、光電子を含むz-high return/escapeへ使用します。このclosureは単調なambient線形応答だけを表し、
光電子space charge、非線形sheath、virtual cathode、space-charge-limited / inverse sheath、trapped populationは
適用外です。

このclosureを`same_batch`およびenabledな負電荷`photo_raycast` speciesと組み合わせると、
loaderは内部の`coupling.update_mode="implicit_mean"`を自動選択します。利用者はこのraw値を入力しません。
局所的なtracked放出・interface到達前の再吸収・ambient吸収による要素間の差は、既存の3D batch traceと
element depositで進めます。interfaceを横切った光電子の実測外向き電流
$J_{pe,\mathrm{out}}^n$を1D平均領域のsourceとします。最初のtraceでstageされた全surface charge deltaから
$J_{\mathrm{pending}}^n$を測り、

$$
J_{\mathrm{other}}^n
=J_{\mathrm{pending}}^n-J_{e,\mathrm{tracked}}^n-J_{i,\mathrm{tracked}}^n-J_{pe,\mathrm{out}}^n
$$

と定義します。これにより、追加species、他面からの流入、別open面へのescape、未解決粒子のcounterchargeを
平均更新から落としません。平面平均電位とcell総電荷は、1 batchの時間幅$\Delta t$について

$$
C_A\frac{\phi_I^{n+1}-\phi_I^n}{\Delta t}
=J_{e,\mathrm{tracked}}^n+J_{i,\mathrm{tracked}}^n
+J_{\mathrm{other}}^n
+J_{pe,\mathrm{out}}^n f_{\mathrm{esc}}(\phi_I^{n+1})+J_{\mathrm{ext}}
$$

をbackward Eulerで解きます。$C_A=\epsilon_0/\lambda_D$は`e_bottom_zero`、
$C_A=2\epsilon_0/\lambda_D$は`symmetric_vacuum`です。ambient電流は、そのbatchで実際に表面へ
吸収された一意なz-high `reservoir_face`電子・イオン種から測定し、それ以外の実測surface changeは
$J_{\mathrm{other}}^n$に残します。$J_{pe,\mathrm{out}}^n$は、そのbatchで実際にinterfaceを通過した
tracked光電子から測定します。`emit_current_density_a_m2`はtracked 3D放出の粒子weightを決めますが、
その設定値を独立な平均sourceとして再利用しません。設定温度が定めるMaxwellian
barrier通過率$f_{\mathrm{esc}}$を$\phi_I^{n+1}$で評価します。このscalar backward Euler解が
$Q^{n+1}$、$\phi_I^{n+1}$、連続Maxwellianのescape / return総量を決めるsource of truthです。
このescape率は放出driftを含まないため、対象の負電荷`photo_raycast` speciesは
`normal_drift_speed=0`を要求し、非零値はfail closedで拒否します。
有限個のparticle sampleでこの平均総量を置き換えません。最初のtraceでstageした全要素deltaを
$\Delta Q_{e,\mathrm{pending}}$とすると、解析return電荷の絶対値は

$$
Q_{\mathrm{pending}}^{n+1}
=Q^n+\sum_e\Delta Q_{e,\mathrm{pending}},\qquad
R_{\mathrm{analytic}}=Q_{\mathrm{pending}}^{n+1}-Q^{n+1}
$$

です。外部追加電流が0なら
$R_{\mathrm{analytic}}=A\Delta t\,J_{pe,\mathrm{out}}^n[1-f_{\mathrm{esc}}(\phi_I^{n+1})]$
に一致します。$R_{\mathrm{analytic}}$が0未満またはtracked outward source総量を超える場合は停止します。

このscalar解は、速い平面平均帯電を同じbatch内で閉じます。解析return率を放出元ごとのoutward chargeへ掛けて
平均候補を作り、batch-startの非零モードを固定したまま全要素候補電荷でmean-only snapshotを1回更新します。
その後、各実crossingの位置・速度・macro weightを保ち、既存の`kinetic_1d_profile_return`写像で離散profileと
far Robin tailをenergy-resolvedに1回追跡します。turning後のlocal軌道がz-highを再びcrossした場合も同じdriverを使い、
terminalなsurface absorptionまたはescapeへ到達するまで追跡します。ただし、このray分類は軌道と着地点の分布標本であり、
平均return総量の再推定には使いません。

terminalに再吸収されたrayから、放出元elementと実hit elementの未正規化分布を作ります。そのsample return総量を
$R_{\mathrm{sample}}$とすると、$R_{\mathrm{analytic}}>0$では全return sampleへ同じ

$$
\alpha=\frac{R_{\mathrm{analytic}}}{R_{\mathrm{sample}}}
$$

を掛けます。$\alpha$は確率ではなくsampling weightなので、$[0,1]$へclampせず1を超え得ます。
負電荷macro粒子のchargeを$q_p<0$、sample weightを$w$、放出元を$s$、着地点を$d$とすると、
正規化後の帰還分は

$$
\Delta Q_s^{\mathrm{source}}=-\alpha q_pw,\qquad
\Delta Q_d^{\mathrm{return}}=\alpha q_pw,\qquad
\sum_e\left(\Delta Q_e^{\mathrm{source}}+\Delta Q_e^{\mathrm{return}}\right)=0
$$

です。放出元counterchargeと実着地点depositを1つの零和transactionとしてcommitし、escape分のsource counterchargeは
表面総電荷に残します。$R_{\mathrm{analytic}}=0$ではreturn weightを0とし、
$R_{\mathrm{analytic}}>0$なのにterminal再吸収sampleがない場合、source/returnが許容値内で相殺しない場合、
または許可されたlocal追跡内にterminal outcomeへ到達しない場合はfail closedで停止します。
outer return後にz-high以外のlocal open faceへ出たshadowも、解析的な上向きescapeへ再分類せずfail closedで停止します。
別要素への残差配分や粒子分類を使ったscalar解の再反復は行いません。

ここで陰的に分離するのは、interface fieldを決めるcell総電荷のscalar modeです。残りの3D element分布には
$k_\parallel\ne0$だけでなく、総和ゼロの高さ依存した平面平均双極子も含まれます。したがってこの実装を
数学的な全$k=0$成分と全$k\ne0$成分の完全分離とは解釈しません。
同じbatch内ではscalar解が速い平均総電荷を決め、正規化した実hit depositは総電荷を変えずに局所分布だけを再配置します。
この局所再分配が非零モードへ反映されるのは通常のbatch commit後です。inner fixed-pointを追加せず、
平均帯電と局所patchの異なるtimescaleをこの2段階で分けます。

この経路は正の `batch_duration`、`outer_update_stride=1`、ちょうど 1 つの負電荷 `photo_raycast` species、
そのspeciesの`normal_drift_speed=0`、および`deposit_opposite_charge_on_emit=true`を要求し、
outer queueとの併用を拒否します。
公開 TOML の `photoelectron_density_model` は省略し、内部で `none` に解決します。update mode や return kernel の
公開 key は追加せず、OBJ / templateの別や特定meshを要求しません。

同じbatch末尾でouter profileを強制更新します。charge ledgerでscalar closureのanalytic escape / return総量へ
置換するのは、z-highを外向きにcrossしてdeferredされた光電子成分だけです。interface到達前のlocal absorption、
z-lowなど別のopen面からのescape、unresolved discardはtracked値を保持します。したがって
`escaped_to_infinity`はtrackedな別open面escapeとanalytic z-high escapeの和、
`absorbed_on_surface`はtrackedなlocal absorptionとanalytic z-high return後absorptionの和です。
ray terminal分類はdeferred成分の
位置分布の標本であり、charge総量のsource of truthではありません。
`interface_outward_gross` / `interface_returned_gross`はterminal量ではなくinterface crossingのgross量です。
ray crossingはanalytic terminal総量へ正規化したweightで集計し、
z-high deferred成分の符号付きanalytic escape電荷を
$Q_{\mathrm{esc,z-high}}^{\mathrm{analytic}}$とすると、
`interface_returned_gross` = `interface_outward_gross` - $Q_{\mathrm{esc,z-high}}^{\mathrm{analytic}}$を保ちます。
別open面からのtracked escapeも含むterminal総量`escaped_to_infinity`を、この式へそのまま代入しません。
return後に再びz-highを横切るrayはoutwardへ、さらにturningして戻るrayはreturnedへ追加されます。
標準出力の`BEACH implicit-mean`行は、零和確認の`transaction_residual_C`、scalar解の
`mean_solver_iterations`、標本上の`sample_escape_fraction`、解析return総量への`return_weight_scale`を記録します。

`summary.txt`は`implicit_mean`について、今回の実行で最後に完了したbatchのanalytic weight適用後のreturn outer
excursionだけを2つの非負診断値へ集約します。restart後に1 batchも進めないno-op実行では、この2 keyはsolver
stateではなくrun-localな派生診断なので出力しません。return excursion $j$の正の電荷magnitudeを$W_j>0$、
outer往復時間を$\tau_j$、
水平面積を$A$、batch幅を$\Delta t$とすると、

$$
\bar{\tau}_{\mathrm{ret},Q}
=\frac{\sum_jW_j\tau_j}{\sum_jW_j},
\qquad
\widehat{\sigma}_{\mathrm{PE,ret}}
=\frac{\sum_jW_j\tau_j}{A\Delta t}
$$

です。`implicit_mean_last_returned_outer_flight_time_mean_s`は$\bar{\tau}_{\mathrm{ret},Q}$、
`implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2`はLittleの法則による
準定常shadow推定$\widehat{\sigma}_{\mathrm{PE,ret}}$を記録します。return excursionがなければ両方を0とします。
$W_j$はanalytic normalization後の正のcharge magnitudeであり、後者は実queueやcharge ledgerに保持された
outer cloud stockではありません。
`escaped_to_infinity` outcome自体には滞在時間を加えませんが、そのrayが最終的にescapeする前に完了したreturn
excursionは集計します。これらは全batchの累積値または最大値ではありません。

`outer_integrated_charge_per_area_C_m2`が非零の場合、

$$
\chi_{\mathrm{PE,shadow}}
=\frac{\widehat{\sigma}_{\mathrm{PE,ret}}}
       {\left|\texttt{outer\_integrated\_charge\_per\_area\_C\_m2}\right|}
$$

は、省略したreturning-photoelectron shadow columnのmagnitudeを1D outer profileの積分電荷magnitudeと比較する
解釈用の比です。これを設定key、停止条件、または組み込みの受理閾値にはしません。

この`ambient_linear_debye`経路は、deferred光電子を準定常shadowとしてenergy分類し、outer flight timeと
frozen-field ratioを診断する一方で、flightをglobal simulation timeへ加えずratio超過でも停止しない
same-batch event transactionです。batch間に残る光電子cloud inventory、光電子space chargeを含む非線形1D密度、
virtual cathode、trapped populationは扱いません。UV turn-onの遅延return currentを時間発展させるmodelではなく、
その用途にはこのclosureへ対応するdelayed inventory / queueを別途設計する必要があります。これらが必要な過渡・
強光電子条件では、このclosureの結果を自己無撞着な光電子シースまたは物理的なreturn-current時間履歴と解釈しません。

ambient electron/ion reservoirを持たないUV-only構成には、この無限遠closureを適用しません。
`external_boundary.field.model="none"`でz-highをescapeさせるUV-only計算は、局所放出・再吸収を調べる
有限box過渡controlであり、無限遠で閉じた準中性シースまたは定常解とは解釈しません。

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

photoemitting Zhaoを`same_batch`で使う強PE経路は、内部の`implicit_mean`を選び、
`zhao_floating`でcommitしたA/B/C branchを前batch状態として要求します。最初の3D traceはbatch-startの場で
光電子を進め、z-highを最初に外向き通過したray $r$について

$$
{\cal E}_{n,r}=\frac{1}{2}m_{pe}v_{z,r}^2,\qquad w_r=|\Delta q_r|
$$

を採取します。各rankの$({\cal E}_{n,r},w_r)$は`MPI_Gatherv`でrootへ集めます。rootはenergyを降順に
stable sortし、浮動小数点比較で完全に等しいenergyだけを1 groupへまとめ、charge weightを加算して累積電荷$C_k$を作ります。
設定bin、histogram補間、smoothingは使いません。このためrootの一時memoryと通信量はinterfaceを通過した
光電子ray数に比例し、現行実装は分散CDFではありません。

最初のtraceでstageされた全surface chargeから全interface光電子source chargeを除いた値を$Q_{\rm base}$とし、
候補charge $Q$に対するZhao profile barrierを

$$
B(Q)=\max\left[
0,\ \max_z q_{pe}\{\phi(z;Q)-\phi_I(Q)\},\
q_{pe}\{\phi_\infty-\phi_I(Q)\}
\right]
$$

と定義します。$q_{pe}<0$であり、最大値は有限gridの全profileとinfinity endpointを走査するため、
Type Aのvirtual cathodeに対応するpotential minimumも含みます。実測CDFが与えるescape charge
$C_{\rm esc}(B)$は${\cal E}_{n,r}>B$のgroupを全量含み、turning-point equalityはreturn側とします。
barrierが1つのenergy groupを横切る場合だけ、そのgroupへ$0\le\theta\le1$のfractional escape weightを
割り当てます。rootは

$$
Q=Q_{\rm base}+C_{\rm esc}[B(Q)]
$$

を満たす一意な一般化根です。interface fieldは`e_bottom_zero`で
$Q=\epsilon_0A E_I$、`symmetric_vacuum`で$Q=2\epsilon_0A E_I$として各候補へ写します。

energy groupを${\cal E}_1>\cdots>{\cal E}_M$、累積chargeを$C_k$とすると、
$Q_k=Q_{\rm base}+C_k$です。solverは最終source上の共通rootを作り、このrootから$Q_0$と$Q_M$へ
connected pathを1本ずつ追跡します。この2本が全候補charge区間を覆い、各adaptive accepted pointの
barrier勾配とsecantを検査します。次のorder predicateは、検査した$B(Q)$の非減少性とenergyの降順性から
falseからtrueへ1回だけ変化します。

$$
P_k=
\begin{cases}
[B(Q_k)\ge{\cal E}_{k+1}], & 0\le k<M,\\
[B(Q_M)\ge{\cal E}_M], & k=M.
\end{cases}
$$

solverはfirst-true indexを二分探索するため、pure root選択に必要なconnected candidate solve数は$O(\log M)$です。
$P_0$がtrueなら$k=0$のall-return、$P_M$がfalseなら$k=M$のall-escapeです。
内部のfirst-true index $k$では、$B(Q_k)<{\cal E}_k$ならpure root、
$B(Q_k)\ge{\cal E}_k$なら$[Q_{k-1},Q_k]$でfractional marginal rootを解きます。
marginal二分では$B(Q_{\rm low})<{\cal E}_k\le B(Q_{\rm high})$を維持し、turning-point equalityを
return側に置くため最終的にupper endpointを採用します。

各candidateのZhao rootは、共通rootを始点とする連結parameter path

$$
\mathbf G\!\left(\mathbf y;E_I(\lambda),n_{pe,0}(\lambda)\right)=\mathbf 0,\qquad
E_I(\lambda)=(1-\lambda)E_{I,0}+\lambda E_{I,1},\qquad
n_{pe,0}(\lambda)=(1-\lambda)n_{pe,0}^{(0)}+\lambda n_{pe,0}^{(1)}
$$

上で追跡します。$\mathbf y$は固定したA/B/C branchのlog-encoded root変数です。Type A/Bは$E_I>0$、
Type Cは$E_I<0$を要求し、このchartでregularでない$E_I=0$も拒否します。adaptive pseudo-arclengthは、
局所補正距離、root jump、接線方向、Jacobian rankを検査します。$\lambda$接線の消失・反転または
固定parameter Jacobianのrank低下をfoldとしてtarget前にfail closedし、targetは$\lambda=1$のevent correctorで
着地します。これは近傍rootへのjumpを制限する有限精度のlocal guardであり、任意に近い別sheetが存在しないことの
数理的証明ではありません。この内部経路に対応する新しいTOML keyはありません。

最終sourceを固定した候補pathでは、各adaptive accepted pointで$B$のprescribed chargeに対する接線勾配も評価し、
負の勾配または逆向きのbarrier secantを検出すると停止します。endpoint path、二分探索のmidpoint、
marginal二分でもcharge順のbarrier bracketを検査します。これはaccepted point間をinterval arithmeticで
包含する連続区間上の数理的証明ではありません。最終候補では実測CDFからescape chargeを再計算し、
$Q-Q_{\rm base}-C_{\rm esc}=0$を要求します。rank低下、fold、barrier勾配反転、order predicateのbracket不整合、
またはmarginal energyを挟めない場合はfail closedとします。

surfaceの`photo_raycast.emit_current_density_a_m2`は放出ray数とmacro weightを決めます。一方、各batchの
Zhao source amplitudeは実測interface outward current

$$
J_{pe,I}^{\rm meas}=\frac{\sum_r w_r}{A\Delta t},\qquad
s_{\rm eff}=\frac{J_{pe,I}^{\rm meas}}
{|q_{pe}|n_{\rm ref}\sin(\alpha)\sqrt{2T_{pe}/m_{pe}}/(2\sqrt{\pi})}
$$

から解き、設定した`photoelectron_source_scale`を同じbatchのinterface sourceとして再利用しません。
設定値は`zhao_floating`の初期branch anchorを構成し、実測値はその後のresolved source scaleとして
outer stateとcheckpointへcommitされます。したがってsurface emissionとinterface到達後のsource normalizationは
別の責務です。

実測sourceが前batch値と異なる場合は、上記の連結pathで前batch rootから最終$s_{\rm eff}$上の共通anchorへ進みます。
その後の全charge候補はsourceを最終値に固定し、共通anchorから同じcommit済みA/B/C root sheet上を追跡します。
branch labelの変更、bootstrap branchへの後退、disconnected rootへのjump、別closureへのfallbackは許しません。
連続stepでは$\phi_0$、$\phi_{\min}$、無限遠electron密度の変化を$0.25$以下の無次元幅に制限します。
さらに1つのfrozen cohortを同じbatch内で使える範囲として、

$$
\frac{|q_e\Delta\phi_I|}{T_e}\le0.25,\quad
\frac{|\Delta B_{e,\mathrm{in}}|}{T_e}\le0.25,\quad
\frac{|q_i\Delta\phi_I|}{E_{i,n}}\le0.25,\quad
\left|\log\frac{n_{e,\infty}^{\mathrm{new}}}{n_{e,\infty}^{\mathrm{old}}}\right|\le0.25,\quad
|\log(s_{\rm eff}/s_{\rm old})|\le0.25,\quad
\frac{\lambda_{D,pe}|\Delta E_I|}{T_{pe}/|q_{pe}|}\le0.25,\quad
\frac{|\Delta B|}{T_{pe}}\le0.25
$$

を要求します。違反時はrayを新しい場で取り直したり旧profileへfallbackしたりせず停止します。
ここで$B_{e,\mathrm{in}}=\max(0,q_e\phi_I,q_e\phi_{\min})$、
$E_{i,n}=\max(T_i,m_i u_{i,n}^2/2)$です。最初の4項は、無限遠からinterfaceへのambient reservoir写像で
凍結される絶対電位差、Type Aを含む電子cutoff、ion drift energy、Zhaoが解く電子流入密度をそれぞれ監視します。
一方、光電子はprofileへ共通の電位移動ではなく、profile内のbarrier $B$とfield変形へ応答するため、
最後の2項は$T_{pe}$と$\lambda_{D,pe}$で測ります。同じ$\Delta\phi_I$を$T_{pe}$でも重複判定しません。

このCDFはbatchごとの生の実測標本なので、収束判定では総ray数だけでなくinterface到達数と
barrierより上のescape tailに残る有効標本数を確認する。source、field、barrierの変化が小さいまま
ambient-potential trust guardだけが統計的に失敗する場合もguardを緩めず`rays_per_batch`を増やし、
少なくとも2倍のray数で$\phi_I$、escape率、表面総電荷の一致を確認する。
ray数を増やしても同じguard超過が収束して残る場合は1 batchの物理更新幅が大きすぎるため、guardを緩めず
`batch_duration`と対応するbatch当たり放出量を小さくし、時間刻み収束を確認する。

受理後は各実測energy groupのescape / return weightをそのままsource elementと実際の再吸収elementへ使い、
ambient経路の共通`return_weight_scale`で再正規化しません。各rayは原則1回のoutward crossingと最大1回の
returnを持つ必要があり、return weightを持つrayの未吸収、escape rayのterminal不一致、再越境を
別rayへのweight付替えで隠さずfail closedとします。
標準出力の`BEACH implicit-zhao`行は、受理したbranch、$\phi_{\min}$、`n_e_inf_m3`、full-profile barrier、
resolved source scale、marginal group、empirical charge residual、recross / terminal mismatchを記録します。

この経路で実測するのはinterface source amplitudeと法線energy CDFです。Zhaoのfree/reflected ambient electron、
free/captured photoelectron、cold ionのdensity closureは解析式のままであり、接線速度分布、任意のtrapped VDF、
衝突、磁化、時間依存Vlasov--PoissonまたはPICを解きません。Type Aの単一virtual-cathode minimumと
Zhao式に含まれるreflected / captured populationだけが適用範囲です。

Zhao profileの物理scaleはphotoelectron温度$T_{pe}$と、$T_{pe}$および$n_{ref}$から導出した$\lambda_{D,pe}$です。
`outer_plasma.debye_length`と`thermal_voltage`はZhao root/profileへ入らず、現時点では`interface_eta_gap`、横方向phi/field比、
local-charge推定などsplit-interface適用性診断のreference inputとして残ります。tracked
`photo_raycast.emit_current_density_a_m2`は、$T_{pe}$をJへ換算して$v_{th,pe}=\sqrt{2T_{pe}/m_{pe}}$とし、有効平面のraw source
$s_{UV}|q_{pe}|n_{ref}\sin(\alpha)v_{th,pe}/(2\sqrt{\pi})$と、`implicit_mean`以外では
1%以内で一致しなければなりません。analytic raw currentは
tracked sourceの整合性検査とcurrent-density診断に使いますが、root、surface charge、ledgerへ別途加えず、tracked放出と
再吸収だけが後二者を更新します。`implicit_mean`では上記の実測$s_{\rm eff}$が解析populationのsource amplitudeを
置換し、surface currentとの1%一致を要求しません。過渡population scale $\eta$はphotoelectron密度、無限遠準中性、Sagdeev項を
scaleしますが、current診断のraw photoelectron emission-current項はscaleせず、full tracked sourceのまま使います。
`implicit_mean`以外では`ray_direction`やrough surfaceからinterfaceへ到達するVDFをZhao outer populationへ
自己無撞着に接続しません。`implicit_mean`でも実測接続はsource amplitudeと法線energy CDFに限られ、
`ray_direction`と`sheath_alpha_deg`はそれぞれ照射rayによる放出面samplingと解析density closureを独立に指定します。
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

`sim.field_normalization` で場計算内部の長さを正規化できます。`"si"` が既定で従来どおり、`"box"` は最大 box 幅、`"mesh"` は mesh bbox 最大幅、`"length"` は `sim.field_length_scale` を長さ基準 `L0` とします。direct/treecode/FMM の Coulomb kernel は座標と periodic cell を `L0` で割った無次元距離で評価し、電場で `k_coulomb/L0^2`、電位で `k_coulomb/L0` を掛けて SI に戻します。入力ファイルと出力 CSV は SI 単位のままです。

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
- outer transfer候補stepでは、Boris端点と整合する二次軌道でz-highと候補x/y面の交差時刻を再評価し、実際に最初のeventを選ぶ。周期・反射面が先なら作用後の残り時間を再積分し、z-highが先ならouterへ渡す。同時刻の周期・反射面は横方向作用を先に合成する。別open面も同時に含む場合はownerが一意でないためfail closedとする
- 二次補正は、候補終点がbox外側にあるため検出されたcrossingだけを対象とする。候補終点がbox内へ戻る途中の一時的越境は探索せず、mesh hitは従来のchord判定を維持する
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

- 通常実行では `batch_count` accepted batchを実行し、再開実行では処理済みaccepted batch数が
  `batch_count` に達するまで実行
- 各粒子は `max_step` で打ち切り

補足:

- `tol_rel` は `rel_change` の監視値として出力されますが、早期終了には使いません。
- 適応的な非零モード進行では、棄却trialはbatch数と`simulated_time_s`を進めません。

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

`summary.txt`は`batches`に加えて`simulated_time_s`、
`periodic2_max_nonzero_mode_potential_step_V`、`adaptive_nonzero_mode_rejected_trials`、
`adaptive_nonzero_mode_last_batch_duration_s`、
`adaptive_nonzero_mode_last_potential_step_V`、`adaptive_nonzero_mode_omp_threads`を記録します。
棄却trialは履歴CSVとcharge ledgerへ記録しません。
棄却数は共通ladderの総数であり、$k\ne0$上限超過と回復可能な`implicit_mean` $k=0$ failureの
内訳は標準出力のreject reasonで区別します。

MPI 実行時はRNGを`rng_state_rankNNNNN.txt`、Zhao過渡queueを`outer_event_queue_rankNNNNN.csv`としてrank別に保存します。
マクロ粒子残差は全rankで共有する状態なので、rootが単一の`macro_residuals.csv`を保存します。

### 8.2 再開（`output.resume = true`）

再開時は既定で `output.dir` から以下を読み込みます。`output.restart_from` を指定した場合は、そのディレクトリから checkpoint を読み込み、新しい出力は引き続き `output.dir` に書きます。

- 必須: `summary.txt`, `charges.csv`, `rng_state.txt`（MPI 時は `rng_state_rankNNNNN.txt`）
- 任意: `macro_residuals.csv`（MPI 時も単一の global ファイル）
- schema v2/v3/v4では電荷収支出力時の`charge_ledger.csv`
- queue有効時はserialの`outer_event_queue.csv`またはMPI全rankの`outer_event_queue_rankNNNNN.csv`

`sim.batch_count` は累積のaccepted batch到達数です。例えば checkpoint が `batches=100` のとき
`batch_count=150` で再開すると、追加で50 accepted batchだけ実行します。
`simulated_time_s`、累積棄却trial数、最後のaccepted trial幅・電位変化も`summary.txt`から復元します。
`batch_count` が checkpoint の処理済みバッチ数より小さい場合は停止します。MPI 実行時の再開では、前回と同一の `mpi_world_size` が必要です。
`output.resume=true` で必須 checkpoint が存在しない場合は新規実行へフォールバックせず停止します。`summary.txt` の統計値、`charges.csv` の電荷、`macro_residuals.csv` の残差は resume 時に有限性と基本範囲を検証します。
新規出力の`summary.txt`は`checkpoint_schema_version=4`とmodel / ordered mesh / ordered species fingerprintを持ちます。schema v4はschema v3のreadyなouter held stateに加え、Zhao過渡closureのpopulation fraction、column target/residual、queue inventoryを復元します。queue本体はserialの`outer_event_queue.csv`またはMPI rank別ファイルにactive event、terminal outcome、due時刻、`next_event_id`を保存します。queue有効時はschema、rank、world size、完了batchの不一致や欠落ファイルをfail closedで拒否します。schema v3は`outer_plasma_profile.csv`に`z, phi, E, rho`を保存し、outer solverのstatus、反復数、residual、積分電荷、species別電流とともにheld stateを復元します。schema v2は3 fingerprintを照合した上でread-only migrationとして受理しますが、旧3列outer profileは初期値にだけ使い、次のouter refreshでroot solveを強制します。廃止された重心source modelを使ったcheckpointは再開できません。
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

### 11.1 電位・電場計算

- simulator が出力した `mesh_potential.csv` / `potential_history.csv` を、要素中心電位の正本として読み込む
- 任意点の電場・電位と object interaction は、三角形 geometry を渡す native `triangle_p0` field kernel で評価する
- 削除済みの重心 source modelをPython側で再構成して、simulatorと同じ電場・電位として扱わない

### 11.2 Coulomb 力と periodic2 対応

- `calc_coulomb`: native P0 panel kernelとtarget面積積分によるメッシュグループ間の Coulomb 力/トルク計算
- `periodic2` パラメータ: Fortran 側の `sim.field_bc_mode="periodic2"` に対応し、finite imageをnative kernelで評価する。`cached_kneq0`は物理的zero modeを合成する専用API以外では拒否する
- `periodic2=None` では出力ディレクトリ近傍の `beach.toml` から `field_bc_mode` を自動判定する。設定が見つからないtotal-field再計算はfree spaceへfallbackしない

### 11.3 電気力線追跡

- `trace_field_lines`: シード点から RK4 積分で電気力線を追跡。`direction` で順方向 / 逆方向 / 両方向を選択可能
- `plot_field_lines_3d`: 力線を 3D 描画し、三角形メッシュを電荷密度で着色してオーバーレイ

### 11.4 Python 側の制限

- Python 側の電場・電位計算も、三角形 geometry を渡す native `triangle_p0` field kernel を使う
- `cached_kneq0` は非零モードだけを返す低水準 operator であり、total field として扱うには物理的 zero mode を同じ規約で合成する
