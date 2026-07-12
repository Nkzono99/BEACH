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
- 静電場（既定は要素重心の点電荷近似 + softening。`field.element_kernel="triangle_p0"` では要素総電荷を三角形上の一定面密度として積分する厳密 free-space direct 核を使用）
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

1. 設定に従ってそのバッチの粒子群を生成
2. 現在の要素電荷に基づいて場ソルバをリフレッシュ (`field_solver%refresh(mesh)`)
3. 各粒子を `max_step` まで前進（OpenMP スレッド並列）
4. 各ステップで
   - 同一時刻の状態 `x0, v0` から予測中点 `x_mid = x0 + 0.5*v0*dt` を計算
   - 場評価点だけをsolverのprimitive target boxへ写像する（周期軸はwrap、非周期軸はbox面へclamp）。軌道候補座標は写像しない
   - 境界要素電場 `E(x_mid)` を1回評価し、一様外部電場 `sim.e0` を1回加える
   - 中点場を使う Boris push と台形位置更新で `x1, v1` を計算
   - `x1` がbox内部なら `x0 -> x1` のmesh collisionを1回探索
   - box faceへ到達する場合は、full chordのmesh hit parameterと最初のface event parameterを比較して最早順序を決める
   - periodic2のfull-chord queryがbox外区間でrange limitに達した場合は、最初のface eventまでに制限して再照会する
   - reflect/periodic後は残り時間を同じBoris規約で再積分し、そのchordのmesh hitを探索（1 outer stepにつき最大8 box event）
   - 衝突時: 粒子を消滅し `q_particle * w_particle` をスレッド別バッファ `dq_thread(elem_idx, tid)` へ加算
   - 残り時間中に9回目のbox eventへ到達し、それ以前にmesh hitがなければ、状態をcommitせず `dt` 縮小を要求する明示的failureとする
5. バッチ終了時に要素電荷差分をコミット: 全スレッドの `dq_thread` を合算し、`photo_emission_dq` を加算した後、MPI allreduce を行い `mesh%q_elem` に反映
6. `rel_change = ||dq|| / max(||q||, q_floor)` を更新
7. 統計と履歴を更新

`photo_raycast` で `deposit_opposite_charge_on_emit=true` の場合は、放出元要素に `-q_particle * w_hit` も加算します。
`photo_escape_model="boltzmann_cutoff"` では、PE escape 係数を掛けた実効重み `w_eff` を粒子重みと放出元要素の逆符号電荷の両方に使います。

## 5. 物理モデル

### 5.1 電場

Fortran 本体の電場計算は次式です（要素重心点電荷近似）:

- `E(r) = k * Σ_j q_j * (r - c_j) / (|r - c_j|^2 + softening^2)^(3/2)`

ここで `c_j` は要素 `j` の重心です。
`field_solver="treecode"` のときはこの核を遠方で monopole 近似し、近傍は direct 和を使います。  
`field_solver="fmm"` のときは simulator 非依存の Coulomb FMM コアを使い、source octree、optional target tree、Cartesian tensor による multipole/local 展開、近傍 direct 和で電場を評価します。現行 adapter の内部既定次数は 4 です。詳しくは `docs/Algorithms.md` の FMM コア詳細を参照してください。

`field.element_kernel="triangle_p0"` は `field_solver="direct" | "fmm" | "auto"`、`softening=0`、`surface_model="insulator"` に限定します。direct は厳密 free-space panel 和、FMM は全頂点を含む topology、近傍の厳密 panel 核、三角形 monomial の厳密 P2M を使います。auto は `tree_min_nelem` 未満で direct、以上で FMM を選びます。treecode と point-source `m2l_root_oracle` は拒否します。各 OBJ または template に `surface_side="normal_plus" | "normal_minus" | "outward_closed"` が必要です。`q_elem` は要素総電荷 [C]、面密度は `q_elem/area` です。面上電位は連続、法線電場は `sigma/epsilon0` だけ跳び、重心電位と principal-value 電場を自己項として用います。非対応 solver へ点電荷 fallback はしません。

`sim.field_bc_mode="periodic2"` かつ `field_solver="fmm"` では、`bc_low/high` が `periodic` の2軸を周期軸として扱います（第三軸は開放）。  
近傍画像和は `sim.field_periodic_image_layers = N` に対して各周期軸 `[-N, N]` を評価します。`periodic2` の遠方補正の既定は `field_periodic_far_correction="none"`（`sim` table）です。`auto` は互換用に受理され、`none` に正規化されます。`none` は explicit image shell だけを評価する有限画像近似であり、完全な周期遠方場を与えるものではありません。`m2l_root_oracle` は `k=0` と charged-wall closure を含む明示 opt-in の高コスト診断 backend で、production physics には使いません。

`field_periodic_far_correction="cached_kneq0"` は production 用の無限 periodic2 非零モード backend です。runtime が加算する有限画像 kernel を `K_shell(N)` とすると、cache は滑らかな full-periodic Ewald residual を root-local operator として保持します。charge refresh 時に source 高さ分布から対称 `k=0` state を一度構築し、各 eval で O(log n) で差し引くため、runtime total は代数的に `K_periodic,k!=0` になります。Ewald all-source 和は cache miss 時の operator generation にだけ使い、particle eval hot path では使いません。物理的な `k=0` は `exclude_k0` / `e_bottom_zero` の zero-mode provider が snapshot 内で一度だけ加えます。したがって non-neutral cell も暗黙の charged walls ではなく、この明示的な zero-mode boundary condition で閉じます。cache fingerprint は周期長、FMM order、画像/Ewald 層、source/target topology、softening、generator version、tolerance、real kind、build version を含みます。MPI では rank 0 だけが lock 下で検証・生成・atomic publish し、operator を全 rank へ broadcast します。
`tree_theta`/`tree_leaf_max` を未指定の場合は、`periodic2` でも通常の自動推定値を使います。現行実装の推定値は `nelem < 1500` で `theta=0.40`, `leaf_max=12`、`1500 <= nelem < 10000` で `0.50` / `16`、`10000 <= nelem < 50000` で `0.58` / `20`、`50000 <= nelem` で `0.65` / `24` です。

`periodic2.nonzero_mode_backend="panel_spectral_reference"` は、P0 panelのFourier `k!=0` 成分、triangle-heightの厳密`k=0`成分、線形Debye outer plasmaを合成する小規模 correctness referenceです。この経路だけは`field_solver="direct"`を用い、`zero_mode_policy="exclude_k0"`、`lower_boundary_model="e_bottom_zero"`、x/y periodic・z open、`e0=0`を必須とします。有限image shellや`charged_walls`とは混用しません。interface面の`k!=0`減衰、gap、局所平均plasma電荷推定、線形性を実測し、設定閾値を超えた場合は`not_applicable`として停止します。外部状態は`outer_update_stride`とともにcheckpointされ、restart後も更新位相を保存します。

`coupling.particle_transfer_mode="electrostatic_1d_instant_return"`では、z-high面を唯一のouter particle interfaceとします。無限遠reservoirの法線VDFはLiouville/flux保存と法線エネルギー保存でinterfaceへ写像します。外向き粒子は同じ線形Debye profile上でinfinityへ到達可能ならescape、turning pointを持つ場合は法線速度を反転してlocalへ返します。return位置のx/yには`v_t*tau_outer`を加えて周期wrapし、残り`dt`を既存stepperで再積分します。この写像はouter flightをglobal simulation timeへ加算せず、`tau_outer/field_evolution_timescale`が上限を超える場合は停止します。persistent outer queueは未対応で、`outer_queue_enabled=true`を拒否します。`b0=0`のみを許可し、磁化outer orbitを近似しません。

`outer_plasma.model="unified_linear_response"` と
`coupling.particle_transfer_mode="electrostatic_3d_explicit_orbit"` の組合せでは、z-high ownership面を
出た粒子を、unified zero modeとscreened nonzero tailを合成した同じ3D静電場中で固定刻み
velocity-Verletにより追跡します。ownership面へ戻る粒子はx/yを周期wrapしてlocal stepperへ返し、
unified grid上端のfar planeを外向きに通過する粒子はinfinity escapeとします。全エネルギー相対誤差、
outer flight time、frozen-field ratioを検査し、step上限到達をdiscardしません。persistent queueが未実装の
ため該当軌道は停止します。外部磁場を無視する条件は未決なので`b0=0`のみを許可します。

`outer_plasma.photoelectron_closure="individual_return"`では、`photo_raycast`粒子がz-high interfaceを外向き通過した時点で、法線運動エネルギーbinごとのsigned charge、全運動エネルギー、接線運動量、個数をMPI-globalに集計します。個別粒子は上記instant-return写像だけで帰還させ、統計量から別粒子を再注入しません。各batchの統計は`previous_batch`として次batchへ渡し、累積統計とともに`photoelectron_histogram.csv`へcheckpointします。charge-conserving modeでは全`photo_raycast` speciesに`deposit_opposite_charge_on_emit=true`を要求し、legacy `photo_escape_model`との併用を拒否します。`statistical_return`は帰還位置・遅延・deposit則が未仕様なので使用不可です。放出signed chargeと`photoelectron_ambient_charge_scale`の比が`max_photoelectron_charge_ratio`を超えると、ambient-only線形closureを適用外として停止し、silent fallbackしません。

`outer_plasma.model="kinetic_1d"`は、z-highの負・正`reservoir_face` speciesを無限遠の
electron half-Maxwellian / cold drifting ion VDFとして用い、伸長1D格子上のPoisson方程式を
interface Neumann条件と遠方Robin条件で解きます。初版は単調・無衝突・非磁化分枝に限定し、
ionにはkinetic Bohm入口条件を課します。`photoelectron_closure="kinetic_mean"`は負電荷
`photo_raycast` speciesの放出fluxからoutgoing/returning平均密度を構成します。解状態は
`converged`、`not_applicable`、`no_physical_solution`、`numerical_failure`を区別し、線形モデルへ
silent fallbackしません。profileは`outer_plasma_profile.csv`へ保存し、restart時のNewton初期値に使います。

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
- 複数候補がある場合は最小 `t`（最初の衝突）を採用
- `sim.field_bc_mode="periodic2"` では、mesh は primitive cell の base element のみ保持しつつ、衝突判定だけ periodic image を implicit に探索する
- periodic2 collision 用 mesh は runtime で canonical unwrapped 形へ平行移動して collision grid を再構築する。raw 頂点は periodic 軸で box 外を含んでよい
- periodic2 hit の `elem_idx` は常に base element index を返し、吸収・堆積先も primitive cell 側へ集約する

### 5.4 ボックス境界

- `open`: 粒子を消滅（`escaped_boundary`）
- `reflect`: 法線成分反転
- `periodic`: 反対側へラップ
- additive な `find_first_boundary_event` は、box 内の始点から候補終点までの最初の交差 fraction と、corner/edge で同時に交差する全 face を bit mask で返す
- additive な `apply_escape_reflect_periodic_event` は同時 face を軸順序に依存せず一括適用し、reflect/periodic 後の位置を境界から1 ULP内側へ置く。非有限値、不正な box/face、event/config 不一致は state を変更せず明示 status を返す
- production particle loop はcandidate生成とmesh queryを先行し、box crossing時だけevent resolverへ進む。候補終点がstrictなbox内部なら追加event geometryを行わず、場評価1回・collision query 1回のfast pathとなる
- reflect/periodic crossingだけ残り時間を最大8回再積分する。各eventはmeshとの最早順序を保って処理し、9回目のbox eventまでにmesh hitがなければ `particle_step_multiple_box_events` でfail closedとする。上限なしのevent loopやadaptive substepは行わない
- `open_boundary_model="potential_barrier"` は既存の単一面scalar energy式をevent位置・補間速度で使うlegacy/experimental扱いとする。複数open faceへの一般化は行わずfail closedとし、物理モデルはshared potential snapshotとともに後段で再設計する
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
- `reservoir_potential_model="infinity_barrier"` 時は注入面平均電位を使って法線速度下限を補正
- `sheath_injection_model` が有効な場合、最初の負電荷 `reservoir_face` species は共有シース解に基づく `n_swe_inf` と `vmin_normal` で上書きされる
- シース 1D 座標の基準平面は共有 `inject_face` の法線方向で定義し、`sim.sheath_reference_coordinate` があればその座標を、未指定なら対応 box face の座標を使う
- `sim.sheath_reference_coordinate` を明示した Zhao モデルでは、基準平面から reservoir 境界までの局所 `phi(z)` を使って electron reservoir の有効密度・cutoff と ion reservoir の局所密度・法線ドリフトを更新し、シースによる VDF 変形を近似する

### 6.3 `photo_raycast`

- 指定面から `rays_per_batch` 本を照射
- 各レイは最初の命中要素からのみ放出（命中しなければ放出なし）
- 1ヒット重み:
  - `w_hit = J_perp * A_perp * batch_duration / (|q| * rays_per_batch)`
  - MPI 実行時は `rays_per_batch` の代わりに `global_rays_per_batch`（全 rank 合計）を使用
- `photo_escape_model="boltzmann_cutoff"` の場合:
  - 放出元要素の自己寄与を除いた中心電位から `barrier = max(phi_emit - phi_infty, 0)` を計算
  - `escape_factor = exp(-|q_particle| * barrier / (k_B * T_PE))`
  - `w_eff = w_hit * escape_factor` とし、PE粒子重みと `deposit_opposite_charge_on_emit` の放出元電荷 bookkeeping に同じ `w_eff` を使う
  - これは戻りPEを即時中和として扱う reduced closure であり、個別PEの再吸収面は追跡しない
- `sim.field_bc_mode="periodic2"` で periodic image に命中した場合も、放出位置は primary cell に wrap した hit 座標を使う
- `sheath_injection_model` が Zhao 系のとき、最初の負電荷 `photo_raycast` species の `emit_current_density_a_m2` は Zhao の自由光電子電流へ上書きされ、法線速度 cutoff も分枝に応じて適用される

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
- `charge_ledger.csv`（species 別の signed charge flux と粒子数）
- `performance_profile.csv`（`BEACH_PROFILE=1` 環境変数設定時）

`mesh_triangles.csv` は要素ごとの `mesh_id` を含み、`mesh_sources.csv` で `mesh_id` と元メッシュ設定を対応付けます。

MPI 実行時は RNG のみ `rng_state_rankNNNNN.txt` として rank 別に保存します。マクロ粒子残差は
全 rank で共有する状態なので、root が単一の `macro_residuals.csv` を保存します。

### 8.2 再開（`output.resume = true`）

再開時は既定で `output.dir` から以下を読み込みます。`output.restart_from` を指定した場合は、そのディレクトリから checkpoint を読み込み、新しい出力は引き続き `output.dir` に書きます。

- 必須: `summary.txt`, `charges.csv`, `rng_state.txt`（MPI 時は `rng_state_rankNNNNN.txt`）
- 任意: `macro_residuals.csv`（MPI 時も単一の global ファイル）
- schema v2/v3 では台帳出力時の `charge_ledger.csv`

`sim.batch_count` は累積の到達バッチ数です。例えば checkpoint が `batches=100` のとき `batch_count=150` で再開すると、追加で50バッチだけ実行します。`batch_count` が checkpoint の処理済みバッチ数より小さい場合は停止します。MPI 実行時の再開では、前回と同一の `mpi_world_size` が必要です。
`output.resume=true` で必須 checkpoint が存在しない場合は新規実行へフォールバックせず停止します。`summary.txt` の統計値、`charges.csv` の電荷、`macro_residuals.csv` の残差は resume 時に有限性と基本範囲を検証します。
新規出力の `summary.txt` は `checkpoint_schema_version=3` と model / ordered mesh / ordered species fingerprint を持ちます。schema v3 は `outer_plasma_profile.csv` に `z, phi, E, rho` を保存し、outer solver の status、反復数、residual、積分電荷、species 別電流とともに held state を完全復元します。schema v2 は3 fingerprint を照合した上で read-only migration として受理しますが、旧3列 outer profile は初期値にだけ使い、次の outer refresh で root solve を強制します。それより古い schema は Phase 0 で実装済みの legacy point-source model に限って読み込めます。
`charge_ledger_residual_C` は surface / flight / unresolved stock の差分と外部 flux から作る transactional 保存残差です。species 間相殺を避けた `discarded_unresolved_abs_C` は別診断であり、残差が 0 でも max-step discard が物理的に許容されることを意味しません。
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
- pure Python の電位再構成は `m2l_root_oracle` / `cached_kneq0` operator を再現せず、explicit image shell のみ。native `FieldKernel` は `cached_kneq0` を利用できる
