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
- 静電場（要素重心の点電荷近似 + softening、`sim.field_solver` による direct/treecode/fmm 切替。`treecode`/`fmm`/`auto` では `tree_theta`/`tree_leaf_max` を要素数から自動推定し、明示指定があれば優先。場境界は `sim.field_bc_mode` で指定）
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
   - 境界要素電場 `E(x_mid)` を1回評価し、一様外部電場 `sim.e0` を1回加える
   - 中点場を使う Boris push と台形位置更新で `x1, v1` を計算
   - `x1` がbox内部なら `x0 -> x1` のmesh collisionを1回探索
   - box faceへ到達する場合は、full chordのmesh hit parameterと最初のface event parameterを比較して最早順序を決める
   - periodic2のfull-chord queryがbox外区間でrange limitに達した場合は、最初のface eventまでに制限して再照会する
   - reflect/periodic後は残り時間を同じBoris規約で一度だけ再積分し、そのchordのmesh hitを探索
   - 衝突時: 粒子を消滅し `q_particle * w_particle` をスレッド別バッファ `dq_thread(elem_idx, tid)` へ加算
   - 残り時間中に2回目のbox eventへ到達し、それ以前にmesh hitがなければ、状態をcommitせず `dt` 縮小を要求する明示的failureとする
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

`sim.field_bc_mode="periodic2"` かつ `field_solver="fmm"` では、`bc_low/high` が `periodic` の2軸を周期軸として扱います（第三軸は開放）。  
近傍画像和は `sim.field_periodic_image_layers = N` に対して各周期軸 `[-N, N]` を評価します。`periodic2` の遠方補正の既定は `field_periodic_far_correction="none"`（`sim` table）です。`auto` は互換用に受理され、adapter と standalone FMM core の両方で `none` に正規化されます。`none` は explicit image shell だけを評価する有限画像近似であり、完全な周期遠方場を与えるものではありません。`m2l_root_oracle` は明示 opt-in の高コスト診断 backend です。build 時だけ exact periodic Ewald residual を oracle として使い、proxy/check 点から root local 演算子へ fit しますが、production exact periodic physics の保証ではありません。runtime では `local(:,root) += T_root * multipole(:,root)` の形で root local へ注入され、tree 外 fallback では同じ Ewald correction を直接足します。非中性ケースでは、slab 外評価に対して `charged_walls` に対応する total-charge 補正を追加します。2 枚の補償壁の場は slab 内では相殺されるため、粒子前進に使う in-box field は従来どおり periodic pair field と一致します。`field_periodic_ewald_alpha` は `m2l_root_oracle` の Ewald 分解パラメータ、`field_periodic_ewald_layers` は real/reciprocal の打切り深さとして使います。将来の解析的 M2L 遠方補正を既定化する場合は Stage 5 の versioned migration として導入し、現在の `auto` から暗黙に有効化しません。
`tree_theta`/`tree_leaf_max` を未指定の場合は、`periodic2` でも通常の自動推定値を使います。現行実装の推定値は `nelem < 1500` で `theta=0.40`, `leaf_max=12`、`1500 <= nelem < 10000` で `0.50` / `16`、`10000 <= nelem < 50000` で `0.58` / `20`、`50000 <= nelem` で `0.65` / `24` です。

`sim.field_normalization` で場計算内部の長さを正規化できます。`"si"` が既定で従来どおり、`"box"` は最大 box 幅、`"mesh"` は mesh bbox 最大幅、`"length"` は `sim.field_length_scale` を長さ基準 `L0` とします。direct/treecode/FMM の Coulomb kernel は座標・softening・periodic cell を `L0` で割った無次元距離で評価し、電場で `k_coulomb/L0^2`、電位で `k_coulomb/L0` を掛けて SI に戻します。入力ファイルと出力 CSV は SI 単位のままです。

### 5.2 粒子前進

- Boris 法（`E`, `B`）
- `B` は `sim.b0` の一様場
- public な粒子入力 `x, v` と出力 `x_new, v_new` は同一時刻の状態であり、half-step staggered 状態ではない
- production の空間電場は予測中点 `x_mid = x + 0.5*v*dt` で1回評価し、`sim.e0` はその評価結果へ1回だけ加える
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
- reflect/periodic crossingだけ残り時間を一度再積分する。2回目のbox eventまでにmesh hitがなければ `particle_step_multiple_box_events` でfail closedとし、任意回数のevent loopやadaptive substepは行わない
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
- `performance_profile.csv`（`BEACH_PROFILE=1` 環境変数設定時）

`mesh_triangles.csv` は要素ごとの `mesh_id` を含み、`mesh_sources.csv` で `mesh_id` と元メッシュ設定を対応付けます。

MPI 実行時は RNG のみ `rng_state_rankNNNNN.txt` として rank 別に保存します。マクロ粒子残差は
全 rank で共有する状態なので、root が単一の `macro_residuals.csv` を保存します。

### 8.2 再開（`output.resume = true`）

再開時は既定で `output.dir` から以下を読み込みます。`output.restart_from` を指定した場合は、そのディレクトリから checkpoint を読み込み、新しい出力は引き続き `output.dir` に書きます。

- 必須: `summary.txt`, `charges.csv`, `rng_state.txt`（MPI 時は `rng_state_rankNNNNN.txt`）
- 任意: `macro_residuals.csv`（MPI 時も単一の global ファイル）

`sim.batch_count` は累積の到達バッチ数です。例えば checkpoint が `batches=100` のとき `batch_count=150` で再開すると、追加で50バッチだけ実行します。`batch_count` が checkpoint の処理済みバッチ数より小さい場合は停止します。MPI 実行時の再開では、前回と同一の `mpi_world_size` が必要です。
`output.resume=true` で必須 checkpoint が存在しない場合は新規実行へフォールバックせず停止します。`summary.txt` の統計値、`charges.csv` の電荷、`macro_residuals.csv` の残差は resume 時に有限性と基本範囲を検証します。
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
- periodic2 の Ewald 遠方補正（`m2l_root_oracle`）は Python 側では再現しない。explicit image shell のみ
