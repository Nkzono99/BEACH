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
- ボックス境界条件（open / reflect / periodic）
- 粒子注入モード
  - `volume_seed`
  - `reservoir_face`
  - `photo_raycast`
- チェックポイント再開（`resume`）

### 2.2 現行で未実装・予約

- 導体再分配モデル
- 表面導電/拡散モデル
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
   - 電場 `E(x)` を評価
   - Boris push で `x1, v1` を計算
   - `x0 -> x1` の最初の衝突要素を探索
   - 衝突時: 粒子を消滅し `q_particle * w_particle` をスレッド別バッファ `dq_thread(elem_idx, tid)` へ加算
   - 非衝突時: ボックス境界を適用（open/reflect/periodic）
5. バッチ終了時に要素電荷差分をコミット: 全スレッドの `dq_thread` を合算し、`photo_emission_dq` を加算した後、MPI allreduce を行い `mesh%q_elem` に反映
6. `rel_change = ||dq|| / max(||q||, q_floor)` を更新
7. 統計と履歴を更新

`photo_raycast` で `deposit_opposite_charge_on_emit=true` の場合は、放出元要素に `-q_particle * w_hit` も加算します。

## 5. 物理モデル

### 5.1 電場

Fortran 本体の電場計算は次式です（要素重心点電荷近似）:

- `E(r) = k * Σ_j q_j * (r - c_j) / (|r - c_j|^2 + softening^2)^(3/2)`

ここで `c_j` は要素 `j` の重心です。
`field_solver="treecode"` のときはこの核を遠方で monopole 近似し、近傍は direct 和を使います。  
`field_solver="fmm"` のときは simulator 非依存の Coulomb FMM コアを使い、source octree、optional target tree、Cartesian tensor による multipole/local 展開、近傍 direct 和で電場を評価します。現行 adapter の内部既定次数は 4 です。詳しくは `docs/fortran_fmm_core.md` を参照してください。

`sim.field_bc_mode="periodic2"` かつ `field_solver="fmm"` では、`bc_low/high` が `periodic` の2軸を周期軸として扱います（第三軸は開放）。  
近傍画像和は `sim.field_periodic_image_layers = N` に対して各周期軸 `[-N, N]` を評価します。`periodic2` の遠方補正は `auto` を既定とし、`m2l_root_oracle` をサポートします。`sim.field_periodic_far_correction="auto"` は `field_solver="fmm"` かつ `field_bc_mode="periodic2"` のとき内部的に `m2l_root_oracle` に正規化され、`none` は遠方補正を無効化する独立の値であり、`auto` の alias ではありません。`m2l_root_oracle` は build 時だけ exact periodic Ewald residual を oracle として使い、proxy/check 点から root local 演算子へ fit します。runtime では `local(:,root) += T_root * multipole(:,root)` の形で root local へ注入され、tree 外 fallback では同じ exact Ewald correction を直接足します。非中性ケースでは、slab 外評価に対して `charged_walls` に対応する total-charge 補正を追加します。2 枚の補償壁の場は slab 内では相殺されるため、粒子前進に使う in-box field は従来どおり periodic pair field と一致します。`field_periodic_ewald_alpha` は `m2l_root_oracle` の Ewald 分解パラメータ、`field_periodic_ewald_layers` は real/reciprocal の打切り深さとして使います。
`tree_theta`/`tree_leaf_max` を未指定の場合は、`periodic2` でも通常の自動推定値を使います。現行実装の推定値は `nelem < 1500` で `theta=0.40`, `leaf_max=12`、`1500 <= nelem < 10000` で `0.50` / `16`、`10000 <= nelem < 50000` で `0.58` / `20`、`50000 <= nelem` で `0.65` / `24` です。

### 5.2 粒子前進

- Boris 法（`E`, `B`）
- `B` は `sim.b0` の一様場

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

## 6. 注入モード

### 6.1 `volume_seed`

- 各 species で `npcls_per_step` 個を毎バッチ生成
- 総生成数は `batch_count * Σ npcls_per_step`

### 6.2 `reservoir_face`

- 上流ドリフト Maxwell 分布から流入フラックス `gamma_in` を計算
- `n_macro_expected = gamma_in * A * batch_duration / w_particle`
- 残差繰越つきで `floor` して今バッチのマクロ粒子数を決定
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
- `sim.field_bc_mode="periodic2"` で periodic image に命中した場合も、放出位置は primary cell に wrap した hit 座標を使う
- `sheath_injection_model` が Zhao 系のとき、最初の負電荷 `photo_raycast` species の `emit_current_density_a_m2` は Zhao の自由光電子電流へ上書きされ、法線速度 cutoff も分枝に応じて適用される

## 7. 実行制御と停止条件

現行実装の停止条件は次の通りです。

- ループは常に `batch_count` バッチ実行
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

MPI 実行時はランク別ファイルが生成されます: `rng_state_rankNNNNN.txt`, `macro_residuals_rankNNNNN.csv`。

### 8.2 再開（`output.resume = true`）

再開時は同一 `output.dir` から以下を読み込みます。

- 必須: `summary.txt`, `charges.csv`, `rng_state.txt`（MPI 時は `rng_state_rankNNNNN.txt`）
- 任意: `macro_residuals.csv`（MPI 時は `macro_residuals_rankNNNNN.csv`）

`sim.batch_count` は「今回追加するバッチ数」です。MPI 実行時の再開では、前回と同一の `mpi_world_size` が必要です。

## 9. 設計方針

- v0.x は insulator accumulation を正規仕様とする
- 拡張点は維持しつつ、現行利用者向けには実装済み挙動を優先して文書化する
- 設定追加・削除時は `docs/fortran_parameter_file.md` を同時更新する

## 10. 実行運用（推奨）

通常利用では、`pip install beach-bem` で導入した `beach` コマンドの利用を推奨します。

- 推奨導線: `pip install beach-bem` で導入し、`beach [config.toml]` で実行
- 開発版導線: `python -m pip install "git+https://github.com/Nkzono99/BEACH.git"` で導入
- 開発時の導線: `python -m pip install -e .` + `make` を利用し、必要に応じて `fpm run -- ...` で直接実行

## 11. Python 後処理レイヤ

Python パッケージ `beach` は Fortran 出力を読み込み、後処理・可視化を行います。
詳細な API リファレンスは `docs/python_postprocess_api.md` を参照してください。

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
