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
- 静電場（要素重心の点電荷近似 + softening、`sim.field_solver` による direct/treecode/fmm 切替。`auto` は要素数で方式とtreecode設定を自動選択）
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
- Hybrid 近傍場（`use_hybrid`, `r_switch_factor`, `n_sub`, `softening_factor` は予約キー）

## 3. データモデル

### 3.1 メッシュ（`mesh_type`）

- 要素数 `nelem`
- 三角形頂点 `v0`, `v1`, `v2`
- 要素重心 `centers`（`center_x/y/z`）
- 要素法線 `normals`
- 要素電荷 `q_elem`
- 衝突高速化用 AABB / grid キャッシュ

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
- `field_time_s`, `push_time_s`, `collision_time_s`

## 4. 1バッチの計算手順

1. 設定に従ってそのバッチの粒子群を生成
2. 各粒子を `max_step` まで前進
3. 各ステップで
- 電場 `E(x)` を評価
- Boris push で `x1, v1` を計算
- `x0 -> x1` の最初の衝突要素を探索
- 衝突時: 粒子を消滅し `q_particle * w_particle` を該当要素へ加算
- 非衝突時: ボックス境界を適用（open/reflect/periodic）
4. バッチ終了時に要素電荷差分をコミット
5. `rel_change = ||dq|| / max(||q||, q_floor)` を更新
6. 統計と履歴を更新

`photo_raycast` で `deposit_opposite_charge_on_emit=true` の場合は、放出元要素に `-q_particle * w_hit` も加算します。

## 5. 物理モデル

### 5.1 電場

Fortran 本体の電場計算は次式です（要素重心点電荷近似）:

- `E(r) = k * Σ_j q_j * (r - c_j) / (|r - c_j|^2 + softening^2)^(3/2)`

ここで `c_j` は要素 `j` の重心です。
`field_solver="treecode"` のときはこの核を遠方で monopole 近似し、近傍は direct 和を使います。  
`field_solver="fmm"` のときは octree の葉ごとに遠方ノード寄与を局所展開（1次）へ事前集約し、近傍は direct 和を使います。

### 5.2 粒子前進

- Boris 法（`E`, `B`）
- `B` は `sim.b0` の一様場

### 5.3 衝突判定

- 線分 `x0 -> x1` と三角形群の交差
- 複数候補がある場合は最小 `t`（最初の衝突）を採用

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

### 6.3 `photo_raycast`

- 指定面から `rays_per_batch` 本を照射
- 各レイは最初の命中要素からのみ放出（命中しなければ放出なし）
- 1ヒット重み:
  - `w_hit = J_perp * A_perp * batch_duration / (|q| * rays_per_batch)`

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
- `rng_state.txt`
- `macro_residuals.csv`

`mesh_triangles.csv` は要素ごとの `mesh_id` を含み、`mesh_sources.csv` で `mesh_id` と元メッシュ設定を対応付けます。

### 8.2 再開（`output.resume = true`）

再開時は同一 `output.dir` から以下を読み込みます。

- 必須: `summary.txt`, `charges.csv`, `rng_state.txt`
- 任意: `macro_residuals.csv`

`sim.batch_count` は「今回追加するバッチ数」です。

## 9. 設計方針

- v0.x は insulator accumulation を正規仕様とする
- 拡張点は維持しつつ、現行利用者向けには実装済み挙動を優先して文書化する
- 設定追加・削除時は `docs/fortran_parameter_file.md` を同時更新する

## 10. 実行運用（推奨）

通常利用では、`pip install git+...` で導入した `beach` コマンドの利用を推奨します。

- 推奨導線: `python -m pip install "git+https://github.com/Nkzono99/BEACH.git"` で導入し、`beach [config.toml]` で実行
- 開発時の導線: `python -m pip install -e .` + `make` を利用し、必要に応じて `fpm run -- ...` で直接実行
