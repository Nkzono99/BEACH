# BEMテストパーティクルシミュレータ 仕様書（v0.1 / Pythonプロトタイプ）

## 1. 目的
BEM（Boundary Element Method）で離散化された境界（三角形要素）上の電荷分布から静電場を計算し、その電場（＋任意の磁場）中でテスト粒子を追跡する。粒子が境界に衝突した場合、衝突要素へ電荷を付与し、一定単位（npcls_per_step）ごとに電荷分布を確定（commit）して電場を更新する反復計算を行う。

本バージョン（v0.1）は **絶縁体モード（蓄積のみ）** を実装対象とし、将来的に導体・有限導電率へ拡張可能な設計を採用する。ただし、導体BEM連立解法等の実装は v0.1 のスコープ外とする。

> 注：Offset（境界外側δオフセット参照）は **採用しない**。

---

## 2. スコープ
### 2.1 実装対象（Must）
- 境界：三角形要素メッシュ（`bem_list`）
- 電場：境界要素電荷からのクーロン場（静電）
- 粒子運動：Boris法（Eのみ、E+Bどちらも対応）
- 衝突判定：粒子線分ステップ vs 三角形要素（最初に当たった要素）
- 相互作用：吸収（衝突で粒子消滅、表面へ電荷付与）
- 電荷更新：`npcls_per_step` 単位でバッチコミット、コミット後に電場キャッシュ更新
- 停止条件：相対変化率 < tol もしくは step >= max_step

### 2.2 非対象（Won’t / v0.1）
- 導体（再分配）BEM連立解
- 有限導電率（表面拡散）モデル
- 二次電子放出（SEE）、鏡面反射、拡散反射（ただし追加可能な拡張点は用意）
- 高速化（FMM / BVH本格導入等）は最適化フェーズで実施

---

## 3. 用語・前提
- **要素電荷** `q_j [C]`：三角形要素 j が保持する総電荷
- **面電荷密度** `σ_j [C/m^2]`：必要なら `σ_j = q_j / A_j`
- **superparticle weight** `w`：粒子が代表する実粒子数（電荷付与量に掛ける）
- 単位系：SI（m, s, kg, C, V）

---

## 4. データモデル
### 4.1 境界要素（BEMElement）
- `v0, v1, v2: ndarray(3,)`：三角形頂点座標
- `q: float`：要素電荷 [C]
- 派生量（前計算推奨）
  - `centroid: ndarray(3,)`
  - `area: float`
  - `normal: ndarray(3,)`
  - `h_elem: float = sqrt(area)`（代表長さ）

### 4.2 メッシュ（BEMMesh）
- `elements: list[BEMElement]`
- 前計算キャッシュ
  - `centers: (N,3)`、`areas: (N,)`、`normals: (N,3)`
  - `bb_min/max: (N,3)`（三角形AABB、衝突のbroad-phase用）
- API
  - `charges() -> (N,)`
  - `add_charges(dq: (N,))`

### 4.3 粒子（Particle）
- `x: ndarray(3,)` 位置 [m]
- `v: ndarray(3,)` 速度 [m/s]
- `q: float` 電荷 [C]
- `m: float` 質量 [kg]
- `w: float` 重み（既定 1.0）
- `alive: bool`

### 4.4 衝突情報（HitInfo）
- `mesh_id: int`
- `elem_idx: int`
- `t: float`（線分パラメータ t∈[0,1]）
- `pos: ndarray(3,)` 衝突点
- `normal: ndarray(3,)` 衝突要素法線

---

## 5. 外部インターフェース
### 5.1 電場計算API（要求仕様）
```python
calc_electric_field(x, y, z, bem_list) -> (Ex, Ey, Ez)
```
- 入力 `x,y,z`：float または numpy配列（ブロードキャスト可能）
- `bem_list`：`list[BEMMesh]`
- 出力：同形状の `Ex,Ey,Ez`

### 5.2 拡張インターフェース
#### (A) 相互作用モデル（InteractionModel）
- `on_collision(particle, hitinfo) -> (alive, dq_to_surface)`
  - v0.1：吸収のみ（alive=False, dq=q*w）
  - 将来：鏡面/拡散/SEE をここに閉じ込める

#### (B) 表面電荷モデル（SurfaceChargeModel）
- `deposit(hitinfo, dq)`：衝突電荷のバッファリング
- `commit_batch(bem_list)`：バッチ確定（絶縁体は単純加算、将来は再分配/拡散へ差替）

#### (C) 磁場モデル（MagneticFieldModel）
- `B(r, t) -> ndarray(3,)`
  - 既定：`ZeroB`
  - オプション：`UniformB` 等

---

## 6. アルゴリズム仕様
### 6.1 粒子プッシャ（Boris）
1. `E = E(x)`、`B = B(x,t)`
2. Boris更新で `x_{n+1}, v_{n+1}` を計算

> EのみでもBoris形式を維持し、B追加時に実装を変更しない。

### 6.2 衝突判定
- 1ステップを線分 `x_n -> x_{n+1}` とみなし、
  - broad-phase：線分AABB vs 三角形AABB で候補抽出
  - narrow-phase：Möller–Trumbore（segment版）で交差判定
- 複数要素に当たる場合は **最小t（最初の衝突）** を採用
- 複数メッシュがある場合も同様に最小tで一意化

### 6.3 電場計算（Hybrid）
#### 基本方針
- **遠方（Far）**：各要素電荷を重心点電荷として近似
- **近傍（Near）**：近傍要素のみ、三角形を細分したサブ要素中心点電荷和で近似（準直接）
- 合成は差分方式：
  - `E = E_far(all) + Σ_{near}(E_near(subdiv) - E_near(centroid))`

#### 近傍判定
- 代表長さ `h_ref`（例：メッシュの `median(h_elem)`）
- `r_switch = r_switch_factor * h_ref`
- 近傍要素：`||center_j - x|| < r_switch`（必要に応じて要素サイズマージン加算）

#### ソフトニング
- 点電荷近似の発散緩和として `softening = softening_factor * h_ref`
- v0.1 は物理厳密性より数値安定性優先（将来は解析積分や自己項処理へ移行可能）

### 6.4 電荷付与とバッチ更新
- 粒子衝突時：
  - `InteractionModel.on_collision` を呼び、`dq_to_surface` を得る
  - `SurfaceChargeModel.deposit(hit, dq)` へ加算
- バッチ更新（npcls_per_step単位）：
  - `commit_batch` で各要素電荷へ反映
  - 電場キャッシュ（重心・電荷配列など）を更新（rebuild）
  - 相対変化率を計算し停止判定へ

> v0.1 では **npcls_per_step = 「処理した粒子数」** を既定とする。
> 「粒子タイムステップ数」で区切りたい場合はカウンタ位置を変更する。

---

## 7. 反復ループ仕様（擬似コード）
```text
initialize bem_list charges q_j = 0 (or given)
initialize SurfaceChargeModel buffers
initialize Field cache

for each particle p:
  for step in 1..max_step:
    E = field.E(p.x)
    B = B_model.B(p.x, t)
    x_new, v_new = boris(p.x, p.v, E, B, dt)

    hit = first_intersection(p.x, x_new, bem_list)
    if hit exists:
      alive, dq = interaction.on_collision(p, hit)
      surface.deposit(hit, dq)
      p.alive = alive
      if not alive: break
      else: (reflection etc future)
    else:
      p.x, p.v = x_new, v_new

  processed_particles += 1
  batch_count += 1

  if batch_count >= npcls_per_step:
     norm_dq, norm_q = surface.commit_batch(bem_list)
     field.rebuild_cache()
     rel = norm_dq / max(norm_q, q_floor)
     if rel < tol_rel: stop
     batch_count = 0
```

---

## 8. 設定パラメータ（SimConfig）
- 時間積分
  - `dt [s]`
  - `max_step`（粒子の最大ステップ数）
- バッチ
  - `npcls_per_step`（バッチ確定の粒子数単位）
- 停止
  - `tol_rel`（相対変化率閾値）
  - `q_floor`（ゼロ割防止）
- 電場（Hybrid）
  - `use_hybrid: bool`
  - `r_switch_factor: float`
  - `n_sub: int`（近傍三角形細分レベル）
  - `softening_factor: float`

---

## 9. 出力
- 統計（例）
  - `processed_particles`
  - `absorbed`
  - `batches`
  - `last_rel_change`
- 最終状態
  - 各メッシュ要素電荷 `q_j`
- オプション（将来）
  - 粒子軌跡ログ（メモリ/IOコストに注意）
  - 境界電位評価（v0.1はスコープ外、後続で追加）

---

## 10. エラーハンドリング・数値安定性
- 幾何
  - 退化三角形（面積≈0）は事前に除外または例外
  - 法線ゼロの場合の処理（規定：ゼロベクトル、衝突応答は吸収なので影響限定）
- 数値
  - 近傍での発散：softeningで抑制（v0.1）
  - dt過大で壁貫通：線分衝突で検出は可能だが、精度向上には dt 管理が必要
- 性能
  - v0.1は far が O(Nelem) / evaluation
  - 要素数が増える場合は以下に移行：
    - Barnes–Hut / Octree による far 加速
    - BVH による衝突と近傍検索加速
    - 影響係数行列（評価点固定なら有効）

---

## 11. テスト方針（最小）
- 幾何・衝突
  - 既知三角形への垂直入射で衝突点・tが一致
  - 複数候補で最小tが選ばれる
- 電場
  - 単一要素に電荷を置いたとき、遠方で点電荷1/r^2に漸近
  - Hybrid有無で近傍の値が連続的に変化（発散しない）
- 保存量（簡易）
  - 吸収のみの場合、表面へ付与された総電荷 = 粒子電荷総和（w含む）

---

## 12. ロードマップ（拡張方針）
- v0.2：InteractionModel 拡張（鏡面/拡散、SEEの骨格）
- v0.3：性能最適化（Octree/BVH）
- v0.4：導体モード（総電荷拘束＋再分配、BEM連立解の導入）
- v0.5：有限導電率（表面拡散モデル：2D格子化 or グラフラプラシアン）

---
以上
