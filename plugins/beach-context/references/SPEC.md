# BEACH 仕様書（現行 Fortran 実装）

## 1. 目的

BEACH は、三角形境界要素上の電荷蓄積とテスト粒子追跡を行うシミュレータです。

- 境界は三角形メッシュ
- 粒子運動は Boris push
- 境界衝突時は吸収して要素へ電荷を堆積
- `batch_count` 個の accepted batch を処理し、統計と履歴を更新

計算の主系は Fortran（`src/`, `app/`）で、Python（`beach/`）は設定補助、後処理、可視化を担います。

## 2. スコープ

### 2.1 実装済み

- template / OBJ の三角形メッシュ
- 三角形上の一定面密度を使う `triangle_p0` 静電場
- direct / treecode / FMM field solver
- free-space と 2 軸周期 `periodic2`
- 一様外部電場 `e0` と磁場 `b0`
- Boris 法による粒子更新
- 線分と三角形の最初の交差判定
- insulator accumulation
- free-space に限った mesh 単位の浮遊 conductor 再配分
- open / reflect / redistributed_reflect / periodic の box 境界
- `volume_seed` / `plane_source` / `photo_raycast` と deprecated `reservoir_face`
- species 別の simulation boundary reservoir 流入
- 境界 reservoir の流入補正
- closed photoelectron の局所反射と neutral-return closure
- checkpoint 再開

### 2.2 未実装・予約

- 外部プラズマ／外部シースの自己無撞着ソルバと粒子連成
- 表面導電・拡散モデル
- periodic2 と conductor の併用
- 誘電体分極・誘電境界条件
- 反射や二次電子放出を含む一般化衝突モデル

## 3. データモデル

### 3.1 メッシュ

`mesh_type` は三角形頂点、重心、法線、面積、要素電荷 `q_elem`、mesh ID、および衝突検索用
AABB/grid cache を保持します。OBJ の scale、rotation、offset はこの順に適用します。

要素電荷は総電荷 [C] であり、場計算時に面積で割って一定面密度へ変換します。

### 3.2 粒子

`particles_soa` は位置、速度、電荷、質量、macro weight、生存フラグ、species ID、
`photo_raycast` の放出元要素を保持します。

### 3.3 統計

`sim_stats` は processed / absorbed / escaped / max-step survivor、accepted batch 数、累積物理時間、
適応 batch の棄却数と直近診断、`last_rel_change` を保持します。

`tol_rel` は監視・出力値であり、早期停止条件ではありません。

## 4. 1 batch の計算手順

1. 現在の要素電荷から静電 snapshot を構築する
2. boundary reservoir inflow、`plane_source`、`reservoir_face`、`volume_seed`、`photo_raycast` の粒子を生成する
3. 各粒子を最大 `max_step` 回まで進める
4. mesh hit、box escape、max-step survivor を分類する
5. 吸収電荷と放出元の逆符号電荷を thread-local buffer へ集計する
6. closed PE が有効なら、解決済み帰還先へ neutral-return 補正を適用する
7. MPI reduction 後の候補電荷を検証して commit する
8. accepted batch の統計、履歴、charge ledger、累積時間を更新する

`periodic2.max_nonzero_mode_potential_step > 0` の場合は、設定された `batch_duration` を最大幅として
`h, h/2, h/4, ...` を試します。候補電荷による cached $k\ne0$ 電位変化が上限以下の trial だけを
commit します。棄却時は RNG と macro 粒子数残差を trial 前へ戻し、統計・履歴・ledger を更新しません。

## 5. 場モデル

### 5.1 free-space

`field_solver="direct" | "treecode" | "fmm" | "auto"` を選べます。

- direct: 厳密な free-space panel 和
- treecode: near は厳密 panel、far は monopole
- FMM: near は厳密 panel、far は三角形 moment を使う FMM
- auto: 要素数に応じて direct または FMM

面上電位は連続で、法線電場は `sigma/epsilon0` だけ跳びます。重心電位と principal-value 電場を
自己項として使います。`surface_side` は各 mesh/template で明示します。

### 5.2 periodic2

`[field_boundary] mode="periodic2"` は `[domain] periodic_axes=["x", "y"]` を要求します。
mesh は primitive cell の base element だけを保持し、場と衝突判定が周期 image を扱います。

`field_periodic_far_correction`:

- `none`: 指定した有限 image shell のみ
- `auto`: 互換用に受理し、`none` に正規化
- `cached_kneq0`: 無限周期の非零 Fourier mode を versioned operator cache から評価

型付き `[periodic2]` 設定では、小規模な検証計算に限り
`nonzero_mode_backend="panel_spectral_reference"` と direct solver の組合せも利用できます。
これは非零 Fourier mode の収束参照であり、外部プラズマモデルには依存しません。

`field_periodic_far_correction="none"` は有限画像近似であり、無限周期解ではありません。
`m2l_root_oracle` は削除済みです。production の無限周期計算には `cached_kneq0` を使用します。
cache fingerprint は generator version を含み、物理的 zero mode の評価は高さ breakpoint の
二分探索により 1 点あたり O(log n) です。

`cached_kneq0` の物理的 $k=0$ は `periodic2.zero_mode_policy="exclude_k0"` と
`lower_boundary_model="symmetric_vacuum" | "e_bottom_zero"` で別に構築し、場へ一度だけ加えます。
外部プラズマ応答は合成しません。

### 5.3 領域・粒子境界・reservoir

公開 TOML は topology、場、外向き粒子作用、外部 reservoir 条件、species 別流入を分離します。

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape" # または "potential_barrier"

[reservoir]
inflow_model = "source_vdf" # または "infinity_barrier"
phi_infty = 0.0
face_potential_grid_n = 5

[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0

[particles.species.boundary_inflow]
z_high = "reservoir"
```

- `boundary_inflow`: 非周期 box 面全体から外部 reservoir VDF を流入させる
- `source_vdf`: 指定した boundary VDF を補正せず境界から注入する
- `infinity_barrier`: `phi_infty` と流入面平均電位の scalar barrier で法線 VDF を補正する
- `escape`: open 面到達で粒子を消滅する
- `potential_barrier`: event 位置の局所電位と法線運動エネルギーから反射／escape を判定する

`boundary_inflow` は外向き粒子作用を上書きしません。周期面への reservoir 流入、および有効な外向き作用が
`open` でない流入面は拒否します。
同じ species では `source_mode="volume_seed"` とだけ併用でき、`plane_source`、`photo_raycast`、
`reservoir_face` との併用は fail closed です。

外部プラズマ profile、外部領域の particle transport、delayed return queue は現行スコープ外です。
削除済みの `[outer_plasma]`、`[coupling]`、`[external_boundary]` は unknown input として拒否します。

## 6. 粒子前進・衝突・box 境界

public な粒子状態は同一時刻の `x, v` です。空間電場は予測中点
`x_mid = x + 0.5*v*dt` で評価し、Boris の速度更新後に
`x_new = x + 0.5*(v + v_new)*dt` とします。

mesh 衝突は Möller–Trumbore 法で最小 `t` を選びます。AABB と、既定では一様 grid + 3D-DDA で
候補を絞ります。非有限値や進行不能は hit なしへ黙って変換せず、明示 status で停止します。

`domain.periodic_axes` は全粒子種に共通の topology です。非周期面の粒子 action は
global `[particle_boundary]` と species 別 `[particles.species.boundary]` で指定します。

- `open`: escape
- `reflect`: 法線速度を反転し、接線速度と接線位置を維持
- `redistributed_reflect`: `reflect` と同じ速度作用に加え、event 面内の位置を一様再配置
- `periodic`: 反対側へ wrap

単一面の `redistributed_reflect` は面内 2 軸を box span の両端 guard を除く範囲から再標本化します。
edge / corner の同時 event では、
event mask に含まれない軸だけを再標本化し、event 軸の座標は内側 guard に置きます。mesh hit と最初の
box event は chord parameter で順序付けます。reflect / redistributed_reflect / periodic 後の残り時間は
同じ Boris 規約で再積分し、1 local step あたり最大 8 box event を許します。

species 別の6面値は `inherit | open | reflect | redistributed_reflect` です。`inherit` はglobal粒子 actionを使います。
周期面はspeciesから上書きできません。

## 7. 粒子 source と境界流入

### 7.1 `volume_seed`

各 species で `npcls_per_step` 個を batch ごとに生成します。

`boundary_inflow` を持つ species では `npcls_per_step=0` を許容します。

### 7.2 `plane_source`

`pos_low` と `pos_high` で box 内部の axis-aligned 矩形面を定義し、`source_normal` 方向へ一方向 flux を
生成します。矩形面はちょうど 1 軸で zero thickness、残る 2 軸で正の面積を持ちます。法線座標は
box 境界より厳密に内側で、接線範囲は box 内に置きます。
`source_normal` は zero-thickness 軸に沿う非ゼロベクトルで、実装内部で単位ベクトルへ正規化します。

Maxwell reservoir または速度 grid の flux、面積、batch duration から macro 粒子数を決めます。
`reservoir.inflow_model`、`phi_infty`、`face_potential_grid_n` は内部平面へ適用しません。
任意の内部面では外部 plasma 側と上流基準電位の対応が一意でないため、流入側の potential/barrier 補正は
新しい設定では simulation boundary に結び付いた `boundary_inflow` が所有します。非推奨の
`reservoir_face` は既存 case の数値挙動を保つため、従来の補正を互換動作として維持します。

### 7.3 species 別 boundary reservoir 流入

`[particles.species.boundary_inflow]` の 6 面キーへ `reservoir` を指定し、非周期 box 面全体から流入させます。
複数面を指定でき、macro 粒子の端数は species と face の組ごとに batch 間で繰り越します。
MPI では global 個数を rank へ分配します。

流入 flux と面積、batch duration から macro 粒子数を決めます。`target_macro_particles_per_batch` 指定時は
macro weight を自動解決します。

`reservoir.inflow_model="infinity_barrier"` では各流入面の `N x N` 電位標本を使い、
`reservoir.phi_infty` とのエネルギー差から法線速度 cutoff と到達速度を同じ写像で補正します。

### 7.4 `reservoir_face`（deprecated）

既存 case の `inject_face` と `pos_low` / `pos_high` による box 面開口の挙動を維持します。
`boundary_inflow` または `plane_source` へ暗黙変換しません。外部 plasma 条件の新規 case は
`boundary_inflow`、内部矩形面は `plane_source` を使います。

### 7.5 `photo_raycast`

指定方向へ `rays_per_batch` 本を飛ばし、各 ray の最初の mesh hit からだけ放出します。
`deposit_opposite_charge_on_emit=true` では放出元へ逆符号電荷を加えます。

closed PE は次の組合せです。

- 負電荷 `photo_raycast`
- `deposit_opposite_charge_on_emit=true`
- species の `inject_face` に対する有効粒子 action が `reflect` または `redistributed_reflect`
- `surface_charge_closure="neutral_return"`

放出電荷を $S<0$、解決済み再吸収を $R<0$、max-step 未解決を $U<0$ とし、`S=R+U` を確認した上で、
帰還先への deposit を `S/R` 倍します。これにより PE species の表面総電荷増分を 0 にしつつ、
帰還先分布による表面内再分配を保持します。

未解決率 `U/S` が 0.05 を超える、`R=0`、実 escape、soft discard、terminal 分類不整合のいずれかでは
補正せず fail closed とします。

## 8. 実行制御

- 新規実行は `sim.batch_count` 個の accepted batch を処理する
- 再開は checkpoint の batch 数から `sim.batch_count` まで処理する
- 各粒子の前進上限は `sim.max_step`
- `tol_rel` は早期停止に使わない

## 9. 出力と再開

主要出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `mesh_potential.csv`（設定時）
- `charge_history.csv` / `potential_history.csv` / `top_reference_history.csv`（設定時）
- `charge_ledger.csv`（ledger がある場合）
- `rng_state.txt`、MPI では `rng_state_rankNNNNN.txt`
- `macro_residuals.csv`
- `performance_profile.csv`（`BEACH_PROFILE=1`）

出力ファイルの条件は `schemas/beach.output-manifest.json` を正本とします。

再開時の必須ファイルは `summary.txt`、`charges.csv`、対応する RNG state です。
`macro_residuals.csv` と `charge_ledger.csv` は状態が記録されている場合に復元します。

`output.checkpoint_stride > 0` では accepted batch の commit 後だけ定期 checkpoint を作ります。
`checkpoints/slot0` と `slot1` を交互に使い、全 rank の状態を書き終えてから `checkpoint_latest.txt` を
原子的に切り替えます。再開時は直下の最終出力と active slot のうち、必須ファイルが揃い
`batches` が最大のものを選びます。`checkpoint_stride=0` でも正常終了時の最終 checkpoint は出力します。

`summary.txt` の checkpoint schema と model / ordered mesh / ordered species fingerprint を照合します。
schema v6 の `macro_residuals.csv` は `species_idx,face,residual` を持ち、`face=0` は従来 source、
`1..6` は boundary face です。旧 2 列形式は読み込み互換です。
globalまたはspecies境界のいずれかで`redistributed_reflect`を使う場合だけ、model fingerprintへ
`sim.rng_seed`と乱数契約識別子`redistributed_reflect_rng_v1`を含めます。通常境界だけの既存fingerprintは変更しません。
必須ファイルの欠落、world size の不一致、非有限値、species 数や mesh 要素数の不一致は
新規実行へ fallback せず停止します。

## 10. 設計方針

- v1.0 の基本系は insulator accumulation
- correctness を優先し、性能機能は明示設定で有効化する
- 実装されていない物理モデルを設定互換性のためだけに保持しない
- public API の変更時は schema、例、日英ドキュメント、テストを同時更新する

## 11. Python 後処理

Python package `beach` は Fortran 出力を読み込み、電位・電場・Coulomb 力・可視化を提供します。
任意点の再評価は三角形 geometry を渡す native `triangle_p0` kernel を使います。

`cached_kneq0` は非零 mode だけの低水準 operator であり、total field として使うには Fortran と同じ
物理的 zero mode の合成が必要です。設定が見つからない total-field 再計算は free-space へ暗黙 fallback しません。
