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
- open / reflect / periodic の box 境界
- `volume_seed` / `reservoir_face` / `photo_raycast`
- 局所 reservoir の流入補正
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
2. `reservoir_face`、`volume_seed`、`photo_raycast` の粒子を生成する
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

`field_bc_mode="periodic2"` はちょうど 2 軸を周期、残り 1 軸を open とします。mesh は primitive cell の
base element だけを保持し、場と衝突判定が周期 image を扱います。

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

### 5.3 外部境界契約

公開 TOML の `[external_boundary]` は次だけを受理します。

```toml
[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "source_vdf" # または "infinity_barrier"

[external_boundary.ordinary_open]
model = "escape" # または "potential_barrier"
```

- `source_vdf`: 指定した boundary VDF をそのまま局所注入する
- `infinity_barrier`: `phi_infty` と注入面平均電位の scalar barrier で法線 VDF を補正する
- `escape`: open 面到達で粒子を消滅する
- `potential_barrier`: event 位置の局所電位と法線運動エネルギーから反射／escape を判定する

外部プラズマ profile、外部領域の particle transport、delayed return queue は現行スコープ外です。
削除済みの `[outer_plasma]`、`[coupling]`、旧 selector key は unknown input として拒否します。

## 6. 粒子前進・衝突・box 境界

public な粒子状態は同一時刻の `x, v` です。空間電場は予測中点
`x_mid = x + 0.5*v*dt` で評価し、Boris の速度更新後に
`x_new = x + 0.5*(v + v_new)*dt` とします。

mesh 衝突は Möller–Trumbore 法で最小 `t` を選びます。AABB と、既定では一様 grid + 3D-DDA で
候補を絞ります。非有限値や進行不能は hit なしへ黙って変換せず、明示 status で停止します。

box action:

- `open`: escape
- `reflect`: 法線速度を反転
- `periodic`: 反対側へ wrap

mesh hit と最初の box event は chord parameter で順序付けます。reflect / periodic 後の残り時間は
同じ Boris 規約で再積分し、1 local step あたり最大 8 box event を許します。

`[[particles.species]].z_high_boundary="reflect"` は、その species の z-high だけを局所反射へ置き換えます。
`sim.use_box=true`、global z-high が open、外部粒子輸送なしの場合だけ有効です。

## 7. 注入モード

### 7.1 `volume_seed`

各 species で `npcls_per_step` 個を batch ごとに生成します。

### 7.2 `reservoir_face`

流入 flux と面積、batch duration から macro 粒子数を決めます。`target_macro_particles_per_batch` 指定時は
macro weight を自動解決します。端数は batch 間で繰り越し、MPI では global 個数を rank へ分配します。

`inflow_model="infinity_barrier"` では注入面の `N x N` 電位標本を使い、
`phi_infty` とのエネルギー差から法線速度 cutoff と到達速度を同じ写像で補正します。

### 7.3 `photo_raycast`

指定方向へ `rays_per_batch` 本を飛ばし、各 ray の最初の mesh hit からだけ放出します。
`deposit_opposite_charge_on_emit=true` では放出元へ逆符号電荷を加えます。

closed PE は次の組合せです。

- 負電荷 `photo_raycast`
- `deposit_opposite_charge_on_emit=true`
- `z_high_boundary="reflect"`
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

`summary.txt` の checkpoint schema と model / ordered mesh / ordered species fingerprint を照合します。
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
