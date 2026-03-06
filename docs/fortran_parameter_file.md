# Fortranパラメータファイル仕様（TOML）

Fortran 実行時は、`fpm run ... -- path/to/config.toml` で設定を読み込めます。引数なし実行でもカレントディレクトリの `beach.toml` を自動読込します。

- 例: `examples/beach.toml`
- 既定値は `src/bem_app_config.f90` の `default_app_config` で定義
- 互換目的で `fortran_config.toml` も自動読込対象ですが、`beach.toml` へ移行してください

---

## 1. 最小例

```toml
[sim]
dt = 1.0e-9
rng_seed = 12345
batch_count = 4
max_step = 400

[[particles.species]]
npcls_per_step = 64

[mesh]
mode = "template"
n_templates = 1

[[mesh.templates]]
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 1
ny = 1
center = [0.0, 0.0, 0.0]

[output]
write_files = true
dir = "outputs/latest"
```

---

## 2. セクション一覧

- `[sim]`: 時間積分・バッチ回数・電磁場パラメータ
- `[[particles.species]]`: 粒子種ごとの注入条件（1件以上必須）
- `[mesh]`: メッシュ入力モード
- `[[mesh.templates]]`: テンプレート形状定義（複数可）
- `[output]`: ファイル出力設定

---

## 3. パラメータ詳細

### 3.1 `[sim]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `dt` | float | `1.0e-9` | 時間刻み[s] |
| `rng_seed` | int | `12345` | 粒子注入乱数のシード |
| `batch_count` | int | `1` | 1回の実行で進めるバッチ数 |
| `batch_duration` | float | `0.0` | 1バッチが表す物理時間[s]（`reservoir_face` 使用時に必要。`target_npcls_species1` と同時指定不可） |
| `target_npcls_species1` | int | `0` | `species[1]` の目標マクロ粒子数/バッチ。指定時は `batch_duration` を自動決定（`species[1]` は `reservoir_face` 必須） |
| `max_step` | int | `400` | 粒子あたり最大ステップ |
| `tol_rel` | float | `1.0e-8` | 相対変化の監視値（早期終了には使わない） |
| `q_floor` | float | 実装既定値 | ゼロ割防止下限 |
| `softening` | float | `1.0e-6` | 点電荷近似のソフトニング |
| `use_hybrid` | bool | 実装既定値 | Hybrid場計算の有効化 |
| `r_switch_factor` | float | 実装既定値 | 近傍切替半径係数 |
| `n_sub` | int | 実装既定値 | 近傍細分レベル |
| `softening_factor` | float | 実装既定値 | `h_ref` 基準ソフトニング係数 |
| `b0` | float[3] | `[0,0,0]` | 一様磁場[T] |
| `reservoir_potential_model` | string | `"none"` | `reservoir_face` の電位補正モデル（`none` / `infinity_barrier`） |
| `phi_infty` | float | `0.0` | 無限遠基準電位 [V]（`infinity_barrier` で使用） |
| `injection_face_phi_grid_n` | int | `3` | 注入開口面平均電位の格子分割数（`N x N`） |
| `use_box` | bool | `false` | シミュレーションボックス境界を有効化 |
| `box_min` | float[3] | `[-1,-1,-1]` | ボックス下限座標[m] |
| `box_max` | float[3] | `[1,1,1]` | ボックス上限座標[m] |
| `bc_x_low` | string | `"open"` | x低側境界 (`open`/`reflect`/`periodic`) |
| `bc_x_high` | string | `"open"` | x高側境界 (`open`/`reflect`/`periodic`) |
| `bc_y_low` | string | `"open"` | y低側境界 (`open`/`reflect`/`periodic`) |
| `bc_y_high` | string | `"open"` | y高側境界 (`open`/`reflect`/`periodic`) |
| `bc_z_low` | string | `"open"` | z低側境界 (`open`/`reflect`/`periodic`) |
| `bc_z_high` | string | `"open"` | z高側境界 (`open`/`reflect`/`periodic`) |

### 3.2 `[[particles.species]]`

単一種でも `[[particles.species]]` を1件以上定義してください。
各エントリは独立した粒子注入条件です。未指定キーは実装既定値で補完されます。
`source_mode="volume_seed"` では従来どおり `npcls_per_step` 個を毎バッチ生成します。
`source_mode="reservoir_face"` では、上流プラズマ条件から物理流入量を計算し、`w_particle` をマクロ粒子重みとして `npcls_per_step` を内部で自動決定します。

| キー | 型 | 説明 |
|---|---|---|
| `enabled` | bool | 有効/無効（既定: `true`） |
| `source_mode` | string | `volume_seed`（既定） / `reservoir_face` |
| `npcls_per_step` | int | `volume_seed` 時の1バッチ生成粒子数 |
| `number_density_cm3` | float | `reservoir_face` 時の上流密度[`cm^-3`] |
| `number_density_m3` | float | `reservoir_face` 時の上流密度[`m^-3`]（`number_density_cm3` と排他） |
| `q_particle` | float | 粒子電荷[C] |
| `m_particle` | float | 粒子質量[kg] |
| `w_particle` | float | superparticle 重み（`reservoir_face` では 1マクロ粒子が代表する実粒子数） |
| `pos_low` | float[3] | 初期位置下限[m] |
| `pos_high` | float[3] | 初期位置上限[m] |
| `drift_velocity` | float[3] | ドリフト速度[m/s] |
| `temperature_k` | float | 温度[K] |
| `temperature_ev` | float | 温度[eV]（`temperature_k` と排他） |
| `inject_face` | string | `reservoir_face` 時の注入面（`x_low` / `x_high` / `y_low` / `y_high` / `z_low` / `z_high`） |

粒子数の計算は以下です。

- `volume_seed` の1バッチあたり総粒子数 = `sum(enabled volume_seed species npcls_per_step)`
- `volume_seed` の総投入粒子数 = `batch_count * sum(enabled volume_seed species npcls_per_step)`
- `reservoir_face` の期待マクロ粒子数（種 `s`）:
  - `n_macro_expected_s = gamma_in_s * A_s * batch_duration / w_particle_s`
  - `gamma_in_s = n_s * (sigma_s * phi(alpha_s) + u_n_s * Phi(alpha_s))`
  - `alpha_s = u_n_s / sigma_s`
  - `sigma_s = sqrt(k_B * T_s / m_s)`
  - `u_n_s = drift_velocity_s ・ inward_normal_s`
- 実際のバッチ注入数（残差繰越込み）:
  - `macro_budget_s = residual_s + n_macro_expected_s`
  - `n_macro_s = floor(max(macro_budget_s, 0))`
  - `residual_s <- macro_budget_s - n_macro_s`
- 上式は `compute_macro_particles_for_batch` の実装と一致します。

`sim.reservoir_potential_model = "infinity_barrier"` を有効化した場合、`reservoir_face` の流入判定と法線速度サンプルに電位差補正が入ります。

- `phi_face`: 注入開口面上 `N x N`（`N = injection_face_phi_grid_n`）格子点で評価した平均電位
- `delta_phi_s = phi_face - phi_infty`
- `barrier_s = 2 * q_s * delta_phi_s / m_s`
- `v_inf_min_s = sqrt(max(barrier_s, 0))`
- `gamma_in_s` は `v_inf >= v_inf_min_s` の切断付き flux-weighted 積分で評価
- サンプルした `v_inf` から注入面法線速度 `v_face = sqrt(max(v_inf^2 - barrier_s, 0))` を与える

`sim.target_npcls_species1` の `batch_duration` 自動決定式は従来どおり（未補正係数）で、`infinity_barrier` 併用時は近似になります。

`target_npcls_species1` を使う場合、`species[1]` を基準に `batch_duration` を次式で自動決定します。

- `batch_duration = target_npcls_species1 * w_1 / (gamma_in_1 * A_1)`
- `species[1]` は `enabled = true` かつ `source_mode = "reservoir_face"` が必須です。
- `batch_duration` と `target_npcls_species1` の同時指定はエラーです。
- 他の `reservoir_face` 種は同じ `batch_duration` で注入されるため、長期平均ではフラックス比が保持されます（整数化誤差は残差繰越で吸収）。

`reservoir_face` では次が必須です。

- `sim.use_box = true`
- `sim.batch_duration > 0`
- `number_density_cm3` または `number_density_m3`
- `inject_face`
- `w_particle > 0`

`target_npcls_species1` を使う場合の追加要件:

- `sim.target_npcls_species1 > 0`
- `sim.batch_duration` を同時に指定しない
- `particles.species[1]` が有効かつ `source_mode = "reservoir_face"`

`pos_low` / `pos_high` は注入面上の矩形開口を表します。法線方向の座標は `inject_face` で指定した箱境界と一致している必要があります。

例: 5/cc, 10eV, 400km/s の上流電子を上面 (`z_high`) から流入させる設定

```toml
[sim]
batch_duration = 2.0e-7
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
w_particle = 1.0e5
inject_face = "z_high"
pos_low = [0.0, 0.0, 4.0]
pos_high = [1.0, 1.0, 4.0]
drift_velocity = [0.0, 0.0, -4.0e5]
```

`target_npcls_species1` で `batch_duration` を自動決定する例:

```toml
[sim]
target_npcls_species1 = 300
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
w_particle = 1.0e5
inject_face = "z_high"
pos_low = [0.0, 0.0, 4.0]
pos_high = [1.0, 1.0, 4.0]
drift_velocity = [0.0, 0.0, -4.0e5]
```

旧キー `n_particles` は廃止されています。
旧 `[particles]` 直下の既定値キーも廃止されました。`rng_seed` は `[sim]` に置き、粒子物性は各 `[[particles.species]]` に記述してください。

### 3.3 `[mesh]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `mode` | string | `"auto"` | `auto` / `obj` / `template` |
| `obj_path` | string | `"examples/simple_plate.obj"` | OBJ パス |
| `n_templates` | int | `1` | 使用する template エントリ数 |

### 3.4 `[[mesh.templates]]`

共通キー:

- `enabled` (bool): 有効/無効
- `kind` (string): `plane`, `box`, `cylinder`, `sphere`
- `center` (float[3]): 中心座標

`kind` ごとの主キー:

- `plane`: `size_x`, `size_y`, `nx`, `ny`
- `box`: `size`, `nx`, `ny`, `nz`
- `cylinder`: `radius`, `height`, `n_theta`, `n_z`, `cap`
- `sphere`: `radius`, `n_lon`, `n_lat`

### 3.5 `[output]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `write_files` | bool | `true` | 結果ファイルを書き出すか |
| `dir` | string | `"outputs/latest"` | 出力先ディレクトリ |
| `history_stride` | int | `1` | 何バッチごとに `charge_history.csv` へ履歴を書き出すか（1=毎バッチ） |
| `resume` | bool | `false` | `dir` に既存の `summary.txt` / `charges.csv` / `rng_state.txt`（あれば `macro_residuals.csv` も） があれば、その続きから再開する |

---

出力ディレクトリには以下のファイルが生成されます。

- `summary.txt`: 集計統計（`escaped_boundary` / `survived_max_step` を含む）
- `charges.csv`: 最終要素電荷
- `mesh_triangles.csv`: 要素三角形頂点と最終電荷
- `charge_history.csv`: 指定した `history_stride` 間隔で逐次書き出される要素電荷履歴（時間発展）
- `rng_state.txt`: 次回 `resume = true` で再開するための Fortran 乱数状態
- `macro_residuals.csv`: `reservoir_face` の自動粒子数計算で残った端数状態

## 4. 注意点

- `[[particles.species]]` は1件以上必須です。
- `sim.batch_count` は1回の実行で追加するバッチ数です。`tol_rel` を下回っても早期終了しません。
- `rng_seed` は `[sim]` に置きます。
- `[[particles.species]]` を使うと、粒子は種ごとラウンドロビンで混在して初期化されます。
- `[[mesh.templates]]` は同じ `kind` を複数回指定可能です。
- `n_templates` は先頭から読み込む最大件数として扱われます。
- `mode="auto"` では `obj_path` が存在すれば OBJ、なければ template を使用します。
- `output.resume = true` の場合、同じ `output.dir` にある前回のチェックポイントを読み込み、今回の `sim.batch_count` 分だけ追加でバッチを進めます。
- 再開時は `charge_history.csv` に追記します。既存履歴も残したいので、通常は同じ `output.dir` を使ってください。

---

## 5. 典型的な実行

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
python examples/inspect_fortran_output.py outputs/latest
```
