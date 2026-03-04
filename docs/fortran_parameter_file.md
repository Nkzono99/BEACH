# Fortranパラメータファイル仕様（TOML）

Fortran 実行時は、`fpm run ... -- path/to/config.toml` で設定を読み込めます。

- 例: `examples/fortran_config.toml`
- 既定値は `src/bem_app_config.f90` の `default_app_config` で定義

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
| `batch_count` | int | `1` | 固定で実行するバッチ数 |
| `max_step` | int | `400` | 粒子あたり最大ステップ |
| `tol_rel` | float | `1.0e-8` | 相対変化の監視値（早期終了には使わない） |
| `q_floor` | float | 実装既定値 | ゼロ割防止下限 |
| `softening` | float | `1.0e-6` | 点電荷近似のソフトニング |
| `use_hybrid` | bool | 実装既定値 | Hybrid場計算の有効化 |
| `r_switch_factor` | float | 実装既定値 | 近傍切替半径係数 |
| `n_sub` | int | 実装既定値 | 近傍細分レベル |
| `softening_factor` | float | 実装既定値 | `h_ref` 基準ソフトニング係数 |
| `b0` | float[3] | `[0,0,0]` | 一様磁場[T] |
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
各バッチでは、全 species について `npcls_per_step` 個ずつ毎回生成されます。
生成順はラウンドロビン（species1→species2→...）なので、1バッチ内に複数種が混在します。
各エントリは独立した粒子注入条件です。未指定キーは実装既定値で補完されます。

| キー | 型 | 説明 |
|---|---|---|
| `enabled` | bool | 有効/無効（既定: `true`） |
| `npcls_per_step` | int | その粒子種が1バッチで生成する粒子数 |
| `q_particle` | float | 粒子電荷[C] |
| `m_particle` | float | 粒子質量[kg] |
| `w_particle` | float | superparticle 重み |
| `pos_low` | float[3] | 初期位置下限[m] |
| `pos_high` | float[3] | 初期位置上限[m] |
| `drift_velocity` | float[3] | ドリフト速度[m/s] |
| `temperature_k` | float | 温度[K] |

粒子数の計算は以下です。

- 1バッチあたり総粒子数 = `sum(enabled species npcls_per_step)`
- 総投入粒子数 = `batch_count * sum(enabled species npcls_per_step)`

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

---

出力ディレクトリには以下のファイルが生成されます。

- `summary.txt`: 集計統計（`escaped_boundary` / `survived_max_step` を含む）
- `charges.csv`: 最終要素電荷
- `mesh_triangles.csv`: 要素三角形頂点と最終電荷
- `charge_history.csv`: 指定した `history_stride` 間隔で逐次書き出される要素電荷履歴（時間発展）

## 4. 注意点

- `[[particles.species]]` は1件以上必須です。
- `sim.batch_count` は固定実行回数です。`tol_rel` を下回っても早期終了しません。
- `rng_seed` は `[sim]` に置きます。
- `[[particles.species]]` を使うと、粒子は種ごとラウンドロビンで混在して初期化されます。
- `[[mesh.templates]]` は同じ `kind` を複数回指定可能です。
- `n_templates` は先頭から読み込む最大件数として扱われます。
- `mode="auto"` では `obj_path` が存在すれば OBJ、なければ template を使用します。

---

## 5. 典型的な実行

```bash
fpm run --profile release --flag "-fopenmp" -- examples/fortran_config.toml
python examples/inspect_fortran_output.py outputs/latest
```
