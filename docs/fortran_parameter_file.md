# Fortranパラメータファイル仕様（TOML）

Fortran 実行時は、`fpm run ... -- path/to/config.toml` で設定を読み込めます。

- 例: `examples/fortran_config.toml`
- 既定値は `src/bem_app_config.f90` の `default_app_config` で定義

---

## 1. 最小例

```toml
[sim]
dt = 1.0e-9
npcls_per_step = 64
max_step = 400

[particles]
n_particles = 256

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

- `[sim]`: 時間積分・停止条件・電磁場パラメータ
- `[particles]`: 注入粒子条件
- `[mesh]`: メッシュ入力モード
- `[[mesh.templates]]`: テンプレート形状定義（複数可）
- `[output]`: ファイル出力設定

---

## 3. パラメータ詳細

### 3.1 `[sim]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `dt` | float | `1.0e-9` | 時間刻み[s] |
| `npcls_per_step` | int | `64` | バッチ確定粒子数 |
| `max_step` | int | `400` | 粒子あたり最大ステップ |
| `tol_rel` | float | `1.0e-8` | 相対変化停止閾値 |
| `q_floor` | float | 実装既定値 | ゼロ割防止下限 |
| `softening` | float | `1.0e-6` | 点電荷近似のソフトニング |
| `use_hybrid` | bool | 実装既定値 | Hybrid場計算の有効化 |
| `r_switch_factor` | float | 実装既定値 | 近傍切替半径係数 |
| `n_sub` | int | 実装既定値 | 近傍細分レベル |
| `softening_factor` | float | 実装既定値 | `h_ref` 基準ソフトニング係数 |
| `b0` | float[3] | `[0,0,0]` | 一様磁場[T] |

### 3.2 `[particles]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `rng_seed` | int | `12345` | 乱数シード |
| `n_particles` | int | `256` | 粒子数 |
| `q_particle` | float | `-1.602176634e-19` | 粒子電荷[C] |
| `m_particle` | float | `9.10938356e-31` | 粒子質量[kg] |
| `w_particle` | float | `1.0` | superparticle 重み |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | 初期位置下限[m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | 初期位置上限[m] |
| `drift_velocity` | float[3] | `[0,0,-8.0e5]` | ドリフト速度[m/s] |
| `temperature_k` | float | `2.0e4` | 温度[K] |

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

---

出力ディレクトリには以下のファイルが生成されます。

- `summary.txt`: 集計統計
- `charges.csv`: 最終要素電荷
- `mesh_triangles.csv`: 要素三角形頂点と最終電荷
- `charge_history.csv`: バッチごとの要素電荷履歴（時間発展）

## 4. 注意点

- `[[mesh.templates]]` は同じ `kind` を複数回指定可能です。
- `n_templates` は先頭から読み込む最大件数として扱われます。
- `mode="auto"` では `obj_path` が存在すれば OBJ、なければ template を使用します。

---

## 5. 典型的な実行

```bash
fpm run --profile release --flag "-fopenmp" -- examples/fortran_config.toml
python examples/inspect_fortran_output.py outputs/latest
```
