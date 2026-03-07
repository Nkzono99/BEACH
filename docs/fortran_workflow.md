# Fortran 中心ワークフロー（現行推奨）

このリポジトリは **Fortran 実行系が主**、Python は後処理・可視化を担当します。  
ここでは、初見ユーザーがそのまま実行できる最短手順を先に示します。

## 1. 事前準備

### 1.1 Fortran 実行

```bash
fpm run --profile release --flag "-fopenmp"
```

- 引数なし実行: カレントの `beach.toml` を自動読込
- 設定ファイルを指定する場合:

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

- スレッド数指定:

```bash
OMP_NUM_THREADS=8 fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

- MPI + OpenMP ハイブリッド実行:

```bash
FPM_FC=mpiifort \
fpm run --profile release --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 4" -- examples/beach.toml
```

- `fpm install` を使う場合の推奨ターゲット:

```bash
make install-auto      # ホスト名で profile 自動判定（camphor* -> camphor, それ以外 -> generic）
make install-generic
make install-camphor
```

### 1.2 Python 後処理

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

## 2. 標準実行フロー

1. `beach.toml` を用意（推奨は `docs/fortran_parameter_file.md` の最小例）
2. Fortran を実行
3. `output.dir` に出力されたファイルを確認
4. Python CLI または `Beach` API で可視化

出力される主なファイル:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `charge_history.csv`（`history_stride > 0` のとき）
- `rng_state.txt`
- `macro_residuals.csv`

MPI実行（`world_size > 1`）では乱数状態・残差はrank別に保存されます。

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals_rank00000.csv`, `macro_residuals_rank00001.csv`, ...

## 3. 実行前の負荷見積もり（推奨）

`reservoir_face` / `photo_raycast` を使う場合、バッチあたり粒子数が動的に決まるため事前見積もりを推奨します。

```bash
beach-estimate-workload examples/beach.toml --threads 8
```

- `resolved_batch_duration=...` が表示されます（`batch_duration_step` 指定時）。
- `target_macro_particles_per_batch` を使う場合、この解決済み `batch_duration` で `w_particle` が算出されます。

残差を考慮して続きから見積もる場合:

```bash
beach-estimate-workload examples/beach.toml \
  --threads 8 \
  --macro-residuals outputs/latest/macro_residuals.csv
```

## 4. 再開実行（resume）

設定例:

```toml
[output]
dir = "outputs/latest"
resume = true
```

同じ `output.dir` で再実行すると、`summary.txt` / `charges.csv` / RNG状態ファイルを読み込んで続きから計算します。  
`sim.batch_count` は「今回追加するバッチ数」です。

MPI実行での再開は、`summary.txt` に記録された `mpi_world_size` と現在のrank数が一致している必要があります。

## 5. Python 後処理

### 5.1 CLI

```bash
beach-inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png \
  --potential-self-term area-equivalent

beach-animate-history outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200
```

`quantity=potential` のアニメーション例:

```bash
beach-animate-history outputs/latest \
  --quantity potential \
  --save-gif outputs/latest/potential_history.gif \
  --potential-self-term area-equivalent \
  --total-frames 200
```

### 5.2 Python API

```python
from beach import Beach

beach = Beach("outputs/latest")
print(beach.result.absorbed, beach.result.escaped)

# 引数省略時は outputs/latest を読む
beach_default = Beach()

beach.plot_bar()
beach.plot_mesh()
beach.plot_potential()
beach.animate_mesh("outputs/latest/charge_history.gif", quantity="charge", total_frames=200)
```

- `animate_mesh(output_path=None)` は保存せず `FuncAnimation` を返します。

### 5.3 互換ラッパー

次のスクリプトは CLI の薄いラッパーです。

- `python examples/inspect_fortran_output.py ...`
- `python examples/animate_fortran_history.py ...`
- `python examples/estimate_fortran_workload.py ...`

新規運用では `beach-*` CLI を優先してください。

## 6. MPIテスト（専用）

MPI経路のみを確認する場合:

```bash
FPM_FC=mpiifort \
fpm test --target test_mpi_hybrid \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 2"
```

## 7. 実装挙動で誤解しやすい点

- 実行は `sim.batch_count` 分だけ進みます。
- `sim.tol_rel` は出力監視値で、現行実装では早期終了条件に使いません。
- Fortran 本体の電場は、要素重心への点電荷近似 + `sim.softening` です。
- `sim.use_hybrid` / `r_switch_factor` / `n_sub` / `softening_factor` は現状予約キーです。

## 8. 運用ルール（推奨）

- 物理モデルやランタイム挙動の変更は Fortran 側を正とする。
- 出力フォーマット変更時は `beach/fortran_results.py` と CLI を同時更新する。
- 設定キーを追加・削除したら `docs/fortran_parameter_file.md` を必ず更新する。
