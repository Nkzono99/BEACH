# Fortran 中心ワークフロー（現行推奨）

このプロジェクトは **Fortran 実行系が主**、Python は後処理・可視化を担当します。  
通常利用の推奨運用は、`pip install git+...` で導入した `beach` コマンドを使う方式です。

## 1. 利用者向けセットアップ（推奨）

### 1.1 ツール確認

```bash
make --version
gfortran --version
fpm --version
python --version
```

### 1.2 Git URL から一括インストール（推奨）

```bash
python -m pip install -U pip setuptools wheel
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

`pip install` 時に `make install` が実行され、Python CLI と Fortran 実行バイナリが同時導入されます。
pip 経由のビルドでは既定で `INSTALL_PROFILE=auto` を使い、失敗時は既定で `generic` にフォールバックします。
フォールバックを無効化する場合は `BEACH_PIP_FALLBACK_GENERIC=0` を指定します。

```bash
export PATH="$HOME/.local/bin:$PATH"
```

## 2. 開発に携わる場合

### 2.1 Python 側（編集可能インストール）

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

### 2.2 Fortran 実行系（`make`）

```bash
make
```

必要に応じてプロファイルを明示します。

```bash
make install-generic
make install-camphor
```

### 2.3 `fpm` 直接実行（開発向け）

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

## 3. 実行フロー

1. `beach.toml` を用意（仕様は `docs/fortran_parameter_file.md`）
2. `beach` でシミュレーション実行
3. `output.dir` の出力ファイルを確認
4. Python CLI または `Beach` API で可視化

## 4. 実行コマンド

### 4.1 推奨: `beach`

```bash
beach examples/beach.toml
```

引数なし実行では、カレントディレクトリの `beach.toml` を自動読込します。

### 4.2 スレッド数指定

```bash
OMP_NUM_THREADS=8 beach examples/beach.toml
```

### 4.3 MPI + OpenMP

MPI ビルドをインストール後、ランチャー経由で `beach` を起動します。

```bash
mpirun -n 4 beach examples/beach.toml
```

## 5. 出力ファイル

主な出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `charge_history.csv`（`history_stride > 0`）
- `rng_state.txt`
- `macro_residuals.csv`

MPI実行（`world_size > 1`）では乱数状態・残差は rank 別です。

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals_rank00000.csv`, `macro_residuals_rank00001.csv`, ...

## 6. 実行前の負荷見積もり（推奨）

`reservoir_face` / `photo_raycast` を使う場合、バッチあたり粒子数が動的に決まるため見積もりを推奨します。

```bash
beach-estimate-workload examples/beach.toml --threads 8
```

残差を考慮した rank 局所見積もり:

```bash
beach-estimate-workload examples/beach.toml \
  --threads 8 \
  --mpi-ranks 4 \
  --mpi-rank 0 \
  --macro-residuals outputs/latest/macro_residuals_rank00000.csv
```

## 7. 再開実行（resume）

```toml
[output]
dir = "outputs/latest"
resume = true
```

同じ `output.dir` で `beach` を再実行すると、`summary.txt` / `charges.csv` / RNG状態を読み込んで続きから計算します。  
`sim.batch_count` は「今回追加するバッチ数」です。

MPI再開時は `summary.txt` の `mpi_world_size` と現在の rank 数が一致している必要があります。

## 8. Python 後処理

### 8.1 CLI

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

### 8.2 Python API

```python
from beach import Beach

beach = Beach("outputs/latest")
print(beach.result.absorbed, beach.result.escaped)

beach.plot_bar()
beach.plot_mesh()
beach.plot_potential()
beach.animate_mesh("outputs/latest/charge_history.gif", quantity="charge", total_frames=200)
```

## 9. MPI経路テスト（開発向け）

```bash
FPM_FC=mpiifort \
fpm test --target test_mpi_hybrid \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 2"
```

## 10. 実装挙動で誤解しやすい点

- 実行は `sim.batch_count` 分だけ進みます。
- `sim.tol_rel` は監視値で、現行実装では早期終了条件に使いません。
- Fortran 本体の電場は要素重心点電荷近似 + `sim.softening` です。
- `sim.use_hybrid` / `r_switch_factor` / `n_sub` / `softening_factor` は予約キーです。

camphor向けのMPIジョブ例は `examples/job_scripts/camphor_mpi_hybrid_job.sh` を参照してください。
