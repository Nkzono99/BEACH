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

### 4.4 性能計測つき実行

粗いフェーズ計測だけを有効にする場合:

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach examples/beach.toml
```

粒子ループ内の `field` / `push` / `collision` まで詳細計測する場合:

```bash
BEACH_PROFILE=1 BEACH_PROFILE_DETAIL=1 OMP_NUM_THREADS=8 beach examples/beach.toml
```

`BEACH_PROFILE_DETAIL=1` は時刻取得回数が増えるため、粗粒度計測より遅くなります。  
詳細計測の `particle_field_eval` / `particle_push` / `particle_collision` と `summary.txt` の
`field_time_s` / `push_time_s` / `collision_time_s` は、他のフェーズ計測と比較しやすいよう
thread 内訳の和ではなく wall-clock ベースで集計されます。  
スケーリング比較には `performance_profile.csv` の `simulation_total` 行にある `rank_max_s` を使うのが推奨です。

可視化例:

```bash
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

## 5. 出力ファイル

主な出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv`（`history_stride > 0`）
- `performance_profile.csv`（`BEACH_PROFILE=1`）
- `rng_state.txt`
- `macro_residuals.csv`

`mesh_triangles.csv` には要素ごとの `mesh_id` 列が含まれ、`mesh_sources.csv` で `mesh_id` と元メッシュ設定（template kind / 要素数）を対応付けます。

MPI実行（`world_size > 1`）では乱数状態・残差は rank 別です。

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals_rank00000.csv`, `macro_residuals_rank00001.csv`, ...

## 6. 実行前の負荷見積もり（推奨）

`reservoir_face` / `photo_raycast` を使う場合、バッチあたり粒子数が動的に決まるため見積もりを推奨します。

```bash
beachx workload examples/beach.toml --threads 8
```

残差を考慮した rank 局所見積もり:

```bash
beachx workload examples/beach.toml \
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
beachx inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png \
  --potential-self-term area-equivalent

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200

beachx slices outputs/latest \
  --grid-n 200 \
  --vmin -20 --vmax 20 \
  --save outputs/latest/potential_slices.png
```

旧 alias の `beach-inspect` / `beach-animate-history` / `beach-plot-potential-slices` /
`beach-estimate-workload` / `beach-plot-performance-profile` も当面は利用できますが、非推奨です。

### 8.2 Python API

```python
from beach import Beach

beach = Beach("outputs/latest")
print(beach.result.absorbed, beach.result.escaped)
# history は常に遅延読込
history_step10 = beach.result.history_at(10)  # 特定 batch の要素電荷を取得
if beach.result.history is not None:
    print(beach.result.history.batch_indices)  # 利用可能な batch 一覧

beach.plot_bar()
beach.plot_mesh()
beach.plot_potential()
beach.plot_potential_slices(
    box_min=[0.0, 0.0, 0.0],
    box_max=[1.0, 1.0, 10.0],
    grid_n=200,
    vmin=-20.0,
    vmax=20.0,
)
beach.animate_mesh("outputs/latest/charge_history.gif", quantity="charge", total_frames=200)

mesh1 = beach.get_mesh(1)
mesh2, mesh3 = beach.get_mesh(2, 3)
mesh1_step10 = beach.get_mesh(1, step=10)
charge_step10 = beach.get_mesh_charge(1, step=10)

interaction = beach.calc_coulomb(target=[mesh1, mesh2], source=[mesh3], step=10)
print(interaction.force_on_a_N, interaction.torque_on_a_Nm)
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
