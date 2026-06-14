title: Fortran 中心ワークフロー（現行推奨）

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
make check
make run CONFIG=examples/beach.toml
```

開発中の標準確認は `make check` です。`BEACH_VERSION_MODE=dev` を使って
Fortran 側へ渡す version 文字列を `1.2.0-dev` のように固定するため、git hash が変わっても
fpm の compile-flag hash が変わらず、差分コンパイルを再利用できます。

`make build` と `make install` は既定で git hash 付き version を埋め込みます。必要なら version mode を明示します。

```bash
make build VERSION_MODE=dev
make build VERSION_MODE=plain
make build VERSION_MODE=git
```

インストールプロファイルは必要に応じて明示します。

```bash
make install-generic
make install-camphor
```

### 2.3 `fpm` 直接実行（低レベル確認向け）

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

通常の開発では `build.sh` 経由の `make run` / `make check` を優先してください。
`build.sh` が `__BEACH_VERSION__` と `__BEACH_VERSION_MODE__` を安定した形で渡します。

### 2.4 テスト

```bash
make test-l0      # L0: static/schema/build check
make test         # L1: normal development loop
make test-l2      # L2: contract/integration
make test-l3      # L3: heavy/release gate
make test-heavy   # heavy Fortran targets only
make test-full    # unfiltered fpm test
```

BEACH のテストは開発ループ向けに階層化しています。

- L0: `git diff --check`、JSON schema parse check、`make check`
- L1: L0 + Python tests + 軽量 Fortran test targets（`make test` / `make test-l1`）
- L2: L1 + contract/integration targets（C field-kernel contract など）
- L3: L2 + heavy FMM targets / full fpm suite（release gate / nightly / main 統合前）

`make test-fortran` は軽量 Fortran target の alias です。重い FMM 系
（`test_dynamics_fmm`, `test_coulomb_fmm_core_basic`, `test_coulomb_fmm_core_periodic`,
`test_periodic2_flat_oracle_diag`）は通常の `make test` から外し、`make test-l3` /
`make test-heavy` / `make test-fortran-heavy` / `make test-full` で明示実行します。

個別 target だけ確認する場合は次を使います。

```bash
FPM_ACTION=test ./build.sh --target test_version
```

KUDPC のログインノード上では、`make test*` / `fpm test` や同等の build/test を直接実行せず、
`tssrun` または `sbatch` で計算ノードに投入してください。

## 3. 実行フロー

通常は、`case.toml` から `beach.toml` を生成して実行します。preset 合成と高水準記法の詳細は
[beachx config / preset / 高水準記法ガイド](config_workflow.html) を参照してください。

1. `case.toml` を用意し、必要なら preset を追加・編集する
2. `beachx config render` で `beach.toml` を生成する
3. `beach` でシミュレーション実行
4. `output.dir` の出力ファイルを確認
5. Python CLI または `Beach` API で可視化

最終 `beach.toml` を手で管理する場合は、仕様を [Fortran パラメータファイル仕様](fortran_parameter_file.html) で確認してください。

### 3.1 `case.toml` から実行する最短例

```bash
mkdir run_periodic2
cd run_periodic2
beachx config init
beachx config render
beach beach.toml
```

### 3.2 `beach.toml` を直接使う場合

1. `beach.toml` を用意（仕様は [Fortran パラメータファイル仕様](fortran_parameter_file.html)）
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
- `mesh_potential.csv`（`write_mesh_potential = true`）
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv`（`history_stride > 0`）
- `potential_history.csv`（`write_potential_history = true` かつ `history_stride > 0`）
- `performance_profile.csv`（`BEACH_PROFILE=1`）
- `rng_state.txt`
- `macro_residuals.csv`

`mesh_triangles.csv` には要素ごとの `mesh_id` 列が含まれ、`mesh_sources.csv` で `mesh_id` と
元メッシュ設定（template kind / surface model / epsilon_r / 要素数）を対応付けます。
`conductor` は `field_bc_mode = "free"` の浮遊導体として等電位再配分され、
`dielectric` は現行ではメタデータのみなので、
含まれる場合は `summary.txt` に注意書きも出力します。
`mesh_potential.csv` を有効にすると、同じ要素順で centroid 電位 [V] も保存されます。

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
`sim.batch_count` は累積の到達バッチ数です。例えば既存 checkpoint が `batches=100` のとき `batch_count=150` で再開すると、追加で50バッチだけ実行します。`batch_count` が既存の処理済みバッチ数より小さい場合は停止します。

MPI再開時は `summary.txt` の `mpi_world_size` と現在の rank 数が一致している必要があります。

## 8. Python 後処理

### 8.1 CLI

```bash
beachx inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png \
  --potential-self-term area-equivalent

# sim.field_bc_mode = "periodic2" の mesh を周期セルへ寄せて描く
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --save-potential-mesh outputs/latest/potential_mesh_periodic.png \
  --apply-periodic2-mesh

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200

beachx slices outputs/latest \
  --grid-n 200 \
  --vmin -20 --vmax 20 \
  --save outputs/latest/potential_slices.png

beachx coulomb outputs/latest \
  --component z \
  --save outputs/latest/coulomb_force_z.png

beachx mobility outputs/latest \
  --density-kg-m3 2500 \
  --mu-static 0.4 \
  --save-csv outputs/latest/mobility_summary.csv
```

`beachx coulomb` は、近傍の `beach.toml` が見つかれば `mesh.templates` から object kind と順序を読み取り、既定では全 object を target 軸に並べて可視化します。特定 kind だけに絞る場合は `--target-kinds sphere` のように指定します。
`beachx mobility` は、既定で `plane` を support とみなし、それ以外の object を対象に合力・合トルクと `lift_ratio` / `slide_ratio` / `roll_ratio` を CSV 化します。質量由来の指標は `--density-kg-m3` と `beach.toml` の幾何情報が必要です。

旧 alias の `beach-inspect` / `beach-animate-history` / `beach-plot-coulomb-force-matrix` /
`beach-plot-potential-slices` / `beach-estimate-workload` / `beach-plot-performance-profile`
も当面は利用できますが、非推奨です。

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
beach.plot_mesh(apply_periodic2_mesh=True)
beach.plot_potential(apply_periodic2_mesh=True)
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

fig_force, ax_force = beach.plot_coulomb_force_matrix(
    component="z",
)
fig_force.savefig("outputs/latest/coulomb_force_z.png", dpi=150)

mobility = beach.analyze_coulomb_mobility(
    density_kg_m3=2500.0,
    mu_static=0.4,
)
for record in mobility.records:
    print(record.label, record.lift_ratio, record.slide_ratio)
```

## 9. MPI経路テスト（開発向け）

```bash
FPM_FC=mpiifort \
fpm test --target test_mpi_hybrid \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 2"
```

## 10. 実装挙動で誤解しやすい点

- 通常実行は `sim.batch_count` 分だけ進みます。再開実行では checkpoint の処理済みバッチ数から `sim.batch_count` に達するまで進みます。
- `sim.tol_rel` は監視値で、現行実装では早期終了条件に使いません。
- Fortran 本体の電場は要素重心点電荷近似 + `sim.softening` です。

camphor向けのMPIジョブ例は `examples/job_scripts/camphor_mpi_hybrid_job.sh` を参照してください。
