title: 開発・運用ワークフロー

Lang: [日本語](Workflow.md) | [English](Workflow.en.md)

# 開発・運用ワークフロー

BEACHのソースを編集し、局所実行からHPC検証まで進めるためのタスクガイドです。
Fortranがシミュレーションを実行し、Pythonが設定支援、後処理、可視化を担当します。
通常の利用だけが目的なら、[インストール](Installation.html)と[実行する](Execution.html)を参照してください。

## インストール済み環境でcaseを実行する

**前提:** `beach-bem`をインストールし、`beach`と`beachx`がPATHにあります。

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem
```

開発版を直接試す場合:

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

**操作:**

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

**期待する出力:** 実行後に`outputs/latest/summary.txt`、`charges.csv`、mesh情報が作成されます。

**解釈:** 正常終了は設定、build、実行経路の動作を示します。定常到達や物理妥当性は保証しません。

**次の選択:** caseの構成変更は[シミュレーションケースを設計する](ConfigurationRecipes.html)、
出力の意味は[出力ファイルを調べる](OutputGuide.html)、妥当性確認は[結果を検証する](ValidationGuide.html)へ進みます。

## 開発環境を作る

**前提:** Python、`make`、Fortran compiler、`fpm`を用意します。

```bash
make --version
gfortran --version
fpm --version
python --version
```

**操作:**

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
make check
```

**期待する出力:** Python packageが編集可能状態で入り、Fortran sourceがcompileされます。

**解釈:** `make check`は開発用の軽量build確認です。`BEACH_VERSION_MODE=dev`を使うため、
git hashだけが変わってもfpmの差分compileを再利用できます。

**次の選択:**

```bash
make run CONFIG=examples/beach.toml
make build VERSION_MODE=plain
make build VERSION_MODE=git
make install-generic
```

通常は`build.sh`経由の`make run`と`make check`を使います。低レベル確認だけに次を使います。

```bash
fpm run --target beach --profile release --flag "-fopenmp" -- examples/beach.toml
```

## 変更をテストする

**前提:** 変更した範囲と、必要なテストtierを決めます。複数の`fpm test`を並行実行しないでください。

**操作:**

```bash
make test-l0      # static/schema/build
make test         # L1: 通常の開発ループ
make test-l2      # L2: contract/integration
make test-l3      # L3: heavy FMMを含む累積検証
```

個別のFortran target:

```bash
FPM_ACTION=test ./build.sh --target test_version
```

**期待する出力:** 各commandが非zero statusなしで完了し、対象testがpassします。

**解釈:**

- L0: `git diff --check`、JSON schema、`make check`
- L1: L0 + Python tests + 軽量Fortran targets
- L2: L1 + C/kernel contractなどのintegration targets
- L3: L2 + heavy FMM targets

重い`test_dynamics_fmm`と`test_coulomb_fmm_core_basic`はL1に含まれません。

**次の選択:**

```bash
make test-heavy
make test-fortran-far-correction
make test-fortran-benchmark
make test-field-kernel-cache
make test-full
```

`test-field-kernel-cache`はnative cache/plane-oracle receipt用のopt-in gateで、L1/L2/L3には含まれません。
release判断は[物理リリースの検証](PhysicsReleaseVerification.html)に従います。

## OpenMPまたはMPIで実行する

**前提:** OpenMP build、またはMPI compilerで作成したMPI buildを用意します。

**操作:**

```bash
OMP_NUM_THREADS=8 beach beach.toml
mpirun -n 4 beach examples/beach.toml
```

粗いフェーズ計測を加える場合:

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach examples/beach.toml
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

**期待する出力:** 通常のsimulation出力に加え、profile有効時は`performance_profile.csv`が作成されます。

**解釈:** scaling比較には`simulation_total`行の`rank_max_s`を使います。MPI再開ではcheckpointの
`mpi_world_size`と現在のrank数を一致させます。

**次の選択:** hybrid経路だけを確認する場合:

```bash
FPM_FC=mpiifx \
fpm test --target test_mpi_hybrid \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 2"
```

Intel `ifx`でOpenMP buildを行う場合は`OPENMP_FLAG=-qopenmp`を使います。

## KUDPCで実行・テストする

**前提:** `hostname`、`module list`、利用可能なら`spartition`と`qgroup`でhostと割当を確認します。

**操作:** login nodeでは編集、軽いログ確認、`make check`、job投入・監視までに留めます。
`make test*`、`fpm test`、長時間実行、benchmarkは計算node上で実行します。

- 短い対話確認: `tssrun`
- batch実行: `sbatch` job内の`srun`

SysAで1 rankあたり112 OpenMP threadsを使う例:

```bash
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
srun beach beach.toml
```

**期待する出力:** Slurm allocation内でsimulationまたはtestが完了し、job logと通常の出力が残ります。

**解釈:** thread配置は性能比較と再現性のため固定します。login nodeでの正常動作を計算nodeの性能とみなしません。

**次の選択:** KUDPCの環境選択とjob構成は、repositoryのKUDPC pluginと
`examples/job_scripts/camphor_mpi_hybrid_job.sh`を参照します。

## 負荷を見積もってから実行する

**前提:** `boundary_inflow`、`plane_source`、deprecatedな`reservoir_face`、`photo_raycast`を使うcaseでは、
batchあたり粒子数が設定から決まります。

**操作:**

```bash
beachx workload examples/beach.toml --threads 8
```

MPI rankとcheckpoint残差を含める場合:

```bash
beachx workload examples/beach.toml \
  --threads 8 \
  --mpi-ranks 4 \
  --mpi-rank 0 \
  --macro-residuals outputs/latest/macro_residuals.csv
```

**期待する出力:** rank局所の`total_particles`と全rankの`global_total_particles`が表示されます。

**解釈:** これはworkloadの事前見積もりであり、walltime予測や物理妥当性の判定ではありません。

**次の選択:** 見積もりが大きい場合はmesh、macro粒子数、`batch_count`を小さくしたsmoke caseから始めます。

## checkpointから再開する

**前提:** 読み込み元に`summary.txt`、`charges.csv`、RNG状態などのcheckpoint一式があります。

**操作:**

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../parent_run/outputs/latest"
```

**期待する出力:** checkpointを読み込み、新しい出力を`outputs/continuation`へ保存します。

**解釈:** `sim.batch_count`は累積到達batch数です。checkpointが`batches=100`で
`batch_count=150`なら50 batchを追加します。既存の処理済みbatch数より小さい値は受理されません。

**次の選択:** 同じdirectoryへ続けて書く場合は`restart_from`を省略し、`dir`をcheckpoint directoryにします。

## 変更後に確認する実装契約

- 通常実行は`sim.batch_count`で指定したbatch数まで進みます。
- `sim.tol_rel`は監視・出力値であり、早期停止条件ではありません。
- box geometryと周期軸は`[domain]`、場closureは`[field_boundary]`で指定します。
- global粒子面は`[particle_boundary]`、species overrideは`[particles.species.boundary]`で指定します。
- 外部reservoirの流入modelと`phi_infty`は`[reservoir]`に置き、流入面は`[particles.species.boundary_inflow]`で選びます。
- v1.0の標準surface modelはinsulator accumulationです。
- 境界reservoir + closed PEは有限box内のreduced closureであり、外部領域を自己無撞着に解くmodelではありません。
- 実行成功、数値収束、物理妥当性は別々に確認します。

公開API、設定、出力を変更した場合は、日英ドキュメント、examples、schema、対応testも同じ変更で更新してください。
