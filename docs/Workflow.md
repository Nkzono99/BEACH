title: 開発・運用ワークフロー

Lang: [日本語](Workflow.md) | [English](Workflow.en.md)

# 開発・運用ワークフロー

このプロジェクトでは、Fortranがシミュレーションを実行し、Pythonが後処理と可視化を担当します。
通常は、`pip install beach-bem`で導入した`beach`コマンドを使って実行します。

通常のcaseの実行方法は[実行する](Execution.html)にまとめています。ソースから開発する場合のsetup、test、
MPI/OpenMP、HPC運用は以下で説明します。

## 1. 利用者向けセットアップ（推奨）

### 1.1 ツール確認

```bash
make --version
gfortran --version
fpm --version
python --version
```

### 1.2 PyPI から一括インストール（推奨）

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem
```

`pip install`は`make install`を実行し、Python CLIとFortran実行バイナリを同時に導入します。
pip経由のbuildは既定で`INSTALL_PROFILE=auto`を使い、失敗した場合は`generic`にfallbackします。
フォールバックを無効化する場合は `BEACH_PIP_FALLBACK_GENERIC=0` を指定します。

```bash
export PATH="$HOME/.local/bin:$PATH"
```

開発版を直接試す場合は Git URL から導入できます。

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
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

開発中の標準的な確認には`make check`を使います。`BEACH_VERSION_MODE=dev`により、
Fortranに渡すversion文字列は`1.5.0-dev`のように固定されます。git hashが変わってもfpmのcompile-flag hashは変わらず、
差分コンパイルを再利用できます。

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
make test-l3      # L3: cumulative L0-L3 verification
make test-physics-release  # HPC: minimal release correctness + MPI manifest
make test-heavy   # heavy Fortran targets only
make test-fortran-far-correction  # oracle far-correction correctness
make test-fortran-far-correction-diagnostics  # assertion-free diagnostics
make test-fortran-benchmark  # release-profile runtime benchmark
make test-field-kernel-cache  # opt-in native cache/plane-oracle receipt gate
make test-full    # unfiltered fpm test
```

BEACH のテストは開発ループ向けに階層化しています。

- L0: `git diff --check`、JSON schema parse check、`make check`
- L1: L0 + Python tests + 軽量 Fortran test targets（`make test` / `make test-l1`）
- L2: L1 + contract/integration targets（C field-kernel contract など）
- L3: L2 + heavy FMM targets（release gate / nightly / main 統合前）

`make test-fortran`は、軽量なFortran targetを実行するaliasです。重いFMM系の
`test_dynamics_fmm`と`test_coulomb_fmm_core_basic`は通常の`make test`には含まれません。
`make test-l3` / `make test-heavy` / `make test-fortran-heavy` / `make test-full`で明示的に実行します。

`m2l_root_oracle`のcorrectnessは計算量が大きいため、`make test-fortran-far-correction`で実行します。
数値assertを持たない`test_periodic2_flat_oracle_diag`は、`make test-fortran-far-correction-diagnostics`で実行します。
速度比較には、release profileの`make test-fortran-benchmark`を使います。
実際のメッシュで`point`と`triangle_p0`の場評価コストを比較する場合は、比較する各設定に対して
次を実行します。`FIELD_KERNEL_BENCHMARK_TARGET_COUNT`は領域内と面近傍でそれぞれ
評価する点数です。出力はメッシュ構築、
ソルバ初期化、電荷更新、電場・電位評価を分けたCSV形式です。

```bash
make benchmark-field-kernel CONFIG=beach.toml FIELD_KERNEL_BENCHMARK_TARGET_COUNT=2048
```

Intel `ifx`では`OPENMP_FLAG=-qopenmp`も指定します。

`make test-field-kernel-cache`はshared kernelをbuildし、その絶対pathをnative periodic plane-oracle receipt testに渡します。
長時間実行するopt-in gateであり、L1/L2/L3と`make test-physics-release`には含まれません。
Intel `ifx` / `mpiifx` の tiered test は、既知の配列一時オブジェクトごとに巨大な stack trace を
出す `arg_temp_created` check だけを既定で抑制します。bounds など他の debug check は維持されます。
一時配列診断そのものを調べる場合は、例えば
`FORTRAN_TEST_FLAGS="-qopenmp" make test-fortran-heavy FPM_FC=mpiifx` のように明示上書きします。

個別 target だけ確認する場合は次を使います。

```bash
FPM_ACTION=test ./build.sh --target test_version
```

KUDPC のログインノード上では、`make test*` / `fpm test` や同等の build/test を直接実行せず、
`tssrun` または `sbatch` で計算ノードに投入してください。

## 3. 実行フロー

通常は、`beach.toml`を編集し、そのまま`beach`に渡します。作成・検査方法は
[`beach.toml`を作成・検証する](Configuration.html)にまとめています。

1. `beach.toml` を用意する
2. `beachx lint beach.toml` で設定を確認する
3. `beach beach.toml` でシミュレーション実行
4. `output.dir` の出力ファイルを確認
5. Python CLI または `Beach` API で可視化

`beach.toml` の仕様は [Fortran パラメータファイル仕様](Parameters.html) で確認してください。

### 3.1 最短例

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
```

### 3.2 `beach.toml` を直接使う場合

1. `beach.toml` を用意（仕様は [Fortran パラメータファイル仕様](Parameters.html)）
2. 必要ならbox基準の座標・配置パラメータを使う（Fortran parserが実座標へ変換）
3. `beach beach.toml` でシミュレーション実行
4. `output.dir` の出力ファイルを確認
5. Python CLI または `Beach` API で可視化

## 4. 実行コマンド

### 4.1 推奨: `beach`

```bash
beach beach.toml
```

引数なし実行では、カレントディレクトリの `beach.toml` を自動読込します。

### 4.2 スレッド数指定

```bash
OMP_NUM_THREADS=8 beach beach.toml
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

`mesh_triangles.csv`には、要素ごとの`mesh_id`列があります。`mesh_sources.csv`は、各`mesh_id`を
元のmesh設定（template kind / surface model / epsilon_r / 要素数）に対応付けます。
`conductor`は`field_bc_mode = "free"`の浮遊導体として扱われ、等電位になるように電荷が再配分されます。
`dielectric`は現行実装ではメタデータのみであり、含まれている場合は`summary.txt`に注意書きを出力します。
`mesh_potential.csv` を有効にすると、同じ要素順で centroid 電位 [V] も保存されます。

MPI実行（`world_size > 1`）では乱数状態だけが rank 別です。reservoir の端数は全 rank で共有するため、
root rank が単一の `macro_residuals.csv` を保存します。

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals.csv`

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
  --macro-residuals outputs/latest/macro_residuals.csv
```

`total_particles`は選択したrankの見積もり、`global_total_particles`は全rankの合計です。reservoirについては、
`local_reservoir_particles`がrankへの分配後、`global_reservoir_particles`が分配前の粒子数を示します。

## 7. 再開実行（resume）

```toml
[output]
dir = "outputs/latest"
resume = true
```

同じ`output.dir`で`beach`を再実行すると、`summary.txt` / `charges.csv` / RNG状態を読み込み、続きから計算します。
読み込み元 checkpoint と新しい出力先を分ける場合は、`restart_from` を使います。

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../parent_run/outputs/latest"
```

この場合、checkpointは`restart_from`から読み込みます。新しい`summary.txt`、`charges.csv`、履歴、RNG状態は`dir`に書き出します。

`sim.batch_count`は、累積の到達batch数です。例えば、既存のcheckpointが`batches=100`のときに
`batch_count=150`で再開すると、追加で50 batchを実行します。`batch_count`が既存の処理済みbatch数より小さい場合は停止します。

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

# 周期メッシュを n 周期分複製して描く（1 なら 3x3 = 9 コピー）
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_tiled.png \
  --periodic2-repeat 1 \
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

# Fortran FMM core と同じ field kernel で object ごとの総電荷・合力・合トルクを出す
make build-kernel
beachx kernel-forces outputs/latest \
  --save-csv outputs/latest/object_forces_kernel.csv

# central-cell primary 自己場だけを除き、周期画像を保持した離脱経路を出す
beachx object-detachment outputs/latest \
  --config beach.toml \
  --target-mesh-id 6 \
  --periodic-model infinite-physical \
  --z-max-m 2.0e-4 \
  --z-points 65 \
  --mass-kg 2.0e-12 \
  --gravity-m-s2 9.80665 \
  --output-dir outputs/latest/object_detachment
```

`beachx coulomb`は、近傍の`beach.toml`が見つかると、`mesh.templates`からobject kindと順序を読み取ります。
既定では全objectをtarget軸に並べて可視化します。特定のkindに絞る場合は、`--target-kinds sphere`のように指定します。

`beachx mobility`は、既定で`plane`をsupportとみなし、それ以外のobjectを解析対象にします。
合力、合トルク、`lift_ratio` / `slide_ratio` / `roll_ratio`をCSVに出力します。質量由来の指標には、
`--density-kg-m3`と`beach.toml`の幾何情報が必要です。

`beachx kernel-forces`は、`libbeach_field_kernel`を介してFortran FMM coreをPythonから呼び出します。
`beach.toml`の`sim.softening` / `sim.field_bc_mode` / periodic2 / tree設定を使い、objectごとのnet forceを計算します。
共有libraryは`make build-kernel`で`build/libbeach_field_kernel.so`に生成できます。別の場所に置く場合は`--library`または
`BEACH_FIELD_KERNEL_LIB`を指定します。設定ファイルが出力directoryの近傍にない場合は、`--config path/to/beach.toml`を指定します。

`kernel-forces`は、target objectの周期画像も除外する旧`exclude_target_lattice`診断です。
`object-detachment`は、central primaryだけを除外する`exclude_primary_keep_images`を使います。
凍結sourceに対する瞬時wrench、鉛直経路、仕事、重力と有限range付着を含むfrom-rest barrierを4つの成果物に出力します。
`configured`はrunの設定をそのまま使い、`infinite-physical`はx/y periodicのcached `k != 0`と`E_bottom=0` zero modeを使います。
CLI 既定重力は月面の `1.62 m/s^2` ですが、上例は地上の `9.80665 m/s^2` を明示しています。
cached operatorの初回生成や、多数のpath点の評価には時間がかかる場合があります。KUDPCではlogin nodeで実行せず、
計算nodeに投入してください。cold cacheの生成専用ジョブは、SysAの`p=1:t=112:c=112`を基準とします。

既存のsimulation allocationがrankを持つ場合、そのrankはcache生成にも使われます。2026-07-12の旧レゴリス入力では、
`1x112=47.0 s`、`2x112=36.7 s`、`4x112=31.5 s`、`6x112=30.3 s`でした。cold buildのために4--6 rankまで増やしても、
coreの利用効率は低い結果です。同じfingerprintのwarm runではcacheを再生成しません。

CLIが正常終了すれば、成果物の作成は完了しています。物理的な妥当性は、`path.status`、仕事/電位差、quadrature、shell/cache、
経路上端への感度を使って別途確認します。非中性periodic cellで得た有限高さのspeedは、無限遠へのescape speedではありません。

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
FPM_FC=mpiifx \
fpm test --target test_mpi_hybrid \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 2"
```

KUDPCのIntel環境では、release MPI buildの既定compilerとして`mpiifx`を使います。
`MPI_FC=mpiifort`を指定すればclassic compilerも使えますが、同一の3000 batch fixtureで大幅な速度低下が観測されています。
そのため、互換性の確認以外のproduction runでは`mpiifx`を使います。

SysAで1 rankあたり多数のOpenMP threadを使うproduction runでは、thread配置を再現できるように次の値を指定します。

```bash
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
srun beach beach.toml
```

この指定は、性能比較時のthread配置を固定するためのものです。実測した300 batch fixtureでは、bindの有無による差は0.4%未満でした。

## 10. 実装挙動の確認項目

- 通常実行は `sim.batch_count` 分だけ進みます。再開実行では checkpoint の処理済みバッチ数から `sim.batch_count` に達するまで進みます。
- `sim.tol_rel` は監視値で、現行実装では早期終了条件に使いません。
- Fortran 本体の要素核は `field.element_kernel` で選びます。互換既定の `point` は要素重心点電荷 +
  `sim.softening`、`triangle_p0` は要素総電荷を三角形上の一定面密度として積分し、`softening=0` を必須とします。

camphor向けのMPIジョブ例は `examples/job_scripts/camphor_mpi_hybrid_job.sh` にあります。
`test-physics-release`は、収束出力に必要なL1 subset、L3 heavy、far-correction correctness、MPI電荷収支、
MPI periodic-cache gateを順に実行します。portable CIで実行済みのL2全体は繰り返しません。

既定の`build/physics-release/manifest.txt`には、commit、dirty state、host、compiler、各gateのstatus、経過時間、最大RSSを保存します。
同じdirectoryの`convergence.csv`には、mesh、dt、FMM order、outer gridなどの収束値を保存します。
`test_l3-target-timings.csv`と`far_correction-target-timings.csv`には、Fortran targetごとのprofile、status、経過秒を保存します。

KUDPCのlogin nodeではこのtestを実行できません。Slurm allocation内では、MPI payloadのrunnerに`srun`を使います。
manifestの出力先は`PHYSICS_RELEASE_MANIFEST=/path/to/manifest.txt`で変更できます。
各検証targetと収束基準は[Physics release verification](PhysicsReleaseVerification.md)にまとめています。
