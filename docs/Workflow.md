title: 開発ワークフロー

Lang: [日本語](Workflow.md) | [English](Workflow.en.md)

# 開発ワークフロー

BEACH のソースを変更し、変更範囲に応じたテストを選ぶための contributor guide です。
通常利用のインストール、case 実行、OpenMP / MPI 実行、負荷見積もり、再開は
[インストール](Installation.html)と[実行する](Execution.html)を参照してください。
コードを初めて読む場合は、先に[開発者向けアーキテクチャ](Architecture.html)で実行フローと
state の所有者を確認してください。

## 開発環境を作る

**前提:** Python、`make`、Fortran compiler、`fpm` を用意します。

```bash
make --version
gfortran --version
fpm --version
python --version
```

**操作:** checkout の root で editable install と軽量 build check を実行します。

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
make check
```

**期待する結果:** Python package が checkout を参照し、Fortran source が compile されます。
`make check` は `BEACH_VERSION_MODE=dev` を使うため、Git hash だけが変わった場合も fpm の
差分 compile を再利用できます。

通常は `build.sh` 経由の Make target を使います。個別の Fortran test target は次の形で実行します。

```bash
FPM_ACTION=test ./build.sh --target test_particle_stepper
```

複数の `fpm test` または `build.sh` test target を並行実行しないでください。同じ `build/` directory を
共有するため、compile artifact が競合します。

## 変更からテストを選ぶ

最初に直接関係する test を実行し、その後で少なくとも表の gate まで確認します。複数 subsystem にまたがる変更では、
該当する行を組み合わせます。公開契約や数値結果へ影響する変更では、直接 test だけで完了としません。

| 変更範囲 | 直接確認する test | 最低 gate / 追加確認 |
| --- | --- | --- |
| documentation、navigation、日英対応 | `pytest -q tests/python/test_docs_sync.py tests/python/test_documentation_contracts.py` | `make test-l1`。site 構成を変えた場合は `python tools/sync_starlight_docs.py` の後に `npm --prefix docs-site run check` |
| TOML、schema、parser、既定値 | `test_app_config_parser`、`test_physics_config_types`、`tests/python/test_config_schema.py`、`tests/python/test_config_cli.py` | `make test-l1` と `make schema-check` |
| mesh template、OBJ import、panel geometry | `test_templates_importers_runtime`、`test_panel_geometry_near`、`test_panel_kernel` | `make test-l1`。panel FMM に影響する場合は `make test-l3` |
| 粒子 source、reservoir、光電子注入 | `test_injection_sampling`、`test_reservoir_injection`、`test_external_field_velocity_grid` | `make test-l1` |
| Boris、collision、box boundary、particle event | `test_particle_stepper`、`test_boundary`、`test_dynamics_basic` | `make test-l1`。MPI / OpenMP 経路も変える場合は `make test-mpi` |
| surface model、charge closure、charge ledger | `test_surface_models`、`test_surface_current_model`、`test_charge_ledger`、`test_simulator` | `make test-l1`。matching-plane は `test_matching_plane_simulator` も確認 |
| field snapshot、Direct / Treecode | `test_electrostatic_snapshot`、`test_dynamics_field_solver`、`test_panel_kernel` | `make test-l1` |
| FMM、periodic2、zero / nonzero mode | 関連する `test_coulomb_fmm_*`、`test_periodic_*`、`test_dynamics_fmm` | `make test-l3`。遠方補正は `make test-fortran-far-correction` |
| C ABI / native field kernel | `test_field_kernel_c`、`test_periodic_zero_mode_c` | `make test-l2`。cache receipt は `make test-field-kernel-cache` |
| output、checkpoint、fingerprint、restart | `test_output_writer_io`、`test_output_writer_potential`、`test_restart`、`test_model_fingerprint` | `make test-l1` と関連する Python reader test |
| Python reader、解析、CLI | 変更 package に対応する `tests/python/test_*.py` | `make test-l1` |

Fortran test 名は `fpm.toml`、tier の所属は `Makefile` が実行上の正本です。subsystem、source、直接 test、
関連 document の対応は[開発者向けアーキテクチャ](Architecture.html)にあります。

## Test tier を使い分ける

```bash
make test-l0      # static / schema / build
make test         # L1: Python + 軽量 Fortran（test-l1 の alias）
make test-l2      # L1 + C / kernel contract
make test-l3      # L2 + heavy FMM / panel
```

- L0 は `git diff --check`、source text、JSON schema、`make check` を確認します。
- L1 は L0 に全 Python test と通常の Fortran test target を加えます。
- L2 は C ABI と periodic zero-mode C contract を加えます。
- L3 は `test_dynamics_fmm`、FMM core、panel near-correction などの heavy target を加えます。

次の gate は tier へ常時含めず、変更内容または release 判断に応じて明示的に実行します。

```bash
make test-heavy
make test-fortran-far-correction
make test-field-kernel-cache
make test-mpi
make test-fortran-benchmark
make test-physics-release
make test-full
```

`test-field-kernel-cache` は native cache / plane-oracle receipt 用です。性能比較は correctness test と分け、
release profile の `make test-fortran-benchmark` を使います。物理 release の判定と収束 fixture は
[物理リリースの検証](PhysicsReleaseVerification.html)に従います。

handoff 前に、Fortran を変更した場合は `make fmt-check-fortran`、Python を変更した場合は `ruff check .` を実行し、
最後に `git diff --check` で whitespace error を確認します。これらの format / lint check が必要な変更では、
test tier の通過だけで代用しません。

## KUDPC で開発テストを実行する

KUDPC では test payload を実行する前に `hostname`、`module list`、利用可能なら `spartition` と `qgroup` を確認し、
host role、active `Sys*` module、割当を確定します。

login node（`camphor*`、`laurel*`、`cinnamon*`、`gardenia*`）では、編集、差分確認、短い log 読み取り、
`make check`、job 投入・監視までに留めます。`make test*`、`pytest`、`fpm test`、MPI / OpenMP test、
benchmark は直接実行しません。

- 短い開発 test は `tssrun` で計算 node の割当を取得して実行する。
- tier test、MPI test、heavy / release gate は `sbatch` job 内の `srun` で実行する。
- job 内でも複数の `fpm test` を並行実行しない。
- checkout、入力、cache、log は計算 node から見える `/home`、`/LARGE0`、`/LARGE1`、または合意した
  `/FAST` path に置く。login node の `/TMP` を前提にしない。

production simulation の job 設計、thread / rank 構成、負荷見積もりは[実行する](Execution.html)の責務です。

## 公開契約を変更するとき

公開契約を変えた場合は、実装だけでなく対応する正本と consumer を同じ変更で更新します。

| 変更した契約 | 同時に確認・更新するもの |
| --- | --- |
| Fortran の simulation behavior | `SPEC.md`、対応する model / numerical-method page、回帰 test、必要なら収束 test |
| TOML table、key、型、既定値、制約 | `schemas/beach.schema.json` と配布 copy、Fortran parser、Python validator、`Parameters.md` / `.en.md`、examples |
| CLI または Python / Fortran API | help、公開 signature、examples、日英 API document、consumer test |
| output file、column、生成条件 | `schemas/beach.output-manifest.json`、Fortran writer、Python reader、`OutputGuide.md` / `.en.md`、restart compatibility |
| checkpoint state または fingerprint | schema version / compatibility rule、writer、loader、restart test、`SPEC.md` の再開契約 |
| 数値手法または物理 model | 適用範囲、fail-closed 条件、Direct / oracle 比較、数値収束、`ValidationGuide.md` / `.en.md` |

すべての文書変更で、日本語版と英語版の command、identifier、制約、警告、期待する結果を対応させます。
`sim.tol_rel` は現行実装では監視・出力値であり、早期停止条件として説明しません。実行成功、test 通過、
数値収束、物理妥当性は別々に報告してください。

正本の責務分担は[開発者向けアーキテクチャ](Architecture.html#正本の責務を区別する)、user case の妥当性確認は
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。
