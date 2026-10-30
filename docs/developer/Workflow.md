title: 開発ワークフロー

# 開発ワークフロー

BEACH のソースを変更し、変更範囲に応じたテストを選ぶための contributor guide です。
通常利用のインストール、case 実行、OpenMP / MPI 実行、負荷見積もり、再開は
[インストール](../start/Installation.md)と[実行する](../guide/Execution.md)を参照してください。
コードを初めて読む場合は、先に[開発者向けアーキテクチャ](Architecture.md)で実行フローと
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
python -m pip install -e '.[dev]' --no-build-isolation
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

### ソース配置変更後のビルドキャッシュ

`build.sh` とそれを使う Make target は、ソースの追加・削除・移動やコンパイラの差し替えを検出すると、
対象コンパイラの fpm キャッシュを削除して再構築します。初回のビルドも同じ処理を行い、古い module や archive の混入を防ぎます。
コンパイラの実体パス・version と、`src/`、`app/`、`tests/fortran/`、`benchmarks/fortran/` にある
コンパイル対象の相対パス一覧を `build/.beach-cache-<compiler>.context` に記録します。

削除対象は `build/<compiler>_<16-hex-digits>/` 内の object、module、archive、実行ファイルです。
他のコンパイラのキャッシュ、依存ライブラリの checkout、log、生成ドキュメントは残ります。
ソース本文の編集、Git commit、profile・flag・target の変更ではこの一括削除を行わず、fpm の通常の増分ビルドと
キャッシュ選択に任せます。直接 `fpm` を実行した場合は、この自動検出の対象外です。

`make build-kernel` は同じコンパイラ・profile・flag で `build.sh --list` が示す archive を選び、共有ライブラリへリンクします。

fpm は、依存先が同じ実行内で再コンパイルされない限り、ソースが変わっていない対象を再コンパイルしません。そのため途中で失敗・中断したビルドの後は、更新済み module に依存する対象が古い型レイアウトのまま残ることがあります。`build.sh` は実行中に `build/.fpm-run-unfinished` を置き、前回分が残っていればそれ以降に書かれた `*.digest` を消してから fpm を呼ぶので、次のビルドで該当対象と依存先が再コンパイルされます。直接 `fpm` を実行した場合はこの保護の対象外です。

古い成果物の混入が疑われる場合は、一度だけ再構築してから通常のテストを実行します。

```bash
BEACH_REBUILD=1 make check
make test-l2
```

再構築時は `[build.sh] Resetting ... caches:` が表示されます。`BEACH_REBUILD=1` は最初のコマンドだけに指定します。
`BEACH_REBUILD=1 make test-l2` のように複数 target のテストへ渡すと、各 target で再構築されます。

## 変更からテストを選ぶ

最初に直接関係する test を実行し、その後で少なくとも表の gate まで確認します。複数 subsystem にまたがる変更では、
該当する行を組み合わせます。公開契約や数値結果へ影響する変更では、直接 test だけで完了としません。

| 変更範囲 | 直接確認する test | 最低 gate / 追加確認 |
| --- | --- | --- |
| documentation、navigation、日英対応 | `pytest -q tests/python/repo/test_docs_sync.py tests/python/repo/test_documentation_contracts.py` | `make test-l1`。site 構成を変えた場合は `python tools/sync_starlight_docs.py` の後に `npm --prefix docs-site run check` |
| TOML、schema、parser、既定値 | `test_app_config_parser`、`test_physics_config_types`、`tests/python/config/`、`tests/python/cli/test_config_commands.py`、`make test-config-contract` | `make test-l2` と `make schema-check` |
| mesh template、OBJ import、panel geometry | `test_templates_importers_runtime`、`test_panel_geometry_near`、`test_panel_kernel` | `make test-l1`。panel FMM に影響する場合は `make test-l2` |
| 粒子 source、reservoir、光電子注入、batch 生成 | `test_injection_sampling`、`test_reservoir_injection`、`test_particle_pool` | `make test-l1` |
| Boris、collision、box boundary、particle event | `test_particle_stepper`、`test_boundary`、`test_dynamics_basic`、`test_step_growth` | `make test-l1`。MPI / OpenMP 経路も変える場合は `make test-mpi` |
| surface model、charge closure、charge ledger | `test_surface_models`、`test_surface_current_model`、`test_charge_ledger`、`test_surface_charge_closure` | `make test-l1`。周期導体は `test_surface_models_periodic`、外部根の更新は `test_zhao_outflow_refresh` も確認 |
| batch loop、統計、停止時の診断文、soft discard | `test_simulator_loop`、`test_simulator_failure_context`、`test_soft_discard` | `make test-l1` |
| field snapshot、Direct / Treecode | `test_electrostatic_snapshot`、`test_dynamics_field_solver`、`test_panel_kernel` | `make test-l1` |
| FMM、periodic2、zero / nonzero mode | 関連する `test_coulomb_fmm_*`、`test_periodic_*`、`test_dynamics_fmm` | `make test-l2`。遠方補正は `make test-fortran-far-correction` |
| C ABI / native field kernel | `test_field_kernel_c`、`test_periodic_zero_mode_c`、`make test-python-native` | `make test-l2`。cache と周期解析解は `make test-field-kernel-cache` |
| output、checkpoint、fingerprint、restart | `test_output_writer_io`、`test_restart` | `make test-l1` と関連する Python reader test |
| Python reader、解析、CLI | 変更 package に対応する `tests/python/<package>/`（`io`、`analysis`、`native`、`plot`、`config`、`cli`、`tooling` など） | `make test-l1` |

挙動を変えないリファクタでは、計算ノードで変更前と変更後の checkout それぞれで `tools/golden_regression.py run OUT_DIR` を実行し、
`tools/golden_regression.py compare BASE_DIR NEW_DIR` で出力を比べます。examples を縮小した3ケース（free-space direct、
periodic2 FMM と光電子、零電流の外部根更新）を release・OpenMP 1 thread・固定 seed で走らせ、`summary.txt`、`charges.csv`、
各 history を、build 情報や時間の行を除いたテキストの完全一致で比較します。演算順が変わる変更では `--rtol 1e-14` のように
許容を指定すると数値を相対差で比べます。どちらでも片側にしかない case・ファイル、有限値から NaN/Inf への変化
（非有限値は同じ表記どうしだけ一致）は失敗です。
`--rtol` は有限の非負値に限り、NaN・Inf・負値は比較を開始せず終了コード 2 で拒否します。

Fortran test 名は `fpm.toml`、tier の所属は `Makefile` が実行上の正本です。subsystem、source、直接 test、
関連 document の対応は[開発者向けアーキテクチャ](Architecture.md)にあります。

## Test tier を使い分ける

```bash
make test-l0      # static / schema / build
make test         # L1: Python + 軽量 Fortran（test-l1 の alias）
make test-l2      # L1 + C ABI / free-space FMM / Python native kernel + Python / Fortran 設定契約
```

- L0 は `git diff --check`、source text、`make python-lint`（`ruff check .` と段階導入中の format check）、`mypy`、JSON schema、`make check` を確認します。lint と型検査の設定は `pyproject.toml` にあります。
- L1 は L0 に通常の Python test、検証 workflow の代表 3 ケース、軽量 Fortran target を加えます。
- L2 は C ABI、periodic zero-mode、free-space FMM の精度比較（`test_dynamics_fmm`、`test_dynamics_panel_fmm`、`test_coulomb_fmm_core_basic`、`test_coulomb_fmm_core_panel`）、共有 kernel を build した Python/native 数値比較、Fortran の設定読込と lint CLI を加えます。
- L1 の Fortran target は 1 本数秒に収めます。数秒を超える test は L2 へ、`cached_kneq0` operator を作る test は `make test-fortran-far-correction` へ置きます。

Python の gate は file ではなく pytest marker で選びます。native の成果物を使う test は fixture で受け取り、
`tests/python/conftest.py` が fixture から marker を付けます。

| marker | 付く test | 実行する gate |
| --- | --- | --- |
| `native` | `field_kernel_library` fixture（共有 kernel）を使う | `make test-python-native` |
| `config_native` | `native_checker` fixture（`beach --check-config`）を使う | `make test-config-contract` |
| `field_kernel_cache` | 実際の periodic cache を生成する（明示的に付ける） | `make test-field-kernel-cache` |
| `workflow` / `smoke` | 検証 workflow の全ケース / 代表 3 ケース | `make test-validation-workflow` / `make test-python` |

`make test-python` は上の 3 marker と `smoke` でない `workflow` を除いて実行するので、`build/` に古い kernel が残っていても native test を選びません。
`make test-validation-workflow` は stage、改変検知、投入、receipt、解析の全ケースを実行します。
CI は L2 とこの専用 gate を実行します。`pytest -q` はすべてを収集します。共有 kernel がない native test は skip し、
`field_kernel_cache` は `BEACH_RUN_FIELD_KERNEL_CACHE_TESTS=1` のときだけ実行します。
実物の build-origin 比較では、同じ clean commit から作った実行ファイルと共有 kernel を
`BEACH_VALIDATION_BINARY`、`BEACH_FIELD_KERNEL_LIB` で指定します。

次の gate は tier へ常時含めず、変更内容または release 判断に応じて明示的に実行します。

```bash
make test-fortran-far-correction
make test-field-kernel-cache
make test-validation-workflow
make test-mpi
make test-fortran-benchmark
make test-physics-release
make test-full
```

`test-field-kernel-cache` は cache 再利用、有限 shell と無限周期の比較、一様・cosine 平面の解析解、oracle receipt を検証します。
mesh 解像度が異なっても同じ box-cell 演算子は再利用できます。各評価の fingerprint、path、SHA256 を照合します。
`BEACH_FIELD_KERNEL_LIB` で指定した共有 kernel が存在しない場合、native test は失敗します。
性能比較は release profile の `make test-fortran-benchmark` で測定し、L1 の数値一致から分けます。物理 release の判定と収束 fixture は
[物理リリースの検証](PhysicsReleaseVerification.md)に従います。

handoff 前に、Fortran を変更した場合は `make fmt-check-fortran`、Python を変更した場合は `make python-lint` を実行し、
最後に `git diff --check` で whitespace error を確認します。これらの format / lint check が必要な変更では、
test tier の通過だけで代用しません。

### Python の format と import 整列

Ruff の formatter と import 整列（rule `I`）は、大きな整形差分が他の変更と混ざらないよう path 単位で段階導入しています。
対象 path は `Makefile` の `PYTHON_FORMAT_PATHS` と、`pyproject.toml` の `[tool.ruff.lint.per-file-ignores]` にある
否定パターン（対象外の path で `I` を無視する）の 2 か所に同じ並びで書きます。

```bash
make python-format        # 対象 path の import を整列し、ruff format を適用する
make python-format-check  # 対象 path の format を検査する（make python-lint と L0 が実行する）
```

対象外の path を巻き込まないよう、段階導入が終わるまで `ruff format .` は使いません。

| stage | 対象 | 状態 |
| --- | --- | --- |
| 1 | `beach/native`、`beach/plot`、`beach/tooling`、`beach/facade.py` | 適用済み |
| 2 | `beach/io`、`beach/config`、`beach/cli`、`beach/analysis`、`beach/__init__.py`（package 全体） | 未適用 |
| 3 | `tests`、`tools`、`examples`、`docs` の Python | 未適用 |

次の stage は、対象 path を触る作業 branch が少ない時期に、次の順で進めます。

1. 上の 2 か所に path を追加し、`make python-format` を実行する。
2. 整形だけの commit を作り、コード変更とは別にする。
3. stage 3 を終えたら path の列挙をやめ、`select` の `I` を全体に効かせて per-file-ignores を削除し、
   `python-format-check` を `ruff format --check .` にする。

format の結果は Ruff の minor version ごとの style に依存するため、`pyproject.toml` の dev 依存で Ruff を minor の範囲に固定しています。
範囲を上げる場合は、同じ変更で対象 path を整形し直します。

formatter の導入後に追加する lint rule の候補は次のとおりです。rule ごとに修正と設定を 1 つの変更にまとめます。

- `UP035`（`typing` からの deprecated import）と `UP037`（不要な引用符付き annotation）。自動修正で直せます。
- `B905`（`strict=` のない `zip`）。`strict=True` は長さの不一致で例外を出すので、呼び出しごとに長さが揃う前提を確認して指定します。

`make python-type-check` は、`beach.io.summary` と設定の共有エラー・schema 読み込みを `mypy --strict` 相当で検査します。
解析・描画などの未移行部分は検査対象に含めず、型注釈を整えた module から `tool.mypy.files` に追加します。
結果の読み込み・解析・描画の回帰テストは、それぞれ `tests/python/io/`、`tests/python/analysis/`、`tests/python/plot/` にあります。

## 設定の契約を確認する

通常運用では実行前に `beachx lint` を通します。lint は同梱 schema の検査の後に `beach --check-config` を呼ぶので、
利用者が Fortran の preflight を別に実行する必要はありません。

構造の正本は `schemas/beach.schema.json` です。編集後は配布 copy を同期し、設定検証の共通ケースを確認します。

```bash
make schema-sync
make schema-check
make test-config-contract
```

`schema-sync` は Python package と BEACH context plugin の 2 copy を正本から更新します。
`schema-check` は JSON の構文、3 ファイルの byte 単位の一致、公開 16 table の schema と Fortran
読取実装のキー集合の一致を確認します。キーの照合は宣言上の検査であり、値域や組合せは実行 test で確認します。
`test-config-contract` は現在の Fortran `beach` を build し、代表ケースごとに `beachx lint` が Fortran の preflight の
採否をそのまま返すことを確認します。preflight だけで simulation output を作った場合も失敗します。
この gate は L2 に含まれます。通常の Python test だけでは Fortran の読み込みを確認したことにはなりません。

検証の責務は 2 層です。キー名・型・列挙値・値域は `schemas/beach.schema.json` が持ち、Python の lint が検査します。
組合せ制約と派生値は Fortran の領域別 preflight (`_preflight_sim` / `_preflight_particles` / `_preflight_surface` /
`_finalize`) が持ちます。同じ制約を Python に書きません。新しい組合せ制約は Fortran の preflight に置き、
`tests/python/config/test_config_contract.py` の代表ケースと、必要なら `test_app_config_parser` の拒否ケースを追加します。
不正な組合せの網羅は保証しません。研究での発生頻度、結果への影響、保守コストから追加検証を選びます。
公開キーの読取は parser、座標・配置の展開は authoring、mesh・粒子の構築は runtime configuration が担当します。

## KUDPC で開発テストを実行する

KUDPC では test payload を実行する前に `hostname`、`module list`、利用可能なら `spartition` と `qgroup` を確認し、
host role、active `Sys*` module、割当を確定します。

login node（`camphor*`、`laurel*`、`cinnamon*`、`gardenia*`）では、編集、差分確認、短い log 読み取り、
`make check`、job 投入・監視までに留めます。`make test*`、`pytest`、`fpm test`、MPI / OpenMP test、
benchmark は直接実行しません。

- 短い開発 test は `tssrun` で計算 node の割当を取得して実行する。
- tier test、MPI test、heavy / release gate は `sbatch` job 内の `srun` で実行する。
- job 内でも複数の `fpm test` を並行実行しない。
- `build.sh` は fpm を `OMP_PROC_BIND=false` で起動する。fpm 自身が OpenMP program のため、KUDPC の job が既定で設定する
  `OMP_PROC_BIND=close` のままでは fpm の main thread が 1 core に固定され、fpm が起動する test・`beach`・MPI launcher が
  その 1 core の affinity を継承する。test payload の環境では `OMP_PROC_BIND=false` が見える。
- checkout、入力、cache、log は計算 node から見える `/home`、`/LARGE0`、`/LARGE1`、または合意した
  `/FAST` path に置く。login node の `/TMP` を前提にしない。

production simulation の job 設計、thread / rank 構成、負荷見積もりは[実行する](../guide/Execution.md)の責務です。

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
`output.diagnostics.charge_rel_change_threshold` は現行実装では監視・出力値であり、早期停止条件として説明しません。実行成功、test 通過、
数値収束、物理妥当性は別々に報告してください。

正本の責務分担は[開発者向けアーキテクチャ](Architecture.md#正本の責務を区別する)、user case の妥当性確認は
[計算結果の妥当性確認](../guide/ValidationGuide.md)を参照してください。
