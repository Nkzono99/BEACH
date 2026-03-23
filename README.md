# BEACH（BEM + Accumulated CHarge）

[![GitHub Pages](https://img.shields.io/website?url=https%3A%2F%2Fnkzono99.github.io%2FBEACH%2F&up_message=GitHub%20Pages&down_message=Pages%20down)](https://nkzono99.github.io/BEACH/)
[![Fortran Docs](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml)
[![Fortran Format](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml)
[![PyPI version](https://img.shields.io/pypi/v/beach-bem)](https://pypi.org/project/beach-bem/)
[![Python versions](https://img.shields.io/pypi/pyversions/beach-bem)](https://pypi.org/project/beach-bem/)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)

公開ドキュメント: [GitHub Pages](https://nkzono99.github.io/BEACH/)

BEACH は、**BEM（境界要素法）ベースの表面帯電 + テスト粒子追跡シミュレータ**です。  
計算本体は Fortran（`fpm` 管理）、Python は後処理・可視化を担当します。

v0.x の主対象は **insulator accumulation（絶縁体への電荷蓄積）** です。

## 推奨ワークフロー

通常利用では、**`pip install beach-bem` で導入して `beach` を実行する運用**を推奨します。  
`pip install -e` や `make` は、開発に参加する場合の手順として扱います。

## 1. 利用者向けセットアップ

### 1.1 前提ツール

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

開発版を直接試したい場合は、Git URL からも導入できます。

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

このインストールでは、ビルド時に `make install` が実行され、Python 側と Fortran 実行バイナリ（`beach`）が同時に導入されます。
ユーザーインストール時は `~/.local/bin` に入るため、必要なら PATH を追加してください。
pip 経由のビルドでは既定で `INSTALL_PROFILE=auto` を使います。必要なら `INSTALL_PROFILE=camphor` などを環境変数で上書きできます。
`auto` 失敗時は既定で `generic` へフォールバックします。無効化したい場合は `BEACH_PIP_FALLBACK_GENERIC=0` を指定してください。

```bash
export PATH="$HOME/.local/bin:$PATH"
```

## 2. 実行

### 2.0 `case.toml` から `beach.toml` を生成する

`beachx config` を使うと、preset の組合せとローカル override から実行用 `beach.toml` を生成できます。
日常的な編集対象は `case.toml`、Fortran 実行系が読むのは生成後の `beach.toml` です。

```bash
mkdir run_periodic2
cd run_periodic2

beachx config init
beachx config validate
beachx config render
beach beach.toml
```

`beachx config init` の既定値は、次の built-in preset を組み合わせた `case.toml` を作ります。

- `sim/periodic2_fmm`
- `species/solarwind_electron`
- `species/solarwind_ion`
- `mesh/plane_basic`
- `output/standard`

preset の探索順は `./.beachx/presets/` → `~/.config/beachx/presets/` → package 同梱 preset です。
研究室ローカルやプロジェクトローカルの設定は、この preset 置き場へ `sim/...` や `mesh/...` として TOML 断片を追加して再利用できます。
サンプルの `case.toml` は [`examples/periodic2_basic/case.toml`](examples/periodic2_basic/case.toml) にあります。
生成される `case.toml` / preset / `beach.toml` にはそれぞれ対応する `#:schema ...` directive を自動で付けています。

user preset を CLI から作る場合は、たとえば次のように使えます。

```bash
beachx preset new sim/lab/periodic2_fast --from sim/periodic2_fmm
beachx preset new output/project/debug --local
beachx preset save sim/lab/current_run --section sim
beachx preset save species/lab/ion --section species --index 2 --rendered
beachx preset edit sim/lab/periodic2_fast
beachx preset list
beachx preset show sim/lab/periodic2_fast
```

### 2.1 推奨: `beach` コマンド

```bash
beach examples/beach.toml
```

引数なしの場合は、カレントディレクトリの `beach.toml` を自動読込します。

### 2.2 出力確認

```bash
ls outputs/latest
beachx inspect outputs/latest
```

主な出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `charge_history.csv`（`history_stride > 0` のとき）

### 2.3 可視化

```bash
beachx inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png

# sim.field_bc_mode = "periodic2" の mesh を周期セルへ寄せて描く場合
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

`beachx coulomb` は、`beach.toml` の `mesh.templates` が見つかれば object の kind と順序をそこから読み取り、既定では全 object を target 軸に並べて行列を描きます。特定 kind だけに絞りたいときは `--target-kinds sphere` のように指定できます。
`beachx mobility` は、既定で `plane` を support とみなし、それ以外の object について合力・合トルクと `lift_ratio` / `slide_ratio` / `roll_ratio` を CSV に書き出します。質量由来の指標は `--density-kg-m3` と `beach.toml` の幾何情報が必要です。

旧 alias の `beach-inspect` / `beach-animate-history` / `beach-plot-coulomb-force-matrix` /
`beach-plot-potential-slices` なども当面は使えますが、今後は `beachx ...` を推奨します。

Python から `mesh source` ごとの面積重み付き箱ひげ図を作る例:

```python
from beach import Beach

run = Beach("outputs/latest")

# 電荷 [C] の面積重み付き箱ひげ図（mesh sourceごと）
fig_q, ax_q = run.plot_mesh_source_boxplot(
    quantity="charge",
    step=-1,          # 最新 history。None なら charges.csv を使う
    showfliers=False,
)
fig_q.savefig("outputs/latest/charge_box_by_source.png", dpi=150)

# periodic2 mesh を周期セルへ寄せて表示
run.plot_mesh(apply_periodic2_mesh=True)
run.plot_potential(apply_periodic2_mesh=True)

# ポテンシャル [V] の面積重み付き箱ひげ図（mesh sourceごと）
fig_phi, ax_phi = run.plot_mesh_source_boxplot(
    quantity="potential",
    softening=0.0,
    self_term="area_equivalent",
)
fig_phi.savefig("outputs/latest/potential_box_by_source.png", dpi=150)

# object ごとの Coulomb 力行列（beach.toml があれば plane/sphere などを自動ラベル化）
fig_f, ax_f = run.plot_coulomb_force_matrix(
    component="z",
)
fig_f.savefig("outputs/latest/coulomb_force_z.png", dpi=150)

# 移動可能性の簡易評価（既定では plane 以外を対象）
mobility = run.analyze_coulomb_mobility(
    density_kg_m3=2500.0,
    mu_static=0.4,
)
for record in mobility.records:
    print(record.label, record.lift_ratio, record.slide_ratio, record.roll_ratio)
```

## 3. 運用でよく使うコマンド

### 3.1 負荷見積もり

```bash
beachx workload examples/beach.toml --threads 8
```

### 3.2 再開実行

```toml
[output]
dir = "outputs/latest"
resume = true
```

同じ `output.dir` で `beach` を再実行すると続きから計算します。

## 4. 開発に携わる場合

ローカル開発では、Python と Fortran を分けて運用するのが便利です。

### 4.1 Python 側（編集可能インストール）

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

### 4.2 Fortran 実行系（`make`）

```bash
make
```

`make` のデフォルトターゲットは `install` で、`BUILD_PROFILE=auto`（ホスト自動判定）を使います。
必要に応じて明示指定できます。

```bash
make install-generic
make install-camphor
```

### 4.3 `fpm` 直接実行（開発向け）

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

MPI + OpenMP（`USE_MPI` 有効化）:

```bash
FPM_FC=mpiifort \
fpm run --profile release --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 4" -- examples/beach.toml
```

### 4.4 テスト

```bash
pytest -q
fpm test
```

### 4.5 Fortran 整形

Fortran の `*.f90` / `*.F90` は `fprettify -i 2` を標準 formatter とします。
生成バックアップの `*.i90` は整形対象外です。

ローカルで hook を有効にするには:

```bash
python -m pip install pre-commit
make install-hooks
```

手動整形と CI 相当の確認:

```bash
make fmt-fortran
make fmt-check-fortran
```

GitHub Actions でも同じ `pre-commit` 設定を `--all-files` で実行し、整形漏れを検知します。

### 4.6 Fortran ドキュメント生成

GitHub Pages 向けの Fortran ドキュメントは FORD で生成します。
モジュール API、`use` 依存グラフ、`docs/` 配下の運用ドキュメントをまとめて公開できます。

```bash
python -m pip install -r docs/requirements.txt
# Ubuntu / Debian 系
sudo apt-get install -y graphviz

make docs-fortran
```

出力先は `build/ford-docs/` です。
`main` への push では GitHub Actions の `Fortran Docs` workflow が同じ手順で GitHub Pages を更新します。

## 5. プロジェクト構成

- [`src/`](src/), [`app/`](app/): Fortran 本体
- [`tests/fortran/`](tests/fortran/): Fortran テスト
- [`beach/`](beach/): Python 後処理ライブラリ
- [`tests/python/`](tests/python/): Python テスト
- [`examples/`](examples/): 設定例・補助スクリプト
- [`docs/`](docs/): 運用ドキュメント
- [`SPEC.md`](SPEC.md): 現行 Fortran 実装仕様

## 6. ドキュメント案内

- Fortran 実行フロー: [`docs/fortran_workflow.md`](docs/fortran_workflow.md)
- `beach.toml` 仕様: [`docs/fortran_parameter_file.md`](docs/fortran_parameter_file.md)
- `beach.toml` JSON Schema: [`schemas/beach.schema.json`](schemas/beach.schema.json)
- Coulomb FMM コア仕様: [`docs/fortran_fmm_core.md`](docs/fortran_fmm_core.md)
- 自動生成の依存関係マップ: `make docs-fortran` 実行後の `docs/fortran_dependency_map.md`
- 実装仕様（source of truth）: [`SPEC.md`](SPEC.md)
- 設定サンプル: [`examples/beach.toml`](examples/beach.toml)
