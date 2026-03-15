# BEACH（BEM + Accumulated CHarge）

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
```

`beachx coulomb` は、`beach.toml` の `mesh.templates` が見つかれば object の kind と順序をそこから読み取り、sphere があれば既定で sphere 群を target にして行列を描きます。

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
    target_kinds=("sphere",),  # 省略時は sphere があれば自動選択
)
fig_f.savefig("outputs/latest/coulomb_force_z.png", dpi=150)
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
- Coulomb FMM コア仕様: [`docs/fortran_fmm_core.md`](docs/fortran_fmm_core.md)
- 実装仕様（source of truth）: [`SPEC.md`](SPEC.md)
- 設定サンプル: [`examples/beach.toml`](examples/beach.toml)
