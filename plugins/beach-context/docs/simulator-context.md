# BEACH simulator context

BEACH は BEM による表面帯電とテスト粒子追跡を組み合わせた simulator です。計算本体は Fortran (`src/`, `app/`) で、Python package (`beach/`) は設定補助、出力読み込み、解析、可視化を担当します。

現行 v1 系の中心仕様:

- 三角形境界要素に電荷を保持する。
- 現在の要素電荷から静電場を評価する。
- 粒子は Boris pusher で前進する。
- 線分と三角形の最初の交差を衝突として扱う。
- 衝突粒子は吸収され、`q_particle * w_particle` が hit element に堆積する。
- batch ごとに電荷差分、統計、履歴を commit する。
- surface model は insulator accumulation が標準。conductor/resistive model は予約・拡張点。

設定入口:

- 普段は `beach.toml` を直接編集し、必要に応じて `beachx config render` で高水準記法を最終キーへ展開する。
- Fortran 実行系が直接読むのは `beach.toml`。
- `beachx lint beach.toml` は TOML parse、JSON Schema、高水準記法、既知制約を含む事前検査。
- 最終キー仕様は `references/fortran_parameter_file.md`。

主要な出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv`
- `potential_history.csv`
- `mesh_potential.csv`
- `rng_state.txt`
- `macro_residuals.csv`
- `performance_profile.csv` (`BEACH_PROFILE=1`)

解析入口:

- CLI: `beachx inspect`, `beachx animate`, `beachx slices`, `beachx coulomb`, `beachx mobility`, `beachx kernel-forces`
- Python: `from beach import Beach`
- Fortran field kernel: `make build-kernel` 後に `beachx kernel-forces`

実行上の注意:

- `sim.batch_count` は今回実行する batch 数。
- `sim.max_step` は粒子ごとの上限 step。
- `sim.tol_rel` は監視値であり停止条件ではない。
- `output.resume = true` では同じ `output.dir` の `summary.txt`, `charges.csv`, `rng_state.txt` を使う。MPI では前回と同じ world size が必要。
- `field_bc_mode = "periodic2"` は 2 軸周期境界を扱い、FMM と box boundary 設定の組み合わせに注意する。
