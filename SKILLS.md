# SKILLS.md — BEACH 開発ナレッジ

このセッションで確立した構造知識・ワークフロー・注意点をまとめる。

## モジュール依存レイヤー（リファクタ後）

```
Layer 0  bem_kinds
Layer 1  bem_constants (pi, eps0, qe, k_coulomb, k_boltzmann)
         bem_string_utils (lower_ascii)
Layer 2  bem_types, bem_coulomb_fmm_types
Layer 3  bem_app_config_types
Layer 4  bem_field, bem_pusher, bem_boundary, bem_particles
Layer 5  bem_mesh, bem_importers, bem_templates
Layer 6  bem_collision, bem_injection
Layer 7  bem_config_helpers (resolve_inject_face, resolve_inward_normal, etc.)
         bem_coulomb_fmm_* (internal/)
Layer 8  bem_coulomb_fmm_core
Layer 9  bem_field_solver (+config/eval/tree submodules)
         bem_app_config_parser (+validate/parse_utils submodules)
Layer 10 bem_sheath_model_core, bem_sheath_runtime
Layer 11 bem_app_config_runtime
Layer 12 bem_app_config (ファサード)
Layer 13 bem_output_writer (純粋I/O、物理計算なし)
         bem_performance_profile, bem_restart
Layer 14 bem_simulator (+loop/io/stats submodules)
Layer 15 app/main.f90
```

### 依存ルール
- 上位層は下位層のみを use する（逆方向依存禁止）
- `bem_output_writer` は FMM 内部型に依存しない（mesh 電位計算は `field_solver` 側）
- `bem_sheath_*` は `bem_app_config_parser` に依存しない（`bem_config_helpers` 経由）
- 物理定数 (`pi`, `eps0`, `qe`) は `bem_constants` に一元管理

## テスト構造

### 重量別テストマップ

| 重さ | テスト | 実行時間目安 |
|------|--------|-------------|
| 軽量 | `test_dynamics_basic` | ~1秒 |
| 軽量 | `test_output_writer_io` | ~1秒 |
| 軽量 | `test_boundary` | ~1秒 |
| 軽量 | `test_performance_profile` | ~1秒 |
| 中量 | `test_dynamics_field_solver` | ~数秒 |
| 中量 | `test_coulomb_fmm_core_basic` | ~数秒 |
| 中量 | `test_output_writer_potential` | ~数秒 |
| 中量 | `test_app_config_parser` | ~数秒 |
| 中量 | `test_injection_sampling` | ~数秒 |
| 中量 | `test_reservoir_injection` | ~数秒 |
| 中量 | `test_sheath_model_core` | ~数秒 |
| 中量 | `test_sheath_injection_model` | ~数秒 |
| 中量 | `test_templates_importers_runtime` | ~数秒 |
| 中量 | `test_simulator` | ~数秒 |
| 中量 | `test_restart` | ~数秒 |
| 重量 | `test_dynamics_fmm` | ~30秒+ |
| 重量 | `test_coulomb_fmm_core_periodic` | ~30秒+ |
| 重量 | `test_periodic2_flat_oracle_diag` | ~30秒+ |

### テスト実行パターン

```bash
# 軽量テストだけ素早く回す（変更の健全性チェック）
fpm test --target test_dynamics_basic
fpm test --target test_boundary

# 特定モジュールの変更後
fpm test --target test_output_writer_potential  # output_writer 変更後
fpm test --target test_sheath_model_core        # sheath 変更後
fpm test --target test_app_config_parser        # parser 変更後

# 全テスト（時間がかかる。バックグラウンド推奨）
fpm test
```

### テストランナー API (`test_support`)

```fortran
program test_example
  use test_support, only: test_init, test_begin, test_end, test_summary, &
    assert_true, assert_close_dp
  implicit none

  call test_init(3)             ! 予定テスト数

  call test_begin('basic')      ! [1/3] basic ...
  call assert_true(.true., 'should pass')
  call test_end()               ! PASS

  call test_begin('math')       ! [2/3] math ...
  call assert_close_dp(1.0d0, 1.0d0, 1.0d-12, 'one equals one')
  call test_end()               ! PASS

  call test_summary()           ! === 2 passed, 0 failed (2 assertions) ===
end program
```

失敗時は actual/expected/diff/tol が自動表示される。

## リファクタリング進捗

### 完了済み (Phase 1)

| Phase | 内容 | コミット |
|-------|------|---------|
| 1-A | `lower_ascii` を `bem_string_utils` に集約 | `4181f45` |
| 1-B | sheath 層違反解消 + 物理定数集約 | `88e0c12` |
| 1-C | `bem_output_writer` の FMM 依存切り離し | `338ff8e` |

### 未着手 (Phase 2 — 中規模)

| Phase | 内容 | リスク |
|-------|------|--------|
| 2-A | `bem_app_config_runtime` (900行) の責務分離 → mesh_builder / particle_init | 中 |
| 2-B | `field_solver_type` (~70フィールド) のサブ型分割 | 高 |
| 2-C | Python レガシー CLI ラッパー除去 | 低（破壊的変更） |

### 触らない方がよい箇所

- **`sim_config` 型**: コア全モジュールが依存。Phase 2 完了後に再検討
- **FMM 内部 (`fmm/internal/`)**: 一方向依存で整理済み。field_solver 側で吸収
- **`mesh_type` フィールド変更**: restart チェックポイントのバイナリ互換に影響
- **`app_config` フィールド名変更**: TOML パーサとの整合が崩れるリスク

## 既知の問題

### `test_coulomb_fmm_core_periodic` — periodic2 FMM 精度テスト (Bug #2)

`test_periodic2_field_accuracy` が相対誤差 5% 閾値を超える場合がある。
`sort_real_ascending` の OOB バグ (Bug #1) は修正済みだが、
periodic2 FMM の精度自体が閾値ギリギリの可能性あり。要調査。

## Fortran submodule 操作の注意

- submodule の `module procedure` を外部モジュールに移動する場合、
  親モジュールの `interface` 宣言も同時に削除する必要がある
- submodule は独自の `use` 文を持てる（親モジュールの use とは独立）
- 親モジュールを変更すると全 submodule が再コンパイルされる
