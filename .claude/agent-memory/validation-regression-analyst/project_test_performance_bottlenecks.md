---
name: Fortran test performance bottlenecks
description: Periodic2 Ewald reference solutions in FMM tests are O(nelem * k_count * screen_count) per query point, causing tests to take minutes
type: project
---

FMM/periodic2 テストのボトルネックは `electric_field_at_periodic2_reference` と `add_periodic2_exact_ewald_correction_all_sources` の参照解計算。

- `test_dynamics_fmm`: 528要素球メッシュ, 10サブテスト, テスト6(fallback)とテスト7(image_layers)がそれぞれ5分以上。全体で16分以上（CPU時間ベース）。
- `test_coulomb_fmm_core_basic`: 1分15秒。`state_update_reuse` の periodic2 M2L キャッシュ構築が主因。
- `test_coulomb_fmm_core_periodic`: 1分11秒。各テストで periodic2 plan 構築 + Ewald 参照解。
- `test_output_writer_potential`: 49秒。`fmm_core_mesh_potential` テストで periodic2 solver init。
- `test_periodic2_flat_oracle_diag`: 800要素メッシュ, 2つのsolver + 8点参照解。推定5分以上。

Ewald 参照解の計算量: nelem * (screen_count + k_count + 1) * n_query。screen_count = (2*(nimg+kmax)+1)^2, k_count = (2*kmax+1)^2 - 1。

**Why:** `fpm test` の全体時間が20分以上になり、開発フィードバックループが極端に遅い。
**How to apply:** periodic2/FMM テストの高速化提案時に参照。テストの正当性（参照解の精度）を損なわずに高速化するには、メッシュ要素数の削減またはewald_layersの削減が有効。
