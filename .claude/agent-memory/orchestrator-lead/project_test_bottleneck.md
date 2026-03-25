---
name: FMM Ewald テストボトルネック
description: Fortranテスト全体の95%以上がperiodic2 Ewald参照解計算に支配されている。高速化戦略あり。
type: project
---

`fpm test` 全19ターゲットの合計実行時間は約53分。うち95%以上がFMM/periodic2関連の5テストに集中。

確定タイミング（2026-03-25計測）:
- test_dynamics_fmm: **38分54秒**（最遅、10テスト全パス）
- test_periodic2_flat_oracle_diag: **11分33秒**（800要素flat plane、1テスト全パス）
- test_coulomb_fmm_core_basic: **1分15秒**
- test_coulomb_fmm_core_periodic: **1分11秒**
- test_output_writer_potential: **48秒**

ボトルネック: `add_periodic2_exact_ewald_correction_all_sources` の O(N_src × K²) 計算。

**Why:** テスト参照解にfull Ewald sum（ewald_layers=4, screen_count=121, k_count=80）を使用しており、各query pointでN_src×200回の特殊関数呼出が発生。

**How to apply:** テストメッシュの要素数削減（528→132）+ ewald_layers削減（4→2）で10-15倍高速化可能。テストの科学的正当性は維持される（誤差閾値が5e-2程度と緩い）。
