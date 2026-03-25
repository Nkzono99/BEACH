---
name: test_injection_sampling photo_raycast_periodic2 failure
description: The photo_raycast_periodic2 test case in test_injection_sampling fails - periodic wrapping of emission position is broken
type: project
---

`test_injection_sampling` の 14/16 番目のテスト `photo_raycast_periodic2` が FAIL する。

失敗アサーション: `photo_raycast periodic2 should use wrapped emission position` (行220)
- 期待: emitted position y < 1e-6 (primary cell 内にラップされるべき)
- 実際: y がラップされていない

テストは periodic2 境界条件下で z_high 面から ray を発射し、wrap された位置にある三角形に命中する場合の emission position がプライマリセル内に収まることを検証。

**Why:** periodic2 境界でのphoto emission wrapping に潜在的バグがある可能性。
**How to apply:** `bem_injection` モジュールの `sample_photo_raycast_particles` 内の periodic2 wrapping ロジックを確認する必要がある。
