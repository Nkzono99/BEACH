---
name: test_injection_sampling photo_raycast_periodic2 FAIL
description: periodic2境界条件下のphoto_raycast放出位置ラッピングに関するテスト失敗（pre-existing、main上で再現）
type: project
---

`test_injection_sampling` の `photo_raycast_periodic2` テスト（テスト14/16）が FAIL。

失敗箇所: `tests/fortran/test_injection_sampling.f90` 行220
メッセージ: "photo_raycast periodic2 should use wrapped emission position"
原因: periodic2境界条件下でphoto emissionの放出位置y座標がprimary cell [0,1)内にラップされていない。

**Why:** git commit 7b309be "Fix pre-existing test failures in sort and photo_raycast" で一部修正されたが、このケースはまだ修正されていない。

**How to apply:** periodic2 photo_raycast の emission position wrapping ロジック（`bem_injection.f90`内）を調査し修正が必要。今回の並列化変更とは無関係。
