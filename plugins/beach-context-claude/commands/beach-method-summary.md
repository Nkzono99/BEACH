---
description: "BEACH の手法説明やアルゴリズム要約を作る"
argument-hint: "<audience/topic>"
---

`beach-method-summary` agent を使い、指定された audience と topic に合わせて BEACH の手法説明を作ってください。

優先参照:
- `docs/Algorithms.md`
- `docs/FieldSolvers.md`, `docs/FMM.md`, `docs/FMMCore.md`
- `docs/PeriodicElectrostatics.md`, `docs/ParticleEvents.md`, `docs/SurfaceModels.md`
- `SPEC.md`
- `README.md`
- 必要なら `src/` と `app/`

`docs/Algorithms.md`は索引として使い、absorption-only interaction、insulator accumulation、
direct/treecode/FMM、periodic2、batch loopの範囲を
正確に区別してください。実装確認なしに未実装 model を有効機能として書かないでください。
