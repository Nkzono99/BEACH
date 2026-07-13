title: 粒子追跡と表面電荷蓄積（移動しました）

Lang: [日本語](ParticleChargeLoop.md) | [English](ParticleChargeLoop.en.md)

# 粒子追跡と表面電荷蓄積（移動しました）

この旧統合ページの内容は、batch loopを入口として工程ごとの詳細ページへ移しました。
計算順序は[計算モデルの全体像](Algorithms.html)から確認してください。

| 旧ページで扱っていた内容 | 移動先 |
| --- | --- |
| 粒子生成、reservoir、photo raycast | [粒子源](ParticleSourcesBoundaries.html)、[reservoir注入](ReservoirInjection.html)、[光電子の放出とライフサイクル](PhotoelectronEmission.html) |
| Boris速度更新と位置更新 | [Boris粒子更新](BorisPusher.html) |
| 三角形衝突、box面、periodic境界 | [粒子の衝突・境界イベント](ParticleEvents.html) |
| escape、反射、outer return | [粒子のescapeとreturn](ParticleEscapeReturn.html) |
| 電荷堆積、surface model、MPI commit | [表面電荷更新](SurfaceModels.html) |
| 統計、履歴、restart | [出力の読み方](OutputGuide.html)、[実行する](Execution.html) |
| OpenMP、MPI、性能計測 | [実行する](Execution.html) |
