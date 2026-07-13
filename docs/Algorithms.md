title: 計算モデルの全体像

Lang: [日本語](Algorithms.md) | [English](Algorithms.en.md)

# 計算モデルの全体像

BEACHは、三角形表面に蓄積した電荷が作る電場と、その電場中を運動する荷電粒子を
batch単位で結合する表面帯電シミュレータです。格子PICではなく、境界要素上の電荷を
sourceとして電場と電位を評価します。

このページでは計算全体の構成だけを説明します。個々の方程式、離散化、実装は各詳細ページに
分けています。

## 計算が保持する状態

batchをまたいで保持する主な状態は次のとおりです。

| 状態 | 内容 |
| --- | --- |
| 表面メッシュ | 三角形の頂点、重心、法線、面積、surface model |
| 要素電荷 | 各三角形に蓄積した総電荷`q_elem` |
| 注入状態 | reservoir粒子数の端数など、次batchへ持ち越す値 |
| 統計 | 吸収、escape、未解決粒子、処理済みbatch数 |
| 再開情報 | RNG状態、model/mesh/species fingerprint |

粒子は原則として各batchで生成され、そのbatch内で追跡されます。表面へ吸収された粒子の電荷は
要素電荷へ加算され、次batchの電場に反映されます。

## 1 batchの流れ

```text
現在の表面電荷
      ↓
電場・電位を評価
      ↓
粒子を生成
      ↓
粒子を前進 ── box境界／三角形との最初のeventを判定
      ↓
吸収粒子の電荷差分を集計
      ↓
要素電荷へcommit
      ↓
統計・履歴・再開状態を更新
```

電荷差分は粒子ごとに即時反映せず、batchの最後にまとめてcommitします。そのため、同じbatchの
粒子は同じbatch開始時の表面電荷から作られた場を共有します。詳しい順序は
[batch結合アルゴリズム](BatchAlgorithm.html)を参照してください。

## 物理モデル

### 表面との相互作用

現行の中心モデルは、粒子を表面で吸収し、その電荷を命中した絶縁体要素へ蓄積する
insulator accumulationです。浮遊導体の等電位緩和には対応しますが、`dielectric`の
`epsilon_r`は現行ではmetadataであり、独立した誘電分極境界条件ではありません。

表面衝突後の電荷堆積とsurface modelは[表面帯電モデル](SurfaceModels.html)を参照してください。

### 粒子源と外部境界

粒子はvolume seed、reservoir face、photo raycastなどから生成できます。open boundaryで外へ
向かう粒子は、選択した境界モデルに従ってescape、反射、returnのいずれかになります。

粒子の生成とbox境界は[粒子源と粒子境界](ParticleSourcesBoundaries.html)、シースとreturn closureは
[外部プラズマモデル](OuterPlasmaModels.html)で説明します。

### 周期領域と外部プラズマ

`periodic2`は2軸周期の電場境界です。近傍画像、無限周期の非zero mode、表面平均を表すzero mode、
必要に応じた外部プラズマ応答を区別して構成します。

周期場の分解は[periodic2場計算](PeriodicElectrostatics.html)、外部closureは
[外部プラズマモデル](OuterPlasmaModels.html)を参照してください。

## 数値手法

### 粒子追跡

荷電粒子の速度はBoris法、位置は同時刻状態を使う台形則で更新します。各stepでは移動線分に対する
box境界と三角形メッシュのeventを調べ、最も早いeventだけを採用します。

### 電場・電位

三角形要素の電荷は、重心の点電荷または三角形上の一定面密度として評価します。要素数と境界条件に
応じてdirect、treecode、FMMを選択できます。選び方は
[場ソルバーと境界条件](FieldSolvers.html)を参照してください。

### 時間スケール

`sim.dt`は粒子追跡のstep幅です。`batch_duration`は粒子供給量と表面電荷更新を結ぶ時間幅であり、
役割が異なります。後者の安定性と収束確認は
[`batch_duration`の安定性と定常値](BatchDurationStability.html)を参照してください。

## 文書の読み分け

| 知りたいこと | ページ |
| --- | --- |
| 1 batchの更新順序 | [batch結合アルゴリズム](BatchAlgorithm.html) |
| 表面吸収と電荷蓄積 | [表面帯電モデル](SurfaceModels.html) |
| 粒子生成、reservoir、box境界 | [粒子源と粒子境界](ParticleSourcesBoundaries.html) |
| Boris法と三角形衝突 | [粒子追跡と衝突](ParticleTrackingCollision.html) |
| direct/treecode/FMMの選択 | [場ソルバーと境界条件](FieldSolvers.html) |
| 周期場とzero mode | [periodic2場計算](PeriodicElectrostatics.html) |
| シース、outer plasma、escape/return | [外部プラズマモデル](OuterPlasmaModels.html) |
| FMMの選択と精度確認 | [FMMによる場計算](FMM.html) |
| FMMの展開と内部処理 | [Coulomb FMMコア詳細](FMMCore.html) |
| 離散化と結果の収束確認 | [計算結果の妥当性確認](ValidationGuide.html) |

設定キーを探す場合は[設定パラメータ](Parameters.html)を参照してください。
