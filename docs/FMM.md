title: FMMによる場計算

Lang: [日本語](FMM.md) | [English](FMM.en.md)

# FMMによる場計算

Fast Multipole Method (FMM)は、多数の境界要素が作るCoulomb場を、多数の粒子位置で評価するための
高速化手法です。このページでは利用者がsolverを選択し、精度を確認するために必要な内容を説明します。

## いつ使うか

| ケース | 推奨 |
| --- | --- |
| 小さなmesh、基準解 | `direct` |
| 中規模のfree境界 | `auto`または`treecode` |
| 大規模mesh、多数の評価点 | `fmm` |
| `periodic2` | `fmm`必須 |

FMMは近似法です。新しいmesh、kernel、周期設定では、小ケースのdirect解または独立oracleとの比較を先に
行ってください。

## 計算の構成

FMMはsourceをoctreeへ分け、遠いcell間の相互作用を多重極展開でまとめます。

```text
source charge
   ↓ P2M
multipole expansion
   ↓ M2M
upward tree
   ↓ M2L
local expansion
   ↓ L2L / L2P
target field
```

近接cellは展開せずdirectに評価します。これにより近距離精度を保ちつつ、遠距離相互作用をまとめます。

## planとstate

BEACHは固定mesh geometryとbatchごとに変化する要素電荷を分けます。

- plan: source位置、tree、interaction listなどの幾何依存情報
- state: 現在の`q_elem`から作るmultipole/local係数

通常はplanを初期化時に一度構築し、各batchではstateだけを更新します。mesh geometryや要素数が変わる場合は
planの再構築が必要です。

## source kernel

| `field.element_kernel` | FMMでのsource |
| --- | --- |
| `point` | 要素重心のsoftened point charge |
| `triangle_p0` | 三角形上の一定面密度 |

`point`と`triangle_p0`は異なる離散化です。FMM orderだけでなく、source meshとkernel自体のrefinementも
独立に確認してください。

## 精度を支配するもの

- 展開order
- treeのleaf sizeとwell-separated判定
- point kernelのsoftening
- `triangle_p0`のnear correction
- periodic image layersとfar correction
- source mesh解像度

あるorderでFMM誤差が小さくても、粗いsource meshの離散化誤差が小さいとは限りません。

## periodic2

`periodic2`ではFMMの有限近傍画像に、必要に応じて`cached_kneq0`遠方operatorを加えます。zero modeと
outer responseはFMM coreだけで完結せず、electrostatic snapshot側で合成します。

詳しくは[periodic2場計算](PeriodicElectrostatics.html)を参照してください。

## 推奨する確認

1. 小meshでdirectとFMMを比較する
2. FMM orderを上げて場と電位が収束することを確認する
3. source meshを細かくして結論が収束することを確認する
4. periodic2ではcache cold/warm一致とzero-mode診断を確認する
5. performance比較はrelease buildで行う

内部API、tree構築、Cartesian展開、M2L、cache operatorの実装は
[Coulomb FMMコア内部実装](FMMCore.html)を参照してください。
