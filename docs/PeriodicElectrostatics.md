title: periodic2場計算

Lang: [日本語](PeriodicElectrostatics.md) | [English](PeriodicElectrostatics.en.md)

# periodic2場計算

`periodic2`はx/yを周期、zを開放とする静電場境界です。有限個の周期画像を足すだけではなく、
横方向Fourier modeと開放方向の境界条件を分けて扱います。

## 場の分解

| 成分 | 意味 | 評価方法 |
| --- | --- | --- |
| near images | primary cell近傍の強い局所場 | direct/tree/FMMの有限画像和 |
| `k != 0` far field | x/y方向に変化する無限周期遠方場 | `cached_kneq0` operator |
| surface `k = 0` | 平面平均された表面電荷の場 | 高さ方向の累積電荷 |
| outer response | interfaceより外側のplasma応答 | 選択したouter model |

最終的なelectrostatic snapshotは、これらを契約に従って一度ずつ合成します。

## finite imageとcached operator

`field_periodic_far_correction="none"`では、設定した`periodic_images`内だけを評価します。これは有限画像
モデルであり、画像数を増やさず無限周期場とみなしてはいけません。

`cached_kneq0`では、Ewald2P referenceからnear imageの寄与とsurface zero modeを除いた滑らかな遠方場を、
root-level線形operatorとして事前計算します。同じgeometry、kernel、softening、periodic cell、FMM設定なら
cacheを再利用できます。

cacheは`k != 0`だけを保存します。zero modeをもう一度含めると平均場を二重加算するため、snapshot構築側で
一度だけ合成します。

## surface zero mode

三角形の面積を高さ方向へ投影し、各高さより下にある累積電荷から平面平均場を作ります。非中性cellでは
z遠方に一定場と線形電位が残り得るため、lower/upper boundary closureを明示する必要があります。

zero modeは数値上消してよい成分ではありません。Gauss residualと境界条件を満たす物理成分として扱います。

## 粒子・衝突との関係

場評価点は周期cellへwrapしますが、軌道eventの物理位置は失いません。mesh衝突は必要なperiodic imageを
幾何的に探索します。field image countとcollision image countは同じ設定値を共有するとは限りません。

## 選択の目安

| 目的 | 設定 |
| --- | --- |
| 小さな有限画像比較 | `periodic_images` + far correctionなし |
| 無限周期production | `cached_kneq0` + 明示したzero/outer closure |
| solver検証 | image refinement、cache cold/warm一致、Ewald/oracle比較 |

## 確認する診断

- cache fingerprintとcold/warm結果の一致
- primary/near/far/zeroの二重加算がないこと
- Gauss residualとlower/upper boundary closure
- image数またはreference operatorに対する収束
- non-neutral cellで有限高さの電位差を無限遠escape energyと解釈していないこと

設定の入口は[場ソルバーと境界条件](FieldSolvers.html)、outer closureは
[外部プラズマモデル](OuterPlasmaModels.html)を参照してください。Ewald/operatorの実装詳細は旧統合文書
[periodic2 zero modeと外部プラズマ](PeriodicZeroModeOuterPlasma.html)に残しています。
