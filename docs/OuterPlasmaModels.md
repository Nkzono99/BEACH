title: 外部プラズマモデル

Lang: [日本語](OuterPlasmaModels.md) | [English](OuterPlasmaModels.en.md)

# 外部プラズマモデル

外部プラズマモデルは、有限な粒子計算領域の外側をどの状態量で閉じるかを定めます。
流入VDFだけを補正するclosure、interfaceに接続する1D profile、外部粒子を進める3D fieldは、
reservoir流入と外向き粒子に与えられる情報がそれぞれ異なります。`open`境界は外部状態を持たず、通過した粒子をその場で除去します。

## 外部へ持つ状態量でモデルを分類する

| モデル | 解くもの | particle pusherへ空間場を与えるか |
| --- | --- | --- |
| `infinity_barrier` | face平均電位による流入cutoff | いいえ |
| `floating_no_photo` | electron/ion流入の簡易電流釣合い | いいえ |
| `zhao_*` | 解析sheath closureによる注入VDF補正 | いいえ |
| `linear_debye` | zero-modeの指数応答 | outer returnに使用 |
| `kinetic_1d` | 無限遠VDFとPoisson方程式の1D profile | outer領域で使用 |
| `unified_linear_response` | rough surfaceを含む線形1D response | localからfarまで合成 |

`floating_no_photo`と`zhao_*`が返す量はsource VDFの補正です。[シース注入closure](SheathInjectionClosures.html)で、
branch、cutoff、適用範囲を説明します。

## split modelはinterfaceでlocal場と1D profileを接続する

`linear_debye`と`kinetic_1d`は、meshを含むlocal領域と、ownership interfaceより外側の1D profileを
接続します。表面が作るzero-mode fieldをinterface条件とし、outer側の電位差から流入粒子の加減速と
外向き粒子のescape/returnを決めます。

`kinetic_1d`は元のPoisson residual、単調分枝、ion accessibility、Bohm entry、無限遠準中性をすべて
満たした解だけを受理します。失敗時に別モデルへsilent fallbackしません。

Poisson問題、VDF density closure、Newton solve、受理条件は
[kinetic 1D外部プラズマ](KineticOuterPlasma.html)にまとめています。

## unified modelは表面から遠方まで一つの場を解く

`unified_linear_response`は、rough surfaceの高さ範囲、表面平均source、plasmaが占有できる面積率、
線形Debye応答を同じ1D gridで扱います。surfaceとinterfaceの間に、真空とみなせるsplit windowがない場合の
線形検証に使います。species VDF、Bohm条件、浮遊電流balanceまで解くのは`kinetic_1d`です。

accessible fraction、zero-mode Poisson solve、nonzero-modeのreflection/transmission、線形性gateは
[unified linear response](UnifiedLinearResponse.html)で説明します。

## 同じouter stateでescapeとreturnを決める

| transfer | 処理 |
| --- | --- |
| instant return | 保存エネルギーからturning pointと往復時間を求め、同じsimulation時刻へ写像 |
| kinetic profile return | 収束した離散$\phi(z)$上でturning pointとflight timeを積分 |
| explicit 3D orbit | zero/nonzero field中を固定stepで追跡 |

instant returnは、定常または準定常sheath向けのreduced closureです。outer flightの遅延をglobal timeに加えないため、
立ち上がりやpulseなど、往復時間と同程度の時間で変化する過渡電流には使いません。

interfaceでのenergy判定、linear/kinetic profileのflight time、3D軌道、frozen-field gateは
[粒子のescapeとreturn](ParticleEscapeReturn.html)で説明します。

## 必要な外部物理から構成を選ぶ

- 単純な有限box: outer modelなし、`open_boundary_model="escape"`
- 過去のscalar barrier再現: `infinity_barrier`
- 文献の注入closure比較: `zhao_*`
- 無限周期表面と自己整合1D sheath: `kinetic_1d`
- roughnessを含む線形応答: `unified_linear_response`
- outer軌道の横方向変化が重要: explicit 3D orbit

## 選んだモデルの受理条件を確認する

profile statusだけでなく、Poisson residual、境界条件、電荷積分、current balance、grid refinement、
frozen-field ratioを確認します。物理的な適用範囲はADR 0001/0002に、数値的な確認手順は
[計算結果の妥当性確認](ValidationGuide.html)にまとめています。

粒子の流入側は[reservoir注入](ReservoirInjection.html)、光電子sourceとの関係は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)で説明します。
