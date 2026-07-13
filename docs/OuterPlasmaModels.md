title: 外部プラズマモデル

Lang: [日本語](OuterPlasmaModels.md) | [English](OuterPlasmaModels.en.md)

# 外部プラズマモデル

外部プラズマモデルは、有限な粒子計算領域の外側をどの物理近似で閉じるかを定めます。単なる
`open`境界、注入VDFの補正、1D Poisson solve、3D outer orbitを区別してください。

## モデルの分類

| モデル | 解くもの | particle pusherへ空間場を与えるか |
| --- | --- | --- |
| `infinity_barrier` | face平均電位による流入cutoff | いいえ |
| `floating_no_photo` | electron/ion流入の簡易電流釣合い | いいえ |
| `zhao_*` | 解析sheath closureによる注入VDF補正 | いいえ |
| `linear_debye` | zero-modeの指数応答 | outer returnに使用 |
| `kinetic_1d` | 無限遠VDFとPoisson方程式の1D profile | outer領域で使用 |
| `unified_linear_response` | rough surfaceを含む線形1D response | localからfarまで合成 |

## split outer model

`linear_debye`と`kinetic_1d`は、meshを含むlocal領域と、ownership interfaceより外側の1D profileを
接続します。表面が作るzero-mode fieldをinterface条件とし、outer側の電位差から流入粒子の加減速と
外向き粒子のescape/returnを決めます。

`kinetic_1d`は元のPoisson residual、単調分枝、ion accessibility、Bohm entry、無限遠準中性をすべて
満たした解だけを受理します。失敗時に別モデルへsilent fallbackしません。

## unified linear response

`unified_linear_response`は、rough surfaceの高さ範囲、表面平均source、plasmaが占有できる面積率、
線形Debye応答を同じ1D gridへ入れます。surfaceとinterfaceの間を真空とみなせるsplit windowがない場合の
線形検証に使います。

これは非線形VDF sheath solverではなく、Bohm条件や浮遊電流balanceを解きません。

## 粒子return

| transfer | 処理 |
| --- | --- |
| instant return | 保存エネルギーからturning pointと往復時間を求め、同じsimulation時刻へ写像 |
| kinetic profile return | 収束した離散$\phi(z)$上でturning pointとflight timeを積分 |
| explicit 3D orbit | zero/nonzero field中を固定stepで追跡 |

instant returnは定常・準定常sheathのreduced closureです。outer flightの遅延をglobal timeへ加えないため、
立ち上がり、pulse、往復時間と同程度で変化する過渡電流には使いません。

## 選択の目安

- 単純な有限box: outer modelなし、`open_boundary_model="escape"`
- 過去のscalar barrier再現: `infinity_barrier`
- 文献の注入closure比較: `zhao_*`
- 無限周期表面と自己整合1D sheath: `kinetic_1d`
- roughnessを含む線形応答: `unified_linear_response`
- outer軌道の横方向変化が重要: explicit 3D orbit

## 必須診断

profile statusだけでなく、Poisson residual、境界条件、電荷積分、current balance、grid refinement、
frozen-field ratioを確認します。物理的な適用範囲はADR 0001/0002と
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。

solver式とmodel別contractの詳細は旧統合文書
[periodic2 zero modeと外部プラズマ](PeriodicZeroModeOuterPlasma.html#4-outer-plasma-modelとの接続)、
粒子写像は[外部シースとreservoir粒子境界](SheathReservoirBoundary.html)に残しています。
