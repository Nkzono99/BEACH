title: 粒子源と粒子境界

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# 粒子源と粒子境界

このページでは、粒子を計算領域へ生成する方法と、粒子がbox境界へ到達したときの処理を整理します。
外部プラズマのPoisson modelそのものは[外部プラズマモデル](OuterPlasmaModels.html)を参照してください。

## 粒子源

| `source_mode` | 用途 | 生成位置 |
| --- | --- | --- |
| `volume_seed` | 初期粒子、簡単な軌道試験 | `pos_low`〜`pos_high` |
| `reservoir_face` | box外のreservoirからの流入 | 指定したbox面 |
| `photo_raycast` | 光照射された表面からの放出 | rayが最初に命中した要素 |

### volume seed

`npcls_per_step`個を指定領域へ生成し、drift付きMaxwell分布から速度を与えます。物理fluxから個数を
計算するsourceではないため、定量的な連続流入には`reservoir_face`を使います。

### reservoir face

面積$A$、流入flux$\Gamma_\mathrm{in}$、`batch_duration`、macro粒子重み$w$から期待粒子数を求めます。

$$
N_\mathrm{macro,expected}=\frac{\Gamma_\mathrm{in}A\,\mathrm{batch\_duration}}{w}
$$

端数はspeciesごとの`macro_residual`として次batchへ持ち越します。MPIでもglobal個数と残差を一度だけ決めてから
rankへ分配するため、world sizeで期待流入量が変わらない契約です。

速度はflux-weighted Maxwell分布またはvelocity gridからsampleします。上流から注入面までに電位障壁がある場合、
到達可能な法線速度を選別し、エネルギー保存で面上速度へ写像します。

### photo raycast

注入面からrayを飛ばし、最初に命中した表面要素から粒子を放出します。ray数、投影面積、電流密度、
`batch_duration`からmacro粒子重みを決めます。放出元への逆符号電荷は
[表面帯電モデル](SurfaceModels.html#光電子放出の電荷ledger)を参照してください。

## box境界

| 境界 | 粒子処理 |
| --- | --- |
| `open` | 外部モデルへ渡すか、escapeとして除去 |
| `reflect` | 面法線速度を反転して残り時間を追跡 |
| `periodic` | 反対側へwrapして残り時間を追跡 |

mesh衝突とbox面が同じstepに現れる場合は、軌道上で早いeventを採用します。数値的なevent処理は
[粒子追跡と衝突](ParticleTrackingCollision.html)を参照してください。

## open境界の処理

open面を通過した粒子の処理は、単純なescapeと外部領域へのtransferを区別します。

- `open_boundary_model="escape"`: その場でescapeとして集計
- legacy potential barrier: 通過点と無限遠の電位差からescape/reflectを判定
- 1D instant return: outer profileの保存エネルギーからescape/returnを判定
- 3D explicit orbit: outer領域の3D場で軌道を追跡

return modelは粒子sourceと独立ではありません。同じ電位差やcutoffをreservoir側とreturn側で二重適用しないよう、
非対応の組み合わせはvalidationで停止します。

## 使い分け

- 小さな決定論的試験: `volume_seed`
- 定常plasma流入: `reservoir_face`
- 表面光電子放出: `photo_raycast`
- 単純な有限box: `open_boundary_model="escape"`
- 自己整合outer sheath: 対応するouter modelとreturn modelを組み合わせる

モデル選択は[外部プラズマモデル](OuterPlasmaModels.html)、全入力keyは
[設定パラメータ](Parameters.html)を参照してください。詳細式は旧統合文書
[外部シースとreservoir粒子境界](SheathReservoirBoundary.html)に残しています。
