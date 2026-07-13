title: 粒子追跡と衝突

Lang: [日本語](ParticleTrackingCollision.md) | [English](ParticleTrackingCollision.en.md)

# 粒子追跡と衝突

BEACHは粒子の位置と速度を同じ時刻で保持し、Boris速度更新、台形位置更新、box／mesh event判定を
1 stepとして実行します。

## 同時刻Boris更新

入力を$(\mathbf{x}^n,\mathbf{v}^n)$とし、予測中点で電場を1回評価します。

$$
\mathbf{x}_\mathrm{mid}=\mathbf{x}^n+\frac12\mathbf{v}^n\Delta t,
\qquad
\mathbf{E}_\mathrm{mid}=\mathbf{E}(\mathbf{x}_\mathrm{mid})+\mathbf{E}_0
$$

この場と一様磁場`sim.b0`でBoris回転を行い、$\mathbf{v}^{n+1}$を得ます。位置は台形則で更新します。

$$
\mathbf{x}^{n+1}=\mathbf{x}^n+
\frac12(\mathbf{v}^n+\mathbf{v}^{n+1})\Delta t
$$

public stateはhalf-step速度を保持しません。smoothな規定場では位置・速度とも二次精度、pure Bでは
速度normを丸め誤差まで保存するcontractです。

## eventの順序

更新候補の線分／二次軌道について、三角形衝突とbox面通過の最初のeventを比較します。

| 最初のevent | 処理 |
| --- | --- |
| mesh hit | 表面で吸収し、候補終点は採用しない |
| open face | escapeまたはouter modelへtransfer |
| reflect face | 法線速度を反転し、残り時間を再積分 |
| periodic face | 反対側へwrapし、残り時間を再積分 |

1 step中のbox event数には安全上限があり、有限回でevent処理が進まない場合はfail closedで停止します。

## 三角形衝突

要素数が小さいmeshでは全要素を調べ、大きいmeshではAABB一様gridと3D-DDAで候補を絞ります。
候補三角形との交差はMöller–Trumbore法で判定し、線分parameterが最小の要素を最初の命中として採用します。

退化三角形、非有限座標、停止したDDA、過大なperiodic image範囲は「命中なし」として続行せず、statusを
集約して停止します。

## periodic2 collision

meshはprimary cellの要素だけを保持します。粒子線分が交差し得る周期image範囲を計算し、線分を各imageの
逆shiftでbase meshへ写して衝突を調べます。命中位置は物理image座標とprimary cellへwrapした座標の両方を
保持します。

場の周期画像と衝突の周期画像は目的が異なります。場の計算方法は
[periodic2場計算](PeriodicElectrostatics.html)を参照してください。

## 未解決粒子

粒子が`sim.max_step`へ達しても吸収やescapeへ分類されなかった場合、`survived_max_step`として記録します。
この粒子を暗黙に吸収・escape扱いしないことが電荷ledgerのcontractです。

## 数値確認

- `sim.dt`を半分にして軌道・吸収位置・統計を比較する
- mesh refinementで最初の衝突が安定することを確認する
- `survived_max_step`が結論へ影響していないことを確認する
- periodic seamを跨ぐ軌道を別途確認する

完全な式とstatus一覧は旧統合文書
[粒子追跡と表面電荷蓄積](ParticleChargeLoop.html#8-粒子時間積分-boris速度更新と同時刻状態)に残しています。
