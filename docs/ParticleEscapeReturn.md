title: box 境界での粒子の扱い

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# box 境界での粒子の扱い

粒子が box の面に達したとき、周期の面では反対側へ回り込み、非周期の面では面ごとに選んだ扱いに従います。
この扱いは、粒子をどの粒子源で作ったかによらず共通です。BEACH は box の外のプラズマを解かないので、
面の外へ出た粒子が戻ってくるかどうかは、ここで選ぶ扱いが代わりに決めます。

## 面ごとの扱いを選ぶ

| 扱い | 設定 | 粒子に起きること |
|---|---|---|
| 周期 | `domain.periodic_axes` | 反対側の面へ移り、速度は変えない |
| 開放 | `"open"` | 脱出する。`open_model="potential_barrier"` なら電位障壁で反射か脱出かを判定する |
| 反射 | `"reflect"` | 法線速度を反転する。位置と接線速度は変えない |
| 再配置つき反射 | `"redistributed_reflect"` | 法線速度を反転し、面内の位置を一様に選び直す |

```toml
[particles.boundary]
z_low = "open"
z_high = "open"
open_model = "escape"      # 開放面の扱い。既定値
```

周期の軸は場と粒子に共通の性質なので、`[particles.boundary]` でも粒子種ごとの設定でも変えられません。
省略した非周期の面は開放です。粒子種ごとに扱いを変えるには `[particles.species.boundary]` で
`inherit`（既定値）、`open`、`reflect`、`redistributed_reflect` を選びます。

## 脱出

開放面を横切った粒子は、横切った位置で取り除きます。粒子の電荷は粒子種ごとの脱出電荷
（`charge_ledger.csv` の `escaped_to_infinity_C`）に数え、表面電荷は変えません。

## 電位障壁（`potential_barrier`）

外部プラズマの上流電位 $\phi_\infty$（`particles.reservoir.phi_infty_v`）を決め、開放面から出ていく粒子が
上流まで行けるかを、横切った点の電位 $\phi_b$ で判定します。

$$
\Delta U=q(\phi_\infty-\phi_b),\qquad
K_n=\frac12 m v_n^2
$$

$v_n>0$ は外向きの法線速度です。$\Delta U>0$ かつ $K_n<\Delta U$ なら法線速度を反転して追跡を続け、それ以外は脱出とします。
接線速度は変えません。

```toml
[particles.boundary]
open_model = "potential_barrier"

[particles.reservoir]
phi_infty_v = 0.0
```

$\phi_b$ は、バッチの初めに固定した場（外部電場を含む）で評価します。一様な外部電場には無限遠の電位が決まらないので、
外部電場と併用するときは $\phi_\infty$ をその基準に合わせて決めます。複数の開放面を同時に横切る角では判定が決まらず、停止します。

流入する粒子を同じ上流電位でふるい分けるには、流入の写像を `infinity_barrier` にします（[境界から粒子を流入させる](ReservoirInjection.html#3-流入の写像を選ぶ)）。

## 外部シースの障壁

外部シースとの接続を使うと、上端面を出ていく電子と光電子は、外部シースの零電流根が与える障壁で反射か脱出かを判定します。
判定の式は電位障壁と同じで、障壁の電位と戻す位置の選び方が異なります（[外部シースとの接続](ZhaoStationaryClosure.html#上端面での粒子の出入り)）。

## 反射と再配置つき反射

どちらも法線速度だけを反転し、接線速度を保ちます。`reflect` は位置を変えません。`redistributed_reflect` は、
1 つの面での反射では面内の 2 軸の位置を、両端の余白を除く範囲から一様に選び直します。角や辺で複数の面に同時に達したときは、
達した面に含まれない軸だけを選び直します。

上端面での光電子の反射は、閉じた光電子（[光電子の放出と電荷の閉じ方](PhotoelectronEmission.html#閉じた光電子neutral_return)）で使います。
同時に複数の面に達したときの判定順は[粒子の衝突・境界イベント](ParticleEvents.html)にあります。

## 確認する出力

| 出力 | 見るもの |
|---|---|
| `summary.txt` | `escaped_boundary`（境界から出た粒子数）、`particle_ordinary_open_model`（使われた開放面の扱い） |
| `charge_ledger.csv` | 粒子種ごとの `escaped_count`、`escaped_to_infinity_C` |

`summary.txt` の `escaped` には、step 上限までに行き先が決まらなかった粒子（`survived_max_step`）も含まれます。
実際に境界から出た数は `escaped_boundary` で読みます。

## 適用範囲

- 開放面の外の電場、折り返す位置、飛行時間、空間電荷は解きません。電位障壁は、上流までの電位差だけで判定する近似です。
- 反射は、box の面に置いた鏡です。外部のプラズマやシースを表すものではありません。
