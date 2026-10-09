title: 光電子の放出と電荷の閉じ方

Lang: [日本語](PhotoelectronEmission.md) | [English](PhotoelectronEmission.en.md)

# 光電子の放出と電荷の閉じ方

光電子は、照射された表面から放出され、セル内の場で表面へ戻るか、box の外へ出ていきます。
BEACH は照射の光線を追跡して放出位置を決め（`source.mode="photo_raycast"`）、放出した光電子を他の粒子と同じ場・同じ方法で追跡します。
このページは、放出量と放出速度の決め方、放出元に残す反作用電荷、光電子の電荷の閉じ方（閉じた光電子、固定電流）を説明します。

## 放出から吸収までの流れ

1. box の面に置いた照射の開口から光線を出す。
2. 光線が最初に当たる三角形を探す。
3. 当たった要素から、照射側へ光電子を放出する。
4. 放出元の要素に反作用電荷を記録する。
5. 光電子を他の粒子と同じく追跡する。box の面に達したら、その面の扱いに従う（[box 境界での粒子の扱い](ParticleEscapeReturn.html)）。
6. 放出と吸収の電荷を、バッチ末尾に表面へ確定反映する。

放出と再吸収は同じバッチの中で起こり得ますが、場はバッチの途中で変えません。正味の表面電荷は次のバッチの場から効きます。

## 照射で放出面を決める

`source.inject_face` と `source.pos_low` / `source.pos_high` が、box の面上の矩形の開口を決めます。`source.ray_direction` は
開口から box の内側へ向く必要があり、省略すると面の内向き法線です。開口の面積を $A$、内向き法線を $\mathbf n_\mathrm{in}$、
光線の方向の単位ベクトルを $\hat{\mathbf d}$ とすると、光線に垂直な投影面積は

$$
A_\mathrm{proj}=A\left|\hat{\mathbf d}\cdot\mathbf n_\mathrm{in}\right|
$$

です。光線の始点は開口内で一様に選びます。光線は次の box の面までの区間ごとに最初に当たる三角形を探し、
非周期の面に達すると放出なしで終わり、周期の面に達すると反対側へ回り込んで続きます。周期の場では周期画像を含めて
最初の当たりを探し、放出位置をもとのセルへ戻します。`particles.raycast.max_bounces` 回の回り込みで当たらない光線は、粒子を作りません。

## 光線 1 本の放出電流

放出電流密度を $J_\mathrm{emit}>0$、実粒子の電荷を $q$、全 MPI rank の光線数の合計を $N_\mathrm{ray}$、バッチ幅を
$\Delta t_\mathrm{batch}$ とすると、当たった光線が作るマクロ粒子の重みは

$$
w_\mathrm{hit}
=\frac{J_\mathrm{emit}A_\mathrm{proj}\,\Delta t_\mathrm{batch}}
{|q|N_\mathrm{ray}}
$$

です。外れた光線は粒子を作らないので、影や見かけの面積は当たる割合として放出量に入ります。
`sampling.rays_per_batch` は放出量ではなく標本数で、増やすと $w_\mathrm{hit}$ が小さくなり、放出位置の統計誤差が減ります。
`sampling.weight` と `sampling.target_macro_particles_per_batch` は使いません。

## 放出速度

当たった三角形の法線のうち、入射した光線と反対側を向くものを放出の法線 $\mathbf n_s$ とします。放出位置は、
同じ要素への直後の再衝突を避けるため、当たった点から $\mathbf n_s$ 方向へ $10^{-12}$ m ずらします。
温度から $\sigma=\sqrt{k_\mathrm{B}T/m}$ を求め、局所的な基底 $(\mathbf n_s,\mathbf t_1,\mathbf t_2)$ で速度を選びます。

- 法線速度: drift `source.normal_drift_speed_m_s` を持つ、束で重み付けした半 Maxwell 分布。
- 接線の 2 成分: 平均 0、標準偏差 $\sigma$ の Gauss 分布（$6\sigma$ で打ち切り）。

## 反作用電荷

`source.deposit_opposite_charge_on_emit=true` なら、放出元の要素 $i$ に

$$
\Delta q_{i,\mathrm{emit}}=-q w
$$

を加えます。電子では $q<0$ なので、表面には正の電荷が残ります。光電子が要素 $j$ に吸収されると、通常の吸収として
$+qw$ を $j$ に加えます。同じ要素に戻れば相殺し、別の要素に戻れば表面内で電荷が移ります。

## 光電子の電荷の閉じ方

セルの高さは Debye 長より小さいので、上端面に達した光電子の多くは、実際には外部のシースで押し戻されます。
上端面を開放のままにすると、光電子の脱出を過大に数えます。これを補う方法を、粒子種の `charging.closure` で選びます。

| `charging.closure` | 光電子の扱い | 使う場面 |
|---|---|---|
| 省略 | 追跡のまま。上端面の扱いだけで決まる | 比較の基準 |
| `neutral_return` | 上端面で反射し、光電子の正味の電流を 0 に合わせる | 表面内の再分配を調べる |
| `fixed_current` | 放出と帰還の電流を、外から与えた目標に合わせる | 外部の電流モデルと組み合わせる |

周期表面のケースで、どれを選ぶかは[プラズマ中の周期表面を設定する](PeriodicPlasmaSurface.html)で比べています。

### 閉じた光電子（`neutral_return`）

光電子の粒子種で、照射の面（`source.inject_face`）を反射にし、`neutral_return` を指定します。

```toml
[particles.species.charging]
closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"
```

反射は法線速度だけを反転し、接線速度と位置を保ちます。`"redistributed_reflect"` にすると、速度は同じに反転し、
戻る位置だけを面内で一様に選び直します（[Zimmerman et al. (2016)](https://doi.org/10.1002/2016JE005049) の上端での
光電子の帰還の扱いと同じ考え方です）。

反射しても、1 粒子あたりの step 上限までに戻らない光電子が残ります。`neutral_return` は、1 バッチの光電子の放出電荷
$S<0$ と、表面に吸収された（帰還した）電荷 $R<0$ を全 MPI rank で合計し、帰還先に置く電荷を $S/R$ 倍します。

$$
(-S)+\frac{S}{R}R=0
$$

放出元の反作用電荷と合わせて、光電子による表面の総電荷の変化はちょうど 0 になります。未帰還の光電子は、
同じバッチで帰還した光電子と同じ分布で戻ると近似しています。総電荷を 0 にしても、異なる高さの面へ電荷が移ると、
面平均した鉛直方向の双極子は残ります。

次の場合はバッチを受理せずに停止します。

- 放出があるのに帰還がない。
- 光電子が開放面から脱出した、または複数の box 面を同時に横切る粒子として捨てられた（soft discard）。
- 値が有限でない、または符号が合わない。
- 未帰還の割合が 5% を超えた。

### 固定電流（`fixed_current`）

放出と帰還の電流を、外部のモデルで決めた目標に合わせます。目標は表面帯電への寄与の符号付きで、別々に与えます。

```toml
[particles.species.charging]
closure = "fixed_current"
target_emission_current_a = 4.5e-15
target_absorbed_current_a = -3.7e-15
```

BEACH は、放出元の分布と帰還先の分布をそれぞれ一律の倍率で目標に合わせます。二つの大きな電流の差である正味の電流は
倍率の分母に使わないので、放出と帰還がほとんど打ち消し合う場合にも安定です。上端面は開放にし、
`neutral_return` と併用しません。外部シースとの接続（[Zhao 定常シース](ZhaoStationaryClosure.html)）を使うと、
目標は外部シースの零電流根から自動で決まります。

倍率が安定でも、要素ごとの分布の統計精度は別の問題です。帰還が 1 件しかなければ、帰還の目標の全量がその要素に置かれます。

## 確認する出力

| 出力 | 見るもの |
|---|---|
| `charge_ledger.csv` | 放出・吸収・脱出の電荷と個数。`neutral_return` では `neutral_return_weight_scale` と `neutral_return_unresolved_fraction`、`fixed_current` では `fixed_*_weight_scale` |
| `fixed_current_history.csv` | `fixed_current` の倍率の時間変化（`target_over_tracked`） |

倍率が 1 から大きく離れる計算では、追跡ではなく閉じ方の仮定が電荷の分布を決めています。

## 収束を確かめる

- `sampling.rays_per_batch` を増やし、当たる割合、放出電流、帯電の分布が変わらないことを確かめる。
- 帰還の位置を評価するなら、`particles.tracking.dt_s` を小さくし、`max_steps_per_particle` を増やして変わらないことを確かめる。
- `neutral_return` では、倍率が 1 に近く、未帰還の割合が小さくなるまで、step 上限と box の高さを見直す。

## 適用範囲

- 放出速度は、設定した表面の半 Maxwell 分布です。表面の材質や仕事関数の分布は扱いません。
- 上端面での反射（`neutral_return`）は有限の box の上に鏡を置く近似で、外部のシースや準中性を解きません。
- 絶縁体の表面では、帰還した電荷はその要素に残り、表面に沿って伝導しません。
