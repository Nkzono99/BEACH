title: 10 分チュートリアル

Lang: [日本語](Tutorial.md) | [English](Tutorial.en.md)

# 10 分チュートリアル

このチュートリアルでは、各 batch に 200 個、合計 4000 個のマクロ電子を 20 batch にわたって
絶縁体平面へ入射します。
吸収位置に電荷分布が形成され、その電荷が作る場によって後続 batch の電子が反射され始めるところまでを確認します。

これは BEACH の batch 間フィードバックを見せる教材です。実在環境を再現した研究ケースではなく、
定常状態や物理時間への収束も主張しません。

## 0. インストールを確認する

```bash
beach --version
beachx --help
```

両方のコマンドが実行できれば準備完了です。コマンドが見つからない場合は、先に
[インストール](Installation.html)を確認してください。HPC のログインノードでは直接実行せず、
計算ノードの割当内で以降のコマンドを実行します。

## 1. 設定を作る

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
```

生成される設定値と挙動はリポジトリの
[`examples/tutorial_insulator.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/tutorial_insulator.toml)
と同一です。コメントや配列の改行など、TOML の表記は異なる場合があります。設定全文を暗記する必要は
ありません。まず次の対応だけを確認します。

| 設定 | 値 | このケースでの役割 |
| --- | --- | --- |
| `sim.batch_count` | `20` | 電荷差分の確定反映（commit）と場の再計算を 20 回繰り返す |
| `sim.dt` / `sim.max_step` | `5e-8 s` / `80` | 各電子の軌道を時間積分する |
| `sim.field_solver` / `field_boundary.mode` | `direct` / `free` | 有限平面の自由空間場を Direct BEM で評価する |
| `npcls_per_step` | `200` | 各 batch に 200 個のマクロ電子を生成する |
| `w_particle` | `2e5` | 1 マクロ電子を実電子 200,000 個分として、場の変化を短い実行で見せる |
| `pos_low` / `pos_high` | z = `0.8 m` の矩形 | 平面上方の広い領域から入射位置を標本化する |
| `drift_velocity` | `[0, 0, -1e6] m/s` | 電子を平面へ向ける |
| `mesh.templates` | 12 × 12 cells | 平面を 288 個の三角形要素へ分割する |
| `history_stride` / `write_potential_history` | `1` / `true` | 各 batch の電荷・電位分布を保存する |

1 batch 内の全電子は、batch 開始時の同じ電場を使います。吸収電荷は batch 末尾に一度だけ
確定反映され、次の batch の電場へ反映されます。この順序が BEACH の中心です。

## 2. 実行する

参照結果と同じ乱数列にするため、このチュートリアルでは 1 OpenMP thread を使います。

```bash
OMP_NUM_THREADS=1 beach beach.toml
beachx inspect outputs/tutorial
```

正常終了すると、次の値が表示されます。

```text
processed_particles=4000
absorbed=3720 escaped=280
batches=20
charge_sum=-1.192019e-10
potential_min=-4.671330e+00
potential_max=-2.579807e+00
```

`processed_particles=200 × 20` と `batches=20` は設定から必ず決まります。上の吸収・escape 件数と
電位範囲は、`rng_seed=12345`、1 thread、現行版での参照値です。異なる thread 数や将来の乱数実装では、
分布の細部が変わる場合があります。

## 3. なぜ後半の電子が escape するか

最初の 13 batch では、各 batch の 200 電子がすべて平面へ吸収されます。負電荷が蓄積すると
表面電位はさらに負になり、負の電子には入射を妨げる電場が生じます。14 batch 目以降は一部が
反射され、上側の open 境界から escape します。

最終表面電荷は

$$
3720\times 2\times10^5\times(-e)
=-1.192019\times10^{-10}\ \mathrm{C}
$$

であり、`charge_sum` と一致します。`charges.csv` の `charge_C` は各三角形の総電荷 [C]、
図の色はそれを要素面積で割った表面電荷密度 [C/m²] です。

この free-space の有限表面では Coulomb 電位を遠方で 0 V とする規約を使います。matching plane や
`reservoir.phi_infty` は接続していないため、この例の電位を外部プラズマ電位と読み替えないでください。

## 4. 電荷履歴と最終電位を見る

```bash
beachx inspect outputs/tutorial \
  --save-mesh outputs/tutorial/charge.png

beachx animate outputs/tutorial \
  --quantity charge \
  --save-gif outputs/tutorial/charge_history.gif \
  --total-frames 20
```

![20 batch 後の表面電荷密度](media/tutorial_insulator_charge.png)

最終状態の表面電位は `mesh_potential.csv` に保存済みです。保存値を再計算せず描画するには、
`reference_point=None` を明示します。

```python
from beach import Beach

run = Beach("outputs/tutorial", config_path="beach.toml")
fig, _ = run.plot_potential(reference_point=None)
fig.savefig("outputs/tutorial/potential.png", dpi=150)
```

![20 batch 後の表面電位](media/tutorial_insulator_potential.png)

![batch ごとの表面電荷密度変化](images/tutorial_insulator_charge.gif)

任意点や各 batch の電位を再構成するには native field kernel が必要です。必要になった段階で
[後処理チュートリアル](PostprocessTutorial.html)の追加手順へ進んでください。

## 5. この結果が示すこと・示さないこと

この実行で確認できるのは、粒子生成、軌道積分、衝突と escape、要素別電荷の確定反映、場の再計算、
履歴出力が batch 間で結合していることです。

一方、`w_particle` は後続 batch への影響を短時間で見せる教材値で、`batch_duration=0` の `volume_seed` には
物理的な秒単位の流入率を割り当てていません。`last_rel_change` と `tol_rel` も自動停止条件ではありません。
したがって、この 20 batch を定常帯電や特定環境の 20 秒として解釈しません。

## 6. 次へ進む

| 目的 | 次に読むページ |
| --- | --- |
| 各出力ファイルと保存則を確認する | [出力ファイルを調べる](OutputGuide.html) |
| 20 batch の checkpoint から一度再開する | [実行・再開する](Execution.html#checkpointから一度再開する) |
| 粒子数、形状、粒子源を変更する | [シミュレーションケースを設計する](ConfigurationRecipes.html) |
| batch、粒子、電荷の確定反映を理解する | [BEACH が解く範囲と計算の流れ](Algorithms.html) |
| 研究結果として受理できるか調べる | [計算結果の妥当性確認](ValidationGuide.html) |
| 実行に失敗した | [トラブルシューティング](Troubleshooting.html) |
