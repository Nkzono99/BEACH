title: 粒子源の全体像

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# 粒子源の全体像

このページは、粒子をどこから生成するかに応じて `source_mode` を選び、生成後に共通して適用される処理を確認するための
入口です。各粒子源の数値計算や物理モデルは、対応する詳細ページで説明します。

## 目的に応じて粒子源を選ぶ

各 `[[particles.species]]` に 1 つの `source_mode` を設定します。

| `source_mode` | 粒子数を決める量 | 生成位置 | 適した用途 |
| --- | --- | --- | --- |
| `volume_seed` | `npcls_per_step` | `pos_low` から `pos_high` の範囲 | 初期粒子、軌道試験、指定個数の生成 |
| `reservoir_face` | 流入流束、開口面積、`batch_duration` | 指定したボックス面 | 外部リザーバーからの連続流入 |
| `photo_raycast` | 電流密度、投影面積、`batch_duration`、光線数 | 光線が最初に命中した表面 | 光照射による表面放出 |

粒子数を直接指定するなら `volume_seed`、外部プラズマから面を横切る流入を与えるなら `reservoir_face`、
照射された表面から放出するなら `photo_raycast` を選びます。

## `volume_seed` で指定個数の粒子を作る

`volume_seed` は、各バッチに `npcls_per_step` 個の粒子を生成します。位置は直方体
`[pos_low, pos_high]` 内の一様分布、速度は

$$
\mathbf v=\mathbf u+\sigma\mathbf Z,
\qquad
\sigma=\sqrt{\frac{k_\mathrm{B}T}{m}}
$$

というドリフト付き Maxwell 分布です。`thermal_speed` を指定した場合は、温度から求めた $\sigma$ より優先します。
標準正規変量は各成分で $6\sigma$ に切られます。

この方式は粒子数を直接指定するため、物理的な面流束から粒子数を計算しません。

## `reservoir_face` で外部からの連続流入を作る

`reservoir_face` は、ボックス外に与えた上流分布関数（VDF）を、指定面から内向きに流入する粒子へ変換します。
各バッチの粒子数は物理流束から決まり、法線速度は面を横切る粒子に適した流束重み付き分布から生成されます。

この方式を選んだ後に必要となる計算は、次のページへ分けています。

- 上流 VDF から流入量、マクロ粒子数、注入開口からの初期位置、面到達速度を作る処理:
  [`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)
- `infinity_barrier`、外部プラズマ、流出境界を含む構成の選択:
  [境界・外部領域の構成を選ぶ](OuterPlasmaModels.html)
`reservoir_face` 自体は、外部シースやボックス外の粒子軌道を解くモデルではありません。

## `photo_raycast` で照射面から粒子を放出する

`photo_raycast` は、ボックス面上の照射開口から光線を発射し、ボックス境界条件に従って進めます。
ボックス内で最初に命中した要素から粒子を放出し、要素法線に対する流束重み付き Maxwell 分布から速度を生成します。

放出元へ置く逆符号電荷と、通常の粒子追跡・外部シースによるescape/returnは
[光電子の放出とライフサイクル](PhotoelectronEmission.html)で説明します。

## 生成後は同じ粒子追跡へ入る

粒子源が異なっても、生成後に保持する状態と追跡処理は共通です。主な状態は位置 $\mathbf x$、速度 $\mathbf v$、
実粒子 1 個の電荷 $q$ と質量 $m$、マクロ粒子重み $w$、粒子種 ID です。追跡するマクロ粒子 1 個の電荷は $qw$ です。

| バッチ内の結果 | 処理 |
| --- | --- |
| メッシュへ吸収 | 命中要素へ $qw$ を堆積 |
| 開放面から無限遠へ脱出 | 粒子を除去し、粒子種別の脱出量へ計上 |
| 外部領域から帰還 | 同じ粒子を界面へ戻し、残りのステップを再積分 |
| `max_step` まで生存 | 未解決粒子としてバッチ末尾で破棄・計上 |

粒子の前進は [Boris 粒子更新](BorisPusher.html)、メッシュ衝突とボックス境界の順序は
[粒子の衝突・境界イベント](ParticleEvents.html)、ボックス外の処理は
[粒子の脱出と帰還](ParticleEscapeReturn.html)で説明します。

## 粒子源はバッチ開始時の場を使う

粒子源は、バッチ開始時に電場・電位と外部プラズマ状態を更新した後で評価されます。したがって、リザーバーの速度補正と
光電子の脱出率は、前のバッチまでに確定した表面電荷を参照します。生成された粒子は、そのバッチ内で固定された場の中を進みます。

粒子生成時に表面電荷を変えるのは、`photo_raycast` が放出元へ逆符号電荷を置く場合だけです。この差分も吸収電荷とともに
バッチ末尾で確定され、現在のバッチの場は変更しません。全体の順序は
[計算モデルの全体像](Algorithms.html)を参照してください。

## MPI と再開で生成量を保つ

`reservoir_face` の全体生成数とマクロ粒子端数はルートランクで一度だけ決め、各ランクへ分配します。
`photo_raycast` の `rays_per_batch` も全ランクの合計です。このため、期待される流入量や放出量は MPI ランク数に依存しません。
一方で、乱数列と個々の粒子軌道はランク数によって変わる場合があります。

リザーバーの端数は粒子種ごとの `macro_residual` として `macro_residuals.csv` に保存され、再開時に復元されます。
詳しい端数処理と確認項目は [`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)にあります。

## Code reference

- 粒子分布と光線追跡: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- 粒子源ごとのバッチ生成とマクロ粒子端数: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- 粒子源入力の検証: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- 電荷収支とバッチ追跡: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- マクロ粒子端数のチェックポイント: [`bem_restart.f90`](../src/runtime/bem_restart.f90)
