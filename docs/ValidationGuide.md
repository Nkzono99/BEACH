title: 計算結果の妥当性確認

Lang: [日本語](ValidationGuide.md) | [English](ValidationGuide.en.md)

# 計算結果の妥当性確認

正常終了しただけでは、計算結果が数値的・物理的に妥当とは限りません。次の3つを順番に確認します。

1. 設定した計算が最後まで実行され、粒子と電荷の収支を説明できる。
2. 時間刻みやmesh解像度を変えても、注目する量が十分に安定している。
3. 採用した物理モデルの前提を満たし、結果から述べる結論がその適用範囲を超えていない。

各値がどのファイルにあるかは[出力ファイルを調べる](OutputGuide.html)にまとめています。
このページでは、その値を使って計算を受理できるか判定します。

## 先に判定基準を決める

収束確認を始める前に、何をもって「結果が変わらない」とするかを決めます。

- 比較する主要量: 表面の総電荷・電位分布、吸収・流出flux、または研究上の結論に直接使う量
- 許容差: 主要量に許容する絶対差または相対差
- 評価区間: 最終値、定常化後の時間平均、最大値など
- 統計誤差の扱い: 乱数seedを変えた計算、または同じ条件のensembleから求めるばらつき

判定基準を結果を見た後で都合よく変えないことが重要です。BEACHに全ケース共通の収束許容差はないため、
計算の目的に応じて基準を明示します。

## 1. 実行が完了したことを確認する

```bash
beachx inspect outputs/latest
```

少なくとも次を確認します。

| 確認先 | 判定 |
| --- | --- |
| processの終了code | `0` |
| `summary.txt` | `batches == sim.batch_count` |
| `charges.csv` | 最終的な要素電荷が保存されている |
| 必要な履歴 | 設定したstrideで欠落なく保存されている |
| restartした計算 | model・mesh・species fingerprintが一致している |

ここで確認できるのは、指定した計算が完了したことだけです。物理量の収束は以降で確認します。

## 2. 粒子数と電荷の収支を確認する

`summary.txt`の粒子数を、設定したsourceと境界条件から説明できるか確認します。

| 項目 | 見方 |
| --- | --- |
| `absorbed` | 表面に到達して吸収された粒子 |
| `escaped_boundary` | box境界またはouter modelから流出した粒子 |
| `survived_max_step` | `sim.max_step`までに吸収・流出が確定しなかった粒子 |

`survived_max_step`は吸収にも流出にも数えられない未解決粒子です。結論に影響する量なら、
`sim.max_step`、`sim.dt`、boxの大きさを見直し、十分に小さくなることを確認します。

次に`charge_ledger.csv`と`summary.txt`の電荷収支を確認します。

- `charge_ledger_residual_C`: 注入、放出、表面吸収、無限遠への流出を含む保存残差
- `charge_ledger_discarded_unresolved_abs_C`: 粒子種別間で相殺しない未解決電荷の絶対値和

保存残差が小さくても、未解決電荷が大きければ計算を受理できません。残差と未解決量は必ず別々に評価します。

## 3. 時間発展と統計的なばらつきを確認する

`charge_history.csv`と、必要なら`potential_history.csv`を使い、定常化したように見えるかを確認します。
最終点だけでなく、評価区間について次を見ます。

- 系統的な増加・減少が残っていないか
- batchごとの振動や突発的な変化がないか
- 時間平均に対するばらつきが、先に決めた許容差より小さいか
- `sim.batch_count`を増やしても平均値や結論が変わらないか

`last_rel_change`と`sim.tol_rel`は監視用です。現行実装では早期停止条件ではなく、
`last_rel_change < tol_rel`だけで定常状態とは判定しません。

Monte Carloノイズと時間離散化の影響は分けて確認します。

| 変更するもの | 主に確認できる影響 |
| --- | --- |
| `sim.rng_seed` | 乱数に起因するばらつき |
| マクロ粒子数または粒子重み | Monte Carloノイズ |
| `sim.batch_count` | 計算時間が十分か |
| `sim.batch_duration`または`sim.batch_duration_step` | batch末尾の表面電荷更新の安定性 |

`batch_duration`は基準値の0.5倍と2倍を目安に比較します。詳しい考え方は
[`batch_duration`の安定性と定常値](BatchDurationStability.html)にまとめています。

## 4. 数値解像度への依存性を確認する

基準ケースを複製し、原則として一度に1つの設定だけを変えます。各ケースで同じ主要量を比較し、
変更による差が先に決めた許容差より小さくなることを確認します。

| 収束軸 | 比較例 | 確認する誤差 |
| --- | --- | --- |
| 粒子時間刻み | `sim.dt`を1/2にする | 軌道、衝突位置、吸収・流出量 |
| 粒子追跡長 | `sim.max_step`を増やす | 未解決粒子による打切り |
| 表面mesh | 三角形を細分化する | 電荷・電位の空間離散化 |
| 場ソルバ | 小規模ケースでDirectと比較する | Treecode/FMM近似 |
| 有限周期画像 | `field_periodic_image_layers`を増やす | 画像和の打切り |
| outer model | grid点数や表面samplingを増やす | 外部profileとgeometry sampling |

meshを固定したまま経路積分や画像層だけを収束させても、mesh離散化誤差まで評価したことにはなりません。
どの収束軸を確認し、どれを未評価のまま残したかを区別します。

## 5. 選んだ物理モデルの診断を確認する

共通の収束確認に加え、使用したモデルに対応する診断を確認します。

| 構成 | 追加で確認するもの |
| --- | --- |
| 有限画像の`periodic2` | `field_periodic_image_layers`依存性。有限画像の結果を無限周期へ外挿しない |
| `cached_kneq0` | cache fingerprint、cold/warm結果の一致、zero mode、Gauss residual |
| `infinity_barrier` | reservoir面の平均電位、画像層依存性、面内電位ばらつきの警告 |
| `potential_barrier` | open面の通過点電位と`sim.phi_infty`、法線運動エネルギーによる反射・流出判定 |
| `kinetic_1d` | solver status、Poisson residual、単調分枝、Bohm条件 |
| `unified_linear_response` | accessibility refinement、線形性、Gauss則の整合性、outer軌道のenergy/frozen-field error |
| 光電子 | 放出・帰還・流出の電荷収支、ambient charge ratio、energy histogramの範囲 |

各診断の定義と適用範囲は、[有限画像構成](FinitePeriodicConfiguration.html)、
[外部プラズマモデル](OuterPlasmaModels.html)、[粒子のescapeとreturn](ParticleEscapeReturn.html)にあります。

## 6. 物理的な結論の範囲を確認する

最後に、数値的に安定した結果が目的の物理的主張を支えているかを確認します。

- 比較するケース間で、意図したモデル以外に変更した入力を列挙する。
- 境界条件、粒子source、表面モデル、場の評価方法を明記する。
- 有限boxの結果を無限遠、有限時間の結果を定常状態、有限画像和を無限周期として扱わない。
- 数値誤差、Monte Carloのばらつき、モデルの適用限界を結論の有効桁や誤差幅に反映する。

最低限、使用した設定ファイル、出力directory、主要量と判定基準、各収束ケースとの差を残します。
コード自体のリリース判定に使う小規模fixtureとHPC検証は、ユーザーケースの妥当性確認とは別に
[物理リリースの検証](PhysicsReleaseVerification.html)で扱います。
