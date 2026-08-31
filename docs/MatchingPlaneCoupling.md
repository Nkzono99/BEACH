title: matching-plane で外部シースを接続する

Lang: [日本語](MatchingPlaneCoupling.md) | [English](MatchingPlaneCoupling.en.md)

# matching-plane で外部シースを接続する

`surface_current_model.model="matching_plane_quasistatic"` は、BEACH の box 上端を外部 1D シースとの
matching plane（整合面）にします。表面電荷が変わるたびに、外部応答から上端電位、ambient 粒子の流入束、
光電子（PE）の return 障壁を更新できます。

これは外部シースの時間発展を解く機能ではありません。BEACH 内の 3D 軌道計算と、外部の準定常な低次元応答を
1 accepted batch ごとに整合させる境界 closure です。全入力キーは[入力パラメータ](Parameters.html#surface_current_model-外部シースclosure)、
出力名は[出力形式リファレンス](OutputReference.html#matching_plane_quasistatic)を正本とします。

このページは、通常の利用判断から設定、実行、診断までを扱います。応答 CSV、`implicit_zero_mode`、固定点の式は
[matching-plane 数値・応答表リファレンス](MatchingPlaneReference.html)で検索できます。

## 1. この model を使うか決める

外部プラズマをどう表現したいかで model を選びます。

| 必要な表現 | 選ぶ model |
|---|---|
| 表面電荷に応じて外部応答も batch ごとに変えたい | `matching_plane_quasistatic` |
| 定常シースの零電流根から得た電流と障壁を run 中固定したい | [`zhao_stationary`](ZhaoStationaryClosure.html) |
| 外部シース closure を使わず、BEACH 内の場と粒子だけを扱う | `none` |

matching-plane では、外部応答の取得方法をさらに選びます。

| backend | 適する用途 | 主な制約 |
|---|---|---|
| `response_backend="table"`（既定） | 独立に検証した Zhao / 1D PIC sweep を監査可能な固定 snapshot として使う | query は安価。応答 CSV と有限な軸範囲が必要 |
| `response_backend="zhao_online"` | CSV を作らず組み込み Zhao を直接使う | `implicit_zero_mode` は選択 branch 内で終点を探索。平面・無衝突・非磁化で、nonlinear solve の費用がある |

`examples/matching_plane_response_synthetic.csv` は table 経路の smoke test 専用です。物理計算には使えません。

### PE の有無

| 構成 | BEACH が扱うもの |
|---|---|
| PE なし | ambient electron / ion の流入と外部電位 gauge。PE feedback、return、escape は 0 |
| PE あり | 上記に加え、整合面を出る PE の束・平均法線 energy、外部障壁による return / escape |

PE なしでも外部シースとの接続は有効です。`photoelectron_species` と `photo_raycast` species を省略するだけで、
整合面電位と ambient 流入束は引き続き response backend が決めます。

### stationary Zhao と online Zhao を区別する

名前は似ていますが、解いている条件が異なります。

| | `zhao_stationary` | matching-plane の `zhao_online` |
|---|---|---|
| 拘束 | 壁面での零電流 $J=0$ | 現在の $D_H/\epsilon_0$ を整合面電場として指定 |
| 更新 | run 開始時に 1 回 | 固定点反復の query ごと |
| PE なし | Type C | Type C には固定されない。$D_H=0$ は縮退した Type B |
| 電流 | species 別の固定 target | 応答束と粒子追跡で得た実測堆積 |

したがって、PE を省略したことだけを理由に online Zhao の結果を Type C と解釈してはいけません。
既定の `zhao_root_selection="require_unique"` では、`zhao_branch="auto"` は現在の $D_H$ に適合する
一意な物理解を探し、一意性を確認できなければ停止します。

## 2. 最小構成を作る

x / y 周期・z open の `periodic2` case を基にします。全 mesh 頂点を
`domain.box_max` の z 成分より厳密に下へ置き、`sim.e0` と `sim.b0` は 0 にします。ambient electron と ion は
z-high reservoir から流入させます。

`response_backend="zhao_online"` を PE ありで使う最小の model 選択は次のとおりです。

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.05, 0.0, 0.0]
```

`coupling_atol` の第 2 成分は PE 平均法線 energy の許容値 [eV] です。上の 0.05 eV は例であり、
ray 数と macro 粒子数を変えた収束試験から決めます。完全な case は
`examples/periodic2_matching_plane_zhao_online.toml` にあります。

検証済みの応答表を使う場合は model 部分を次のように変えます。

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "table"
response_table_path = "matching_plane_response.csv"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.0, 0.0, 0.0]
```

`response_table_path` は設定ファイルのディレクトリから解決します。完全な配線例は
`examples/periodic2_matching_plane_quasistatic.toml` です。

PE なしの online case は PE role を指定しません。

```toml
[periodic2]
lower_boundary_model = "e_bottom_zero"

[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "electron"
ion_species = "ion"
implicit_zero_mode = true
```

この経路には `response_table_path` も事前の `matching_query.csv` も不要です。完全な CSV 不要例は
`examples/periodic2_matching_plane_zhao_implicit.toml` です。table snapshot を使う場合だけ query grid と
`beach-zhao-response` を別途使います。

implicit 化だけでは Zhao branch を選びません。既定では `auto` が現在の seed で一意な物理解を保証できない場合に停止します。
強い PE では `a` / `b` / `c` を別々に走査してから、検証した branch を明示してください。

複数根を物理量で選ぶ場合は、online Zhao に限り次を指定できます。

```toml
zhao_root_selection = "minimum_energy"
```

BEACH は multistart 探索で検出した候補について、表面から無限遠までの profile から

$$
U=-\frac{\epsilon_0}{2}\int_0^\infty E^2\,dx
$$

を計算し、最も低い $U$ を選びます。明示 branch ではその branch 内、`auto` では検証できた A / B / C 間の比較です。
別 branch の数値失敗で候補集合を確定できない場合や、最小値が相対 $10^{-6}$ 以内で縮退する場合は停止します。
この判定は [Mishra et al. (2023)](https://academic.oup.com/mnras/article/520/1/233/6987684) の
sheath potential-energy 比較に基づきますが、有限 multistart が全数学根を列挙する保証や時間依存安定性の証明ではありません。

table で PE を省略する場合は、応答表の PE flux / energy 軸も 0 の singleton にします。species、境界、
`periodic2` の完全な入力条件は[入力パラメータ](Parameters.html#matching-plane-quasistatic-closure)で確認してください。

まず、同梱の PE あり `response_backend="zhao_online"` 例をそのまま確認します。この例は 4 batch の
配線確認用であり、研究条件の妥当性を示すものではありません。repository root で次を実行します。

```bash
beachx lint examples/periodic2_matching_plane_zhao_online.toml
beach examples/periodic2_matching_plane_zhao_online.toml
beach-inspect outputs/periodic2_matching_plane_zhao_online
```

正常終了すると `outputs/periodic2_matching_plane_zhao_online/` に少なくとも `summary.txt`、`charges.csv`、
`matching_plane_history.csv` が生成されます。

`summary.txt` の `batches=4` と `matching_plane_state_valid=T` を確認してください。これは 4 accepted batch と
固定点の成立を示しますが、外部シースの物理妥当性や粒子数への収束は示しません。

`beachx lint` は TOML と既知の組合せを検査しますが、応答 CSV の内容は読みません。`response_backend="table"` では、
table の header、直積格子、整合面高度を `beach` の起動時に検査します。

## 3. 1 accepted batch で起こること

1. **整合面の状態を測る。** 現在の表面電荷と下側境界から、整合面直下の平均変位 $D_H$ を求めます。
2. **外部応答を得る。** backend は $D_H$ と外向き feedback から、整合面電位 $\Phi_H$、electron / ion の
   inward flux と access potential、PE barrier を返します。online implicit では、各 feedback 反復の内側で
   後退 Euler 終点を解き直します。
3. **同じ batch を追跡する。** $\Phi_H$ を `periodic2` の zero-mode gauge に使い、粒子を追跡して外向き moment と
   PE の return / escape を測ります。
4. **固定点を確認する。** 測定値が仮定した feedback と一致しなければ緩和して再実行します。各 trial は同じ
   batch 開始 RNG state と macro 粒子端数から再生され、未収束 trial は状態を変更しません。
5. **1 回だけ commit する。** 収束した trial だけが表面電荷、RNG、ledger、history、外部 state を更新します。
   adaptive $k\ne0$ 条件が trial を棄却した場合は、batch 幅を半分にして外部 state も巻き戻します。

この反復により Monte Carlo の乱数差ではなく、同じ粒子写像に対する外部応答の整合性を評価します。

### 0 V reservoir と PE return

外部 model の上流 plasma potential を 0 V とし、`matching_potential_v` と 3 種類の access / barrier potential は
同じ gauge を使います。これらの potential だけへ任意の定数を加えることはできません。

PE が z-high を外向きに横切ると、BEACH は局所電位、法線運動 energy、外部 PE barrier を比較します。障壁を越えない
PE は z-high で鏡面反射して return、越える PE は escape です。したがって PE return は考慮されますが、外部の
turning point までの距離や飛行時間は解かず、境界上の即時反射へ縮約しています。

ambient electron / ion の外向き粒子は局所反射しません。table は外部で戻る成分を含む総 inward flux を返せます。
online Zhao v1 は ambient 外向き feedback を transparent とするため、その外部 return が重要な case には使えません。

online Zhao の branch と barrier は次の関係です。

| branch | $\Phi_H$ | electron access / PE barrier |
|---|---:|---:|
| Type A | 正 | $\phi_m<0$ |
| Type B | 正（$D_H=0$ では 0） | 0 V |
| Type C | 負 | 0 V |

`summary.txt` の `surface_current_model_zhao_branch=auto` は選択方針であり、各 query の実 branch 名ではありません。
online 応答では、accepted state の $\Phi_H$ と access / barrier potential から上表の branch を判別します。

$H$ は interface、zero-mode gauge、PE moment の測定面を固定します。online Zhao は平面・並進対称なので、
$H$ の絶対座標を Sagdeev 方程式の距離 parameter には使わず、壁面から $H$ までの 1D profile は解きません。

## 4. 出力から成否を判断する

`output.history_stride>0` なら `matching_plane_history.csv` が生成されます。まず次を確認します。

| 確認すること | 出力 | 判定 |
|---|---|---|
| accepted state がある | `matching_plane_state_valid` | `summary.txt` で `T` |
| 固定点が収束した | `matching_plane_residual` | `surface_current_model_coupling_rtol` 以下 |
| 反復に余裕がある | `matching_plane_iterations` | 上限へ張り付かず、条件変更でも安定 |
| PE の分類が閉じる | outward / return / escape flux | $\Gamma_{pe}^{out}\simeq\Gamma_{pe}^{return}+\Gamma_{pe}^{escape}$ |
| table の由来を識別できる | response path / content fingerprint | production data と一致 |
| 電位と帯電が安定する | $D_H$、$\Phi_H$、mesh charge / potential | batch 幅・粒子数・mesh を変えて許容差内 |

accepted state の全 17 列、summary receipt、時刻の意味は
[出力形式リファレンス](OutputReference.html#matching_plane_quasistatic)にまとめています。

### 停止したときの見方

| 症状 | 主な原因 | 次に試すこと |
|---|---|---|
| response preflight 失敗 | path、header、$H$、直積格子の不一致 | [CSV 契約](MatchingPlaneReference.html#table-backend-の応答-csv-v1)と `domain.box_max` の z 成分を照合する |
| table query が範囲外 | active 軸の sweep が過渡状態を覆っていない | 外挿せず、物理的に検証した範囲で表を再生成する |
| 固定点が反復上限に到達 | 粒子 noise、強すぎる feedback、狭すぎる許容値 | ray / macro 粒子数を増やし、残差が減るなら緩和係数や上限を調整する |
| online Zhao に物理解がない、または曖昧 | $D_H$ と branch の不整合、複数根、数値失敗 | `a` / `b` / `c` を個別に scan。必要なら検証後に `minimum_energy` を使う |
| table implicit root を bracket できない | 応答表内に backward-Euler 終点がない | [`implicit_zero_mode` の契約](MatchingPlaneReference.html#implicit_zero_mode)に沿って $D_H$ 範囲を見直すか `batch_duration` を小さくする |
| online implicit root を bracket できない | Zhao branch が終わるか、幾何拡張または signed natural-scale scan で符号変化がない | branch と初期電荷を確認し、必要なら `batch_duration` を小さくする |
| soft discard 上限に到達 | 周期境界 event の未解決粒子が誤差 budget を超えた | [soft discard の停止条件](ParticleEvents.html#境界通過後の残り時間を進める)に従い、batch ごとの burst、累積率、絶対電荷を調べる |

run の正常終了が示すのは、backend 評価と数値的な固定点収束です。外部シースの物理妥当性、matching-plane 高度への
不変性、Monte Carlo 収束までは証明しません。

## 5. 受理される構成と適用限界

matching-plane は平均場や粒子 channel の二重計上を防ぐため、次の構成に限定されます。

| 項目 | 必要条件 |
|---|---|
| box / 場 | x / y periodic、z open、`field_boundary.mode="periodic2"`、`sim.e0=sim.b0=[0,0,0]` |
| periodic split | `cached_kneq0` または `panel_spectral_reference`、`exclude_k0`、`symmetric_vacuum` または `e_bottom_zero` |
| reservoir / open 面 | `[reservoir].inflow_model="source_vdf"`、`ordinary_open_model="escape"` |
| ambient species | electron と ion だけを role に指定。`volume_seed`、`npcls_per_step=0`、z-high reservoir 流入 |
| PE species | 任意。負電荷の `photo_raycast`、z-high 注入、放出反作用あり |
| surface closure | 全 role が `explicit`。手動 `fixed_current` target と `neutral_return` は使わない |
| event policy | `abort`、または[率・件数猶予・絶対電荷で制限した `soft_discard`](ParticleEvents.html#境界通過後の残り時間を進める) |

`reference_area_m2` と stationary Zhao 専用の source key は指定しません。面積は domain の x-y 面積、$H$ は
`domain.box_max` の z 成分、更新間隔は 1 accepted batch から決まります。online Zhao 固有の単価電荷、温度、密度、drift の
制約は[入力パラメータ](Parameters.html#matching-plane-quasistatic-closure)を参照してください。

この model は次を解きません。

- 外部の 6D VDF、粒子 inventory、flight time、遅延 return queue
- 衝突、磁化された return、外部シースの過渡
- BEACH 領域内の volume plasma charge
- online Zhao v1 における ambient 外向き population の外部 return

online Zhao は PE の束と平均法線 energy を、その 2 moment を再現する half-Maxwellian へ縮約します。
そのため、高 energy tail は保持しません。

`auto` の複数根判定は有限個の multistart で見つけた根の比較であり、数学的な root isolation ではありません。
branch 別の検証では `a` / `b` / `c` を明示して scan します。

各 query は物理的に stateless で、前の root からの continuation は行いません。解けない場合も明示 branch や
backend を暗黙に切り替えません。

これらが主要効果なら、独立した 1D--3D kinetic coupling または full PIC で検証します。応答表の grid、固定点許容値、
matching-plane 高度を変える検証項目は[数値・応答表リファレンス](MatchingPlaneReference.html#収束と適用性を検証する)にまとめています。
