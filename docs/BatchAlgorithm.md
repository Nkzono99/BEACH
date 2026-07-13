title: batch結合アルゴリズム

Lang: [日本語](BatchAlgorithm.md) | [English](BatchAlgorithm.en.md)

# batch結合アルゴリズム

BEACHは、表面電荷が作る場と粒子輸送をbatch単位で結合します。このページでは初期化から
電荷commitまでの更新順序と、順序から生じる数値的な意味を説明します。

## 初期化

実行開始時に、概ね次の順序で状態を構築します。

1. `beach.toml`を読み、高水準記法を正規化して物理制約を検査する
2. 三角形メッシュとsurface modelを構築する
3. 粒子speciesと注入状態を初期化する
4. 場ソルバーと周期operatorを準備する
5. 新規実行またはcheckpointから表面電荷、統計、RNG状態を設定する
6. 出力先と履歴writerを準備する

再開時はmodel、mesh、speciesのfingerprintがcheckpointと一致しなければ停止します。詳しくは
[実行する](Execution.html#再開実行)を参照してください。

## batch loop

設定した`sim.batch_count`へ到達するまで、次の処理を繰り返します。

```text
q_elem(batch start)
        │
        ├─ 場ソルバーのstateを更新
        ├─ 必要なpotential／boundary profileを更新
        └─ particle sourceを準備
                    ↓
              粒子を生成・追跡
                    ↓
       absorbed / escaped / unresolvedを集計
                    ↓
             delta_q_elemを集約
                    ↓
       q_elem ← q_elem + delta_q_elem
                    ↓
       surface model処理・統計・履歴を更新
```

### 1. 場の更新

batch開始時の`q_elem`から、direct、treecode、FMMいずれかの場ソルバーstateを更新します。
`periodic2`では近傍画像、cached nonzero operator、zero mode、outer modelを選択した契約に従って
組み合わせます。

### 2. 粒子生成

speciesごとにvolume seed、reservoir、photoelectronなどのsourceを評価します。reservoirの
非整数macro-particle数は残差として次batchへ持ち越します。

### 3. 粒子追跡

各粒子は最大`sim.max_step`まで進みます。各stepでは場を評価して速度と位置を更新し、box境界または
三角形との最初のeventを確定します。

- 表面へ衝突: 吸収し、命中要素の`delta_q_elem`へ電荷を加える
- open boundaryから脱出: boundary modelに従ってescapeまたはreturnを判定する
- `max_step`へ到達: `survived_max_step`として記録し、吸収やescapeへ読み替えない

粒子積分と衝突判定は[粒子追跡と衝突](ParticleTrackingCollision.html)を参照してください。

### 4. 電荷commit

粒子追跡中は`q_elem`を変更せず、命中要素ごとの`delta_q_elem`を集約します。全粒子の処理後に
差分をcommitするため、同じbatchの粒子同士は、そのbatch内で堆積した電荷を介して相互作用しません。

insulatorでは命中電荷をその要素に保持します。floating conductorを使う場合は、commit後に
導体要素間で電荷を再配分して等電位条件を近似します。

### 5. 統計と履歴

commit後の状態からbatch統計と履歴を更新します。`tol_rel`は表面電荷変化のmonitoring metricであり、
現行実装では早期停止条件ではありません。通常実行は必ず設定した`sim.batch_count`まで進みます。

## 並列実行で保つ契約

OpenMPでは粒子処理を並列化し、thread局所の統計と電荷差分をbatch末尾で集約します。MPIでは粒子を
rankへ分配し、要素電荷差分と統計をcollectiveで統合します。

並列方式にかかわらず、以下を維持します。

- batch開始時の表面電荷から場を作る
- 粒子処理中に共有`q_elem`を更新しない
- commit後に全rankが同じ表面電荷を持つ
- 電荷ledgerで注入、吸収、escape、未解決を区別する

## 時間刻みの区別

| 値 | 役割 |
| --- | --- |
| `sim.dt` | 粒子軌道の1 step |
| `sim.max_step` | 1粒子を追跡する最大step数 |
| `batch_duration` | 1 batchが代表する物理時間と粒子供給量 |
| `sim.batch_count` | 実行する累積batch数 |

`dt`を小さくする確認と`batch_duration`を変える確認は別の収束軸です。後者は
[`batch_duration`の安定性と定常値](BatchDurationStability.html)を参照してください。
