title: 計算モデルの全体像

Lang: [日本語](Algorithms.md) | [English](Algorithms.en.md)

# 計算モデルの全体像

BEACHは、表面電荷が作る電場と、その電場中を運動する荷電粒子をbatch単位で結合します。
三角形要素上の電荷から電場と電位を評価し、粒子が衝突した要素へ新たな電荷を蓄積することで、
表面帯電の時間発展を計算します。

## n batchの計算フロー

`n = sim.batch_count`です。再開時はcheckpointに記録された次のbatchから処理を再開し、
完了したbatchの総数が`n`に達するまで同じ流れを繰り返します。

```mermaid
flowchart TD
    start(["初期化 / checkpoint読込"])
    more{"batch i < n ?"}

    subgraph batch["batch i"]
        direction TB
        subgraph prepare["snapshot / source"]
            direction LR
            p1["1. snapshot更新"] --> p2["2. 粒子生成"]
        end

        subgraph particle_loop["粒子step loop"]
            direction LR
            p3["3. 1 step前進"] --> p4["4. 衝突・境界処理"] --> tracking{"続ける?"}
            tracking -- "次step" --> p3
        end

        subgraph batch_end["batch末尾"]
            direction LR
            p5["5. 結果集約"] --> p6["6. 電荷commit"] --> p7["7. 統計・履歴"]
        end

        prepare --> particle_loop --> batch_end
    end

    finish(["最終出力 / checkpoint"])
    start --> prepare
    more -- "はい" --> prepare
    batch_end --> more
    more -- "いいえ" --> finish
```

同じbatchの粒子は、手順1で確定した同じfield snapshotを共有します。手順6で更新した表面電荷が
場へ反映されるのは、次batchの手順1です。

### 1. 場とouter modelのsnapshotを更新

前batchまでにcommitされた`q_elem`から、粒子追跡中に固定する電場と電位を作ります。direct、treecode、
FMMの選択は[場の評価](FieldSolvers.html)、周期和は[periodic2静電場](PeriodicElectrostatics.html)、
外部領域との結合は[外部プラズマモデル](OuterPlasmaModels.html)を参照してください。

### 2. batch粒子を生成

`volume_seed`、`reservoir_face`、`photo_raycast`から、そのbatchで追跡する粒子を作ります。reservoir粒子の
個数は流入fluxと`batch_duration`から決まり、光電子はrayが最初にhitした表面から放出されます。全体は
[粒子源](ParticleSourcesBoundaries.html)、詳細は[reservoir注入](ReservoirInjection.html)と
[光電子の放出とライフサイクル](PhotoelectronEmission.html)を参照してください。

### 3. 粒子を1 step前進

固定snapshot中の電場と任意の一様磁場を使い、速度をBoris法、位置を同時刻状態の台形則で更新します。この時点では
候補軌道を作るだけで、step内の最終状態は衝突・境界処理後に確定します。式と更新順は
[Boris粒子更新](BorisPusher.html)を参照してください。

### 4. 最初の衝突・境界イベントを判定・処理

候補軌道上の三角形hitとbox面交差を比較し、最初に起きるものを選びます。その位置で吸収、reflect、
periodic wrap、escape、outer returnのいずれかを処理します。粒子が生存し、stepの残り時間があれば手順3へ
戻ります。判定順は[粒子の衝突・境界イベント](ParticleEvents.html)、外部領域での処理は
[粒子のescapeとreturn](ParticleEscapeReturn.html)を参照してください。

### 5. batch結果を集約

表面hitによる電荷差分、光電子放出の反作用電荷、粒子outcome、outer-interface診断を集約します。
OpenMPでは衝突電荷をthread-localに保持し、必要な量をMPI rank間で集約します。並列実行の構成は
[実行する](Execution.html)を参照してください。

### 6. 表面電荷をcommit

全粒子の電荷差分を`q_elem`へ一度だけ加え、必要なら浮遊導体を総電荷一定のまま等電位化します。絶縁体帯電、
光電子放出の符号、conductor処理、charge ledgerは[表面電荷更新](SurfaceModels.html)を参照してください。

### 7. 統計と履歴状態を更新

吸収、escape、`max_step`まで生存した粒子の数と`tol_rel`を更新し、設定したstrideで電荷・電位履歴を保存します。
`tol_rel`は監視値であり、早期終了条件ではありません。出力の意味は[出力の読み方](OutputGuide.html)、再開方法は
[実行する](Execution.html)を参照してください。`n` batch完了後に最終結果とcheckpointを書きます。

## batch間で引き継ぐもの

| 状態 | 次batchでの役割 |
| --- | --- |
| 要素電荷`q_elem` | 手順1のfield sourceになる |
| reservoir端数、RNG、outer state | 次の粒子生成やouter refreshを継続する |
| 統計と履歴 | 累積結果とrestart位置を保持する |
| model / mesh / species fingerprint | checkpointとの互換性を検証する |

`sim.dt`は手順3の粒子step幅、`batch_duration`は手順2の粒子供給量と手順6の電荷更新を結ぶ時間幅です。
`batch_duration`依存性は[`batch_duration`の安定性と定常値](BatchDurationStability.html)で確認してください。

周期境界を含む計算全体の組み方は、[periodic2有限画像構成](FinitePeriodicConfiguration.html)と
[periodic2無限周期＋outer plasma構成](InfinitePeriodicOuterConfiguration.html)にまとめています。設定キーは
[設定パラメータ](Parameters.html)、離散化と結果の収束確認は[計算結果の妥当性確認](ValidationGuide.html)を参照してください。
