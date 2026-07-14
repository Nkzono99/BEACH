title: 実行する

Lang: [日本語](Execution.md) | [English](Execution.en.md)

# 実行する

設定済みの`beach.toml`は、入力検査、実行、出力確認の順に扱います。計算規模に応じて、
実行前に負荷も見積もってください。既存のcheckpointからは、同じ設定を使って計算を再開できます。
未インストールの場合は[インストール](Installation.html)、最初のケースを作る場合は
[10分チュートリアル](Tutorial.html)から始めてください。

## 基本フロー

```bash
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

1. `beachx lint`でTOML、JSON Schema、高水準記法、既知の制約を検査します。
2. `beach`に設定ファイルを渡してシミュレーションを実行します。
3. `output.dir`に作成された結果を`beachx inspect`で確認します。

引数を省略した`beach`は、カレントディレクトリの`beach.toml`を読み込みます。

## 並列実行

OpenMPのスレッド数は環境変数で指定します。

```bash
OMP_NUM_THREADS=8 beach beach.toml
```

MPIビルドを利用する場合はランチャーから起動します。

```bash
mpirun -n 4 beach beach.toml
```

MPIとOpenMPの組み合わせ、コンパイラ設定、KUDPCでの実行方法は
[開発環境とテスト](Workflow.html)にまとめています。

OpenMPでは、粒子indexを`dynamic, 1`で分配します。これにより、粒子ごとの追跡step数の違いから生じる
負荷の偏りを抑えます。衝突による電荷変化は`dq_thread(nelem, nth)`にthread-localで集計し、最後に結合します。

MPIでは粒子の生成と追跡をrank間で分担し、各rankが同じmeshと`q_elem`を保持します。batch末尾に
`dq(nelem)`と終了状態別の粒子数をallreduceするため、commit後の表面電荷と統計は全rankで一致します。
通常の結果、履歴、全rankで共有するmacro residualはroot rankが書き出し、RNG stateはrankごとに保存します。

## 実行前に負荷を見積もる

`reservoir_face`や`photo_raycast`ではバッチごとの粒子数が動的に決まるため、実行前の
見積もりを推奨します。

```bash
beachx workload beach.toml --threads 8
```

MPI rankを含めて見積もる場合:

```bash
beachx workload beach.toml \
  --threads 8 \
  --mpi-ranks 4 \
  --mpi-rank 0
```

## 性能計測

粗いフェーズ計測を有効にする場合は`BEACH_PROFILE=1`を設定します。

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach beach.toml
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

スケーリング比較には`performance_profile.csv`の`simulation_total`行にある`rank_max_s`を使います。
初期化、batch準備、場のrefresh、粒子追跡、電荷commit、MPI集約、統計・履歴更新、
結果・checkpoint出力を計測します。CSV上では`load_or_init`、`field_solver_init`、`prepare_batch`、
`field_refresh`、`particle_batch`、`commit_charge`、`mpi_reduce`、`stats_update`、`history_write`、
`write_results`、`write_checkpoint`として区別されます。

## 再開実行

同じ出力ディレクトリから再開する場合:

```toml
[output]
dir = "outputs/latest"
resume = true
```

checkpointと新しい出力先を分ける場合は`restart_from`を指定します。

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../parent_run/outputs/latest"
```

`sim.batch_count`は追加バッチ数ではなく、累積の到達バッチ数です。既存checkpointが
`batches=100`で`batch_count=150`なら、追加で50バッチ実行します。MPI再開では保存時と現在の
rank数が一致する必要があります。

## 実行後に確認すること

- processの終了codeが0である
- `summary.txt`の`batches`が`sim.batch_count`と一致する
- `absorbed`、`escaped_boundary`、`survived_max_step`の内訳を確認する
- `tol_rel`を自動停止条件として扱わない

ファイルの意味は[出力の読み方](OutputGuide.html)、図の作成は
[後処理チュートリアル](PostprocessTutorial.html)、結果の判定は
[計算結果の妥当性確認](ValidationGuide.html)へ進んでください。
