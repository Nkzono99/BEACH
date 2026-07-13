title: 実行する

Lang: [日本語](Execution.md) | [English](Execution.en.md)

# 実行する

このページでは、設定済みの`beach.toml`を検査して実行し、必要に応じて負荷見積もりや
再開を行うまでの流れを説明します。インストール前の場合は[インストール](Installation.html)、
最初のケースを作る場合は[10分チュートリアル](Tutorial.html)から始めてください。

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
[開発環境とテスト](Workflow.html)を参照してください。

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
