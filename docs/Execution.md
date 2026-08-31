title: 実行・再開する

Lang: [日本語](Execution.md) | [English](Execution.en.md)

# 実行・再開する

このページは、設定済みの `beach.toml` を安全に実行し、必要なら checkpoint から続きを計算する手順です。
入力検査、負荷見積もり、実行、出力確認の順に進めます。再開時は元の出力を残し、累積 `batch_count` を
増やした別設定へ引き継ぎます。

未インストールの場合は[インストール](Installation.html)、最初のケースを作る場合は
[10 分チュートリアル](Tutorial.html)から始めてください。

KUDPC の login node では `beach` や `mpirun` を直接実行しません。
他の HPC 環境でもサイトの login-node 運用規則に従ってください。
ローカル環境、または計算ノードの割当内で以下のコマンドを実行してください。
KUDPC では短い確認を `tssrun`、長い計算を `sbatch` 内の `srun` で実行します。job script の例は
[`examples/job_scripts/`](https://github.com/Nkzono99/BEACH/tree/main/examples/job_scripts)を参照してください。

## 基本フロー

次は公式入門ケースの `output.dir="outputs/tutorial"` を使う例です。
別のケースでは、最後の引数をその `output.dir` に置き換えてください。

```bash
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/tutorial
```

1. `beachx lint`でTOML、JSON Schema、座標・配置パラメータの組合せ、既知の制約を検査します。
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

MPI build、launcher、compiler module は利用する環境に合わせてください。MPI / OpenMP 内部の state 所有と
集約方法を変更する開発者は[ランタイムのアーキテクチャ](Architecture.html)を参照してください。

## 実行前に負荷を見積もる

`boundary_inflow`、`plane_source`、deprecatedな`reservoir_face`、`photo_raycast`ではバッチごとの粒子数が動的に決まるため、実行前の
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
beachx profile outputs/tutorial/performance_profile.csv \
  --save outputs/tutorial/performance_profile.png
```

スケーリング比較には `performance_profile.csv` の `simulation_total` 行にある `rank_max_s` を使います。
初期化、場の更新、粒子追跡、電荷更新、MPI、出力を別の phase として比較できます。

## checkpointから一度再開する

正常終了時の最終出力は、`checkpoint_stride=0` でも checkpoint として再開できます。
最初の再開確認では、元の `outputs/tutorial` を残し、続きを別のディレクトリへ書く方法を推奨します。

### 公式入門ケースを 1 バッチ延長する

**前提:** [10 分チュートリアル](Tutorial.html)が完了し、作業ディレクトリに `beach.toml` と
`outputs/tutorial` があることを確認します。完了済みバッチ数を読みます。

```bash
grep '^batches=' outputs/tutorial/summary.txt
```

公式入門ケースでは `batches=20` と表示されます。一般には、この完了済みバッチ数を `B` とします。

元の設定を残して再開用の設定を作ります。

```bash
cp beach.toml resume.toml
```

`resume.toml` の既存の `[sim]` と `[output]` にある該当値を変更します。
次は公式入門ケースを 20 batch から 21 batch へ延長する完全な例です。

```toml
[sim]
batch_count = 21

[output]
write_files = true
dir = "outputs/resumed"
resume = true
restart_from = "outputs/tutorial"
```

`resume.toml` の他の設定は変更しません。`restart_from` は checkpoint の読み込み元、`dir` は新しい出力先です。
チュートリアルと同じ作業ディレクトリで、入力検査、実行、出力確認を順に行います。

```bash
beachx lint resume.toml
beach resume.toml
beachx inspect outputs/resumed
```

合格時は、`beach` と `beachx inspect` の出力に次が含まれます。

```text
resuming_from_batches=20
...
batches=21 ...
```

別のケースでも、完了済みバッチ数が `B` なら `batch_count` を整数 `B+1` に変更し、
`resuming_from_batches=B` と `batches=B+1` を確認します。
`sim.batch_count` は追加バッチ数ではなく、累積の到達バッチ数です。

### 長時間実行の定期 checkpoint

accepted batch 数で定期 checkpoint を有効にできます。

```toml
[output]
dir = "outputs/tutorial"
checkpoint_stride = 1000
```

この例は 1000 accepted batch ごとに再開状態を二重 slot へ保存します。`0` は定期保存を無効にしますが、
正常終了時の最終 checkpoint は引き続き出力します。

### 同じ出力先へ続ける場合と MPI 再開

元の出力ディレクトリへ続けて書く場合は、`restart_from` を省略し、`dir` を checkpoint のディレクトリにします。
元の最終出力をそのまま比較に残したい場合は、上の別ディレクトリ手順を使ってください。

```toml
[output]
write_files = true
dir = "outputs/tutorial"
resume = true
```

一般に、checkpoint が `batches=100` で新しい `batch_count=150` なら 50 バッチを追加します。
MPI 再開では、checkpoint の `mpi_world_size` と現在の rank 数を一致させます。
適応的な $k\ne0$ 進行では、同じ実行内の再試行で全 MPI rank の実 OpenMP team size を揃えます。
restart 前後で同じ team size を使う必要はなく、再開後の実 team size を新しい診断値として記録します。

必須ファイルと checkpoint 選択の契約は
[再開に使うファイル](OutputReference.html#再開に使うファイル)が正本です。
再開が拒否された場合は[トラブルシューティング](Troubleshooting.html#checkpointから再開できない)で
必須ファイル、fingerprint、MPI world size、累積 `batch_count` を確認してください。

## 実行後に確認すること

- processの終了codeが0である
- `summary.txt`の`batches`が`sim.batch_count`と一致する
- `absorbed`、`escaped_boundary`、`survived_max_step`の内訳を確認する
- `tol_rel`を自動停止条件として扱わない

ファイルの意味は[出力ファイルを調べる](OutputGuide.html)、図の作成は
[後処理チュートリアル](PostprocessTutorial.html)、結果の判定は
[計算結果の妥当性確認](ValidationGuide.html)へ進んでください。
