title: 実行・再開する

# 実行・再開する

このページは、設定済みの `beach.toml` を安全に実行し、必要なら checkpoint から続きを計算する手順です。
入力検査、負荷見積もり、実行、出力確認の順に進めます。再開時は元の出力を残し、累積 `batch_count` を
増やした別設定へ引き継ぎます。

未インストールの場合は[インストール](../start/Installation.md)、最初のケースを作る場合は
[入門チュートリアル](../start/Tutorial.md)から始めてください。

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

1. 実行前に `beachx lint` で TOML、JSON Schema、座標・配置パラメータの組合せ、既知の制約を検査し、`status=ok` を確認します。
2. `beach` に設定ファイルを渡してシミュレーションを実行します。
3. `output.dir` に作成された結果を `beachx inspect` で確認します。

設定検査で確認する範囲と外部データの扱いは[設定を検証する](Configuration.md)を参照してください。

引数を省略した`beach`は、カレントディレクトリの`beach.toml`を読み込みます。

## 並列実行

OpenMPのスレッド数は環境変数で指定します。

```bash
OMP_NUM_THREADS=8 beach beach.toml
```

`beach` は `OMP_SCHEDULE` が未指定なら、通常の粒子追跡を `dynamic,1`（1 粒子ずつ動的に割り当て）で実行します。
粒子ごとに飛行時間が異なっても、空いたスレッドが次の粒子を処理します。明示した `OMP_SCHEDULE` を優先するので、
固定割り当てで比較したい場合は次のように指定します。

```bash
OMP_NUM_THREADS=8 OMP_SCHEDULE=static beach beach.toml
```

適応的な $k\ne0$ 進行は、同じバッチの再試行を再現するため、粒子追跡中だけ `static` を使います。
並列実行では電荷の加算順序が変わり、丸め誤差によって個々の粒子の結果が変わる場合があります。

MPIビルドを利用する場合はランチャーから起動します。

```bash
mpirun -n 4 beach beach.toml
```

MPI build、launcher、compiler module は利用する環境に合わせてください。MPI / OpenMP 内部の state 所有と
集約方法を変更する開発者は[ランタイムのアーキテクチャ](../developer/Architecture.md)を参照してください。

バッチの粒子は種・source・面ごとの global slot に並び、各 slot は (seed, batch, slot) の PCG64 counter-based 乱数で
生成されます。既定（`run.batch.particle_distribution = "static"`）は slot を rank に連続区間で割り振るので、
rank 数を変えても同じ粒子集合になります。rank 間で粒子追跡の所要時間が偏るときは `"shared"` にすると、全 rank が
1 つの slot pool を `particle_pool_chunk` 個ずつ取り合い（rank 0 のカウンタへの fetch-and-add）、負荷が動的に均されます。
`summary.txt` の `particle_pool_claim_seconds` で claim の待ち時間を確認できます。各 rank は次の chunk を先に要求して
粒子追跡と通信を重ねるので、ノードをまたぐ実行でも既定の 64 から 256 程度の chunk なら claim の待ち時間はバッチ時間の
1% 未満です。

## 実行前に負荷を見積もる

`inflow`、`photo_raycast`ではバッチごとの粒子数が動的に決まるため、実行前の
見積もりを推奨します。

```bash
beachx lint beach.toml
beachx workload beach.toml --threads 8
```

`workload` は粒子数の計算に必要な流束・重み・バッチ幅を解決します。設定の組合せは `lint` で検査します。
光電子は `rays_per_batch` による放出粒子数の上限を数え、光線と表面の交差や放出後の軌道は計算しません。
Python の見積もり自体には numpy や native field kernel は不要です。

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

正常終了時の最終出力は、`output.checkpoint.stride_batches=0` でも checkpoint として再開できます。
最初の再開確認では、元の `outputs/tutorial` を残し、続きを別のディレクトリへ書く方法を推奨します。

### 公式入門ケースを 1 バッチ延長する

**前提:** [入門チュートリアル](../start/Tutorial.md)が完了し、作業ディレクトリに `beach.toml` と
`outputs/tutorial` があることを確認します。完了済みバッチ数を読みます。

```bash
grep '^batches=' outputs/tutorial/summary.txt
```

公式入門ケースでは `batches=150` と表示されます。一般には、この完了済みバッチ数を `B` とします。

元の設定を残して再開用の設定を作ります。

```bash
cp beach.toml resume.toml
```

`resume.toml` の既存の `[run]` と `[output]` にある該当値を変更します。
次は公式入門ケースを 150 batch から 151 batch へ延長する完全な例です。

```toml
[run]
batch_count = 151

[output]
enabled = true
dir = "outputs/resumed"

[run.restart]
from = "outputs/tutorial"
```

`resume.toml` の他の設定は変更しません。`run.restart.from` は checkpoint の読み込み元、`dir` は新しい出力先です。
チュートリアルと同じ作業ディレクトリで、入力検査、実行、出力確認を順に行います。

```bash
beachx lint resume.toml
beach resume.toml
beachx inspect outputs/resumed
```

合格時は、`beach` と `beachx inspect` の出力に次が含まれます。

```text
resuming_from_batches=150
...
batches=151 ...
```

別のケースでも、完了済みバッチ数が `B` なら `batch_count` を整数 `B+1` に変更し、
`resuming_from_batches=B` と `batches=B+1` を確認します。
`run.batch_count` は追加バッチ数ではなく、累積の到達バッチ数です。

### 長時間実行の定期 checkpoint

accepted batch 数で定期 checkpoint を有効にできます。

```toml
[output]
dir = "outputs/tutorial"

[output.checkpoint]
stride_batches = 1000
```

この例は 1000 accepted batch ごとに再開状態を二重 slot へ保存します。`0` は定期保存を無効にしますが、
正常終了時の最終 checkpoint は引き続き出力します。

### 同じ出力先へ続ける場合と MPI 再開

元の出力ディレクトリへ続けて書く場合は、`run.restart.from` を省略し、`dir` を checkpoint のディレクトリにします。
元の最終出力をそのまま比較に残したい場合は、上の別ディレクトリ手順を使ってください。

```toml
[output]
enabled = true
dir = "outputs/tutorial"

[run.restart]
```

一般に、checkpoint が `batches=100` で新しい `batch_count=150` なら 50 バッチを追加します。
古い slot から再開した場合も、history はその checkpoint の batch まで切り戻してから追記します。
`[run.restart]` を書かずに既存の出力先で実行すると新規実行として扱い、前の run の checkpoint と history を
消します（[出力先の世代](../reference/OutputReference.md#出力先の世代)）。
MPI 再開の rank 数は保存時と違ってかまいません。乱数はすべて (seed, batch, 粒子 slot, ...) で決まる PCG64 の
counter-based 乱数で、checkpoint は乱数状態を持ちません。同じ seed なら rank 数によらず同じ粒子集合を追跡します。
適応的な $k\ne0$ 進行では、同じ実行内の再試行で全 MPI rank の実 OpenMP team size を揃えます。
restart 前後で同じ team size を使う必要はなく、再開後の実 team size を新しい診断値として記録します。

必須ファイルと checkpoint 選択の契約は
[再開に使うファイル](../reference/OutputReference.md#再開に使うファイル)が正本です。
再開が拒否された場合は[トラブルシューティング](Troubleshooting.md#checkpointから再開できない)で
必須ファイル、mesh 識別子、累積 `batch_count` を確認してください。

## 実行後に確認すること

- processの終了codeが0である
- `summary.txt`の`batches`が`run.batch_count`と一致する
- `absorbed`、`escaped_boundary`、`survived_max_step`の内訳を確認する
- `output.diagnostics.charge_rel_change_threshold`を自動停止条件として扱わない

ファイルの意味は[出力ファイルを調べる](OutputGuide.md)、図の作成は
[後処理チュートリアル](PostprocessTutorial.md)、結果の判定は
[計算結果の妥当性確認](ValidationGuide.md)へ進んでください。
