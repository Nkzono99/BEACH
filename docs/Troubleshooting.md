title: トラブルシューティング

Lang: [日本語](Troubleshooting.md) | [English](Troubleshooting.en.md)

# トラブルシューティング

最初に症状に近い節を選び、確認コマンドの結果と「期待する状態」を比較してください。
エラーを回避するために未対応モデルを強制したり、数値許容値だけを緩めたりしないでください。

KUDPC の login node では `beach`、`mpirun`、長時間の解析を直接実行しません。
以下の実行コマンドは、ローカル環境または計算ノードの割当内で実行します。

以下の一般例では出力先を `outputs/latest` とします。実際の `beach.toml` にある `output.dir` へ
読み替えてください。公式入門ケースの出力先は `outputs/tutorial` です。

## `beach` または `beachx` が見つからない

**確認コマンド:** 現在の Python、build tool、Fortran compiler、BEACH コマンドがどこから読まれるかを確認します。

```bash
command -v python
command -v make
command -v gfortran
command -v beach
command -v beachx
python -m site --user-base
beach --version
beachx --help
```

`gfortran` は例です。`FC` または `FPM_FC` で別の Fortran compiler を選んだ場合は、そのコマンドを確認します。
checkout を `--no-build-isolation` でインストールする場合、または `make` を直接使う場合だけ、
`command -v fpm` も確認します。

**期待:** 必要なコマンドごとに実行ファイルの絶対パスが表示され、`beach --version` と `beachx --help` が終了コード `0` になります。

**安全な対処:**

- `beach` / `beachx` だけがなければ、[インストール](Installation.html)に従って現在の Python 環境へ再導入します。
- user base 配下へ入っている場合は、表示された user base の `bin` を `PATH` へ追加します。
- `make` または compiler がなければ、OS package または HPC サイトの compiler module を用意してから再インストールします。
- 別の Python 環境の `pip` と混在させず、`python -m pip ...` を使います。

## `beachx lint` が失敗する

**確認コマンド:** 実行に使うのと同じファイルを検査します。

```bash
beachx lint beach.toml
```

**期待:** TOML、JSON Schema、BEACH 制約を通過し、最後に `status=ok` が表示されます。

**安全な対処:** 最初に表示された未知キー、型、値域、または排他条件のエラーを 1 件ずつ修正します。
[設定のよくある失敗](Configuration.html#4-よくある失敗)と[入力パラメータ](Parameters.html)でキーの所属、型、単位を確認します。

## lint は通るが `beach` が途中で停止する

**確認コマンド:** 標準出力と標準エラーを保存し、終了コードと最後の診断を確認します。

```bash
beach beach.toml > beach-run.log 2>&1
run_status=$?
printf 'exit_code=%s\n' "$run_status"
tail -n 50 beach-run.log
```

**期待:** `exit_code=0` になり、ログの最後に実行統計と `results written to ...` が表示されます。

**安全な対処:** `ERROR STOP` またはその直前の最初の診断から修正します。
Fortran の物理制約検査は Python schema より厳しい場合があるため、lint 成功だけで実行可能とは限りません。
互換表や必須条件に合わせて設定を修正し、エラー回避のためだけに許容値を緩めないでください。

## 出力ディレクトリがない、または `beachx inspect` が読めない

**確認コマンド:** 実行時の作業ディレクトリ、`[output]` の設定、必須の最終出力を確認します。

```bash
pwd
grep -A 10 '^\[output\]' beach.toml
ls -l outputs/latest/summary.txt outputs/latest/charges.csv
beachx inspect outputs/latest
```

**期待:** `write_files=true` で、`dir` の下に `summary.txt` と `charges.csv` があり、`beachx inspect` が終了コード `0` になります。

**安全な対処:** 相対パスの `output.dir` は `beach` を起動した作業ディレクトリから確認します。
先に `beach` の終了コードとログを修正し、空の `summary.txt` や `charges.csv` を手作業で作らないでください。

## 実行は終わるが粒子統計が期待と合わない

**確認コマンド:** 実行完了と各粒子の終了状態を分けて読みます。

```bash
beachx inspect outputs/latest
grep -E '^(processed_particles|absorbed|escaped|escaped_boundary|survived_max_step|multiple_box_events_soft_discarded|batches)=' \
  outputs/latest/summary.txt
```

**期待:** `batches` が `sim.batch_count` と一致し、粒子数の内訳を source、吸収、境界 escape、
`survived_max_step`、および未解決の soft discard から説明できます。
`survived_max_step` には全ケース共通の合格値はありません。

**安全な対処:** `survived_max_step` が目的量に影響する場合は、`dt`、`max_step`、box サイズ、注入速度を
1 項目ずつ変えて依存性を確認します。未解決粒子を後から吸収または escape として数え直さないでください。
受理判定は[計算結果の妥当性確認](ValidationGuide.html)に従います。

## 履歴が空、または大きすぎる

**確認コマンド:** 保存間隔と実際の行数を確認します。

```bash
grep -A 10 '^\[output\]' beach.toml
ls -lh outputs/latest/*history.csv
wc -l outputs/latest/*history.csv
```

**期待:** `history_stride > 0` なら `charge_history.csv` が生成されます。
`potential_history.csv` は `write_potential_history=true` かつ `history_stride > 0` のときだけ生成されます。

**安全な対処:** 空な場合は `history_stride`、`batch_count`、実際の `output.dir` を確認します。
大きすぎる場合は `history_stride` を増やし、電位履歴が必要なケースだけ `write_potential_history` を有効にします。

## checkpointから再開できない

**確認コマンド:** まず `restart_from` が指すディレクトリの直下を確認します。
`restart_from` を省略した同一ディレクトリ再開では `output.dir` を確認します。
次は serial 出力 `outputs/latest` の例です。

```bash
checkpoint_dir=outputs/latest
grep -A 12 '^\[sim\]' resume.toml
grep -A 12 '^\[output\]' resume.toml
ls -l "$checkpoint_dir/summary.txt" \
  "$checkpoint_dir/charges.csv" \
  "$checkpoint_dir/checkpoint_complete.txt"
ls -l "$checkpoint_dir"/rng_state*.txt
grep -E '^(checkpoint_schema_version|batches|model_fingerprint|mesh_fingerprint|species_fingerprint|mpi_world_size)=' \
  "$checkpoint_dir/summary.txt"
sed -n '1,80p' "$checkpoint_dir/checkpoint_complete.txt"
```

**期待:** 現行 schema の checkpoint では次の条件をすべて満たします。

- `summary.txt`、`charges.csv`、serial の `rng_state.txt` または MPI の全 `rng_state_rankNNNNN.txt` がある
- schema v8 以降の `checkpoint_complete.txt` が `state=complete` で、その `batches` と `mpi_world_size` が `summary.txt` と一致する
- manifest が宣言する `macro_residuals.csv` と `charge_ledger.csv` がある
- ordered mesh fingerprint が一致する。model / species fingerprint の不一致は warning 付きで継続できる
- MPI 再開では保存時の `mpi_world_size` と現在の rank 数が一致する
- `resume.toml` で `output.write_files=true`、`output.resume=true` である。`output.restart_from` を指定した場合は
  確認中の checkpoint を指し、`output.dir` は意図した新しい出力先である。省略した場合は `output.dir` 自体が再開元である
- 新しい `sim.batch_count` が保存済み `batches` 以上である。1 バッチ以上進めるならより大きい

**安全な対処:**

- `restart_from` は checkpoint の読み込み元だけを指定し、新しい結果は別の `output.dir` へ書いて元出力を保存します。
- 元の物理・数値設定を変えて再開する場合は warning を記録し、出力先を分けて変更前後の連続性を確認します。
- 異なる世代のファイルを手作業で混ぜたり、`checkpoint_complete.txt` を書き換えたりしません。
- `resume=false` にして同じ `output.dir` へ新規実行することを、再開エラーの回避策にしません。
- 定期 slot がある場合、BEACH は直下の最終出力と `checkpoints/slot0`、`slot1` から完全で最新の候補を自動選択します。

修正後は[実行と再開](Execution.html#checkpointから一度再開する)の手順で別出力先へ再開し、次を確認します。

```bash
beachx lint resume.toml
beach resume.toml > resume.log 2>&1
grep '^resuming_from_batches=' resume.log
beachx inspect outputs/resumed
```

`resuming_from_batches` が使用した checkpoint の `batches` と一致し、新しい出力の `batches` が
`resume.toml` の `sim.batch_count` と一致すれば復旧完了です。
必須ファイルと自動選択の正本は[再開に使うファイル](OutputReference.html#再開に使うファイル)です。

## 表面モデルまたは periodic 構成が拒否される

**確認コマンド:** まず入力検査で拒否された組合せを確認します。

```bash
beachx lint beach.toml
```

lint は通るが実行時に拒否される場合だけ、実行時ログを採取します。

```bash
beach beach.toml > beach-run.log 2>&1
run_status=$?
printf 'exit_code=%s\n' "$run_status"
tail -n 50 beach-run.log
```

**期待:** [入力パラメータ](Parameters.html)と[場の評価](FieldSolvers.html)の互換条件を満たします。
現行実装で dielectric polarization は未実装です。

**安全な対処:** サポートされた surface model、field solver、field boundary の組合せへ設定を戻します。
未実装の物理を別のキーや数値許容値で代用しないでください。

## FMM far correction が遅い

**確認コマンド:** 計算ノード上で workload とフェーズ別時間を確認します。

```bash
beachx workload beach.toml --threads 8
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach beach.toml
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

**期待:** `cached_kneq0` を使う繰返し実行では warm cache が再利用され、cold operator 生成を毎回行いません。

**安全な対処:** cache fingerprint と cache directory を確認し、互換性のある warm cache を再利用します。
速度のために物理的な zero mode や遠方補正モデルを無断で変えず、変更時は[計算結果の妥当性確認](ValidationGuide.html)をやり直します。

## 解決しない問題を報告する

問題を再現できる場合は、次を添えて GitHub issue に報告してください。

- 使用した設定ファイルと最小 mesh
- `beach --version`、compiler / MPI、rank / thread 数
- 実行コマンド、終了コード、エラー全文
- 再開問題では `summary.txt` と `checkpoint_complete.txt`。RNG state の内容自体は不要
- 入力やパスに含まれる秘密情報を除いた再現手順
