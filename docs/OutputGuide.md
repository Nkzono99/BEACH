title: 出力ファイルを調べる

# 出力ファイルを調べる

Lang: [日本語](OutputGuide.md) | [English](OutputGuide.en.md)

> **質問:** 公式チュートリアルは正常に終わり、表面電荷は保存則と整合していますか。
>
> **一文での回答:** まず `beachx inspect` で実行件数を確認し、次に `charges.csv` と
> `charge_ledger.csv` で最終表面電荷と粒子の行き先を照合します。

読了後には、公式ケースの実行完了を判定し、最終電荷、mesh、species 別電荷収支、履歴、電位を
順に確認できます。モデル固有の receipt、完全な列定義、checkpoint schema は
[出力形式リファレンス](OutputReference.html)に分離しています。

## 公式入門ケースの結果を確認する

[10 分チュートリアル](Tutorial.html)を実行したディレクトリで、`beach.toml` の
`output.dir` に指定された結果を調べます。公式入門ケースでは `outputs/tutorial` です。

```bash
beachx inspect outputs/tutorial
```

チュートリアルどおり 1 OpenMP thread で実行すると、少なくとも次の行が含まれます。

```text
directory=outputs/tutorial
processed_particles=4000
absorbed=3720 escaped=280
batches=20 last_rel_change=...
charge_sum=-1.192019e-10
```

最小合格条件は次のとおりです。

- `beach` と `beachx inspect` の終了コードが `0`
- `batches=20` が `beach.toml` の `sim.batch_count=20` と一致する
- `processed_particles=4000` が 1 batch あたり 200 粒子 × 20 batch と一致する
- `outputs/tutorial/summary.txt` と `outputs/tutorial/charges.csv` が存在する

`absorbed=3720`、`escaped=280`、`charge_sum=-1.192019e-10` は、現行版、`rng_seed=12345`、
1 OpenMP thread での参照値です。thread 数や乱数実装が異なる場合は、乱数列に依存する値まで
同一であるとは仮定しません。

この判定が示すのは実行完了だけです。`last_rel_change` と `tol_rel` は自動停止条件ではなく、
数値収束や物理的妥当性は[計算結果の妥当性確認](ValidationGuide.html)で別に判定します。
コマンドが失敗する、または最小条件を満たさない場合は
[トラブルシューティング](Troubleshooting.html)へ進んでください。

## 最初に見るファイル

| 確認したいこと | 最初に見る場所 |
| --- | --- |
| 実行件数、batch 数、最終変化量 | `summary.txt` と `beachx inspect` |
| 最終的な要素電荷 | `charges.csv` |
| 三角形 geometry と mesh ID | `mesh_triangles.csv` |
| mesh ID と入力 mesh の対応 | `mesh_sources.csv` |

### `summary.txt`

`summary.txt` は、実行統計と解決済み設定を `key=value` 形式で記録します。最初は
`processed_particles`、`absorbed`、`escaped`、`batches`、`last_rel_change` を読みます。
モデル固有の key を近くの `beach.toml` から推測せず、実行時の receipt として読む場合は
[構成固有の出力](OutputReference.html#構成固有の値を探す)を参照してください。

### `charges.csv`

`charges.csv` は最終状態の `elem_idx,charge_C` を持ちます。`charge_C` は各三角形要素の総電荷 [C] で、
表面電荷密度ではありません。公式ケースでは 288 行あり、総和は `beachx inspect` の
`charge_sum=-1.192019e-10` と一致します。

```bash
head -n 3 outputs/tutorial/charges.csv
```

### mesh ファイル

`mesh_triangles.csv` は三角形頂点、要素電荷、`mesh_id` を持ちます。`mesh_sources.csv` は各 `mesh_id` を
入力 mesh、template、表面モデルへ対応付けます。公式ケースは 288 三角形からなる 1 個の plane mesh です。

```bash
head -n 3 outputs/tutorial/mesh_triangles.csv
head -n 2 outputs/tutorial/mesh_sources.csv
```

## charge ledger

`charge_ledger.csv` は、各 species の注入、表面吸収、escape、未解決 discard を signed charge と count で
集計します。公式ケースでは電子 1 species だけなので、次の短い確認で粒子数、電荷、最終表面電荷を
同時に照合できます。

```bash
python - <<'PY'
import csv
from math import fsum
from pathlib import Path

out = Path("outputs/tutorial")
with (out / "charge_ledger.csv").open(newline="", encoding="utf-8") as f:
    row = next(csv.DictReader(f))
with (out / "charges.csv").open(newline="", encoding="utf-8") as f:
    surface = fsum(float(item["charge_C"]) for item in csv.DictReader(f))

terminal_count = sum(int(row[name]) for name in (
    "absorbed_count", "escaped_count", "discarded_unresolved_count"
))
terminal_charge = fsum(float(row[name]) for name in (
    "absorbed_on_surface_C", "escaped_to_infinity_C", "discarded_unresolved_C"
))
print(f"counts: injected={row['injected_count']} terminal={terminal_count}")
print(f"charge_C: injected={float(row['injected_from_remote_C']):.12e} "
      f"terminal={terminal_charge:.12e}")
print(f"surface_C: charges={surface:.12e} "
      f"absorbed={float(row['absorbed_on_surface_C']):.12e}")
PY

grep -E '^(charge_ledger_surface_charge_after_C|charge_ledger_residual_C)=' \
  outputs/tutorial/summary.txt
```

現行版の参照結果は次の関係を満たします。

```text
counts: injected=4000 terminal=4000
charge_C: injected=-1.281741307200e-10 terminal=-1.281741307200e-10
surface_C: charges=-1.192019415696e-10 absorbed=-1.192019415696e-10
charge_ledger_surface_charge_after_C= -1.1920194156960000E-10
charge_ledger_residual_C= 5.5497729723797089E-25
```

このケースでは、注入粒子は吸収、escape、未解決 discard のいずれかへ進みます。表面放出と外部補正が
ないため、吸収電荷は `charges.csv` の総和になります。`charge_ledger_residual_C` は丸め誤差程度で 0 に
近いことを確認しますが、小さい保存残差だけでは統計収束や物理妥当性を証明できません。
補正を含む一般の保存則と全列は[charge ledger リファレンス](OutputReference.html#charge-ledger)を参照してください。

## 履歴

公式ケースは各 batch の電荷と電位を保存します。

| ファイル | 読み取るもの |
| --- | --- |
| `charge_history.csv` | batch ごとの各要素電荷と `rel_change` |
| `potential_history.csv` | batch ごとの各要素電位 |
| `top_reference_history.csv` | box の z-high 面における電位統計 |

```bash
head -n 3 outputs/tutorial/charge_history.csv
head -n 3 outputs/tutorial/potential_history.csv

beachx animate outputs/tutorial \
  --quantity charge \
  --save-gif outputs/tutorial/charge_history.gif \
  --total-frames 20
```

公式ケースでは batch 1 から 20 までの 20 snapshot があり、負電荷が蓄積する様子を確認できます。
電位履歴の再構成には native field kernel が必要です。必要な場合は
[後処理チュートリアル](PostprocessTutorial.html)の追加手順を使います。
生成条件と `matching_plane_history.csv` の完全な列は
[履歴リファレンス](OutputReference.html#履歴)を参照してください。

## mesh 電位

`mesh_potential.csv` は最終状態における各要素重心の電位 [V] です。公式ケースには 288 行あり、
1 thread の参照範囲は `potential_min=-4.671330e+00`、`potential_max=-2.579807e+00` です。
`potential_history.csv` は batch ごとの履歴、`mesh_potential.csv` は最終状態という違いがあります。

電位の基準、periodic2、field reconstruction の条件は
[mesh 電位リファレンス](OutputReference.html#mesh-電位)にまとめています。

## 再開に使うファイル

再開するときは checkpoint 内のファイルを手作業で結合せず、`output.restart_from` に読み込み元を指定します。
20 batch の公式結果から 21 batch 目へ進む操作は
[checkpoint から一度再開する](Execution.html#checkpointから一度再開する)に従ってください。

必須ファイル、periodic slot の選択、schema 互換性を調べる場合だけ
[checkpoint 出力リファレンス](OutputReference.html#再開に使うファイル)を参照します。

## 構成固有の値を探す

通常の公式ケースには、以下の詳細は必要ありません。使用したモデルまたは診断に対応する行だけを選びます。

| 調べたいもの | 参照先 |
| --- | --- |
| 全ファイルの生成条件 | [ファイル生成条件](OutputReference.html#ファイル生成条件) |
| field solver、periodic2、粒子境界の解決結果 | [構成固有の receipt](OutputReference.html#構成固有の値を探す) |
| `zhao_stationary` の signed 電流と補正 | [`zhao_stationary`](OutputReference.html#zhao_stationary) |
| matching-plane の accepted state | [`matching_plane_quasistatic`](OutputReference.html#matching_plane_quasistatic) |
| `charge_ledger.csv` の全列 | [charge ledger](OutputReference.html#charge-ledger) |
| adaptive periodic2 の trial | [適応 batch の診断](OutputReference.html#適応-batch-の診断) |
| checkpoint schema と必須ファイル | [再開に使うファイル](OutputReference.html#再開に使うファイル) |
| Python reader の対応属性 | [Python から読む](OutputReference.html#python-から読む) |

### 適応 batch の診断

`periodic2.max_nonzero_mode_potential_step > 0` の receipt と、履歴に含まれない棄却 trial は
[適応 batch の診断リファレンス](OutputReference.html#適応-batch-の診断)で確認します。

## 次へ進む

| 目的 | 次に読むページ |
| --- | --- |
| 電荷・電位を可視化する | [後処理チュートリアル](PostprocessTutorial.html) |
| checkpoint から再開する | [実行・再開する](Execution.html#checkpointから一度再開する) |
| 研究結果として受理できるか調べる | [計算結果の妥当性確認](ValidationGuide.html) |
| 値の不一致や欠落を診断する | [トラブルシューティング](Troubleshooting.html) |
| 出力契約を検索する | [出力形式リファレンス](OutputReference.html) |

### Python から読む

Python での最初の読込みと可視化は[後処理チュートリアル](PostprocessTutorial.html)、
class と関数の完全な仕様は [Python 後処理 API](PythonPostprocessAPI.html)を参照してください。
