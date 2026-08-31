title: `batch_duration` をどう決めるか

Lang: [日本語](BatchDurationStability.md) | [English](BatchDurationStability.en.md)

# `batch_duration` をどう決めるか

このページは、「1 batch の物理時間をどこまで大きくしてよいか」に答える利用ガイドです。
単独の run から安全な値を決めることはできません。固定幅を 1/2 倍、1 倍、2 倍に変え、
同じ物理時刻で結果が変わらない範囲を探すのが基本です。

読了後は、固定幅の比較を実行し、幅を採用するか、`cached_kneq0` の適応進行を使うかを判断できます。

> `sim.tol_rel` は監視・出力値です。現行実装は `tol_rel` を早期停止条件に使わず、
> `sim.batch_count` で指定した accepted batch 数まで実行します。

## まず選ぶ経路

| 目的 | 使う経路 |
| --- | --- |
| 一般のケースで時間幅依存性を調べる | 固定幅の 1/2 倍、1 倍、2 倍を比較する |
| `cached_kneq0` で 1 batch の局所電位変化を制限する | `max_nonzero_mode_potential_step` を使う |
| $k=0$ を含む全体の安定性や大域精度を確認する | 固定幅比較を行う。適応進行だけでは判定しない |

`sim.batch_duration` は 1 batch の物理時間と、表面電荷を更新する時間幅です。
`sim.batch_duration_step` を使う場合、解決される幅は

$$
\texttt{batch\_duration}=\texttt{dt}\times\texttt{batch\_duration\_step}
$$

です。2 つのキーは同時に指定できません。`sim.dt` は粒子を前進させる時間刻みであり、
表面電荷をまとめて反映する `batch_duration` とは役割が異なります。

`boundary_inflow`、`plane_source`、`reservoir_face`、`photo_raycast`、および
`surface_charge_closure="fixed_current"` は、正の解決済み `batch_duration` を必要とします。
完全な入力条件は[入力パラメータリファレンス](Parameters.html)で確認してください。

## 固定幅を比較する

### 比較条件を揃える

- `write_files=true`、`history_stride > 0` とし、`summary.txt` と `charge_history.csv` を保存します。
- メッシュ、粒子分布、乱数 seed、OpenMP thread 数、MPI rank 数を揃えます。
- 3 run の `output.dir` を分けます。
- 幅だけでなく `batch_count` も変え、同じ `simulated_time_s` 付近で比較します。

基準幅を $h$、基準 batch 数を $N$ とすると、次の組合せを使います。

| run | `batch_duration` | `batch_count` | `output.dir` |
| --- | ---: | ---: | --- |
| half | $h/2$ | $2N$ | `outputs/batch-half` |
| reference | $h$ | $N$ | `outputs/batch-reference` |
| double | $2h$ | $N/2$ | `outputs/batch-double` |

$N$ は偶数にします。例えば $h=1.0\times10^{-7}$ s、$N=100$ なら、half、reference、double の
設定値はそれぞれ `(5.0e-8, 200)`、`(1.0e-7, 100)`、`(2.0e-7, 50)` です。
これは編集方法を示す例であり、推奨する物理値ではありません。
`batch_duration_step` を使う場合は `dt` を固定し、step 値を同じ 1/2 倍、1 倍、2 倍にします。

### 設定を作って実行する

元の設定を 3 つに複製し、上表に対応する `[sim]` と `[output]` の値だけを変更します。

```bash
cp beach.toml batch-half.toml
cp beach.toml batch-reference.toml
cp beach.toml batch-double.toml
```

入力検査に通ってから、ローカル環境または計算ノードで実行します。

```bash
beachx lint batch-half.toml
beachx lint batch-reference.toml
beachx lint batch-double.toml

beach batch-half.toml
beach batch-reference.toml
beach batch-double.toml
```

実行環境の選び方は[実行する](Execution.html)を参照してください。

### 期待する出力を確認する

各 run に `summary.txt`、`charges.csv`、`charge_history.csv` が必要です。
次のコマンドで、完了 batch 数と物理終了時刻を確認します。

```bash
for output_dir in outputs/batch-half outputs/batch-reference outputs/batch-double; do
  beachx inspect "$output_dir"
  grep -E '^(batches|last_rel_change|simulated_time_s)=' "$output_dir/summary.txt"
done
```

期待する状態は次のとおりです。

- 3 run が終了コード `0` で完了する
- `batches` が各設定の `sim.batch_count` と一致する
- `simulated_time_s` が 3 run で同じ物理終了時刻になる
- `charge_history.csv` の `batch` 列を各 run の幅で物理時刻へ換算し、対応する要素電荷を比較できる

正常終了は、数値安定性や定常到達を証明しません。`last_rel_change` も停止判定ではなく、
履歴と最終状態を比較するための診断値です。

### 採用する幅を決める

最終表面電荷分布、総電荷、局所電位幅、absorbed / escaped 数、履歴の振動を比較します。
「一致」の許容値は、研究で必要な目的量の精度から先に決めてください。

| 観測 | 判断 |
| --- | --- |
| half と reference が許容値内で一致する | reference は実用上の候補。half を検証基準として残す |
| reference または double で振動・発散する | `batch_duration` を下げる |
| 幅を変えると最終電荷が系統的に変わる | より小さい幅でもう一度比較する |
| 履歴が Monte Carlo ノイズに埋もれる | `w_particle` または `target_macro_particles_per_batch` を先に調整する |
| 終了時にも変化が続く | 物理終了時刻を揃えたまま `batch_count` を増やす |

この比較は Richardson 外挿ではなく、step-size sensitivity check です。誤差の冪乗則や
特定の収束次数は仮定しません。理論的な理由は
[`batch_duration` の理論的背景](BatchDurationTheory.html)で説明します。

## 適応的な $k\ne0$ 進行を使う

### 適用できるケース

この経路は既存の periodic2 ケースに対する高度な選択肢で、次をすべて要求します。

- `[periodic2].nonzero_mode_backend = "cached_kneq0"`
- 正の `sim.batch_duration`、または正の値へ解決される `sim.batch_duration_step`
- time-scaled な `boundary_inflow`、`plane_source`、`reservoir_face`、`photo_raycast`
- reservoir 流入と `plane_source` では、固定 `w_particle` ではなく `target_macro_particles_per_batch`

正の `npcls_per_step` を持つ `volume_seed` はこの経路で使えません。

### 上限を設定して実行する

元の periodic2 設定を複製し、`adaptive.toml` の既存の `[periodic2]` と `[output]` にある
該当値を変更します。

```bash
cp beach.toml adaptive.toml
```

```toml
[periodic2]
nonzero_mode_backend = "cached_kneq0"
max_nonzero_mode_potential_step = 1.0e-2 # V

[output]
dir = "outputs/adaptive"
```

`1.0e-2` V は入力例です。許容できる局所電位変化に合わせて決め、上限を半分にした比較も行います。

```bash
beachx lint adaptive.toml
beach adaptive.toml
beachx inspect outputs/adaptive
```

解決済み `sim.batch_duration` を $h_0$ とすると、BEACH は各 accepted batch で
$h_0,h_0/2,h_0/4,\ldots$ を試します。候補電荷が作る $k\ne0$ 電位変化を全 panel 重心で評価し、
最大絶対値が上限以下となる最初の幅を受理します。

### 期待する出力を確認する

```bash
grep -E '^(batches|simulated_time_s|periodic2_max_nonzero_mode_potential_step_V|adaptive_nonzero_mode_rejected_trials|adaptive_nonzero_mode_last_batch_duration_s|adaptive_nonzero_mode_last_potential_step_V|adaptive_nonzero_mode_omp_threads)=' \
  outputs/adaptive/summary.txt
```

- `adaptive_nonzero_mode_last_batch_duration_s` は最後に受理した幅です。
- `adaptive_nonzero_mode_last_potential_step_V` は最後に受理した $k\ne0$ 電位変化で、設定上限以下になります。
- `adaptive_nonzero_mode_rejected_trials=0` でも異常ではありません。上限を満たすまで半減した累積回数を示します。
- `simulated_time_s` は受理した幅の累積であり、一般には
  `batch_count * batch_duration` と一致しません。

棄却した再試行では、RNG と macro 粒子数残差を元へ戻し、統計、履歴、charge ledger を更新しません。
1 batch で 24 回半減しても上限を満たさない場合は停止します。この場合は最大幅を下げるだけでなく、
場モデル、粒子統計、電荷変化が意図どおりかも確認してください。

### 上限を採用するか判断する

1. `max_nonzero_mode_potential_step` を 1/2 にした run を作ります。
2. 同じ `simulated_time_s` 付近で表面電荷、局所電位幅、総電荷、粒子統計を比較します。
3. 設定した許容値内で一致すれば、大きい方の上限を実用候補にします。
4. 固定幅 control では、このキーを省略するか `0` にします。

`max_nonzero_mode_potential_step` は、batch 中に $k\ne0$ 場を固定するための局所電位 trust bound です。
局所打切り誤差、大域精度の次数、$k=0$ 更新の安定性は保証しません。このため、適応進行を使っても
固定幅または上限半減の sensitivity check が必要です。

target-count reservoir と固定 `rays_per_batch` の photo source では、trial 幅を半分にすると
macro 粒子電荷も半分になります。上限半減比較には時間離散化差と Monte Carlo 分散差の両方が含まれるため、
同じ乱数 seed を使い、電荷分布 norm と粒子統計を併記します。

棄却した再試行は、同じ実行内で固定 OpenMP team size を使います。restart 後は別の team size を
使用できるため、restart 前後の bitwise 一致ではなく、電荷分布と accepted width の数値的一致を確認します。

## 関連文書

- [`batch_duration` の理論的背景](BatchDurationTheory.html) — 固定点、線形安定性、適用限界
- [入力パラメータリファレンス](Parameters.html) — キーの型と完全な制約
- [境界 reservoir の流入量と速度サンプリング](ReservoirInjection.html) — 粒子数と重み
- [出力形式リファレンス](OutputReference.html#適応-batch-の診断) — adaptive receipt
- [結果を検証する](ValidationGuide.html) — 数値収束と物理妥当性
- [BEACH の計算サイクル](Algorithms.html) — batch 中に場を固定する更新順序
