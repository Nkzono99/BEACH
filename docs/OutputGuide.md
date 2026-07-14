title: 出力の読み方

Lang: [日本語](OutputGuide.md) | [English](OutputGuide.en.md)

# 出力の読み方

初回実行後は、`outputs/latest/`の完了状態、電荷収支、要素電荷、履歴の順に確認します。
詳細な後処理コマンドは [後処理チュートリアル](PostprocessTutorial.html)、全パラメータは [入力パラメータリファレンス](Parameters.html) から確認できます。

## まず確認すること

最初に、計算が定義されたbatchまで実行され、必要な出力を書き出したことを確認します。完了を確認した後、
[計算結果の妥当性確認](ValidationGuide.html)に従って物理モデルと離散化を評価します。

1. `outputs/latest/summary.txt` がある。
2. `summary.txt` の `batches` が `sim.batch_count` に到達している。
3. `charges.csv` があり、要素ごとの最終電荷が出ている。
4. `beachx inspect outputs/latest` がエラーなく概要を表示する。

```bash
beachx inspect outputs/latest
```

`output.dir` を変更している場合は、`outputs/latest` を設定した出力ディレクトリに置き換えてください。

## 主要ファイル

| ファイル | いつ出るか | まず見る内容 |
| --- | --- | --- |
| `summary.txt` | 常時 | batch 数、吸収数、脱出数、最後の相対変化、MPI rank 数 |
| `outer_plasma_profile.csv` | `kinetic_1d` / `unified_linear_response`でouter stateが有効 | 収束した outer grid の `z, phi, E, rho`。条件付きcheckpoint |
| `photoelectron_histogram.csv` | `photoelectron_closure="individual_return"` | 前batchと累積のphotoelectron histogram。条件付きcheckpoint |
| `charges.csv` | 常時 | 要素ごとの最終電荷 |
| `mesh_triangles.csv` | 常時 | 三角形頂点、要素 ID、`mesh_id` |
| `mesh_sources.csv` | 常時 | OBJまたはtemplateの`mesh_id`、source kind、surface model、要素数 |
| `mesh_potential.csv` | `output.write_mesh_potential=true` | 最終時刻の要素重心電位 |
| `charge_history.csv` | `output.history_stride > 0` | batch ごとの要素電荷履歴 |
| `potential_history.csv` | `write_potential_history=true` かつ `history_stride > 0` | batch ごとの要素重心電位履歴 |
| `performance_profile.csv` | `BEACH_PROFILE=1` | phase 別実行時間 |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | 常時 | serialまたはMPI rank別の再開用乱数状態 |
| `macro_residuals.csv` | マクロ粒子残差stateが有効 | MPIでも1個だけ書くglobalな端数マクロ粒子数 |
| `charge_ledger.csv` | 常時 | 粒子種別の電荷収支と粒子数 |

生成条件と再開時の役割は、機械可読な[`beach.output-manifest.json`](../schemas/beach.output-manifest.json)を正本とします。

## 粒子種別の電荷収支

`charge_ledger.csv`は、粒子種別ごとの電荷と個数の流入・移送・流出をまとめた電荷収支fileです。
主な項目は次のとおりです。

| 項目 | 意味 |
| --- | --- |
| `injected_from_remote` | `volume_seed`または`reservoir_face`から入った量 |
| `emitted_from_surface` | `photo_raycast`で表面から出たtracked量 |
| `absorbed_on_surface` | meshへ吸収された量 |
| `escaped_to_infinity` | open/outer modelで無限遠へ脱出した量 |
| `discarded_unresolved` | `max_step`で未解決のまま破棄した量 |
| `interface_outward_gross` / `returned_gross` | outer ownership面を横切った往路・復路の総量 |

電荷保存残差は、surface、local flight、outer flight、unresolved stockのbatch前後差に、remote injection、escape、
discardを合わせて計算します。表面放出と表面吸収はsurface/flight stock間の内部移送なので、独立した外部sourceとして
二重には数えません。

粒子種別間で正負が相殺し得る`charge_ledger_residual_C`とは別に、
`charge_ledger_discarded_unresolved_abs_C`は$\sum_s|Q_{s,\mathrm{discard}}|$を示します。
保存残差が小さくても、未解決のまま破棄した電荷が大きい計算は受理しません。

## 成功と注意の読み分け

`summary.txt` で最初に確認する量は次のとおりです。

| 項目 | 見方 |
| --- | --- |
| `batches` | 通常実行では `sim.batch_count` と一致していれば完了 |
| `absorbed` | 表面に吸収された粒子数。帯電が進んでいるかを見る主指標 |
| `escaped` | open boundary から出た粒子数。注入・境界条件の確認に使う |
| `survived_max_step` | `sim.max_step` まで生存した粒子数。多い場合は `dt`、箱、注入条件を見直す |
| `last_rel_change` | 最終 batch の電荷変化監視値。現行実装では早期停止条件ではない |
| `charge_ledger_residual_C` | 電荷保存残差。0 でも unresolved discard は別途確認する |
| `charge_ledger_discarded_unresolved_abs_C` | species 間で相殺しない max-step discard 電荷の絶対値和 |

### モデル固有の診断

`field_source_model`と`field_kernel_id`は、出力の計算に使ったfield kernelを示します。
`triangle_p0_exact_tree_near`は、全頂点を含むnode半径、解析的なpanel near評価、monopole far評価を使うTreecodeです。
`triangle_p0_exact_p2m_near`は、全頂点のtopology、解析的なpanel near評価、厳密なpanel P2Mを使うFMMです。
`FieldKernel.from_result`は`triangle_p0`をpanel C ABIにdispatchします。その他のPython側のpotential/field/force/field-line estimatorは
point sourceのみに対応するため、`triangle_p0`の結果では停止します。

split periodic2では、`summary.txt`にinterface potential/normal field、`eta_phi_kneq0`、`eta_field_kneq0`、
`eta_gap`、`eta_local_charge`、Gauss residual、outer積分電荷、最後にouterを更新したbatchを保存します。
これらは物理modelの適用性を判定する診断値であると同時に、restart状態の一部でもあります。
そのため、これらの値が欠けたsplit checkpointからは再開できません。

`unified_linear_response`では、`outer_accessible_fraction_min/max`と
`outer_accessible_fraction_refinement_error`も確認します。後者は、rough surfaceの高さ標本数を
各周期軸で2倍に増やしたときの、accessible fractionの最大差です。この値が
`outer_plasma.accessible_fraction_tolerance`を超えると、初期化時に停止します。

`cached_kneq0` では次のcache診断を確認します。

| summary key | cold miss | warm reuse |
| --- | --- | --- |
| `periodic2_cache_hit` | `F` | `T` |
| `periodic2_operator_build_count` | root rankで`1` | `0` |
| `periodic2_cache_fingerprint` | 生成identity | 再利用identityと一致 |
| `periodic2_cache_path` | 公開先 | 読み込み元 |

cold buildでは、target sliceをMPI rank間に分配し、proxy RHSを各rankのOpenMP thread間に分配します。
cache I/Oはroot rankだけが担当します。operator本体は再生成できるため、checkpointには含めません。
cache directoryとgeneration toleranceは`summary.txt`に保存します。

particle transferを有効にすると、`charge_ledger.csv`の`interface_outward_gross_C`と
`interface_returned_gross_C`にinterfaceを通過する往路・復路の電荷量を記録します。これらを保存残差に二重加算することはありません。
`summary.txt`の`max_outer_flight_time_s`、`max_outer_frozen_field_ratio`、
`max_outer_energy_relative_error`は、run全体からMPI集約した最大値です。

`conductor` や `dielectric` を含む mesh では、`summary.txt` に注意書きが出る場合があります。
`dielectric` は現行実装ではメタデータであり、誘電体境界条件を解くモデルではありません。

## 履歴を使う

時間発展を見るには、`output.history_stride` を正の値にします。

```toml
[output]
history_stride = 1
write_potential_history = true
```

`charge_history.csv`は

$$
(\texttt{stats.batches}-1)\bmod\texttt{history\_stride}=0
$$

となるbatchで書かれ、batch 1は常に含まれます。`write_potential_history=true`なら、同じstrideで現在の
`q_elem`から電場・電位を更新し、要素重心の`potential_history.csv`も書きます。

電位履歴は追加の場評価を行うため、要素数や履歴頻度が大きい場合は計算コストが増えます。
最初は `history_stride = 1` で小ケースを確認し、大きい計算では間引いてください。

## 次に実行するコマンド

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif
```

Python から読む場合:

```python
from beach import Beach

b = Beach("outputs/latest")
print(b.result.absorbed, b.result.escaped)
b.plot_mesh()
b.plot_potential()
```

## 再開実行の出力

`output.resume=true`では、設定と保存済みstateに応じて次のファイルを読み込みます。

| checkpointファイル | 再開時の扱い |
| --- | --- |
| `summary.txt`, `charges.csv` | 常に必須 |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | serialでは前者、MPIでは全rank分が必須 |
| `macro_residuals.csv` | 存在する場合にglobal残差を復元。MPIでも単一ファイル。旧rank別ファイルは拒否 |
| `charge_ledger.csv` | summaryにledger checkpoint metadataがある場合に必要 |
| `outer_plasma_profile.csv` | 保存済みouter stateがreadyで、`kinetic_1d` / `unified_linear_response`を再開するときに必須 |
| `photoelectron_histogram.csv` | `photoelectron_closure="individual_return"`を再開するときに必須 |

schema v3ではmodel / mesh / species fingerprintを照合し、outer solverのprofile/stateを完全に復元します。

schema v2の3列outer profileも読み込めます。ただし、held stateとしては使わず、次のrefreshでouter profileを解き直します。
`output.restart_from`を指定すると、checkpointは`restart_from`から読み込み、新しい出力は`output.dir`に書き出します。
