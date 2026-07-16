title: 出力の読み方

Lang: [日本語](OutputGuide.md) | [English](OutputGuide.en.md)

# 出力の読み方

初回実行後は、`outputs/latest/` の完了状態、電荷収支、要素電荷、履歴の順に確認します。
詳細な後処理コマンドは [後処理チュートリアル](PostprocessTutorial.html)、全パラメータは
[入力パラメータリファレンス](Parameters.html) から確認できます。

## まず確認すること

このページの確認手順は、`output.write_files = true` を前提とします。`false` の場合、BEACH は当該実行の
`summary.txt` や `charges.csv` などを `output.dir` に作成・更新しないため、`beachx inspect` でその実行を確認できません。
既存ファイルが残っていても過去の結果です。この場合は、端末に表示される実行概要だけを確認できます。

最初に、計算が設定したバッチ数まで実行され、必要な出力を書き出したことを確認します。完了を確認した後、
[計算結果の妥当性確認](ValidationGuide.html) に従って物理モデルと離散化を評価します。

1. `outputs/latest/summary.txt` がある。
2. `summary.txt` の `batches` が `sim.batch_count` に到達している。
3. `charges.csv` があり、要素ごとの最終電荷が出ている。
4. `beachx inspect outputs/latest` がエラーなく概要を表示する。

```bash
beachx inspect outputs/latest
```

`output.dir` を変更している場合は、`outputs/latest` を設定した出力先に置き換えてください。

## 主要ファイル

| ファイル | いつ出るか | まず見る内容 |
| --- | --- | --- |
| `summary.txt` | 常時 | バッチ数、吸収数、脱出数、最後の相対変化、MPI ランク数 |
| `outer_plasma_profile.csv` | `kinetic_1d` / `unified_linear_response` で外部状態が有効 | 収束した外部グリッドの `z, phi, E, rho`。条件付きチェックポイント |
| `photoelectron_histogram.csv` | `outer_plasma.photoelectron_histogram_enabled=true` でヒストグラム状態が準備済み | 前バッチと累積の光電子ヒストグラム。条件付きチェックポイント |
| `charges.csv` | 常時 | 要素ごとの最終電荷 |
| `mesh_triangles.csv` | 常時 | 三角形頂点、要素 ID、`mesh_id` |
| `mesh_sources.csv` | OBJ またはテンプレートメッシュ | `mesh_id`、`source_kind`、`surface_model`、要素数 |
| `mesh_potential.csv` | `output.write_mesh_potential=true` | 最終時刻の要素重心電位 |
| `charge_history.csv` | `output.history_stride > 0` | バッチごとの要素電荷履歴 |
| `potential_history.csv` | `output.write_potential_history=true` かつ `output.history_stride > 0` | バッチごとの要素重心電位履歴 |
| `performance_profile.csv` | `BEACH_PROFILE=1` | フェーズ別実行時間 |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | 常時 | シリアルまたは MPI ランク別の再開用乱数状態 |
| `macro_residuals.csv` | マクロ粒子残差状態が有効 | MPI でも 1 個だけ書く全体の端数マクロ粒子数 |
| `charge_ledger.csv` | 電荷台帳状態が有効 | 粒子種別の電荷収支と粒子数 |

表中のファイルは、すべて `output.write_files = true` の場合に限って生成されます。追加の生成条件と再開時の役割は、
機械可読な [`beach.output-manifest.json`](../schemas/beach.output-manifest.json) を正本とします。

## 粒子種別の電荷収支

`charge_ledger.csv` は、粒子種別ごとの電荷と個数の流入・移送・流出をまとめた電荷収支ファイルです。
主な項目は次のとおりです。

| 項目 | 意味 |
| --- | --- |
| `injected_from_remote_C` | `volume_seed` または `reservoir_face` から入った電荷量 |
| `emitted_from_surface_C` | `photo_raycast` で表面から出た追跡対象の電荷量 |
| `absorbed_on_surface_C` | メッシュへ吸収された電荷量 |
| `escaped_to_infinity_C` | 開放境界または外部モデルから無限遠へ脱出した電荷量 |
| `discarded_unresolved_C` | `max_step` 到達時に生存したまま破棄した電荷量 |
| `interface_outward_gross_C` / `interface_returned_gross_C` | 外部領域との界面を横切った往路・復路の総電荷量 |

電荷保存残差は、表面、局所飛行、外部飛行、未解決粒子の電荷ストックのバッチ前後差に、外部からの注入、
無限遠への脱出、未解決粒子の破棄を合わせて計算します。表面放出と表面吸収は、表面と飛行中粒子の間の内部移送なので、独立した外部源として
二重には数えません。

粒子種別間で正負が相殺し得る `charge_ledger_residual_C` とは別に、
`charge_ledger_discarded_unresolved_abs_C` は $\sum_s|Q_{s,\mathrm{discard}}|$ を示します。
保存残差が小さくても、未解決のまま破棄した電荷が大きい計算は受理しません。

## 成功と注意の読み分け

`summary.txt` で最初に確認する量は次のとおりです。

| 項目 | 見方 |
| --- | --- |
| `batches` | 通常実行では `sim.batch_count` と一致していれば完了 |
| `absorbed` | 表面に吸収されたマクロ粒子のイベント数。電荷の符号や粒子重みは反映しない |
| `escaped` | 開放境界から出た粒子数。注入・境界条件の確認に使う |
| `survived_max_step` | `sim.max_step` まで生存した粒子数。多い場合は `dt`、箱、注入条件を見直す |
| `last_rel_change` | 最終バッチの電荷変化監視値。現行実装では早期停止条件ではない |
| `charge_ledger_residual_C` | 電荷保存残差。0 でも未解決粒子の破棄は別途確認する |
| `charge_ledger_discarded_unresolved_abs_C` | 粒子種間で相殺しない、最大ステップ到達時の破棄電荷の絶対値和 |

`absorbed` だけでは帯電量を判断できません。符号・重み付き吸収電荷は `absorbed_on_surface_C`、表面電荷の正味変化は
`charge_ledger_surface_charge_before_C` と `charge_ledger_surface_charge_after_C` の差、最終分布は `charges.csv` で確認します。

### モデル固有の診断

`field_source_model` と `field_kernel_id` は、出力の計算に使った場カーネルを示します。
`triangle_p0_exact_tree_near` は、全頂点を含むノード半径、解析的なパネル近傍評価、単極子遠方評価を使う Treecode です。
`triangle_p0_exact_p2m_near` は、全頂点のトポロジー、解析的なパネル近傍評価、厳密なパネル P2M を使う FMM です。
`FieldKernel.from_result` は `triangle_p0` を panel C ABI へ振り分けます。その他の Python 側の電位・電場・力・電気力線の推定器は
点電荷源（point source）のみに対応するため、`triangle_p0` の結果では停止します。

split periodic2 構成では、`summary.txt` に `interface_potential_V`、`interface_normal_field_V_m`、
`interface_eta_phi_kneq0`、`interface_eta_field_kneq0`、`interface_eta_gap`、`interface_eta_local_charge`、
`gauss_residual_C`、`outer_integrated_charge_C`、`last_outer_update_batch` を保存します。
これらは、界面の電位・法線電場、ガウス則の残差、外部領域の積分電荷、更新時点を記録し、
物理モデルの適用性を判定する診断値であると同時に、再開状態の一部でもあります。
そのため、これらの値が欠けた split チェックポイントからは再開できません。

光電子histogram stateがreadyな場合は、次の`summary.txt`キーも確認します。

| キー | 内容 |
| --- | --- |
| `photoelectron_histogram_bins` | 法線運動エネルギーbin数 |
| `photoelectron_histogram_energy_max_J` | histogram上端 [J] |
| `photoelectron_last_completed_batch` | histogramへcommit済みの最終batch |
| `photoelectron_cumulative_signed_charge_C` | z-high outward interface crossingの累積signed charge [C] |
| `photoelectron_cumulative_kinetic_energy_J` | 累積全運動エネルギー [J] |
| `photoelectron_cumulative_count` | 累積crossing数 |
| `photoelectron_previous_signed_current_A` | 前batchのsigned chargeを`batch_duration`で割った電流 [A] |
| `photoelectron_previous_charge_ratio` | 前batchのoutward crossing chargeとambient charge scaleの比 |
| `photoelectron_max_charge_ratio` | 設定した適用性上限 |
| `photoelectron_linear_applicability_status` | 正常完了したbatchでは`applicable` |

`photoelectron_previous_charge_ratio`が上限を超えたbatchは停止するため、適用外の状態を正常な出力として継続しません。

`unified_linear_response` では、`outer_accessible_fraction_min`、`outer_accessible_fraction_max`、
`outer_accessible_fraction_refinement_error` も確認します。後者は、粗い表面の高さ標本数を
各周期軸で 2 倍に増やしたときの、アクセス可能率（accessible fraction）の最大差です。この値が
`outer_plasma.accessible_fraction_tolerance` を超えると、初期化時に停止します。

`cached_kneq0` では次のキャッシュ診断を確認します。

| `summary.txt` のキー | 初回構築 | 再利用 |
| --- | --- | --- |
| `periodic2_cache_hit` | `F` | `T` |
| `periodic2_operator_build_count` | ルートランクで `1` | `0` |
| `periodic2_cache_fingerprint` | 生成した識別子 | 再利用した識別子と一致 |
| `periodic2_cache_path` | 公開先 | 読み込み元 |

コールドビルドでは、対象スライスを MPI ランク間に分配し、proxy RHS を各ランクの OpenMP スレッド間に分配します。
キャッシュ I/O はルートランクだけが担当します。演算子本体は再生成できるため、チェックポイントには含めません。
キャッシュディレクトリと生成許容差は `summary.txt` に保存します。

粒子移送を有効にすると、`charge_ledger.csv` の `interface_outward_gross_C` と
`interface_returned_gross_C` に界面を通過する往路・復路の電荷量を記録します。これらを保存残差に二重加算することはありません。
`summary.txt` の `max_outer_flight_time_s`、`max_outer_frozen_field_ratio`、
`max_outer_energy_relative_error` は、実行全体から MPI 集約した最大値です。

`dielectric` を含むメッシュでは、`summary.txt` に `surface_model_dielectric_elem_count` と
`surface_model_note=metadata_only_dielectric_present` を出力します。`dielectric` は現行実装ではメタデータであり、
誘電体境界条件を解くモデルではありません。`conductor` だけを含む場合は、この注意書きを出力しません。

## 履歴を使う

時間発展を見るには、`output.history_stride` を正の値にします。

```toml
[output]
history_stride = 1
write_potential_history = true
```

`charge_history.csv` は

$$
(\texttt{stats.batches}-1)\bmod\texttt{history\_stride}=0
$$

となるバッチで書かれ、バッチ 1 は常に含まれます。`output.write_potential_history=true` なら、同じストライドで現在の
`q_elem` から電場・電位を更新し、要素重心の `potential_history.csv` も書きます。

電位履歴は追加の場評価を行うため、要素数や履歴頻度が大きい場合は計算コストが増えます。
最初は `history_stride = 1` で小ケースを確認し、大きい計算では間引いてください。

## 次に実行するコマンド

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh.png

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif
```

`field_source_model="point"` の結果に限り、Python 側で電位を再構成してメッシュへ描画できます。

```bash
beachx inspect outputs/latest \
  --save-potential-mesh outputs/latest/potential_mesh.png
```

`triangle_p0` の結果では `--recompute-potential`、`--save-potential-mesh`、`--show` が停止します。
要素重心電位が必要な場合は、実行前に
`output.write_mesh_potential = true` を設定し、BEACH が書き出す `mesh_potential.csv` を使用してください。

Python から読む場合:

```python
from beach import Beach

b = Beach("outputs/latest")
print(b.result.absorbed, b.result.escaped)
b.plot_mesh()

if b.result.field_source_model == "point":
    b.plot_potential()
```

`Beach.plot_potential()` も点電荷源専用です。

## 再開実行の出力

`output.resume=true` では、設定と保存済み状態に応じて次のファイルを読み込みます。

| チェックポイントファイル | 再開時の扱い |
| --- | --- |
| `summary.txt`, `charges.csv` | 常に必須 |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | シリアルでは前者、MPI では全ランク分が必須 |
| `macro_residuals.csv` | 存在する場合に全体の残差を復元。MPI でも単一ファイル。旧ランク別ファイルは拒否 |
| `charge_ledger.csv` | `summary.txt` に台帳のチェックポイントメタデータがある場合に必要 |
| `outer_plasma_profile.csv` | 保存済み外部状態が準備済みで、`kinetic_1d` / `unified_linear_response` を再開するときに必須 |
| `photoelectron_histogram.csv` | `outer_plasma.photoelectron_histogram_enabled=true` の状態を再開するときに必須 |

スキーマ v3 ではモデル、メッシュ、粒子種のフィンガープリントを照合し、外部ソルバーのプロファイルと状態を完全に復元します。

スキーマ v2 の 3 列外部プロファイルも読み込めます。ただし、保持状態としては使わず、次の更新時に外部プロファイルを解き直します。
`output.restart_from` を指定すると、チェックポイントは `restart_from` から読み込み、新しい出力は `output.dir` に書き出します。
