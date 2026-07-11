title: 計算結果の妥当性確認

Lang: [日本語](ValidationGuide.md) | [English](ValidationGuide.en.md)

# 計算結果の妥当性確認

正常終了は、物理的・数値的に妥当であることを意味しません。確認を3段階に分けます。

## レベル1: 実行が完了した

- processの終了codeが0
- `summary.txt`、`charges.csv`、必要な履歴が存在する
- `batches == sim.batch_count`
- `beachx inspect`が読める
- restart時はmodel/mesh/species fingerprintが一致する

## レベル2: 数値的に qualification できる

- `absorbed`、`escaped_boundary`、`survived_max_step`の内訳を物理条件から説明できる
- `charge_ledger_residual_C`が丸め誤差範囲で、`discarded_unresolved`を別に確認した
- `tol_rel`を収束停止条件と誤解せず、履歴が十分な長さを持つ
- `sim.dt`を半分にして主要量が変わらない
- mesh解像度、FMM order/tolerance、outer gridを上げて主要量が収束する
- `batch_duration`を0.5倍/2倍にして結論が安定する
- stochastic caseはseedまたはensemble依存を確認する

ここでいう qualification は、宣言した離散化・収束基準を満たすという意味です。
process の終了、CSV の生成、または一つの `status="converged"` だけでは成立しません。

## レベル3: 物理的な結論を支持できる

- 比較する case 間で、意図した物理 model 以外の入力差を列挙した
- 境界条件、self interaction、surface trace、source/target の移動規則を説明できる
- 有限 box、有限時間、有限 image shell から無限遠・定常・無限周期の結論へ外挿していない
- 数値誤差、stochastic uncertainty、model uncertainty を結論の精度に反映した

## model固有の確認

| model | 必須診断 |
| --- | --- |
| periodic2 cached | cache fingerprint、cold/warm一致、zero-mode/Gauss residual |
| unified outer | accessibility refinement、linearity、outer energy/frozen-field error |
| kinetic outer | solver status、Poisson residual、Bohm/branch applicability |
| photoelectron | emission/return ledger、ambient charge ratio、histogram範囲 |
| object detachment | primary-only self exclusion、PV trace、work/potential一致、quadrature、finite-shell/cache、from-rest barrier |

## 周期 object 離脱解析の追加 gate

1. `configured` と `infinite_physical` を区別し、どちらを結論に使ったか記録します。
   `configured` は run の finite/cached 設定を再現し、`infinite_physical` は cached
   `k != 0` と `E_bottom=0` の物理 zero mode を組み合わせます。
2. self policy が `exclude_primary_keep_images` であることを確認します。旧
   `kernel-forces` の `exclude_target_lattice` や、電位再構成の `area_equivalent` と
   混同しません。
3. triangle source では order 3/7 または mesh refinement で面積積分の依存性を確認します。
   object mechanics は PV trace、粒子 pusher は片側 plus trace です。cached metadata の
   `cached_kneq0_trace_correction` は `periodic_kneq0` に適用済みなので再加算しません。
4. 凍結 source path で `integral(F_z dh)` と `U_env(0)-U_env(h)` が宣言 tolerance 内で
   一致し、`path.status="converged"` であることを確認します。外部の凍結場に対する
   `U_env=sum(q phi_env)` には係数 `1/2` を付けません。
5. finite shell は native canonical-unwrapped source 表現を使い、raw symmetric と
   `E_bottom=0` corrected の両方を保存します。隣接 shell の raw 増分だけでは false
   convergence が生じたため、`force_tail_proxy_N` / `work_tail_proxy_J` と、
   `infinite_physical` 時の `reference_force_error_N` /
   `reference_work_error_J` を組み合わせます。`increment_converged` はこの combined gate
   であり、2回連続で真になるまで採用しません。`status="not_converged"` では
   `selected_image_layers=None`, `selected_path=None` が正しい結果です。
6. `evaluate_release()` では endpoint energy だけでなく、有限 range adhesion と重力を
   含む全経路の `barrier_free_from_rest` を確認します。
7. 非中性 x/y 周期 cell では遠方の一定場・線形電位が残り得ます。有限高さの work/speed
   を無限遠への escape energy/speed として報告しません。
8. cached 無限周期実装は opt-in 解析 oracle でも確認します。一様な非中性 triangle plane
   (`sigma=Q/A`) では `E_bottom=0` に対して below `E_z=0`、above
   `E_z=sigma/epsilon0`、surface PV `E_z=sigma/(2 epsilon0)`、面圧
   `sigma^2/(2 epsilon0)` と全 object 力 `Q^2/(2 epsilon0 A)` を確認します。中性な
   `sigma_0 cos(kx)` sheet では `k != 0` の場・電位振幅が
   `exp(-k |z-z0|)` で減衰し、triangle mesh refinement で解析解誤差が減ることを確認します。

CLI が exit code 0 で4成果物を書いたことはレベル1です。`path.status` や shell 収束は
レベル2の一部であり、上記の model 選択と感度解析まで完了して初めてレベル3の主張を
検討できます。

release用の小規模基準は[Physics release verification](PhysicsReleaseVerification.html)にあります。
本番caseでは同じ収束軸を自分の観測量に対して評価してください。
