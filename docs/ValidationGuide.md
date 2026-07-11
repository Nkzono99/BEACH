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
8. cached 無限周期実装は解析解 oracle でも確認します。軽量な pytest oracle は
   `BEACH_RUN_FIELD_KERNEL_CACHE_TESTS=1` の opt-in ですが、本番 strict chain では解析 job が
   `probe-periodic-oracles` を `analyze --require-complete` より先に必ず実行し、同じ staged
   library に結び付いた write-once receipt を要求します。本番入力と同じ
   `point` source / `softened_point` kernel を主 oracle とし、`triangle_p0` は補助 oracle として
   同じ library で評価します。
9. 一様な非中性 plane (`sigma=Q/A`) では、`E_bottom=0` closure に対して below
   `E_z=0`、surface PV `E_z=sigma/(2 epsilon0)`、above `E_z=sigma/epsilon0`、面圧
   `sigma^2/(2 epsilon0)` と全 object 力 `Q^2/(2 epsilon0 A)` を確認します。potential は sheet
   高さ `z0` で gauge を固定し、`phi=0` on/below、
   `phi=-sigma (z-z0)/epsilon0` above とします。surface PV は field trace であり、potential は
   この gauge で連続です。この力は closure 依存の Maxwell traction であり、自由空間の普遍的な
   self-force ではありません。
10. production `point` oracle は各 `z=-0.25, 0, +0.25` で24x24の x/y sampleを空間平均し、
   source grid 4x4 と8x8を評価します。off-surface modulation RMS は4から8で減少し、fine grid
   で `0.12 V/m` 以下、point-centroid wrench は両 grid で解析値から12%以内、grid間 force 差は
   1%以内を要求します。force・transverse force・torque の誤差はそれぞれ4から8で減少します。
   primary-free subtraction residual も減少して fine grid で1%以内でなければなりません。
   decomposition は other force、external force、`total_external-target_periodic_images` をそれぞれ
   `1e-12` 以下とし、報告する component residual がこの3成分と primary-free residual の最大値に
   一致することも確認します。
   別の単一 primary check は primary self-force の除去と softened self potential
   energy `-K q^2/a` の減算を相対 `1e-11` で検証します。これは
   `exclude_primary_keep_images`、すなわち primary だけを除いて周期像を保持する契約です。補助
   `triangle_p0` oracle は Gauss-Duffy order 3/7 の wrench 差を1%以内にします。
11. 中性な `sigma_0 cos(kx)` sheet では、4x4から8x8へ field と potential の解析解誤差が
   それぞれ減少し、fine grid で各8%以内、振幅が `exp(-k |z-z0|)` で減衰することを確認します。
   production point source は両 grid で charge-neutrality ratio `<=1e-12` を要求し、対になった
   `+z/-z` sample で tangential field は even、normal field は odd、potential は even として、
   field/potential の parity 誤差を別々に8%以内にします。さらに `a/L=1e-6` の softened-point
   micro-oracle は `r/a=0,1,2,3` の4点で analytic softened field/potential と direct evaluator、
   および ordinary/direct の一致を相対 `1e-11` で確認し、normalized self-field は
   `32 epsilon_machine` 以下とします。これは局所 kernel 契約の検査であり、periodic closure の
   代替ではありません。一様面の12%と cosine の8%は production ABI/cache 経路の smoke gate
   であり、object path の0.5%や finite-shell の1%という収束基準の代わりではなく、8%/12%を
   物理精度として主張しません。

## 旧 run・finite・infinite の SysA 比較手順

`tools/periodic_object_validation.py` は、既存 archive を基準に、同じ物理入力の
`finite_configured` と `infinite_physical` を比較するための fail-closed harness です。
検証 root は repository の外に置き、既存内容のない directory を指定します。

```bash
current_sys="$(module -t list 2>&1 | grep -E '^Sys(A|B|C|CL|G)/' | head -n 1 || true)"
if [ -n "${current_sys}" ] && [ "${current_sys}" != "SysA/2022" ]; then
  module switch "${current_sys}" SysA/2022
elif [ -z "${current_sys}" ]; then
  module load SysA/2022
fi
module load intel/2023.2 intelmpi/2023.2

python3.11 tools/periodic_object_validation.py stage \
  --archive-run /path/to/RUN \
  --validation-root /LARGE1/.../beach-periodic-object-force-validation \
  --binary /path/to/clean-build/beach \
  --library /path/to/clean-build/libbeach_field_kernel.so \
  --require-clean-source
python3.11 tools/periodic_object_validation.py verify-inputs \
  --validation-root /LARGE1/.../beach-periodic-object-force-validation
bash /LARGE1/.../beach-periodic-object-force-validation/submit/submit_chain.sh
```

`stage` は executable、Python source、解析 tool、native kernel library を hash 付きで snapshot
します。運用上は、同じ clean commit から SysA/Intel 環境で executable と library を build して
から `stage` に渡します。`stage --require-clean-source` は executable の `--build-info` と library の
C ABI build-info を読み、version/mode/full source SHA/`SHA:clean` が相互に一致し、staged source
commit と一致する場合だけ受理します。`verify-inputs` は staged artifact から同じ情報を再取得し、
各 simulation の `summary.txt` と plane-oracle receipt も同じ build origin に固定されます。
`analyze --require-complete` はこの build origin を無条件に要求するため、production stage で
`--require-clean-source` を省略した manifest は strict qualification を通りません。
compiler identity 自体は binary に埋め込まず、SysA/Intel の module log と build log を別の実行証拠として
保持します。

生成される DAG は smoke（cache prime を含む）、finite 140000、finite 280000、infinite
140000、infinite 280000、analysis の6 jobです。各 model は140000 batch まで新規実行し、
検証済み checkpoint から別 directory の 280000 batch segment へ restartします。最後の解析
job は有限・無限の両 280000 job に `afterok` で依存し、必須の production plane oracle の後に
`analyze --require-complete` を実行します。`submit_chain.sh` は各投入直後に atomic 更新する永続
`job_ids.json` journal と排他的 lock を使い、部分投入後を含む同一 chain の再投入を拒否します。
strict 解析は6個の一意な job ID、全 job の SysA/Intel module/hash log、先行5 job の
source commit・resource・exit code 0 の status、および実行中 analysis job の ID を再検証します。

`verify-run` は schema、geometry、charge/particle ledger、全 MPI rank checkpoint、履歴、
`mesh_sources.csv` / `mesh_potential.csv`、restart fingerprint、cache fingerprint と
cache file hash を再検証します。初回成功時に output tree の全 regular file を含む
write-once execution receipt を `provenance/verified/` に作り、以後は上書きせず現在値を
その receipt と比較します。restart parent と cache prime の receipt も依存関係として
固定します。cached run と後処理 evaluator は、cache fingerprint、path、file hash に加えて
検証済み cache-prime receipt の hash に結び付け、別 cache の黙示的な再利用を認めません。
最終解析では stage 時の archive input、manifest、source snapshot、case graph、
binary/library を再検証し、7 case 全てに既存 execution receipt があることを要求します。
新 run がない archive-only preflight では `--require-complete` を付けず、`missing` /
`not_evaluated` を結論に使いません。

比較の意味は固定されています。`archive -> new finite` は version、restart、runtime の
drift を含む再現性比較であり、境界条件だけの差ではありません。境界 model の効果として
解釈できるのは、同じ新 executable と入力で実行した `new finite -> new infinite` です。
`vdw_work` の速度換算は初期接着力と全接着仕事を保存する等価な有限 range 定数力を使い、
元の `1/s^2` 障壁形状を再現するものではありません。主結果は constrained translation の
`eta_translation=1`、旧解析の `0.5` は感度値として別に保存します。

旧解析が出した `Fz>0` 自体も別に監査します。archive の
`force_timeseries.csv`、`moving_sphere_force_curves.csv`、
`moving_sphere_release_summary.csv`、`moving_sphere_model_summary.json` を hash と exact schema で固定し、
共通の batch/mesh `(149001,7)`, `(180001,6)`, `(279001,6)`, `(279001,7)` だけを比較します。
旧重心間 direct、local-pair proxy、moving-sphere の `z=0` / `z=2R` の力と仕事を、同じ archived
charge を current native finite evaluator で再評価した other objects、自身の周期像、total に分解して
`legacy_estimator_comparison.csv` へ保存します。280000 は旧表にないため補間しません。旧 estimator は
target 全体を self source から除き、current evaluator は primary のみを除いて target periodic images を
保持するので、数値の近さは合否条件にしません。strict gate は coverage、有限値、charge/radius/batch 対応、
および旧 curve と summary の endpoint 内部整合です。

解析は全 case の共通保存 batch で per-object charge と element charge の L1/L2/Linf 差を
計算します。高価な force/path は mesh 7 の 149001、mesh 6 の 180001、両者の 279001 と
final 280000 に限定します。`snapshot_manifest.csv` は history と final を分離し、charge
vector と source file の hash を持ちます。`comparison_matrix.csv` は archive version drift、
同一 charge に対する field closure 差、共通 infinite evaluator での charging-history 差、
end-to-end 差を別の `comparison_kind` として保存します。
strict comparison artifact contract はこの4種類だけを許し、各種類に force、endpoint work、
minimum available energy、from-rest barrier、endpoint reachability を要求します。参照 snapshot
は全て解決可能で、frozen-field 比較は同一 charge snapshot、実効 far-correction の組は宣言どおりで、
構造的な `missing` / `invalid` / `not_evaluated` 行がないことを確認します。

archive/new finite/new infinite の mesh ID と順序は厳密一致を要求し、triangle 座標は
`max(1e-18 m, 64 epsilon Lbox)` の明示許容差で比較します。さらに new finite/new infinite の
`field_source_model`、`field_kernel_id` と、`mesh_sources.csv` の source/template kind、surface
model、`epsilon_r`、element count を staged input および相互で厳密一致させます。旧 archive に
同じ metadata がない場合は、その欠落を new run の契約を弱める理由にしません。cached evaluator
については実効 model に加えて cache hit、build count、fingerprint、path、file hash と
cache-prime receipt hash も `object_wrench.csv` に保存します。

`analyze --require-complete` は strict input、receipt、oracle、geometry の検証と artifact 生成を
一時 directory で行います。publish 前に失敗した場合は一時 directory を除去し、`analysis/` が
未作成または空のままなら修正後に同じ validation root で retry できます。

`numerical_qualification_for_local_frozen_model` は exact 30 path/wrench key、6 shell group、path/work
収束を含みます。finite-shell の相対 tolerance は1%、adaptive path は0.5%（force absolute
floor `1e-12 N`、work absolute floor `1e-18 J`）です。さらに force floor が peak force の
0.5%以下、瞬時離脱力 margin が `1e-12 N` を超え、barrier 判定境界からの energy margin が
energy tolerance の10%を超えることを要求します。これらの resolution gate を満たさない
barrier/speed は、積分自体が収束していても主張しません。この資格は局所 frozen-charge model
に限られます。path または shell が未収束なら barrier/speed は
`not_claimed_unqualified` であり、上方一定場が残る非中性系の `0..2R` work/speed は局所
frozen-field 指標であって無限遠 escape 量ではありません。

strict 解析は次の exact 14 artifact を一時 directory に作ります: `run_summary.csv`、
`charge_history_pair.csv`、`particle_ledger_pair.csv`、`mesh_potential_pair.csv`、
`snapshot_manifest.csv`、`object_wrench.csv`、`object_path_curves.csv`、
`object_path_summary.csv`、`finite_shell_convergence.csv`、`comparison_matrix.csv`、
`legacy_estimator_comparison.csv`、`analysis_summary.json`、`review_ja.md`、`artifacts.json`。
file set、非空、size/hash manifest を
検証して fsync してから `analysis/` へ atomic publish します。strict success には publish 後に
`numerical_qualification_for_local_frozen_model.status="qualified"` と
`comparison_artifact_contract.status="complete"` の両方を要求します。どちらかが失敗すれば CLI は
非ゼロですが、publish 済み診断は上書きしません。publish 後に strict 解析を再実行するには
新しい validation root を使います。

exact 14 artifact が完全に読めることはレベル1の実行証拠です。strict CLI の exit code 0 は
上記の構造・oracle・comparison・局所数値 gate を全て通過したことを示しますが、それでも
レベル2の局所 model qualification です。model 選択と感度解析まで完了して初めてレベル3の
主張を検討できます。

release用の小規模基準は[Physics release verification](PhysicsReleaseVerification.html)にあります。
本番caseでは同じ収束軸を自分の観測量に対して評価してください。
