title: FMMを使う

Lang: [日本語](FMM.md) | [English](FMM.en.md)

# FMMを使う

このページは、BEACH の利用者が FMM を選び、計算領域を設定し、Direct との誤差と実行時間を確認するための
利用手順です。FMM は、多数の三角形要素が作る遠方場を近似して、同じ mesh 上で多数の粒子位置を評価する計算を
高速化します。小規模計算では固定費が上回るため、Direct との精度比較と実測時間の両方で採用を決めてください。

読了後は、free-space FMM の最小設定、精度調整、Direct 比較、性能測定を一つのケースで実行できます。
展開式、内部データ構造、更新アルゴリズム、並列実装は
[Coulomb FMMコア内部実装](FMMCore.html)を正本とします。

## FMMを選ぶか判断する

| 条件 | 選択 |
| --- | --- |
| 要素数が多く、1 batch で多数の粒子 step を評価する | FMM を候補にし、release build で Direct と実測比較する |
| 無限 2 軸周期場を本番計算で使う | FMM と `cached_kneq0` を使う |
| 小規模ケース、基準解、近似誤差の確認 | Direct を使う |
| free-space で要素数による自動選択を使いたい | `field_solver="auto"` を使い、解決された solver を出力で確認する |

FMM の有利・不利は要素数だけでなく、粒子数、追跡 step 数、評価点の分布にも依存します。
Treecode を含む solver 全体の選択は[場の評価](FieldSolvers.html)を参照してください。

## 最小設定で自由空間FMMを使う

既に lint が通るケースの `[sim]`、`[domain]`、`[field_boundary]`、`[output]` を次のように設定します。
box の値は例なので、実際の計算領域を使ってください。

```toml
[sim]
field_solver = "fmm"

[domain]
box_min = [-0.5, -0.5, -0.1]
box_max = [ 0.5,  0.5,  1.0]
periodic_axes = []

[field_boundary]
mode = "free"

[output]
write_files = true
dir = "outputs/fmm"
```

完全なケースは [`examples/panel_fmm.toml`](../examples/panel_fmm.toml) にあります。この example は
短い確認用に `write_files=false` なので、結果を調べる場合は上の `[output]` 設定へ変更してください。
入力検査、実行、solver の確認を順に行います。

```bash
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/fmm
grep '^field_reconstruction_resolved_field_solver=' outputs/fmm/summary.txt
```

lint の末尾が `status=ok`、実行の終了コードが `0` になり、最後のコマンドが
`field_reconstruction_resolved_field_solver=fmm` を表示すれば FMM 経路で完了しています。
これは実行完了の確認であり、近似誤差や物理的妥当性の合格判定ではありません。

## domainで評価範囲を覆う

`[domain]` は粒子領域だけでなく、FMM が高速評価する target tree の範囲も定めます。
粒子が到達し得る位置をすべて box 内に入れてください。要素中心の電位や履歴を評価する場合は、mesh も
box 内に収めると範囲外評価を避けられます。

box は `box_min` と `box_max`、または `box_origin` と正の `box_size` のどちらか一組で指定します。
free 場で範囲外の点を評価すると、BEACH は全要素の Direct 和へ fallback します。結果は得られますが、頻発すると
FMM の速度上の利点が失われます。`cached_kneq0` を使う periodic2 計算では範囲外 target を受理しません。

## 精度を調整する

最初は `tree_theta` と `tree_leaf_max` を省略し、要素数から解決される値を基準にします。
この自動調整は、明示的に `field_solver="fmm"` を選んだ場合にも適用されます。解決値は次で確認できます。

```bash
grep -E '^field_reconstruction_(tree_theta|tree_leaf_max|fmm_expansion_order)=' \
  outputs/fmm/summary.txt
```

| 設定 | 調整時の見方 |
| --- | --- |
| `tree_theta` | 小さくすると遠方近似の採用条件が厳しくなり、一般に高精度・低速になる。`0 < theta <= 1` |
| `tree_leaf_max` | 木の深さと近傍 Direct 計算量のバランスを変える。`1` 以上で複数値を比較する |
| `field_normalization` | 内部座標の数値スケールだけを変える。物理単位や FMM 近似の合格条件には使わない |

現行 simulator の展開次数は 4 固定で、入力から変更できません。次の順で調整します。

1. 自動解決値で基準実行する。
2. `tree_theta` を小さくした実行を追加する。
3. `tree_leaf_max` も変え、目的の観測量が安定する範囲を探す。
4. Direct との差が許容範囲内で、かつ実測時間が短い組合せを採用する。

要素数別の自動解決値と全キーの型・既定値は[入力パラメータ](Parameters.html#場ソルバ)を参照してください。

<a id="source-kernel"></a>

## Directと比較して精度を確認する

Direct と FMM は、どちらも要素電荷を同じ triangle P0 source として扱います。そのため、同じ mesh での差は
主に FMM の近似誤差を表します。mesh 離散化誤差は、別に mesh を細分化して確認してください。

自由空間の縮小ケースを複製し、solver と出力先だけを分けます。最初の solver 誤差比較では
`sim.batch_count=1` とし、初期表面電荷が同じ状態から得た要素電荷を使います。

| 設定ファイル | `sim.field_solver` | `output.dir` |
| --- | --- | --- |
| `direct.toml` | `"direct"` | `"outputs/direct"` |
| `fmm.toml` | `"fmm"` | `"outputs/fmm"` |

mesh、粒子、乱数 seed、thread 数、`dt`、`field_normalization`、出力条件は一致させます。

```bash
beachx lint direct.toml
beachx lint fmm.toml
OMP_NUM_THREADS=1 beach direct.toml
OMP_NUM_THREADS=1 beach fmm.toml
beachx inspect outputs/direct
beachx inspect outputs/fmm
```

同じ評価点での電場と電位は、[電場計算](PythonPostprocessAPI.html#6-電場計算)と
[電位再構成](PythonPostprocessAPI.html#4-電位再構成)を使って比較できます。
基準場がほぼ 0 の点では相対誤差が不安定になるため、絶対誤差または代表場で正規化した誤差も使ってください。
solver 誤差が合格した後に代表的な複数 batch へ延長し、吸収・escape 数、要素別の最終電荷、
`charge_ledger.csv` も比較します。許容値は研究目的に応じて先に決め、
[計算結果の妥当性確認](ValidationGuide.html)に従って判定します。

panel 面上ちょうどの電場は FMM と Direct で trace の規約が異なるため、通常の点ごとの一致判定に使いません。
periodic2 の Direct 基準は solver 名だけを変更して作れません。
[periodic2 静電場](PeriodicElectrostatics.html)の split reference を使用してください。

## 性能を測る

精度条件を満たした後、release build と実運用に近い粒子数で測定します。Direct と FMM で MPI rank、OpenMP thread、
mesh、粒子、出力条件を揃え、ローカル環境または計算ノードの割当内で実行してください。

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach direct.toml
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach fmm.toml
beachx profile outputs/direct/performance_profile.csv \
  --save outputs/direct/performance_profile.png
beachx profile outputs/fmm/performance_profile.csv \
  --save outputs/fmm/performance_profile.png
```

全体比較には `simulation_total` 行の `rank_max_s` を使います。FMM の固定費は `field_solver_init`、
batch ごとの場更新は `field_refresh`、粒子評価を含む追跡は `particle_batch` で確認できます。
小ケースでは初期化の固定費が Direct を上回る場合があります。

`write_potential_history` などの出力条件も計算量を変えるため、比較間で一致させます。periodic2 の
`cached_kneq0` は cold cache と warm cache を分けて測定してください。計測項目の定義は
[実行時の性能計測](Execution.html#性能計測)を参照してください。

## periodic2でFMMを使う

`periodic2` の通常の本番経路は FMM です。無限周期の非零 mode には `cached_kneq0` を使い、
物理的 zero mode を `[periodic2]` で明示します。`field_periodic_far_correction="none"` は有限画像近似であり、
無限周期解ではありません。

完全な構成は [`examples/periodic2_cached_panel.toml`](../examples/periodic2_cached_panel.toml)を参照してください。
成分の選択と検証は[periodic2 静電場](PeriodicElectrostatics.html)、cache の生成と再利用は
[periodic2 遠方補正](PeriodicFarCorrection.html)を正本とします。

## 現行制約を確認する

- FMM が扱う source は Coulomb の triangle P0 要素で、別の kernel は選べない。
- simulator の展開次数は 4 固定で、次数収束は入力から実行できない。
- 対応する場境界は `free` と `periodic2`。完全な互換表は[場の評価](FieldSolvers.html#solverと場境界の互換表)に従う。
- mesh geometry は実行中に固定される。電荷だけが batch ごとに更新される。
- `field_normalization` を変更しても、出力される電場と電位は SI 単位へ戻される。
- FMM 自体は外部 plasma / sheath を解かない。matching-plane response は FMM の外側で合成される。

FMM の数式、内部 API、geometry と電荷状態の分離、範囲外 fallback、OpenMP 実装を調べる開発者は
[Coulomb FMMコア内部実装](FMMCore.html)へ進んでください。
