title: 外部 1D1V kinetic oracle を検証に使う

Lang: [日本語](OuterKineticOracle.md) | [English](OuterKineticOracle.en.md)

# 外部 1D1V kinetic oracle を検証に使う

`beach-kinetic-response` は、matching-plane の入力を固定し、外側の無衝突 1D1V Vlasov--Poisson 系が
定常または統計的定常へ到達するかを調べる研究用の offline tool です。BEACH runtime と時間発展を連成せず、
認証できた完全直積の応答だけを既存の `response_backend="table"` 用 CSV へ変換します。

この tool の正常終了は物理的な妥当性を意味しません。各 query の `classification`、格子収束、遠方境界、保存則を
確認してから table へ変換します。

## 1. Reference sweep を実行する

repository root で次を実行します。出力先は空または未作成でなければなりません。

```bash
beach-kinetic-response \
  examples/outer_kinetic_reference.toml \
  examples/outer_kinetic_k1_k6.csv \
  outputs/outer-kinetic-atlas
```

この例は Zhao の一意解、複数根、no-root 領域を含む診断 sweep です。全点が認証されることを期待する
production fixture ではありません。完了すると次を生成します。

| 出力 | 用途 |
|---|---|
| `kinetic_response_raw.csv` | 応答、分類、stationarity、遠方境界、保存則の query 別結果 |
| `kinetic_response_manifest.json` | 格子、速度範囲、CFL、許容値、commit、source hash、raw CSV hash |
| `kinetic_response_profiles/` | 電位・電場・電荷密度 profile と時系列 |

## 2. 分類を読む

| `classification` | table への使用 |
|---|---|
| `steady` | 格子・領域収束後に使用可能 |
| `stationary_average` | 自己相関時間と SEM を確認し、converter で明示許可した場合だけ使用可能 |
| `unresolved_transient` | 使用不可。計算時間または averaging window を見直す |
| `secular` | 使用不可。任意時刻の平均で置き換えない |
| `far_boundary_not_converged` | 使用不可。まず外部長と遠方の物理 model を検証する |
| `numerical_failure` | 使用不可。`failure_reason` と速度領域・保存則を確認する |

solver は $E(H)=D_H/\epsilon_0$ と $\Phi(L)=0$ を境界条件にします。`E(L)=0` は追加拘束ではなく、
外部長が十分かを判定する診断です。ambient の外向き入力は v1 では 0 のみ受理します。
raw CSV の `far_field_v_m` と `far_charge_imbalance` は、最後の averaging window における絶対値の最大です。
符号と profile は対応する NPZ を確認します。

`ambient_electron.number_density_m3` と `ambient_ion.number_density_m3` は、遠方境界から入る drifting
Maxwellian 全体の振幅です。境界で実際に与えるのは $v<0$ の半空間だけなので、この値は計算された遠方密度と
同じとは限りません。v1 は electron source 振幅を `E(L)=0` へ自動調整しません。Zhao との overlap 比較では
Zhao root が求めた ambient electron source density と揃え、Zhao が解けない点では source density 自体も独立な
収束軸として走査します。固定値が不整合なままの far-boundary failure を、外部長不足と解釈してはいけません。

## 3. 認証済み subset を table にする

converter は raw CSV の hash、全 query の分類、重複、完全 Cartesian product を検査します。欠損点を補間せず、
既定では `steady` 以外を拒否します。

```bash
beach-kinetic-table \
  outputs/outer-kinetic-atlas/kinetic_response_raw.csv \
  outputs/matching_plane_response.csv \
  --range displacement_c_m2=-1e-11:1e-11
```

manifest は raw CSV と同じディレクトリから自動検出します。別の場所にある場合は `--manifest` を指定します。
`stationary_average` を採用する場合だけ `--allow-stationary-average` を明示します。

生成した CSV は既存の
[matching-plane table backend](MatchingPlaneCoupling.html#2-最小構成を作る)で読み込めます。

## 4. Production 候補にする前の比較

少なくとも次を個別に変えます。

- `nz` と各 species の `nv`
- velocity の上下限
- `z_length_m`
- `cfl`
- `max_time_s` と averaging window

Zhao の一意解がある点では、$\Phi_H$、ambient inward flux、PE return / escape flux を Zhao と比較します。
その後、同じ BEACH case を `zhao_online` と kinetic-derived table で実行し、表面電荷、gap 電位、吸収電流を
比較します。最後に matching-plane 高度を 3 点以上変え、主要量が収束しなければ `no_matching_overlap` と判断します。

この oracle は無衝突・非磁化・1Dで、PEを束と平均法線 energyから half-Maxwellianへ再構成します。強いPE beamで
遠方準中性が成立しない場合、外部長の延長だけで認証を通してはいけません。衝突、磁場、幾何拡散、energy-resolved
VDFが必要かを別 model で検証します。
