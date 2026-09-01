title: matching-plane 数値・応答表リファレンス

Lang: [日本語](MatchingPlaneReference.md) | [English](MatchingPlaneReference.en.md)

# matching-plane 数値・応答表リファレンス

`surface_current_model.model="matching_plane_quasistatic"` の応答 CSV、面平均電荷の陰的更新、固定点収束条件を
調べるためのリファレンスです。model の選択、最初の 4 batch 実行、出力診断は先に
[matching-plane で外部シースを接続する](MatchingPlaneCoupling.html)を参照してください。

## 目的別索引

| 調べるもの | 節 |
|---|---|
| 秒スケールの batch 幅で面平均電荷だけを陰的に更新する | [`implicit_zero_mode`](#implicit_zero_mode) |
| 応答表の exact header、11 列、単位、直積格子 | [Table backend の応答 CSV v1](#table-backend-の応答-csv-v1) |
| `beach-zhao-response` で応答表を作る | [Table backend用の応答表を作る](#table-backend用の応答表を作る) |
| 固定点の受理式と緩和式 | [固定点の数値契約](#固定点の数値契約) |
| grid、batch 幅、matching-plane 高度を収束させる | [収束と適用性を検証する](#収束と適用性を検証する) |

## `implicit_zero_mode`

秒スケールの `batch_duration` で面平均電流の陽的更新が硬い場合、`implicit_zero_mode=true` で
面平均 $D_H$ だけを後退 Euler 更新できます。table と online Zhao の両方を選べます。

### 設定契約

| backend | $D_H$ の探索範囲 | feedback |
|---|---|---|
| `table` | CSV の $D_H$ 軸内。2 node 以上必要 | PE ありは正の PE flux / energy singleton、PE なしは両方 0。ambient outward は 0 の singleton |
| `zhao_online` | 選択 branch 内で現在値から探索 | PE moment は粒子固定点の各反復で更新。ambient outward は transparent |

どちらも `periodic2.lower_boundary_model="e_bottom_zero"` が必要です。table は監査済みの有限範囲、online は
CSV を用意せず組み込み Zhao を直接解く経路です。

```toml
[periodic2]
lower_boundary_model = "e_bottom_zero"

[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
implicit_zero_mode = true
```

この online 実行には `response_table_path` も `matching_query.csv` も要りません。query CSV は
`beach-zhao-response` で固定 table snapshot を別途作る場合だけ使います。

### 後退 Euler の終点

BEACH は

$$
D_H^{n+1}=D_H^n+hJ(D_H^{n+1})
$$

の符号変化を bracket して終点を解きます。table は CSV の両端を固定 bracket として二分法を使い、
符号変化がなければ外挿せず停止します。online は guard 付き secant と中点 fallback を使います。

online は前の outer 反復の終点（最初は $D_H^n$）を seed とし、valid な seed では明示更新の変位を自然なシース尺度

$$
D_{ref}=\sqrt{\epsilon_0 n_i e T_e}
$$

以下に抑えた初期幅から、幅を 2 倍ずつ最大 64 回拡張します。seed が明示 Type A / B / C の解領域外なら、
branch と整合する符号を $D_{ref}/32$ 刻み、最大 $8D_{ref}$ まで走査します。数値未保証の gap をまたいだ2点は
bracket とみなしません。これは table を永続的に拡張する処理ではなく、その batch の終点を探す処理です。
Zhao branch が終点より前に終わる、走査範囲または数値範囲を超える、または符号変化を見つけられない場合は
停止します。

signed scan は明示した `a` / `b` / `c` にだけ使います。既定の
`zhao_root_selection="require_unique"` では、`auto` が seed で一意な物理解を保証できない場合、implicit solver は
branch を代わりに選ばず停止します。したがって強い PE の A/B 共存は implicit 化だけでは解消せず、branch 別の
可解性評価が必要です。

`zhao_root_selection="minimum_energy"` では、multistart で検出した各候補について表面から無限遠までの profile から

$$
U=-\frac{\epsilon_0}{2}\int_0^\infty E^2\,dx
$$

を評価し、最小の $U$ を選びます。明示 branch ではその branch 内の根、`auto` では数値的に検証できた A / B / C
候補を比較します。候補 branch の数値失敗で集合を確定できない場合と、最小値が相対 $10^{-6}$ 以内で縮退する場合は
停止します。このエネルギー比較は候補選択規則であり、有限 multistart による全根の列挙や時間依存安定性を保証しません。
最小エネルギー根の切替点では応答が不連続になり得ます。backward-Euler 残差がその不連続をまたいでも通常の零点が
なければ、2根を混合せず数値失敗として停止します。

PE ありでは half-Maxwellian 近似から

$$
\Gamma_{pe}^{escape}(D)=\Gamma_{pe}^{out}
\exp\left[-\frac{\Phi_H(D)-\Phi_{pe,barrier}(D)}
{\langle K_{pe,n}^{out}\rangle}\right]
$$

を求め、$q_{pe}<0$ として

$$
D_H^{n+1}=D_H^n+h\left[
q_e\Gamma_e^{in}(D_H^{n+1})+q_i\Gamma_i^{in}(D_H^{n+1})
-q_{pe}\Gamma_{pe}^{escape}(D_H^{n+1})\right]
$$

の終点を解きます。PE なしでは最後の PE 項を除き、

$$
J=q_e\Gamma_e^{in}+q_i\Gamma_i^{in}
$$

だけを使い、PE target は作りません。

table implicit の PE moment は CSV の singleton 値に固定されます。online implicit は、現在の PE feedback
$X^m$ でこの終点を解き、同じ trial の粒子追跡で得た PE moment を緩和し、次の反復で終点を解き直します。
したがって PE return を含む outer feedback と $D_H^{n+1}$ は入れ子に整合されます。

陰的になるのは $k=0$ の面平均だけです。要素別 $k\ne0$ 分布は batch 開始場から追跡します。したがって、
6 s のような幅を使えるかは、局所電位変化、粒子 sampling、root bracket、物理範囲を別々に検証します。
実務上の比較手順は [`batch_duration` の選択](BatchDurationStability.html)を参照してください。

## Table backend の応答 CSV v1

header より前に整合面高度を 1 回だけ書きます。この値は `domain.box_max` の z 成分と一致させます。

```csv
# matching_plane_z_m=1.0e-3
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,electron_outward_number_flux_m2_s,ion_outward_number_flux_m2_s,matching_potential_v,electron_inward_number_flux_m2_s,ion_inward_number_flux_m2_s,electron_access_potential_v,ion_access_potential_v,photoelectron_barrier_potential_v
```

最初の 5 列が入力軸、後ろの 6 列が response です。

| 列 | 単位 | 意味 |
|---|---|---|
| `displacement_c_m2` | C/m2 | 整合面直下の平均 $D_z$。+z が正 |
| `photoelectron_outward_number_flux_m2_s` | 1/(m2 s) | 整合面へ到達した外向き PE 束 |
| `photoelectron_outward_mean_normal_energy_ev` | eV | 外向き PE の平均法線運動 energy |
| `electron_outward_number_flux_m2_s` | 1/(m2 s) | ambient electron の外向き束 |
| `ion_outward_number_flux_m2_s` | 1/(m2 s) | ion の外向き束 |
| `matching_potential_v` | V | 外部シースが返す $\Phi_H$ |
| `electron_inward_number_flux_m2_s` | 1/(m2 s) | BEACH へ入る electron の総束。外部 return を含められる |
| `ion_inward_number_flux_m2_s` | 1/(m2 s) | BEACH へ入る ion の総束。外部 return を含められる |
| `electron_access_potential_v` | V | electron reservoir から整合面への access bottleneck |
| `ion_access_potential_v` | V | ion reservoir から整合面への access bottleneck |
| `photoelectron_barrier_potential_v` | V | PE が外向きに越える最大外部 barrier |

### 表の格子と値の契約

- 5 入力軸の完全な Cartesian product を持ち、重複点と欠損点がない。行順は任意。
- flux、PE 平均 energy、出力 flux は非負で、すべての値が有限。
- 2 node 以上の feedback 軸は初期評価のため 0 を含む。BEACH は範囲外へ外挿しない。
- 外部 model が依存しない feedback 軸は singleton にする。singleton は任意の有限 query を受理し、その依存性を無効化する。
- 4 potential 列は上流 0 V の同じ gauge を使う。
- 数値 token は十進実数だけとし、`/`、`2*0`、空欄などの Fortran list-directed 制御記法を使わない。

補間は最大 32 corner の多重線形補間です。読込メモリは行数に比例し、MPI では各 rank が表を保持します。

### Table backend用の応答表を作る

production 表は、同じ $H$、上流分布、符号規約を使う独立した Zhao / 1D PIC sweep から作ります。入力する PE moment は
壁面放出量ではなく、整合面を実際に通過した束と法線 energy です。非単調電位では、外部 profile 全体の最大 barrier を使います。
表生成 code、上流条件、solver version、単位変換も production data と一緒に保存してください。

組み込み online Zhao を事前評価して table 形式の snapshot を作る場合は、次を実行します。

```console
beach-zhao-response \
  examples/periodic2_matching_plane_zhao_online.toml \
  examples/matching_plane_zhao_query_grid.csv \
  response.csv
```

設定ファイルは完全な `response_backend="zhao_online"` matching case とし、`response_table_path` は指定しません。
query CSV は空行と `#` comment を許し、最初の非 comment 行を次の exact header にします。全値は有限、flux と
PE energy は非負です。

```csv
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,electron_outward_number_flux_m2_s,ion_outward_number_flux_m2_s
```

v1 generator では PE flux 軸に 0 を含め、PE energy は singleton にします。正の PE flux node がある場合は
energy も正にします。transparent な ambient outward 2 軸は 0 の singleton です。

5 軸の完全な直積を与え、全 query が解けた場合だけ 11 列の `response.csv` を書きます。sample grid は固定 3 eV の
配線確認用で、production 範囲ではありません。PE energy 依存を持つ production 表は、独立した外部 solver から
直接生成してください。

`beach-zhao-response` が見つからない場合は、[この文書に対応する現行版](Installation.html#このドキュメントと一致する版をインストール)
をインストールしてください。生成表を使う run は別の設定ファイルで `response_backend="table"` と
`response_table_path="response.csv"` を指定します。online backend を直接使う run には応答表は不要です。

### Zhao の解領域を調べる

`beach-zhao-atlas` は、指定した matching-plane moment に対して Zhao A/B/C を独立に評価します。
simulation 用の branch を選ぶ前に、多重解、物理解なし、solver の数値失敗を分けて調べるための offline 診断です。
応答表は作らず、BEACH runtime の設定も変更しません。

```console
beach-zhao-atlas \
  examples/periodic2_matching_plane_zhao_online.toml \
  query_grid.csv \
  atlas.csv
```

設定ファイルは `response_backend="zhao_online"` の完全な matching case とします。`zhao_branch` と
`zhao_root_selection` の設定値にかかわらず、atlas は A/B/C をそれぞれ `require_unique` で評価します。
query CSV は空行と `#` comment を許し、最初の非 comment 行を次の exact header にします。完全な直積は不要で、
調べたい点だけを並べられます。

```csv
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev
```

`atlas.csv` は query ごとに A/B/C の 3 行を持ちます。`status` は次の意味です。

| `status` | 判定 |
|---|---|
| `ok` | その branch の一意な物理解を solver が認証した |
| `no_physical_solution` | 現行 solver がその branch を物理的に不適格と判定した |
| `numerical_failure` | 存在・不存在を数値的に判定できなかった |
| `ambiguous_within_branch` | 同じ branch 内で複数根が残った |
| `invalid_input` | flux や PE energy の入力契約に違反した |

複数の branch が `ok` なら `multiple` です。1 branch だけが `ok` でも、ほかに
`numerical_failure` または `ambiguous_within_branch` があれば一意性は未認証です。全 branch が
`no_physical_solution` の場合だけ `no_root` とし、数値失敗を `no_root` に含めないでください。
この判定は現行の有限個の初期値と profile 検査による solver certificate であり、数学的な不存在証明ではありません。

## 固定点の数値契約

feedback vector の成分順は

$$
X=(\Gamma_{pe}^{out},\langle K_{z,pe}\rangle^{out},\Gamma_e^{out},\Gamma_i^{out})
$$

です。active 成分 $j$ の backend scale を $s_j$、相対許容値を $r$、絶対許容値を $a_j$ とすると、

$$
|X_{raw,j}^{m+1}-X_j^m|\le\max(r s_j,a_j)
$$

を全成分で満たした trial を受理します。未収束時は `coupling_relaxation` を $\alpha$ として

$$
X^{m+1}=(1-\alpha)X^m+\alpha X_{raw}^{m+1}
$$

と更新します。inactive 成分は判定から除外し、その `coupling_atol` は 0 にします。

backend scale と inactive 成分は次のように決まります。

| backend | $s_j$ | inactive 成分 |
|---|---|---|
| `table` | 対応する active feedback 軸の最大値と最小値の差 | singleton feedback 軸 |
| `zhao_online` | Zhao model が定める基準 flux または基準 energy | transparent な ambient electron / ion outward 軸 |

$\Delta_j=X_{raw,j}^{m+1}-X_j^m$ とすると、出力 `matching_plane_residual` は

$$
\max_j \rho_j,\qquad
\rho_j=
\begin{cases}
r|\Delta_j|/a_j, & a_j>r s_j,\\
|\Delta_j|/s_j, & a_j\le r s_j.
\end{cases}
$$

です。この正規化により、絶対許容値が支配する成分があっても、収束した trial では
`matching_plane_residual <= coupling_rtol` になります。

history の応答列は accepted trial の $X^m$ で評価した値、feedback 列は同じ trial で観測した
$X_{raw}^{m+1}$ です。収束後に、実行していない緩和更新 $X^{m+1}$ を加えて記録することはありません。

`coupling_max_iterations` までに収束しなくても、feedback と応答が有限なら最終 trial を warning 付きで commit します。
`matching_plane_residual > coupling_rtol` と最大反復回数が、その batch の未収束 receipt になります。次 batch は観測した
feedback から再開します。table の active 軸が範囲外の場合、online solve が失敗した場合、または非有限値が出た場合は、
有効な trial がないため停止します。state と残差の出力契約は
[出力形式リファレンス](OutputReference.html#matching_plane_quasistatic)を参照してください。

## 収束と適用性を検証する

1. `coupling_rtol`、`coupling_atol`、緩和係数、粒子数を変えて accepted observables を比較する。
2. table の grid 解像度と範囲、または online Zhao の明示 branch を独立に変える。
3. `batch_duration`、mesh、periodic cell を収束させる。
4. $H$ を外部 model との overlap region 内で動かし、grain charge、gap potential、PE escape fraction の不変性を調べる。

$H$ 依存性が小さいことは、この連成に固有の中心的な検証です。実行完了、数値収束、物理的妥当性は別々に判定してください。
