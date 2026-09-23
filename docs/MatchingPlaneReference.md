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
| online Zhao の複数根・root family 選択を調べる | [`zhao_root_selection`](#zhao_root_selection) |
| 応答表の exact header、11 列、単位、直積格子 | [Table backend の応答 CSV v1](#table-backend-の応答-csv-v1) |
| `beach-zhao-response` で応答表を作る | [Table backend用の応答表を作る](#table-backend用の応答表を作る) |
| 固定点の受理式と緩和式 | [固定点の数値契約](#固定点の数値契約) |
| grid、batch 幅、matching-plane 高度を収束させる | [収束と適用性を検証する](#収束と適用性を検証する) |

## Zhao Type B の電子密度

Type B の ambient 電子は、上流分布をエネルギー保存で輸送し、正電位へ加速される粒子の速度下限を残して積分します。
零ドリフトの場合は、BEACH の規格化で

$$
\hat n_{e,f}=\frac{\hat n_{e,\infty}}{2}
\operatorname{erfcx}\!\left(\sqrt{\hat\phi/\tau}\right)
$$

です。$\hat\phi=e(\Phi-\Phi_\infty)/(k_BT_{ph})$、$\tau=T_e/T_{ph}$、
密度は基準 PE 密度で規格化します。有限ドリフトでは、この式の erfc 引数を単にずらさず、
[上流 VDF の軌道積分](#ambient-vdf-とシース密度の整合性)を使います。
密度修正前に作った応答表は `beach-zhao-response` で再生成してください。table reader は既存 CSV の数値をそのまま使います。

## `photoelectron_closure`

online Zhao の `moment_matched_half_maxwellian`（既定）は従来どおり H の PE 外向き束と平均法線エネルギーを使います。
`energy_spectrum` は H の法線エネルギー別外向き流束 $F_H(K)=d\Gamma_H/dK$ を使います。
後者は `model="matching_plane_quasistatic"`、`response_backend="zhao_online"`、有効な
`photoelectron_species` が必要です。5 入力の table、stationary Zhao、PE なしには指定できません。
ambient electron / ion の密度式と A/B/C の branch 制約は変えません。

### 分布と電流

H の外向き通過を外部障壁との比較前に記録し、macro weight を面積と batch 時間で割って束へ直します。
再通過も数えます。H より内側で表面に戻った PE はこの束に含めず、外部から戻った PE は既存の
H での反射を経て内側を再追跡します。設定の表面放出電流を H の束に足し戻しません。

$K$ を eV、$B=\max(0,\Phi_H-\Phi_{pe,barrier})$ を同じ単位の障壁エネルギーとすると、

$$
\Gamma_{escape}=\int_B^\infty F_H(K)\,dK,\qquad
\Gamma_{return,H}=\Gamma_H-\Gamma_{escape}.
$$

外部の密度にも同じ分布を使い、到達する各エネルギーの束を局所速度で割ります。
Type A の極小までの区間と Type B では戻り population も加え、極小から遠方の区間では
透過 population だけを使います。bin ごとの密度とその電位積分は区分一定の $F_H$ に対して
解析的に評価します。Sagdeev 積分・遠方中性・profile の E² を満たす根を既存の有限探索で求めます。
根の未検出は一般的な不存在の証明ではありません。

`implicit_zero_mode=true` の PE escape target も上式を使い、平均エネルギーの指数式と混在させません。
表面の emission target は設定電流、全 return target は表面放出束からこの escape を引いた量です。
外部での return と表面まで戻った総量は別の量として扱います。

### 格子・反復・再開

`photoelectron_spectrum_bins_per_decade=N`（既定 32、正の int32）に対し、bin 境界は

$$
K_j=T_{pe,config}\left(10^{j/N}-1\right)\quad[j=0,1,\ldots]
$$

です。$T_{pe,config}$ は設定 PE 温度 [eV] で、測定エネルギーを覆うまで格子を延ばします。
各 bin の積分束 $f_j$ を保存し、その区間の $F_H=f_j/(K_{j+1}-K_j)$ を一定とします。
設定温度は格子の尺度であり、測定分布をその温度の Maxwell 分布へ戻す操作ではありません。
spectrum mode の平均エネルギーはこの bin 表現から計算するため、既存の標本平均とは離散化誤差だけ異なります。

初回、または spectrum を持たない旧 checkpoint からの再開時は、利用できる PE モーメントから
Maxwell 初期分布を作ります。各 trial では分布全体を `coupling_relaxation` で緩和し、
既存のモーメント収束条件に加えて $\sum_j|f_{j,observed}-f_{j,input}|$ が PE 束の許容値
$\max(\mathtt{coupling\_rtol}\,s_\Gamma,\mathtt{coupling\_atol}[1])$ 以下であることを要求します。
ここで `[1]` は第 1 成分を指します。有限な未収束 trial の warning と受理方針は従来どおりです。

観測分布と応答へ渡した分布は `matching_plane_spectrum_history.csv` に別々に保存します。
再開用には現在の checkpoint の summary に両分布・格子・応答入力を追加し、既存 scalar history の列は変えません。
保存済み格子と設定の PE 温度・bin 分解能が同じ場合は、非 Maxwell の形を含めて保存分布を復元します。
異なる場合は警告し、保存済みの束・平均エネルギーを持つ Maxwell 分布で初期値を作り直して H で再計測します。
これは再開用の初期値であり、以後の反復は測定分布を使います。既定の moment mode でも診断格子の変更だけでは
停止しません。出力列と units は[出力形式](OutputReference.html#pe-spectrum-の観測と応答入力)を参照してください。

分解能、ray 数、batch 幅、H の位置を変えて収束・収支を確認してください。これは外部の平面・無衝突・非磁化 1D
近似です。法線エネルギー以外の相関、外部の遅延 return、BEACH 内部の PE 体積空間電荷は追加していません。

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
符号変化がなければ外挿せず停止します。

online の Type B は、`zhao_branch="b"`、または `auto` で正の内向き電子 drift により Type B だけが物理的に
許される場合、境界電位 $\Phi_H$ を未知数にして陰的終点を直接解きます。各電位で上流中性条件から ambient 密度を
求め、Sagdeev 積分から $D_H(\Phi_H)$ を計算します。`energy_spectrum` では分布の bin 境界と各 bin の四分点も
探索するため、$D_H$ 方向に狭い解区間を粗い変位刻みで飛び越える問題を避けます。
有限探索は全数学根の列挙を保証しません。`continuation` に有効な seed があれば一意な最近傍根を選び、
初期 seed がない場合は一意根を要求します。複数の陰的終点を初期値の検出順で選びません。

その他の online 条件では $D_H$ を探索し、guard 付き secant と中点 fallback を使います。
前の outer 反復の終点（最初は $D_H^n$）を seed とし、valid な seed では明示更新の変位を自然なシース尺度

$$
D_{ref}=\sqrt{\epsilon_0 n_i e T_e}
$$

以下に抑えた初期幅から、幅を 2 倍ずつ最大 64 回拡張します。seed が明示 Type A / B / C の解領域外なら、
branch と整合する符号を $D_{ref}/32$ 刻み、最大 $8D_{ref}$ まで走査します。数値未保証の gap をまたいだ2点は
bracket とみなしません。これは table を永続的に拡張する処理ではなく、その batch の終点を探す処理です。
`auto` + `continuation` は初期点が未解決なら両符号を走査し、見つかった有効根をその区間内の探索 seed に使います。
この局所 seed は未受理の間は次 batch の accepted state に保存しません。
Zhao branch が終点より前に終わる、走査範囲または数値範囲を超える、または符号変化を見つけられない場合は
停止します。

有限な両端で符号変化を確認した後、丸め誤差で残差許容値に届かなかった場合は、残差が小さい方の端点を
採用して警告します。根を bracket できない場合、応答評価が非有限な場合、物理解がない場合は従来どおり停止します。

強い PE の A/B 共存は implicit 化だけでは解消せず、branch 別の可解性評価が必要です。

### `zhao_root_selection`

| 値 | 適用範囲 | 規則 |
|---|---|---|
| `require_unique` | online Zhao | query ごとに一意な物理解を要求する。`auto` が一意性を確認できなければ branch を代わりに選ばず停止 |
| `minimum_energy` | online Zhao | multistart で検出した候補から全シース電位エネルギーが最小の根を選択 |
| `continuation` | online Zhao、`implicit_zero_mode=true` | 最後に受理した endpoint を seed に局所追跡し、fallback では設定 branch 内で seed に最も近い検出根だけを受理 |

既定値は `require_unique` です。`continuation` は履歴依存の opt-in で、`auto` / `a` / `b` / `c` に使えます。
明示更新、table、stationary Zhao には使えません。

各 branch の multistart は最大 8 個の初期値を使います。電位候補は PE 流束分布の中央値、PE 平均エネルギー、
$\epsilon_0 E_H^2/(en_i)$、電子温度、冷イオンの運動エネルギー上限から定め、ambient 密度の初期値は
各候補の上流中性条件から求めます。固定の V や m$^{-3}$ に依存した初期値は使いません。

`zhao_root_selection="minimum_energy"` では、multistart で検出した各候補について表面から無限遠までの profile から

$$
U=-\frac{\epsilon_0}{2}\int_0^\infty E^2\,dx
$$

を評価し、最小の $U$ を選びます。明示 branch ではその branch 内の根、`auto` では数値的に検証できた A / B / C
候補を比較します。候補 branch の数値失敗で集合を確定できない場合と、最小値が相対 $10^{-6}$ 以内で縮退する場合は
停止します。この比較は [Mishra et al. (2023)](https://academic.oup.com/mnras/article/520/1/233/6987684) の
sheath potential-energy 比較に基づく候補選択規則であり、有限 multistart による全根の列挙や時間依存安定性を保証しません。
最小エネルギー根の切替点では応答が不連続になり得ます。backward-Euler 残差がその不連続をまたいでも通常の零点が
なければ、2根を混合せず数値失敗として停止します。

`continuation` は新規 run の初回に有限 multistart を行い、`require_unique` と同じ一意性条件を要求します。
エネルギー順位で初期根を選びません。batch trial は最後の受理済み endpoint から始め、trial 内で有効な
endpoint が得られた後は、その根を次の feedback 反復の Newton seed にします。初回 batch も同じ規則です。
Newton が収束しない、root を復号できない、profile 検査に通らない、または候補が大きく離れた場合だけ
full multistart へ戻ります。局所 Newton の受理判定は、PE 温度尺度で規格化した境界電位差、経路最低電位差、
ambient 密度比の対数の最大値

$$
d=\max\left(\frac{|\Delta\Phi_H|}{T_{pe}},\frac{|\Delta\Phi_{min}|}{T_{pe}},
\left|\log\frac{n_{e,\infty}}{n_{e,\infty}^{seed}}\right|\right)
$$

が 0.25 以下であることを要求します。経路最低電位は Type A で $\phi_m$、Type B で 0、Type C で $\Phi_H$ です。
full multistart でも同じ距離で最も近い根を選びます。最近傍距離を $d_1$、2 番目を $d_2$ としたとき、
$|d_2-d_1|\le10^{-6}\max(1,d_1)$ なら数値的に区別できないため、初期 guess の順では選ばず曖昧状態として停止します。
それ以外は、最近傍 root を距離 0.25 の内外にかかわらず受理します。0.25 は局所 Newton をそのまま使う fast path の
上限であり、full multistart 後の root family に対する物理的な距離上限ではありません。full multistart で root を
検出できない場合、または探索や profile 検査が数値的に失敗した場合は停止します。最後の有効 root の直後で解なしまたは
数値失敗となった implicit probe には二分を試し、branch 終端付近の root を粗い走査で飛び越えないようにします。
明示した branch は切り替えません。`auto` では検証済み A / B / C 候補が対象です。

この方法は pseudo-arclength continuation ではありません。full multistart と branch-boundary subdivision は局所 Newton が
失った root を再取得する手段ですが、同じ物理 family の保持、fold の位置、fold の通過可能性は証明しません。
初回と fallback の multistart にも、
有限個の初期値では全数学根を列挙できないという制限があります。

implicit root 探索で棄却した probe、固定点の未受理 trial、adaptive batch の棄却 trial は候補 root を accepted
continuation state に commit しません。trial 内の継続 seed は trial の棄却時に破棄し、accepted endpoint だけを
次の batch の seed にします。restart では保存済み
accepted response から seed の再構成を試み、再構成できない場合は初回の一意根探索へ戻ります。保存状態と
receipt は[出力形式リファレンス](OutputReference.html#matching_plane_quasistatic)を参照してください。

PE ありで moment closure または table を使う場合は half-Maxwellian 近似から

$$
\Gamma_{pe}^{escape}(D)=\Gamma_{pe}^{out}
\exp\left[-\frac{\max(0,\Phi_H(D)-\Phi_{pe,barrier}(D))}
{\langle K_{pe,n}^{out}\rangle}\right]
$$

を求めます。`energy_spectrum` では上の指数式を使わず、測定分布の障壁以上を積分します。
どちらも $q_{pe}<0$ として

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

## Ambient VDF とシース密度の整合性

ambient electron の密度は、流入束と同じ上流 drifting Maxwellian をエネルギー保存で写像します。
電位を V、$T_e$ を eV の数値で表し、$u=v_d/\sqrt{2eT_e/m_e}$ を内向き drift とすると、
Type A の電位極小を通過する population は

$$
\frac{n_e(\phi)}{n_{e,\infty}}=
\frac{1}{\sqrt\pi}\int_{\sqrt{-\phi_m/T_e}}^\infty
\frac{s\,e^{-(s-u)^2}}{\sqrt{s^2+\phi/T_e}}\,ds
$$

です。$s$ は上流の内向き法線速度を $\sqrt{2eT_e/m_e}$ で割った値で、$\phi_m<0$、$\phi\ge\phi_m$ とします。
Type B では下限を 0 にし、反射する区間では同じ上流分布の反射 population を加えます。
有限 drift を Boltzmann 因子とずらした erfc の積に置き換えません。PE 密度も放出点からの軌道と
反射・透過の速度範囲を保ちます。Type A は $\phi_m<\min(0,\Phi_H)$ を要求し、$\Phi_H<0$ も許容します。

中性条件・Sagdeev 積分条件の代数根だけでは物理解にならないため、全 profile の $E^2\ge0$ を検査します。
Type B では無限遠近傍の電荷密度の $\sqrt{\phi}$ 項も検査し、有限個の profile 点では見落とす負の $E^2$ を
棄却します。この係数には障壁位置の PE 流束密度を使い、bin 境界では両側の値を区別します。
低速 ambient 電子が完全反射する Type A / C に正の内向き drift を与えると、厳密な中性・零電場の無限遠条件の下では
無限遠近傍でこの条件を満たしません。有限上流境界は別の境界値問題であり、現在の online closure には含めません。

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
表には accepted endpoint の履歴がないため、`zhao_root_selection="continuation"` を指定した生成は拒否されます。
表生成には `require_unique` または `minimum_energy` を使います。
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
