# BEACH 仕様書（現行 Fortran 実装）

## 1. 目的

BEACH は、三角形境界要素上の電荷蓄積とテスト粒子追跡を行うシミュレータです。

- 境界は三角形メッシュ
- 粒子運動は Boris push
- 境界衝突時は吸収して要素へ電荷を堆積
- `batch_count` 個の accepted batch を処理し、統計と履歴を更新

計算の主系は Fortran（`src/`, `app/`）で、Python（`beach/`）は設定補助、後処理、可視化を担います。

## 2. スコープ

### 2.1 実装済み

- template / OBJ の三角形メッシュ
- 三角形上の一定面密度を使う `triangle_p0` 静電場
- direct / treecode / FMM field solver
- free-space と 2 軸周期 `periodic2`
- 一様外部電場 `e0` と磁場 `b0`
- Boris 法による粒子更新
- 線分と三角形の最初の交差判定
- insulator accumulation
- free-space に限った mesh 単位の浮遊 conductor 再配分
- open / reflect / redistributed_reflect / periodic の box 境界
- `volume_seed` / `plane_source` / `photo_raycast` と deprecated `reservoir_face`
- species 別の simulation boundary reservoir 流入
- 境界 reservoir の流入補正
- closed photoelectron の局所反射と neutral-return closure
- 応答表または組み込み `zhao_online` による matching-plane 準定常外部シース連成
- checkpoint 再開

### 2.2 未実装・予約

- BEACH 内で解く full-VDF / 1D PIC / time-dependent outer-sheath solver と外部領域の粒子 transport
- 表面導電・拡散モデル
- periodic2 と conductor の併用
- 誘電体分極・誘電境界条件
- 反射や二次電子放出を含む一般化衝突モデル

## 3. データモデル

### 3.1 メッシュ

`mesh_type` は三角形頂点、重心、法線、面積、要素電荷 `q_elem`、mesh ID、および衝突検索用
AABB/grid cache を保持します。OBJ の scale、rotation、offset はこの順に適用します。

要素電荷は総電荷 [C] であり、場計算時に面積で割って一定面密度へ変換します。

### 3.2 粒子

`particles_soa` は位置、速度、電荷、質量、macro weight、生存フラグ、species ID、
`photo_raycast` の放出元要素を保持します。

### 3.3 統計

`sim_stats` は processed / absorbed / escaped / max-step survivor、accepted batch 数、累積物理時間、
適応 batch の棄却数と直近診断、`last_rel_change` を保持します。

`tol_rel` は監視・出力値であり、早期停止条件ではありません。

## 4. 1 batch の計算手順

1. 現在の要素電荷から静電 snapshot を構築する
2. boundary reservoir inflow、`plane_source`、`reservoir_face`、`volume_seed`、`photo_raycast` の粒子を生成する
3. 各粒子を最大 `max_step` 回まで進める
4. mesh hit、box escape、max-step survivor を分類する
5. 吸収電荷と放出元の逆符号電荷を thread-local buffer へ集計する
6. closed PE が有効なら、解決済み帰還先へ neutral-return 補正を適用する
7. MPI reduction 後の候補電荷を検証して commit する
8. accepted batch の統計、履歴、charge ledger、累積時間を更新する

`periodic2.max_nonzero_mode_potential_step > 0` の場合は、設定された `batch_duration` を最大幅として
`h, h/2, h/4, ...` を試します。候補電荷による cached $k\ne0$ 電位変化が上限以下の trial だけを
commit します。棄却時は RNG と macro 粒子数残差を trial 前へ戻し、統計・履歴・ledger を更新しません。

## 5. 場モデル

### 5.1 free-space

`field_solver="direct" | "treecode" | "fmm" | "auto"` を選べます。

- direct: 厳密な free-space panel 和
- treecode: near は厳密 panel、far は monopole
- FMM: near は厳密 panel、far は三角形 moment を使う FMM
- auto: 要素数に応じて direct または FMM

面上電位は連続で、法線電場は `sigma/epsilon0` だけ跳びます。重心電位と principal-value 電場を
自己項として使います。`surface_side` は各 mesh/template で明示します。

### 5.2 periodic2

`[field_boundary] mode="periodic2"` は `[domain] periodic_axes=["x", "y"]` を要求します。
mesh は primitive cell の base element だけを保持し、場と衝突判定が周期 image を扱います。

`field_periodic_far_correction`:

- `none`: 指定した有限 image shell のみ
- `auto`: 互換用に受理し、`none` に正規化
- `cached_kneq0`: 無限周期の非零 Fourier mode を versioned operator cache から評価

型付き `[periodic2]` 設定では、小規模な検証計算に限り
`nonzero_mode_backend="panel_spectral_reference"` と direct solver の組合せも利用できます。
これは非零 Fourier mode の収束参照であり、外部プラズマモデルには依存しません。

`field_periodic_far_correction="none"` は有限画像近似であり、無限周期解ではありません。
`m2l_root_oracle` は削除済みです。production の無限周期計算には `cached_kneq0` を使用します。
cache fingerprint は generator version を含み、物理的 zero mode の評価は高さ breakpoint の
二分探索により 1 点あたり O(log n) です。

`cached_kneq0` の物理的 $k=0$ は `periodic2.zero_mode_policy="exclude_k0"` と
`lower_boundary_model="symmetric_vacuum" | "e_bottom_zero"` で別に構築し、場へ一度だけ加えます。
このzero-mode構築自体は外部プラズマ応答を合成しません。matching-plane modelを選んだ場合、その応答は
別の連成層で適用します。

### 5.3 領域・粒子境界・reservoir

公開 TOML は topology、場、外向き粒子作用、外部 reservoir 条件、species 別流入を分離します。

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape" # または "potential_barrier"

[reservoir]
inflow_model = "source_vdf" # または "infinity_barrier"
phi_infty = 0.0
face_potential_grid_n = 5

[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0

[particles.species.boundary_inflow]
z_high = "reservoir"
```

- `boundary_inflow`: 非周期 box 面全体から外部 reservoir VDF を流入させる
- `source_vdf`: 指定した boundary VDF を補正せず境界から注入する
- `infinity_barrier`: `phi_infty` と流入面平均電位の scalar barrier で法線 VDF を補正する
- `escape`: open 面到達で粒子を消滅する
- `potential_barrier`: event 位置の局所電位と法線運動エネルギーから反射／escape を判定する

`boundary_inflow` は外向き粒子作用を上書きしません。周期面への reservoir 流入、および有効な外向き作用が
`open` でない流入面は拒否します。
同じ species では `source_mode="volume_seed"` とだけ併用でき、`plane_source`、`photo_raycast`、
`reservoir_face` との併用は fail closed です。

外部プラズマ profile、外部領域の particle transport、delayed return queue は現行スコープ外です。
削除済みの `[outer_plasma]`、`[coupling]`、`[external_boundary]` は unknown input として拒否します。

## 6. 粒子前進・衝突・box 境界

public な粒子状態は同一時刻の `x, v` です。空間電場は予測中点
`x_mid = x + 0.5*v*dt` で評価し、Boris の速度更新後に
`x_new = x + 0.5*(v + v_new)*dt` とします。

mesh 衝突は Möller–Trumbore 法で最小 `t` を選びます。AABB と、既定では一様 grid + 3D-DDA で
候補を絞ります。非有限値や進行不能は hit なしへ黙って変換せず、明示 status で停止します。

`domain.periodic_axes` は全粒子種に共通の topology です。非周期面の粒子 action は
global `[particle_boundary]` と species 別 `[particles.species.boundary]` で指定します。

- `open`: escape
- `reflect`: 法線速度を反転し、接線速度と接線位置を維持
- `redistributed_reflect`: `reflect` と同じ速度作用に加え、event 面内の位置を一様再配置
- `periodic`: 反対側へ wrap

単一面の `redistributed_reflect` は面内 2 軸を box span の両端 guard を除く範囲から再標本化します。
edge / corner の同時 event では、
event mask に含まれない軸だけを再標本化し、event 軸の座標は内側 guard に置きます。mesh hit と最初の
box event は chord parameter で順序付けます。reflect / redistributed_reflect / periodic 後の残り時間は
同じ Boris 規約で再積分し、1 local step あたり最大 8 box event を許します。

species 別の6面値は `inherit | open | reflect | redistributed_reflect` です。`inherit` はglobal粒子 actionを使います。
周期面はspeciesから上書きできません。

## 7. 粒子 source と境界流入

### 7.1 `volume_seed`

各 species で `npcls_per_step` 個を batch ごとに生成します。

`boundary_inflow` を持つ species では `npcls_per_step=0` を許容します。

### 7.2 `plane_source`

`pos_low` と `pos_high` で box 内部の axis-aligned 矩形面を定義し、`source_normal` 方向へ一方向 flux を
生成します。矩形面はちょうど 1 軸で zero thickness、残る 2 軸で正の面積を持ちます。法線座標は
box 境界より厳密に内側で、接線範囲は box 内に置きます。
`source_normal` は zero-thickness 軸に沿う非ゼロベクトルで、実装内部で単位ベクトルへ正規化します。

Maxwell reservoir または速度 grid の flux、面積、batch duration から macro 粒子数を決めます。
`reservoir.inflow_model`、`phi_infty`、`face_potential_grid_n` は内部平面へ適用しません。
任意の内部面では外部 plasma 側と上流基準電位の対応が一意でないため、流入側の potential/barrier 補正は
新しい設定では simulation boundary に結び付いた `boundary_inflow` が所有します。非推奨の
`reservoir_face` は既存 case の数値挙動を保つため、従来の補正を互換動作として維持します。

### 7.3 species 別 boundary reservoir 流入

`[particles.species.boundary_inflow]` の 6 面キーへ `reservoir` を指定し、非周期 box 面全体から流入させます。
複数面を指定でき、macro 粒子の端数は species と face の組ごとに batch 間で繰り越します。
MPI では global 個数を rank へ分配します。

流入 flux と面積、batch duration から macro 粒子数を決めます。`target_macro_particles_per_batch` 指定時は
macro weight を自動解決します。

`reservoir.inflow_model="infinity_barrier"` では各流入面の `N x N` 電位標本を使い、
`reservoir.phi_infty` とのエネルギー差から法線速度 cutoff と到達速度を同じ写像で補正します。

### 7.4 `reservoir_face`（deprecated）

既存 case の `inject_face` と `pos_low` / `pos_high` による box 面開口の挙動を維持します。
`boundary_inflow` または `plane_source` へ暗黙変換しません。外部 plasma 条件の新規 case は
`boundary_inflow`、内部矩形面は `plane_source` を使います。

### 7.5 `photo_raycast`

指定方向へ `rays_per_batch` 本を飛ばし、各 ray の最初の mesh hit からだけ放出します。
`deposit_opposite_charge_on_emit=true` では放出元へ逆符号電荷を加えます。

closed PE は次の組合せです。

- 負電荷 `photo_raycast`
- `deposit_opposite_charge_on_emit=true`
- species の `inject_face` に対する有効粒子 action が `reflect` または `redistributed_reflect`
- `surface_charge_closure="neutral_return"`

放出電荷を $S<0$、解決済み再吸収を $R<0$、max-step 未解決を $U<0$ とし、`S=R+U` を確認した上で、
帰還先への deposit を `S/R` 倍します。これにより PE species の表面総電荷増分を 0 にしつつ、
帰還先分布による表面内再分配を保持します。

未解決率 `U/S` が 0.05 を超える、`R=0`、実 escape、soft discard、terminal 分類不整合のいずれかでは
補正せず fail closed とします。

### 7.6 species 別 fixed-current closure

`surface_charge_closure="fixed_current"` は、species の吸収 channel と `photo_raycast` の放出反作用 channel を
独立に補正します。signed target charge は accepted trial の batch 幅 $\Delta t$ に対して

$$
Q_{s,\mathrm{abs}}^{\mathrm{target}}=I_{s,\mathrm{abs}}^{\mathrm{target}}\Delta t,
\qquad
Q_{s,\mathrm{emit}}^{\mathrm{target}}=I_{s,\mathrm{emit}}^{\mathrm{target}}\Delta t
$$

です。raw channel 電荷を $R$ とすると、各要素の raw deposit を同じ $Q^{\mathrm{target}}/R$ で倍率化し、空間分布を
保持します。`target_absorbed_current_a` の符号は `q_particle` と一致し、
`target_emission_current_a` の符号は `q_particle` と逆です。非ゼロ target に対して $R=0$、非有限倍率、負倍率なら
commit 前に fail closed とします。

この倍率化が保持するのは raw Monte Carlo 標本の経験分布であり、母分布の統計精度ではありません。raw hit が
1 件なら target 全量をその hit 要素へ割り当てます。必要標本数や許容倍率は mesh 解像度と評価量の許容誤差に依存するため、
BEACH は恣意的な count / scale 閾値では停止しません。`charge_ledger.csv` の raw charge、`absorbed_count` / `emitted_count`、
`fixed_*_weight_scale`を確認し、粒子数・ray 数・batch 幅・乱数 seed を変えた要素別電荷分布の収束で妥当性を判定します。
ledger の保存残差が小さいことは電荷収支を検証しますが、この統計収束を保証しません。

PE の emission と return は別 channel のまま扱い、net current を倍率分母に使いません。`fixed_current` と
`neutral_return` は同じ species で排他です。外部 return VDF の注入と top reflection / `neutral_return` の併用は
同じ return current を二重計上するため、構成上使用しません。

### 7.7 自動表面電流model

トップレベル`[surface_current_model]`は、設定型、model dispatch、model固有solver、species channelとkinetic境界写像への割当を分離します。
未指定または`model="none"`はtargetを生成せず、speciesに記述した手動targetを使います。初期実装の
`model="zhao_stationary"`は、Zhao A/B/Cの平面・無衝突・非磁化シースについて零電流定常根をrun開始時に一度解きます。
新しい電流modelは同じdispatch resultへspecies別の吸収・放出targetと診断値を返すことで追加します。

Zhao modelはambient electronとcold ionを参照し、PE有効時はphotoelectronも明示的に参照します。各speciesは
`surface_charge_closure="fixed_current"`を要求し、手動targetとの併用を禁止します。単価電荷、electron/PEの同一質量、
z-highからの内向きambient流入、負電荷`photo_raycast`の放出反作用、PEのopenなz-high境界、$T_i\le0.1T_e$を
fail-closedに検証します。
非磁化closureなので`sim.b0`はゼロを要求します。Zhao固有の0 V reservoirと速度写像を使うため、genericな
`reservoir.inflow_model="infinity_barrier"`との併用も拒否します。
電流密度から電流へ変換する面積は`reference_area_m2`、省略時はdomainのx-y面積です。

`photoelectron_source_scale=0.0`は光電子なしの契約です。この場合はPE固有設定とPE speciesを省略し、
`zhao_branch`は`"auto"`または`"c"`を指定します。Zhao Type Cの$J_e+J_i=0$を解いて
electron/ionの吸収targetとkinetic境界写像だけを生成します。
PE emission、return、escape targetは生成しません。

$n_{pe,0}=s_{UV}n_{pe,ref}\sin\alpha$とし、解いたPE emissionを$J_{emit}>0$、escapeを$J_{escape}>0$とすると、

$$
J_{return}=J_{escape}-J_{emit}\le0,
\qquad
J_e+J_i+J_{escape}=0
$$

です。$AJ_e$、$AJ_i$、$AJ_{return}$を各吸収channel、$AJ_{emit}$をPE放出反作用channelへ渡します。
さらに、外向きPE粒子電流$-AJ_{escape}$を外部escape targetとして記録します。escape targetは表面要素へ
depositせず、rawな境界escape電荷との差を外部境界closureの補正としてledgerへ記録します。したがって、適用後は

$$
AJ_{return}+AJ_{emit}-AJ_{escape}=0,
\qquad
AJ_e+AJ_i+AJ_{return}+AJ_{emit}=0
$$

を独立に検証できます。PEの大きなemission、return、escapeを別channelで補正し、差であるnet PE電流を
scaleの分母に使いません。raw軌道統計は上書きせず、raw / target / appliedを別々に出力します。
Zhao の零電流根と budget residual は総電流 closure を検証するだけで、raw kinetic map の空間的な統計精度を
保証しません。各 target channel には上記 fixed-current の収束確認を適用します。

Zhaoの定常電位は`zhao_barrier_v1` kinetic contractとしてz-highへも適用します。ambientの設定Maxwell VDFを
無限遠電位0 Vのreservoir分布とし、Type A electronは$\phi_m$、Type B/C electronとionは0 Vを外部access
bottleneckとして使います。各batchのz-high面平均電位を$\phi_f$とすると、流入tailは外部bottleneckと
$\phi_f$の両方へ到達できる粒子から生成し、法線速度を0 Vから$\phi_f$までエネルギー保存で写像します。

PE放出速度は従来どおり設定した表面half-Maxwellianです。PE、ambient electron、ionがopenなz-highを外向きに
横切るとき、粒子位置の局所電位から残りのZhao barrierを評価します。electron/PEのbarrier電位はType Aで
$\phi_m$、Type B/Cで0 V、ionは0 Vです。法線運動エネルギーが不足する粒子はz-highで鏡面反射してreturnへ、
十分な粒子だけをescapeへ分類します。固定電流targetはその後に各channel総量を正規化するため、kinetic写像は
rawなreturn/escape分類と空間分布を決め、Zhao電流は総電流収支を決めます。

このstationary modelはbox外の電場・空間電荷・Debye shielding・return軌道・遅延を解かず、run中の表面電位に応じて
targetを再計算しません。z-high反射は外部turning pointまでの距離・飛行時間を省略するadiabaticな境界closureです。
外部シースの過渡解ではなく、BEACHの軌道追跡から得る空間分布へ固定総電流を与えるclosureです。

### 7.8 matching-plane 準定常連成

`model="matching_plane_quasistatic"`は、`domain.box_max`のz成分をmatching plane $H$ とし、外部1Dシースを
非線形境界演算子としてBEACHへ接続します。`response_backend="table"`（既定）は事前計算済み応答表、
`response_backend="zhao_online"`は有限$H$のcharge-driven Zhao A/B/C準定常solveを使います。
全mesh頂点は$H$より厳密に下へ置きます。BEACHは
全triangleの実表面電荷を保持し、$H$より下の
$k=0$と$k\neq0$を解きます。外部Zhao/PIC場をmicrodomainへ重ねず、Zhaoの定常壁電位または表面電荷も追加しません。
speciesの設定は`surface_charge_closure="explicit"`のままです。既定の陽的更新はraw trajectory depositを使います。
後述の`implicit_zero_mode=true`だけは、粒子追跡が与えた要素別分布を保ちながらchannel総量を陰的終点へ正規化します。

このmodelはexplicitな`periodic2`構成だけに対応します。nonzero backendは`cached_kneq0`または
`panel_spectral_reference`、zero-mode policyは`exclude_k0`、lower boundaryは`e_bottom_zero`または
`symmetric_vacuum`とします。x/yはperiodic、z-low/z-highはopen、`sim.e0`と`sim.b0`はゼロです。
enabled speciesはambient electron、ion、および任意のphotoelectron roleだけを別speciesとして明示し、前2者はz-high
reservoir流入、PEを指定する場合は負電荷の`photo_raycast`かつopenなz-highを使います。generic `infinity_barrier`、手動fixed-current target、
`reference_area_m2`は併用しません。面積はdomainのx-y面積、$H$はbox上端、更新間隔は
1 accepted batchから導出し、重複parameterを公開しません。multiple-box-event policyは`abort`または
累積率で制限した `soft_discard` とします。soft discard の累積件数を $D$、accepted batch で
処理した累積 macro particle 数を $P$、累積絶対 macro charge を $Q$ とすると、commit 前の停止条件は

$$
\left(D>G\ \text{and}\ \frac{D}{P}>f_{\mathrm{limit}}\right)
$$

です。`multiple_box_events_soft_discard_count_grace` の既定値 $G=1000$ は累積件数の単独上限ではなく、
率判定を開始する件数猶予です。`multiple_box_events_soft_discard_fraction_limit` の既定値は $10^{-6}$、
制約は $0<f_{\mathrm{limit}}\le1$ で、$G\ge0$ とします。いずれの閾値も等値では停止しません。
`multiple_box_events_soft_discard_abs_charge_limit` は停止条件ではなく、累積絶対電荷の初回超過を知らせる
警告閾値です。
`summary.txt` と checkpoint には `multiple_box_events_soft_discarded`、
`multiple_box_events_soft_discarded_abs_charge_C`、および $D/P$ から導出した
`multiple_box_events_soft_discard_fraction` を残します。累積率は長い正常履歴によって後半の burst を
希釈しうるため、監査では batch ごとの集約 log も確認します。
`multiple_box_events_retry_backend="upper_panel_fourier"`は、通常の`cached_kneq0`場を変更せず、
`multiple_box_events`となった 1 step だけを元の状態から再試行します。再試行の$k\neq0$場は triangle P0 電荷の
有限 Fourier 展開を全 mesh 頂点より上で因子化し、既存の$k=0$場と外部一様場を合成します。potential-barrier
境界電位も同じ展開とgaugeで評価します。全評価点が成立域に
入らない場合または再試行も失敗する場合は、元の status に対して設定済み policy を適用します。

`response_backend="table"`は`response_table_path`を必須とし、`zhao_branch`を含むZhao固有keyを拒否します。
response CSV v1は、headerより前に一意な`# matching_plane_z_m=<finite>`を持ち、その値を$H$と照合します。
入力5列は

$$
(D_H,\ \Gamma_{pe}^{out},\ \langle K_{pe,n}^{out}\rangle,\
 \Gamma_e^{out},\ \Gamma_i^{out})
$$

であり、出力6列は

$$
(\Phi_H,\ \Gamma_e^{in},\ \Gamma_i^{in},\
 \Phi_{e,access},\ \Phi_{i,access},\ \Phi_{pe,barrier})
$$

です。列名と単位は`docs/MatchingPlaneReference.md`を正本とします。rowは5入力軸の完全Cartesian productで、
重複・欠損・非有限値を拒否します。flux、PE平均法線energy、出力fluxは非負です。2 node以上のfeedback軸2--5は
初期評価のためゼロを含みます。2 node以上の軸はclosed range内で最大32 cornerの多重線形補間を行い、外挿しません。
singleton軸はnodeが0以外でもよく、そのfeedback依存を意図的に無効化し、任意のfinite queryを係数0で受理します。
tableはprocess内でpathごとのimmutable snapshotとして読み、canonical axis/value列をmodel fingerprintへ含めます。
potential 4列は同じgaugeを使い、外部modelの上流reservoirを0 Vとします。inward VDFはこの0 Vから
access potentialと$\Phi_H$へ写像するため、potential列だけの定数shiftは同値ではありません。

`implicit_zero_mode=true`は、硬い面平均帯電だけを後退Eulerで更新します。table / `zhao_online`の両backendで
使用でき、共通して`lower_boundary_model="e_bottom_zero"`を要求します。tableは2 node以上の$D_H$軸と
singletonのfeedback軸2--5を要求します。singleton参照値はPEありでは
$\Gamma_{pe}^{out}>0$、$\langle K_{pe,n}^{out}\rangle>0$、PEなしではこの2値を両方0とし、どちらも
$\Gamma_e^{out}=\Gamma_i^{out}=0$とします。onlineはresponse/query CSVを要求せず、現在のfeedbackを使います。
PEありではhalf-Maxwellian近似から

$$
\Gamma_{pe}^{escape}(D)=\Gamma_{pe}^{out}
\exp\left[-\frac{\Phi_H(D)-\Phi_{pe,barrier}(D)}
{\langle K_{pe,n}^{out}\rangle}\right]
$$

を求め、

$$
D_H^{n+1}=D_H^n+h\left[q_e\Gamma_e^{in}(D_H^{n+1})
+q_i\Gamma_i^{in}(D_H^{n+1})-q_{pe}\Gamma_{pe}^{escape}(D_H^{n+1})\right]
$$

を解きます。PEなしでは最後のPE項を除きます。tableは$D_H$範囲の両端でbracketし、二分法で解きます。
両端が根を挟まない場合は外挿せず停止します。onlineは前のouter反復の終点（最初は$D_H^n$）をseedとし、
validなら明示終点変位を$D_{ref}=\sqrt{\epsilon_0n_i eT_e}$以下に抑えた初期幅から最大64回まで2倍にします。
seedが明示A/B/C branchの解領域外なら、branchと整合する符号を$D_{ref}/32$刻み、最大$8D_{ref}$まで走査し、
未保証区間をまたがない隣接valid点だけでbracketします。Zhao branch境界を越えたprobeは最後のvalid点との間を
縮小探索し、responseを外挿しません。bracket後はguard付きsecantと中点fallbackで解きます。有限なbracketを
最後まで縮小しても残差が許容値の8倍を超える場合は、残差が小さい方の有限端点をwarning付きで受理します。
signed scanは明示A/B/C branchだけに適用します。既定の`zhao_root_selection="require_unique"`では、`auto`が
seedで一意な物理解を返さない場合は、探索中にbranchを選ばずfail closedとします。したがってimplicit化だけでは
branch多重性を解消せず、強いPEではA/B/Cの事前scanが必要です。
局所軌道はbatch開始時の表面電荷から計算し、陰的終点のresponseを流入VDF、PE barrier、matching gaugeへ使います。
ambient吸収の総量は陰的応答へ、PEありではPE放出を設定した表面放出電流へ、PE returnを
「表面放出flux - 外部escape flux」へ正規化し、要素別のraw分布は維持します。陰的終点との整合性は PE あり・なしと
強い正負相殺を含む回帰 test で検証します。runtime は mesh 電荷の補償和から得た有限な commit 済み $Q/A$ を
次 batch の canonical な$D_H$とします。PE なしでは
PE target を生成しません。
これは$k=0$の時間刻み安定化であり、$k\ne0$の局所電位変化、backendの物理範囲、粒子samplingに対する
`batch_duration`の上限を除去しません。

`response_backend="zhao_online"`は`response_table_path`を禁止し、`zhao_branch="auto" / "a" / "b" / "c"`と
`zhao_root_selection="require_unique" / "minimum_energy" / "continuation"`を受理します。`continuation`は
`zhao_branch="a"`かつ`implicit_zero_mode=true`に限定します。table backendとstationary Zhaoは
`zhao_root_selection`を拒否します。各queryで$E_H=D_H/\epsilon_0$を境界条件とし、上流0 V・零電場へ接続する有限$H$の
Sagdeev A/B/C rootを解きます。これは壁面の零電流根ではなく、零電流条件を課さないcharge-driven responseです。
既定の`require_unique`では、`auto`が複数の物理解を検出した場合、または数値失敗により一意なbranchを
確認できない場合はfail closedとします。`minimum_energy`では検出候補の表面から無限遠までについて

$$
U=-\frac{\epsilon_0}{2}\int_0^\infty E^2\,dx
$$

を評価し、明示branchではbranch内、`auto`では検証済みA/B/C候補間で最小の$U$を選びます。候補branchの
数値失敗により集合を確定できない場合、または最小値が相対$10^{-6}$以内で縮退する場合はfail closedとします。
v1の複数根検出は有限個のmultistartから得た収束根のcluster判定であり、数学的なroot isolationではありません。
電位エネルギー比較も時間依存安定性の証明ではありません。

`continuation`の初回queryは`minimum_energy`と同じmultistartでType A rootを選びます。以後は最後にacceptedとなった
endpointの$(\phi_0,\phi_m,n_{e,\infty})$をNewton seedにします。候補とseedをType Aの対数未知数へ写像したときの
最大成分差が0.25以下なら局所Newtonの根を受理します。Newton失敗、rootのdecode失敗、profile検証失敗、または
この距離を超える場合だけfull multistartへ戻り、検出したType A rootのうちseedに最も近いものを調べます。
最近傍距離を$d_1$、2番目を$d_2$としたとき、$d_1\le0.25$かつ
$|d_2-d_1|\le10^{-6}\max(1,d_1)$なら、guessの検出順では選ばず曖昧状態として停止します。

最近傍rootも距離0.25を超える場合、Zhao evaluatorは解なしやfamily ambiguityではなく専用の
step-too-large statusを返します。implicit solverは有効なrootと棄却probeの間を二分し、距離条件を満たす
中間rootから探索を続けます。有効rootの直後で解なしまたは数値失敗となったprobeも、branch終端を粗く
飛び越えないよう同じ二分を試します。二分しても中間rootを再取得できない場合だけ停止します。
この規則はA/B/Cを暗黙に切り替えず、pseudo-arclength continuationのようにfoldの位置や通過を保証もしません。
`beach-zhao-response`にはaccepted endpointがないため、history-dependentな`continuation`を指定した表生成は拒否します。
最小エネルギー根の切替でresponseが不連続になり、backward-Euler残差が零点を持たず不連続だけをまたぐ場合は、
根を補間または混合せずnumerical failureとします。Newton multistartの検出順は物理的なroot IDとして公開しません。
branch別の物理検証では`a` / `b` / `c`を明示してparameter scanします。
ここで$H$は外部半無限領域のinterface原点、zero-mode gauge、PE moment測定面を固定します。平面・並進対称の
online closureでは$H$の絶対座標をSagdeev方程式の数値parameterにせず、壁面から$H$までの距離拘束は解きません。
外向きPE number fluxと平均法線energyは、その2 momentを再現するhalf-Maxwellianへ写像します。PE fluxが0なら
PE populationは0のまま、PE speciesがない場合はambient electron温度を数値scaleのfallbackに使います。

online MVPはambient electron / ionの外向きfeedbackをtransparentとして扱い、外部profile、戻りflux、応答値へ
反映しません。`require_unique`と`minimum_energy`の各queryはstatelessです。`continuation`はaccepted endpointのrootだけを
次batchのseedとして保持し、棄却したimplicit probe、固定点trial、adaptive trialをcommitしません。restartでは既存の
accepted responseからseedを再構成し、再構成できなければ初回multistartへ戻ります。outer inventoryとflight-time queueは
どのpolicyも保存しません。
設定したbranch policyで解が存在しない、branch制約を満たさない、Sagdeev積分が実数にならない、または非線形solveが
収束しない場合は停止します。明示したbranchやbackendを暗黙に切り替えません。stationary Zhaoの`solar_elevation_deg`、
`photoelectron_ref_density_m3`、`photoelectron_source_scale`はonline入力ではありません。

online implicitのPEありでは、各outer feedback反復の現在値$X^m$で後退Euler終点を解き、同じtrialの粒子追跡で
得たPE momentを緩和してから次の終点を解き直します。PE fluxが0のtrialでは平均energyは未定義なので、現在の
canonical energyを保ち、escape fluxを0とします。このnested反復によりPE returnと$D_H^{n+1}$を整合させます。

online Zhaoは平面・無衝突・非磁化、全roleの単価電荷、$T_e>0$、$0\le T_i\le0.1T_e$、
正の無限遠ion密度、ambient electron / ionの正の内向きdrift
（`drift_velocity`のz成分は負）を要求します。設定検査は
これらを満たさないcaseを拒否します。PE指定時はambient electronとPEの同一質量および$T_{pe}>0$も要求します。

各batch trialでは、表面総電荷とlower boundaryから

$$
D_H=D_b+Q_{cell}/A
$$

を得ます。accepted済みouter feedback $X^0$から、選択したresponse backendの評価、
$\Phi_H$を指定したinner field、
応答flux/barrierを使う粒子追跡、実際の$H$外向きmoment測定を反復します。同じbatch開始RNG stateと
macro-particle端数を毎反復で再生し、Monte Carlo写像を反復間で変えません。新規runの通常modeはゼロから
開始します。implicit tableはsingleton参照値、implicit onlineのPEありは設定した放出number fluxと$T_{pe}$、
PEなしはゼロを初期値とします。raw response $X_{raw}^{m+1}$は

$$
X^{m+1}=(1-\alpha)X^m+\alpha X_{raw}^{m+1},
\qquad 0<\alpha\le1
$$

で緩和します。feedback vectorと`coupling_atol`の成分順は、PE外向きnumber flux [m^-2 s^-1]、
PE外向き平均法線energy [eV]、ambient electron外向きnumber flux [m^-2 s^-1]、
ion外向きnumber flux [m^-2 s^-1]です。active軸$j$のbackend scaleを$s_j$、
`coupling_rtol`を$r$、`coupling_atol`の成分を$a_j$とすると、収束条件は

$$
\left|X_{raw,j}^{m+1}-X_j^m\right|\le\max(r s_j,a_j)
$$

です。`coupling_atol`は有限な非負4-vectorで、既定値`[0, 0, 0, 0]`は従来の相対許容値だけの判定を保ちます。
inactive軸の$a_j$は0でなければならず、非零値は設定検査またはprovider初期化で拒否します。
出力する`matching_plane_residual`は、$a_j>r s_j$の成分を$r|\Delta_j|/a_j$、それ以外を
$|\Delta_j|/s_j$として最大値を取るため、絶対許容値を使っても収束したtrialでは
`matching_plane_residual <= coupling_rtol`を保ちます。
tableはactive feedback軸のspanを$s_j$とし、singleton軸を除外します。onlineはZhao modelの基準flux / energyから
$s_j$を導出し、transparentなambient electron / ion outward軸を除外します。tableのactive軸が補間範囲外、
online solveが失敗した場合、または非有限値を含む場合はfail closedとします。
有限なtrialが`coupling_max_iterations`までに収束しない場合はwarningを出し、残差と最大反復回数をreceiptとして
最終trialをcommitします。次batchのouter stateは観測feedbackから開始します。
adaptive batch-durationで棄却したtrialはouter stateもbatch開始値へrollbackします。

z-highの外向きeventはspecies別にmacro weightを掛けて集約し、PEについては外向き数、法線energy、外部barrierでの
return数、escape数を独立に保持します。$\Gamma_{pe}^{out}=\Gamma_{pe}^{return}+\Gamma_{pe}^{escape}$を診断し、
outer state、反復回数、残差をhistoryとcheckpoint schema v9へ保存します。restartは保存したfeedbackから反復を
再開します。tableでは応答内容、onlineではZhao設定をfingerprintへ含めますが、不一致はwarning付きの
changed-condition continuationとして許可します。
onlineの自動bracket点は永続tableやmodel stateではなく、restart後に同じZhao contractから再評価します。

このmodelは準定常・無衝突・非磁化の低次元closureです。完全6D VDF、外部flight time、遅延return queue、
外部過渡、BEACH領域内volume plasma chargeは解きません。online Zhaoもfull VDF、1D PIC、time-dependent outer sheathの
代替ではありません。production tableは独立に検証したZhao/1D PIC sweepから作成し、どちらのbackendでも
matching-plane高度$H$をoverlap region内で変えたときの主要量の不変性を連成検証とします。

## 8. 実行制御

- 新規実行は `sim.batch_count` 個の accepted batch を処理する
- 再開は checkpoint の batch 数から `sim.batch_count` まで処理する
- 各粒子の前進上限は `sim.max_step`
- `tol_rel` は早期停止に使わない

## 9. 出力と再開

主要出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `mesh_potential.csv`（設定時）
- `charge_history.csv` / `potential_history.csv` / `top_reference_history.csv`（設定時）
- `matching_plane_history.csv`（matching-plane連成かつ履歴出力時）
- `charge_ledger.csv`（ledger がある場合）
- `rng_state.txt`、MPI では `rng_state_rankNNNNN.txt`
- `macro_residuals.csv`
- `checkpoint_complete.txt`
- `performance_profile.csv`（`BEACH_PROFILE=1`）

出力ファイルの条件は `schemas/beach.output-manifest.json` を正本とします。

再開時の必須ファイルは `summary.txt`、`charges.csv`、対応する RNG state です。
schema v8 以降は `checkpoint_complete.txt` 自体も必須で、manifest が宣言する
`macro_residuals.csv` と `charge_ledger.csv` も必須です。
summary に ledger metadata があれば、schema の世代によらず `charge_ledger.csv` の欠落を拒否します。

`output.checkpoint_stride > 0` では accepted batch の commit 後だけ定期 checkpoint を作ります。
`checkpoints/slot0` と `slot1` を交互に使い、全 rank の状態を書き終えて `checkpoint_complete.txt` を完了状態へ
原子的に切り替えてから、`checkpoint_latest.txt` を原子的に切り替えます。再開時は直下の最終出力と両 slot を検査し、
必須ファイルが揃う load 可能な checkpoint のうち `batches` が最大のものを選びます。
`checkpoint_latest.txt` が欠落、破損、または古い場合も、完了 manifest を持つ slot は回収対象です。
`checkpoint_stride=0` でも正常終了時の最終 checkpoint は出力します。

`summary.txt` の checkpoint schema と model / ordered mesh / ordered species fingerprint を照合します。
ordered mesh fingerprintの不一致は要素電荷を安全に対応付けられないため停止します。model / species fingerprintの
不一致はprovenance warningとして報告し、保存状態の配列形状が読込可能なら現在の条件で継続します。
schema v6 の `macro_residuals.csv` は `species_idx,face,residual` を持ち、`face=0` は従来 source、
`1..6` は boundary face です。旧 2 列形式は読み込み互換です。
schema v9はmatching-planeのaccepted feedback、potential、return/escape flux、反復receiptを`summary.txt`へ保存します。
non-matching modelではこれらを無効値として保持し、v8以前のload可能なcheckpointとの互換性を維持します。
globalまたはspecies境界のいずれかで`redistributed_reflect`を使う場合だけ、model fingerprintへ
`sim.rng_seed`と乱数契約識別子`redistributed_reflect_rng_v1`を含めます。また、境界event速度をchord方向かつ
予測中点電場の離散workと整合させる契約と、表面注入を未照会飛行なしで1 ULP内側から開始する契約にも
version tagを持たせます。
tree solverの`tree_theta`と`tree_leaf_max`は値だけでなく明示指定の有無もfingerprintへ含め、条件変更をsummaryと
warningから追跡可能にします。fingerprintは同一条件の証明とprovenanceに使い、mesh以外は再開禁止には使いません。
必須ファイルの欠落、world size の不一致、非有限値、保存ファイルが要求する配列形状や mesh 要素数の不一致は
新規実行へ fallback せず停止します。

## 10. 設計方針

- v1.0 の基本系は insulator accumulation
- correctness を優先し、性能機能は明示設定で有効化する
- 実装されていない物理モデルを設定互換性のためだけに保持しない
- public API の変更時は schema、例、日英ドキュメント、テストを同時更新する

## 11. Python 後処理

Python package `beach` は Fortran 出力を読み込み、電位・電場・Coulomb 力・可視化を提供します。
任意点の再評価は三角形 geometry を渡す native `triangle_p0` kernel を使います。
`summary.txt`のfield-reconstruction schema v2は、実際に解決された`direct` / `treecode` / `fmm`と
固定FMM展開次数を記録します。自動再構成は`direct`を非周期exact-directへ、`fmm`をreceiptの展開次数を
使うFMMへ対応させます。`treecode`は高水準の電場・電位・力APIでFMMへ置換せず停止します。
resolved `direct`の電場・電位・力はexact-directで評価し、uniform `E0`を一度だけ合成します。

`cached_kneq0` は非零 mode だけの低水準 operator であり、total field として使うには Fortran と同じ
物理的 zero mode の合成が必要です。設定が見つからない total-field 再計算は free-space へ暗黙 fallback しません。
