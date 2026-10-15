# BEACH 仕様書（現行 Fortran 実装）

## 1. 目的

BEACH は、三角形境界要素上の電荷蓄積とテスト粒子追跡を行うシミュレータです。

- 境界は三角形メッシュ
- 粒子運動は Boris push
- 境界衝突時は吸収して要素へ電荷を堆積
- `batch_count` 個の accepted batch を処理し、統計と履歴を更新

計算の主系は Fortran（`src/`, `app/`）で、Python（`beach/`）は設定補助、後処理、可視化を担います。

## 2. スコープ

### 2.1 サポートする構成

BEACH がサポートする製品構成は、2 軸周期セルに外部 1-D シースを接続した照射実験の 1 種類です。
判断の経緯は [docs/adr/0002-irradiation-cell-scope.md](docs/adr/0002-irradiation-cell-scope.md) にあります。

| 要素 | 製品構成 |
| --- | --- |
| 領域 | `[domain]` box、`periodic_axes = ["x", "y"]` |
| 場 | `periodic2`、FMM、`cached_kneq0`、split zero mode (`exclude_k0`)。direct は `panel_spectral_reference` の検証用 |
| 表面 | insulator accumulation、`mesh_id` ごとの浮遊 conductor 再配分 |
| 粒子境界 | z は `open`。`reflect` / `redistributed_reflect` は closed photoelectron 用 |
| 粒子源 | `volume_seed`、boundary reservoir 流入、`photo_raycast`、衝突時の `secondary_emission` |
| closure | `explicit`、`fixed_current`、`neutral_return`、`zhao_stationary` と outflow refresh |

参照構成として、free-space の direct / FMM と入門ケース (free + direct + `volume_seed` + plane template) の
ビルドとテストを維持します。参照構成に新しい機能は加えません。それ以外の組合せは lint と Fortran が拒否します。

### 2.2 実装済み

- template / OBJ の三角形メッシュ
- 三角形上の一定面密度を使う `triangle_p0` 静電場
- direct / treecode / FMM field solver
- free-space と 2 軸周期 `periodic2`
- 一様外部電場 `fields.external.electric_v_m` と磁場 `fields.external.magnetic_t`
- Boris 法による粒子更新
- 線分と三角形の最初の交差判定
- insulator accumulation
- mesh 単位の浮遊 conductor 再配分 (free-space と periodic2 の snapshot 電位で等電位化)
- open (脱出) / reflect / redistributed_reflect / periodic の box 境界と Zhao の種別 kinetic 障壁
- `volume_seed` / `photo_raycast` / 吸収時の `secondary_emission`
- species 別の simulation boundary reservoir 流入
- 境界 reservoir の流入補正 (`source_vdf` / `infinity_barrier`)
- closed photoelectron の局所反射と neutral-return closure
- 外部 1-D シースの Zhao 零電流根による表面電流 closure と、観測 PE + SEE 流出による外部根の準定常更新
- checkpoint 再開

### 2.3 未実装・予約

- BEACH 内で解く full-VDF / 1D PIC / time-dependent outer-sheath solver と外部領域の粒子 transport
- 表面導電・拡散モデル
- 誘電体分極・誘電境界条件
- 反射、材質ごとの二次電子応答を含む一般化衝突モデル

### 2.4 設定入力の検証

公開設定は `beach.toml` です。未知の table/key、基本型の不一致、非有限の実数・配列成分、
Fortran の格納先の長さを超える文字列は読み込み時に拒否します。文字列を暗黙に切り詰めません。
整数は符号付き 32 bit の `-2147483648..2147483647` に制限し、範囲外は変換前に拒否します。
無効化した species も基本値の検査対象です。

検証の正本は 2 つです。キー名、型、列挙値、値域は `schemas/beach.schema.json` が持ちます。
組合せの制約と派生値は Fortran の preflight が持ち、Python には同じ制約を書きません。
`beachx lint beach.toml` は TOML の読込、schema の検査、座標・配置の正規化の後に `beach --check-config` を実行し、
Fortran の診断をそのまま表示します。通常運用では `beachx lint` に成功してから `beach beach.toml` を実行します。

`beach --check-config beach.toml` は通常実行と同じ既定値・読取・正規化と preflight を適用します。
成功時は終了コード 0 と `status=ok`、検出した入力エラーには非ゼロの終了コードを返します。
MPI 初期化、mesh 構築、外部データの内容読込、simulation、結果・checkpoint 出力は行いません。

Python の `load_config_file`、`normalize_config_document` は schema と authoring の正規化だけを行います。
contract test は、代表ケースごとに `beachx lint` が Fortran の採否をそのまま返すことを確認します。
`beachx lint --schema` は同梱契約に制約を追加し、正規化前後に検査します。同梱契約を緩める指定には使えません。

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

`output.diagnostics.charge_rel_change_threshold` は監視・出力値であり、早期停止条件ではありません。

## 4. 1 batch の計算手順

1. 現在の要素電荷から静電 snapshot を構築する
2. boundary reservoir inflow、`volume_seed`、`photo_raycast` の粒子を生成する
3. 各粒子を最大 `particles.tracking.max_steps_per_particle` 回まで進める
4. mesh hit、box escape、max-step survivor を分類する
5. 吸収電荷と放出元の逆符号電荷を thread-local buffer へ集計する。SEE は世代ごとに追加追跡する
6. closed PE が有効なら、解決済み帰還先へ neutral-return 補正を適用する
7. MPI reduction 後の候補電荷を検証して commit する
8. accepted batch の統計、履歴、charge ledger、累積時間を更新する

`run.batch.adaptive.max_nonzero_mode_potential_step_v > 0` の場合は、設定された `run.batch.duration_s` を最大幅として
`h, h/2, h/4, ...` を試します。候補電荷による cached $k\ne0$ 電位変化が上限以下の trial だけを
commit します。棄却時は RNG と macro 粒子数残差を trial 前へ戻し、統計・履歴・ledger を更新しません。

## 5. 場モデル

### 5.1 free-space

`fields.solver.method="direct" | "treecode" | "fmm" | "auto"` を選べます。

- direct: 厳密な free-space panel 和
- treecode: near は厳密 panel、far は monopole
- FMM: near は厳密 panel、far は三角形 moment を使う FMM
- auto: 要素数に応じて direct または FMM

面上電位は連続で、法線電場は `sigma/epsilon0` だけ跳びます。重心電位と principal-value 電場を
自己項として使います。`surface_side` は各 mesh/template で明示します。

### 5.2 periodic2

`[fields] boundary="periodic2"` は `[domain] periodic_axes=["x", "y"]` を要求します。
mesh は primitive cell の base element だけを保持し、場と衝突判定が周期 image を扱います。

`fields.periodic.backend`:

- `finite_images`: 指定した有限 image shell のみ
- `cached_kneq0`: 無限周期の非零 Fourier mode を versioned operator cache から評価
- `panel_spectral_reference`: 小規模な direct solver 検証で非零 Fourier mode を直接評価

型付き `[fields.periodic]` 設定では、小規模な検証計算に限り
`fields.periodic.backend="panel_spectral_reference"` と direct solver の組合せも利用できます。
これは非零 Fourier mode の収束参照であり、外部プラズマモデルには依存しません。

`fields.periodic.backend="finite_images"` は有限画像近似であり、無限周期解ではありません。
`m2l_root_oracle` は削除済みです。production の無限周期計算には `cached_kneq0` を使用します。
cache fingerprint は generator version を含み、物理的 zero mode の評価は高さ breakpoint の
二分探索により 1 点あたり O(log n) です。

平均場の高さ分布は保持し、全電荷源の上側では累積電荷を総電荷で厳密に固定します。

`cached_kneq0` の遠方補正は、線源の木の深さ 2 の節点ごとの multipole から深さ 2 の target anchor の local 展開への
operator を当てはめます。有限画像殻の面平均は無限面電荷の平均場と滑らかな残り $s$ の和で、$s$ は閉じた式で
求めて線源ごとの高さの表から実行時に引き、operator には面平均 0 の残差だけを当てはめます。こうして非零 mode の
面平均を 0 に保ちます。

`fields.periodic.backend="cached_kneq0"` または `"panel_spectral_reference"` では、物理的 $k=0$ を
`fields.periodic.lower_boundary_model="symmetric_vacuum" | "e_bottom_zero"` で別に構築し、場へ一度だけ加えます。
このzero-mode構築自体は外部プラズマ応答を合成しません。外部シースclosureを選んだ場合、その壁電位は
z-high面の電位基準として別の層で与えます。

### 5.3 領域・粒子境界・reservoir

公開 TOML は `run / domain / mesh / particles / fields / sheath / output` の7グループ形式だけを受け付けます。
設定キーは以下も含めて公開 TOML のパスで表します。
再開・粒子源・周期場 backend に暗黙の優先順位を作りません。

公開 TOML は topology、場、外向き粒子作用、外部 reservoir 条件、species 別流入を分離します。

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]
periodic_axes = ["x", "y"]

[fields]
boundary = "periodic2"

[particles.boundary]
z_low = "open"
z_high = "open"

[particles.reservoir]
inflow_model = "source_vdf" # または "infinity_barrier"
phi_infty_v = 0.0
face_potential_grid_n = 5

[[particles.species]]

[particles.species.distribution]
number_density_cm3 = 5.0
temperature_ev = 10.0

[particles.species.inflow]
z_high = "reservoir"
```

- `inflow`: 非周期 box 面全体から外部 reservoir VDF を流入させる
- `source_vdf`: 指定した boundary VDF を補正せず境界から注入する
- `infinity_barrier`: `particles.reservoir.phi_infty_v` と流入面平均電位の scalar barrier で法線 VDF を補正する
- 通常の open 面到達で粒子を消滅する。Zhao の種別 kinetic 障壁がある面では外部根の電位で帰還／脱出を判定する

`inflow` は外向き粒子作用を上書きしません。周期面への reservoir 流入、および有効な外向き作用が
`open` でない流入面は拒否します。
境界流入だけの species では `source.mode` と `sampling.volume_macro_particles_per_batch` を省略します。体積生成は行いません。Maxwell 分布で体積生成を追加する場合だけ
`source.mode="volume_seed"` と正の `sampling.volume_macro_particles_per_batch` を明示します。
同じ species では解決後の `source.mode="volume_seed"` とだけ併用でき、`photo_raycast` との併用は拒否します。

外部プラズマ profile、外部領域の particle transport、delayed return queue は現行スコープ外です。
削除済みの `[outer_plasma]`、`[coupling]`、`[external_boundary]` は unknown input として拒否します。

## 6. 粒子前進・衝突・box 境界

public な粒子状態は同一時刻の `x, v` です。空間電場は予測中点
`x_mid = x + 0.5*v*particles.tracking.dt_s` で評価し、Boris の速度更新後に
`x_new = x + 0.5*(v + v_new)*particles.tracking.dt_s` とします。

mesh 衝突は Möller–Trumbore 法で最小 `t` を選びます。AABB と、既定では一様 grid + 3D-DDA で
候補を絞ります。非有限値や進行不能は hit なしへ黙って変換せず、明示 status で停止します。

`domain.periodic_axes` は全粒子種に共通の topology です。非周期面の粒子 action は
global `[particles.boundary]` と species 別 `[particles.species.boundary]` で指定します。

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

各 species で `sampling.volume_macro_particles_per_batch` 個を batch ごとに生成します。

`inflow` を持つ species では `sampling.volume_macro_particles_per_batch=0` を許容します。

### 7.2 species 別 boundary reservoir 流入

`[particles.species.inflow]` の 6 面キーへ `reservoir` を指定し、非周期 box 面全体から流入させます。
複数面を指定でき、macro 粒子の端数は species と face の組ごとに batch 間で繰り越します。
MPI では global 個数を rank へ分配します。

流入 flux と面積、batch duration から macro 粒子数を決めます。`sampling.target_macro_particles_per_batch` 指定時は
macro weight を自動解決します。

`particles.reservoir.inflow_model="infinity_barrier"` では各流入面の `N x N` 電位標本を使い、
`particles.reservoir.phi_infty_v` とのエネルギー差から法線速度 cutoff と到達速度を同じ写像で補正します。

### 7.3 `photo_raycast`

指定方向へ `sampling.rays_per_batch` 本を飛ばし、各 ray の最初の mesh hit からだけ放出します。
`source.deposit_opposite_charge_on_emit=true` では放出元へ逆符号電荷を加えます。

closed PE は次の組合せです。

- 負電荷 `photo_raycast`
- `source.deposit_opposite_charge_on_emit=true`
- species の `source.inject_face` に対する有効粒子 action が `reflect` または `redistributed_reflect`
- `charging.closure="neutral_return"`

放出電荷を $S<0$、解決済み再吸収を $R<0$、max-step 未解決を $U<0$ とし、`S=R+U` を確認した上で、
帰還先への deposit を `S/R` 倍します。これにより PE species の表面総電荷増分を 0 にしつつ、
帰還先分布による表面内再分配を保持します。

未解決率 `U/S` が 0.05 を超える、`R=0`、実 escape、soft discard、terminal 分類不整合のいずれかでは
補正せず fail closed とします。

### 7.4 species 別 fixed-current closure

`charging.closure="fixed_current"` は、species の吸収 channel と `photo_raycast` の放出反作用 channel を
独立に補正します。signed target charge は accepted trial の batch 幅 $\Delta t$ に対して

$$
Q_{s,\mathrm{abs}}^{\mathrm{target}}=I_{s,\mathrm{abs}}^{\mathrm{target}}\Delta t,
\qquad
Q_{s,\mathrm{emit}}^{\mathrm{target}}=I_{s,\mathrm{emit}}^{\mathrm{target}}\Delta t
$$

です。raw channel 電荷を $R$ とすると、各要素の raw deposit を同じ $Q^{\mathrm{target}}/R$ で倍率化し、空間分布を
保持します。`charging.target_absorbed_current_a` の符号は `charge_c` と一致し、
`charging.target_emission_current_a` の符号は `charge_c` と逆です。非ゼロ target に対して $R=0$、非有限倍率、負倍率なら
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

`[sheath]`は、設定型、model dispatch、model固有solver、species channelとkinetic境界写像への割当を分離します。
未指定または`sheath.closure="none"`はtargetを生成せず、speciesに記述した手動targetを使います。初期実装の
`sheath.closure="zero_current"`は、Zhao A/B/Cの平面・無衝突・非磁化シースについて零電流定常根をrun開始時に解きます。
`outflow_refresh_batches=N>0`では、z-highで観測したPE流出を外部シースの放出源として、accepted batch $N$個ごとに根を解き直します。
根は外部ライブラリsheath-model（`fpm.toml`でcommit固定）で解き、プロファイル全域の$E^2\ge0$、中性で電場0の上流への接続、
イオン流の通過を満たす根だけを採用します。成立する根がなければ起動時に停止し、外部根の更新では前回の根を保ちます。
電子が内向きにdriftし反射電子を持つA/Cの根は、上流のごく近くで$E^2<0$になるため既定では採用しません。
`sheath.zhao.upstream_band_tolerance`$=f>0$では、$E^2<0$が上流に接する帯に限られ、その電位幅がAでは$f|\phi_m|$、Cでは$f|\phi_0|$以下の根を採用し、
幅を`surface_current_model_upstream_negative_band_V`に記録します。`sheath.zhao.upstream_field_floor_v_m`$=E_{floor}>0$では、帯の内側の場の不足$\sqrt{-\min E^2}$が$E_{floor}$以下の根を採用し、不足を`surface_current_model_upstream_field_deficit_V_m`に記録します。これは上流境界を「場が外部レベルまで落ちた点」に置くことに相当し、太陽風の運動電場$|\mathbf{u}\times\mathbf{B}|$が目安です。両方を指定すれば両方を満たす根だけを採用します。設定例は電子にイオンと同じdriftを与え、$E_{floor}$で受理します。
新しい電流modelは同じdispatch resultへspecies別の吸収・放出targetと診断値を返すことで追加します。

Zhao modelはambient electronとionを参照し、PE有効時はphotoelectronも明示的に参照します。ionの温度と`sheath.zhao.ion_pressure_factor`は外部根のイオン流体（Sana & Mishra 2026）に渡します。各speciesは
`charging.closure="fixed_current"`を要求し、手動targetとの併用を禁止します。単価電荷、electron/PEの同一質量、
z-high reservoirからのambient流入（electron driftの内向き成分は0以上、ionは正。electronのdriftがあるとType A/Cは成立しないため0を推奨）、負電荷`photo_raycast`の放出反作用、PEのopenなz-high境界、イオン入口の超音速条件（$m_i u_z^2/e>$`sheath.zhao.ion_pressure_factor`$\cdot T_i$）を
fail-closedに検証します。PEの`source.emit_current_density_a_m2`（レイ垂直面基準）× ray仰角のsinが外部源から5%以上、またはray仰角が太陽高度から1度以上ずれる設定は警告します。
非磁化closureなので`fields.external.magnetic_t`はゼロを要求します。Zhao固有の0 V reservoirと速度写像を使うため、genericな
`particles.reservoir.inflow_model="infinity_barrier"`との併用も拒否します。
電流密度から電流へ変換する面積は`reference_area_m2`、省略時はdomainのx-y面積です。

`sheath.photoelectrons.source_scale=0.0`は光電子なしの契約です。この場合はPE固有設定とPE speciesを省略し、
`sheath.zhao.branch`は`"auto"`または`"c"`を指定します。Zhao Type Cの$J_e+J_i=0$を解いて
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
周期k=0成分を持つ明示的なsplit `[fields.periodic]`では、z-high面の水平平均電位を根の壁電位$\phi_0$へ固定します。
障壁と流入写像は上流0 Vからの電位差であり、z-high面は外部1-Dシースから見た壁だからです。このため
`zhao_stationary`は`fields.boundary="periodic2"`、x/y周期、明示的なsplit `[fields.periodic]`（`cached_kneq0`または
`panel_spectral_reference`、`exclude_k0`）を要求し、split zero modeを持たない場ソルバーは設定エラーにします。

PE放出速度は従来どおり設定した表面half-Maxwellianです。PE、ambient electron、ionがopenなz-highを外向きに
横切るとき、粒子位置の局所電位から残りのZhao barrierを評価します。electron/PEのbarrier電位はType Aで
$\phi_m$、Type B/Cで0 V、ionは0 Vです。法線運動エネルギーが不足する粒子は法線速度を反転してreturnへ、
十分な粒子だけをescapeへ分類します。接線速度は保ち、x/yがともに周期ならreturn位置をz-high面内で一様に
選び直します。外部シース内の横移動はセル幅よりはるかに大きいためです。一様乱数は粒子eventのcounterから作り、
OpenMP実行順序に依存しません。選び直した位置では、横切り位置との局所電位差だけ法線運動エネルギーを補正して全エネルギーを保ち、
補正後に届かない低速粒子は横切り位置へ戻します。x/yが周期でなければ横切り位置へ戻します。固定電流targetはその後に各channel総量を正規化するため、kinetic写像は
rawなreturn/escape分類と空間分布を決め、Zhao電流は総電流収支を決めます。

`outflow_refresh_batches=N>0` は上の場の条件に加えて PE または SEE の放出電子グループを要求します。以下の PE 規格化に SEE の重み由来の束を合算します（7.1 節）。
$N$個のaccepted batchで、z-highを外向きに横切ったPE数$N_{out}$、その法線運動エネルギー和、表面放出数
$N_{emit}$をraw weightで合計し、$\eta=N_{out}/N_{emit}$と平均法線エネルギー$\bar K$を得ます。$N_{out}$は法線運動
エネルギーのbin（$[0,1\,\mathrm{meV})$と1 meV〜1 keVの1桁24個の対数bin、1 keV以上は最後のbin）別にも数えます。
外部根は放出源を、全体が$\eta\Gamma_{emit}$になるよう比例配分したこのbinスペクトル（bin内で$d\Gamma/dK$一定）として
解き直します。流出は放出温度に近い本体とセル内で加速された裾からなり、平均は裾で決まるため、同じ平均の
half-Maxwellianへは縮約しません。表面放出$J_{emit}$は設定から決まる値のままで、
$J_{return}=J_{escape}-J_{emit}$はセル内再吸収と外部returnの和です。`sheath.zhao.branch`が明示ならそのbranchだけ、
`"auto"`なら直前に採用したbranchを先に試します。根が得られない窓、または流出0の窓は警告して直前の根を保ちます。
採用した外部状態は`matching_plane_*`のstatsへ保存し、`matching_plane_history.csv`へ、放出源のスペクトルは
`matching_plane_spectrum_history.csv`へ出力します。checkpointからの再開では保存した放出源（スペクトル、無ければ
half-Maxwellian）から根を再構成し、途中の窓は再開時点から数え直します。

このstationary modelはbox外の電場・空間電荷・Debye shielding・return軌道・遅延を解かず、
`outflow_refresh_batches=0`ではrun中の表面電位に応じてtargetを再計算しません。z-high反射は外部turning pointまでの距離・飛行時間を省略するadiabaticな境界closureです。
外部シースの過渡解ではなく、BEACHの軌道追跡から得る空間分布へ固定総電流を与えるclosureです。

### 7.1 二次電子放出

入射種の `particles.species.surface_interaction.secondary_emission` で `enabled`、`model`、`yield_max`、
`energy_max_ev`、`const_yield`、`target_species` を指定します。既定は無効、有効時は放出先必須です。
Whipple は正の最大収率・最大エネルギー、constant は 0 以上の収率を要求します。
入射エネルギー・入射側法線との角度から収率を求め、角度の cosine を 0.1 に制限します。
`delta*abs(q_in/q_see)` を整数部 + Bernoulli 端数で放出し、重みは入射粒子と同じです。
位置は衝突点 + 法線方向 1e-12 m、周期座標を canonicalize します。速度は PE と共有する半 Maxwell サンプラーを
放出先温度・drift 0 で使います。吸収電荷の後に放出元へ必ず `-q_see*w_in` を置きます。

放出先の `source.mode="secondary_emission"` は負電荷・正質量・0 以上の温度を持ち、独立した
`sampling.*`、`inflow.*`、放出電流密度、光線数、`charging.closure` を指定できません。
thread ごとの buffer に生成粒子を積み、OpenMP ループ後に親・放出順で追加して同一バッチ・同一凍結場で追跡します。
一次粒子の世代を 0 とし、`particles.secondary_emission.max_generations`（既定 3、0 以上）に達した命中は
さらに放出せず、放出予定個数を記録します。追加分も processed/fresh 粒子、ledger、上端スペクトルへ含めます。
粒子状態は次のバッチへ持ち越さず、生成した rank に留まります。counter-based RNG は実行順に依存しません。

外部シースの `sheath.species.secondary_electron` は SEE を PE と同じ z-high 運動論的障壁と面内一様帰還へ接続し、
吸収・放出・脱出の固定目標は付けません。外部根の専用 SEE 項は作りません。更新窓の PE 流出は設定放出束 × 透過率、
SEE は `sum(w_out)/(reference_area*sum(run.batch.duration_s))` へ規格化し、bin と平均エネルギーを束で合算します。
PE なしの初期根は源 0 の Type C、SEE は最初の更新から入ります。
Phase 1 は暫定で、PE の脱出目標を窓の PE 束割合で配分し、帰還目標を脱出 − PE 放出にします。
SEE は完全に追跡し、追跡脱出と外部根の SEE 配分との差を `secondary_budget_residual_a_m2` へ出します。
モデルと適用範囲は [docs/models/SecondaryElectronEmission.md](docs/models/SecondaryElectronEmission.md) を参照してください。

## 8. 実行制御

- 新規実行は `run.batch_count` 個の accepted batch を処理する
- 再開は checkpoint の batch 数から `run.batch_count` まで処理する
- 各粒子の前進上限は `particles.tracking.max_steps_per_particle`
- `output.diagnostics.charge_rel_change_threshold` は早期停止に使わない

## 9. 出力と再開

主要出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `mesh_potential.csv`（設定時）
- `charge_history.csv` / `potential_history.csv` / `top_reference_history.csv`（設定時）
- `matching_plane_history.csv` / `matching_plane_spectrum_history.csv`（`zhao_stationary`の外部根更新かつ履歴出力時）
- `fixed_current_history.csv`（固定電流 closure があり履歴出力時。species・channel 別の追跡電荷と目標電荷）
- `charge_ledger.csv`（ledger がある場合）
- `secondary_emission_history.csv`（SEE 有効な入射種がある場合。受理バッチごと）
- `rng_state.txt`、MPI では `rng_state_rankNNNNN.txt`
- `macro_residuals.csv`
- `checkpoint_complete.txt`
- `performance_profile.csv`（`BEACH_PROFILE=1`）

出力ファイルの条件は `schemas/beach.output-manifest.json` を正本とします。

再開時の必須ファイルは `summary.txt`、`charges.csv`、対応する RNG state です。
schema v8 以降は `checkpoint_complete.txt` 自体も必須で、manifest が宣言する
`macro_residuals.csv` と `charge_ledger.csv` も必須です。
summary に ledger metadata があれば、schema の世代によらず `charge_ledger.csv` の欠落を拒否します。

`output.checkpoint.stride_batches > 0` では accepted batch の commit 後だけ定期 checkpoint を作ります。
`checkpoints/slot0` と `slot1` を交互に使い、全 rank の状態を書き終えて `checkpoint_complete.txt` を完了状態へ
原子的に切り替えてから、`checkpoint_latest.txt` を原子的に切り替えます。再開時は直下の最終出力と両 slot を検査し、
必須ファイルが揃う load 可能な checkpoint のうち `batches` が最大のものを選びます。
`checkpoint_latest.txt` が欠落、破損、または古い場合も、完了 manifest を持つ slot は回収対象です。
`output.checkpoint.stride_batches=0` でも正常終了時の最終 checkpoint は出力します。

`summary.txt` の checkpoint schema と ordered mesh fingerprint を照合します。
メッシュの不一致は保存電荷を同じ要素へ対応付けられないため停止します。モデル・粒子種の設定はハッシュ化せず、
保存状態の配列形状が読込可能なら現在の設定で継続します。
schema v6 の `macro_residuals.csv` は `species_idx,face,residual` を持ち、`face=0` は従来 source、
`1..6` は boundary face です。旧 2 列形式は読み込み互換です。
schema v9は`matching_plane_*`の外部状態（放出源、壁電位、return/escape flux、採用根数と相対変化）を`summary.txt`へ
保存します。外部根を更新しないmodelではこれらを無効値として保持し、v8以前のload可能なcheckpointとの互換性を維持します。
現行 schema v10 は model / species / response-content fingerprint の出力と照合を廃止します。
v2〜v9 も引き続き読み、旧 fingerprint は無視します。メッシュ識別子の計算法と状態配列の形式は変更しません。
以前の実行条件は設定ファイルと summary の具体的な値で確認してください。
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
