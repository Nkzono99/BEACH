title: matching-plane 準定常連成を使う

Lang: [日本語](MatchingPlaneCoupling.md) | [English](MatchingPlaneCoupling.en.md)

# matching-plane 準定常連成を使う

`surface_current_model.model="matching_plane_quasistatic"` は、BEACH の上端
`z = domain.box_max` の z 成分を外部 1D シースとの整合面にします。BEACH は全表面電荷と整合面より下の
`k=0` / `k!=0` 場を保持し、外部モデルの事前計算済み応答表から整合面電位、ambient 流入束、
およびアクセス・return 障壁を求めます。各 batch の粒子軌道が返す外向き flux と応答表を固定点反復し、
収束した trial だけを電荷更新へ commit します。

この機能は Zhao の定常壁電位を run 中に再評価するものではありません。応答表は Zhao または 1D PIC の
独立した parameter sweep から作成し、その外部モデルと適用範囲を別途検証する必要があります。
`examples/matching_plane_response_synthetic.csv` は実行経路を確認するための合成データであり、物理計算には使えません。

## 最小構成

**前提:** x / y 周期、z open の `periodic2` case を用意し、全 mesh 頂点を上端より厳密に下へ置きます。
外部一様場と一様磁場はゼロにし、electron、ion、photoelectron の 3 species を別々に指定します。

**操作:** `[surface_current_model]` を追加します。

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_table_path = "matching_plane_response.csv"
electron_species = "electron"
ion_species = "ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

完全な合成 smoke case は `examples/periodic2_matching_plane_quasistatic.toml` にあります。
`response_table_path` は `beach.toml` のあるディレクトリからの相対 path として解決されます。

**期待する出力:** 設定検査が model 構成を受理し、実行開始時に応答表の構文と整合面高度が検証されます。
`summary.txt` には最後に accepted された連成 state に加え、応答表の content fingerprint、3 role、反復設定が
記録されます。`output.history_stride > 0` なら、
該当する accepted batch ごとの応答値、4 つの feedback moment、return / escape flux、反復回数、
正規化残差を `matching_plane_history.csv` で追跡できます。

`beachx lint` は response file の内容を読みません。

**解釈:** run の正常終了は表の構文、補間範囲、および固定点収束を確認します。外部シースの物理的妥当性、
matching-plane 高度への不変性、Monte Carlo 収束は別の検証です。

## 構成上の制約

この model は二重計上を避けるため、次の構成だけを受理します。

- `field_boundary.mode="periodic2"`
- `periodic2.nonzero_mode_backend="cached_kneq0"` または `"panel_spectral_reference"`
- `periodic2.zero_mode_policy="exclude_k0"`
- `periodic2.lower_boundary_model="symmetric_vacuum"` または `"e_bottom_zero"`
- x / y の粒子境界は periodic、z-low / z-high は open
- `sim.e0 = [0, 0, 0]`、`sim.b0 = [0, 0, 0]`
- `[reservoir].inflow_model="source_vdf"`（generic な `infinity_barrier` は無効）
- enabled species は 3 つの role だけとし、相異なる各 role の `surface_charge_closure="explicit"`
- ambient electron / ion は `source_mode="volume_seed"`、`npcls_per_step=0` とし、
  `boundary_inflow` は z-high の `reservoir` だけ
- photoelectron は負電荷の `photo_raycast`
- photoelectron の z-high 外向き境界は open
- `particle_boundary.ordinary_open_model="escape"`
- `sim.multiple_box_events_policy="abort"`

`reference_area_m2`、Zhao 専用キー、手動 `fixed_current` target は併用できません。整合面積は
`domain` の x-y 面積、整合面高度は `domain.box_max` の z 成分、更新間隔は 1 accepted batch に固定されています。
これらを別 parameter にしないことで、互いに不整合な値を設定できないようにしています。

## 応答 CSV v1

header より前に、整合面高度を 1 回だけ指定します。値は有限で、`domain.box_max` の z 成分と一致しなければなりません。

```csv
# matching_plane_z_m=1.0e-3
displacement_c_m2,photoelectron_outward_number_flux_m2_s,photoelectron_outward_mean_normal_energy_ev,electron_outward_number_flux_m2_s,ion_outward_number_flux_m2_s,matching_potential_v,electron_inward_number_flux_m2_s,ion_inward_number_flux_m2_s,electron_access_potential_v,ion_access_potential_v,photoelectron_barrier_potential_v
```

最初の 5 列が入力軸、後ろの 6 列が応答です。

3 つの potential 出力と `matching_potential_v` は同じ gauge を使い、外部モデルの上流 reservoir を
0 V とします。BEACH は electron / ion の inward VDF をこの 0 V reservoir から access potential と
整合面電位へ写像するため、potential 列だけへ任意の定数を加えることはできません。

| 列 | 単位 | 意味 |
| --- | --- | --- |
| `displacement_c_m2` | C/m2 | BEACH の表面電荷と下端 closure から得る整合面直下の平均 $D_z$。+z を正とする |
| `photoelectron_outward_number_flux_m2_s` | 1/(m2 s) | 整合面へ到達した外向き PE number flux |
| `photoelectron_outward_mean_normal_energy_ev` | eV | 外向き PE の平均法線運動エネルギー |
| `electron_outward_number_flux_m2_s` | 1/(m2 s) | ambient electron の外向き number flux |
| `ion_outward_number_flux_m2_s` | 1/(m2 s) | ion の外向き number flux |
| `matching_potential_v` | V | 外部シースが返す整合面電位 $\Phi_H$ |
| `electron_inward_number_flux_m2_s` | 1/(m2 s) | 外部で折り返す成分も含む、BEACH へ入る electron の総 number flux |
| `ion_inward_number_flux_m2_s` | 1/(m2 s) | 外部で折り返す成分も含む、BEACH へ入る ion の総 number flux |
| `electron_access_potential_v` | V | electron reservoir から整合面へ入る VDF の access bottleneck |
| `ion_access_potential_v` | V | ion reservoir から整合面へ入る VDF の access bottleneck |
| `photoelectron_barrier_potential_v` | V | 外向き PE が越える最大外部 potential barrier |

入力 5 軸の完全な直積格子を 1 行ずつ記述します。行順は任意ですが、重複点と欠損点は許されません。
読込メモリは行数に比例し、MPI では各 rank が表を保持します。不要な軸 node を増やす前に必要メモリを見積もってください。
BEACH は最大 32 corner の多重線形補間を行い、複数 node を持つ軸では外挿せず範囲外を停止します。
数値 token は十進実数だけを受理し、Fortran の `/`、`2*0`、空欄などの list-directed 制御記法は拒否します。

flux、平均エネルギー、および出力 flux は非負です。feedback 軸 2--5 は初期状態を評価できるよう 0 を含めます。
ある入力への依存を外部 sweep で扱わない場合、その軸を 1 node にできます。1 node 軸は任意の有限 query を受理し、
その feedback を明示的に無効化します。これは clamp や暗黙の外挿ではありません。

## 固定点反復

accepted 済みの feedback を $X^0$ とし、各 trial で次を反復します。

1. 現在の表面電荷から $D_H$ を計算し、$(D_H, X^m)$ で応答表を補間する。
2. 応答の $\Phi_H$ を `periodic2` の zero-mode gauge に設定し、electron / ion の総 inward flux と access map、
   PE の外向き barrier を適用する。ambient の外向き粒子は局所反射せず、次の外部応答へ渡す。
3. 同じ batch 開始 RNG state と macro-particle 端数から粒子 trial を再生し、外向き moment $X_{\mathrm{raw}}^{m+1}$ を測る。
4. 応答表の各 feedback 軸幅で正規化した $X^m$ と $X_{\mathrm{raw}}^{m+1}$ の最大残差が
   `coupling_rtol` 以下なら、この trial を収束済みとして受理する。
5. 未収束なら `coupling_relaxation` を $\alpha$ として
   $X^{m+1}=(1-\alpha)X^m+\alpha X_{\mathrm{raw}}^{m+1}$ を計算し、次の replay へ進む。

履歴の応答値は受理した trial の $X^m$ で評価した値、feedback 列は同じ trial で実測した
$X_{\mathrm{raw}}^{m+1}$ です。両者の差は記録された残差以下であり、収束後に未実行の緩和 step は加えません。

反復中は表面電荷、RNG、注入端数、ledger、history を commit しません。adaptive batch-duration trial が棄却された場合も、
outer state を含めて batch 開始状態へ戻します。`coupling_max_iterations` までに収束しない場合や補間範囲を外れた場合は、
未収束値で継続せず停止します。

## 応答表を作る

1D Zhao または 1D PIC 側で、同じ `z=H`、上流分布、符号規約を使って入力 5 軸を sweep します。
壁面での放出 flux ではなく、整合面を実際に上向き通過した PE flux と法線 energy を入力にします。
非単調 potential では、局所的な電位差ではなく外部 profile 全体の最大 barrier を出力します。

表には少なくとも、想定する初期状態、定常状態、および反復の過渡点を余裕をもって含めます。依存性を扱わない軸は
singleton にし、5 軸を一律に細分化して表を膨張させないでください。ただし格子を広げることは
外部モデルの適用範囲を広げません。表生成 code、上流条件、solver version、単位変換を production data と一緒に保存してください。

## 妥当性を確認する

**数値検証:** 次を独立に変え、主要量が許容誤差内で安定することを確認します。

- `coupling_rtol`、`coupling_relaxation`、`coupling_max_iterations`
- response grid の刻みと範囲
- `sim.batch_duration` と macro-particle 数
- mesh / periodic cell の解像度

**連成検証:** 外部モデルを同じ物理条件で再生成し、overlap region 内で `H` を変えます。grain charge、gap potential、
PE escape fraction がほぼ変わらないことが matching-plane 連成の中心的な検証です。PE は整合面で
`outward = return + escape` となることも確認します。

**適用限界:** この準定常 model は 6D VDF、outer flight-time queue、衝突、磁化された return、外部シースの過渡、
BEACH 領域内の volume plasma charge を解きません。これらが結果を支配する場合は、1D--3D kinetic coupling または
full PIC を検証用・production 用に使います。
