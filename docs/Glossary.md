title: 用語集

Lang: [日本語](Glossary.md) | [English](Glossary.en.md)

# 用語集

BEACH のドキュメントで使う用語と、対応する英語、設定キー・出力名の対応です。
本文ではここの日本語を使い、設定や出力を探すときは識別子の列を使います。
設定キーは新形式（7 グループ）の名前です。

## 計算の進め方

| 用語 | 英語 | 識別子 | 意味 |
|---|---|---|---|
| バッチ | batch | `run.batch_count` | 場を固定して粒子を追跡し、電荷差分を表面へ確定反映するまでの 1 周期 |
| バッチ幅 | batch duration | `run.batch.duration_s` | 1 バッチが表す物理時間。時間に比例する粒子源の注入量と、表面電荷を更新する間隔を決める |
| 確定反映 | commit | — | バッチ末尾に、集めた電荷差分を表面電荷へ一度だけ加えること |
| 時間刻み | time step | `particles.tracking.dt_s` | 粒子軌道を 1 回進める時間 |
| 適応進行 | adaptive progression | `run.batch.adaptive.*` | 面内変動成分の電位変化を上限内に抑えるよう、バッチを短い幅で再試行すること |
| 再開 | restart | `run.restart` | checkpoint の状態から計算を続けること |

## 粒子

| 用語 | 英語 | 識別子 | 意味 |
|---|---|---|---|
| 粒子種 | species | `particles.species[].species_key` | 電荷・質量・分布・粒子源を共有する粒子の集まり |
| マクロ粒子 | macro-particle | — | 実粒子の集まりを代表する計算粒子 |
| 重み | weight | `sampling.weight` | 1 マクロ粒子が代表する実粒子数 |
| 粒子源 | particle source | `source.mode` | 粒子を生成する方法（`volume_seed`、`plane_source`、`photo_raycast`） |
| 境界流入 | boundary inflow | `particles.species[].inflow` | box 外の外部プラズマから、非周期の box 面を通して粒子を入れること |
| 外部プラズマ（reservoir） | reservoir | `particles.reservoir` | box 外にあると仮定する、分布が既知のプラズマ |
| 流入の写像 | inflow mapping | `particles.reservoir.inflow_model` | 外部プラズマの分布を、流入面での分布に置き換える方法（`source_vdf`、`infinity_barrier`） |
| 吸収 | absorption | `absorbed_*` | 粒子が表面の三角形に当たり、電荷を残して消えること |
| 脱出 | escape | `escaped_*` | 粒子が開放面（`open`）から出て、戻らないこと |
| 反射 | reflection | `reflect`、`redistributed_reflect` | box 面で法線速度を反転して追跡を続けること |
| 帰還 | return | — | 表面から放出した粒子が、再び表面に吸収されること |
| 未解決粒子 | unresolved particle | `discarded_unresolved_*` | 1 粒子あたりの step 上限までに吸収も脱出もしなかった粒子。バッチ末尾に破棄する |
| 光電子 | photoelectron（PE） | `source.mode = "photo_raycast"` | 照射された表面から放出する電子 |

## 表面と電荷

| 用語 | 英語 | 識別子 | 意味 |
|---|---|---|---|
| 表面モデル | surface model | `surface_model` | 吸収した電荷を表面でどう保持するか（`insulator`、`conductor`） |
| 反作用電荷 | reaction charge | `source.deposit_opposite_charge_on_emit` | 粒子を放出した要素に残す、放出粒子と逆符号の電荷 |
| 電荷の閉じ方 | charge closure | `charging.closure` | 粒子種の電荷を、追跡結果に加えて何に合わせるか（`neutral_return`、`fixed_current`） |
| 閉じた光電子 | closed photoelectrons | `charging.closure = "neutral_return"` | 光電子を上端で反射し、正味の光電子電流を 0 に合わせる構成 |
| 目標電流・目標電荷 | target current / charge | `charging.target_*_current_a`、`fixed_*_target_charge_C` | 固定電流で合わせる電流と、それにバッチ幅を掛けた電荷 |
| 倍率 | weight scale | `*_weight_scale`、`target_over_tracked` | 目標電荷と追跡で得た電荷の比。追跡した分布をこの倍率で拡大・縮小して表面へ配る |
| 電荷台帳 | charge ledger | `charge_ledger.csv` | 粒子種ごとの注入・吸収・脱出・補正の電荷と個数の集計 |

## 外部シース

| 用語 | 英語 | 識別子 | 意味 |
|---|---|---|---|
| 外部シースとの接続 | sheath closure | `sheath.closure = "zero_current"` | box 外に平面定常シースを仮定し、その零電流根から上端の電流・障壁・電位基準を与えること |
| 零電流根 | zero-current root | `sheath.zhao.branch` | 表面への正味電流が 0 になる外部シースの定常解（Zhao Type A / B / C） |
| 壁電位 | wall potential（$\phi_0$） | `matching_plane_history.csv` の `phi_H_V` | 外部シースから見た表面（BEACH の上端面）の電位。上流のプラズマを 0 V とする |
| 電位極小 | potential minimum（$\phi_m$） | — | Type A で、上流と壁の間にできる電位の最小値。電子と光電子の障壁になる |
| 上流電子密度 | upstream electron density（$n_{e,\infty}$） | — | 外部シースの根が与える、無限遠の電子 Maxwell 分布の密度 |
| 外部根の更新 | outflow refresh | `sheath.coupling.outflow_refresh_batches` | 上端面を出る光電子の観測値から、外部シースの根を定期的に解き直すこと |
| 上端面 | top face（z-high） | `z_high` | box の上側の面。外部シースと接続する面 |
| 電位障壁 | potential barrier | `particles.boundary.open_model = "potential_barrier"` | 開放面を出る粒子を、上流電位との差で反射か脱出かに分ける判定 |

## 周期場

| 用語 | 英語 | 識別子 | 意味 |
|---|---|---|---|
| 場ソルバ | field solver | `fields.solver.method` | 表面電荷から場を計算する方法（`direct`、`treecode`、`fmm`） |
| 場境界 | field boundary | `fields.boundary` | 場をどの空間で解くか（`free`、`periodic2`） |
| 有限画像 | finite images | `fields.periodic.image_layers` | 周期セルの周囲に陽に並べる複製の層 |
| 面内変動成分 | nonzero mode（$k\ne0$） | `fields.periodic.backend` | x/y 方向に変化する場の成分 |
| 面平均成分 | zero mode（$k=0$） | `fields.periodic.lower_boundary_model` | x/y 方向に平均した場の成分。各高さより下の総電荷と下側境界条件で決まる |
| 遠方補正 | far correction | `fields.periodic.backend = "cached_kneq0"` | 有限画像の外側の周期複製が作る面内変動成分を、事前に作った演算子で加えること |

## 出力

| 用語 | 英語 | 識別子 | 意味 |
|---|---|---|---|
| 実行記録 | run summary | `summary.txt` | 実行統計と、解決した設定を `key=value` で書いたファイル |
| 履歴 | history | `*_history.csv` | バッチごとに書く時系列 |
