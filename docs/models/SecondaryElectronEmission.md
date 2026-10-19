title: 二次電子放出

# 二次電子放出

二次電子放出（SEE）は、表面へ吸収された粒子のエネルギーと入射角から電子を生成するモデルです。
入射種に放出条件、放出先の種に温度を指定すると、BEACH は二次電子を同じバッチの凍結場で追跡します。
表面の反作用電荷と、帰還・脱出を含む収支を確認できます。

## 何を解くか

入射粒子の質量を $m$、速度を $\mathbf v$ とし、エネルギー $E=m|\mathbf v|^2/(2e)$ を eV で求めます。
入射側を向く単位法線を $\hat{\mathbf n}$ とし、$c=\max(0.1,-\mathbf v\cdot\hat{\mathbf n}/|\mathbf v|)$ を使います。
`whipple` の収率は

$$
\delta = \frac{1.114\delta_{max}}{c}\left(\frac{E_{max}}{E}\right)^{0.35}
\left[1-\exp\left(-2.28c\left(\frac{E}{E_{max}}\right)^{1.35}\right)\right]
$$

です（Whipple 1981; Marchand et al. 2014, 式 (2)）。$E\le0$ では収率を 0 にします。
`constant` は入射エネルギー・角度に依存せず `const_yield` を使い、イオンにも指定できます。

## 仮定

放出マクロ粒子数は $\delta|q_{in}/q_{see}|$ の整数部と、端数を確率とする Bernoulli の和です。
各二次電子は入射粒子の重み $w_{in}$ を受け継ぎます。速度は光電子と共有するサンプラーで生成し、
法線成分は束で重み付けした半 Maxwell 分布、接線成分は平均 0 の Gauss 分布（$6\sigma$ 打ち切り）です。
法線方向の drift は 0、温度は放出先の種の設定を使います。

## BEACH での扱い

入射粒子を吸収して $q_{in}w_{in}$ を衝突要素へ置いた後、二次電子を衝突点から法線方向へ $10^{-12}$ m ずらして生成します。
周期座標はセル内へ戻します。二次電子 1 個ごとに同じ要素へ $-q_{see}w_{in}$ を必ず置き、
放出電荷と反作用電荷を相殺します。光電子の放出電荷チャネルを共有します。

生成粒子は世代ごとに追加して同じ凍結場で追跡し、吸収、脱出、step 上限の判定へ含めます。
二次電子自身の再吸収にも放出条件を指定できます。一次粒子を世代 0 とし、設定した世代上限に達した粒子は
さらに放出せず、放出するはずだった個数を診断へ記録します。その分の反作用電荷は置きません。
乱数は seed とバッチ・rank・親粒子・step・世代・放出順から決まり、thread 数や実行順序には依存しません。

外部シースを使う場合は SEE 種を `sheath.species.secondary_electron` に指定します。SEE は光電子と同じ
上端障壁（Type A は $\phi_m$、B/C は 0 V）と、周期セルの面内一様帰還を使います。
外部根に専用の SEE 項は加えず、PE と SEE の上端流出スペクトルを合算します。
規格化と固定電流の配分は [外部シースの定常解](ZhaoStationaryClosure.md#外部根を更新する)を参照してください。
初期根は PE の表面半 Maxwell 源だけで解き、PE なしなら源 0 の Type C から始めます。SEE は最初の更新から入ります。

## 設定

世代上限は `particles.secondary_emission.max_generations`（既定 `3`、0 以上）で指定し、一次粒子を世代 0 とします。

既存の入射種へ、次の差分を加えます。

```toml
[particles.species.surface_interaction.secondary_emission]
enabled = true
model = "whipple"
yield_max = 1.5
energy_max_ev = 500.0
target_species = "secondary_electron"
```

放出先は独立した粒子源を持たない種にします。

```toml
[[particles.species]]
species_key = "secondary_electron"
charge_c = -1.602176634e-19
mass_kg = 9.1093837139e-31

[particles.species.source]
mode = "secondary_emission"

[particles.species.distribution]
temperature_ev = 2.0

[particles.species.boundary]
z_high = "open"
```

`sampling.*`、`inflow.*`、放出電流密度、光線数、`charging.closure` は省略します。
全キーの既定値・制約は [入力パラメータ](../reference/Parameters.md)を参照してください。
完全な設定は [`zero_current_secondary.toml`](../../examples/zero_current_secondary.toml) と
夜側の [`zero_current_no_photo_secondary.toml`](../../examples/zero_current_no_photo_secondary.toml) にあります。

## 確認する出力

`secondary_emission_history.csv` は受理バッチごとに入射種の命中数、放出個数・電荷、平均収率、世代上限で抑えた個数を記録します。
`summary.txt` は累積値を、`charge_ledger.csv` は SEE 種の放出・再吸収・脱出を記録します。
外部根を更新する場合は `matching_plane_history.csv` の `secondary_share` と
`secondary_budget_residual_a_m2` を確認します。列の単位と定義は [出力リファレンス](../reference/OutputReference.md)にあります。

## 適用範囲

材質や `mesh_id` によるパラメータの違い、弾性・非弾性反射、表面導電は含みません。
入射角の下限 0.1 は暫定の発散回避であり、斜入射での材質応答を保証しません。
外部根の脱出束は PE と SEE の合算に対する値です。PE があれば SEE は完全に追跡し、各バッチで追跡した SEE の脱出電流を
差し引いた残りを PE の脱出・帰還目標にします。PE がなければ SEE の帰還・脱出を根の帰還束・脱出束に固定し、放出の
反作用電荷を根に与えた外側源（追跡した H 外向き束の窓平均）に揃えます。放出量そのものは追跡で決まり、根は予測しません。
表面収支はどちらも毎バッチ閉じ、`secondary_budget_residual_a_m2` は閉じられなかった分（PE あり: SEE の脱出が根の
脱出束を超えた分、PE なし: SEE が放出されなかったバッチの根の脱出分）だけを記録します。
SEE の放出個数は追跡した衝突から決まり、入射種の `fixed_current` の倍率は掛かりません。
`fixed_current` / `neutral_return` を SEE 種へ指定すると設定エラーです。
粒子はバッチ内で完結し、次のバッチへ飛行状態を持ち越しません。

負荷見積もりは独立した注入源だけを数え、衝突時の SEE 生成を除外します。
SEE 有効時は `secondary_emission_scope=excluded_collision_generated_particles` と表示します。
