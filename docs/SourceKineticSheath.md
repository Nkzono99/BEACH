title: 境界から到達する粒子のシース応答

Lang: [日本語](SourceKineticSheath.md) | [English](SourceKineticSheath.en.md)

# 境界から到達する粒子のシース応答

matching-plane の外部密度を粒子軌道から求める場合は、`density_model="source_kinetic"` を選びます。
このページではモデルの条件、根の分類、診断 CSV を確認できます。時間依存の安定枝を選ぶ機能ではありません。

## 設定と互換性

既存の [matching-plane 設定](MatchingPlaneCoupling.html)の online backend に、次を指定します。
実行可能な PE なしの例は `examples/periodic2_matching_plane_source_kinetic.toml` です。
`beach examples/periodic2_matching_plane_source_kinetic.toml` は 2 バッチを実行し、
`outputs/periodic2_matching_plane_source_kinetic/` に summary と履歴を出力します。

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
density_model = "source_kinetic"
zhao_branch = "auto"
zhao_root_selection = "require_unique"
```

`density_model` の既定値 `zhao_legacy` は、公開済みの BEACH 実装の密度式・根選択を維持します。
新しい source-connected モデルを今後の主モデル候補として検証しますが、既存入力の物理モデルを暗黙に変更しません。
この互換動作に削除予定はありません。stationary Zhao と table backend は `density_model` を受理しません。

初期実装は、1D・無衝突・非磁化、ドリフトなしの ambient Maxwell 電子、half-Maxwell 放出電子、
反射しない冷たいイオンです。ambient electron の `drift_velocity=[0,0,0]`、ion の `temperature_ev=0` と
内向き drift が必要です。非零電子 drift に旧式の補正係数を掛ける代用はしません。
境界につながらない捕捉軌道は空です。ambient の外向き feedback は従来どおり transparent です。

`zhao_branch` は `auto`、`a`、`n`、`b`、`c`、`b0` を受理します。
全候補の探索・検証後に Type で絞り、検出根がちょうど一つのときだけ応答を使います。
複数根・未解決では応答を作れないため停止し、別 Type、従来密度式、有限境界へ切り替えません。
`minimum_energy` と `continuation` は新モデルでは未対応です。有限探索で一根を検出しても数学的な一意性は保証しません。

## Zhao 原著との対応

Zhao の電子モデル全体を Boltzmann 密度モデルと呼ぶのは不正確です。
[Zhao et al. (2020), §III.A、式 (3)–(6)](https://scholarworks.indianapolis.iu.edu/server/api/core/bitstreams/e20ede43-d66d-4b73-b89c-3d3192672188/content)
は Maxwellian 分布を速度領域ごとに積分し、到達・反射・帰還成分を区別しています。
零 drift の式 (3) は $n_{e,f}=(A_e/2)e^\phi\operatorname{erfc}(\sqrt{\phi-\phi_m})$ です。
Type B に対して経路の最低電位 $\phi_m=0$ を入れると、下記の source-connected 式になります。
これは原式の零 drift・Type B への適用であり、「Zhao 原著を Boltzmann モデルから置き換えた」という意味ではありません。

既存 BEACH の `evaluate_zhao_density_hat` の Type B は、上式の位置依存の速度下限を持たず、
零 drift では $(A_e/2)e^\phi$ になります。この差を source モードで修正しています。
A/C は既存式にも速度領域の cutoff があり、一律に Boltzmann 型ではありません。
`zhao_legacy` は既存 BEACH の再現用という名前で、原著への忠実性を示す名前ではありません。

[Zhao (2022) 学位論文](https://scholarsmine.mst.edu/doctoral_dissertations/3176/)の §4.3.1、式 (4.29) と
付録の式 (1) でも同じ速度下限付き積分を確認でき、最低電位は経路全体の最低電位と定義されています。2021 年 IEEE TPS 版
[DOI:10.1109/TPS.2021.3110946](https://doi.org/10.1109/TPS.2021.3110946) は本文を取得できていないため、
その版の Type B 固有の式を確認済みとはしません。非零 drift の正確な境界 VDF 輸送は別途検証が必要です。

## 規格化と密度

$\phi=e(\Phi-\Phi_{out})/(k_BT_e)$、$x=(z-H)/\lambda_{De}$、$E=-\phi'$、
$\tau=T_{ph}/T_e$、$G=\Gamma_{ph,H}^{out}/[n_i\sqrt{k_BT_e/m_e}]$ を使います。
旧 Zhao 内部の温度比は逆数です。指定電場は SI の $D_H/\epsilon_0$ を $T_e[\mathrm{eV}]/\lambda_{De}$ で割ります。
放出温度は整合面で測った平均法線 energy、$G$ は整合面を横切る外向き束です。

Type B の ambient 密度は、加速された入射粒子の速度下限を残した

$$
n_e(\phi)=\frac{A_e}{2}\operatorname{erfcx}(\sqrt\phi)
$$

です。$\phi_*$ を経路最低電位、$q=G\exp[-(\phi_H-\phi_*)/\tau]$ とすると、
A/N の極小の内側・外側では

$$
n_{e,in}=\frac{A_e}{2}e^{\phi_*}\operatorname{erfcx}(\sqrt{\phi-\phi_*}),\qquad
n_{e,out}=A_e e^\phi-n_{e,in},
$$
$$
n_{ph,out}=\sqrt{\frac\pi{2\tau}}q\operatorname{erfcx}(\sqrt{(\phi-\phi_*)/\tau}),\qquad
n_{ph,in}=2\sqrt{\frac\pi{2\tau}}G e^{(\phi-\phi_H)/\tau}-n_{ph,out}.
$$

B は内側の式で $\phi_*=0$、C は外側の式で $\phi_*=\phi_H$ とします。
イオンは $n_i=(1-2\phi/M^2)^{-1/2}$ で、全経路の $\max\phi<M^2/2$ を要求します。
イオン反射限界の clipping はしません。

$\phi(\infty)=E(\infty)=\rho(\infty)=0$ と指定 $E_H$ を課し、ambient 振幅 $A_e$ を中性条件から解きます。
固定された外部 VDF 振幅の問題とは異なります。電流 $J=M/\sqrt{m_i/m_e}+q-A_e e^{\phi_*}/\sqrt{2\pi}$ は出力です。
$J=0$ を追加条件にしません。準中性、第一積分、全経路で標本化した $E^2$ と電荷密度の極値、
極小曲率、イオン反射余裕、放出障壁、外端の到達可能性を検査します。

| Type | 電位 | 最低電位 |
|---|---|---|
| B | 正の整合面から 0 へ単調減少 | 0 |
| C | 負の整合面から 0 へ単調増加 | $\phi_H$ |
| A | $\phi_H>0>\phi_m$、内部極小あり | $\phi_m$ |
| N | $\phi_m<\phi_H\le0$、内部極小あり | $\phi_m$ |
| B0 | 全域 $\phi=E=0$ | 0 |

零電場では B0 に加えて非一様 B/C の境界極値も探索します。B0 の条件は
$A_e=2[1-\sqrt{\pi/(2\tau)}G]>0$ です。A/N の electron access と PE barrier は $\Phi_m$、B/C/B0 は従来の 0 V gauge です。

## 全候補の診断

source モデルの設定を `beach-zhao-atlas` に渡すと、旧 A/B/C 各一行形式に代わり、検出候補ごとに一行を書きます。
入力の 3 列は [既存 atlas](MatchingPlaneReference.html) と同じです。

```bash
beach-zhao-atlas source.toml queries.csv source-atlas.csv
```

`classification` は `solutions`、`analytically_excluded`、`search_unresolved`、`numerical_exception`、`invalid_input` です。
`numerical_failures` は候補評価中の非有限演算の回数です。正常な根を検出した場合もこの回数を残します。
未対応の drift や根選択 policy による運用上の除外は設定時にエラーとなり、解析的な不存在と混同しません。
`root_count` は採択した全 Type の根数で、`accepted=F` の棄却候補と理由も保存します。
根がない query も一行残し、電位を NaN にします。`phi_h_te` と `phi_min_te` は $T_e/e$ 単位で、
物理的最低電位を保存します。旧単調枝の placeholder とは異なります。
振幅、逃走束、流入束、電流、各残差、探索範囲、障壁の丸め指標を併記します。

探索は深さ $10^{-24}$ から $10^8$、基本 481 点の有限走査です。浅い領域と深い領域を追加し、
符号反転、標本化した極値、物理領域境界を刻み直します。`search_complete=T` は解析的除外にだけ付けます。
新モデルの $M\ge1,E_H>0$ では $E_H^2<2\sqrt{2\pi\tau}G$ が必要です。
Type B には $q\ge\sqrt{2/\pi}\tau/(1+\sqrt\tau)$ も必要ですが、この制約を他 Type の不存在証明には使いません。

`deep_root=T` は深さ $10^4$ 超、`barrier_resolution_limited=T` は
$\epsilon_{mach}(|\phi_H|+|\phi_m|)/(\phi_H-\phi_m)>10^{-7}$ を示します。
深い根を実用的な存在域の拡大と数えず、電位差に $T_e$ を掛けて非相対論近似も評価してください。
atlas は応答表ではありません。`beach-zhao-response` は選択後の完全直積だけを書き、欠損や曖昧点を補完しません。
生成表は `# density_model=...` を保持します。online の summary も密度モデルを記録します。

## 検証範囲

独立配布物 `lunar-sheath-solver-v0.1.0.zip` の 86 テストと、BEACH の Fortran 回帰を別に実行しました。
配布物の 1,211 条件では全 739 根の数・Type が一致し、電位・振幅・束・電流の尺度調整後の最大差は
$3.6\times10^{-13}$ でした。B/A/N の独立速度積分、A/N の Poisson 第一積分、B/N 共存、
B0 と非一様零電場 C、深さ打切り、SI 変換を回帰対象にしています。
再現手順と確認範囲は [検証例の README](../examples/source_kinetic_validation/README.md) にあります。

PE なしの BEACH 帯電では、source online と同じモデルから作った 129 行の table が、
整合面電位・変位について 0.2% 以内で一致しました。これはこの小規模 fixture の結合回帰であり、
PE を含む帯電の格子・時間刻み収束を保証しません。

共存例の各 $A_e$ を指定した既存 kinetic oracle の 14 条件では、3 条件が `steady`、
11 条件が `far_boundary_not_converged` でした。B の速度格子細分化と N の外部長延長で静的電位に近付きますが、
この試験だけでは時間依存の安定枝・半無限極限・格子収束は確定できません。未収束結果を応答表へ昇格させていません。

これは静的根探索の検証です。[時間依存 kinetic oracle](OuterKineticOracle.html) との比較には、各根の
$A_e$ に外部 electron source を合わせる必要があります。格子、速度範囲、冷たいイオン極限、外部長、時間の収束、
さらに BEACH の帯電・PE return/escape の収束は別の確認事項です。静的根の順序や小さい残差から時間安定性を判定しません。
