title: 粒子更新

Lang: [日本語](ParticleTrackingCollision.md) | [English](ParticleTrackingCollision.en.md)

# 粒子更新

BEACHは粒子の位置と速度を同じ時刻で保持します。1 stepごとに次時刻の候補を作り、
meshとの衝突とbox境界の通過のうち、軌道上で最初に起こるものまでを確定します。

```text
(xⁿ, vⁿ)
    │ 予測中点で電場を1回評価
    ▼
Boris速度更新 + 台形位置更新
    │
    ▼
(xⁿ⁺¹, vⁿ⁺¹) の候補
    │ mesh hitとbox面通過を比較
    ├─ meshが先 ────── 吸収
    ├─ open面が先 ──── escapeまたはouter interfaceへ渡す
    └─ reflect/periodic ─ 残り時間を再積分
```

| 工程 | 詳細 |
| --- | --- |
| 電場標本、Boris回転、位置更新 | [Boris粒子更新](BorisPusher.html) |
| 三角形衝突、box面、periodic image、発生順序 | [粒子の衝突・境界イベント](ParticleEvents.html) |
| 衝突・境界処理後の吸収と要素電荷差分 | [表面電荷更新](SurfaceModels.html) |
| open面のescape、return、outer transfer | [粒子源](ParticleSourcesBoundaries.html)、[外部プラズマモデル](OuterPlasmaModels.html) |

## 1 stepで確定する状態

stepの結果は、通常の次時刻状態、表面吸収、box escape、outer interface crossing、または未完了statusの
いずれかです。reflectとperiodicは粒子を生存させ、境界処理後の残り時間を同じ更新法で進めます。

mesh hitではhit位置と要素indexを確定し、候補終点を採用しません。吸収粒子の電荷はその場でfieldへ入れず、
thread-localな要素電荷差分へ加え、batch末尾にcommitします。

## step数上限

粒子は吸収またはescapeするまで最大`sim.max_step`回進めます。上限に達しても分類されなかった粒子は
`survived_max_step`として記録し、暗黙に吸収やescapeへ数えません。

`sim.dt`は1回の粒子更新幅、`sim.max_step * sim.dt`は1粒子を追跡できる最大時間です。
一方、`batch_duration`は粒子供給量と表面電荷更新を結ぶ時間幅で、粒子step幅とは別です。

## 最初に確認するもの

- `sim.dt`を半分にして軌道、命中要素、吸収数が安定するか
- `survived_max_step`が結論に影響する割合になっていないか
- periodic seam、corner、reflect後の残り時間を跨ぐ軌道が意図どおりか
- 粒子の吸収・escape・未解決数と電荷収支が一致するか

設定値は[設定パラメータ](Parameters.html)、出力での分類は[出力の読み方](OutputGuide.html)から確認できます。

## Code reference

- step候補と衝突・境界イベントの順序: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- 粒子を追跡するbatch loop: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- step回帰テスト: [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)
