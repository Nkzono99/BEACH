title: 粒子更新

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
    ├─ open面が先 ──── escapeまたはscalar-barrier反射
    └─ reflect/periodic ─ 残り時間を再積分
```

| 工程 | 詳細 |
| --- | --- |
| 電場標本、Boris回転、位置更新 | [Boris粒子更新](BorisPusher.md) |
| 三角形衝突、box面、periodic image、発生順序 | [粒子の衝突・境界イベント](../developer/ParticleEvents.md) |
| 衝突・境界処理後の吸収と要素電荷差分 | [表面はどう帯電するか](SurfaceModels.md) |
| open面のescape、scalar-barrier反射 | [粒子源](../guide/ParticleSourcesBoundaries.md)、[粒子のescapeとreturn](ParticleEscapeReturn.md) |

## 1 stepで確定する状態

step の結果は、次時刻状態、表面吸収、box escape、または未完了 status のいずれかです。reflect、periodic、
scalar-barrier 反射は粒子を生存させ、境界処理後の残り時間を同じ更新法で進めます。

mesh hitではhit位置と要素indexを確定し、候補終点を採用しません。吸収粒子の電荷はその場でfieldへ入れず、
thread-localな要素電荷差分へ加え、batch末尾にcommitします。

## step数上限

粒子は吸収またはescapeするまで最大`particles.tracking.max_steps_per_particle`回進めます。上限に達しても分類されなかった粒子は
`survived_max_step`として記録し、暗黙に吸収やescapeへ数えません。

`particles.tracking.dt_s`は1回の粒子更新幅、`particles.tracking.max_steps_per_particle * particles.tracking.dt_s`は1粒子を追跡できる最大時間です。
一方、`run.batch.duration_s`は粒子供給量と表面電荷更新を結ぶ時間幅で、粒子step幅とは別です。

## 最初に確認するもの

- `particles.tracking.dt_s`を半分にして軌道、命中要素、吸収数が安定するか
- `survived_max_step`が結論に影響する割合になっていないか
- periodic seam、corner、reflect後の残り時間を跨ぐ軌道が意図どおりか
- 粒子の吸収・escape・未解決数と電荷収支が一致するか

設定値は[設定パラメータ](../reference/Parameters.md)、出力での分類は[出力ファイルを調べる](../guide/OutputGuide.md)から確認できます。

## Code reference

- step候補と衝突・境界イベントの順序: [`bem_particle_stepper.f90`](../../src/runtime/simulator/bem_particle_stepper.f90)
- 粒子を追跡するbatch loop: [`bem_simulator_loop.f90`](../../src/runtime/simulator/bem_simulator_loop.f90)
- step回帰テスト: [`test_particle_stepper.f90`](../../tests/fortran/test_particle_stepper.f90)
