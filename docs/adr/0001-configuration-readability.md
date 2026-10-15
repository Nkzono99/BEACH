# パラメータ設計の見直し

状態: 置き換え済み。入力形式の改訂は 2.0 の 7 グループ形式 (`a2e521d`) で実施した。

対象は、現在の `beach.toml`、Python の設定検査、Fortran の読み込み処理です。
設定の判断順を「物理条件 → 数値精度 → 出力」に揃えます。既存形式で解消できる混乱は
利用ガイドとリファレンスへ反映し、形式の変更が必要なものは以下の設計課題として区別します。

## 今回の判断

公開キー、既定値、物理モデルは維持します。名前の変更だけでは粒子数・重み・時間幅の関係は
簡単にならず、既存ケースの移行も必要になります。新規ケースでは、使う機能に必要な項目だけを
記述し、全パラメータを並べた設定を出発点にしません。

| 現状の混乱 | 判断と文書への反映 |
| --- | --- |
| `sim.dt` と `batch_duration` の区別が難しい | 軌道刻みと電荷更新幅を別々に選ぶ。まず秒単位の `batch_duration` を使い、`batch_duration_step` は両者の比を指定したい場合に使う |
| `npcls_per_step` は実際には batch 当たりの個数 | ガイドで時間単位を明示する。粒子数を倍にするとき、固定した投入量の比較では重みを半分にする |
| `[sim]` に時間、solver、診断、例外処理が集中している | リファレンスを用途別に分け、通常の時間設定の途中へ診断や破棄設定を挟まない |
| `tol_rel` が収束停止条件に見える | 監視値の項へ分ける。収束確認は観測量・時間幅・標本数の比較で行う |
| 流束源の密度と標本数が同列に並ぶ | 密度・温度・流束を先に決める。`target_macro_particles_per_batch` と `w_particle` は、その物理供給量の標本化方法として一方を選ぶ |
| 刻みの比較で追跡時間や終了時刻も変わる | `dt` と `max_step`、`batch_duration` と `batch_count` の組合せを明記する |

操作は[ケース設計](../ConfigurationRecipes.md)、全キーの契約は[入力リファレンス](../Parameters.md)、
結果の採用条件は[妥当性確認](../ValidationGuide.md)に置きます。このメモを読まなくてもケースを作れる構成にします。

## 入力形式を改訂する場合の候補

以下は **未実装の設計案** です。現在の TOML に新しいセクションや別名を追加するものではありません。
一つの判断を一つの場所で表現し、同じ値への入口を増やさない方針とします。

| 課題 | 改訂案 | 維持すべき意味 |
| --- | --- | --- |
| `npcls_per_step` の時間単位 | 次の形式改訂で「batch 当たり」を名前に含める | `volume_seed` の個数と、流束源の標本数 target は異なるため統合しない |
| 大きな `[sim]` | 時間進行、場の評価、境界イベント制御を責務別の table に分ける | 同じ物理量を複数の table から上書きさせない |
| 境界流入のみでも `volume_seed` と `npcls_per_step=0` が必要 | 領域内の生成と境界流入を独立に有効化する構造にする | 流入は生成、`particle_boundary` は外向きの作用。この二つは統合しない |
| `target_macro_particles_per_batch=-1` は species 1 への暗黙参照 | 重みの共有先を安定した `species_key` で指定する | 並べ替えで共有先が変わらないこと、循環参照を許さないこと |
| `coupling_atol` は異なる単位を持つ 4 成分配列 | PE 流束、PE energy、electron 流束、ion 流束を名前付きで表す | 各単位と、inactive 成分を判定から除外する規則を保つ |
| `multiple_box_events_soft_discard_abs_charge_limit` は警告閾値 | 改訂時に警告を表す名前へ揃える | 最初の閾値超過で警告し、絶対電荷だけでは停止しない |

## 形式改訂の検証範囲

実装時は、現在の設定を読む経路、移行後の設定、両者が生成する物理条件を一緒に確認します。
単なる名前変更なら、同じ乱数条件で注入量・電荷収支・結果が一致することが必要です。
互換入力を残す場合は、対象リリースと廃止条件を明記します。

変更対象は Python と Fortran の parser、配布用を含む JSON Schema、設定例、CLI の生成設定、
リファレンス、必要に応じて model fingerprint と restart の契約です。
schema だけ、または Python の lint だけが新形式を受け付ける状態では完了としません。

実装を調べる入口:

- [Python 設定処理](../../beach/config/core.py)
- [Fortran parser](../../src/config/app_config_parser/bem_app_config_parser.f90)
- [Fortran 入力の組合せ検査](../../src/config/app_config_parser/bem_app_config_parser_finalize.f90)
- [流束とマクロ粒子数](../../src/particles/injection/bem_injection_flux.f90)
