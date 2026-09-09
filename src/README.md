# Fortran source layout

`src/` は fpm のライブラリソースです。モジュール名（`bem_*`）は互換維持のため変更せず、責務ごとに配置しています。

- `core/`: 基本型・定数
- `mesh/`: メッシュ初期化、テンプレート生成、外部メッシュ読込。`panel/` は三角形の幾何・面の向き・求積点と重み
- `physics/`: 物理モデルとその数値解法。電場評価、Boris pusher、衝突判定、ボックス境界、シースを含む
- `particles/`: 粒子 SoA 初期化、注入サンプリング
- `config/`: 設定型・TOML の基本値読み込み・領域別パーサ・正規化・実行前検証
- `runtime/`: 設定からの実行データ構築、リスタート入出力、シミュレータ本体、シースとの連成
- `tools/`: 応答表生成や分岐診断など、通常のシミュレーションとは別に使う Fortran ツール

`config/` は入力を `app_config` へ読み、派生値・参照を解決し、実行に必要な最小条件を確認します。
`bem_config_toml` は基本型・有限値・文字列長、`app_config_parser/` は table ごとの読取と領域別 preflight、
`bem_app_config_authoring*` は authoring 型と座標・配置の展開を担当します。
メッシュや粒子、境界電位を設定から構築する処理は `runtime/configuration/` に置きます。
詳細な列挙値・値域・組合せ診断は実行前の `beachx lint` が担当します。
開発・診断用の `beach --check-config beach.toml` は読取・正規化・残されたモデル成立条件だけを確認し、lint を代替しません。

電場の合成と評価は `physics/field_solver/` にまとめています。
`periodic/` は平面平均の zero mode と Fourier 評価、`panel/` は面電荷の Coulomb 積分、
`fmm/` は FMM の木・展開・周期遠方補正を担当します。zero mode は周期電場の平均成分を解く処理です。
snapshot がこれを非ゼロ成分と一様場へ加え、粒子が参照する場を作ります。

`mesh/panel/` の三角形求積は電荷や Coulomb kernel を扱いません。
幾何・離散化を mesh 側、幾何から物理量を求める積分を field solver 側に置くため、
メッシュ構築は電場評価に依存しません。

シースの式と解法は `physics/sheath/zhao/`、入出力の契約は `physics/sheath/` に置きます。
設定・MPI・バッチ状態との接続は `runtime/sheath/`、応答表の読み込みと補間はその `table/`、
オフライン生成器は `tools/sheath/` が担当します。物理モデルは設定型や MPI に依存しません。

公開モジュールは呼び出し口を保ち、粒子生成、出力、再開などの詳細を責務別の submodule に委譲します。
各入口と実装の対応は[開発者向けアーキテクチャ](../docs/Architecture.md)を参照してください。

注意: `*.i90` の生成物を `src/` 直下に置くと、fpm が重複ソースとして検出して同名モジュールを二重にコンパイルする場合があります。生成物は `src/` 外へ出力してください。
