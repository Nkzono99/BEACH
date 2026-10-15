# ドキュメント構成の作り直し

作成日: 2026-10-15。基準: `main` の `6059366`。状態: 採用。段階 1〜3 (原稿のサイト掲載、原稿の差分の解消、入門ケースの作り直し) と、`docs/` のファイル配置の整理を実装済み。

公開サイト (GitHub Pages) のドキュメントは 43 ページ × 日英 2 言語、約 2 万行に増え、読みにくくなった。
一方、`docs/manuscript/` の技術原稿は物理モデルと数値手法を一つの流れで説明できている。
本書は、理論を原稿へ一本化し、サイトを「物理的に意味のあるケース」と「使い方」に絞る判断を記録する。

## 1. 現状の問題

| 問題 | 具体例 |
| --- | --- |
| 理論が細切れで原稿と二重管理 | 原稿 10 章とほぼ同じ内容を、サイトでは約 20 ページ (Algorithms、SurfaceModels、Photoelectron、Zhao×2、Periodic×2、BatchDuration×2、Direct/Treecode/FMM/FMMCore、Boris、ParticleEvents など) に分割している |
| 入口が重複 | `index` に「最初の 30 分」「目的から選ぶ」「ドキュメントの構成」の 3 つの読み順がある |
| 実装の語彙が利用者向け本文に出る | commit、receipt、ledger、teacher、zero mode などが説明なしに出る。日英混在の文が多い |
| ケース設計が設定キーの順 | `ConfigurationRecipes` はチュートリアルを複製し、形状 → 表面 → 粒子源 → box → solver → 時間の順にキーを変える。環境とスケールから数値条件を決める流れがない |
| 入門ケースに物理的な意味がない | `run.batch.duration_s=0`、`weight=2e5` で実時間を持たず、本文で「研究ケースではない」と繰り返し断っている |
| examples が機能の組合せ表 | `zero_current_no_photo_secondary`、`panel_fmm`、`cached_periodic` など。シナリオとして読めない |
| 原稿とサイト・SPEC の不一致 | 原稿 7 章は浮遊導体を「2 軸周期場とは併用しない」と書くが、SPEC と ADR 0002 は periodic2 でも使えるとする |

## 2. 判断

| 項目 | 判断 |
| --- | --- |
| 正本の言語 | 日本語だけを正本にする。英語ページ (`*.en.md`) は削除し、構成が固まってから翻訳し直す。旧英語版は git 履歴に残る |
| 理論の正本 | `docs/manuscript/` (LaTeX) を唯一の正本にする。サイトには原稿から生成した章ページと PDF を載せ、手書きの理論ページは廃止する |
| サイトの役割 | 原稿が扱わない「使い方」(原稿の理論をどの設定で使い、どの出力で確かめるか) とケース集に絞る |
| ケースの範囲 | ADR 0002 の製品構成 (periodic2 照射セル) を中心にする。free-space は入門ケースだけ |
| 入門ケース | 実時間と電流密度を持つ「5 eV の電子ビームで浮遊導体板を帯電させる」に作り直し、飽和電位 ≈ −E_b/e、容量、帯電時間を手計算の見積もりと照合する。最後に絶縁体へ替えて違いを見る (8 節) |
| 開発者向けページ | サイト末尾の折りたたみ節へ移す。数値実装の詳細ページはここに集める |

## 3. 新しい構成

| 節 | ページ | 内容 |
| --- | --- | --- |
| はじめる | BEACH とは | 何を解き何を解かないか、計算サイクルの図 1 枚、読む順 1 つ |
| | インストール | 現行のまま |
| | 入門: 電子ビームで絶縁体を帯電させる | 物理条件付きの 10 分チュートリアル |
| ケース集 | 日照下の平坦な月面 | 太陽風 + 光電子の periodic2 セル。外部シース接続ありで、表面電位を 1-D シース解と照合する |
| | 凹凸のある表面の不均一帯電 | closepack 形状、太陽高度を変えて影と帯電の関係を見る |
| | 表面粒子に働く力と離脱 | 帯電後の力、離脱経路 |
| ケースを設計する | 設計の流れ | 環境 (n、T、v、光電流) → 特性スケール (λ_D、帯電時定数、セル寸法との比) → 形状 → 上端の閉じ方 → 数値条件 → 検証 |
| | 形状と表面モデル | mesh template、OBJ、insulator / conductor |
| | 太陽風・光電子・二次電子を入れる | 境界流入、光電子、二次電子、`volume_seed` |
| | 上端の閉じ方と外部シース | 開放、閉じた光電子、スカラー補正、外部シース接続 |
| | 場の計算 | periodic2、FMM、cache、Direct との比較 |
| | 時間幅と粒子数 | `dt_s`、`batch.duration_s`、adaptive、粒子数と重み |
| | 設定ファイルの作成と検証 | 現行の Configuration |
| 実行して調べる | 実行・HPC・再開 / 出力を読む / 可視化 / 結果を検証する / 問題を診断する | 現行ページを整理。検証の理論は原稿 9 章へ |
| 理論 (原稿) | 目次と PDF、第 1〜10 章 | `docs/manuscript/` から生成 |
| リファレンス | 入力パラメータ / 出力形式 / 用語集 / Python API / Fortran API | 現行のまま |
| 開発者向け (折りたたみ) | アーキテクチャ / 開発とテスト / 物理リリース検証 / 数値実装の詳細 / Fortran 依存関係 | |

## 4. 現行ページの行き先

「設定・出力」は利用者向けの使い方ページへ、「理論」は原稿の該当章へ移す。原稿にない理論は、
削除する前に原稿へ追記する (5 節)。

| 現行ページ | 行き先 |
| --- | --- |
| `index` | 書き直し |
| `Installation` | 据え置き |
| `Tutorial` | 入門ケースとして作り直し |
| `ConfigurationRecipes` | 「設計の流れ」として書き直し。§2–3 は形状ページ、§7 は時間幅ページへ |
| `ParticleSourcesBoundaries`、`ReservoirInjection` | 「太陽風・光電子・二次電子を入れる」へ統合 |
| `PeriodicPlasmaSurface` | 共通構成は「日照下の平坦な月面」、方式の選択は「上端の閉じ方」へ |
| `FieldSolvers`、`FMM` | 「場の計算」へ統合 |
| `Configuration`、`Execution`、`OutputGuide`、`PostprocessTutorial`、`Troubleshooting` | 据え置き。用語と相互リンクを整理 |
| `ObjectForcesDetachment` | 「表面粒子に働く力と離脱」ケースへ |
| `ValidationGuide` | 手順は残し、理論は原稿 9 章へ |
| `Algorithms` | 廃止。図 1 枚を `index` へ、本文は原稿 2 章 |
| `SurfaceModels` | 廃止。設定と出力は形状ページへ、理論は原稿 7 章 |
| `PhotoelectronEmission` | 廃止。設定と出力は粒子ページへ、理論は原稿 6・7 章 |
| `SecondaryElectronEmission` | 廃止。設定と出力は粒子ページへ。理論は原稿に追記 (5 節) |
| `ZhaoStationaryClosure`、`ZhaoUpstreamBand` | 廃止。設定と確認は上端ページへ、理論は原稿 8 章 |
| `ParticleEscapeReturn` | 廃止。上端ページへ、理論は原稿 5 章 |
| `PeriodicElectrostatics` | 廃止。設定と収束確認は場のページへ、理論は原稿 4 章 |
| `BatchDurationStability` | 時間幅ページへ統合 |
| `BatchDurationTheory` | 廃止。原稿 7 章 |
| `BorisPusher`、`ParticleTrackingCollision` | 廃止。`dt_s` の選び方は時間幅ページへ、理論は原稿 5 章 |
| `DirectSolver`、`Treecode`、`ParticleEvents`、`SurfaceChargeNumerics`、`PeriodicFarCorrection`、`FMMCore` | 開発者向け「数値実装の詳細」へ。Code reference を中心に残す |
| `Parameters`、`OutputReference`、`Glossary`、`PythonPostprocessAPI` | 据え置き |
| `Architecture`、`Workflow`、`PhysicsReleaseVerification`、`FortranDependencyMap` | 開発者向け節へ |
| `agent-user-guide` | 新しい構成に合わせてルーティングを更新 |

## 5. 原稿側で先に直すこと

ページを廃止する前に、そのページの理論が原稿に含まれることを確認する。2026-10-15 に廃止予定の 11 ページを原稿と
突き合わせ、原稿になかった理論を次のように取り込んだ (段階 2)。設定キー、出力列、確認手順は原稿へ入れず、
段階 4 の使い方ページへ移す。

| 廃止するページ | 原稿の該当箇所 | 段階 2 で原稿へ追加したもの |
| --- | --- | --- |
| `SecondaryElectronEmission` | 6.5 節 (新設) | 収率の式、生成数の標本化、重みと速度、反作用電荷、世代上限、乱数の決め方、適用範囲。文献 Marchand2014 を追加 |
| `PhotoelectronEmission` | 6.4 節、7.2 節 | 周期軸での ray の回り込み、放出の法線と速度分布の詳細、閉じた光電子で残る鉛直双極子、固定電流の統計の注意 |
| `SurfaceModels` | 7.1 節、7.3 節 | 浮遊導体を場の評価法から作る影響行列で定式化し、periodic2 でも使えるよう修正 (旧版は「2 軸周期場とは併用しない」で実装と不一致)。LU 分解の再利用、外部シースの基準との関係 |
| `BatchDurationTheory` | 7.4 節 | 物理的な時間尺度 ($\omega_{pe}^{-1}$、帯電時間) と、Monte Carlo 誤差との結合 |
| `BorisPusher`、`ParticleTrackingCollision` | 5.1〜5.3 節 | 時間刻みの選び方、event の速さの性質と退化 event の停止 |
| `ParticleEscapeReturn` | 5.3 節、8.8.4 節 | 再配置反射の辺と角での扱い、二次電子の外部障壁 |
| `ZhaoStationaryClosure` | 8.8.4〜8.8.5 節 | 二次電子を含む外部根の更新と暫定の配分、零モードを分ける周期場が必要な理由と試験計算 |
| `ZhaoUpstreamBand` | 8.5.2 節、8.8.3 節 | 帯を埋める物理 (旋回混合、揺らぎ、不安定性、有限の天体) の見積もり |
| `PeriodicElectrostatics` | 4.4〜4.6 節 | 追加なし (既に含まれていた) |
| `Algorithms` | 2 章 | 追加なし (既に含まれていた) |

## 6. 原稿をサイトへ載せる方法

`tools/manuscript_to_starlight.py` を追加し、`make docs-pages` の前段で実行する。

1. 既存の `docs/manuscript/Makefile` で PDF を作る (XeLaTeX + BibTeX)。生成された `.aux` と `.bbl` を以降の番号付けに使う。
2. 章ごとに LaTeX を前処理する。`\ref` / `\eqref` を `.aux` の番号付きリンクへ、`\cite` を `.bbl` の番号付きリンクへ置き換え、
   ラベル付きの式に `\tag{}` を付け、`\input{figures/...}` を SVG 画像と図番号付きキャプションへ置き換える。
3. pandoc 3 (`pypandoc-binary` で同梱版を固定) で Markdown に変換し、Starlight の content へ書き出す。
4. TikZ 図は図ごとに単独でコンパイルし、`pdftocairo -svg` で SVG にする。
5. PDF を `public/manuscript/beach_models.pdf` として配置し、目次ページからダウンロードできるようにする。

CI には TeX Live (xetex、CJK)、Noto CJK、TeX Gyre、poppler-utils を追加する。生成物は commit しない。

## 7. 進め方

| 段階 | 内容 | 完了条件 |
| --- | --- | --- |
| 1 | 原稿をサイトへ載せる (6 節)。英語 locale を外す | 第 1〜10 章が式番号・図・文献つきで表示され、PDF を取得できる |
| 2 | 原稿の差分を直す (5 節) | 廃止予定ページの理論がすべて原稿にある |
| 3 | 入門ケースを作り直す (実装済み) | 計算ノードで実行し、飽和電位を見積もりと照合した図と数値を載せる |
| 4 | 「ケースを設計する」節を新設し、旧ページを統合・廃止する | 4 節の表のとおりにナビゲーションを置き換え、リンク切れがない |
| 5 | ケース集 (平坦な月面、凹凸面、離脱) | 各ケースを実行し、解析解または 1-D シース解との照合を載せる |
| 6 | `index` と `agent-user-guide` を書き直す | |
| 後日 | 英語版 | 日本語版の構成が固まってから翻訳する |

## 8. 入門ケースの選び方 (段階 3)

当初は絶縁体板を予定したが、絶縁体では電荷が当たった場所に残るため、板の電位が一様にならない。中央に溜まった電荷が
電子を外へそらし、端が $-E_b/e$ に届く前に吸収が止まるので、「飽和電位 ≈ $-E_b/e$」という一つの数で照合できない。
浮遊導体板なら板全体が等電位になり、次の 3 点を手計算と比べられる。

| 見積もり | 計算 (`examples/tutorial_electron_beam.toml`) |
| --- | --- |
| 板と放出面の電位差が $-E_b/e=-5$ V で吸収が止まる | $-5.38$ V で吸収が 0 (エネルギー幅 0.32 eV の分だけ深い) |
| 容量 = 同面積の円板 $8\epsilon_0r$ = 0.80 pF | 0.796 pF |
| 帯電時間 約 0.5 μs 以上 | 約 0.85 μs で $-5$ V に達する |

絶縁体は、同じビームで電位の分布と吸収の止まり方が変わる例として、チュートリアルの後半で扱う。
ビームの電子密度は約 $10^5\ \mathrm{m^{-3}}$ で、空間電荷による電位は $10^{-7}$ V 程度なので、BEACH のテスト粒子近似が成り立つ。
結果は OpenMP の thread 数によらず一致する (1 thread で約 8 分、8 thread で約 1 分)。

## 9. `docs/` のファイル配置

`docs/` 直下に並んでいたページを、役割ごとのディレクトリへ移した。ページ間のリンクはファイルからの相対パス
(`../reference/Parameters.md#run`) で書き、GitHub 上でそのまま開けるようにした。サイト生成
(`tools/sync_starlight_docs.py`) がリンク先を解決してサイトの URL に変換する。

| ディレクトリ | 内容 |
| --- | --- |
| `docs/` 直下 | `index.md`、`agent-user-guide.md` |
| `start/` | インストール、入門チュートリアル |
| `guide/` | ケースを作る・実行する・結果を調べる手順 |
| `models/` | モデルの説明。理論は原稿へ移したので、段階 4 で使い方を `guide/` へ統合して廃止する |
| `reference/` | 入力パラメータ、出力形式、用語集、Python API |
| `developer/` | アーキテクチャ、開発とテスト、数値実装の詳細、生成した依存関係マップ |
| `manuscript/` | 技術原稿 (理論の正本) |
| `adr/` | 設計判断の記録 |
| `media/`、`images/` | 図 |

段階 5 のケース集は `cases/` に置く。

