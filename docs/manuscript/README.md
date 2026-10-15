# BEACH の物理モデルと数値アルゴリズム

Lang: [日本語](README.md) | [English](README.en.md)

BEACH の物理モデル、数値手法、先行研究をまとめた、日本語の論文・技術書形式の LaTeX 初稿です。
[beach_models.tex](beach_models.tex) が入口です。章別の本文は `sections/`、文献は
[references.bib](references.bib) にあります。日本語の概要と英語の abstract、目次、数式、
20 点の図、図目次、参考文献を含みます。

## 内容

1. 表面帯電、月面・無大気天体、先行研究
2. モデル分担、支配方程式、バッチ計算
3. P0 三角形要素、解析積分、自己項・片側電場
4. Direct、treecode、Cartesian FMM、2 軸周期場と Ewald・零モード
5. 同時刻 Boris 更新、衝突、箱境界・周期画像
6. 境界流入、速度分布、照射 ray と光電子放出、吸収時の二次電子放出
7. 電荷収支、帰還 closure、浮遊導体、バッチ幅の安定性
8. 外部シースのモデルと上端面への接続（仮定、軌道と粒子集団、密度・流束の導出、Poisson 方程式の第一積分と Type A、B、C、閉じ方、物理解の条件、数値アルゴリズム、数値例、外部根の条件と型の選択、電子 drift と場の下限、上端面への写像、外部根の更新）
9. 検証、収束、Monte Carlo 誤差、再現性
10. 力・離脱の評価、論文に記載する方法と限界、原稿作成に関する記載

## 説明図

計算サイクル図と、次の 19 点を本文の該当箇所に含みます。
図は TikZ で描いたベクトル図で、拡大しても文字と線が鮮明に表示されます。
作図ソースは `figures/`、色と線の共通設定は [figures/styles.tex](figures/styles.tex) です。
図を編集して `make` を実行すると PDF に反映されます。追加の作図ソフトや shell escape は不要です。

| 章 | 図の内容 | 作図ソース |
| --- | --- | --- |
| 2 | 計算領域とモデル分担 | [model_domain.tex](figures/model_domain.tex) |
| 3 | 三角形の評価幾何と片側電場 | [panel_geometry.tex](figures/panel_geometry.tex) |
| 4 | FMM の集約・変換・共有 | [fmm_translations.tex](figures/fmm_translations.tex) |
| 4 | 周期画像と零モードの分担 | [periodic_split.tex](figures/periodic_split.tex) |
| 5 | Boris 更新と最初の幾何イベント | [particle_step.tex](figures/particle_step.tex) |
| 6 | 体積分布と境界流束分布 | [flux_distribution.tex](figures/flux_distribution.tex) |
| 6 | 照射による遮蔽と光電子の電荷収支 | [photoemission.tex](figures/photoemission.tex) |
| 7 | バッチ幅による線形緩和の安定性 | [batch_stability.tex](figures/batch_stability.tex) |
| 8 | 電子のエネルギー図と軌道の分類（数値解） | [sheath_model_populations.tex](figures/sheath_model_populations.tex) |
| 8 | 電荷密度と第一積分（数値解） | [sheath_model_first_integral.tex](figures/sheath_model_first_integral.tex) |
| 8 | 外部シースの Type A、B、C | [sheath_types.tex](figures/sheath_types.tex) |
| 8 | 上流の帯と場の不足（数値解） | [sheath_model_upstream_band.tex](figures/sheath_model_upstream_band.tex) |
| 8 | 零電流根の求解手順 | [sheath_model_algorithm.tex](figures/sheath_model_algorithm.tex) |
| 8 | Type A、B、C の電位と電場（数値解） | [sheath_model_profiles.tex](figures/sheath_model_profiles.tex) |
| 8 | 粒子集団ごとの密度（数値解） | [sheath_model_densities.tex](figures/sheath_model_densities.tex) |
| 8 | 電子の速度分布（数値解） | [sheath_model_vdf.tex](figures/sheath_model_vdf.tex) |
| 8 | 太陽高度に対する零電流根（数値解） | [sheath_model_alpha_sweep.tex](figures/sheath_model_alpha_sweep.tex) |
| 8 | 物理解が見つかる型の分布（数値解） | [sheath_model_type_map.tex](figures/sheath_model_type_map.tex) |
| 8 | 外部シースと上端面の接続 | [sheath_mapping.tex](figures/sheath_mapping.tex) |

模式図の形状・軌道・電位曲線は数値解ではありません。流束分布とバッチ幅の図は、本文の式から
導いた解析例です。表で「数値解」とした第 8 章の図は、外部シースモデルを解いた結果です。
各キャプションに区別を示し、根拠となる節・式、元データ、図の適用限界を
[figures/figure_trace.yaml](figures/figure_trace.yaml) に記録しています。

### 第 8 章の数値解の図を作り直す

数値解の図は、計算と作図の 2 段階で作ります。生成した TikZ ソースと元データ（`figures/data/sheath_model/`）は
リポジトリに含めるので、`make` だけなら Python は不要です。計算条件や外部シースモデルを変えたときだけ、
次を計算ノードで実行します。計算には sheath-model の Python パッケージ（NumPy、SciPy）が必要です。

```bash
cd /path/to/BEACH/docs/manuscript
# 数値解を CSV/JSON へ（32 並列で約 9 分。大半は型の分布の格子）
PYTHONPATH=/path/to/sheath-model python tools/sheath_model_data.py figures/data/sheath_model --workers 32
# CSV から figures/sheath_model_*.tex を生成（NumPy だけを使う）
python tools/sheath_model_figures.py figures/data/sheath_model figures
make
```

`provenance.json` に、使った sheath-model の commit と共通の物理条件を記録します。
`sheath_model_algorithm.tex` は手で描いた図で、生成の対象ではありません。

## 書き方の方針

本文は、プラズマ物理と数値計算の一般的な知識を前提に、モデルとアルゴリズムを原稿の中で完結させます。
入力キー、ソースのパス、出力ファイル名、テスト名、commit、計算機環境は本文に書きません。
手法の由来は文献で示します。実装との対応は次の表で管理します。

| 章 | 主な実装・文書 |
| --- | --- |
| 2 計算サイクル | `src/runtime/simulator/bem_simulator_loop.f90`、`docs/models/Algorithms.md` |
| 3 三角形場 | `src/physics/field_solver/panel/bem_panel_kernel.f90`、`src/physics/field_solver/panel/bem_panel_self_terms.f90`、`docs/developer/DirectSolver.md` |
| 4 treecode・FMM | `docs/developer/Treecode.md`、`src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`、`docs/developer/FMMCore.md` |
| 4 周期遠方場 | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`、`docs/developer/PeriodicFarCorrection.md` |
| 4 零モード | `src/physics/field_solver/periodic/bem_periodic_zero_mode_plan.f90`、`src/physics/field_solver/periodic/bem_periodic_zero_mode_eval.f90`、`docs/models/PeriodicElectrostatics.md` |
| 5 粒子積分・衝突 | `src/physics/bem_pusher.f90`、`src/physics/bem_collision.f90`、`src/runtime/simulator/bem_particle_stepper.f90`、`docs/models/BorisPusher.md`、`docs/developer/ParticleEvents.md` |
| 6 流入・電子放出 | `src/particles/injection/bem_injection_flux.f90`、`src/particles/injection/bem_photoelectron_injection.f90`、`src/particles/bem_secondary_emission.f90`、`src/runtime/simulator/bem_simulator_particles.f90`、`docs/guide/ParticleSourcesBoundaries.md`、`docs/models/PhotoelectronEmission.md` |
| 7 電荷更新・導体 | `src/physics/bem_surface_models_conductor.f90`、`docs/developer/SurfaceChargeNumerics.md`、`docs/models/BatchDurationTheory.md` |
| 8 外部シースモデル | sheath-model（`fpm.toml` で commit を固定）の `docs/kinetic-model.md`、`docs/algorithm.md`、`docs/spectral-api.md`、`src/internal/sheath_model_admissibility.f90`、`src/sheath_model_photoelectrons.f90`、`sheath_model/_physics.py`、`sheath_model/_profile.py`、`sheath_model/search.py` |
| 8 上端面への接続 | `src/runtime/sheath/bem_surface_current_model.f90`、`src/runtime/sheath/bem_zhao_outflow_refresh.f90`、`docs/models/ZhaoStationaryClosure.md`、`docs/models/ZhaoUpstreamBand.md` |
| 9 検証 | `docs/guide/ValidationGuide.md`、`docs/developer/PhysicsReleaseVerification.md`、`docs/guide/BatchDurationStability.md` |
| 10 力・離脱 | `docs/guide/ObjectForcesDetachment.md`、`docs/reference/PythonPostprocessAPI.md` |

## PDF の作成

XeLaTeX、BibTeX、`xeCJK`、標準的な LaTeX パッケージ、TeX Gyre と Noto CJK フォントが必要です。
本文は Noto Serif CJK JP、見出しは Noto Sans CJK JP、コードは Noto Sans Mono CJK JP を使います。
別のフォントを使う場合は `beach_models.tex` のフォント指定を変更してください。
TeX の shell escape は不要です。

ローカル PC または計算ノードの割当内で、次を実行します。

```bash
cd /path/to/BEACH/docs/manuscript
make
```

生成物は `build/beach_models.pdf` です。補助ファイルも `build/` にまとめ、Git の追跡対象にしません。
LaTeX → BibTeX → LaTeX 2 回の順で引用と相互参照を解決します。

KUDPC のログインノードでは組版を直接実行せず、利用可能な Sys module、`spartition`、`qgroup` を確認して、
計算ノードへ投入します。確認済みのキューに置き換えてください。

```bash
tssrun -p <queue> -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /path/to/BEACH/docs/manuscript && make'
```

### この環境での組版確認

2026 年 10 月 15 日、System C の計算ノードで 20 点の図を含む 50 ページの A4 PDF を生成しました。
未解決の引用・参照、行のはみ出し、BibTeX の警告はありません。
既存の TeX Live 2018 に不足していた `xeCJK` と LaTeX3 バンドルは、同世代の TeX Live アーカイブから
`BEACH/build/manuscript-tools/texmf2018/` に一時展開して使いました。既存の TeX インストールは変更していません。
この一時領域は Git の追跡対象外です。通常の TeX 環境では必要パッケージを導入して上の `make` を使います。

ログインノードの Sys module と異なる system へ投入するときは、Sys module を切り替えてキューを確認します。
次は現在の module が SysG の場合です。

```bash
module switch SysG SysC
module list
spartition
qgroup
tssrun -p gr20001c -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH/docs/manuscript && TEXINPUTS=/LARGE0/gr20001/b36291/Github/BEACH/build/manuscript-tools/texmf2018/tex//: make'
```

## ドキュメントサイトへの掲載

GitHub Pages のサイトには、この原稿から生成した章ごとのページ（「理論（技術原稿）」の節）と PDF を載せます。
生成は [tools/manuscript_to_starlight.py](../../tools/manuscript_to_starlight.py) が行い、`make docs-pages` の中で実行されます。
章ページの節・式・図・表・文献の番号は、PDF の `.aux` と照合してから書き出すので、PDF と一致します。
本文は pandoc 3（`docs/requirements.txt` の `pypandoc-binary`）で Markdown に変換し、TikZ の図は 1 点ずつ
XeLaTeX でコンパイルして `pdftocairo` で SVG にします。生成物は Git の追跡対象にしません。

原稿を書くときは、次を守るとサイト側の変換が通ります。

- 図は `figure` 環境に `tikzpicture` を 1 つだけ置き、`\caption` と `\label` を 1 つずつ付ける
- 式の番号は `equation` と `align` だけで付ける。番号を付けない行には `\nonumber` を使う
- 参照は `\ref`、`\eqref`、引用は `\citep` を使う

KUDPC では計算ノードに `pdftocairo` がないため、2 段階で実行します。図の PDF 化は計算ノードで、
SVG への変換はログインノードで行います。2 回目の実行は、内容が変わっていない図の再コンパイルを省略します。

```bash
# 計算ノード: PDF、章ページ、図の PDF
tssrun -p <queue> -t 0:30:0 --rsc p=1:t=4:c=4 bash -lc \
  'cd /path/to/BEACH && TEXINPUTS=$PWD/build/manuscript-tools/texmf2018/tex//: .venv/bin/python tools/manuscript_to_starlight.py'
# ログインノード: SVG への変換だけ
.venv/bin/python tools/manuscript_to_starlight.py --skip-pdf
```

## 原稿の根拠と更新

2026 年 10 月 5 日の作業ツリー（基点 `57561eb2fe795d84e06e5b6abaf29ff15c1d8070`）を対象に作成しました。
2026 年 10 月 13 日に第 8 章を実装（commit `7e56bb5`、sheath-model `b4a2521`）に合わせて改訂し、
同日に全章から実装固有の記述を除いて、モデルとアルゴリズムの記述として書き直しました。
このとき、第 4 章の周期遠方作用素（深さ 2 の節点間の作用素と零モードの補正 $s$、$S$）と、
第 5 章の削除済みの連成経路の記述を現行の実装に合わせました。
2026 年 10 月 14 日に、外部シースのモデルそのものを説明する第 8 章を追加しました（sheath-model `26569c1`。
物理の実装は `b4a2521` と同じで、差は `.gitignore` だけです）。仮定、軌道と粒子集団、密度・流束の導出、
第一積分、閉じ方、物理解の条件、数値アルゴリズムと数値例をここへ移しました。
2026 年 10 月 15 日に、旧第 8 章の外部シースの定常 closure を第 8 章の最後の節「BEACH の上端面への接続」へ統合しました。
重複していた上流の帯の模式図は、第 8 章の数値図に置き換えて削除しました。
上流の帯の密度の展開（式 `eq:hlogh`）は、遅い電子が両方向に存在することによる係数 2 を積分の前に明記しました。
同日、サイトの理論ページを原稿へ一本化するため（ADR 0003）、廃止するページにだけあった理論を取り込みました。
第 6 章に吸収時の二次電子放出の節と周期軸での照射 ray の回り込みを、第 7 章に periodic2 でも使える浮遊導体の定式化
（場の評価法から作る影響行列）、closure の統計の注意、バッチ幅の物理的な時間尺度を追加しました。
第 5 章には時間刻みの選び方を、第 8 章には二次電子を含む外部根の更新、零モードを分ける周期場が必要な理由、
上流の帯を埋める物理の見積もりを追加しました。
実装を変更したときは、上の表の対応する章を確認してください。

文献の DOI、著者、題名、掲載情報は Crossref、出版社・研究機関の公開資料、arXiv で照合しました。
本文または abstract による確認であり、全論文の全文精読を完了したとする原稿ではありません。
この初稿は Codex と Claude の支援で作成しています。著者・所属、貢献、研究費、利益相反は確定時に記入してください。
