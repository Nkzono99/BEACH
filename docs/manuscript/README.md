# BEACH の物理モデルと数値アルゴリズム

Lang: [日本語](README.md) | [English](README.en.md)

BEACH の物理モデル、数値手法、先行研究をまとめた、日本語の論文・技術書形式の LaTeX 初稿です。
[beach_models.tex](beach_models.tex) が入口です。章別の本文は `sections/`、文献は
[references.bib](references.bib) にあります。日本語の概要と英語の abstract、目次、数式、
計算サイクル図、参考文献を含みます。

## 内容

1. 表面帯電、月面・無大気天体、先行研究
2. モデル分担、支配方程式、バッチ計算
3. P0 三角形要素、解析積分、自己項・片側電場
4. Direct、treecode、Cartesian FMM、2 軸周期場と Ewald・零モード
5. 同時刻 Boris 更新、衝突、箱境界・周期画像
6. 境界流入、速度分布、照射 ray と光電子放出
7. 電荷収支、帰還 closure、浮遊導体、バッチ幅の安定性
8. 外部シースの定常 closure（零電流根、イオン流体、上流の帯と場の下限、上端面への写像、外部根の更新）
9. 検証、収束、Monte Carlo 誤差、再現性
10. 力・離脱の評価、論文に記載する方法と限界、原稿作成に関する記載

## 書き方の方針

本文は、プラズマ物理と数値計算の一般的な知識を前提に、モデルとアルゴリズムを原稿の中で完結させます。
入力キー、ソースのパス、出力ファイル名、テスト名、commit、計算機環境は本文に書きません。
手法の由来は文献で示します。実装との対応は次の表で管理します。

| 章 | 主な実装・文書 |
| --- | --- |
| 2 計算サイクル | `src/runtime/simulator/bem_simulator_loop.f90`、`docs/Algorithms.md` |
| 3 三角形場 | `src/physics/field_solver/panel/bem_panel_kernel.f90`、`src/physics/field_solver/panel/bem_panel_self_terms.f90`、`docs/DirectSolver.md` |
| 4 treecode・FMM | `docs/Treecode.md`、`src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`、`docs/FMMCore.md` |
| 4 周期遠方場 | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`、`docs/PeriodicFarCorrection.md` |
| 4 零モード | `src/physics/field_solver/periodic/bem_periodic_zero_mode_plan.f90`、`src/physics/field_solver/periodic/bem_periodic_zero_mode_eval.f90`、`docs/PeriodicElectrostatics.md` |
| 5 粒子積分・衝突 | `src/physics/bem_pusher.f90`、`src/physics/bem_collision.f90`、`src/runtime/simulator/bem_particle_stepper.f90`、`docs/BorisPusher.md`、`docs/ParticleEvents.md` |
| 6 流入・光電子 | `src/particles/injection/bem_injection_flux.f90`、`src/particles/injection/bem_photoelectron_injection.f90`、`docs/ParticleSourcesBoundaries.md`、`docs/PhotoelectronEmission.md` |
| 7 電荷更新・導体 | `src/physics/bem_surface_models_conductor.f90`、`docs/SurfaceChargeNumerics.md`、`docs/BatchDurationTheory.md` |
| 8 外部シース | `src/runtime/sheath/bem_surface_current_model.f90`、`src/runtime/sheath/bem_zhao_outflow_refresh.f90`、`docs/ZhaoStationaryClosure.md`、`docs/ZhaoUpstreamBand.md`、sheath-model（`fpm.toml` で commit を固定） |
| 9 検証 | `docs/ValidationGuide.md`、`docs/PhysicsReleaseVerification.md`、`docs/BatchDurationStability.md` |
| 10 力・離脱 | `docs/ObjectForcesDetachment.md`、`docs/PythonPostprocessAPI.md` |

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

2026 年 10 月 13 日、System C の計算ノードで 30 ページの A4 PDF を生成しました。
未解決の引用・参照、行のはみ出し、BibTeX の警告はありません。
既存の TeX Live 2018 に不足していた `xeCJK` と LaTeX3 バンドルは、同世代の TeX Live アーカイブから
`BEACH/build/manuscript-tools/texmf2018/` に一時展開して使いました。既存の TeX インストールは変更していません。
この一時領域は Git の追跡対象外です。通常の TeX 環境では必要パッケージを導入して上の `make` を使います。

ログインノードの Sys module と異なる system へ投入するときは、`SLURM_CONF` で対象を切り替えます。

```bash
SLURM_CONF=/etc/slurm/sysC/slurm.conf tssrun -p gr20001c -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH/docs/manuscript && TEXINPUTS=/LARGE0/gr20001/b36291/Github/BEACH/build/manuscript-tools/texmf2018/tex//: make'
```

## 原稿の根拠と更新

2026 年 10 月 5 日の作業ツリー（基点 `57561eb2fe795d84e06e5b6abaf29ff15c1d8070`）を対象に作成しました。
2026 年 10 月 13 日に第 8 章を実装（commit `7e56bb5`、sheath-model `b4a2521`）に合わせて改訂し、
同日に全章から実装固有の記述を除いて、モデルとアルゴリズムの記述として書き直しました。
このとき、第 4 章の周期遠方作用素（深さ 2 の節点間の作用素と零モードの補正 $s$、$S$）と、
第 5 章の削除済みの連成経路の記述を現行の実装に合わせました。
実装を変更したときは、上の表の対応する章を確認してください。

文献の DOI、著者、題名、掲載情報は Crossref、出版社・研究機関の公開資料、arXiv で照合しました。
本文または abstract による確認であり、全論文の全文精読を完了したとする原稿ではありません。
この初稿は Codex と Claude の支援で作成しています。著者・所属、貢献、研究費、利益相反は確定時に記入してください。
