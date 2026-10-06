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
6. reservoir 流入、速度分布、照射 ray と光電子放出
7. 電荷収支、帰還 closure、浮遊導体、バッチ幅の安定性
8. stationary Zhao、matching-plane、固定点・陰的零モード
9. 検証、収束、Monte Carlo 誤差、再現性
10. `examples/beach.toml` の読み方、力・離脱後処理、研究への接続

付録に実装ファイルとの対応と文献の確認範囲をまとめています。

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

2026 年 10 月 6 日、System B の計算ノードで改訂後の 31 ページの A4 PDF を生成しました（Job 24860895）。
引用 12 件、章・数式の参照、実装ファイルのパスを照合し、未解決参照と行のはみ出しはありません。
既存の TeX Live 2018 に不足していた `xeCJK` と LaTeX3 バンドルは、同世代の TeX Live アーカイブから
`BEACH/build/manuscript-tools/texmf2018/` に一時展開して確認しました。既存の TeX インストールは変更していません。
この一時領域は Git の追跡対象外です。通常の TeX 環境では必要パッケージを導入して上の `make` を使います。
同日の Job 24860896 で、更新した設定例 7 件・境界流入 species 14 件について、省略形と既定値の明示形の
両方が設定検査を通ることを確認しました。既存の境界流入・schema の関連テスト 3 件も成功しています。

確認に用いた組版コマンドは次のとおりです。投入前に `module switch SysG/2022 SysB` と
`module list`、`spartition`、`qgroup` で対象を確認しました。

```bash
tssrun -p gr20001b -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH/docs/manuscript && TEXINPUTS=/LARGE0/gr20001/b36291/Github/BEACH/build/manuscript-tools/texmf2018/tex//: make'
```

## 原稿の根拠と更新

2026 年 10 月 5 日のローカル作業ツリーを対象とします。基点は
`57561eb2fe795d84e06e5b6abaf29ff15c1d8070` で、確認時の未コミットの周期場・零モード修正を含みます。
各章末に仕様・実装の参照を示しました。新規シミュレーション結果や速度比較は生成していません。
2026 年 10 月 6 日に設定例と粒子源の説明を改訂し、境界流入だけの場合は `source_mode` と
`npcls_per_step` を省略する書き方に揃えました。
新しい実装と一致するように本文を更新してから研究結果に引用してください。

文献の DOI、著者、題名、掲載情報は Crossref と出版社・研究機関の公開資料で照合しました。
本文または abstract による確認であり、全論文の全文精読を完了したとする原稿ではありません。
この初稿は Codex の支援で作成しています。著者・所属、貢献、研究費、利益相反は確定時に記入してください。
