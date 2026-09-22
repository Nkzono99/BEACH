# Source-connected モデルの再現検証

Lang: [日本語](README.md) | [English](README.en.md)

この例は #38 の独立配布物と Fortran の静的解を比較し、既存 kinetic oracle で時間依存の確認を行います。
根の一致、有限探索の解集合、動的な収束を別々に評価できます。モデルと CSV の契約は
[専用リファレンス](../../docs/SourceKineticSheath.md)を参照してください。

## 実行

BEACH 開発環境、GNU Fortran、NumPy、Matplotlib が必要です。独立 Python の実行には SciPy と
外部入力 `lunar-sheath-solver-v0.1.0.zip` が必要です。これらの検証用依存を BEACH runtime に追加していません。
以下をリポジトリ直下の計算ノードで実行します。KUDPC のログインノードでは `tssrun` の割当内へ移してください。

```bash
mkdir -p build/source-validation
gfortran -O2 -J build/source-validation -I build/source-validation \
  src/core/bem_kinds.f90 src/core/bem_constants.f90 \
  src/physics/sheath/bem_source_kinetic_sheath.f90 \
  examples/source_kinetic_validation/scan.f90 -o build/source-validation/scan
python examples/source_kinetic_validation/compare_handoff.py \
  _handoff/lunar-sheath-solver-v0.1.0.zip build/source-validation/scan \
  build/source-validation/archive
unzip _handoff/lunar-sheath-solver-v0.1.0.zip -d build/source-validation/reference
PYTHONPATH=build/source-validation/reference/lunar-sheath-solver/src \
  python -m unittest discover -s build/source-validation/reference/lunar-sheath-solver/tests
PYTHONPATH=build/source-validation/reference/lunar-sheath-solver/src \
  python examples/source_kinetic_validation/multiaxis.py \
  build/source-validation/scan build/source-validation/multiaxis
PYTHONPATH=. python examples/source_kinetic_validation/oracle_probe.py build/source-validation/oracle
```

`compare_handoff.py` は 1,211 条件の `native_roots.csv`、差分・checksum の `comparison.json`、
E–G 断面の `existence.png` を出力します。`multiaxis.py` は M–E、τ–G と A fixture の 181 条件を
拡張ソルバー `solve_extended` と比較し、481 / 961 点の走査、全根 CSV、比較 JSON、2 枚の図を保存します。
図のセルは標本点であり、精密な分岐曲線ではありません。どちらも不一致があれば終了コード 1 です。

`oracle_probe.py` は各根の $A_e$ に電子源を合わせ、空間格子、速度格子、外部長、イオン温度、
積分時間、速度範囲を個別に変えます。14 条件の完全な oracle 出力と集計 `comparison.json` を保存します。
既存 oracle は温かいイオンなので、$T_i=0.12$ / $0.03$ eV の比較は冷たい極限への確認です。
終了コード 0 は全 probe の実行完了を表し、全条件の `steady` 判定を意味しません。

## 2026-09-22 の検証結果

独立 ZIP の SHA-256:
`7c3a2e885c5df9bf816c0a54c4facafcb35d76381be12b6a966d1200bbfd11f4`。

- 独立 Python: 86 テスト成功。
- 配布物: 1,211 条件・739 根で状態、根数、Type が一致。尺度調整後の根の最大差 $3.6\times10^{-13}$。
- E–G の 980 条件: 採択 600、解析的除外 255、未解決 125、検出 624 根。
- M–E、τ–G、A fixture の追加 181 条件も拡張 Python と一致（最大差 $2.1\times10^{-13}$）。
  481 点から 961 点への細分化で根数・Type・分類の変化はありませんでした。
- BEACH の公開 main (`c3fc093`) を基準にした最終 L1: Fortran 全 target が成功、Python は 843 件成功・49 件 skip。
  KUDPC の計算ノードで `CC=gcc` を指定して実行し、文書の目次・同梱リファレンスの整合も確認しました。
- source の独立速度積分・Poisson 回帰 7 cases / 60 assertions、online/table の PE なし帯電回帰、
  2 バッチの実行例が成功。

時間依存比較は $M=10,\tau=0.2,G=0.3,E_H=0.1,T_e=12$ eV です。
静的電位は B が 0.274641 V、N が −1.595019 V です。

| 変更軸 | B 電位 [V] | B 判定 | N 電位 [V] | N 判定 |
|---|---:|---|---:|---|
| 基準: nz=32, nv=128, L=3λ, 6 ion transit | 0.303067 | far boundary 未収束 | −1.385123 | far boundary 未収束 |
| nz=64 | 0.305228 | far boundary 未収束 | −1.390159 | far boundary 未収束 |
| nv=256 | 0.282565 | steady | −1.388728 | far boundary 未収束 |
| L=6λ, nz=64 | 0.298639 | far boundary 未収束 | −1.574772 | steady |
| Ti=0.03 eV, nv=256 | 0.282565 | steady | −1.388725 | far boundary 未収束 |
| 12 ion transit | 0.306552 | far boundary 未収束 | −1.385124 | far boundary 未収束 |
| vmax=9ve, 同じ Δv | 0.303135 | far boundary 未収束 | −1.385129 | far boundary 未収束 |

静的根の安定性、全軸の同時収束、PE を含む BEACH 帯電の table/online 収束は未確定です。
14 条件の最大電荷収支残差は $2.01\times10^{-14}$、Gauss 残差は $1.01\times10^{-27}$ C/m²、
速度境界損失率は $5.46\times10^{-4}$ でした。保存則が良好でも遠方境界・格子の収束とは区別します。
根の順序から安定枝を選ばず、未収束 oracle 結果を応答表に変換していません。
この実装だけで #38 の動的検証を完了したとは扱いません。
