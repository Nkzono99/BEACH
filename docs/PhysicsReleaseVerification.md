# 物理リリースの検証

新しいperiodic2・outer plasma modelを定量的な議論に使う前に、portable CIとHPC用のリリース検証を
両方実行します。

```bash
make test-l2
make test-physics-release
```

HPC側では、収束データを生成する3本のL1 testに続けて、L3 heavy test、far correctionの正しさ、
MPI電荷収支、MPIによるperiodic cacheの同時生成を検証します。portable CIで実行済みのL2全体は重ねて実行しません。

`manifest.txt`には、各検証工程の経過時間とGNU timeが報告した最大RSSを記録します。最大RSSが既定の8 GiBを超えると、
リリース検証は失敗します。上限は`BEACH_RELEASE_MAX_RSS_KB`で変更できます。

`convergence.csv`は、test logの`BEACH_CONVERGENCE`行から毎回生成します。必要な収束軸が1つでも欠けると、
リリース検証は失敗します。Fortran targetごとの実行時間は、`test_l3-target-timings.csv`と
`far_correction-target-timings.csv`に記録します。

## 基準となる収束表

次の表は、2026-07-11にIntel 2023.2とSystem Bのcompute nodeで取得した、正しさを確認する小規模fixtureの値です。
実運用caseでは、この表とは別にcase固有の収束確認を行います。

| 軸 | 設定 | metric 1 | metric 2 | 判定 |
| --- | --- | ---: | ---: | --- |
| same-time Boris dt | 0.25 / 0.125 / 0.0625 | 5.1600e-3 / 1.2990e-3 | 3.2533e-4 | 予測中点場 + 台形位置の二次収束 |
| panel FMM order | 2 / 3 / 4 / 5 | field: 2.3566e-2 / 2.0761e-3 / 2.3681e-4 / 2.8845e-5 | potential(order 5): 4.7828e-6 | order 5で2e-3未満 |
| rough panel mesh | n=4 / 8 / 16、n=32基準 | field: 1.0645e-2 / 3.0336e-3 / 6.3047e-4 | potential: 1.0805e-2 / 2.9001e-3 / 5.9489e-4 | refinementで減少 |
| rough outer grid | 65 → 129 | field差 3.2521e-6 | potential差 4.3034e-7 | 5e-5未満 |
| rough accessibility | 8x8 → 16x16 | 3.1250e-2 | tolerance 1.0e-1 | tolerance以下 |
| outer orbit dt | 0.1 → 0.05 | energy error 1.2000e-3 | 3.8835e-4 | dt半減で減少 |

periodic far correctionの数値検証は`make test-fortran-far-correction`で実行します。
無限周期operatorは`test_periodic2_infinite_operator`、cold/warm cacheは`test_periodic2_operator_cache`、
MPIによる同時生成は`test_periodic2_operator_cache_mpi`で検証します。速度比較はdebug buildで行う正しさの検証とは分け、
release profileの`make test-fortran-benchmark`で実行します。

実運用caseでは、少なくともmesh、`sim.dt`、FMM order/tolerance、outer grid、height samplingを変えてください。
そのうえで、表面電位、吸収/escape flux、電荷収支、主要な結論が収束することを確認します。
