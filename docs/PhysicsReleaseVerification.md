# Physics release verification

新しい periodic2 / outer-plasma model を定量利用する前に、portable CI と HPC release gate の
両方を通します。

```bash
make test-l2
make test-physics-release
```

HPC gate は L3、far-correction、MPI ledger、MPI periodic-cache concurrency を逐次実行します。
`manifest.txt` は各 gate の elapsed time と GNU time の最大RSSを記録し、既定の8 GiB予算を超えると
失敗します。予算は `BEACH_RELEASE_MAX_RSS_KB` で明示変更できます。
`convergence.csv` は test log の `BEACH_CONVERGENCE` 行から毎回生成され、必要な収束軸が一つでも
欠けると release gate は失敗します。

## Reference convergence table

2026-07-11、Intel 2023.2、System B compute nodeで取得した小規模 correctness fixture の値です。
これは実運用caseの収束確認を置き換えません。

| 軸 | 設定 | metric 1 | metric 2 | 判定 |
| --- | --- | ---: | ---: | --- |
| Boris dt | 0.25 / 0.125 / 0.0625 | 5.1600e-3 / 1.2990e-3 | 3.2533e-4 | 二次収束 |
| panel FMM order | 2 / 3 / 4 / 5 | field: 2.3566e-2 / 2.0761e-3 / 2.3681e-4 / 2.8845e-5 | potential(order 5): 4.7828e-6 | order 5で2e-3未満 |
| rough panel mesh | n=4 / 8 / 16、n=32基準 | field: 1.0645e-2 / 3.0336e-3 / 6.3047e-4 | potential: 1.0805e-2 / 2.9001e-3 / 5.9489e-4 | refinementで減少 |
| rough outer grid | 65 → 129 | field差 3.2521e-6 | potential差 4.3034e-7 | 5e-5未満 |
| rough accessibility | 8x8 → 16x16 | 3.1250e-2 | tolerance 1.0e-1 | tolerance以下 |
| outer orbit dt | 0.1 → 0.05 | energy error 1.2000e-3 | 3.8835e-4 | dt半減で減少 |

flat surface は `test_periodic2_flat_oracle_diag`、無限周期 operator は
`test_periodic2_infinite_operator`、cold/warm cache は `test_periodic2_operator_cache`、
MPI同時生成は `test_periodic2_operator_cache_mpi` が担当します。

実運用caseでは少なくとも mesh、`sim.dt`、FMM order/tolerance、outer grid、height samplingを変え、
表面電位、吸収/escape flux、電荷収支、主要な結論が収束することを別途確認します。
