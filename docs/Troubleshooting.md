title: トラブルシューティング

Lang: [日本語](Troubleshooting.md) | [English](Troubleshooting.en.md)

# トラブルシューティング

| 症状 | 確認と対処 |
| --- | --- |
| `gfortran` / `fpm`がない | [インストール](Installation.html)に従ってversionを確認。HPCではcompiler moduleをload |
| `beach` / `beachx`がない | `python -m site --user-base`を確認し、user baseの`bin`を`PATH`へ追加 |
| lintは通るが実行時に停止 | errorに表示されるmodel combinationを確認。Fortranのphysics validationはPython schemaより厳しい場合がある |
| output directoryがない | `[output] write_files=true`、`dir`、processのexit codeを確認 |
| `survived_max_step`が多い | `dt`、`max_step`、boxサイズ、注入速度を見直す。未解決粒子を吸収/escape扱いしない |
| 履歴が空 | `history_stride > 0`、`batch_count`、output pathを確認 |
| 履歴が大きすぎる | `history_stride`を増やし、必要時だけpotential historyを有効化 |
| restartが拒否される | model/mesh/species fingerprint、`restart_from`、累積`batch_count`を確認 |
| conductor/dielectric/periodicの組合せエラー | [入力パラメータ](Parameters.html)のsupport matrixを確認。fallbackは行わない |
| FMM far correctionが遅い | `cached_kneq0`のwarm cacheを使い、cold operator生成の頻度を確認 |

問題を再現できる場合は、使用したconfig、BEACH version、compiler/MPI、rank/thread数、errorの全文、
再現用の最小meshを添えてGitHub issueに報告してください。物理的に適用できないmodelを、
数値許容値だけを緩めて実行しないでください。
