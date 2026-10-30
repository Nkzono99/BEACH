# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added

- ドキュメントに、計算サイクル、粒子源・境界、周期場、放出モデル、再開と後処理を説明する日本語の概念図を 24 点追加した。

### Changed

- 公式チュートリアルと `beachx config init` が生成する設定を、プラズマ中の浮遊金属板の帯電に替えた
  （`examples/tutorial_floating_plate.toml`）。密度 1 cm⁻³、温度 1 eV の電子と陽子を z-high 面から境界 reservoir
  （`source_vdf`）で流入させ、x/y 周期（`periodic2`、`cached_kneq0`）の浮遊導体平面が、電子と陽子の電流が釣り合う
  浮遊電位（流入面に対して −3.76 V）まで帯電する 10000 batch（5 s）を追う。`examples/tutorial_floating_analysis.py` は
  電位差、板の電荷、粒子種ごとの吸収数を電流の釣り合いの式と比べる。電子ビームのケース
  （`examples/tutorial_electron_beam.toml`、`examples/tutorial_beam_analysis.py`）は削除した。
- 初期 Zhao 根と mesh 識別子（`mesh_fingerprint`）を run 開始時に一度だけ求め、定期 checkpoint と最終出力の
  `summary.txt` で再利用するようにした（#60）。
  同梱の zero_current 例では checkpoint 1 回あたり約 0.4 s の Zhao 求解が、
  20 万要素の mesh では約 4.6 s の識別子計算が不要になる。`summary.txt` の内容は変わらない。
- `test_mpi_hybrid` の rank 数変更 restart 試験に、1 rank で書いた checkpoint を全 rank で継続する方向と、
  非整数の reservoir 残差と batch 1 の電荷が batch 2 の軌道を曲げる fixture を追加した（#75）。
- 周期物体検証ツールの電荷履歴の batch 合計を、他の集計と同じ `math.fsum` にした（#59）。
- `run.batch.particle_distribution = "shared"`
  で各 rank が global slot 数分の粒子配列・フラグ・slot 表を持たなくなり、
  claim した分だけを保持するようにした（不足時は claim の合間に 1.5 倍以上へ拡張）（#60）。
  4M 粒子/batch・2 rank で rank あたりのピークメモリが約 620 MB から約 350 MB に減る。
  slot 順・乱数・結果は変わらない。
- static 分配の batch 生成で一時配列を経由せず粒子を直接書き込み、
  二次電子の世代追加で全粒子を毎回写し直さないようにした（4M 粒子/batch・1 rank でピークメモリ約 950 MB → 約 600 MB）（#58, #60）。
- counter-based 乱数で同じ鍵の連続 draw を 1 つの PCG64 生成器から引くようにした（値は従来と同一）（#60）。
- 固定電流 closure を使わない実行で raw 電荷の集計と MPI 合計を省き、
  吸収 channel の補正を種ごとの全粒子走査なしで行うようにした。
  reservoir 面の電位格子を layout ごとに面あたり 1 回だけ評価するようにした（#60）。
- 実行終了時の標準出力で全要素の電荷（`mesh charges=`）を出さず、
  `mesh_charge_sum_C`・`mesh_charge_min_C`・`mesh_charge_max_C` を表示するようにした（#60）。
  要素ごとの電荷は従来どおり `charges.csv` にある。
- periodic2 の upper Fourier 切替高さの自動選択で、
  proxy 電荷 pattern ごとの source 和を 1 回だけ作るようにした（従来は候補 16 高さごと）。
  mesh ごとの中央値 pattern は要素を mesh id 順に 1 回並べて求め、mesh 数に比例した全要素走査をなくした。
  選択される高さと mismatch は変わらない（#60）。
- periodic zero-mode の state 更新は確保済み配列を再利用し、
  gauge 高さの電荷原始関数を更新時に 1 回だけ計算するようにした（評価点ごとの区間探索をなくした）（#60）。
- FMM plan 構築で、読まれていなかった遠方 node 一覧を作らず、周期 wrap で動かない panel の幾何を作り直さず、
  octree 分割の作業配列を固定長にした。finite_images の seam 探索のソートを O(N log N) にした。
  結果は変わらない（#60）。
- cached_kneq0 operator の checksum を target block ごとに計算し、
  operator 全体の byte 複製をなくした（checksum 値と cache 形式は同じ）。
  MPI を使わない run では operator の flatten/reshape 複製をしない（#60）。
- periodic2 の run 開始時、
  canonical cell への平行移動が要らない mesh では幾何と collision grid を作り直さないようにした（公開の
  `prepare_periodic2_collision_mesh` は従来どおり常に再構築する）（#60）。
- 再直交化 QR の実装を `bem_regularized_qr` の 1 つに統一し、
  FMM plan と Ewald 補助データの reset を型の既定値で行うようにした（#58）。
- `ObjectProbe.vertical_path` の適応細分化で、評価済みの高さ（格子点と分割しなかった区間の中点）
  を再評価しないようにした。1 高さで評価点上限を超える target は評価点方向にも分割する。
  高さごとの位置エネルギーは評価のまとめ方に依存しない行ごとの和で求め、従来値との差は最大 2 ULP（#60）。
- `Beach` の転送メソッドを別の `step` で呼んでも、
  同じ周期条件・FMM 設定の kernel が開いていれば電荷を一時的に入れ替えて使い、
  step ごとに kernel を増やさないようにした（#60）。
- `animate_history_mesh` は選んだフレームだけを読み、履歴全体を展開しない。`compute_potential_slices`
  は 3 面の全点配列を作らず chunk ごとに評価する（#60）。
- 設定 1 回の読込で schema validator を 1 回だけ構築するようにした。release 判定の仕事積分表、
  boxplot の重み付き四分位数、runops の stderr 末尾の読み込みもそれぞれ 1 回にした（#60）。
- `FieldKernel.source_count` を追加した。`Beach.get_mesh_charge` は `mesh_triangles.csv` がなくても電荷を返す（#60）。
- `Beach` のメソッドが ParamSpec により元の関数の引数と戻り値の型を型チェッカーと補完に示すようにした（#58）。
- Python の formatter（`ruff format`）と import 整列を path 単位で段階導入し、第 1 段階として
  `beach/native`・`beach/plot`・`beach/tooling`・`beach/facade.py` を整形した（#58）。`make python-lint`（L0）
  がこれらの format を検査し、`make python-format` で整形できる。dev 依存の Ruff は minor の範囲（0.17 系）
  に固定した。
- 周期 object 検証ツールの job 構成（script、job 名、walltime、担当 case）
  と object 評価 plan（charge-history batch と target mesh）を 1 か所の定義から導出するようにした（#58）。
  生成される job script と検証内容は変わらない。

### Fixed

- `beach --check-config`（Fortran preflight）が、JSON schema だけが拒否していた
  `particles.tracking.max_steps_per_particle = 0`、`output.enabled = false` での `[run.restart]`、`inflow` を持たない
  `volume_seed` 種の `sampling.target_macro_particles_per_batch` を受理していた。
  preflight でも公開キー名を示して停止する。
- periodic2 FMM の点評価で、周期下端の半 ulp 未満だけ下の点が `modulo` の丸めで上端 `lo+L` に写り、
  leaf 探索で再度折り返されて反対端の leaf が選ばれていた。この点の電場・電位は 1 周期離れた位置の展開で評価され、
  桁違いに誤っていた。評価点と source を `[lo, lo+L)` へ一度だけ折り返す（#60）。
- OBJ メッシュの読込で 1024 文字を超える行を切り詰めず、行末まで読むようにした。
  長い面行が短い別の面として取り込まれることがなくなった（#59）。
- Zhao 外部シースの outflow refresh で、解けた根の closure が棄却されたとき、
  外部根を前回の受理状態のまま残すようにした（#59）。これまでは外部根だけが進み、
  次回の refresh と restart の再構成がずれていた。
- 周期物体検証ツールの stage / probe-periodic-oracles が途中で失敗したとき、その試行の生成物を消すようにした（#59）。
  receipt のない作業域が次の試行を止めることがなくなった。

## [4.1.0] - 2026-10-30

### Added

- `kind = "sphere"` に `tessellation = "icosphere"` と `subdivisions` を追加した。正二十面体を再帰分割した
  `20 * 4**subdivisions` 要素の球面で、要素面積の最大/最小比は約 1.3 に収まる。`beachx model close-pack` は
  `--sphere-subdivisions` で icosphere の球を生成する。
- `kind = "disk"`, `"annulus"`, `"cylinder"`（キャップ）に `tessellation = "uniform"` を追加した。各リングの点数を
  半径に比例させ、要素面積の最大/最小比を分割数によらず約 1.5 に保つ。円柱のキャップは `n_r` で半径方向に分割できる。
  既定の `"polar"` は従来の極座標格子（要素面積が半径に比例）のまま。
- `fixed_current_history.csv` に、全 species の表面正味電荷を closure 前後で記録する `surface_net` 行（`species_idx = 0`）を
  追加した。`zero_current` の run で closure 後の正味が 0 にならない batch があると警告する。

### Changed

- `kind = "sphere"` の既定の分割を icosphere（`subdivisions = 2`、320 要素）にした。`n_lon` または `n_lat` を書いた球は
  従来どおり UV 球になる。どちらも書いていない球は、従来の UV 球（24×12、528 要素）から icosphere に変わる。
  従来のメッシュを保つには `tessellation = "uv"` を書く。
- 通常の開発テストから周期検証 workflow の全ケースを分離し、代表 3 ケースを残した。
  全ケースは `make test-validation-workflow` と CI で実行する。
- L2 で共有 field kernel を build して Python/native 数値比較を実行し、lint CLI の native 設定契約も確認する。
  `make test-field-kernel-cache` は周期解析解と有限 shell の検証も実行する。
- 周期導体のテストを遠方補正 gate に、Laplace 偏微分の速度比較を release benchmark に移した。
  同じ schema 値の再比較と文書書式・内部呼出し回数の固定を削減し、C header の compile probe を compiler 間で利用できる形式にした。
- periodic2 image-layer の FMM 回帰を panel FMM と同じ相対誤差 0.2% で確認する。
- L3 tier を L2 に統合した。`make test-l3`、`make test-heavy`、`make test-fortran-heavy` を削除し、
  free-space FMM の精度比較は `make test-l2` で実行する。`make test-fortran-release-correctness` は収束データを出す
  `test_dynamics_basic` と `test_coulomb_fmm_core_panel` だけを実行し、release manifest は schema 3 で
  `release_correctness.*` と `release_correctness-target-timings.csv` を記録する。
- cached operator を作らない periodic target 4 本（`test_periodic_nonzero_reference`、`test_periodic_nonzero_upper_vacuum`、
  `test_periodic_shell_mean`、`test_coulomb_fmm_core_periodic`）を遠方補正 gate から L1 に移し、
  `test_panel_near_correction` を `test_panel_geometry_near` に統合した。
- L1 の重い 3 target を同じ oracle と許容誤差のまま短縮した。auto solver の要素数表は公開 pure 関数で直接検査し、
  外部根の更新 test は平板の一様場で結果が変わらない時間刻みに広げ、sheath-model 自身が検査する根の棄却 2 件と
  平均エネルギー近似の根探索を Zhao 根 test から除いた。
- `test_simulator` を `test_simulator_loop`、`test_surface_charge_closure`、`test_simulator_failure_context`、
  `test_soft_discard` に分け、粒子 batch の生成は `test_particle_pool`、適応 sub-step は `test_particle_stepper` と
  `test_step_growth` に移した。cached snapshot の adaptive batch は `test_periodic2_adaptive_batch` に分けた。
  `test_field_kernel_cache_c` と `test_output_writer_potential` は、固有の検査を `test_periodic2_cached_snapshot` と
  `test_output_writer_io` に移して削除した。`make test-field-kernel-cache` の Fortran 部分は `test_periodic2_cached_snapshot` を実行する。
- 複数 layer で同じ契約を検査していた Fortran test（panel FMM と direct の比較、free box dual tree、round-robin 順、
  redistributed reflect の counter 要求、top-reference history の書式、上端 Neumann 場、ABI 番号の literal、
  examples 固有の値）を、契約を持つ layer の 1 か所に寄せた。子 process probe、cache の後片付け、直接 panel 和の
  oracle は共通 helper（`test_support`、`test_panel_oracle`）にまとめた。
- Python の gate は file 列挙ではなく pytest marker で選ぶ。native kernel と Fortran 設定 checker を使う test は
  fixture から `native` / `config_native` が付き、`make test-python` は native test を選ばない。そのため
  `build/` に残った古い共有 kernel を使わない。periodic cache を生成する test は `field_kernel_cache` で選ぶ。
  これにより、gate 間で約 180 件あった二重実行がなくなった。
- `tests/python/` を `io`、`analysis`、`native`、`facade`、`plot`、`config`、`cli`、`tooling`、`repo`、`workflow` の
  directory に分けた。周期検証 tool の test は工程別の 8 file に分けた。Fortran parser test と同じ preflight 規則を
  Python 側で再検査していた case は、lint が Fortran の判定と診断文をそのまま返すことを示す代表 case に絞った。
- `dev` extra に `runops>=2.0,<3` を追加し、runops adapter の test を CI で実行する。
- periodic2 zero-mode plan の高さのソートと区切り探索を O(N log N) にした（10000 panel で約 6 s から 0.01 s）。
  `mesh.parts` の結合は各 part を一度だけ copy する。box target FMM の target tree は深さの上限から決まる
  4681 node だけを確保する。MPI の配列集約は `MPI_IN_PLACE` で作業配列を使わない。
- `ObjectForcePath.evaluate_release` は全区間の評価点をまとめて積分し、鉛直経路の box 検査は変位の最小・最大だけを
  使う。結果は変わらない。
- Python の解析・描画 test を、形状ではなく値で検査するようにした（#76）。対象は電位 slice の全格子と座標対応、
  Coulomb force matrix の各成分、boxplot の加重四分位と描画量、GIF の frame 内容、`far_correction`
  正規化後の kernel option、history 全体での kernel 1 個の再利用である。
- `beachx lint` の設定契約 test を整理した（#77）。schema が拒否する入力は native checker なしで schema 診断を照合し、
  preflight が拒否する入力は Fortran 診断を照合する。
  batch duration の overflow/underflow ケースは schema で先に拒否されていたため、native 検査に届く入力に直した。
  `probe_library` の返却値と charge-ledger 自己 probe の期待診断も検査する。
- `surface_side = "outward_closed"` の辺の照合を全辺対全辺の探索から許容差格子の hash に置き換え、
  part ごとに自身の要素だけを処理するようにした（#62）。8192 三角形の閉曲面で 11.8 s → 0.009 s。判定結果、エラー、
  vacuum normal は変わらない。
- native field kernel（`FieldKernel`）は、
  非周期 plan の FMM 木・演算子・展開を最初の FMM 評価まで作らないようにした（#63）。
  exact-direct だけを使う評価（`resolved_field_solver = "direct"`、ObjectProbe の primary kernel）
  では構築と電荷更新が大幅に軽くなり、32768 panel ではピークメモリが 4.1 GB から 17 MB に減る。結果は bit 単位で同じ。
  同じ handle の評価を複数の host thread から同時に呼ぶ場合は、呼出し側で直列化する必要がある。
- 数値テストの検出力を強化した（#74）。NaN/Inf を相対誤差比較より前に拒否し、
  panel kernel テストの標本を overflow しない乱数へ置き換え、
  非零原点 FMM box で box 内 target の E と phi を直接総和と比較し、
  Python 側の許容差を物理スケールに合わせて明示する。
- runtime テストで粒子状態・空間分配・出力値を検査するようにした（#75）。legacy sampler の分布モーメント、
  Whipple 収率曲線、append の全フィールド、fixed-current 補正の要素分配、テンプレートの寸法と cap、
  履歴 CSV の値と列順、reservoir の絶対粒子数が対象。
- 並列条件と再開を検証するようにした（#75）。OpenMP replay が 4 スレッドで走ったことを確認し、`make test-mpi`
  は要求した rank 数（`BEACH_TEST_MPI_RANKS`）と一致しなければ失敗する。
  多 rank の checkpoint を 1 rank で継続した結果を中断なしの実行と比較し、
  adaptive sub-step の着地点を細かい固定 step の参照と比較する。

### Removed

- Fortran の `particle_step_result` から `field_eval_count` と `collision_query_count` を削除した。src からは読まれず、
  test が内部の呼出し回数を固定するためだけに使っていた。
- Fortran の `find_first_hit` は `status` を必須にし、`status` を省略したときに不完全な衝突照会で停止する経路を削除した。
  production の呼出し側はすべて `status` を受け取って扱っている。`compute_periodic_shift_bounds` は `bem_collision` の
  内部手続きにした。
- 3.8.0 で batch の粒子が `generate_batch_slot` から作られるようになってから呼ばれていなかった Fortran の batch sampler
  `init_random_beam_particles`、`sample_reservoir_face_particles`、`sample_photo_raycast_particles` と、それを包む
  `sample_species_state` / `sample_photo_species_state` を削除した。`sample_flux_weighted_normal_component` の
  `uniform` と `sample_photo_emission_velocity` の `rng_seed` / `rng_counter` は必須にした。注入面と photo ray の
  sampling は counter-based の `sample_reservoir_face_particle_counter`、`sample_photo_ray_slot`、`photo_ray_macro_weight`
  で検査する。
- Fortran の `init_particle_batch_from_config` と `build_batch_injection_layout` は `source_plan` を必須にした。
  省略時に呼出しの中で plan を作る経路は test だけが使っていた。
- charge ledger の累積時に粒子数の非負を再検査する停止を削除した。粒子数は src の加算でしか増えず、
  checkpoint から読む値は resume の読込みが検査している。

### Fixed

- 反射電子を外部源に含めた `zero_current` の固定電流 closure で、追跡した反射電子の脱出が根の脱出束を超えた
  batch に、超過分の正の正味電荷が残っていた（#93）。PE の脱出目標は根の脱出束で 0 に止まるのに、電子の吸収目標へは
  追跡した脱出の全量を足し戻していた。電子へは放出電子側が実際に差し引いた量だけを足し戻す。
- 正方形の `plate_hole` で `n_theta` が 4 の倍数のとき（既定の 1×1、`n_theta = 48` を含む）、角に当たる放射線の
  交点が丸め誤差で角とずれ、面積が丸め誤差程度の三角形が角に生成されていた。角に当たる交点を角に一致させ、
  その要素を生成しないようにした。該当する設定では要素数が減る（1×1、`n_theta = 48`、`n_r = 4` では 386 から 384）。
- `cylinder` の側面最上段の頂点の z 座標を上面キャップと同じ式で計算し、側面と上面キャップが外周の頂点を
  ビット一致で共有するようにした（従来は丸め誤差で 1 ulp ずれる場合があった）。
- 再開時の `charges.csv`、`macro_residuals.csv`、`charge_ledger.csv` を Fortran の自由形式入力で読んでいたため、
  空欄、`/`、繰り返し指定（`3*0.0` など）や次行の値を受理し、前行の値を再利用する場合があった。各行が所定数の
  数値だけからなることを確認してから変換する。空白、`D` 指数、旧 2 列の残差ファイルは従来どおり読める。
- Python の `charge_ledger.csv` reader が、旧版で存在しない列と、存在する列の空欄・欠けた行を区別せず既定値で
  補っていた。既定値は列がない場合だけに使い、空欄や列数の合わない行はエラーにする。
- `FieldKernelOptions` の `box_min` / `box_max` に 3 成分でない配列を渡すと、そのまま共有ライブラリへ渡っていた。
  有限な 3 成分であることと、両方が指定されていることを確認する。
- 粒子種ごとの `sampling.volume_macro_particles_per_batch` が範囲内でも、合計が 32 bit 整数を超えると誤った
  粒子数になっていた。合計が上限を超える設定は停止する。
- `beachx coulomb --json` の `force_matrix` が source 行 × target 列の行列を target 名と行で対応させていたため、
  source と target の数が異なると値が誤対応・欠落していた。target ごとに全 source の力を報告する。
- `beachx lint --beach ./checker` の相対パスが `PATH` 上で探されていた。起動ディレクトリから解決する。
- `tools/golden_regression.py --rtol nan` が数値差を常に合格させていた。NaN・Inf・負の `--rtol` は終了コード 2 で拒否する。
- native preflight（`beach --check-config` と通常実行）が `volume_seed` 粒子種の `mass_kg > 0`、`sampling.weight > 0`、
  有限で非負の温度、各成分 `source.pos_low <= source.pos_high` を生成前に診断するようにした（#68）。
  これまで reservoir 流入のない `volume_seed` は検査を通らず、質量は最初の batch で停止し、負の温度は 0 に丸められていた。
  温度 0 と点状・平面状の生成領域は引き続き受け入れ、slot 生成は検査済みの値を再検査・補正せずに使う。
- 同じ `output.dir` を再利用する新規実行が、前の run の history・`mesh_potential.csv`・`performance_profile.csv` と
  定期 slot を残し、解析が前の run の history を、再開が前の run の slot を選ぶ問題を直した（#79）。新規実行は書き始める前に、
  前の run の最終出力と両 slot から `summary.txt` と `checkpoint_complete.txt`、`checkpoint_latest.txt` を削除して
  再開候補から外し、既知の名前の任意出力だけを削除する。
- history より古い checkpoint から再開すると、history の未来行の後へ追記して batch が重複・逆転する問題を直した（#80）。
  再開時に各 history から checkpoint の `batches` より後の行と改行で終わらない末尾の行を原子的に除き、
  0 byte の history には header を書く。定期 checkpoint を書く batch では完了 manifest より先に history を flush する。
- `build.sh` が fpm を `OMP_PROC_BIND=false` で起動するようにした。KUDPC の job 既定の `OMP_PROC_BIND=close` では
  fpm が起動する test と `beach` が 1 core の affinity を継承し、OpenMP test が数倍遅くなっていた。
- `test_upper_fourier_switch` が固定の cache directory を再利用して 2 回目以降に cold build を検査しなくなる問題を直した。
- opt-in 周期力テストの zero-mode 設定と、平面 oracle receipt の backend 期待値を現行の明示的な周期場契約に合わせた。
- 平面 oracle が mesh 解像度ごとに異なる cache identity を要求する誤検知を修正した。
  同じ box-cell 演算子の再利用を認め、評価ごとの path・fingerprint・SHA256 の照合を維持する。
- 近接場の有限 shell テストで経路の標本点を 9 点から 33 点に増やし、1% の仕事誤差条件を満たす積分解像度で比較する。
- `beachx workload` の rank 別粒子数を本体の `static` 分配に合わせた。種・source・reservoir 面ごとに
  round-robin に並ぶ global slot の総数を一度だけ分割した連続区間から数え（以前は source・面ごとに分割して合計していた）、
  `shared` 分配では均等負荷の仮定であることを `rank_count_basis` に示す（#64）。
- `tools/golden_regression.py compare` が有限値から NaN/Inf への変化、新しい側だけにある case・ファイル、
  `--rtol 0` での表記差（`1.0` → `1.00` など）を成功扱いしていた問題を直した。`--rtol 0` は volatile 行を除いた
  テキストの完全一致を要求し、`--rtol` が正のときは数値を相対差で比べて非有限値は同じ表記どうしだけ一致とする（#67）。
- runops adapter が、run の `input/` にある OBJ を絶対 path（または同じ実体への symlink）で指す mesh part を
  自分自身へ copy しようとして `SameFileError` で止まる問題を直した。同名で内容の異なる OBJ は従来どおり
  `part<index>-<name>` へ配置する（#70）。
- `ClosePackSpec.validate` が `radius`、`floor_z`、`box_height`、`top_clearance` の NaN/Inf を受理していた問題を直し、
  有限であることを同じ検査で要求する（#71）。
- `ObjectProbe.vertical_path` の適応細分化が、区間ごとの判定を通ると返す grid の経路全体の work/potential 一致を
  確かめずに終了し、細分化の余力を残したまま `work_potential_mismatch` で未収束を返していた問題を直した。
  経路全体の一致を loop 内の受理条件に含め、不一致なら区間幅に比例して許容を配分し、超えた区間を細分化する（#81）。
- physics release gate が、列の欠けた行・非有限の metric・空の acceptance を含む `convergence.csv`
  でも成功 manifest を作る問題を直した（#77）。各行に 6 列を要求し、`metric_1`・`metric_2` は有限値、`metric_3`
  は空欄か有限値、`acceptance` は空でない値とする。満たさない行があれば `convergence.status=failed_malformed_row`
  で失敗する。

## [4.0.0] - 2026-10-23

Python API と CLI の破壊的変更を含む。移行に必要な変更点を以下に示す。

### Changed

- 外部境界の基底契約を設定の preflight で一度解決して保持し、粒子生成・run 制御・summary で共有する。
  粒子源計画は整数の種別を保持し、バッチの slot 構築で `source_mode` を解釈し直さない。
  Fortran から手組みする `app_config` も `boundary_contract` に解決済みの値を持たせる。
  `reservoir_face_velocity_correction` の任意引数 `boundary_contract` は削除した。
- 粒子追跡・衝突・境界 event の設定を `particle_domain_config`、photo ray の設定を `photo_raycast_config` に分けた。
  `sim_config` はそれらを継承して既存の属性を保ち、runtime は physics に必要な親 component だけを渡す。
- `workload` は、粒子数の計算に不要な光電子の放出条件、粒子境界の作用、一様外部場の組合せを再検証しない。
  光電子は ray 数の上限を見積もる。設定の検証には `beachx lint` を使う。流束・重み・バッチ幅の計算と、
  計算が未定義または曖昧になる入力の検出は見積もり側が引き続き担当する。
- **BREAKING**: Python package を依存の層ごとのサブパッケージに分けた。`beach.io`（出力 reader・データ型・`summary.txt`
  parser、numpy のみ）、`beach.native`（field kernel と zero mode の ctypes binding）、`beach.analysis`（電位・電場・力・object
  解析）、`beach.plot`（matplotlib の図とアニメーション）、`beach.tooling`（workload・close-pack・runops adapter、numpy 不要）、
  `beach.facade`（`Beach`）。`beach.fortran_results` と import 互換の
  `beach/cli_*.py` 6 本は互換層を残さず削除した。サブパッケージから import していたコードは、`from beach import ...` または
  新しい module（例: `beach.summary` は `beach.io.summary`、`beach.workload` は `beach.tooling.workload`、
  `beach.fortran_results.periodic.Periodic2Config` は `beach.native.periodic2.Periodic2Config`）へ書き換える。runops の entry point は
  `beach.tooling.runops_adapter:Adapter` になった。
- **BREAKING**: `FortranRunResult`、`load_fortran_result`、`list_fortran_runs` を `RunResult`、`load_run`、`list_runs` に改名した。
  `RunResult` の統計・形状・状態記録・診断・場の出自は、それぞれ `summary`、`mesh`、`receipts`、`diagnostics`、`field` に分けた。
  例として `result.triangles` は `result.mesh.triangles`、`result.batches` は `result.summary.batches` になる。
  `FortranChargeHistory` は `BatchHistory` に一般化し、電荷と電位の履歴は `result.history.charge` と
  `result.history.potential` で得る。旧名の alias と旧属性の property は残さない。保存済みの出力形式は変更しない。
- **BREAKING**: `FieldKernel.from_result(...)` を削除した。native binding が run の解析層に依存しないよう、run から kernel を作る処理は
  `beach.analysis.field_kernel.field_kernel_from_result(...)` に移した。`Beach` からは `Beach.kernel(...)` を使う（下記）。
- **BREAKING**: `Beach` を、1 つの出力ディレクトリの結果・`RunContext`・native field kernel を持つセッションにした。
  `Beach(output_dir, *, config_path=None, library_path=None)` で作り、`run.context` が `beach.toml` の解決を一手に持つ。
  `Beach.field_kernel(...)` は alias を残さず `Beach.kernel(*, step=-1, periodic2=None, theta=None, leaf_max=None, order=None)`
  に置き換えた。`kernel()` は同じ電荷と plan に対して 1 つの kernel を cache し、セッションが所有する（呼び出し側は閉じない）。
  `close()`、`reload()`、`with Beach(...) as run:` の終了で cache 中の kernel を閉じる。`Beach.config_path` 属性は削除した。
  解決済みの path は `run.context.config_path` で得る。
- **BREAKING**: `RunContext.from_value` は、`Beach` から `config_path` 属性を推測せず、`Beach.context` をそのまま使う。
  `Beach` の代わりに独自 object を渡すコードは `result` に加えて `context`（`RunContext`）を持たせる。
- **BREAKING**: `Beach` の委譲メソッドは、委譲先のモジュール関数の signature・既定値・docstring をそのまま持つ
  （`help(Beach.plot_potential)` で確認できる）。既定値はモジュール関数だけが持つ。この結果、委譲先の引数はすべて
  メソッドでも使える（例: `Beach.calc_coulomb` の `quadrature_order`・`config_path`・`library_path`・`kernel`）。
  メソッドは `library_path` を省略するとセッションの値を使い、kernel を評価する呼び出しには cached `Beach.kernel()` を渡す。
  `kernel=`、`config_path=`、セッションと異なる `library_path=` を明示した呼び出しは、その呼び出しだけの kernel を構築する。
  保存済み `mesh_potential.csv` を使う `compute_potential()` / `plot_potential()`、`quantity="charge"` の `animate_mesh()`、
  `workers > 1` の `plot_coulomb_force_matrix()` は kernel を構築しない。
- runops adapter の `summarize` は `summary.txt` を `beach.io.summary.parse_summary_text` で読む。重複 key や空 key は黙って
  上書きせず `ValueError` になる。
- **BREAKING**: `reference_point` の既定を、電位の計算と描画で共通の `None`（絶対電位）にした。
  `plot_potential_mesh`、`plot_potential_slices`、`animate_history_mesh(quantity="potential")` と `Beach.plot_potential`、
  `Beach.plot_potential_slices`、`Beach.animate_mesh` は、これまで species 1 の注入面中心を基準にした電位差を既定で描いていた。
  `compute_potential_*` は以前から `None` が既定なので、同じ引数で計算と描画の値が一致する。`mesh_potential.csv` があれば
  potential mesh の描画は kernel を使わずに保存値を使う。従来の図は `reference_point="species1_injection_center"` で得られる。
  `beachx inspect` の potential plot と `--recompute-potential` の要約、`beachx slices`、`beachx animate --quantity potential`
  も絶対電位になる。
- **BREAKING**: 1.0.0 で非推奨にした console script 6 本（`beach-inspect`、`beach-animate-history`、
  `beach-plot-coulomb-force-matrix`、`beach-plot-potential-slices`、`beach-estimate-workload`、
  `beach-plot-performance-profile`）と `beach/cli/legacy.py` を削除した。CLI は `beachx` だけになる
  （`beachx inspect`、`beachx animate`、`beachx coulomb`、`beachx slices`、`beachx workload`、`beachx profile`）。
  旧 script 名を残すためだけの `examples/generate_closepack_config.py` も削除した。`beachx model close-pack` を使う。
  生成した TOML の先頭コメントは `# Generated by beachx model close-pack` になる。
- **BREAKING**: `beachx` の保存オプションを、1 ファイルは `--save PATH`、決まった名前の複数ファイルは
  `--save-dir DIR` に統一した。旧名の alias は残さない。

  | コマンド | 旧 | 新 |
  |---|---|---|
  | `inspect` | `--save-bar A --save-mesh B --save-potential-mesh C` | `--save-dir DIR`（`charge_bar.png`、`charge_mesh.png`、`potential_mesh.png` を書く） |
  | `animate` | `--save-gif PATH` | `--save PATH` |
  | `mobility`、`kernel-forces` | `--save-csv PATH` | `--save PATH` |
  | `object-detachment` | `--output-dir DIR` | `--save-dir DIR` |

  `inspect --save-dir` は `--show` と同じ 3 枚を書く。`mesh_potential.csv` が無い出力では電位 mesh のために電位を計算する。
- **BREAKING**: `beachx` の標準出力を `key=value` 行に統一した。スカラーは 1 行 1 組（`absorbed=7` と `escaped=3` は
  別の行）、浮動小数点数は `%.6e`、真偽値は `true` / `false`、値なしは `null`、配列は空白なしの JSON で書く。表は 1 行
  1 レコードの `field=value` 列で、`batch_details=` や `top_regions_by_<metric>=` のような見出し行と字下げ行は無くなった。
  保存先は `saved=<path>`、複数なら `saved_<name>=<path>` で、`saved_gif=`、`saved_csv=`、`saved_dir=`、
  `saved_potential_slices=`、`saved_bar=`、`saved_mesh=` は無くなった。`inspect` の adaptive 幅は `summary.txt` と同じ
  key（`periodic2_max_nonzero_mode_potential_step_V` など）で 1 行ずつ出す。`mesh_triangles.csv` が無いときの
  「mesh plot を省略」の通知は標準エラーに出す。
- **BREAKING**: `beachx` のオプション値の検査は引数の解析時に行い、`argument --view-elev: must be finite` の形の usage error
  （終了コード 2）になる。library 関数が同じ検査をする値（`mobility` の密度・摩擦係数、`animate --fps` など）は CLI で
  検査せず、library の `ValueError` の文言で終了する。どのコマンドでも library の `ValueError` と `FieldKernelError` は
  traceback ではなく文言で終了する。
- `beachx coulomb`、`beachx mobility`、`beachx kernel-forces` の `--step` は `object-detachment` と同じく `final`
  （`charges.csv`）を受け付ける。
- **BREAKING**: `beach.cli` のサブコマンド module は `configure(parser)` と `run(...)` だけを持ち、各 module の `main`、
  `build_parser`、`add_subparser`、`COMMAND_NAME`、`LEGACY_COMMAND_NAME` を削除した。コマンドの表と dispatch は
  `beach/cli/_registry.py`、共通オプションは `beach/cli/_options.py` が持つ。Python から CLI を呼ぶ場合は
  `beach.cli.main.main(["inspect", ...])` を使う。`beach.cli.config` と `beach.cli.model` は `config_init`、
  `config_diff` と registry の入れ子に置き換え、`beach.cli.lint.run_lint` は `run` に改名した。
  `beach.cli.estimate_fortran_workload` は `beach.tooling.workload` の名前を re-export しない。
- `summary.txt` のスカラーと `RunResult` の属性の対応を `SUMMARY_FIELDS` にまとめた。
  `tools/generate_summary_reference.py` は同じ表から出力リファレンスの対応表を更新する。
- `beachx` の解析コマンドは `Beach` セッションを使い、`inspect --recompute-potential` の計算結果を要約と図で共有する。
- `examples/tutorial_beam_analysis.py` は、生 CSV の列番号に依存せず `load_run` と履歴の公開 API を使う。

### Added

- runops adapter に `progress(run_dir)` を追加した（runops >= 2.5 が `runo sync` で呼ぶ）。manifest の job の
  stdout から現在の batch と `run.batch_count` を返し、再開した run は先頭の `resuming_from_batches=` を `start` として
  返す。runops はこれを時刻付きで保存し、再開点から数えた速度と残り時間を `runo status` / `list` に出す。
- wheel の agent context（`beach.agent_context`）に `plugins/beach-context/docs/` を含めた。skill が参照する
  `../../docs/*.md` が runops の展開先（`beach-docs`）で読めるようになる。
- `beachx inspect`、`mobility`、`kernel-forces`、`coulomb`、`object-detachment`、`workload` に `--json` を追加した。
  `key=value` 行と同じ内容を JSON 1 オブジェクトで標準出力に出す。表は配列（`objects`、`batch_details`、
  `mesh_sources`、`force_matrix`）、非有限値は文字列 `"inf"` / `"-inf"` / `"nan"` になる。`object-detachment` は
  保存先と要約値だけを出し、詳細は従来どおり `summary.json` に書く。`coulomb` は target ごとの力行列も出す。

- `compute_potential_mesh`、`compute_potential_points`、`compute_potential_slices`、`compute_electric_field_points`、
  `trace_field_lines`、`calc_coulomb`、`calc_object_forces_kernel`、`plot_potential_mesh`、`plot_potential_slices`、
  `plot_field_lines_3d`、`plot_coulomb_force_matrix`、`animate_history_mesh` に `kernel: FieldKernel | None = None` を追加した。
  `None` は従来どおり呼び出しごとに kernel を構築して閉じる。渡した kernel は同じ run の panel 数と全電場 plan を確認したうえで
  その plan のまま評価し、閉じない（`periodic2`・`theta`・`leaf_max`・`order`・`library_path` は使わない）。`step` を持つ関数と
  アニメーションは呼び出し中だけ電荷を差し替えて元に戻す。`calc_coulomb` は全 mesh の kernel で source 以外の電荷を 0 にして
  評価し、一様外部場を除く。Coulomb 力行列は worker process 間で kernel を共有できないため、`workers == 1` のときだけ使う。
- `compute_potential_mesh` と `plot_potential_mesh` に `source: "auto" | "saved" | "kernel" = "auto"` を追加した。`"saved"` は
  `mesh_potential.csv` だけを返し、`"kernel"` は保存値があっても再計算する。`"auto"` は従来の選択（保存値があり、
  `periodic2`・`reference_point`・`kernel` の指定が無ければ保存値）。`plot_mesh_source_boxplot(quantity="potential")` の
  履歴 step は、result の保存値を消す代わりにこの `source="kernel"` で再計算する。
- `FieldKernel.options`、`FieldKernel.source_charges`、`FieldKernel.closed` の読み取り専用 property を追加した。
- `BeachScene.from_result(..., library_path=None)` と `BeachScene.library_path` を追加した。scene の kernel メソッドは
  `library_path` を省略するとこの値を使う。`Beach.scene()` はセッションの `library_path` を渡す。
- 公開解析・描画関数の docstring の冒頭に、native kernel を構築する（heavy）か保存済みの値だけを使う（light）かを書いた。
- `top_reference_history.csv`、`fixed_current_history.csv`、`matching_plane_spectrum_history.csv` と
  `performance_profile.csv` の公開 reader を追加し、`RunResult` からも参照できるようにした。
- Sphinx autodoc と napoleon で Python API リファレンスを生成する。`make docs-python` で単独生成し、
  `make docs-pages` は Fortran API と並べて `python/` に同梱する。Sphinx は開発用の `docs` extra に置く。
- `dev` extra と `make python-type-check` を追加した。L0 は軽量な summary・設定処理3ファイルを mypy の strict 設定で検査する。
  結果のテストは IO・解析・描画に分け、公開描画関数の戻り値に図と axes の型を付けた。

### Fixed

- `make build-kernel` が PCG64 の C++ 部品に `-fPIC` を渡しておらず、共有ライブラリをリンクできなかった。
  C・C++ にも位置独立コードを指定し、C++ runtime をリンクする。
- 参照点の iterator を保存値の選択と計算で二重に消費していた問題を修正した。モジュール関数と `Beach` で
  `source="saved"` の引数制約もそろえた。
- `FieldKernel` は入力電荷をコピーして保持する。呼び出し側が元の配列を書き換えても、借用 kernel の復元状態は変わらない。
  `Beach` の kernel cache は、解決済みの設定と同じ電荷を持つ呼び出しで共有し、委譲メソッドは指定した履歴 step を使う。
- モジュール関数に `Beach(..., config_path=...)` を渡したとき、`analyze_coulomb_mobility` と `plot_coulomb_force_matrix` の
  内部の Coulomb 計算、`plot_mesh_source_boxplot(quantity="potential")`、`BeachScene.from_result` が明示した `beach.toml` を
  使わず近傍の自動探索に戻っていた。内部の呼び出しは run の `RunContext` を引き継ぐ。
- `python -m beach.cli.main` は何も実行せずに終了していた。`beachx` と同じ CLI を実行する。
- `animate_history_mesh` と `Beach.animate_mesh`（`beachx animate`）は `view_elev` / `view_azim` を省略すると matplotlib の
  既定（仰角 30°、方位角 -60°）で描いていた。docstring と CLI help のとおり、他の mesh 描画と同じ 24° / -58° を使う。
  描画の既定の view 角度と colormap は `beach/plot/_defaults.py` の定数だけが持ち、関数・`Beach`・CLI の
  既定値と help はそこから取る。

## [3.8.0] - 2026-10-23

### Changed

- 乱数をすべて PCG64（[pcg-cpp-fortran-bind](https://github.com/Nkzono99/pcg-cpp-fortran-bind)、DRIFT・MPIEMSES3D と同じ
  revision）に切り替えた。粒子注入・境界イベント・二次電子は (seed, batch, slot, ...) を Murmur3 で 128 bit の鍵にした
  counter-based 乱数、単体の sampling API は seed と thread から開く逐次 PCG64 stream を使う。Fortran 組み込みの
  `random_number` と自前の 31 bit ハッシュ乱数は使わない。同じ seed でも 3.7.x とは粒子列が変わる。
- バッチの粒子は種・source・面ごとの global slot に並べて生成する。既定の `run.batch.particle_distribution = "static"` は
  slot を rank に連続区間で割り振り、粒子ごとの境界・二次電子乱数も global slot で決まるので、同じ seed ならバッチの
  粒子集合は rank 数に依らない。
- checkpoint schema v11 は乱数状態ファイル（`rng_state.txt` / `rng_state_rankNNNNN.txt`）を書かず、再開時に MPI
  rank 数の一致を要求しない。v10 以前の checkpoint はこれまでどおり読み、乱数状態ファイルは存在だけを確認する。
- MPI 版は `MPI_Init_thread(MPI_THREAD_SERIALIZED)` で初期化する（`shared` の claim を OpenMP critical 区間から呼ぶため）。
- build に C++ コンパイラが必要になった（fpm が `icpx` / `g++` などを使う）。

### Added

- `particles.tracking.step_growth`（既定 `false`）、`step_growth_max`（既定 5）、`step_growth_energy_fraction`（既定 0.01）。
  場の変化が小さい領域で粒子の step を `dt_s` の 2^n 倍まで伸ばす。受理した step の中点電場の変化が
  `substep_field_change_max` の 4 分の 1 未満（または変化 × 変位の仕事が 1e-4 eV と運動エネルギーの
  `step_growth_energy_fraction` 倍の大きい方以下）なら次の step を 2 倍にし、大きければ既存の半分化で戻す。
  `substep_field_change_max` が 0 なら 0.05 を自動で使う。本番相当の入力では step の 4 分の 3 が粒子層から 20 µm 以上
  離れた弱い場にあり、粒子追跡が 4.4 倍速くなった（最長の捕捉光電子は 20 万 push から 3 千 push、455 ms から 43 ms）。
  吸収・脱出数と表面電荷は実行間の揺らぎの範囲で一致した。`summary.txt` に `adaptive_grown_steps=` と
  `adaptive_step_growth_max=` を記録する。
- `run.batch.particle_distribution`（既定 `"static"`）と `run.batch.particle_pool_chunk`（既定 64）。`"shared"` で全 rank
  共有の slot pool を rank 0 のカウンタ（MPI RMA fetch-and-add、次の chunk を先行要求して粒子追跡と重ねる）から
  chunk ずつ取り合い、rank 内では thread が 1 slot ずつ受け取る。rank 間の粒子追跡負荷を動的に均す。`summary.txt` に
  `particle_distribution=`、`particle_pool_claims=`、`particle_pool_claim_seconds=` を記録する。
- `fields.solver.near_field_precision`（既定 `"double"`）。`"single"` で FMM 近接場の葉ごとの panel 表を葉中心基準の
  単精度で持ち、8 レーン SIMD の核で和を取る。本番相当の入力でその和が約 1.5 倍速く、場の相対誤差は 10⁻⁶ 程度で
  展開次数 4 の FMM 打ち切り（2×10⁻⁴）の 100 分の 1。平面上の点は倍精度の核に落ちる。`summary.txt` に
  `field_near_field_precision=` を記録する。
- `bem_counter_rng`（PCG64 counter-based の `counter_uniform`、`counter_key`、`counter_child_id`）と逐次 PCG64 stream の
  `bem_random_stream`（`seed_random_stream`、`random_uniform`）。
- `bem_injection` に counter-based の単一粒子生成 `sample_reservoir_face_particle_counter`、
  `sample_photo_ray_slot`、`cast_photo_ray`、`photo_ray_macro_weight` を追加し、`bem_app_config` に
  バッチ注入 layout `batch_injection_layout_type` / `build_batch_injection_layout` / `generate_batch_slot` を追加。

## [3.7.1] - 2026-10-22

### Changed

- CLI `beach` の起動時に `OMP_SCHEDULE` が未指定なら `dynamic,1` を選び、通常の粒子追跡を 1 粒子ずつ
  空いた OpenMP スレッドへ割り当てる。明示した `OMP_SCHEDULE` は優先し、適応的な $k\ne0$ 進行では
  再試行の再現性のために従来どおり粒子追跡中だけ `static` を使う。

## [3.7.0] - 2026-10-22

### Fixed

- camphor のビルドプロファイル（`env/camphor.env`）で作った MPI 版が、どの設定でもバッチループの先頭で
  `forrtl: severe (154): array index out of bounds` で停止していた。ifx 2023.2.4 が `-xHost`（Xeon Max の AVX-512）で
  8 byte 境界の stack 一時領域へ `vmovaps %xmm16` を出す誤コンパイルで、少なくとも 3.5.0 から起きていた。
  プロファイルの命令セットを `-xCORE-AVX2` にした。ifx 2024.0 と 2025.0.4 でも同じ非整列ストアが出て停止する
  （2026-10-22 に計測）ので、これらでも `-xHost` は使わない。

### Changed

- `beach --check-config` が、これまで実行中のバッチで初めて停止していた設定を preflight で拒否する: `run.batch_count <= 0`、
  負の `run.batch.duration_s`、photo_raycast 種の finite `[domain]` なし・`run.batch.duration_s <= 0`・質量 <= 0・温度 < 0。
  エラーの文言は変わらない。
- FMM の近接場で、像シフトした panel を評価するたびに 328 byte の geometry を複製して一般の核へ渡していたのをやめ、
  保存済み geometry をその場で読み、頂点だけ平行移動し、頂点距離を辺ごとに二度ではなく一度だけ求める核
  （`panel_potential_field_translated`）を追加した。葉の近接 panel の総和は核モジュール内の 1 つのループ
  （`panel_field_sum_translated` / `panel_potential_sum_translated`）で行う。演算順は同じなので電場・電位は
  bit 単位で一致し、テストで逐次評価と比較している。粒子ステップの大半を占める近接和が速くなる。

- FMM の各葉の近接 panel を、像シフトを掛けた成分配列（頂点 x/y/z、法線、辺長、辺の外向き法線、面積）として
  plan に保存し（`near_vx` など、`build_near_panel_tables`）、点評価は葉の区間をそのまま読む。近接和は 64 panel ずつ
  SIMD ループで評価し、合算は従来どおりリスト順に行う。値の演算順は変えていないので結果は bit 単位で同じ。
  本番相当の入力（3,392 panel、4 MPI×16 スレッド）で 1 batch が 1.33 s から 0.69 s になった。
- `update_state` と `snapshot%refresh`、`field_solver%refresh` に省略可能な `mpi` 引数を追加した。batch ループが渡す
  MPI コンテキストで M2L の target node を組数がほぼ等しい連続ブロックに分け、各 rank が自分のブロックだけ計算して
  局所展開を allgather で揃える（`bem_mpi` に in-place の allgatherv 補助を追加）。各 node の展開は従来と同じループで
  1 つの rank が求めるので、全 rank で複製していた結果と bit 単位で同じ。`test_mpi_hybrid` で確認する。

- opt-in の `env/camphor-fast.env`（`make install-camphor-fast`）を追加した。`-fp-model fast=1` で近接核の `log`/`atan2` が
  4 レーンの SVML になり、20 バッチの測定（4 MPI rank × 16 スレッド）で 1 batch が約 1.3 倍速い（0.62 s → 0.48 s）。
  短い単一スレッド比較では電荷・履歴に相対 10⁻¹³ 程度の差があった。粒子の衝突・脱出結果も変わり得るため、
  この値を他の入力や長時間実行の誤差上限とはしない。テストと golden 比較は従来の `camphor.env` で行う。
  fpm は `--flag` の後に自分のプロファイル
  フラグを足して `-fp-model` を上書きするため、`build.sh` に `FPM_NO_PROFILE=1`（`FFLAGS` だけを fpm に渡す）を追加した。
- SIMD の近接和で、電荷ゼロの panel と面上の評価点を非特異な経路へ分けた。使わないレーンでの特異式の評価による
  浮動小数点例外を避ける。ビルドのプロファイル選択には macOS 標準 Bash にも対応する構文を使う。

## [3.6.1] - 2026-10-22

### Fixed

- periodic2 の zero mode（k=0）で、傾斜 panel の累積電荷を絶対 z の二次・三次の冪で積分していたため、薄い高さ分布を
  原点から離すと有効桁を失い、電位が符号ごと狂っていた（高さ差 1 µm の三角形を z=1 m に置くと −2.26 V が +122 V）。
  区間ごとの累積電荷を区間左端を原点とする局所座標で持ち、panel の寄与を panel 自身の高さ幅だけで表すようにした。
  平行移動しても一次モーメントと Gauss 則が保たれる（#41）。
- 浮遊導体の LU 分解は L の既存列ごと行を交換するのに、解く側は右辺の行交換と前進代入を交互に行っていたため、途中で
  行交換が起きると電荷再配分後も同じ導体内に電位差が残っていた。右辺の置換を全て先に適用する（#46）。
- `run.batch.adaptive.max_nonzero_mode_potential_step_v` を書くと、値が無効を意味する 0 でも、明示した
  `fields.periodic.backend="finite_images"`（または既定の finite_images）が `panel_spectral_reference` / `exclude_k0` /
  `e_bottom_zero` に黙って切り替わっていた（Fortran parser と Python の正規化の両方）。adaptive 上限は backend を選ばない。
  0 は何もせず、正の値は従来どおり `cached_kneq0` 以外では設定エラーになる（#44）。
- PE と SEE があり外部根が Type C（放出電子の帰還 0）のとき、SEE が脱出したバッチで固定電流 closure が存在しない PE 帰還へ
  目標を割り当てて停止していた（#48）。Type C では PE が全て脱出して帰還がないので、PE の帰還目標を 0、脱出目標を PE 放出とし、
  SEE を PE なしの構成と同じく根へ固定する（帰還・脱出を根の値、放出の反作用電荷を根に与えた外側源へ）。SEE の取り分は
  根の脱出 − PE 放出で、最初の更新前は 0。Type A/B では従来どおり PE の帰還が SEE 脱出を差し引く。

## [3.6.0] - 2026-10-22

### Changed

- `sheath.coupling.include_reflected_electrons` の既定を `true` にした。外部根を更新する run（`outflow_refresh_batches>0`）では、
  セル内で反射して z-high を抜けた上流電子を既定で外部源に数える（sheath-model の結合ガイドどおり）。外部根を固定する run では
  効かない。Type C の更新 run では φ0 が +0.06〜0.1 V 動き、3.5.0 までの結果や checkpoint からの続きと一致しなくなる。
  以前の挙動は `include_reflected_electrons = false` で得られる。明示的に `true` と書いた場合だけ、更新間隔 0 を設定エラーにする。
- `sheath.closure="zero_current"` で `secondary_emission` 種が 1 つだけなら、`sheath.species.secondary_electron` を省いても
  その種を SEE の役に当てる。以前は登録しない SEE 種が表面収支に入らず、セル総電荷が random walk していた。複数あって
  指定がなければ設定エラーにする。

### Added

- `sheath.coupling.include_reflected_electrons=true` を PE も SEE もない構成でも使えるようにした。セル内で反射して z-high を
  抜けた上流電子だけが外部根の源になり、源 0 の Type C から最初の更新で入る。根の脱出・帰還は全て反射電子のものなので、
  電子の吸収目標を根の「流入 − 脱出」、脱出目標を根の脱出に固定し、追跡した脱出の揺らぎでセル総電荷が動かないようにする。
  この場合に限り、放出電子の種がなくても `outflow_refresh_batches>0` を受け付ける。

### Fixed

- `build.sh` は fpm が失敗・中断した次の実行で、その間に更新された対象を作り直して依存側へ伝える。fpm 0.12 は同じ実行で
  依存先が再ビルドされない限り digest の一致した対象を飛ばすため、型を変えた後の失敗ビルドで依存側が古い `.mod` と
  派生型レイアウトのまま残り、無関係な箇所で落ちていた。`build.sh` を通さない直接の `fpm` 実行は対象外。

## [3.5.0] - 2026-10-20

### Added

- `sheath.coupling.include_reflected_electrons`（既定 `false`）: セル内で反射して z-high を外向きに抜けた上流電子を
  外部根の放出源に含める。sheath-model の結合ガイドが求める扱いで、各バッチで追跡したその脱出を電子の吸収目標に足し戻し、
  PE（なければ SEE）の脱出・帰還（SEE だけなら放出も）目標から差し引く（closure contract の `reflected_electron_species`、
  `absorbed_target_credits_tracked_escape`）。窓の反射束は `electron_outward_flux_m2_s` に出し、更新ログに
  `reflected_flux_m2_s` と `reflected_mean_normal_energy_eV` を加える。`outflow_refresh_batches>0` が必要。

### Changed

- 上端の流出スペクトルの bin を sheath-model の既定格子（`log_energy_edges` / `log_energy_bin`、145 bin）から取り、
  BEACH 側の bin 定義を削除。束と Maxwell 規格化密度の換算も sheath-model の `emitted_fluxes` / `with_outward_flux` を使う。
  sheath-model の依存を plane coupling helpers（PR #10 の merge commit 2455f44）へ更新。

## [3.4.0] - 2026-10-19

### Changed

- `zhao_stationary` で PE がなく SEE だけの構成でも表面収支を毎バッチ閉じる。根の脱出束・帰還束は全て SEE のものなので、
  PE と同じく SEE の帰還・脱出を根の値に固定し、放出の反作用電荷は根に与えた外側源（追跡した H 外向き束の窓平均）に
  揃える（closure contract の `emission_target_from_tracked_source`）。SEE の放出個数は追跡のまま決まり、根は予測しない。
  最初の更新前の源 0 の根では反作用電荷を置かない。生の帰還がないバッチでは帰還目標を放出へ畳み、SEE が放出されなかった
  バッチの根の脱出分だけを `secondary_budget_residual_a_m2` に残す。3.3.1 までは SEE だけの構成でセル総電荷が SEE 脱出の
  雑音で random walk し、最初のバッチの SEE 脱出がそのまま板の電位に残っていた。

## [3.3.1] - 2026-10-19

### Changed

- `zhao_stationary` で SEE を使うとき、PE の脱出・帰還目標を「外部根の脱出束 − そのバッチで追跡した SEE の脱出電流」にした
  （closure contract の `escape_partner_species`）。外部根の脱出束は PE+SEE 合算の値なので、追跡した SEE の分を差し引けば
  表面の電流収支が毎バッチ厳密に閉じ、セルの総電荷は PE だけのときと同じく固定される。従来は PE の脱出目標を窓の
  H 外向き束の割合で配分していたため、SEE の脱出の統計雑音がセル総電荷を random walk させ（例題の 0.1×0.1×1 mm セルでは
  電子 1 個が板と H の電位差 1.8 mV に相当し、40 バッチで数 V）、系統的なずれも残っていた。`secondary_share` は
  そのバッチで SEE に充てた割合、`secondary_budget_residual_a_m2` は SEE の脱出が根の脱出束を超えて差し引けなかった分
  （通常 0）になる。履歴・summary・checkpoint の列は変えない。PE なしの SEE だけの構成は従来どおり差し引く先がない。
- sheath-model を `db32ee7` へ固定する。文書の再構成と、全入力の背景電子ドリフトの既定値を 0 に揃えた変更を含む。
  BEACH は外部根の入力を設定から全て与えるので、計算結果は変わらない。原稿の参照ファイル一覧も新しい文書名に更新した。
- 3.3.0 はリリース作業の手違いで GitHub 上のタグ名が使用済みになったため欠番。公開した内容はこの版と同じ。

## [3.2.0] - 2026-10-19

### Added

- `particles.tracking.substep_field_change_max` / `substep_max_halvings`: 粒子ごとの適応サブステップ。step 候補の中点電場が
  前回受理した中点電場から指定の割合より大きく変わる step を、上限回数まで半分にして取り直す。帯電した要素の近傍で
  1 step の変位が場の変化スケールに近いときにエネルギー保存を保つ。既定は無効。半分にした step 数と sub-step 総数、
  上限に達した step 数を `summary.txt` に出力する。三層レゴリス（1.98 ps、球が +290〜+700 V）の 100 バッチ継続では、
  H を反射して出る電子の全エネルギー差の rms が dt のまま 73 eV（A）/ 545 eV（C）だったものが、閾値 0.2 で
  0.56 / 2.7 eV、0.1 で 0.22 / 0.69 eV（dt/4 の一様刻みは 0.19 / 1.2 eV）になる。光電子が帯電した球の表面から出る
  この case では大半の粒子が強い場を通るので、計算時間は閾値 0.1 で dt/4 と同程度である。
- runops adapter: `runo new --from RUN --continue` resumes from the source run's checkpoint. The new run links
  `work/restart` to the source output, sets `[run.restart].from`, and records the source in `input/restart_source.txt`.
  `run.batch_count` must exceed the saved `batches`; otherwise `runo new` stops with the number to use.
- runops adapter: `runo new` / `runo sweep` warn about `--set` keys that `beach.schema.json` does not define.

### Fixed

- runops adapter: a run counts as completed only when `summary.txt` reports `batches` equal to `run.batch_count`, so
  a run stopped at the walltime after a periodic checkpoint is `failed`. `runo tidy --purge-work` keeps the final
  checkpoint in `work/latest/` (RNG state, residuals, ledger, `checkpoint_complete.txt`), so a purged run can still be
  continued.

## [3.1.0] - 2026-10-18

### Added

- periodic2 cached 非零モードの上部真空 Fourier 切替。高さの幾何自動選択と次数制御、接続面の電位・法線場 rms、再構築 receipt を追加。
- Zhao 外部根の法線電場を零モードへ渡す `e_top_from_sheath` と、シース側・セル側の上端電場履歴を追加。checkpoint schema は維持。
- 17/19列 matching-plane 履歴と8列 top-reference 履歴を、新しい列へ拡張して再開する互換処理。
- `beach --prepare-periodic-cache beach.toml`: periodic2 の遠方補正 operator を作成・publish して終了する CLI。本番より多い MPI rank で cold build だけを先に済ませられる。

### Changed

- periodic2 `cached_kneq0` の線源セルと評価側 anchor を、軸ごとの分割数を持つ box 格子にした（生成器 v12）。
  どのセルもほぼ立方体になるまで最も幅の広い軸だけを分け、箱の高さが周期長より大きくても横方向の分割は増えない。
  `examples/zero_current_refresh.toml`（0.1×0.1×1 mm）の演算子は 4096 anchor × 4096 セル（5 GB、数時間）から
  60 anchor × 120 セルになり、三層レゴリスも 64 anchor × 512 セルから 20 × 40 になる。anchor は適応木の節点ではなく
  点評価で参照する。v11 の fingerprint は一致せず、演算子を再生成する。
- periodic2 `cached_kneq0` の生成器を v11 に更新。線源を box 固定セルへ分け、相対オフセットと鏡映対称性で Ewald 評価を共有する。
  同じ box と設定なら異なるメッシュで演算子を再利用する。v10 の fingerprint は一致せず、必要な演算子を自動再生成する。
- `cached_kneq0` の演算子は要素を含む線源セルの分だけ作り、cache にないセルだけを追加生成して同じファイルへ書き足す。小さなテスト mesh や上部が空の box で生成が短くなる。
- FMM の M2L 微分表を多項式の反復乗算から漸化式に置き換え、三層レゴリスの起動（cache 読込済み）が 28 スレッドで 98 s から 5.5 s になった。
- `fields.periodic.cache_dir` の未指定時は `$BEACH_CACHE_DIR/periodic2`、`$XDG_CACHE_HOME/beach/periodic2`、
  `$HOME/.cache/beach/periodic2` の順で保存先を解決する。明示指定した保存先は維持する。
- `zhao_stationary` の外部根、根が z-high 面に課す条件（壁電位・上端電場、access/障壁電位、上流貯留部、各 channel の流束）、
  外部根の更新間で引き継ぐ根を sheath-model の結合 API（`sheath_matching_iteration`、`evaluate_matching_conditions`、
  `observed_outflow_source`）から得る。粒子種の対応、参照面積、固定電流 target、PE/SEE の分割、summary による再開は BEACH が持つ。
  sheath-model を放出源 API を改名した `384985f` へ固定する。

### Fixed

- runops adapter: `runo case new` writes the grouped `beach.toml` that `beachx config init` produces (the template was
  still in the pre-3.0 flat layout, which BEACH 3 rejects). A flat `input/beach.toml` now fails at `runo new` with a
  clear message instead of in the queue. OBJ meshes in `[[mesh.parts]]` are referenced from the run directory: a path
  shipped with the case becomes `input/<path>`, and an absolute path is copied into `input/`.
- `import beach` no longer imports numpy and matplotlib. The public names from `beach.fortran_results` load on first
  access, so `beach.config` and the runops adapter work inside runops's own environment.

### Known issues

- 数百 V に帯電した要素の近傍では、既定の時間刻みの指針で選んだ `particles.tracking.dt_s` で電子のエネルギーが保存しない。
  三層レゴリス（球が +290〜+700 V、1.98 ps、1 ステップ約 20 µm）では H を反射して出る電子の約 2 割が 10 eV 以上を
  数値的に獲得し、H での PE 平均法線エネルギーと更新後の壁電位（−0.19 V）も変わる。dt を 1/4 にすると保存する
  （エネルギー差 rms 0.2〜1 eV）。帯電した表面の近傍を扱う case では dt を 1/4 程度にして収束を確認する。

## [3.0.0] - 2026-10-15

### Changed

- Documentation ([ADR 0003](docs/adr/0003-documentation-restructure.md)): the LaTeX manuscript in `docs/manuscript/`
  is the single source for the models and numerical methods. The site publishes it as numbered chapter pages and a PDF
  (`tools/manuscript_to_starlight.py`, run by `make docs-pages`), and the manuscript now covers secondary electron emission,
  the periodic2 floating conductor, time-step and batch-width scales, and the SEE-aware outer-root refresh.
  Japanese is the only documentation language for now, and the pages under `docs/` are grouped into `start/`, `guide/`,
  `models/`, `reference/`, and `developer/` with relative links.
- The official tutorial is a 5 eV electron beam (2.2 uA/cm^2) charging a floating 20 mm plate for 1.5 us
  (`examples/tutorial_electron_beam.toml`, also generated by `beachx config init`). Its page compares the saturation
  potential, capacitance, and charging time with hand estimates using `examples/tutorial_beam_analysis.py`, and
  contrasts an insulating plate.

- 照射セルの壁・周期軸・下端モデルを実行時契約へ集約し、物理結果と保存出力名を保つ。

- Mesh input now uses ordered `[[mesh.parts]]` entries and mixes OBJ files with procedural templates. Each enabled part receives one `mesh_id`, surface policy, and conductor group. Removed `mesh.mode`, `[mesh.obj]`, and `[[mesh.templates]]` without input aliases. Added mixed geometry, preflight contracts, and synchronized examples/docs.

- Fortran の設定エラー・警告を grouped 公開キー名と公開選択値へ統一し、粒子種の添字付き診断と reader 名の混入を検査する回帰テストを追加。

- 入力パラメータリファレンスと日英の利用者向け文書・plugin の設定例を grouped 7 グループの公開キー名へ統一。公開 schema から文書のキー網羅と見出し階層を検査する。

- 公開 TOML 入力を `run / domain / mesh / particles / fields / sheath / output` の grouped 7 グループ形式に統一。Fortran と Python は flat 入力を明確なエラーで拒否する。
- grouped examples を `examples/` 直下へ移動。direct/FMM、入門・単一電子・形状見本の参照構成も grouped に変換して維持。

### Removed

- `examples/tutorial_insulator.toml` (replaced by `examples/tutorial_electron_beam.toml`) and the English `docs/*.en.md` pages.
- `beachx config migrate`、flat 入力の公開 schema (`legacyConfig`)、旧形式 fixture、重複した flat examples。
- `docs/GroupedConfiguration*.md` と plugin 内の移行ページ、ナビゲーションと移行リンク。

### Added
- Secondary electron emission after absorption: Whipple / constant yields, charge-scaled Bernoulli macro counts,
  inherited weights, shared PE velocity sampler and reaction-charge channel, and deterministic per-thread spawn buffers
  tracked generation by generation in the same frozen batch field (`particles.secondary_emission.max_generations`, default 3).
- Optional `sheath.species.secondary_electron`: shared kinetic barrier and PE + SEE outflow spectrum, including SEE-only
  refresh. Provisional Phase 1 allocates only the PE escape target by its flux share; SEE remains fully tracked.
- `secondary_emission_history.csv`, cumulative SEE summary diagnostics, matching-plane `secondary_share` and
  `secondary_budget_residual_a_m2`, paired SEE documentation, two examples, and focused contracts. Released 17-column
  matching-plane histories remain readable and are extended atomically on resume.
- Added [docs/adr/0002-irradiation-cell-scope.md](docs/adr/0002-irradiation-cell-scope.md): BEACH supports one
  product configuration, the `periodic2` irradiation cell with the outer 1-D sheath, and keeps the free-space direct / FMM
  cases and the tutorial as validation references. The record lists the features removed for that scope and the phases
  (one validation source, feature removal, grouped-only input, `[[mesh.parts]]`, and collecting the cell assumptions).
  `SPEC.md` 2.1 states the supported configuration.
- The Fortran preflight now owns every configuration combination that the Python validators used to check: exclusive
  `batch_duration` / `batch_duration_step`, exclusive `e0` / `e0_abs`, `[particle_boundary]` and species boundary overrides
  require a `[domain]` and cannot touch a periodic face, boundary inflow needs an open face and a `volume_seed` species,
  unique `species_key`, nonzero `q_particle`, exactly one density form, the `fixed_current` target signs and
  `neutral_return` prerequisites, tree and `cached_kneq0` parameter ranges when their solver is active,
  the adaptive nonzero-mode step prerequisites, `mesh.surface_side` for OBJ meshes, template hole radii, and the
  photoelectron settings when `photoelectron_source_scale = 0`.
- Added `surface_current_model.zhao_upstream_field_floor_v_m` (grouped `sheath.zhao.upstream_field_floor_v_m`, default 0).
  It accepts drifting Zhao A/C roots when the field deficit $\sqrt{-\min E^2}$ inside the upstream negative-$E^2$ band
  does not exceed this external field level, which places the upstream boundary where the sheath field has fallen to
  the level the solar wind itself carries (the motional field $|\mathbf{u}\times\mathbf{B}|$, about 2 mV/m). The deficit
  is written as `surface_current_model_upstream_field_deficit_V_m`. `zhao_upstream_band_tolerance` remains available;
  when both are set a root must satisfy both. The sheath-model dependency is pinned to `b4a2521`, which adds `upstream_field_floor_v_m`.
- Added `surface_current_model.zhao_ion_pressure_factor` (grouped `sheath.zhao.ion_pressure_factor`, default 1). The
  `zhao_stationary` outer root now passes the ion species temperature and this pressure coefficient to the warm-ion fluid of
  sheath-model (Sana & Mishra 2026) instead of a cold beam. Preflight replaces the `T_i <= 0.1 T_e` stop with the supersonic
  entry condition `m_i u_z^2/e > factor*T_i`. For solar-wind entries (about 1 keV) the root moves by less than 0.1 mV and the
  upstream negative band is unchanged; each root solve costs 40–400 ms.
- `zhao_stationary` preflight warns when the photo_raycast surface emission (`emit_current_density_a_m2` times the sine of
  the ray elevation) differs by more than 5% from the outer photoelectron source, or when the ray elevation differs from
  `solar_elevation_deg` by more than 1 degree. The fixed-current scale absorbed such mismatches silently.

### Changed
- `surface_model="conductor"` now works with `field_boundary.mode="periodic2"`. The floating-conductor solver builds its
  influence matrix from the field snapshot (periodic images and zero mode included) once at initialization, keeps the LU
  factorization, and takes the fixed potential from the snapshot each batch; the adaptive nonzero-mode step judges the
  redistributed candidate charge. The `commit_batch_charge` / `prepare_adaptive_charge_candidate` interfaces take the
  snapshot and the `floating_conductor_solver_type` instead of `external_e` / `field_bc_mode`, and
  `apply_surface_model_charge_relaxation` is removed.
- `beachx lint` runs the packaged JSON Schema and then `beach --check-config`, and prints the Fortran diagnostic when the
  preflight rejects the input (`checks=toml,schema,preflight`). The `beach` executable is found through `--beach`,
  `BEACH_CONFIG_CHECK_EXE`, then `PATH`; lint fails without one. The Python semantic validators
  (`beach/config/_*_validation.py`), `validate_runtime_config`, and `beachx config validate` are removed, so the schema
  owns keys, types, enumerations, and ranges and the Fortran preflight owns combinations and derived values.
  `tests/python/test_config_contract.py` checks that lint returns the preflight verdict for each representative case.
- `surface_current_model.model="zhao_stationary"` (grouped `sheath.closure="zero_current"`) now requires a periodic2 box with
  x/y periodic axes and an explicit split `[periodic2]` table (`cached_kneq0` or `panel_spectral_reference`, `exclude_k0`);
  other field settings are a configuration error instead of a startup warning. Without the split zero mode the top-face mean
  potential could not be fixed to the wall potential $\phi_0$, so the inflow map and the outer barrier used a top face near
  0 V. In an 8x8 flat plate the mean plate potential stayed near 0 V instead of $\phi_0$, the Type C electron scale was 0.53,
  and the Type A photoelectron return and escape scales were 1.99 and 0.30. With the split zero mode the plate stayed within
  0.1 mV of $\phi_0$ and all scales were 1 within the statistical error. The fixed-root Zhao examples
  (`periodic2_zhao_fixed_current.toml`, `periodic2_zhao_no_photo_fixed_current.toml`, and their grouped versions) now use
  `cached_kneq0`, and the stale comment about a 10% band tolerance in the Zhao examples now describes the field floor.
- Restructured [Connecting to the outer sheath](docs/models/ZhaoStationaryClosure.md) around the four conditions the root gives the
  top face, and moved the electron-drift discussion to a new page, [Electron drift and the upstream band](docs/models/ZhaoUpstreamBand.md):
  why the band appears, what happens to the potential profile, the two acceptance criteria with measured deficits, and the
  physics that fills the band in a real plasma.
- The Zhao examples emit photoelectrons at the normal-incidence value of the outer source (2.5446e-6 A/m²) along a ray at the
  configured 60-degree solar elevation, so the fixed-current emission scale is 1 instead of about 0.02. They now accept the
  drifting root with `zhao_upstream_field_floor_v_m = 2.0e-3` instead of `zhao_upstream_band_tolerance = 0.1`.
- Documented that the outer-root refresh retries branches in A → B → C order after the previous branch; the
  solar-elevation ordering applies only to the initial root.

### Removed
- Removed `particle_boundary.ordinary_open_model` (grouped `particles.boundary.open_model`) and the ordinary-open
  `potential_barrier` model, including its `particle_ordinary_open_model` summary receipt. Ordinary open faces now
  always escape. Reservoir `source_vdf` / `infinity_barrier` inflow, `phi_infty`, `face_potential_grid_n`, and the
  Zhao `zhao_barrier_v1` per-species kinetic barriers remain available.
- Removed particle source modes `reservoir_face` and `plane_source`, along with `source_normal` (grouped
  `particles.species.source.normal`). Use boundary reservoir inflow for external-plasma supply; partial boundary
  apertures and internal plane sources are no longer supported. Removed inputs are rejected without aliases.
- Removed the `velocity_distribution = "grid"` CSV velocity distribution: `velocity_grid_path`, `velocity_grid_pdf_kind`,
  `velocity_grid_sampling`, `particle_flux_m2_s`, `current_density_a_m2` (grouped `distribution.model`, `distribution.grid_*`,
  `distribution.particle_flux_m2_s`, `distribution.current_density_a_m2`, `sampling.velocity_grid_sampling`),
  `bem_injection_velocity_grid`, `test_external_field_velocity_grid`, and `examples/velocity_grid.csv`. Flux-driven sources
  use the drifting Maxwellian from the species density, temperature, and drift velocity.

### Fixed
- `beachx slices` reads the current 7-group `beach.toml`; it looked for the legacy `[sim]` table and always failed.
- `bem_surface_current_model` no longer uses `findloc`, which GCC 8 gfortran rejects, so the gfortran build of the outer-root
  spectrum source works again.

## [2.0.0] - 2026-10-13

### Added
- Added `surface_current_model.zhao_upstream_band_tolerance` (grouped `sheath.zhao.upstream_band_tolerance`). With
  inward electron drift, Zhao Type A/C roots have negative $E^2$ in a narrow band next to the upstream potential and were
  always rejected; a positive tolerance accepts them when the band is no wider than that fraction of $|\phi_m|$ (A) or
  $|\phi_0|$ (C), and records the width as `surface_current_model_upstream_negative_band_V`. The default 0 keeps the exact
  condition. Smoothing only the slow upstream electrons (the kinetically unstable cusp of the mirrored drifting
  distribution) removes the band and moves the wall potentials by tens of mV, far less than neglecting the drift
  (about 1 V for Type A, 2.7 V for Type C in the documented case).
- Added `surface_current_model.outflow_refresh_batches` (grouped `sheath.coupling.outflow_refresh_batches`) for
  `zhao_stationary`. Every N accepted batches it re-solves the outer zero-current root with the PE outflow observed at
  z-high as the outer emission source, while the surface emission target and the per-batch floating targets stay
  unchanged. The outflow is passed as its normal-kinetic-energy spectrum (24 logarithmic bins per decade from 1 meV to
  1 keV) scaled to the window transmission, not reduced to a half-Maxwellian with the same mean energy: the outflow can
  be a body near the emission temperature plus a tail of electrons accelerated in the cell, the tail sets the mean, and
  the mean-energy half-Maxwellian then miscounts escaping and slow electrons and can lose the outer root. The accepted
  outer state is written to
  the `matching_plane_*` summary receipts (the spectrum as `matching_plane_photoelectron_spectrum_m2_s`),
  `matching_plane_history.csv`, and `matching_plane_spectrum_history.csv`, and restarts reconstruct it from the
  checkpoint. A checkpoint without the spectrum key resumes with the half-Maxwellian source until the next refresh.
- Added `fixed_current_history.csv`, written every `output.history_stride` batches when a fixed-current closure is
  active. Each row holds the tracked and target signed charge of one species channel (`absorbed`, `emission`, `escape`)
  summed over the write interval, so the closure's scale factor can be followed in time instead of only as the
  run-cumulative value in `charge_ledger.csv`.
- Added a pluggable `[surface_current_model]` closure with a restored Zhao A/B/C stationary zero-current solver. It resolves
  ambient-electron, ion, PE-emission, PE-escape, and PE-return targets independently and applies them through per-species
  `fixed_current` spatial maps.
- Added a config-driven field-kernel runtime benchmark that separates mesh construction, solver initialization, charge refresh, and volume/near-panel P0 field and potential evaluation across solver configurations.

### Fixed
- The periodic2 `cached_kneq0` far correction leaked a plane-mean component into the nonzero-mode field that grew with the
  surface charge (in a charged three-layer regolith about 6 V at the top of the layers and 1.4 V at the top face; the
  offset also made fewer photoelectrons escape than the outer Zhao root). Operator generator v10 fixes it at unchanged
  per-update cost: the proxy multipole rows are normalized before the pseudoinverse QR (the absolute rank threshold
  discarded higher moments in micrometre SI geometry); the exact plane mean of the finite image shell is split off in
  closed form and subtracted from a per-source height table, so the operator fits only a mean-free residual; and the
  operator maps the multipoles of depth-2 source-tree nodes to depth-2 target anchors instead of a single root multipole.
  With one image layer and order 4 the leak drops to 0.03 V or less and the lateral potential error above the layers from
  0.5–1.4 V to 0.01–0.03 V. The first operator build takes longer (about 27 minutes on one 28-thread process for the
  regolith cell, split over MPI ranks) and is cached; older operator caches are rebuilt. The mean field above all sources
  now uses the exact total charge.
- `zhao_stationary` now fixes the z-high plane-mean potential to the outer wall potential $\phi_0$ when a split
  `[periodic2]` zero mode is available. Previously the absolute outer barrier $\phi_m$ was compared with a z-high
  potential tied to the bottom of the charged geometry, so tracked PE and electron returns and the inflow energy map used
  a barrier lower than the root by about $\phi_0$; the fixed-current scaling hid this in the totals. Results of PE
  `zhao_stationary` runs with a split zero mode change. Field backends without a split zero mode warn at startup.
- `zhao_stationary` outer-barrier returns in an x/y periodic cell now re-enter at a counter-based uniform position on the
  z-high plane instead of the crossing point, because the lateral drift in the outer sheath far exceeds the cell width.
  The normal kinetic energy is corrected by the local potential difference between the crossing and return positions
  so that total energy is conserved; a slow particle that cannot reach the drawn position returns at its crossing point.
- The Camphor `ifx` install profile no longer passes unsupported inline/report options or enables the unstable IPO link path, allowing release executables to link reliably with Intel 2023.2.
- Reflected and periodic particle events now place surviving particles a scale-aware distance inside the box, preventing zero-time boundary chatter caused by a subnormal one-ULP offset at zero-valued faces. The existing eight-event safety limit is unchanged.
- FMM mesh-centroid potential output no longer adds an area-equivalent point self term after the analytic triangle-panel self integral has already been evaluated.

### Changed
- `zhao_stationary` now solves its zero-current root with the external sheath-model library (fpm git dependency pinned
  to sheath-model `bebe58d`, which also builds with Intel ifx) instead of the in-tree Zhao core, which is removed.
  sheath-model accepts a root only when the potential profile has $E^2\ge0$ throughout, approaches a neutral zero-field
  upstream state, and does not block the ion flow. With an inward electron drift, Type A and C have negative $E^2$ in a
  narrow band next to the upstream potential, so by default roots that the in-tree core accepted are rejected (the
  previous examples move from Type A to Type B at 60 degrees and have no root at 10 degrees);
  `zhao_upstream_band_tolerance` accepts them. The examples give solar-wind electrons the ion drift with a tolerance of
  0.1; test fixtures that check the exact condition use a zero electron drift. The in-tree core's singular `u=0` Type A
  path is gone. Ions remain a cold beam.
- `zhao_stationary` now injects the z-high ambient-electron inflow at the upstream electron Maxwellian density
  $n_{e,\infty}$ solved by the Zhao root (the density that keeps infinity quasi-neutral after the wall absorbs fast
  electrons), and requires the electron and ion species to share the solar-wind `number_density_*`. Previously the
  electrons were injected at the configured density, so the injected flux fell short of the root's electron current
  (about 17 % in the three-layer regolith runs) and the fixed-current closure rescaled every absorbed electron by
  $n_{e,\infty}/n_i$. Totals are unchanged; the electron sample count and per-batch noise change.
- `zhao_stationary` now accepts a nondrifting ambient-electron reservoir (`drift_velocity` z component 0).
- The supported boundary workflow is now local reservoir injection plus closed
  photoelectron reflection/neutral-return closure. `[external_boundary]` is
  limited to local-source inflow and ordinary open-face behavior.
- The native field-kernel C ABI is now version 2 and accepts triangle vertices directly; the former point-source build signature is not retained.
- FMM expansion order must now be at least one in the core, C/Python field-kernel APIs, and CLI; order zero cannot represent the far/local electric-field expansion and is rejected before plan construction.
- The mesh-input default is now deterministic `mode="template"`. Explicit `mode="auto"` remains available but, like `mode="obj"`, requires an OBJ `surface_side` so either file-selection branch has a complete triangle-panel contract.

### Removed
- Removed the unreleased `matching_plane_quasistatic` coupling (table and `zhao_online` backends, `implicit_zero_mode`,
  `zhao_root_selection`, `coupling_*` keys, fixed-point replay, `beach-zhao-response`, `beach-zhao-atlas`, and the offline
  `beach-kinetic-*` oracle). For a cell much smaller than the Debye length the mean charge is a stiff, noise-dominated
  state; `zhao_stationary` with `outflow_refresh_batches` replaces it by enforcing the floating condition every batch
  and re-solving only the outer root. The grouped `sheath.stationary.*` keys moved to `sheath.photoelectrons.*`,
  `sheath.reference_area_m2`, and `sheath.coupling.outflow_refresh_batches`.
- **BREAKING**: Removed the self-consistent outer-plasma/sheath subsystem,
  including `[outer_plasma]` and `[coupling]`, kinetic/Zhao closures, interface
  particle transfer, delayed return queues, outer checkpoint state and
  diagnostics, and their examples and tests. Removed keys are rejected rather
  than mapped to the local reservoir model.
- **BREAKING**: Removed the public centroid `point` element source, `sim.softening`, and the `[field].element_kernel` selector. All surface models now use the implicit P0 triangle-panel discretization; retained legacy keys fail as unknown input rather than selecting a compatibility path.
- **BREAKING**: Removed the executable `m2l_root_oracle` periodic far-correction backend and its exact-Ewald out-of-box fallback. New configurations are rejected with guidance to use `cached_kneq0`; historical output metadata remains readable.
- **BREAKING**: Removed `photo_escape_model` and its `boltzmann_cutoff`
  reduced-emission closure. Photoelectrons now keep their full emitted weight;
  the closed workflow uses local reflection and neutral-return correction.
- **BREAKING**: Removed the `individual_return` photoelectron-closure setting.

## [1.6.2] - 2026-07-15

### Added
- Opt-in `BEACH_COLLISION_DIAGNOSTICS=1` output for collision-grid traversal failures.

### Fixed
- Collision-grid DDA now classifies subnormal periodic-seam displacements as stationary using a mesh-scale tolerance, preventing false `grid_stalled` failures caused by an infinite axis increment.

## [1.6.1] - 2026-07-15

### Changed
- Package metadata now uses an SPDX license expression and declares the license file through the current `pyproject.toml` fields.

### Fixed
- PyPI source installations now include Fortran benchmark sources referenced by `fpm.toml`.
- Pull-request and release workflows now validate every explicit fpm target in the built source distribution and build a wheel from that artifact before publication.

## [1.6.0] - 2026-07-14

### Added
- Free-space `triangle_p0` Treecode evaluation with exact analytic panel near interactions, all-vertex MAC radii, and monopole far interactions for both electric field and potential.
- A machine-readable output manifest and documentation-contract tests for configuration, output, and restart guidance.

### Changed
- Treecode now accelerates arbitrary-point and mesh-centroid potential evaluation in addition to electric-field evaluation.
- Field-solver documentation now provides a canonical compatibility matrix for solver, kernel, boundary condition, and periodic far-correction combinations.
- Parameter, schema, and agent guidance now document the supported Direct `periodic2` split-reference exception without repeating cross-reference prose.
- Output and restart documentation now describes conditional outer-plasma and photoelectron checkpoints, MPI-global residuals, and OBJ mesh-source paths consistently.
- Starlight development pages now identify main-branch documentation, expose last-updated and edit links, and receive content/type checks in CI.

### Fixed
- Repaired the malformed `[coupling]` parameter table and documented every top-level configuration table, including `coupling.update_mode`.

## [1.5.0] - 2026-07-11

### Added
- `sim.open_boundary_model="potential_barrier"` for open-boundary return/escape decisions based on normal kinetic energy and the local potential barrier.
- Direct continuous triangle-P0 panel potential/field kernels with explicit vacuum-side, self-potential, principal-value, and field-jump contracts.
- Panel-aware FMM topology, near correction, exact triangle multipole moments, and a matching public C/Python kernel path.
- Explicit periodic2 zero-mode projection and panel-spectral nonzero-mode reference evaluation for non-neutral surface charge.
- Linear Debye, nonlinear 1D kinetic, and unified linear-response outer-plasma models with fail-closed applicability and branch diagnostics.
- Conservative particle transfer across the outer interface, including infinity-to-interface reservoir mapping, individual photoelectron return accounting, and explicit 3D electrostatic outer orbits.
- A versioned cached infinite-periodic nonzero-mode operator with checksum validation, corruption recovery, MPI root-only generation, and warm-cache reuse.
- Checkpoint schema 3, which restores complete outer-plasma potential, electric-field, charge-density, residual, and current state while retaining schema-2 read-only migration.
- Physics release gates covering L3, far-correction, two-rank MPI, cache concurrency, memory budgets, and six required convergence categories.
- Bilingual installation, tutorial, validation, troubleshooting, and physical-redesign completion-audit documentation.

### Changed
- `beach-inspect` now uses only precomputed mesh potential for its ordinary summary; `--recompute-potential` explicitly requests the quadratic reconstruction path.
- Periodic far-correction documentation now consistently describes `none` as the default, `auto` as its compatibility alias, and `m2l_root_oracle` as an explicit diagnostic mode.
- `outer_plasma.interface_z` is an ownership/handoff plane in the unified model rather than a truncation of the electrostatic field domain.
- The Starlight documentation navigation now follows installation, first run, output inspection, and numerical-validation workflows.

### Fixed
- Field-solver submodules no longer re-import host-associated FMM types, restoring editable builds with GFortran 13.
- Boris particle candidates now use same-time trapezoidal position updates and sample the boundary-element electric field at the predicted midpoint. The existing `boris_push` signature is unchanged, but corrected trajectories differ from previous releases.
- Particle tracking now orders the first mesh hit against the first box-face event, advances one reflected/periodic remainder, and fails closed when one outer step would cross a second box face.
- Reservoir injection now rounds one global expected macro-particle count before distributing particles across MPI ranks, making count and residual sequences independent of MPI world size.
- MPI checkpoints now store one root-owned global reservoir residual while retaining rank-local RNG state; ambiguous legacy rank-local residual files are rejected on resume.
- The workload estimator now rounds reservoir counts globally before MPI rank splitting and reports both global and selected-rank totals.
- The production particle loop now keeps box-event state off the no-crossing hot path; a representative 800-triangle direct-field benchmark measured 1.13% median overhead versus the pre-event baseline.
- Output-directory creation now uses the POSIX filesystem binding without shell evaluation, including paths with literal shell metacharacters.
- Periodic collision queries now fail closed on image/range limits and propagate deterministic OpenMP/MPI failure context for both particle tracking and `photo_raycast` injection.
- Treecode now descends mixed-sign charge nodes instead of accepting an unreliable near-cancelled monopole approximation.
- Charge-history readers now reject incomplete dense batch snapshots before exposing the batch through indexed or per-step access.

### Removed
- Removed the user-facing output-generation subcommand from `beachx config`; high-level BEACH authoring keys are now accepted directly by the Fortran config loader.

## [1.4.0] - 2026-07-07

### Added
- `output.restart_from` for reading resume checkpoints from a separate directory while writing new outputs to `output.dir`
- `photo_raycast` example species in the default `beachx config init` template
- `photo_escape_model="boltzmann_cutoff"` for `photo_raycast` species to approximate photoelectron return/escape with a reduced Boltzmann closure

### Changed
- README now focuses on quick start and links detailed usage to `docs/`
- `periodic2` far correction now defaults to `none`; `auto` remains a compatibility alias for `none`, and `m2l_root_oracle` is an explicit opt-in diagnostic mode
- Moved `m2l_root_oracle` far-correction diagnostics out of the default heavy Fortran target list and behind `make test-fortran-far-correction`

### Fixed
- Maxwellian velocity sampling now bounds extreme tails to avoid runaway particle steps
- `periodic2` collision search now guards against runaway periodic image enumeration

## [1.3.0] - 2026-06-22

### Added
- `beachx lint <beach.toml>` for TOML parsing, packaged JSON Schema validation, high-level config normalization, and semantic BEACH config checks

### Removed
- **BREAKING**: Removed the preset/case config layer, including `beachx preset`, `case.toml`, preset package data, and case/preset schemas

### Changed
- Fortran config loading now uses `toml-f` for standard TOML syntax, including multiline arrays, dotted keys, inline tables, and literal `#` characters in strings
- Python config normalization now uses `tomli-w` instead of maintaining duplicate in-project TOML writer code
- `beachx config` now creates, validates, normalizes, and diffs direct `beach.toml` files
- BEACH context plugin references now document the direct `beach.toml` workflow

### Fixed
- Serial builds now reject multi-rank launcher configuration before entering MPI initialization

## [1.2.0] - 2026-06-14

### Added
- Codex context plugin for BEACH repository workflows
- Tiered developer test targets (`test-l0`, `test-l1`, `test-l2`, `test-l3`, `test-heavy`, `test-full`)
- Object surface material models, including limited free-space floating conductor support
- Stable development version metadata mode for faster incremental fpm builds

### Changed
- **BREAKING**: `sim.batch_count` is now the cumulative target batch count when `output.resume=true`; resuming from `batches=100` with `batch_count=150` runs 50 more batches
- Resume progress, history, and workload estimates now use remaining batches derived from checkpoint state
- Surface material validation and documentation are aligned with current conductor/dielectric support

### Fixed
- Resume now fails when required checkpoint files are missing instead of silently starting a new run
- Restart, config, CLI, and Python result readers now reject non-finite or invalid numeric values more consistently
- Zero-softening self-field singularities are skipped in direct and tree field paths
- Workload estimator accepts the current normalization and sheath configuration keys

## [1.1.0] - 2026-05-27

### Added
- Internal field normalization controls (`field_normalization`, `field_length_scale`) for direct, treecode, and FMM Coulomb kernels
- External electric-field and velocity-grid injection support
- Field-kernel scene analysis utilities
- Batch-duration stability theory documentation
- Agent user guide for AI-assisted simulation workflows

### Changed
- Periodic mesh visualization can preserve whole mesh objects while wrapping periodic domains
- Coulomb force matrix plotting supports parallel workers for larger mesh-group analyses
- JSON Schema numeric constraints now consistently reject nonzero-only invalid values

### Fixed
- Small-scale collision tolerances for regolith-scale geometries

## [1.0.0] - 2026-03-26

### Added
- `beachx` unified CLI with subcommands: `inspect`, `animate`, `workload`, `slices`, `coulomb`, `mobility`, `profile`, `config`, `preset`, `model`
- Preset-based config workflow with validation and preset management commands
- JSON Schema for `beach.toml` config validation
- OBJ mesh transform support (scale, rotation, offset) with CRLF line ending handling
- `periodic2` support for Coulomb force calculation and 3D electric field line plotting
- Periodic mesh replication for visualization (`periodic2_repeat`)
- Fortran potential history output (`write_potential_history`, `potential_history.csv`)
- OpenMP parallelization for mesh potential computation and FMM evaluation
- `beachx model close-pack` for generating close-packed sphere models
- Sheath injection model (`zhao_auto`) with configurable parameters
- 3D electric field line tracing and plotting (`trace_field_lines`, `plot_field_lines_3d`)
- Coulomb mobility analysis (`analyze_coulomb_mobility`)
- `Beach` facade class for high-level Python API
- Performance profiling output and `beachx profile` visualization
- Comprehensive Fortran test suite (20 test targets)
- Python post-processing API documentation (`docs/python_postprocess_api.md`)
- SPEC.md behavioral specification

### Changed
- **BREAKING**: Old CLI entry points (`beach-inspect`, `beach-animate-history`, etc.) are deprecated in favor of `beachx` subcommands
- `sim_stats` particle counters upgraded from 32-bit to 64-bit integers to prevent overflow in large simulations
- `total_particles_from_config` now validates against integer overflow
- `plot_potential()` default colormap unified to `"viridis"` (was `"jet"` in facade)
- Internal `_select_frame_columns` and `_surface_charge_density` removed from public `__all__`

### Fixed
- Empty CSV handling in `load_fortran_result` for edge cases with zero elements
- Empty `mesh_triangles.csv` no longer causes `IndexError`
- Invalid `m2l_root_trunc` reference removed from example config comments

## [0.8.0] - 2026-03-08

### Added
- Preset-based `beachx config` workflow with high-level spatial config normalization
- JSON Schema for case and preset TOML files
- `beachx preset save` and `beachx preset edit` commands
- High-level config validation

### Changed
- Removed unused reserved sim keys

## [0.7.0] - 2026-02-15

### Added
- `periodic2` m2l_root_oracle far correction (Ewald residual build-time fit)
- Exact periodic2 Ewald far correction
- Mesh view angle controls for Python plotting

### Changed
- Reorganized FMM source layout into subdirectories
- Switched Fortran formatter to 2-space indent

## [0.6.0] - 2026-01-25

### Added
- Plate-hole template and per-cap cylinder controls
- Simulation-box potential slice plotting with configurable vmin/vmax
- Expanded Fortran test cases

### Changed
- Split large Fortran modules into focused submodules (simulator, field_solver, config_parser)
- Organized split modules into subdirectories

## [0.5.0] - 2026-01-12

### Added
- Treecode field solver with auto-tuning heuristics
- Treecode sim config controls (`tree_theta`, `tree_leaf_max`, `tree_min_nelem`)
- Treecode parameter support in workload estimate CLI

### Fixed
- Photo raycast now ignores out-of-box mesh hits

## [0.4.0] - 2025-12-28

### Added
- Mesh-group charge APIs and Coulomb interaction calculation
- `rel_change` displayed in batch progress output

### Changed
- **BREAKING**: Simplified Fortran history API
- Renamed `calc_coulomb` args to `target`/`source`
- Default mesh step selection to latest history

## [0.3.0] - 2025-12-15

### Changed
- **BREAKING**: Removed deprecated config keys

## [0.2.0] - 2025-12-01

### Added
- pip-installable packaging with bundled Fortran binary
- MPI-rank-aware workload estimation
- Generic HPC job script template
