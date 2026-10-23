# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Changed

- **BREAKING**: `reference_point` の既定を、電位の計算と描画で共通の `None`（絶対電位）にした。
  `plot_potential_mesh`、`plot_potential_slices`、`animate_history_mesh(quantity="potential")` と `Beach.plot_potential`、
  `Beach.plot_potential_slices`、`Beach.animate_mesh` は、これまで species 1 の注入面中心を基準にした電位差を既定で描いていた。
  `compute_potential_*` は以前から `None` が既定なので、同じ引数で計算と描画の値が一致する。`mesh_potential.csv` があれば
  potential mesh の描画は kernel を使わずに保存値を使う。従来の図は `reference_point="species1_injection_center"` で得られる。
  `beachx inspect` の potential plot と `--recompute-potential` の要約、`beachx slices`、`beachx animate --quantity potential`
  も絶対電位になる。

### Fixed

- `animate_history_mesh` と `Beach.animate_mesh`（`beachx animate`）は `view_elev` / `view_azim` を省略すると matplotlib の
  既定（仰角 30°、方位角 -60°）で描いていた。docstring と CLI help のとおり、他の mesh 描画と同じ 24° / -58° を使う。
  描画の既定の view 角度と colormap は `beach/fortran_results/_defaults.py` の定数だけが持ち、関数・`Beach`・CLI の
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
