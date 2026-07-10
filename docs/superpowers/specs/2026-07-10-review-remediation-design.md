# Review Remediation Design

## Goal

先行レビューと物理レビューで確認された問題を、途中段階では fail-fast で誤計算を防ぎつつ、最終的には粒子積分、周期場、電位境界、再開、出力が一つの物理・数値契約を共有する設計へ統合する。

段階導入は実装と検証を小さく保つための手段であり、暫定 guard を最終状態として残さない。最終状態は、レビューで指摘された組み合わせを可能な限り再び利用可能にした上で、同じ場、同じゲージ、同じ時間積分規約に従わせる。

## Scope

本設計は次を対象とする。

- `output.dir` の shell injection
- treecode の混合符号電荷に対する不正な monopole 近似
- Boris pusher の位置更新と時刻規約
- box crossing、反射、periodic wrap、mesh collision のイベント順序
- periodic collision の 4096 画像 fail-open と整数変換範囲
- MPI reservoir の rank 同期 burst と global residual
- restart の mesh identity、世代整合、atomic commit、history 整合
- Python の charge/potential history 読み込みと破損検出
- inspect CLI の暗黙 `O(N^2)` 電位再計算
- `periodic2` の有限画像集合、canonical seam、softening、Ewald 安定性
- periodic root operator の電位ゲージと精度契約
- root oracle の QR、MPI 重複、cache、最終的な解析 M2L
- open、reservoir、photo が使う電位と基準面の統一
- photo self term と条件付き escape VDF
- Zhao sheath の境界 closure と prescribed field の整合
- SPEC、schema、docs、examples と実装の同期

現在の v1 系が表面三角形を重心 softened point charge として扱うことは、この改修の基準契約とする。三角形上の連続面電荷積分へ全面移行すること、および体積電荷を自己無撞着に解く PIC/Poisson solver の追加は別プロジェクトとする。ただし、将来その kernel と volume-field provider を差し替えられる境界を作る。

## Design Principles

1. **Correctness before compatibility**: 物理的に未定義、または silent data loss を起こす経路は一時的に拒否する。
2. **Temporary guards have removal criteria**: guard は対応する統合実装と acceptance test が揃ったコミットで解除する。
3. **One field contract**: 電場と電位は同じ source kernel、periodicity、softening、non-neutral closure、gauge を使う。
4. **One time convention**: 粒子状態は同時刻 `(x^n, v^n)` と定義し、境界・衝突判定も同じ step trajectory を使う。
5. **Committed state only**: resume は世代 manifest で commit 済みと確認できる checkpoint だけを読む。
6. **Corruption is an error**: 欠損履歴、重複要素、incomplete collision query を物理値ゼロや no-hit と解釈しない。
7. **Geometry work is reusable**: FMM operator、target topology、potential sampling grid は、電荷状態と分離して再利用する。
8. **Tests define scientific contracts**: solver 同士の一致だけでなく、解析解、収束、保存則、周期不変性を検証する。

## Final Architecture

### 1. Filesystem And Checkpoint Transactions

`bem_filesystem` を追加し、recursive directory creation、atomic rename、必要な同期を POSIX C binding 経由で提供する。ユーザー入力を shell command に連結しない。空文字、`.`、`..`、absolute path、空白、引用符、shell metacharacter を通常の path component として扱う。

Checkpoint は `.checkpoints/generation-<batch>/` に次を一旦 temp 名で書く。

- summary/stats
- charges
- rank 別 RNG state
- global injection residual
- mesh/config fingerprint
- checkpoint schema version

全 rank の書き込み完了後に barrier を行い、root rank が generation manifest を書き、directory と `current` manifest を atomic rename する。loader は `current` が指す generation の全ファイル、schema、batch、MPI metadata、fingerprint を検証してから state を変更する。旧形式 checkpoint は read-only migration path として残し、新形式で次回保存する。

Mesh fingerprint は、primary physical mesh の element order、全頂点、mesh ID、surface model、epsilon metadata を安定した byte serialization にし、SHA-256 で計算する。periodic collision 用の派生 canonical representation は fingerprint に含めない。element order は charge index の意味を決めるため、順序変更も不一致とする。

同一出力 directory へ resume する場合、history の最終 batch が committed checkpoint を越えていれば、commit batch より後の行を temp file へ除外コピーして atomic replace する。header、重複 batch、要素完全性も append 前に検証する。

### 2. Particle Step And Boundary Event Engine

粒子状態を同時刻 `(x^n, v^n)` と定義する。Boris velocity update を位置更新から分離し、step は次を行う。

1. midpoint 位置を予測し、保存的外場を含む `E(x^{n+1/2})` を評価する。
2. Boris rotation で `v^{n+1}` を求める。
3. `x^{n+1} = x^n + 0.5 * (v^n + v^{n+1}) * dt` とする。
4. step trajectory に対し、最初の mesh hit と box face crossing を比較する。
5. 最初の event 時刻まで state を進め、event 処理後に残り時間を同じ規約で再積分する。

`boundary_event_type` は event fraction、同時に交差した face mask、boundary mode を保持する。corner/edge crossing は face mask を一括処理し、軸順序に依存させない。反射は event 点で法線速度を反転し、periodic は対応する対向 face へ移し、open は energy closure を適用する。event loop 上限超過は明示的な numerical failure とする。

Collision query は `complete`, `image_limit`, `index_range`, `invalid_geometry` の status を返す。status を扱わない既存 API 呼び出しでは incomplete を `error stop` とし、simulator の OpenMP region では thread-local status を集約して region 外で停止する。画像 index は変換前に `i64` 範囲を検証する。

Trajectory はまず二次 chord dense output を使い、gyro angle、加速度変化、element representative length に基づく adaptive substep を追加する。解析的な曲線対三角形交差は要求しないが、`dt -> 0` の hit element/time 収束を acceptance condition とする。

### 3. Reservoir Counts And MPI Determinism

Reservoir の物理流束から、species ごとの global expected macro count と global residual を一度だけ更新する。global count は `mpi_split_count` で rank に分配する。これにより、期待値が1未満でも rank 数単位の同時 burst を起こさない。

Global residual は root checkpoint に一つだけ保存し、batch 開始時に broadcast する。rank 別 RNG state は exact same-MPI-size resume のため維持する。MPI size を変えた resume は、global count/state は復元できても乱数列が変わるため、別の explicit non-bitwise resume mode を追加するまで拒否する。

Python workload estimator は同じ global count/residual 規約を実装し、Fortran の batch count と一致させる。

### 4. Output And Python History Model

Fortran の charge history と potential history は、各 snapshot が全要素を一度ずつ含む dense element history と定義する。Python では共通 parser/index 基盤を使い、次を拒否する。

- batch 内の欠損 element
- 重複 element
- element index 範囲外
- batch metadata の不一致
- 非有限 charge/potential/relative change
- batch index の逆行または不正な resume 重複

`FortranRunResult` は charge history と potential history を別フィールドで保持する。potential animation と inspection は Fortran の `potential_history.csv` / `mesh_potential.csv` を優先する。再計算は明示的な `--recompute-potential` でだけ行い、単なる inspect は `O(N^2)` kernel を構築しない。

### 5. Periodic Field Kernel Contract

`periodic2` の source は primary cell 座標と image class を保持する。largest-gap unwrap は tree 構築上の派生表現に限定し、有限画像集合の物理的意味を変更しない。整数周期長だけ source または target を移動しても、field/potential が tolerance 内で不変でなければならない。

v1 の物理 kernel は全距離で次の softened point kernel とする。

```text
G_epsilon(r) = 1 / sqrt(|r|^2 + epsilon^2)
```

Near direct、finite image shell、periodic teacher、field、potential、self policy が同じ `epsilon` を使う。oracle の都合で inner softened shell を unsoftened Coulomb に置き換えない。`epsilon=0` は self interaction を除く Coulomb limit として別の明示的 branch にする。

Ewald reciprocal term は scaled-erfc/log-domain helper を field と potential で共有し、中間 `Inf*0` を生成しない。`alpha`、real cutoff、reciprocal cutoff は固定 heuristic ではなく requested tolerance から選び、推定 truncation error を metadata と performance profile に出す。

Non-neutral 2P は config で次を明示する。

- `error`: total charge が tolerance を超えれば拒否する既定値
- `charged_walls`: box の free-axis 両面に補償 sheet を置く closure

`charged_walls` は wall positions、charge split、上下 reference potential を出力 metadata に記録する。単一の暗黙 `phi_infty` で上下両側を表現しない。

### 6. Periodic M2L And Potential Gauge

移行中の root oracle は validation backend として残すが、次を修正する。

- check matrix を anchor ごとに一度だけ factor し、全 proxy RHS を block solve する
- Ewald teacher を batch 評価し、対称な reciprocal mode を再利用する
- potential と field を同時に fit する
- 全 anchor で同じ global gauge constraint を使う
- geometry/options/algorithm version fingerprint で operator を disk cache する
- MPI root だけが operator を生成し、他 rank へ broadcast する

最終 runtime backend は、2P softened Green function の multi-index derivative を次数 `2p` まで評価し、Cartesian multipole から local への translation operator を直接構築する analytic periodic M2L とする。p=4 では proxy/check sampling の代わりに、8次までの derivative table を使う。root oracle は analytic M2L との cross-check と development diagnostics に限定する。

Potential gauge は finite reference plane と関連付ける。neutral free-axis problem では指定 reference plane 上の cell-average potential をゼロとし、charged-walls では各 side の exterior plateau potential を reference とする。全 anchor の local constant coefficient はこの共通 gauge から導出する。`E = -grad(phi)`、anchor boundary continuity、source order invariance を満たすまで potential closure には使用しない。

### 7. Shared Potential Snapshot

`field_solver` に batched field/potential evaluation を追加し、batch 開始時の charge state から immutable `potential_snapshot` を作る。snapshot は次を保持する。

- open/reservoir face の tile center potential と interpolation metadata
- 全 emitting element の center/hit-point potential
- element self contribution policy
- side 別 finite reference plane potential
- periodic kernel、softening、non-neutral closure、gauge ID
- uniform external field の reference-relative potential

Uniform `e0` は無限遠電位へ変換せず、有限 reference point `r_ref` からの `-e0 dot (r-r_ref)` として potential difference に加える。open、reservoir、photo は snapshot の `delta_phi(point, reference_side)` だけを使い、個別の free-space direct sum を持たない。

Reservoir は face tile ごとに local barrier、流束、energy mapping を計算する。tile 選択確率は `Gamma_i A_i` に比例し、速度はその tile の cutoff 条件付き flux VDF から生成する。面平均 potential を非線形流束へ代入しない。

Photo escape の一般式は次とする。

```text
E_barrier = max(q * (phi_reference - phi_emit), 0)
escape_fraction = Gamma(v_cut) / Gamma(0)
```

Emitter の macro current/weight は escape fraction を反映し、生成速度も同じ cutoff の条件付き VDF から sample する。非ゼロ normal drift は tail integral を数値的に正しく使う。self potential は field kernel と同じ softened source policyを使い、自己 element 全体を除外しない。

### 8. Zhao Sheath Modes

Zhao 系は二つの明示 mode に分ける。

1. `external_closure`: 解析 sheath は計算領域外にあり、box face が finite matching plane となる。Zhao はその面の incoming VDF と reference potential を与える。
2. `prescribed_profile`: `phi_Zhao(z)` と `E_Zhao(z)=-dphi/dz` を保存的外場 provider として pusher、potential snapshot、reservoir、photo の全てに加える。

`external_closure` では reference coordinate が domain interior に入る設定を拒否する。`prescribed_profile` が完成するまで、Zhao と独立 `photo_escape_model` の二重適用を拒否する。完成後は Zhao が返す free/captured population と conditional VDF を唯一の closure とし、別の Boltzmann factor を重ねない。

Type A/B/C の各 branch は、解析密度だけでなく VDF の0次、1次、2次 moment、energy mapping、field-potential derivative を数値積分 oracle と比較する。完全な volume-charge self-consistent mode は、将来の field provider として追加可能にする。

## Temporary Guards And Removal Conditions

段階1で次を fail-fast にする。

| Temporary guard | Removal condition |
|---|---|
| `periodic2` と open/reservoir/photo potential closure | shared snapshot、periodic potential gauge、multi-anchor continuity tests が通る |
| `e0 != 0` と infinity-based closure | finite reference plane と external-field potential difference が入る |
| Zhao と独立 photo cutoff | Zhao mode と単一 conditional VDF contract が入る |
| root oracle と `softening > 0` | softened periodic teacher/analytic M2L の kernel tests が通る |
| non-neutral periodic2 の暗黙 closure | explicit `error` / `charged_walls` config と metadata が入る |

Guard は parser、schema、docs で同じ条件にする。Guard を解除する実装 commit は対応する regression/acceptance tests と docs を必ず含む。

## Delivery Stages And Commit Boundaries

### Stage 1: Immediate Safety And Data Integrity

1. Safe filesystem directory creation
2. Fail-closed periodic collision status
3. Treecode mixed-sign descent
4. Strict Python history validation
5. Inspect precomputed-only default
6. SPEC `auto -> none` synchronization and lint cleanup

### Stage 2: Particle And MPI Correctness

1. Same-time second-order Boris kinematics
2. Reflection folding and earliest boundary event engine
3. Mesh/box event integration and reflected-path collision
4. Global reservoir count/residual and MPI distribution

### Stage 3: Transactional Resume And Output

1. Mesh/config fingerprint
2. Generation checkpoint manifest and atomic filesystem operations
3. MPI commit protocol and history truncation/validation
4. Potential history Python API and animation integration

### Stage 4: Periodic Kernel Correctness

1. Stable scaled-erfc Ewald terms
2. Primary-cell image-class semantics
3. Softened 2P kernel and tolerance-controlled teacher
4. Explicit non-neutral closure
5. Common potential gauge and exact fallback removal

### Stage 5: Periodic Performance And Final M2L

1. Block QR and batch teacher
2. Redundant state update elimination
3. Operator cache and MPI root broadcast
4. Analytic periodic M2L
5. Root oracle demotion to diagnostics

### Stage 6: Unified Potential Closures

1. Batched field/potential API and snapshot
2. Local reservoir flux/VDF
3. Photo self term and conditional VDF
4. Finite reference plane with `e0`
5. Zhao external closure
6. Zhao prescribed profile
7. Removal of all temporary compatibility guards whose conditions are satisfied

Each numbered item is an independently reviewable commit or a short series of commits that keeps its focused test target green. Public API/config changes include schema, docs, examples, Python authoring support, workload estimator, and output metadata in the same stage.

## Test Strategy

### Unit And Contract Tests

- path metacharacters create literal directories without shell side effects
- mixed-sign tree nodes match direct field across cancellation ratios
- constant E displacement, pure B gyro orbit, time reversal, timestep convergence
- reflected-path and corner boundary collision ordering
- collision query never converts incomplete work into no-hit
- global reservoir count sequence is independent of MPI rank count
- checkpoint rejects same-nelem different mesh and uncommitted generations
- history rejects missing/duplicate/non-finite rows
- inspect does not invoke potential recomputation without an explicit flag
- periodic source/target integer-translation invariance
- finite shell, Ewald teacher, FMM field, and potential use one kernel
- high-aspect Ewald evaluation remains finite
- `E = -grad(phi)` and multi-anchor potential continuity
- reservoir tile flux and conditional VDF moments
- photo escape current and velocity distribution use the same cutoff
- Zhao branch density/VDF moments and prescribed field consistency

### Integration Tests

- serial and MPI runs produce the same global count/statistics within stochastic contract
- restart continuation matches uninterrupted execution for supported exact-resume mode
- crash-injected checkpoint leaves the previous generation loadable
- free and periodic2 simulations conserve the intended energy to the integrator tolerance when no absorption occurs
- potential boundary configurations match direct high-accuracy reference cases
- FMM error converges with order/tolerance against direct or spectral Ewald reference

### KUDPC Verification

`laurel*` login nodesでは編集、静的確認、`make check` 相当の軽い build までに留める。Fortran runtime tests、MPI/OpenMP tests、FMM accuracy/performance tests は一つずつ `tssrun` または `sbatch` + `srun` で計算ノードへ送る。同じ `build/` を使う `fpm test` は並列実行しない。

検証は L0 から L3 の順で行う。長時間 far-correction diagnostics は専用 opt-in target とし、通常 L1 へ戻さない。性能改善 commit は wall time、operator build count、cache hit、peak memory の基準値を記録し、精度 test を先に満たす。

## Compatibility And Migration

- 既存 API signature は optional status や追加 procedure で可能な限り維持する。
- silent fail-open と一次精度位置更新はバグとして挙動を変更し、release note に記載する。
- legacy checkpoint は読み込み可能にするが、mesh fingerprint がない場合は明示 opt-in migration と警告を要求する。
- `field_periodic_far_correction="auto"` は移行中 `none` と同義に保ち、analytic backend が既定精度・性能条件を満たした段階で、versioned migration として新しい `auto` policy を導入する。
- Python の欠損 history zero-fill は廃止し、破損ファイルを `ValueError` とする。
- Temporary guard で拒否される既存設定には、エラー文で代替設定と解除予定の実装条件を示す。

## Completion Criteria

本改修は次を全て満たしたときに完了とする。

1. Scope に列挙した各問題が regression test で再現され、修正後に通る。
2. Temporary guard のうち、最終設計で対応可能とした組み合わせが全て再有効化されている。
3. open、reservoir、photo、Zhao prescribed profile が同じ potential snapshot contract を使う。
4. periodic field と potential が同じ softened kernel、closure、gauge を使う。
5. analytic periodic M2L が既定 backend として精度・性能 acceptance を満たし、oracle は diagnostics に限定される。
6. Particle step が同時刻二次規約と earliest-event processing を使う。
7. Restart が mesh identity と committed generation を検証する。
8. Python が charge/potential history を欠損なく読み、inspect が暗黙 quadratic work を行わない。
9. Docs、SPEC、schema、examples、output metadata が実装と一致する。
10. KUDPC 計算ノード上で必要な L0-L3、MPI、far-correction tests が完了し、残る制約が明示されている。
