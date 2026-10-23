!> シミュレーション設定・統計・メッシュ・粒子・衝突情報の主要データ型を定義する。
module bem_types
  use bem_kinds, only: dp, i32, i64
  use bem_outflow_spectrum, only: outflow_spectrum_bins
  implicit none

  integer(i32), parameter :: bc_open = 0_i32
  integer(i32), parameter :: bc_reflect = 1_i32
  integer(i32), parameter :: bc_periodic = 2_i32
  integer(i32), parameter :: bc_redistributed_reflect = 3_i32

  integer(i32), parameter :: surface_model_insulator = 1_i32
  integer(i32), parameter :: surface_model_conductor = 2_i32
  integer(i32), parameter :: surface_model_dielectric = 3_i32

  !> 粒子の追跡・衝突・境界 event が使う解決済みの領域と乱数 seed。
  !! 入力の has_*、外部場の角度指定、solver の選択や batch 制御は含めない。
  type :: particle_domain_config
    integer(i32) :: rng_seed = 12345_i32
    logical :: use_box = .false.
    real(dp) :: box_min(3) = [-1.0d0, -1.0d0, -1.0d0]
    real(dp) :: box_max(3) = [1.0d0, 1.0d0, 1.0d0]
    integer(i32) :: bc_low(3) = [bc_open, bc_open, bc_open]
    integer(i32) :: bc_high(3) = [bc_open, bc_open, bc_open]
    character(len=16) :: field_bc_mode = 'free'
  end type particle_domain_config

  !> photo ray に必要な領域、解決済み batch 幅、反射回数上限。
  type, extends(particle_domain_config) :: photo_raycast_config
    real(dp) :: batch_duration = 0.0d0
    integer(i32) :: raycast_max_bounce = 16_i32
  end type photo_raycast_config

  !> 入力と run 制御の設定。physics へは必要な親 component だけを渡す。
  type, extends(photo_raycast_config) :: sim_config
    real(dp) :: dt = 1.0d-12
    integer(i32) :: batch_count = 1_i32
    logical :: has_batch_duration = .false.
    real(dp) :: batch_duration_step = 0.0d0
    logical :: has_batch_duration_step = .false.
    integer(i32) :: max_step = 100
    !> 粒子ごとの適応サブステップ。step 候補の中点電場が前回受理した中点電場からこの割合より大きく変わり、
    !! その差が step の仕事に換算して 1e-4 eV を超えるときは step を半分にして取り直す。0 で無効。
    real(dp) :: substep_field_change_max = 0.0d0
    !> 1 つの dt を半分にする回数の上限（最小 sub-step は dt/2^n）。
    integer(i32) :: substep_max_halvings = 6_i32
    !> 場の変化が小さい領域で step を `dt_s` の 2^n 倍まで伸ばす（opt-in）。有効時に `substep_field_change_max` が 0 なら 0.05 を使う。
    logical :: step_growth = .false.
    integer(i32) :: step_growth_max = 5_i32
    !> 伸長時の仕事の床値を運動エネルギーのこの割合にする（1e-4 eV との大きい方）。
    real(dp) :: step_growth_energy_fraction = 1.0d-2
    real(dp) :: tol_rel = 1.0d-4
    real(dp) :: q_floor = 1.0d-30
    character(len=16) :: field_solver = 'auto'
    character(len=16) :: field_normalization = 'si'
    real(dp) :: field_length_scale = 1.0d0
    integer(i32) :: field_periodic_image_layers = 1_i32
    character(len=16) :: field_periodic_far_correction = 'none'
    real(dp) :: field_periodic_ewald_alpha = 0.0d0
    integer(i32) :: field_periodic_ewald_layers = 4_i32
    character(len=256) :: field_periodic_cache_dir = ''
    real(dp) :: field_periodic_generation_tolerance = 1.0d-8
    real(dp) :: tree_theta = 0.5d0
    logical :: has_tree_theta = .false.
    integer(i32) :: tree_leaf_max = 16_i32
    !> 'double' (default, bitwise-reproducible) or 'single' for the FMM near-field panel sums.
    character(len=8) :: near_field_precision = 'double'
    logical :: has_tree_leaf_max = .false.
    integer(i32) :: tree_min_nelem = 256_i32
    real(dp) :: e0(3) = 0.0d0
    logical :: has_e0_vector = .false.
    real(dp) :: e0_abs = 0.0d0
    logical :: has_e0_abs = .false.
    real(dp) :: e0_phi_xy_deg = 0.0d0
    logical :: has_e0_phi_xy_deg = .false.
    real(dp) :: e0_phi_z_deg = 0.0d0
    logical :: has_e0_phi_z_deg = .false.
    real(dp) :: b0(3) = 0.0d0
    character(len=32) :: reservoir_potential_model = 'none'
    real(dp) :: phi_infty = 0.0d0
    !> バッチ粒子の rank 配分。`static` は rank ごとに固定分割、`shared` は全 rank 共有の slot pool を chunk 単位で取り合う。
    character(len=8) :: particle_distribution = 'static'
    !> `shared` で 1 回の claim が取る slot 数。
    integer(i32) :: particle_pool_chunk = 64_i32
    character(len=32) :: multiple_box_events_policy = 'abort'
    character(len=32) :: multiple_box_events_retry_backend = 'none'
    integer(i32) :: multiple_box_events_soft_discard_count_grace = 1000_i32
    real(dp) :: multiple_box_events_soft_discard_fraction_limit = 1.0d-6
    real(dp) :: multiple_box_events_soft_discard_abs_charge_limit = 1.0d-12
    integer(i32) :: injection_face_phi_grid_n = 3_i32
  end type sim_config

  !> adaptive batch（非零 mode の電位変化で batch 幅を縮める）の状態。summary の adaptive_nonzero_mode_* キー。
  type :: adaptive_nonzero_mode_stats
    integer(i64) :: rejected_trials = 0_i64
    real(dp) :: last_batch_duration = 0.0_dp
    real(dp) :: last_potential_step = 0.0_dp
    integer(i32) :: omp_threads = 0_i32
  end type adaptive_nonzero_mode_stats

  !> 外部シースの receipt。summary の matching_plane_* キー（top_field_v_m は surface_current_top_field_v_m）。
  type :: matching_plane_stats
    logical :: state_valid = .false.
    real(dp) :: displacement_c_m2 = 0.0_dp
    real(dp) :: phi_v = 0.0_dp
    ! Recomputed from the restored external root and snapshot; not persisted in checkpoints.
    real(dp) :: top_field_v_m = 0.0_dp
    real(dp) :: top_cell_mean_v_m = 0.0_dp
    real(dp) :: response(6) = 0.0_dp
    real(dp) :: feedback(4) = 0.0_dp
    real(dp) :: photoelectron_return_flux_m2_s = 0.0_dp
    real(dp) :: photoelectron_escape_flux_m2_s = 0.0_dp
    !> 外部根の光電子源に使った H 外向きスペクトル [1/(m2 s)]、bin 別。全て 0 は表面放出の半 Maxwell 源。
    real(dp) :: photoelectron_spectrum_m2_s(outflow_spectrum_bins) = 0.0_dp
    integer(i32) :: iterations = 0_i32
    real(dp) :: residual = 0.0_dp
    real(dp) :: secondary_share = 0.0_dp
    real(dp) :: secondary_budget_residual_a_m2 = 0.0_dp
  end type matching_plane_stats

  !> 処理済み粒子数、吸着/脱出数、バッチ数、最終相対変化量を集計する統計型。
  type :: sim_stats
    integer(i64) :: processed_particles = 0_i64
    integer(i64) :: absorbed = 0_i64
    integer(i64) :: escaped = 0_i64
    integer(i64) :: escaped_boundary = 0_i64
    integer(i64) :: survived_max_step = 0_i64
    integer(i64) :: multiple_box_events_retry_attempted = 0_i64
    integer(i64) :: multiple_box_events_retry_resolved = 0_i64
    integer(i64) :: multiple_box_events_soft_discarded = 0_i64
    real(dp) :: multiple_box_events_soft_discarded_abs_charge = 0.0_dp
    !> 適応サブステップ: 半分にした dt step の数、その中で受理した sub-step の総数、半分化の上限に達した step の数。
    integer(i64) :: adaptive_substep_halved_steps = 0_i64
    integer(i64) :: adaptive_substeps = 0_i64
    integer(i64) :: adaptive_substep_depth_limited_steps = 0_i64
    !> step 伸長で `dt_s` より長く進めた push の数。
    integer(i64) :: adaptive_grown_steps = 0_i64
    integer(i32) :: batches = 0
    real(dp) :: last_rel_change = -1.0d0
    real(dp) :: simulated_time = 0.0_dp
    !> 非零 mode の電位変化で batch 幅を縮める adaptive batch の状態。
    type(adaptive_nonzero_mode_stats) :: adaptive_nonzero_mode
    !> shared particle pool の claim 回数と、claim の critical 区間で費やした秒数（全 rank 合計）。
    integer(i64) :: particle_pool_claims = 0_i64
    real(dp) :: particle_pool_claim_seconds = 0.0_dp
    !> 外部シース（Zhao 根と z-high の matching-plane）の receipt。summary と checkpoint に出る。
    type(matching_plane_stats) :: matching_plane
    integer(i64) :: secondary_emission_hits = 0_i64
    integer(i64) :: secondary_emission_count = 0_i64
    integer(i64) :: secondary_generation_capped_count = 0_i64
    real(dp) :: secondary_emitted_charge_c = 0.0_dp
    real(dp) :: secondary_yield_sum = 0.0_dp
  end type sim_stats

  !> 種ごとのマクロ粒子端数を保持し、再開時にも注入期待値を保つ。
  type :: injection_state
    real(dp), allocatable :: macro_residual(:)
    real(dp), allocatable :: boundary_macro_residual(:, :)
  end type injection_state

  !> 三角形頂点と前計算幾何量、要素電荷を保持する境界要素メッシュ型。
  type :: mesh_type
    integer(i32) :: nelem = 0
    real(dp), allocatable :: v0(:, :)
    real(dp), allocatable :: v1(:, :)
    real(dp), allocatable :: v2(:, :)
    real(dp), allocatable :: centers(:, :)
    real(dp), allocatable :: center_x(:)
    real(dp), allocatable :: center_y(:)
    real(dp), allocatable :: center_z(:)
    real(dp), allocatable :: normals(:, :)
    real(dp), allocatable :: bb_min(:, :)
    real(dp), allocatable :: bb_max(:, :)
    real(dp), allocatable :: h_elem(:)
    real(dp), allocatable :: panel_area(:)
    real(dp), allocatable :: panel_moment1(:, :)
    real(dp), allocatable :: panel_moment2(:, :, :)
    real(dp), allocatable :: panel_edge_length(:, :)
    real(dp), allocatable :: panel_edge_outward(:, :, :)
    real(dp), allocatable :: panel_quad_position(:, :, :)
    real(dp), allocatable :: panel_quad_weight(:, :)
    real(dp), allocatable :: q_elem(:)
    integer(i32), allocatable :: elem_mesh_id(:)
    integer(i32), allocatable :: elem_surface_model(:)
    real(dp), allocatable :: elem_epsilon_r(:)
    integer(i32), allocatable :: elem_vacuum_sign(:)
    real(dp), allocatable :: vacuum_normals(:, :)
    real(dp) :: grid_bb_min(3) = 0.0d0
    real(dp) :: grid_bb_max(3) = 0.0d0
    integer(i32) :: grid_ncell(3) = 1_i32
    real(dp) :: grid_inv_cell(3) = 1.0d0
    integer(i32), allocatable :: grid_cell_start(:)
    integer(i32), allocatable :: grid_cell_elem(:)
    logical :: use_collision_grid = .false.
    logical :: periodic2_collision_ready = .false.
    !> prepare_periodic2_collision_mesh が解決した周期軸と周期長。ray query のたびに設定を解き直さない。
    integer(i32) :: periodic2_axes(2) = 0_i32
    real(dp) :: periodic2_len(2) = 0.0d0
  end type mesh_type

  !> 粒子の位置・速度・物性値・生存フラグをSoA形式で保持する型。
  type :: particles_soa
    integer(i32) :: n = 0
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: v(:, :)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: m(:)
    real(dp), allocatable :: w(:)
    integer(i32), allocatable :: species_id(:)
    !> `photo_raycast` の放出元要素。表面起源でない粒子は -1。
    integer(i32), allocatable :: source_element(:)
    logical, allocatable :: alive(:)
  end type particles_soa

  !> 衝突有無、命中要素番号、線分パラメータ、交点座標を返す衝突結果型。
  type :: hit_info
    logical :: has_hit = .false.
    integer(i32) :: elem_idx = -1
    real(dp) :: t = 0.0d0
    real(dp) :: pos(3) = 0.0d0
    integer(i32) :: image_shift(2) = 0_i32
    real(dp) :: pos_wrapped(3) = 0.0d0
  end type hit_info

end module bem_types
