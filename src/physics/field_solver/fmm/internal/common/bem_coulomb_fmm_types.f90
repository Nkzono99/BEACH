!> Coulomb FMM コアで共有する型定義。
module bem_coulomb_fmm_types
  use bem_kinds, only: dp, sp, i32, i64
  use bem_mpi, only: mpi_context
  use bem_panel_geometry, only: panel_geometry_type
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type
  use bem_periodic_shell_mean, only: periodic_shell_mean_plan_type, periodic_shell_mean_state_type, &
                                     reset_periodic_shell_mean_state
  implicit none
  private

  real(dp), parameter, public :: inv_sqrt_pi = 0.56418958354775628695d0

  ! Fixed box-cell geometry shared by preflight, proxy generation and cache identity.
  real(dp), parameter, public :: periodic_operator_proxy_shell_scale = 1.5_dp
  real(dp), parameter, public :: periodic_operator_source_radius_ratio = 0.4_dp
  real(dp), parameter, public :: periodic_operator_target_radius_ratio = 0.5_dp

  public :: periodic_source_proxy_half_size
  public :: periodic2_ewald_data_type
  public :: fmm_options_type
  public :: fmm_plan_type
  public :: fmm_state_type
  public :: reset_periodic2_ewald_data
  public :: reset_fmm_plan
  public :: initialize_fmm_state
  public :: reset_fmm_state

  type :: fmm_options_type
    real(dp) :: theta = 0.5d0
    integer(i32) :: leaf_max = 16_i32
    integer(i32) :: order = 4_i32
    !> Panel sources only: keep the per-leaf near-field tables in single precision (relative to the leaf center)
    !! and sum them with the single-precision SIMD kernel. Off by default; results then differ at the 1e-6 level.
    logical :: near_field_single = .false.
    real(dp) :: softening = 0.0d0
    logical :: use_periodic2 = .false.
    character(len=16) :: periodic_far_correction = 'none'
    integer(i32) :: periodic_axes(2) = 0_i32
    real(dp) :: periodic_len(2) = 0.0d0
    integer(i32) :: periodic_image_layers = 1_i32
    real(dp) :: periodic_ewald_alpha = 0.0d0
    integer(i32) :: periodic_ewald_layers = 4_i32
    character(len=256) :: periodic_cache_dir = ''
    real(dp) :: periodic_generation_tolerance = 1.0d-8
    !> cached_kneq0 の遠方補正が使う source セルと target anchor の box 格子の各軸の分割数。負なら box と周期長から
    !! 自動で選ぶ（anchor は半対角線が短い周期長の 0.5 倍以下、source は 0.4 倍以下になるほぼ立方体のセルで、
    !! source の分割数は anchor の 2^k 倍）。粗いと作成は速いが、粒子層の高さが周期長に近いとき非零モードに
    !! 面平均が漏れる。
    integer(i32) :: periodic_operator_source_cells(3) = -1_i32
    integer(i32) :: periodic_operator_anchor_cells(3) = -1_i32
    real(dp) :: target_box_min(3) = 0.0d0
    real(dp) :: target_box_max(3) = 0.0d0
    !> cached_kneq0 の root operator を rank で分担して作るときの MPI 文脈。呼び出し側（runtime）が渡す。
    !! 既定は単一 rank で、collective は何もしない（C ABI・テスト・非 MPI build）。
    type(mpi_context) :: mpi
  end type fmm_options_type

  type :: periodic2_ewald_data_type
    logical :: ready = .false.
    integer(i32) :: axis1 = 0_i32
    integer(i32) :: axis2 = 0_i32
    integer(i32) :: axis_free = 0_i32
    integer(i32) :: nimg = 0_i32
    integer(i32) :: img_outer = 0_i32
    integer(i32) :: kmax = 0_i32
    real(dp) :: alpha = 0.0d0
    real(dp) :: soft2 = 0.0d0
    real(dp) :: cell_area = 0.0d0
    real(dp) :: k0_pref = 0.0d0
    integer(i32) :: screen_count = 0_i32
    integer(i32) :: inner_count = 0_i32
    integer(i32) :: k_count = 0_i32
    real(dp), allocatable :: screen_shift1(:), screen_shift2(:)
    real(dp), allocatable :: inner_shift1(:), inner_shift2(:)
    real(dp), allocatable :: k1(:), k2(:), kmag(:), karg0(:), kpref1(:), kpref2(:), kprefz(:)
  end type periodic2_ewald_data_type

  type :: fmm_plan_type
    type(fmm_options_type) :: options = fmm_options_type()
    !> Options, `src_pos` and `panel_geometry` are set; enough for exact-direct evaluation.
    logical :: sources_ready = .false.
    !> Tree, interaction lists and operators are built as well (implies `sources_ready`).
    logical :: built = .false.
    integer(i32) :: nsrc = 0_i32
    integer(i32) :: ncoef = 0_i32
    integer(i32) :: nderiv = 0_i32
    integer(i32), allocatable :: alpha(:, :)
    integer(i32), allocatable :: alpha_degree(:)
    real(dp), allocatable :: alpha_factorial(:)
    real(dp), allocatable :: alpha_sign(:)
    integer(i32), allocatable :: alpha_map(:, :, :)
    integer(i32), allocatable :: alpha_plus_axis(:, :)
    integer(i32), allocatable :: deriv_alpha(:, :)
    integer(i32), allocatable :: deriv_degree(:)
    real(dp), allocatable :: deriv_factorial(:)
    integer(i32), allocatable :: deriv_map(:, :, :)
    integer(i32), allocatable :: alpha_beta_deriv_idx(:, :)
    integer(i32) :: eval_term_count = 0_i32
    integer(i32), allocatable :: eval_exp(:, :)
    integer(i32), allocatable :: eval_deriv_idx(:, :)
    real(dp), allocatable :: eval_inv_factorial(:)
    real(dp), allocatable :: src_pos(:, :)
    logical :: panel_source = .false.
    type(panel_geometry_type), allocatable :: panel_geometry(:)
    integer(i32) :: max_node = 0_i32
    integer(i32) :: nnode = 0_i32
    integer(i32) :: node_max_depth = 0_i32
    integer(i32), allocatable :: elem_order(:)
    integer(i32), allocatable :: node_start(:), node_count(:)
    integer(i32), allocatable :: child_count(:), child_idx(:, :), child_octant(:, :)
    integer(i32), allocatable :: node_depth(:)
    integer(i32), allocatable :: node_level_start(:), node_level_nodes(:)
    real(dp), allocatable :: node_center(:, :)
    real(dp), allocatable :: node_half_size(:, :)
    real(dp), allocatable :: node_radius(:)
    logical :: target_tree_ready = .false.
    integer(i32) :: target_max_node = 0_i32
    integer(i32) :: target_nnode = 0_i32
    integer(i32) :: target_node_max_depth = 0_i32
    integer(i32), allocatable :: target_child_count(:), target_child_idx(:, :), target_child_octant(:, :)
    integer(i32), allocatable :: target_node_depth(:)
    integer(i32), allocatable :: target_level_start(:), target_level_nodes(:)
    real(dp), allocatable :: target_node_center(:, :)
    real(dp), allocatable :: target_node_half_size(:, :)
    real(dp), allocatable :: target_node_radius(:)
    integer(i32) :: nsource_leaf = 0_i32
    integer(i32), allocatable :: source_leaf_nodes(:)
    integer(i32) :: nleaf = 0_i32
    integer(i32), allocatable :: leaf_nodes(:)
    integer(i32), allocatable :: leaf_slot_of_node(:)
    integer(i32), allocatable :: near_start(:), near_nodes(:)
    integer(i32), allocatable :: near_source_start(:), near_source_idx(:)
    real(dp), allocatable :: near_source_shift1(:), near_source_shift2(:)
    !> Panel sources only: the near-field panels in `near_source_idx` order as translated component arrays
    !! (vertex x/y/z per vertex, normal, edge lengths, outward edge normals per component, area), so a
    !! point evaluation reads its leaf's panels without gathering them from `panel_geometry`.
    real(dp), allocatable :: near_vx(:, :), near_vy(:, :), near_vz(:, :)
    real(dp), allocatable :: near_nx(:), near_ny(:), near_nz(:)
    real(dp), allocatable :: near_edge(:, :), near_ox(:, :), near_oy(:, :), near_oz(:, :)
    real(dp), allocatable :: near_area(:)
    !> Single-precision counterparts used when options%near_field_single: vertices relative to the leaf center.
    real(sp), allocatable :: near_sp_vx(:, :), near_sp_vy(:, :), near_sp_vz(:, :)
    real(sp), allocatable :: near_sp_nx(:), near_sp_ny(:), near_sp_nz(:)
    real(sp), allocatable :: near_sp_edge(:, :), near_sp_ox(:, :), near_sp_oy(:, :), near_sp_oz(:, :)
    real(sp), allocatable :: near_sp_area(:)
    real(dp), allocatable :: near_leaf_center(:, :)
    integer(i32), allocatable :: far_start(:), far_nodes(:)
    integer(i32) :: m2l_pair_count = 0_i32
    integer(i32) :: m2l_build_count = 0_i32
    integer(i32) :: m2l_visit_count = 0_i32
    integer(i32), allocatable :: m2l_target_nodes(:), m2l_source_nodes(:)
    integer(i32), allocatable :: m2l_shift_idx1(:), m2l_shift_idx2(:)
    integer(i32), allocatable :: m2l_target_start(:), m2l_pair_order(:)
    integer(i32), allocatable :: source_parent_of(:)
    integer(i32), allocatable :: parent_of(:)
    integer(i32), allocatable :: m2m_term_count(:)
    integer(i32), allocatable :: m2m_alpha_list(:, :)
    integer(i32), allocatable :: m2m_delta_list(:, :)
    integer(i32), allocatable :: l2l_term_count(:)
    integer(i32), allocatable :: l2l_gamma_list(:, :)
    integer(i32), allocatable :: l2l_delta_list(:, :)
    real(dp), allocatable :: shift_axis1(:), shift_axis2(:)
    type(periodic2_ewald_data_type) :: periodic_ewald
    logical :: periodic_root_operator_ready = .false.
    !> 遠方補正の (target anchor, source cell) 組の数。anchor は box を各軸 periodic_operator_anchor_cells に分けた
    !! 格子で、適応木とは独立に点評価で参照する。operator は panel を含む source cell だけを持ち、
    !! 組 k = (slot-1)*anchor_count + anchor で並ぶ。slot は periodic_operator_cell_slot(cell)、
    !! 0 のセルは operator を持たない（線源を置けない）。
    integer(i32) :: periodic_root_target_count = 0_i32
    integer(i32) :: periodic_operator_anchor_count = 0_i32
    integer(i32) :: periodic_operator_anchor_cells(3) = 0_i32
    integer(i32) :: periodic_operator_source_cells(3) = 0_i32
    real(dp) :: periodic_anchor_half(3) = 0.0_dp
    real(dp) :: periodic_source_half(3) = 0.0_dp
    real(dp), allocatable :: periodic_anchor_center(:, :)
    integer(i32) :: periodic_operator_built_count = 0_i32
    integer(i32), allocatable :: periodic_operator_cell_slot(:)
    integer(i32), allocatable :: periodic_operator_cells(:)
    !> 今回の初期化で新たに作った source cell の数（cache にあったセルは数えない）。
    integer(i32) :: periodic_operator_new_cell_count = 0_i32
    integer(i32) :: periodic_root_source_count = 0_i32
    integer(i32), allocatable :: periodic_source_cell(:)
    real(dp), allocatable :: periodic_source_center(:, :), periodic_source_p2m_basis(:, :)
    real(dp), allocatable :: periodic_root_operator(:, :, :)
    logical :: periodic_cache_hit = .false.
    integer(i32) :: periodic_operator_build_count = 0_i32
    integer(i32) :: periodic_qr_preparation_count = 0_i32
    integer(i32) :: periodic_operator_local_target_count = 0_i32
    integer(i32) :: periodic_operator_thread_count = 0_i32
    integer(i32) :: periodic_operator_class_count = 0_i32
    integer(i64) :: periodic_operator_ewald_evaluation_count = 0_i64
    character(len=16) :: periodic_cache_fingerprint = ''
    character(len=512) :: periodic_cache_path = ''
    logical :: periodic_k0_ready = .false.
    integer(i32) :: periodic_k0_free_axis = 0_i32
    type(periodic_zero_mode_plan_type) :: periodic_k0_plan
    !> 有限画像殻の面平均のうち、対称な平均場を除いた滑らかな残り（cached_kneq0 で実行時に引く）。
    type(periodic_shell_mean_plan_type) :: periodic_shell_mean
    real(dp), allocatable :: m2l_deriv(:, :)
    real(dp), allocatable :: source_p2m_basis(:, :)
    real(dp), allocatable :: source_shift_monomial(:, :)
    real(dp), allocatable :: target_shift_monomial(:, :)
  end type fmm_plan_type

  type :: fmm_state_type
    logical :: ready = .false.
    integer(i32) :: update_count = 0_i32
    real(dp), pointer :: src_q(:) => null()
    real(dp), pointer :: multipole(:, :) => null()
    real(dp), pointer :: local(:, :) => null()
    integer(i32), pointer :: multipole_active(:) => null()
    integer(i32), pointer :: local_active(:) => null()
    real(dp), allocatable :: periodic_multipole(:, :)
    !> 遠方補正の anchor ごとの局所展開。点評価で、点を含む anchor の展開を葉の局所展開に加える。
    real(dp), allocatable :: anchor_local(:, :)
    integer(i32), allocatable :: anchor_local_active(:)
    type(periodic_zero_mode_state_type) :: periodic_k0_state
    type(periodic_shell_mean_state_type) :: periodic_shell_mean_state
  end type fmm_state_type

contains

  !> Expand all box faces by the same padding, with proxy corner radius = scale * cell half diagonal.
  pure function periodic_source_proxy_half_size(half) result(proxy_half)
    real(dp), intent(in) :: half(3)
    real(dp) :: proxy_half(3), radius, normalized_sum, padding

    radius = norm2(half)
    normalized_sum = sum(half)/radius
    padding = radius*(sqrt(normalized_sum**2 + 3.0_dp*(periodic_operator_proxy_shell_scale**2 - 1.0_dp)) - &
                      normalized_sum)/3.0_dp
    proxy_half = half + padding
  end function periodic_source_proxy_half_size

  !> periodic2 Ewald 用の補助データを初期状態へ戻す。
  !! intent(out) deallocates every allocatable component and applies the type's default values.
  !! @param[out] data リセット対象の補助データ。
  subroutine reset_periodic2_ewald_data(data)
    type(periodic2_ewald_data_type), intent(out) :: data
  end subroutine reset_periodic2_ewald_data

  !> FMM plan を初期状態へ戻す。
  !! The plan owns only allocatable storage (no pointers), so intent(out) deallocates every component
  !! and applies the type's default values; there is no second list of fields to keep in sync.
  !! @param[out] plan リセット対象の FMM plan。
  subroutine reset_fmm_plan(plan)
    type(fmm_plan_type), intent(out) :: plan
  end subroutine reset_fmm_plan

  !> FMM state のポインタ成分を未関連状態へ初期化する。
  !! @param[inout] state 初期化対象の FMM state。
  subroutine initialize_fmm_state(state)
    type(fmm_state_type), intent(inout) :: state

    nullify (state%src_q)
    nullify (state%multipole)
    nullify (state%local)
    nullify (state%multipole_active)
    nullify (state%local_active)
    if (allocated(state%periodic_multipole)) deallocate (state%periodic_multipole)
    if (allocated(state%anchor_local)) deallocate (state%anchor_local)
    if (allocated(state%anchor_local_active)) deallocate (state%anchor_local_active)
    state%ready = .false.
    state%update_count = 0_i32
    state%periodic_k0_state = periodic_zero_mode_state_type()
    call reset_periodic_shell_mean_state(state%periodic_shell_mean_state)
  end subroutine initialize_fmm_state

  !> FMM state を初期状態へ戻す。
  !! @param[inout] state リセット対象の FMM state。
  subroutine reset_fmm_state(state)
    type(fmm_state_type), intent(inout) :: state

    if (associated(state%src_q)) deallocate (state%src_q)
    if (associated(state%multipole)) deallocate (state%multipole)
    if (associated(state%local)) deallocate (state%local)
    if (associated(state%multipole_active)) deallocate (state%multipole_active)
    if (associated(state%local_active)) deallocate (state%local_active)
    if (allocated(state%anchor_local)) deallocate (state%anchor_local)
    if (allocated(state%anchor_local_active)) deallocate (state%anchor_local_active)
    call initialize_fmm_state(state)
  end subroutine reset_fmm_state

end module bem_coulomb_fmm_types
