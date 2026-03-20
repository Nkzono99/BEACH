!> Coulomb FMM コアで共有する型定義。
module bem_coulomb_fmm_types
  use bem_kinds, only: dp, i32
  implicit none
  private

  real(dp), parameter, public :: inv_sqrt_pi = 0.56418958354775628695d0

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
    real(dp) :: softening = 0.0d0
    logical :: use_periodic2 = .false.
    character(len=16) :: periodic_far_correction = 'auto'
    integer(i32) :: periodic_axes(2) = 0_i32
    real(dp) :: periodic_len(2) = 0.0d0
    integer(i32) :: periodic_image_layers = 1_i32
    real(dp) :: periodic_ewald_alpha = 0.0d0
    integer(i32) :: periodic_ewald_layers = 4_i32
    real(dp) :: target_box_min(3) = 0.0d0
    real(dp) :: target_box_max(3) = 0.0d0
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
    type(periodic2_ewald_data_type) :: periodic_ewald = periodic2_ewald_data_type()
    logical :: periodic_root_operator_ready = .false.
    real(dp), allocatable :: periodic_root_operator(:, :)
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
  end type fmm_state_type

contains

  !> periodic2 Ewald 用の補助データを初期状態へ戻す。
  !! @param[inout] data リセット対象の補助データ。
  subroutine reset_periodic2_ewald_data(data)
    type(periodic2_ewald_data_type), intent(inout) :: data

    if (allocated(data%screen_shift1)) deallocate (data%screen_shift1)
    if (allocated(data%screen_shift2)) deallocate (data%screen_shift2)
    if (allocated(data%inner_shift1)) deallocate (data%inner_shift1)
    if (allocated(data%inner_shift2)) deallocate (data%inner_shift2)
    if (allocated(data%k1)) deallocate (data%k1)
    if (allocated(data%k2)) deallocate (data%k2)
    if (allocated(data%kmag)) deallocate (data%kmag)
    if (allocated(data%karg0)) deallocate (data%karg0)
    if (allocated(data%kpref1)) deallocate (data%kpref1)
    if (allocated(data%kpref2)) deallocate (data%kpref2)
    if (allocated(data%kprefz)) deallocate (data%kprefz)
    data%ready = .false.
    data%axis1 = 0_i32
    data%axis2 = 0_i32
    data%axis_free = 0_i32
    data%nimg = 0_i32
    data%img_outer = 0_i32
    data%kmax = 0_i32
    data%alpha = 0.0d0
    data%soft2 = 0.0d0
    data%cell_area = 0.0d0
    data%k0_pref = 0.0d0
    data%screen_count = 0_i32
    data%inner_count = 0_i32
    data%k_count = 0_i32
  end subroutine reset_periodic2_ewald_data

  !> FMM plan を初期状態へ戻す。
  !! @param[inout] plan リセット対象の FMM plan。
  subroutine reset_fmm_plan(plan)
    type(fmm_plan_type), intent(inout) :: plan

    if (allocated(plan%alpha)) deallocate (plan%alpha)
    if (allocated(plan%alpha_degree)) deallocate (plan%alpha_degree)
    if (allocated(plan%alpha_factorial)) deallocate (plan%alpha_factorial)
    if (allocated(plan%alpha_sign)) deallocate (plan%alpha_sign)
    if (allocated(plan%alpha_map)) deallocate (plan%alpha_map)
    if (allocated(plan%alpha_plus_axis)) deallocate (plan%alpha_plus_axis)
    if (allocated(plan%deriv_alpha)) deallocate (plan%deriv_alpha)
    if (allocated(plan%deriv_degree)) deallocate (plan%deriv_degree)
    if (allocated(plan%deriv_factorial)) deallocate (plan%deriv_factorial)
    if (allocated(plan%deriv_map)) deallocate (plan%deriv_map)
    if (allocated(plan%alpha_beta_deriv_idx)) deallocate (plan%alpha_beta_deriv_idx)
    if (allocated(plan%eval_exp)) deallocate (plan%eval_exp)
    if (allocated(plan%eval_deriv_idx)) deallocate (plan%eval_deriv_idx)
    if (allocated(plan%eval_inv_factorial)) deallocate (plan%eval_inv_factorial)
    if (allocated(plan%src_pos)) deallocate (plan%src_pos)
    if (allocated(plan%elem_order)) deallocate (plan%elem_order)
    if (allocated(plan%node_start)) deallocate (plan%node_start)
    if (allocated(plan%node_count)) deallocate (plan%node_count)
    if (allocated(plan%child_count)) deallocate (plan%child_count)
    if (allocated(plan%child_idx)) deallocate (plan%child_idx)
    if (allocated(plan%child_octant)) deallocate (plan%child_octant)
    if (allocated(plan%node_depth)) deallocate (plan%node_depth)
    if (allocated(plan%node_level_start)) deallocate (plan%node_level_start)
    if (allocated(plan%node_level_nodes)) deallocate (plan%node_level_nodes)
    if (allocated(plan%node_center)) deallocate (plan%node_center)
    if (allocated(plan%node_half_size)) deallocate (plan%node_half_size)
    if (allocated(plan%node_radius)) deallocate (plan%node_radius)
    if (allocated(plan%target_child_count)) deallocate (plan%target_child_count)
    if (allocated(plan%target_child_idx)) deallocate (plan%target_child_idx)
    if (allocated(plan%target_child_octant)) deallocate (plan%target_child_octant)
    if (allocated(plan%target_node_depth)) deallocate (plan%target_node_depth)
    if (allocated(plan%target_level_start)) deallocate (plan%target_level_start)
    if (allocated(plan%target_level_nodes)) deallocate (plan%target_level_nodes)
    if (allocated(plan%target_node_center)) deallocate (plan%target_node_center)
    if (allocated(plan%target_node_half_size)) deallocate (plan%target_node_half_size)
    if (allocated(plan%target_node_radius)) deallocate (plan%target_node_radius)
    if (allocated(plan%source_leaf_nodes)) deallocate (plan%source_leaf_nodes)
    if (allocated(plan%leaf_nodes)) deallocate (plan%leaf_nodes)
    if (allocated(plan%leaf_slot_of_node)) deallocate (plan%leaf_slot_of_node)
    if (allocated(plan%near_start)) deallocate (plan%near_start)
    if (allocated(plan%near_nodes)) deallocate (plan%near_nodes)
    if (allocated(plan%near_source_start)) deallocate (plan%near_source_start)
    if (allocated(plan%near_source_idx)) deallocate (plan%near_source_idx)
    if (allocated(plan%near_source_shift1)) deallocate (plan%near_source_shift1)
    if (allocated(plan%near_source_shift2)) deallocate (plan%near_source_shift2)
    if (allocated(plan%far_start)) deallocate (plan%far_start)
    if (allocated(plan%far_nodes)) deallocate (plan%far_nodes)
    if (allocated(plan%m2l_target_nodes)) deallocate (plan%m2l_target_nodes)
    if (allocated(plan%m2l_source_nodes)) deallocate (plan%m2l_source_nodes)
    if (allocated(plan%m2l_shift_idx1)) deallocate (plan%m2l_shift_idx1)
    if (allocated(plan%m2l_shift_idx2)) deallocate (plan%m2l_shift_idx2)
    if (allocated(plan%m2l_target_start)) deallocate (plan%m2l_target_start)
    if (allocated(plan%m2l_pair_order)) deallocate (plan%m2l_pair_order)
    if (allocated(plan%source_parent_of)) deallocate (plan%source_parent_of)
    if (allocated(plan%parent_of)) deallocate (plan%parent_of)
    if (allocated(plan%m2m_term_count)) deallocate (plan%m2m_term_count)
    if (allocated(plan%m2m_alpha_list)) deallocate (plan%m2m_alpha_list)
    if (allocated(plan%m2m_delta_list)) deallocate (plan%m2m_delta_list)
    if (allocated(plan%l2l_term_count)) deallocate (plan%l2l_term_count)
    if (allocated(plan%l2l_gamma_list)) deallocate (plan%l2l_gamma_list)
    if (allocated(plan%l2l_delta_list)) deallocate (plan%l2l_delta_list)
    if (allocated(plan%shift_axis1)) deallocate (plan%shift_axis1)
    if (allocated(plan%shift_axis2)) deallocate (plan%shift_axis2)
    call reset_periodic2_ewald_data(plan%periodic_ewald)
    if (allocated(plan%periodic_root_operator)) deallocate (plan%periodic_root_operator)
    if (allocated(plan%m2l_deriv)) deallocate (plan%m2l_deriv)
    if (allocated(plan%source_p2m_basis)) deallocate (plan%source_p2m_basis)
    if (allocated(plan%source_shift_monomial)) deallocate (plan%source_shift_monomial)
    if (allocated(plan%target_shift_monomial)) deallocate (plan%target_shift_monomial)
    plan%options = fmm_options_type()
    plan%built = .false.
    plan%nsrc = 0_i32
    plan%ncoef = 0_i32
    plan%nderiv = 0_i32
    plan%eval_term_count = 0_i32
    plan%max_node = 0_i32
    plan%nnode = 0_i32
    plan%node_max_depth = 0_i32
    plan%target_tree_ready = .false.
    plan%target_max_node = 0_i32
    plan%target_nnode = 0_i32
    plan%target_node_max_depth = 0_i32
    plan%nsource_leaf = 0_i32
    plan%nleaf = 0_i32
    plan%m2l_pair_count = 0_i32
    plan%m2l_build_count = 0_i32
    plan%m2l_visit_count = 0_i32
    plan%periodic_root_operator_ready = .false.
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
    state%ready = .false.
    state%update_count = 0_i32
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
    call initialize_fmm_state(state)
  end subroutine reset_fmm_state

end module bem_coulomb_fmm_types
