!> Coulomb FMM コアで共有する型定義。
module bem_coulomb_fmm_types
  use bem_kinds, only: dp, i32
  implicit none
  private

  real(dp), parameter, public :: inv_sqrt_pi = 0.56418958354775628695d0

  public :: fmm_options_type
  public :: fmm_plan_type
  public :: fmm_state_type

  type :: fmm_options_type
    real(dp) :: theta = 0.5d0
    integer(i32) :: leaf_max = 16_i32
    integer(i32) :: order = 4_i32
    real(dp) :: softening = 0.0d0
    logical :: use_periodic2 = .false.
    character(len=16) :: periodic_far_correction = 'none'
    integer(i32) :: periodic_axes(2) = 0_i32
    real(dp) :: periodic_len(2) = 0.0d0
    integer(i32) :: periodic_image_layers = 1_i32
    real(dp) :: periodic_ewald_alpha = 0.0d0
    integer(i32) :: periodic_ewald_layers = 4_i32
    real(dp) :: target_box_min(3) = 0.0d0
    real(dp) :: target_box_max(3) = 0.0d0
  end type fmm_options_type

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
    integer(i32), allocatable :: far_start(:), far_nodes(:)
    integer(i32) :: m2l_pair_count = 0_i32
    integer(i32) :: m2l_build_count = 0_i32
    integer(i32) :: m2l_visit_count = 0_i32
    integer(i32), allocatable :: m2l_target_nodes(:), m2l_source_nodes(:)
    integer(i32), allocatable :: m2l_target_start(:), m2l_pair_order(:)
    integer(i32), allocatable :: source_parent_of(:)
    integer(i32), allocatable :: parent_of(:)
    integer(i32), allocatable :: m2m_delta_idx(:, :)
    integer(i32), allocatable :: l2l_delta_idx(:, :)
    real(dp), allocatable :: shift_axis1(:), shift_axis2(:)
    real(dp), allocatable :: m2l_deriv(:, :)
    real(dp), allocatable :: source_shift_monomial(:, :)
    real(dp), allocatable :: target_shift_monomial(:, :)
  end type fmm_plan_type

  type :: fmm_state_type
    logical :: ready = .false.
    integer(i32) :: update_count = 0_i32
    real(dp), allocatable :: src_q(:)
    real(dp), allocatable :: multipole(:, :)
    real(dp), allocatable :: local(:, :)
  end type fmm_state_type

end module bem_coulomb_fmm_types
