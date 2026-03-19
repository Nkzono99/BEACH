!> 粒子位置での電場評価を direct / treecode / fmm で切り替える場ソルバ。
module bem_field_solver
!$ use omp_lib
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config, bc_periodic
  use bem_field, only: electric_field_at
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type
  implicit none
  private

  real(dp), parameter :: charge_cancellation_tol = 1.0d-10

  type :: field_solver_type
    character(len=16) :: mode = 'direct'
    character(len=16) :: field_bc_mode = 'free'
    real(dp) :: softening = 1.0d-6
    real(dp) :: theta = 0.5d0
    integer(i32) :: leaf_max = 16_i32
    integer(i32) :: min_nelem = 256_i32
    logical :: use_periodic2 = .false.
    integer(i32) :: periodic_axes(2) = 0_i32
    real(dp) :: periodic_len(2) = 0.0d0
    integer(i32) :: periodic_image_layers = 1_i32
    character(len=16) :: periodic_far_correction = 'none'
    real(dp) :: periodic_ewald_alpha = 0.0d0
    integer(i32) :: periodic_ewald_layers = 4_i32
    real(dp) :: target_box_min(3) = 0.0d0
    real(dp) :: target_box_max(3) = 0.0d0
    logical :: tree_ready = .false.
    integer(i32) :: nelem = 0_i32
    integer(i32) :: max_node = 0_i32
    integer(i32) :: nnode = 0_i32
    integer(i32), allocatable :: elem_order(:)
    integer(i32), allocatable :: node_start(:), node_count(:)
    integer(i32), allocatable :: child_count(:), child_idx(:, :), child_octant(:, :)
    integer(i32), allocatable :: node_depth(:)
    integer(i32) :: node_max_depth = 0_i32
    integer(i32), allocatable :: node_level_start(:), node_level_nodes(:)
    real(dp), allocatable :: node_center(:, :)
    real(dp), allocatable :: node_half_size(:, :)
    real(dp), allocatable :: node_radius(:)
    real(dp), allocatable :: node_q(:), node_abs_q(:)
    real(dp), allocatable :: node_qx(:), node_qy(:), node_qz(:)
    real(dp), allocatable :: node_charge_center(:, :)
    logical :: fmm_ready = .false.
    integer(i32) :: nleaf = 0_i32
    integer(i32), allocatable :: leaf_nodes(:)
    integer(i32), allocatable :: leaf_slot_of_node(:)
    logical :: target_tree_ready = .false.
    integer(i32) :: target_max_node = 0_i32
    integer(i32) :: target_nnode = 0_i32
    integer(i32), allocatable :: target_child_count(:), target_child_idx(:, :), target_child_octant(:, :)
    integer(i32), allocatable :: target_node_depth(:)
    integer(i32) :: target_node_max_depth = 0_i32
    integer(i32), allocatable :: target_level_start(:), target_level_nodes(:)
    real(dp), allocatable :: target_node_center(:, :)
    real(dp), allocatable :: target_node_half_size(:, :)
    real(dp), allocatable :: target_node_radius(:)
    integer(i32), allocatable :: near_start(:), near_nodes(:)
    integer(i32), allocatable :: far_start(:), far_nodes(:)
    integer(i32) :: fmm_m2l_pair_count = 0_i32
    integer(i32) :: fmm_m2l_build_count = 0_i32
    integer(i32) :: fmm_m2l_visit_count = 0_i32
    integer(i32) :: fmm_near_interaction_count = 0_i32
    integer(i32) :: fmm_far_interaction_count = 0_i32
    integer(i32), allocatable :: fmm_m2l_target_nodes(:), fmm_m2l_source_nodes(:)
    integer(i32), allocatable :: fmm_m2l_target_start(:), fmm_m2l_pair_order(:)
    integer(i32), allocatable :: fmm_parent_of(:)
    real(dp), allocatable :: fmm_node_local_e0(:, :)
    real(dp), allocatable :: fmm_node_local_jac(:, :, :)
    real(dp), allocatable :: fmm_node_local_hess(:, :, :, :)
    real(dp), allocatable :: fmm_shift_axis1(:), fmm_shift_axis2(:)
    real(dp), allocatable :: leaf_far_e0(:, :)
    real(dp), allocatable :: leaf_far_jac(:, :, :)
    real(dp), allocatable :: leaf_far_hess(:, :, :, :)
    logical :: fmm_profile_enabled = .false.
    integer(i32) :: fmm_refresh_count = 0_i32
    real(dp) :: fmm_last_moment_time_s = 0.0d0
    real(dp) :: fmm_last_clear_time_s = 0.0d0
    real(dp) :: fmm_last_m2l_time_s = 0.0d0
    real(dp) :: fmm_last_l2l_time_s = 0.0d0
    real(dp) :: fmm_last_copy_time_s = 0.0d0
    real(dp) :: fmm_last_refresh_time_s = 0.0d0
    real(dp) :: fmm_total_moment_time_s = 0.0d0
    real(dp) :: fmm_total_clear_time_s = 0.0d0
    real(dp) :: fmm_total_m2l_time_s = 0.0d0
    real(dp) :: fmm_total_l2l_time_s = 0.0d0
    real(dp) :: fmm_total_copy_time_s = 0.0d0
    real(dp) :: fmm_total_refresh_time_s = 0.0d0
    logical :: fmm_use_core = .false.
    logical :: fmm_core_ready = .false.
    type(fmm_options_type) :: fmm_core_options = fmm_options_type()
    type(fmm_plan_type) :: fmm_core_plan
    type(fmm_state_type) :: fmm_core_state
  contains
    procedure :: init => init_field_solver
    procedure :: refresh => refresh_field_solver
    procedure :: eval_e => eval_e_field_solver
  end type field_solver_type

  public :: field_solver_type

  interface
    !> 設定とメッシュから電場ソルバを初期化する。
    module subroutine init_field_solver(self, mesh, sim)
      class(field_solver_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      type(sim_config), intent(in) :: sim
    end subroutine init_field_solver

    !> 要素数に応じて treecode の代表パラメータを推定する。
    pure module subroutine estimate_auto_tree_params(nelem, theta, leaf_max)
      integer(i32), intent(in) :: nelem
      real(dp), intent(out) :: theta
      integer(i32), intent(out) :: leaf_max
    end subroutine estimate_auto_tree_params

    !> 現在の要素電荷から treecode/FMM モーメントを再計算する。
    module subroutine refresh_field_solver(self, mesh)
      class(field_solver_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
    end subroutine refresh_field_solver

    !> 観測点 `r` の電場を設定されたソルバで評価する。
    module subroutine eval_e_field_solver(self, mesh, r, e)
      class(field_solver_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: e(3)
    end subroutine eval_e_field_solver

    !> メッシュ重心を octree 分割して treecode トポロジを構築する。
    module subroutine build_tree_topology(self, mesh)
      class(field_solver_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
    end subroutine build_tree_topology

    !> 要素添字区間を1ノードとして登録し、必要なら8分割で子ノードを作る。
    recursive module subroutine build_node(self, mesh, node_idx, start_idx, end_idx, depth)
      class(field_solver_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      integer(i32), intent(in) :: node_idx, start_idx, end_idx, depth
    end subroutine build_node

    !> treecode 判定に使う8分木の octant 添字を返す。
    pure module function octant_index(x, y, z, center) result(oct)
      real(dp), intent(in) :: x, y, z
      real(dp), intent(in) :: center(3)
      integer(i32) :: oct
    end function octant_index

    !> ノードを再帰走査し、葉では direct 総和、遠方は monopole 近似を適用する。
    recursive module subroutine traverse_node(self, mesh, node_idx, rx, ry, rz, soft2, ex, ey, ez)
      class(field_solver_type), intent(in) :: self
      type(mesh_type), intent(in) :: mesh
      integer(i32), intent(in) :: node_idx
      real(dp), intent(in) :: rx, ry, rz, soft2
      real(dp), intent(inout) :: ex, ey, ez
    end subroutine traverse_node

    !> ノード半径・距離と電荷打ち消し度合いで遠方近似を判定する。
    module function accept_node(self, node_idx, rx, ry, rz) result(accept_it)
      class(field_solver_type), intent(in) :: self
      integer(i32), intent(in) :: node_idx
      real(dp), intent(in) :: rx, ry, rz
      logical :: accept_it
    end function accept_node

    !> ノード配列を要求サイズで確保し、未使用要素をゼロ初期化する。
    module subroutine ensure_tree_capacity(self, max_node_needed)
      class(field_solver_type), intent(inout) :: self
      integer(i32), intent(in) :: max_node_needed
    end subroutine ensure_tree_capacity

    !> source tree ノードを深さごとの連続バケットへ並べ替える。
    module subroutine rebuild_source_level_cache(self)
      class(field_solver_type), intent(inout) :: self
    end subroutine rebuild_source_level_cache

    !> treecode/FMM 作業配列を解放する。
    module subroutine reset_tree_storage(self)
      class(field_solver_type), intent(inout) :: self
    end subroutine reset_tree_storage

    !> core FMM plan の可観測メタデータを solver view へ同期する。
    module subroutine sync_core_plan_view(self)
      class(field_solver_type), intent(inout) :: self
    end subroutine sync_core_plan_view

    !> mesh 重心から core FMM 用の source 座標配列 `src_pos(3,n)` を作る。
    module subroutine build_core_source_positions(mesh, src_pos)
      type(mesh_type), intent(in) :: mesh
      real(dp), allocatable, intent(out) :: src_pos(:, :)
    end subroutine build_core_source_positions

    !> core FMM plan 由来の pair 数・near/far 統計値を solver へ同期する。
    module subroutine sync_core_plan_stats(self)
      class(field_solver_type), intent(inout) :: self
    end subroutine sync_core_plan_stats

    !> OpenMP 有効時は `omp_get_wtime`、それ以外は `cpu_time` を返す簡易タイマ。
    module function field_solver_time_seconds() result(time_s)
      real(dp) :: time_s
    end function field_solver_time_seconds

    !> ASCII英字を小文字化する。
    pure module function lower_ascii(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
    end function lower_ascii
  end interface

end module bem_field_solver
