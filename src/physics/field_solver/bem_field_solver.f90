!> 粒子位置での電場評価を direct / treecode / fmm で切り替える場ソルバ。
module bem_field_solver
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config
  use bem_field, only: electric_field_at
  implicit none
  private

  real(dp), parameter :: charge_cancellation_tol = 1.0d-10

  type :: field_solver_type
    character(len=16) :: mode = 'direct'
    real(dp) :: softening = 1.0d-6
    real(dp) :: theta = 0.5d0
    integer(i32) :: leaf_max = 16_i32
    integer(i32) :: min_nelem = 256_i32
    logical :: tree_ready = .false.
    integer(i32) :: nelem = 0_i32
    integer(i32) :: max_node = 0_i32
    integer(i32) :: nnode = 0_i32
    integer(i32), allocatable :: elem_order(:)
    integer(i32), allocatable :: node_start(:), node_count(:)
    integer(i32), allocatable :: child_count(:), child_idx(:, :), child_octant(:, :)
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
    integer(i32), allocatable :: near_start(:), near_nodes(:)
    integer(i32), allocatable :: far_start(:), far_nodes(:)
    real(dp), allocatable :: leaf_far_e0(:, :)
    real(dp), allocatable :: leaf_far_jac(:, :, :)
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
      class(field_solver_type), intent(in) :: self
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
    recursive module subroutine build_node(self, mesh, node_idx, start_idx, end_idx)
      class(field_solver_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      integer(i32), intent(in) :: node_idx, start_idx, end_idx
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

    !> 1つのターゲット葉ノードに対して近傍/遠方ノード候補を再帰収集する。
    recursive module subroutine gather_leaf_interactions(self, target_leaf, source_node, near_buf, near_n, far_buf, far_n)
      class(field_solver_type), intent(in) :: self
      integer(i32), intent(in) :: target_leaf
      integer(i32), intent(in) :: source_node
      integer(i32), intent(inout) :: near_buf(:)
      integer(i32), intent(inout) :: near_n
      integer(i32), intent(inout) :: far_buf(:)
      integer(i32), intent(inout) :: far_n
    end subroutine gather_leaf_interactions

    !> source/target ノード対が遠方近似可能かを判定する。
    module function nodes_well_separated(self, target_node, source_node) result(accept_it)
      class(field_solver_type), intent(in) :: self
      integer(i32), intent(in) :: target_node
      integer(i32), intent(in) :: source_node
      logical :: accept_it
    end function nodes_well_separated

    !> 各葉ノードの near/far 相互作用リストを構築する。
    module subroutine build_fmm_interactions(self)
      class(field_solver_type), intent(inout) :: self
    end subroutine build_fmm_interactions

    !> FMM 遠方相互作用から葉ノード局所展開係数を更新する。
    module subroutine refresh_fmm_locals(self, mesh)
      class(field_solver_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
    end subroutine refresh_fmm_locals

    !> 観測点の属する葉ノードを返す。
    module function locate_target_leaf(self, r) result(node_idx)
      class(field_solver_type), intent(in) :: self
      real(dp), intent(in) :: r(3)
      integer(i32) :: node_idx
    end function locate_target_leaf

    !> FMM 局所展開と近傍 direct 和を用いて観測点の電場を評価する。
    module subroutine eval_e_fmm(self, mesh, r, e)
      class(field_solver_type), intent(in) :: self
      type(mesh_type), intent(in) :: mesh
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: e(3)
    end subroutine eval_e_fmm

    !> treecode/FMM 作業配列を解放する。
    module subroutine reset_tree_storage(self)
      class(field_solver_type), intent(inout) :: self
    end subroutine reset_tree_storage

    !> ASCII英字を小文字化する。
    pure module function lower_ascii(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
    end function lower_ascii
  end interface

end module bem_field_solver
