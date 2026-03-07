!> シミュレーション設定・統計・メッシュ・粒子・衝突情報の主要データ型を定義する。
module bem_types
  use bem_kinds, only: dp, i32
  implicit none

  integer(i32), parameter :: bc_open = 0_i32
  integer(i32), parameter :: bc_reflect = 1_i32
  integer(i32), parameter :: bc_periodic = 2_i32

  !> 時間刻み・収束判定・バッチ回数・外部磁場など実行制御パラメータを保持する。
  type :: sim_config
    real(dp) :: dt = 1.0d-12
    integer(i32) :: rng_seed = 12345_i32
    integer(i32) :: batch_count = 1_i32
    real(dp) :: batch_duration = 0.0d0
    logical :: has_batch_duration = .false.
    real(dp) :: batch_duration_step = 0.0d0
    logical :: has_batch_duration_step = .false.
    integer(i32) :: max_step = 100
    real(dp) :: tol_rel = 1.0d-4
    real(dp) :: q_floor = 1.0d-30
    real(dp) :: softening = 1.0d-4
    logical :: use_hybrid = .true.
    real(dp) :: r_switch_factor = 3.0d0
    integer(i32) :: n_sub = 2
    real(dp) :: softening_factor = 0.1d0
    real(dp) :: b0(3) = 0.0d0
    character(len=32) :: reservoir_potential_model = 'none'
    real(dp) :: phi_infty = 0.0d0
    integer(i32) :: injection_face_phi_grid_n = 3_i32
    integer(i32) :: raycast_max_bounce = 16_i32
    logical :: use_box = .false.
    real(dp) :: box_min(3) = [-1.0d0, -1.0d0, -1.0d0]
    real(dp) :: box_max(3) = [1.0d0, 1.0d0, 1.0d0]
    integer(i32) :: bc_low(3) = [bc_open, bc_open, bc_open]
    integer(i32) :: bc_high(3) = [bc_open, bc_open, bc_open]
  end type sim_config

  !> 処理済み粒子数、吸着/脱出数、バッチ数、最終相対変化量を集計する統計型。
  type :: sim_stats
    integer(i32) :: processed_particles = 0
    integer(i32) :: absorbed = 0
    integer(i32) :: escaped = 0
    integer(i32) :: escaped_boundary = 0
    integer(i32) :: survived_max_step = 0
    integer(i32) :: batches = 0
    real(dp) :: last_rel_change = -1.0d0
    real(dp) :: field_time_s = 0.0d0
    real(dp) :: push_time_s = 0.0d0
    real(dp) :: collision_time_s = 0.0d0
  end type sim_stats

  !> 種ごとのマクロ粒子端数を保持し、再開時にも注入期待値を保つ。
  type :: injection_state
    real(dp), allocatable :: macro_residual(:)
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
    real(dp), allocatable :: q_elem(:)
    real(dp) :: grid_bb_min(3) = 0.0d0
    real(dp) :: grid_bb_max(3) = 0.0d0
    integer(i32) :: grid_ncell(3) = 1_i32
    real(dp) :: grid_inv_cell(3) = 1.0d0
    integer(i32), allocatable :: grid_cell_start(:)
    integer(i32), allocatable :: grid_cell_elem(:)
    logical :: use_collision_grid = .false.
  end type mesh_type

  !> 粒子の位置・速度・物性値・生存フラグをSoA形式で保持する型。
  type :: particles_soa
    integer(i32) :: n = 0
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: v(:, :)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: m(:)
    real(dp), allocatable :: w(:)
    logical, allocatable :: alive(:)
  end type particles_soa

  !> 衝突有無、命中要素番号、線分パラメータ、交点座標を返す衝突結果型。
  type :: hit_info
    logical :: has_hit = .false.
    integer(i32) :: elem_idx = -1
    real(dp) :: t = 0.0d0
    real(dp) :: pos(3) = 0.0d0
  end type hit_info

end module bem_types
