module bem_types
  use bem_kinds, only: dp, i32
  implicit none

  type :: sim_config
    real(dp) :: dt = 1.0d-12
    integer(i32) :: npcls_per_step = 128
    integer(i32) :: max_step = 100
    real(dp) :: tol_rel = 1.0d-4
    real(dp) :: q_floor = 1.0d-30
    real(dp) :: softening = 1.0d-4
  end type sim_config

  type :: sim_stats
    integer(i32) :: processed_particles = 0
    integer(i32) :: absorbed = 0
    integer(i32) :: escaped = 0
    integer(i32) :: batches = 0
    real(dp) :: last_rel_change = -1.0d0
  end type sim_stats

  type :: mesh_type
    integer(i32) :: nelem = 0
    real(dp), allocatable :: v0(:, :)
    real(dp), allocatable :: v1(:, :)
    real(dp), allocatable :: v2(:, :)
    real(dp), allocatable :: centers(:, :)
    real(dp), allocatable :: normals(:, :)
    real(dp), allocatable :: bb_min(:, :)
    real(dp), allocatable :: bb_max(:, :)
    real(dp), allocatable :: h_elem(:)
    real(dp), allocatable :: q_elem(:)
  end type mesh_type

  type :: particles_soa
    integer(i32) :: n = 0
    real(dp), allocatable :: x(:, :)
    real(dp), allocatable :: v(:, :)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: m(:)
    real(dp), allocatable :: w(:)
    logical, allocatable :: alive(:)
  end type particles_soa

  type :: hit_info
    logical :: has_hit = .false.
    integer(i32) :: elem_idx = -1
    real(dp) :: t = 0.0d0
    real(dp) :: pos(3) = 0.0d0
  end type hit_info

end module bem_types
