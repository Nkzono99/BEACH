!> 粒子追跡の失敗 status code と、その名前。
!! collision query（0–4）と粒子 step（1001–）は result の status と失敗 receipt に出る同じ名前空間で、
!! box 境界 event（0–1）は stepper の内部で particle_step_invalid_boundary へ写される。値はログと出力に出るので変えない。
module bem_status_codes
  use bem_kinds, only: i32
  implicit none
  private

  integer(i32), parameter, public :: collision_query_ok = 0_i32
  integer(i32), parameter, public :: collision_query_image_limit = 1_i32
  integer(i32), parameter, public :: collision_query_index_range = 2_i32
  integer(i32), parameter, public :: collision_query_invalid_segment = 3_i32
  integer(i32), parameter, public :: collision_query_grid_stalled = 4_i32

  integer(i32), parameter, public :: boundary_event_ok = 0_i32
  integer(i32), parameter, public :: boundary_event_invalid_geometry = 1_i32

  integer(i32), parameter, public :: particle_step_ok = collision_query_ok
  integer(i32), parameter, public :: particle_step_invalid_boundary = 1001_i32
  integer(i32), parameter, public :: particle_step_multiple_box_events = 1002_i32
  integer(i32), parameter, public :: particle_step_ambiguous_open_corner = 1003_i32

  public :: particle_status_name

contains

  !> collision query または粒子 step の失敗 status の名前。未知の値は 'unknown'。
  pure function particle_status_name(status) result(name)
    integer(i32), intent(in) :: status
    character(len=32) :: name

    select case (status)
    case (collision_query_image_limit)
      name = 'image_limit'
    case (collision_query_index_range)
      name = 'index_range'
    case (collision_query_invalid_segment)
      name = 'invalid_segment'
    case (collision_query_grid_stalled)
      name = 'grid_stalled'
    case (particle_step_invalid_boundary)
      name = 'invalid_boundary'
    case (particle_step_multiple_box_events)
      name = 'multiple_box_events'
    case (particle_step_ambiguous_open_corner)
      name = 'ambiguous_open_corner'
    case default
      name = 'unknown'
    end select
  end function particle_status_name

end module bem_status_codes
