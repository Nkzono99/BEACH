module bem_outer_coupler
  use bem_kinds, only: i32
  use bem_types, only: mesh_type
  use bem_physics_config_types, only: coupling_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  type, public :: outer_coupler_type
    integer(i32) :: update_stride = 1_i32
    integer(i32) :: last_outer_update_batch = -1_i32
  contains
    procedure :: init => init_outer_coupler
    procedure :: refresh => refresh_outer_coupler
  end type outer_coupler_type

contains

  subroutine init_outer_coupler(self, config, last_outer_update_batch)
    class(outer_coupler_type), intent(out) :: self
    type(coupling_config), intent(in) :: config
    integer(i32), intent(in), optional :: last_outer_update_batch

    if (trim(lower_ascii(config%update_mode)) /= 'explicit') then
      error stop 'outer coupler currently requires update_mode=explicit.'
    end if
    if (trim(lower_ascii(config%particle_transfer_mode)) /= 'none' .and. &
        trim(lower_ascii(config%particle_transfer_mode)) /= 'electrostatic_1d_instant_return') then
      error stop 'Unknown outer-coupler particle transfer mode.'
    end if
    if (config%outer_update_stride < 1_i32) error stop 'outer_update_stride must be >= 1.'
    self%update_stride = config%outer_update_stride
    self%last_outer_update_batch = -1_i32
    if (present(last_outer_update_batch)) self%last_outer_update_batch = last_outer_update_batch
  end subroutine init_outer_coupler

  subroutine refresh_outer_coupler( &
    self, snapshot, mesh, batch_index, outer_updated, continuation_stage &
    )
    class(outer_coupler_type), intent(inout) :: self
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: batch_index
    logical, intent(out), optional :: outer_updated
    character(len=*), intent(in), optional :: continuation_stage
    logical :: update_now

    if (batch_index < 1_i32) error stop 'outer coupler batch index must be >= 1.'
    update_now = self%last_outer_update_batch < 0_i32 .or. &
                 mod(batch_index - 1_i32, self%update_stride) == 0_i32
    call snapshot%refresh( &
      mesh, update_outer=update_now, continuation_stage=continuation_stage, &
      continuation_batch=batch_index &
      )
    if (update_now) self%last_outer_update_batch = batch_index
    if (present(outer_updated)) outer_updated = update_now
  end subroutine refresh_outer_coupler

end module bem_outer_coupler
