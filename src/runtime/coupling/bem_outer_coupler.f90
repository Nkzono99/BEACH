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
    logical :: implicit_mean_reuse = .false.
    logical :: snapshot_matches_mesh = .false.
  contains
    procedure :: init => init_outer_coupler
    procedure :: mark_mesh_changed => mark_outer_coupler_mesh_changed
    procedure :: accept_external_update => accept_outer_coupler_external_update
    procedure :: accept_restored_snapshot => accept_outer_coupler_restored_snapshot
    procedure :: refresh => refresh_outer_coupler
  end type outer_coupler_type

contains

  subroutine init_outer_coupler(self, config, last_outer_update_batch)
    class(outer_coupler_type), intent(out) :: self
    type(coupling_config), intent(in) :: config
    integer(i32), intent(in), optional :: last_outer_update_batch

    select case (trim(lower_ascii(config%update_mode)))
    case ('explicit')
      self%implicit_mean_reuse = .false.
    case ('implicit_mean')
      self%implicit_mean_reuse = .true.
    case default
      error stop 'outer coupler received an unsupported update_mode.'
    end select
    if (trim(lower_ascii(config%particle_transfer_mode)) /= 'none' .and. &
        trim(lower_ascii(config%particle_transfer_mode)) /= 'electrostatic_1d_instant_return') then
      error stop 'Unknown outer-coupler particle transfer mode.'
    end if
    if (config%outer_update_stride < 1_i32) error stop 'outer_update_stride must be >= 1.'
    self%update_stride = config%outer_update_stride
    self%last_outer_update_batch = -1_i32
    self%snapshot_matches_mesh = .false.
    if (present(last_outer_update_batch)) self%last_outer_update_batch = last_outer_update_batch
  end subroutine init_outer_coupler

  !> Simulatorがmesh chargeを更新したことを通知し、次のrefreshを省略不可にする。
  subroutine mark_outer_coupler_mesh_changed(self)
    class(outer_coupler_type), intent(inout) :: self

    self%snapshot_matches_mesh = .false.
  end subroutine mark_outer_coupler_mesh_changed

  !> Simulator側でmeshとouter stateを同じbatch時刻へ更新したことを記録する。
  subroutine accept_outer_coupler_external_update(self, batch_index)
    class(outer_coupler_type), intent(inout) :: self
    integer(i32), intent(in) :: batch_index

    if (batch_index < 1_i32) then
      error stop 'outer coupler external update batch index must be >= 1.'
    end if
    self%last_outer_update_batch = batch_index
    self%snapshot_matches_mesh = .true.
  end subroutine accept_outer_coupler_external_update

  !> checkpointから復元したouter stateと現meshが同じ時刻にあることを記録する。
  subroutine accept_outer_coupler_restored_snapshot(self)
    class(outer_coupler_type), intent(inout) :: self

    self%snapshot_matches_mesh = .true.
  end subroutine accept_outer_coupler_restored_snapshot

  subroutine refresh_outer_coupler( &
    self, snapshot, mesh, batch_index, outer_updated, continuation_stage, force_outer_update &
    )
    class(outer_coupler_type), intent(inout) :: self
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: batch_index
    logical, intent(out), optional :: outer_updated
    character(len=*), intent(in), optional :: continuation_stage
    logical, intent(in), optional :: force_outer_update
    logical :: update_now, force_update, reuse_pre_batch

    if (batch_index < 1_i32) error stop 'outer coupler batch index must be >= 1.'
    force_update = .false.
    if (present(force_outer_update)) force_update = force_outer_update
    reuse_pre_batch = .false.
    if (present(continuation_stage)) then
      reuse_pre_batch = self%implicit_mean_reuse .and. self%snapshot_matches_mesh .and. &
                        trim(lower_ascii(continuation_stage)) == 'pre_batch' .and. .not. force_update
    end if
    if (reuse_pre_batch) then
      if (present(outer_updated)) outer_updated = .false.
      return
    end if

    update_now = self%last_outer_update_batch < 0_i32 .or. &
                 mod(batch_index - 1_i32, self%update_stride) == 0_i32
    update_now = update_now .or. force_update
    call snapshot%refresh( &
      mesh, update_outer=update_now, continuation_stage=continuation_stage, &
      continuation_batch=batch_index &
      )
    self%snapshot_matches_mesh = .true.
    if (update_now) self%last_outer_update_batch = batch_index
    if (present(outer_updated)) outer_updated = update_now
  end subroutine refresh_outer_coupler

end module bem_outer_coupler
