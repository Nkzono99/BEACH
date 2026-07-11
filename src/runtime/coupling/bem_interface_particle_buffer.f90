module bem_interface_particle_buffer
  use bem_kinds, only: dp, i32
  use bem_interface_types, only: interface_crossing_type
  implicit none
  private

  type, public :: interface_particle_record_type
    type(interface_crossing_type) :: crossing
    real(dp) :: charge = 0.0_dp
    real(dp) :: mass = 0.0_dp
    real(dp) :: weight = 0.0_dp
    integer(i32) :: species_id = 0_i32
  end type interface_particle_record_type

  type, public :: interface_particle_buffer_type
    type(interface_particle_record_type), allocatable :: records(:)
    integer(i32) :: count = 0_i32
  contains
    procedure :: init => init_interface_particle_buffer
    procedure :: clear => clear_interface_particle_buffer
    procedure :: append => append_interface_particle_buffer
  end type interface_particle_buffer_type

contains

  subroutine init_interface_particle_buffer(self, capacity)
    class(interface_particle_buffer_type), intent(out) :: self
    integer(i32), intent(in), optional :: capacity
    integer(i32) :: initial_capacity

    initial_capacity = 0_i32
    if (present(capacity)) initial_capacity = max(0_i32, capacity)
    allocate (self%records(initial_capacity))
    self%count = 0_i32
  end subroutine init_interface_particle_buffer

  subroutine clear_interface_particle_buffer(self)
    class(interface_particle_buffer_type), intent(inout) :: self

    self%count = 0_i32
  end subroutine clear_interface_particle_buffer

  subroutine append_interface_particle_buffer(self, record)
    class(interface_particle_buffer_type), intent(inout) :: self
    type(interface_particle_record_type), intent(in) :: record
    type(interface_particle_record_type), allocatable :: grown(:)
    integer(i32) :: old_capacity, new_capacity

    if (.not. allocated(self%records)) allocate (self%records(0))
    old_capacity = int(size(self%records), i32)
    if (self%count == old_capacity) then
      new_capacity = max(1_i32, 2_i32*old_capacity)
      allocate (grown(new_capacity))
      if (self%count > 0_i32) grown(1:self%count) = self%records(1:self%count)
      call move_alloc(grown, self%records)
    end if
    self%count = self%count + 1_i32
    self%records(self%count) = record
  end subroutine append_interface_particle_buffer

end module bem_interface_particle_buffer
