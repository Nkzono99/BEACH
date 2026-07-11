module bem_outer_plasma_photoelectron
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  implicit none
  private

  integer(i32), parameter, public :: photoelectron_closure_ok = 0_i32
  integer(i32), parameter, public :: photoelectron_closure_not_applicable = 1_i32
  integer(i32), parameter, public :: photoelectron_closure_invalid = 2_i32

  type, public :: photoelectron_histogram_type
    integer(i32) :: nbins = 0_i32
    real(dp) :: energy_max = 0.0_dp
    real(dp), allocatable :: signed_charge(:)
    real(dp), allocatable :: kinetic_energy(:)
    real(dp), allocatable :: tangential_momentum_x(:)
    real(dp), allocatable :: tangential_momentum_y(:)
    integer(i64), allocatable :: count(:)
  contains
    procedure :: init => init_photoelectron_histogram
    procedure :: reset => reset_photoelectron_histogram
    procedure :: add => add_photoelectron_crossing
    procedure :: merge => merge_photoelectron_histogram
    procedure :: total_signed_charge => total_photoelectron_signed_charge
    procedure :: total_kinetic_energy => total_photoelectron_kinetic_energy
    procedure :: total_count => total_photoelectron_count
  end type photoelectron_histogram_type

  type, public :: photoelectron_coupling_state_type
    logical :: ready = .false.
    integer(i32) :: last_completed_batch = 0_i32
    type(photoelectron_histogram_type) :: previous_batch
    type(photoelectron_histogram_type) :: cumulative
  contains
    procedure :: init => init_photoelectron_coupling_state
    procedure :: begin_batch => begin_photoelectron_batch
    procedure :: commit_batch => commit_photoelectron_batch
  end type photoelectron_coupling_state_type

  public :: validate_photoelectron_linear_applicability

contains

  subroutine init_photoelectron_histogram(self, nbins, energy_max)
    class(photoelectron_histogram_type), intent(out) :: self
    integer(i32), intent(in) :: nbins
    real(dp), intent(in) :: energy_max

    if (nbins < 1_i32 .or. .not. ieee_is_finite(energy_max) .or. energy_max <= 0.0_dp) then
      error stop 'Photoelectron histogram requires positive bins and energy_max.'
    end if
    self%nbins = nbins
    self%energy_max = energy_max
    allocate ( &
      self%signed_charge(nbins), self%kinetic_energy(nbins), self%tangential_momentum_x(nbins), &
      self%tangential_momentum_y(nbins), self%count(nbins) &
      )
    call self%reset()
  end subroutine init_photoelectron_histogram

  subroutine reset_photoelectron_histogram(self)
    class(photoelectron_histogram_type), intent(inout) :: self

    if (.not. allocated(self%signed_charge)) return
    self%signed_charge = 0.0_dp
    self%kinetic_energy = 0.0_dp
    self%tangential_momentum_x = 0.0_dp
    self%tangential_momentum_y = 0.0_dp
    self%count = 0_i64
  end subroutine reset_photoelectron_histogram

  subroutine add_photoelectron_crossing(self, charge, mass, weight, velocity)
    class(photoelectron_histogram_type), intent(inout) :: self
    real(dp), intent(in) :: charge, mass, weight, velocity(3)
    real(dp) :: normal_energy, total_energy, macro_mass
    integer(i32) :: bin

    if (self%nbins < 1_i32 .or. mass <= 0.0_dp .or. weight <= 0.0_dp .or. velocity(3) <= 0.0_dp .or. &
        .not. all(ieee_is_finite([charge, mass, weight, velocity]))) then
      error stop 'Invalid photoelectron interface crossing.'
    end if
    normal_energy = 0.5_dp*mass*velocity(3)*velocity(3)
    bin = min(self%nbins, 1_i32 + int(normal_energy/self%energy_max*real(self%nbins, dp), i32))
    macro_mass = mass*weight
    total_energy = 0.5_dp*macro_mass*sum(velocity*velocity)
    self%signed_charge(bin) = self%signed_charge(bin) + charge*weight
    self%kinetic_energy(bin) = self%kinetic_energy(bin) + total_energy
    self%tangential_momentum_x(bin) = self%tangential_momentum_x(bin) + macro_mass*velocity(1)
    self%tangential_momentum_y(bin) = self%tangential_momentum_y(bin) + macro_mass*velocity(2)
    self%count(bin) = self%count(bin) + 1_i64
  end subroutine add_photoelectron_crossing

  subroutine merge_photoelectron_histogram(self, other)
    class(photoelectron_histogram_type), intent(inout) :: self
    type(photoelectron_histogram_type), intent(in) :: other

    if (self%nbins /= other%nbins .or. self%energy_max /= other%energy_max) then
      error stop 'Photoelectron histogram contracts do not match.'
    end if
    self%signed_charge = self%signed_charge + other%signed_charge
    self%kinetic_energy = self%kinetic_energy + other%kinetic_energy
    self%tangential_momentum_x = self%tangential_momentum_x + other%tangential_momentum_x
    self%tangential_momentum_y = self%tangential_momentum_y + other%tangential_momentum_y
    self%count = self%count + other%count
  end subroutine merge_photoelectron_histogram

  pure real(dp) function total_photoelectron_signed_charge(self) result(total)
    class(photoelectron_histogram_type), intent(in) :: self

    total = sum(self%signed_charge)
  end function total_photoelectron_signed_charge

  pure real(dp) function total_photoelectron_kinetic_energy(self) result(total)
    class(photoelectron_histogram_type), intent(in) :: self

    total = sum(self%kinetic_energy)
  end function total_photoelectron_kinetic_energy

  pure integer(i64) function total_photoelectron_count(self) result(total)
    class(photoelectron_histogram_type), intent(in) :: self

    total = sum(self%count)
  end function total_photoelectron_count

  subroutine init_photoelectron_coupling_state(self, nbins, energy_max)
    class(photoelectron_coupling_state_type), intent(out) :: self
    integer(i32), intent(in) :: nbins
    real(dp), intent(in) :: energy_max

    call self%previous_batch%init(nbins, energy_max)
    call self%cumulative%init(nbins, energy_max)
    self%last_completed_batch = 0_i32
    self%ready = .true.
  end subroutine init_photoelectron_coupling_state

  subroutine begin_photoelectron_batch(self, batch_histogram)
    class(photoelectron_coupling_state_type), intent(in) :: self
    type(photoelectron_histogram_type), intent(out) :: batch_histogram

    if (.not. self%ready) error stop 'Photoelectron coupling state is not initialized.'
    call batch_histogram%init(self%previous_batch%nbins, self%previous_batch%energy_max)
  end subroutine begin_photoelectron_batch

  subroutine commit_photoelectron_batch(self, batch_index, batch_histogram)
    class(photoelectron_coupling_state_type), intent(inout) :: self
    integer(i32), intent(in) :: batch_index
    type(photoelectron_histogram_type), intent(in) :: batch_histogram

    if (.not. self%ready .or. batch_index /= self%last_completed_batch + 1_i32) then
      error stop 'Photoelectron batch sequence is not contiguous.'
    end if
    self%previous_batch = batch_histogram
    call self%cumulative%merge(batch_histogram)
    self%last_completed_batch = batch_index
  end subroutine commit_photoelectron_batch

  subroutine validate_photoelectron_linear_applicability(photoelectron_charge, ambient_charge_scale, max_ratio, status)
    real(dp), intent(in) :: photoelectron_charge, ambient_charge_scale, max_ratio
    integer(i32), intent(out) :: status

    status = photoelectron_closure_invalid
    if (.not. all(ieee_is_finite([photoelectron_charge, ambient_charge_scale, max_ratio])) .or. &
        ambient_charge_scale <= 0.0_dp .or. max_ratio <= 0.0_dp) return
    if (abs(photoelectron_charge)/ambient_charge_scale > max_ratio) then
      status = photoelectron_closure_not_applicable
    else
      status = photoelectron_closure_ok
    end if
  end subroutine validate_photoelectron_linear_applicability

end module bem_outer_plasma_photoelectron
