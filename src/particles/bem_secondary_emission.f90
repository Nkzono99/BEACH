!> 吸収に伴う SEE の収率・丸めと、OpenMP thread ごとの生成粒子 buffer。
module bem_secondary_emission
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_counter_rng, only: counter_uniform
  implicit none
  private
  public :: secondary_emission_yield, secondary_emission_count
  public :: secondary_spawn_type, secondary_spawn_buffer_type

  type :: secondary_spawn_type
    real(dp) :: x(3), v(3), w
    integer(i32) :: species_idx, element, parent, ordinal
  end type secondary_spawn_type

  type :: secondary_spawn_buffer_type
    type(secondary_spawn_type), allocatable :: particles(:)
    integer(i32) :: n = 0_i32
  contains
    procedure :: push
  end type secondary_spawn_buffer_type
contains

  pure real(dp) function secondary_emission_yield(model, energy_ev, cos_theta, yield_max, energy_max_ev, const_yield) &
    result(delta)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: energy_ev, cos_theta, yield_max, energy_max_ev, const_yield
    real(dp) :: cosine, ratio, exponent

    delta = 0.0_dp
    if (model == 'constant') then
      delta = const_yield
    else if (energy_ev > 0.0_dp) then
      cosine = max(0.1_dp, min(1.0_dp, cos_theta))
      ratio = energy_ev/energy_max_ev
      exponent = 2.28_dp*cosine*ratio**1.35_dp
      ! Avoid subtracting two nearly equal numbers at low incident energies.
      if (exponent < 1.0e-5_dp) then
        delta = exponent*(1.0_dp - 0.5_dp*exponent + exponent**2/6.0_dp)
      else
        delta = 1.0_dp - exp(-exponent)
      end if
      delta = 1.114_dp*yield_max/cosine*ratio**(-0.35_dp)*delta
    end if
  end function secondary_emission_yield

  integer(i32) function secondary_emission_count(delta, charge_ratio, seed, counter) result(n)
    real(dp), intent(in) :: delta, charge_ratio
    integer(i32), intent(in) :: seed
    integer(i64), intent(in) :: counter(:)
    real(dp) :: expected

    expected = delta*charge_ratio
    if (.not. ieee_is_finite(expected) .or. expected >= real(huge(n), dp)) then
      error stop 'secondary emission particle count overflowed.'
    end if
    n = int(expected, i32)
    if (counter_uniform(seed, counter, 0_i64) < expected - real(n, dp)) n = n + 1_i32
  end function secondary_emission_count

  subroutine push(self, particle)
    class(secondary_spawn_buffer_type), intent(inout) :: self
    type(secondary_spawn_type), intent(in) :: particle
    type(secondary_spawn_type), allocatable :: grown(:)
    integer :: capacity

    capacity = 0
    if (allocated(self%particles)) capacity = size(self%particles)
    if (self%n == capacity) then
      allocate (grown(max(16, 2*capacity)))
      if (self%n > 0) grown(:self%n) = self%particles(:self%n)
      call move_alloc(grown, self%particles)
    end if
    self%n = self%n + 1_i32
    self%particles(self%n) = particle
  end subroutine push
end module bem_secondary_emission
