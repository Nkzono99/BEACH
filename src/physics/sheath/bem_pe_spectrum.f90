!> Matching-plane PE の外向き数流束を法線エネルギーの区分一定分布として保持する。
module bem_pe_spectrum
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  implicit none
  private

  real(dp), parameter :: log_ten = log(10.0_dp)

  ! flux(j) is the integrated number flux in [edge(j-1), edge(j)], not dGamma/dK.
  ! Energies are normal kinetic energies in eV; potentials below are in volts.
  type, public :: pe_spectrum_type
    real(dp) :: energy_scale_ev = 1.0_dp
    integer(i32) :: bins_per_decade = 32_i32
    real(dp), allocatable :: flux(:)
  contains
    procedure :: clear => spectrum_clear
    procedure :: reset => spectrum_reset
    procedure :: add => spectrum_add
    procedure :: edge => spectrum_edge
    procedure :: total_flux => spectrum_total_flux
    procedure :: mean_energy => spectrum_mean_energy
    procedure :: tail_flux => spectrum_tail_flux
    procedure :: bootstrap_maxwellian => spectrum_bootstrap_maxwellian
    procedure :: combine => spectrum_combine
    procedure :: density => spectrum_density
    procedure :: integrated_density => spectrum_integrated_density
  end type pe_spectrum_type

contains

  subroutine spectrum_clear(self)
    class(pe_spectrum_type), intent(inout) :: self
    if (allocated(self%flux)) deallocate (self%flux)
  end subroutine spectrum_clear

  subroutine spectrum_reset(self)
    class(pe_spectrum_type), intent(inout) :: self
    if (allocated(self%flux)) self%flux = 0.0_dp
  end subroutine spectrum_reset

  subroutine spectrum_add(self, energy_ev, weight)
    class(pe_spectrum_type), intent(inout) :: self
    real(dp), intent(in) :: energy_ev, weight
    integer :: j

    call validate_grid(self)
    if (.not. ieee_is_finite(energy_ev) .or. energy_ev < 0.0_dp .or. &
        .not. ieee_is_finite(weight) .or. weight < 0.0_dp) &
      error stop 'PE spectrum requires finite nonnegative energies and flux weights.'
    if (weight == 0.0_dp) return
    j = energy_bin(self, energy_ev)
    call extend_flux(self, j)
    if (weight > huge(weight) - self%flux(j)) error stop 'PE spectrum bin flux overflow.'
    self%flux(j) = self%flux(j) + weight
  end subroutine spectrum_add

  pure real(dp) function spectrum_edge(self, j) result(value)
    class(pe_spectrum_type), intent(in) :: self
    integer, intent(in) :: j
    real(dp) :: exponent, factor

    exponent = log_ten*real(j, dp)/real(self%bins_per_decade, dp)
    if (exponent < log(huge(value))) then
      factor = exp_minus_one(exponent)
      if (factor > 1.0_dp) then
        if (self%energy_scale_ev > huge(value)/factor) error stop 'PE spectrum energy edge overflow.'
      end if
      value = self%energy_scale_ev*factor
    else
      ! Do not overflow exp(exponent) when a small scale still gives a finite edge.
      exponent = log(self%energy_scale_ev) + exponent
      if (exponent >= log(huge(value))) error stop 'PE spectrum energy edge overflow.'
      value = exp(exponent) - self%energy_scale_ev
    end if
  end function spectrum_edge

  pure real(dp) function spectrum_total_flux(self) result(value)
    class(pe_spectrum_type), intent(in) :: self
    value = 0.0_dp
    if (allocated(self%flux)) value = sum(self%flux)
  end function spectrum_total_flux

  pure real(dp) function spectrum_mean_energy(self) result(value)
    class(pe_spectrum_type), intent(in) :: self
    real(dp) :: total, left, right
    integer :: j

    value = 0.0_dp
    total = self%total_flux()
    if (total <= 0.0_dp) return
    left = 0.0_dp
    do j = 1, size(self%flux)
      right = self%edge(j)
      value = value + (self%flux(j)/total)*(0.5_dp*left + 0.5_dp*right)
      left = right
    end do
  end function spectrum_mean_energy

  pure real(dp) function spectrum_tail_flux(self, barrier_ev) result(value)
    class(pe_spectrum_type), intent(in) :: self
    real(dp), intent(in) :: barrier_ev
    real(dp) :: left, right
    integer :: j

    value = 0.0_dp
    if (.not. allocated(self%flux)) return
    left = 0.0_dp
    do j = 1, size(self%flux)
      right = self%edge(j)
      if (right > barrier_ev) value = value + self%flux(j)*(right - max(left, barrier_ev))/(right - left)
      left = right
    end do
  end function spectrum_tail_flux

  subroutine spectrum_bootstrap_maxwellian(self, total_flux, temperature_ev)
    class(pe_spectrum_type), intent(inout) :: self
    real(dp), intent(in) :: total_flux, temperature_ev
    real(dp) :: left, right
    integer :: j, n

    call validate_grid(self)
    if (.not. ieee_is_finite(total_flux) .or. total_flux < 0.0_dp .or. &
        .not. ieee_is_finite(temperature_ev) .or. temperature_ev <= 0.0_dp) &
      error stop 'Maxwell PE spectrum requires finite nonnegative flux and positive temperature.'
    call self%clear()
    if (total_flux == 0.0_dp) return
    if (temperature_ev > huge(temperature_ev)/40.0_dp) error stop 'Maxwell PE spectrum energy overflow.'
    n = energy_bin(self, 40.0_dp*temperature_ev)
    allocate (self%flux(n))
    left = 0.0_dp
    do j = 1, n
      right = self%edge(j)
      self%flux(j) = total_flux*exp(-left/temperature_ev)*(-exp_minus_one(-(right - left)/temperature_ev))
      left = right
    end do
  end subroutine spectrum_bootstrap_maxwellian

  subroutine spectrum_combine(self, other, factor_self, factor_other)
    class(pe_spectrum_type), intent(inout) :: self
    type(pe_spectrum_type), intent(in) :: other
    real(dp), intent(in) :: factor_self, factor_other
    integer :: n_other

    call validate_grid(self)
    call validate_grid(other)
    if (self%energy_scale_ev /= other%energy_scale_ev .or. self%bins_per_decade /= other%bins_per_decade) &
      error stop 'PE spectra cannot be combined on different energy grids.'
    if (.not. ieee_is_finite(factor_self) .or. factor_self < 0.0_dp .or. &
        .not. ieee_is_finite(factor_other) .or. factor_other < 0.0_dp) &
      error stop 'PE spectrum combination requires finite nonnegative factors.'
    n_other = 0
    if (allocated(other%flux)) n_other = size(other%flux)
    call extend_flux(self, n_other)
    if (allocated(self%flux)) self%flux = factor_self*self%flux
    if (n_other > 0) self%flux(:n_other) = self%flux(:n_other) + factor_other*other%flux
  end subroutine spectrum_combine

  ! n_free includes one outgoing escaping leg; n_returning includes BOTH legs
  ! of each reflected orbit. Their sum is the total PE density at this point.
  pure subroutine spectrum_density(self, phi_v, phi_h_v, phi_min_v, electron_mass_kg, n_free, n_returning, upper_side)
    class(pe_spectrum_type), intent(in) :: self
    real(dp), intent(in) :: phi_v, phi_h_v, phi_min_v, electron_mass_kg
    real(dp), intent(out) :: n_free, n_returning
    logical, intent(in), optional :: upper_side
    real(dp) :: b, barrier, coefficient

    b = phi_h_v - phi_v
    barrier = max(phi_h_v - phi_min_v, 0.0_dp)
    coefficient = sqrt(electron_mass_kg/(2.0_dp*qe))
    n_free = coefficient*power_moment(self, b, barrier, huge(b), .false.)
    n_returning = 0.0_dp
    if (present(upper_side)) then
      if (upper_side) return
    end if
    n_returning = 2.0_dp*coefficient*power_moment(self, b, 0.0_dp, barrier, .false.)
  end subroutine spectrum_density

  ! Exact integral of (n_free+n_returning) dphi from phi_min_v to phi_v, in m^-3 V.
  pure real(dp) function spectrum_integrated_density( &
    self, phi_v, phi_h_v, phi_min_v, electron_mass_kg, upper_side) result(value)
    class(pe_spectrum_type), intent(in) :: self
    real(dp), intent(in) :: phi_v, phi_h_v, phi_min_v, electron_mass_kg
    logical, intent(in), optional :: upper_side
    real(dp) :: b, b_min, barrier, coefficient

    b = phi_h_v - phi_v
    b_min = phi_h_v - phi_min_v
    barrier = max(b_min, 0.0_dp)
    coefficient = sqrt(2.0_dp*electron_mass_kg/qe)
    value = coefficient*sqrt_moment_change(self, b, b_min, barrier)
    if (present(upper_side)) then
      if (upper_side) return
    end if
    value = value + 2.0_dp*coefficient*power_moment(self, b, 0.0_dp, barrier, .true.)
  end function spectrum_integrated_density

  pure real(dp) function power_moment(self, b, lower, upper, square_root) result(value)
    class(pe_spectrum_type), intent(in) :: self
    real(dp), intent(in) :: b, lower, upper
    logical, intent(in) :: square_root
    real(dp) :: left, right, lo, hi, root_lo, root_hi, primitive
    integer :: j

    value = 0.0_dp
    if (.not. allocated(self%flux)) return
    left = 0.0_dp
    do j = 1, size(self%flux)
      right = self%edge(j)
      lo = max(left, lower, b)
      hi = min(right, upper)
      if (hi > lo .and. self%flux(j) > 0.0_dp) then
        root_lo = sqrt(max(lo - b, 0.0_dp))
        root_hi = sqrt(max(hi - b, 0.0_dp))
        primitive = 2.0_dp*(hi - lo)/(root_hi + root_lo)
        if (square_root) primitive = primitive*((hi - b) + root_hi*root_lo + (lo - b))/3.0_dp
        value = value + self%flux(j)*(primitive/(right - left))
      end if
      left = right
      if (left >= upper) exit
    end do
  end function power_moment

  ! Difference of sqrt moments evaluated together. Factoring BOTH the potential
  ! difference and bin width keeps near-minimum Sagdeev integrals accurate.
  pure real(dp) function sqrt_moment_change(self, b, b_reference, lower) result(value)
    class(pe_spectrum_type), intent(in) :: self
    real(dp), intent(in) :: b, b_reference, lower
    real(dp) :: near_b, far_b, delta, left, right, lo, hi, primitive
    real(dp) :: near_lo, near_hi, far_lo, far_hi, sum_lo, sum_hi, r_change
    integer :: j

    value = 0.0_dp
    if (.not. allocated(self%flux) .or. b == b_reference) return
    near_b = min(b, b_reference)
    far_b = max(b, b_reference)
    delta = far_b - near_b
    left = 0.0_dp
    do j = 1, size(self%flux)
      right = self%edge(j)
      if (self%flux(j) > 0.0_dp) then
        ! Energies below far_b contribute only to the nearer square root.
        lo = max(left, lower, near_b)
        hi = min(right, far_b)
        primitive = 0.0_dp
        if (hi > lo) then
          near_lo = sqrt(max(lo - near_b, 0.0_dp))
          near_hi = sqrt(max(hi - near_b, 0.0_dp))
          primitive = (2.0_dp/3.0_dp)*(hi - lo)* &
                      ((hi - near_b) + near_hi*near_lo + (lo - near_b))/(near_hi + near_lo)
        end if
        lo = max(left, lower, far_b)
        hi = right
        if (hi > lo) then
          near_lo = sqrt(max(lo - near_b, 0.0_dp))
          near_hi = sqrt(max(hi - near_b, 0.0_dp))
          far_lo = sqrt(max(lo - far_b, 0.0_dp))
          far_hi = sqrt(max(hi - far_b, 0.0_dp))
          sum_lo = near_lo + far_lo
          sum_hi = near_hi + far_hi
          r_change = (hi - lo)/(near_hi + near_lo)*(1.0_dp - (far_hi/sum_hi)*(far_lo/sum_lo)) + &
                     (hi - lo)/(far_hi + far_lo)*(1.0_dp - (near_hi/sum_hi)*(near_lo/sum_lo))
          primitive = primitive + (2.0_dp/3.0_dp)*delta*r_change
        end if
        value = value + self%flux(j)*(primitive/(right - left))
      end if
      left = right
    end do
    if (b > b_reference) value = -value
  end function sqrt_moment_change

  subroutine validate_grid(self)
    class(pe_spectrum_type), intent(in) :: self
    if (.not. ieee_is_finite(self%energy_scale_ev) .or. self%energy_scale_ev <= 0.0_dp .or. &
        self%bins_per_decade <= 0_i32) error stop 'PE spectrum requires a finite positive energy scale and bin count.'
  end subroutine validate_grid

  integer function energy_bin(self, energy_ev) result(j)
    class(pe_spectrum_type), intent(in) :: self
    real(dp), intent(in) :: energy_ev
    real(dp) :: coordinate

    if (energy_ev > self%energy_scale_ev) then
      coordinate = log(energy_ev) - log(self%energy_scale_ev) + log_one_plus(self%energy_scale_ev/energy_ev)
    else
      coordinate = log_one_plus(energy_ev/self%energy_scale_ev)
    end if
    coordinate = coordinate*real(self%bins_per_decade, dp)/log_ten
    if (coordinate >= real(huge(j) - 1, dp)) error stop 'PE spectrum energy grid size is not representable.'
    j = int(coordinate) + 1
    ! Keep exact-edge samples deterministic despite log/exp rounding.
    do while (energy_ev >= self%edge(j))
      j = j + 1
    end do
    do while (j > 1)
      if (energy_ev >= self%edge(j - 1)) exit
      j = j - 1
    end do
  end function energy_bin

  subroutine extend_flux(self, n)
    class(pe_spectrum_type), intent(inout) :: self
    integer, intent(in) :: n
    real(dp), allocatable :: extended(:)
    integer :: old_n

    if (n <= 0) return
    old_n = 0
    if (allocated(self%flux)) old_n = size(self%flux)
    if (old_n >= n) return
    allocate (extended(n), source=0.0_dp)
    if (old_n > 0) extended(:old_n) = self%flux
    call move_alloc(extended, self%flux)
  end subroutine extend_flux

  pure real(dp) function log_one_plus(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: rounded
    rounded = 1.0_dp + x
    if (rounded == 1.0_dp) then
      value = x
    else
      value = log(rounded)*(x/(rounded - 1.0_dp))
    end if
  end function log_one_plus

  pure real(dp) function exp_minus_one(x) result(value)
    real(dp), intent(in) :: x
    if (abs(x) < 0.5_dp) then
      value = 2.0_dp*exp(0.5_dp*x)*sinh(0.5_dp*x)
    else
      value = exp(x) - 1.0_dp
    end if
  end function exp_minus_one

end module bem_pe_spectrum
