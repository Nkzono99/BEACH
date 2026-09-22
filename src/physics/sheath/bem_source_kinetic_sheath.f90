!> Zero-drift source-connected, semi-infinite 1D sheath.
!! Units: Te/e, electron Debye length, ni(infinity), sqrt(Te/me).
!! tau=Tph/Te, unlike the reciprocal temperature ratio in the legacy Zhao model.
!! The ambient Maxwellian amplitude is determined by neutrality, not prescribed.
!! Current is an output; finite searches never certify uniqueness or completeness.
module bem_source_kinetic_sheath
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use bem_kinds, only: dp
  use bem_constants, only: pi
  implicit none
  private

  integer, parameter, public :: source_kinetic_ok = 0, source_kinetic_invalid = 1
  integer, parameter, public :: source_kinetic_excluded = 2, source_kinetic_unresolved = 3
  integer, parameter, public :: source_kinetic_numerical_failure = 4
  integer, parameter :: nq = 96

  type, public :: source_kinetic_root
    character(len=2) :: branch = ''
    real(dp) :: phi_h = 0, phi_min = 0, amplitude = 0, escaping_flux = 0
    real(dp) :: electron_flux = 0, ion_flux = 0, current = 0
    real(dp) :: neutrality_residual = 0, field_squared_residual = 0, outer_residual = 0
    real(dp) :: minimum_field_squared = 0, minimum_curvature = 0, edge_coefficient = 0
    real(dp) :: barrier_roundoff = 0, ion_turning_margin = 0
    logical :: accepted = .false., boundary_turning = .false.
    character(len=80) :: rejection = ''
  end type

  type, public :: source_kinetic_result
    type(source_kinetic_root), allocatable :: roots(:), rejected(:)
    integer :: status = source_kinetic_unresolved, search_points = 481
    integer :: numerical_failures = 0
    real(dp) :: depth_min = 1.e-24_dp, depth_max = 1.e8_dp
    logical :: search_complete = .false.
    character(len=160) :: message = ''
  end type

  type :: context_type
    real(dp) :: mach, tau, emission, field, mass_ratio
    real(dp) :: x(nq), w(nq)
  end type

  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function
  end interface

  public :: solve_source_kinetic, source_kinetic_density, source_kinetic_field_squared

contains

  subroutine solve_source_kinetic(mach, tau, emission, field, mass_ratio, result, search_points, maximum_depth)
    real(dp), intent(in) :: mach, tau, emission, field, mass_ratio
    type(source_kinetic_result), intent(out) :: result
    integer, intent(in), optional :: search_points
    real(dp), intent(in), optional :: maximum_depth
    type(context_type) :: ctx
    type(source_kinetic_root) :: root
    real(dp), allocatable :: grid(:), found(:), boundaries(:)
    real(dp) :: a, q, lo, hi, k, d, f, edge_flux
    integer :: points, i, mode

    allocate (result%roots(0), result%rejected(0))
    if (present(search_points)) result%search_points = search_points
    if (present(maximum_depth)) result%depth_max = maximum_depth
    result%status = source_kinetic_invalid
    if (.not. all(ieee_is_finite([mach, tau, emission, field, mass_ratio, result%depth_max]))) then
      result%message = 'source_kinetic inputs must be finite.'
      return
    end if
    if (min(mach, tau, mass_ratio, result%depth_max) <= 0 .or. emission < 0 .or. &
        result%search_points < 65 .or. result%depth_max <= result%depth_min) then
      result%message = 'source_kinetic requires positive Mach, temperatures, mass and depth, G>=0, search_points>=65.'
      return
    end if
    call make_context(mach, tau, emission, field, mass_ratio, ctx)
    k = sqrt(pi/(2*tau))
    if (.not. all(ieee_is_finite([mach**2, field**2, k, k*emission, 4*k*tau*emission]))) then
      result%status = source_kinetic_numerical_failure
      result%numerical_failures = 1
      result%message = 'Non-finite arithmetic in source_kinetic normalization.'
      return
    end if
    if (field > 0 .and. mach >= 1 .and. field**2 > nearest(4*k*tau*emission, 1._dp)) then
      result%status = source_kinetic_excluded
      result%search_complete = .true.
      result%message = 'Analytical exclusion: E_H^2 < 2*sqrt(2*pi*tau)*G is necessary for M>=1.'
      return
    end if
    if (field > 0 .and. emission == 0) then
      result%status = source_kinetic_excluded
      result%search_complete = .true.
      result%message = 'No emission: no interior minimum and no accessible positive kinetic-B edge.'
      return
    end if
    points = result%search_points
    if (field == 0 .and. k*emission < 1) then
      root = make_root(ctx, 'B0', 0._dp, 0._dp, 2*(1 - k*emission), emission)
      root%accepted = .true.
      result%roots = [result%roots, root]
    end if

    ! B: neutrality eliminates A, q=G exp(-phi_H/tau). Keep the marginal
    ! kinetic edge q=q_*; the ion-turning endpoint remains strictly excluded.
    if (field >= 0 .and. emission > 0) then
      edge_flux = sqrt(2/pi)*tau/(1 + sqrt(tau))
      lo = max(0._dp, tau*log(k*emission))
      hi = min(mach**2/2, tau*log(emission/edge_flux))
      if (hi > lo) then
        call bounded_grid(lo, hi, points, grid)
        if (hi < mach**2/2) grid = [grid, hi]
        mode = 1
        call find_roots(equation, grid, 1.e-11_dp*max(field**2, 1.e-30_dp), found)
        do i = 1, size(found)
          call state_at(found(i), root)
          call retain(root)
        end do
      end if
    end if

    call depth_grid(result%depth_min, result%depth_max, points, grid)
    if (field <= 0) then
      mode = 2
      call find_roots(equation, grid, 1.e-11_dp*max(field**2, 1.e-30_dp), found)
      do i = 1, size(found)
        call state_at(found(i), root)
        call retain(root)
      end do
    else if (emission > 0) then
      ! A and N share the outer closure; classify only after solving the
      ! inner field. Add physical-domain boundaries before scanning narrow
      ! intervals, without interpolating or filling missing responses.
      mode = 3
      call find_roots(domain_boundary, grid, 0._dp, boundaries)
      call add_boundary_neighbors(grid, boundaries)
      mode = 4
      call find_roots(domain_boundary, grid, 0._dp, boundaries)
      call add_boundary_neighbors(grid, boundaries)
      mode = 3
      call find_roots(equation, grid, 1.e-11_dp*max(field**2, 1.e-30_dp), found)
      do i = 1, size(found)
        call state_at(found(i), root)
        call retain(root)
      end do
    end if
    if (size(result%roots) > 0) then
      result%status = source_kinetic_ok
      result%message = 'All detected admissible roots retained; finite search does not certify completeness or stability.'
    else if (result%numerical_failures > 0) then
      result%status = source_kinetic_numerical_failure
      result%message = 'No admissible root retained; non-finite candidate arithmetic was encountered.'
    else
      result%status = source_kinetic_unresolved
      result%message = 'No admissible root detected in the recorded finite search; this is not a nonexistence certificate.'
    end if

  contains

    subroutine state_at(value, r)
      real(dp), intent(in) :: value
      type(source_kinetic_root), intent(out) :: r
      real(dp) :: ac, bc, ii, ie, ip, gap, wall
      r = source_kinetic_root()
      select case (mode)
      case (1)
        q = emission*exp(-value/tau)
        a = 2*(1 - k*q)
        r = make_root(ctx, 'B', value, 0._dp, a, q)
      case (2)
        call outer_basis(ctx, value, ac, bc, ii, ie, ip, q, a)
        a = (1 - bc*emission)/ac
        r = make_root(ctx, 'C', -value, -value, a, emission)
      case (3)
        call outer_basis(ctx, value, ac, bc, ii, ie, ip, q, a)
        if (min(q, a) <= 0 .or. .not. all(ieee_is_finite([q, a]))) return
        gap = tau*log_ratio(emission, q)
        wall = -value + gap
        if (gap <= 0 .or. wall >= mach**2/2) return
        if (wall > 0) then
          r = make_root(ctx, 'A', wall, -value, a, q)
        else
          r = make_root(ctx, 'N', wall, -value, a, q)
        end if
      end select
    end subroutine

    function equation(value) result(y)
      real(dp), intent(in) :: value
      real(dp) :: y, ac, bc, ii, ie, ip, aq, aa
      type(source_kinetic_root) :: r
      y = nan()
      call state_at(value, r)
      if (r%amplitude <= 0 .or. len_trim(r%branch) == 0) return
      if (mode == 2 .and. field == 0) then
        call outer_basis(ctx, value, ac, bc, ii, ie, ip, aq, aa)
        if (aq > 0 .and. aa > 0) y = aq - emission
      else
        y = endpoint_field_squared(ctx, r) - field**2
        if (.not. ieee_is_finite(y)) result%numerical_failures = result%numerical_failures + 1
      end if
    end function

    function domain_boundary(value) result(y)
      real(dp), intent(in) :: value
      real(dp) :: y, ac, bc, ii, ie, ip, aq, aa
      call outer_basis(ctx, value, ac, bc, ii, ie, ip, aq, aa)
      y = nan()
      if (min(aq, aa) <= 0 .or. .not. all(ieee_is_finite([aq, aa]))) return
      y = tau*log_ratio(emission, aq)
      if (mode == 4) y = y - value - mach**2/2
    end function

    subroutine retain(r)
      type(source_kinetic_root), intent(inout) :: r
      integer :: j
      call validate_root(ctx, r)
      if (r%rejection == 'non-finite candidate arithmetic') result%numerical_failures = result%numerical_failures + 1
      do j = 1, size(result%roots)
        d = max(abs(r%phi_h), abs(r%phi_min), 1.e-30_dp)
        f = max(abs(r%phi_h - result%roots(j)%phi_h), abs(r%phi_min - result%roots(j)%phi_min))
        if (r%branch == result%roots(j)%branch .and. f < 2.e-8_dp*d) return
      end do
      if (r%accepted) then
        result%roots = [result%roots, r]
      else
        result%rejected = [result%rejected, r]
      end if
    end subroutine
  end subroutine solve_source_kinetic

  subroutine make_context(mach, tau, emission, field, mass_ratio, ctx)
    real(dp), intent(in) :: mach, tau, emission, field, mass_ratio
    type(context_type), intent(out) :: ctx
    integer :: i, j, it
    real(dp) :: z, old, p1, p2, p3, derivative
    ctx%mach = mach
    ctx%tau = tau
    ctx%emission = emission
    ctx%field = field
    ctx%mass_ratio = mass_ratio
    do i = 1, nq/2
      z = cos(pi*(i - 0.25_dp)/(nq + 0.5_dp))
      do it = 1, 30
        p1 = 1
        p2 = 0
        do j = 1, nq
          p3 = p2
          p2 = p1
          p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
        end do
        derivative = nq*(z*p1 - p2)/(z*z - 1)
        old = z
        z = old - p1/derivative
        if (abs(z - old) <= 2*epsilon(z)) exit
      end do
      ctx%x(i) = (1 - z)/2
      ctx%x(nq + 1 - i) = (1 + z)/2
      ctx%w(i) = 1/((1 - z*z)*derivative**2)
      ctx%w(nq + 1 - i) = ctx%w(i)
    end do
  end subroutine

  function make_root(ctx, branch, wall, minimum, amplitude, q) result(root)
    type(context_type), intent(in) :: ctx
    character(len=*), intent(in) :: branch
    real(dp), intent(in) :: wall, minimum, amplitude, q
    type(source_kinetic_root) :: root
    root%branch = branch
    root%phi_h = wall
    root%phi_min = minimum
    root%amplitude = amplitude
    root%escaping_flux = q
    root%electron_flux = amplitude*exp(minimum)/sqrt(2*pi)
    root%ion_flux = ctx%mach/sqrt(ctx%mass_ratio)
    root%current = root%ion_flux + q - root%electron_flux
    root%ion_turning_margin = ctx%mach**2/2 - max(wall, 0._dp)
    root%boundary_turning = ctx%field == 0 .and. branch /= 'B0'
  end function

  ! These density and first-integral entry points support independent orbit
  ! quadrature tests and profile receipts; callers provide a validated root.
  subroutine source_kinetic_density(mach, tau, emission, root, phi, lower, ni, ne, npe)
    real(dp), intent(in) :: mach, tau, emission, phi
    type(source_kinetic_root), intent(in) :: root
    logical, intent(in) :: lower
    real(dp), intent(out) :: ni, ne, npe
    real(dp) :: s, sp, free
    ni = (1 - 2*phi/mach**2)**(-0.5_dp)
    s = sqrt(max(0._dp, phi - root%phi_min))
    sp = s/sqrt(tau)
    free = 0.5_dp*exp(root%phi_min)*erfc_scaled(s)
    if (root%branch == 'B' .or. root%branch == 'B0' .or. lower) then
      ne = root%amplitude*free
      npe = sqrt(pi/(2*tau))*(2*emission*exp((phi - root%phi_h)/tau) - root%escaping_flux*erfc_scaled(sp))
    else
      ne = root%amplitude*(exp(phi) - free)
      npe = sqrt(pi/(2*tau))*root%escaping_flux*erfc_scaled(sp)
    end if
  end subroutine

  elemental function erfcx_increment(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y, coefficients(0:14)
    integer :: i
    if (abs(x) >= 0.05_dp) then
      y = erfc_scaled(x) - 1
      return
    end if
    coefficients(0:1) = [1._dp, -2/sqrt(pi)]
    do i = 2, 14
      coefficients(i) = 2*coefficients(i - 2)/i
    end do
    y = coefficients(14)
    do i = 13, 1, -1
      y = y*x + coefficients(i)
    end do
    y = x*y
  end function

  elemental function exp_increment(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y, term
    integer :: i
    if (abs(x) > 0.05_dp) then
      y = exp(x) - 1
      return
    end if
    y = x
    term = x
    do i = 2, 12
      term = term*x/i
      y = y + term
    end do
  end function

  elemental function ion_increment(phi, mach) result(y)
    real(dp), intent(in) :: phi, mach
    real(dp) :: y, r
    r = sqrt(1 - 2*phi/mach**2)
    y = (2*phi/mach**2)/(r*(1 + r))
  end function

  function log_ratio(a, b) result(y)
    real(dp), intent(in) :: a, b
    real(dp) :: y, t, term
    integer :: i
    if (abs(a - b) >= 0.5_dp*b) then
      y = log(a) - log(b)
    else
      t = (a - b)/(a + b)
      y = t
      term = t
      do i = 3, 31, 2
        term = term*t*t
        y = y + term/i
      end do
      y = 2*y
    end if
  end function

  function charge(ctx, r, phi, lower) result(rho)
    type(context_type), intent(in) :: ctx
    type(source_kinetic_root), intent(in) :: r
    real(dp), intent(in) :: phi
    logical, intent(in) :: lower
    real(dp) :: rho, ni, ne, npe, s, s0, dec, dpc, k, increment
    k = sqrt(pi/(2*ctx%tau))
    if (r%branch == 'B' .or. r%branch == 'B0') then
      s = sqrt(max(phi, 0._dp))
      if (phi/ctx%tau < 0.5_dp) then
        increment = r%escaping_flux*exp_increment(phi/ctx%tau)
      else
        increment = ctx%emission*exp((phi - r%phi_h)/ctx%tau) - r%escaping_flux
      end if
      rho = ion_increment(phi, ctx%mach) - 0.5_dp*r%amplitude*erfcx_increment(s) &
            - k*(2*increment - r%escaping_flux*erfcx_increment(s/sqrt(ctx%tau)))
    else if (-r%phi_min < 1.e-3_dp*min(1._dp, ctx%tau) .and. abs(phi) < 1.e-3_dp*min(1._dp, ctx%tau)) then
      s = sqrt(max(phi - r%phi_min, 0._dp))
      s0 = sqrt(-r%phi_min)
      if (lower) then
        dec = 0.5_dp*(exp_increment(phi)*(1 - erf(s)) - erf(s) - erf(s0))
        dpc = k*(2*exp_increment((phi - r%phi_min)/ctx%tau) &
                 - erfcx_increment(s/sqrt(ctx%tau)) - erfcx_increment(s0/sqrt(ctx%tau)))
      else
        dec = 0.5_dp*(exp_increment(phi)*(1 + erf(s)) + erf(s) - erf(s0))
        dpc = k*(erfcx_increment(s/sqrt(ctx%tau)) - erfcx_increment(s0/sqrt(ctx%tau)))
      end if
      rho = ion_increment(phi, ctx%mach) - r%amplitude*dec - r%escaping_flux*dpc
    else
      call source_kinetic_density(ctx%mach, ctx%tau, ctx%emission, r, phi, lower, ni, ne, npe)
      rho = ni - ne - npe
    end if
  end function

  function ion_integral(a, b, mach) result(y)
    real(dp), intent(in) :: a, b, mach
    real(dp) :: y
    y = 2*(b - a)/(sqrt(1 - 2*a/mach**2) + sqrt(1 - 2*b/mach**2))
  end function

  subroutine outer_basis(ctx, depth, a, b, ii, ie, ip, q, amplitude)
    type(context_type), intent(in) :: ctx
    real(dp), intent(in) :: depth
    real(dp), intent(out) :: a, b, ii, ie, ip, q, amplitude
    real(dp) :: s, s0, phi, ec, pc, relative_ec, delta_pc, numerator, denominator, k, weight
    integer :: j
    k = sqrt(pi/(2*ctx%tau))
    s0 = sqrt(depth)
    a = 0.5_dp*(1 + erf(s0))
    b = k*erfc_scaled(s0/sqrt(ctx%tau))
    ii = ion_integral(-depth, 0._dp, ctx%mach)
    if (depth < 0.05_dp*min(1._dp, ctx%tau)) then
      numerator = 0
      denominator = 0
      ie = 0
      ip = 0
      do j = 1, nq
        s = s0*ctx%x(j)
        phi = -depth*(1 - ctx%x(j)**2)
        weight = 2*ctx%x(j)*ctx%w(j)
        ec = 0.5_dp*exp(phi)*(1 + erf(s))
        pc = k*erfc_scaled(s/sqrt(ctx%tau))
        relative_ec = exp_increment(phi) + exp(phi)*(erf(s) - erf(s0))/(1 + erf(s0))
        delta_pc = k*(erfcx_increment(s/sqrt(ctx%tau)) - erfcx_increment(s0/sqrt(ctx%tau)))
        numerator = numerator + weight*(ion_increment(phi, ctx%mach) - relative_ec)
        denominator = denominator + weight*(delta_pc - b*relative_ec)
        ie = ie + depth*weight*ec
        ip = ip + depth*weight*pc
      end do
    else
      ie = 0.5_dp*(1 + erf(s0) - exp(-depth)*(1 + 2*s0/sqrt(pi)))
      ip = k*ctx%tau*(erfc_scaled(s0/sqrt(ctx%tau)) - 1 + 2*s0/sqrt(pi*ctx%tau))
      numerator = (ii - ie/a)/depth
      denominator = (ip - b*ie/a)/depth
    end if
    q = nan()
    amplitude = nan()
    if (denominator == 0) return
    q = numerator/denominator
    amplitude = (1 - b*q)/a
  end subroutine

  function endpoint_field_squared(ctx, r) result(e2)
    type(context_type), intent(in) :: ctx
    type(source_kinetic_root), intent(in) :: r
    real(dp) :: e2, gap, ie, ip, ac, bc, ii, aq, aa, k, s
    k = sqrt(pi/(2*ctx%tau))
    gap = r%phi_h - r%phi_min
    select case (r%branch)
    case ('B', 'A', 'N')
      if (max(abs(r%phi_h), abs(r%phi_min)) < 0.01_dp*min(1._dp, ctx%tau, ctx%mach**2) &
          .or. r%branch == 'N') then
        e2 = field_squared(ctx, r, r%phi_h, .true.)
        return
      end if
      s = sqrt(gap)
      ie = 0.5_dp*r%amplitude*exp(r%phi_min)*(erfc_scaled(s) - 1 + 2*s/sqrt(pi))
      ip = k*ctx%tau*(2*ctx%emission*(-exp_increment(-gap/ctx%tau)) &
                      - r%escaping_flux*(erfc_scaled(s/sqrt(ctx%tau)) - 1 + 2*s/sqrt(pi*ctx%tau)))
      e2 = 2*(ie + ip - ion_integral(r%phi_min, r%phi_h, ctx%mach))
    case ('C')
      if (-r%phi_min < 1.e-3_dp*min(1._dp, ctx%tau)) then
        e2 = field_squared(ctx, r, r%phi_h, .false.)
      else
        call outer_basis(ctx, -r%phi_min, ac, bc, ii, ie, ip, aq, aa)
        e2 = 2*(ii - r%amplitude*ie - r%escaping_flux*ip)
      end if
    case default
      e2 = 0
    end select
  end function

  function source_kinetic_field_squared(mach, tau, emission, root, phi, lower) result(e2)
    real(dp), intent(in) :: mach, tau, emission, phi
    type(source_kinetic_root), intent(in) :: root
    logical, intent(in) :: lower
    real(dp) :: e2
    type(context_type) :: ctx
    call make_context(mach, tau, emission, 0._dp, 1._dp, ctx)
    e2 = field_squared(ctx, root, phi, lower)
  end function

  function field_squared(ctx, r, phi, lower) result(e2)
    type(context_type), intent(in) :: ctx
    type(source_kinetic_root), intent(in) :: r
    real(dp), intent(in) :: phi
    logical, intent(in) :: lower
    real(dp) :: e2, gap, sample, weight, low, top, width, t, ni, ne, npe
    integer :: j
    logical :: from_min, subtract_ion
    from_min = (r%branch == 'A' .or. r%branch == 'N') .and. (lower .or. phi < r%phi_min/2)
    e2 = 0
    if (from_min .or. r%branch == 'B' .or. r%branch == 'B0') then
      gap = phi - r%phi_min
      subtract_ion = phi > ctx%mach**2/4
      do j = 1, nq
        sample = r%phi_min + gap*ctx%x(j)**2
        weight = 2*gap*ctx%x(j)*ctx%w(j)
        if (subtract_ion) then
          call source_kinetic_density(ctx%mach, ctx%tau, ctx%emission, r, sample, lower, ni, ne, npe)
          e2 = e2 + 2*weight*(ne + npe)
        else
          e2 = e2 - 2*weight*charge(ctx, r, sample, lower)
        end if
      end do
      if (subtract_ion) e2 = e2 - 2*ion_integral(r%phi_min, phi, ctx%mach)
    else
      top = sqrt(-r%phi_min)
      low = sqrt(max(0._dp, phi - r%phi_min))
      if (top + low == 0) return
      width = -phi/(top + low)
      do j = 1, nq
        t = low + width*ctx%x(j)
        sample = -width*(1 - ctx%x(j))*(low + top + width*ctx%x(j))
        e2 = e2 + 4*width*t*ctx%w(j)*charge(ctx, r, sample, .false.)
      end do
    end if
  end function

  subroutine validate_root(ctx, r)
    type(context_type), intent(in) :: ctx
    type(source_kinetic_root), intent(inout) :: r
    real(dp) :: ac, bc, ii, ie, ip, aq, aa, gap, lo, hi, value, scale, e2, curvature
    real(dp), allocatable :: grid(:), zeros(:)
    integer :: side, i, sides
    logical :: lower
    r%rejection = 'nonpositive amplitude or ion turning'
    if (r%amplitude <= 0 .or. r%ion_turning_margin <= 0) return
    r%neutrality_residual = charge(ctx, r, 0._dp, .false.)
    r%field_squared_residual = endpoint_field_squared(ctx, r) - ctx%field**2
    gap = r%phi_h - r%phi_min
    if (gap > 0) r%barrier_roundoff = epsilon(gap)*(abs(r%phi_h) + abs(r%phi_min))/gap
    r%rejection = 'emission barrier is unresolved by potential subtraction'
    if (abs(r%escaping_flux - ctx%emission*exp(-gap/ctx%tau)) > 2.e-7_dp*max(1._dp, ctx%emission)) return
    if (r%branch == 'B') then
      r%edge_coefficient = (r%amplitude - 2*sqrt(pi/(2*ctx%tau))*r%escaping_flux/sqrt(ctx%tau))/sqrt(pi)
    else
      value = -r%phi_min
      call outer_basis(ctx, value, ac, bc, ii, ie, ip, aq, aa)
      r%edge_coefficient = 1/ctx%mach**2 - r%amplitude*(ac + exp(-value)/(2*sqrt(pi*value))) &
                           - r%escaping_flux*bc/ctx%tau + r%escaping_flux/(sqrt(2._dp)*ctx%tau*sqrt(value))
    end if
    r%rejection = 'inaccessible neutral edge'
    if (r%edge_coefficient > 2.e-8_dp) return
    sides = 1
    if (r%branch == 'A' .or. r%branch == 'N') then
      sides = 2
      r%minimum_curvature = -charge(ctx, r, r%phi_min, .true.)
      r%outer_residual = -2*(ii - r%amplitude*ie - r%escaping_flux*ip)
      r%rejection = 'nonpositive minimum curvature'
      if (r%minimum_curvature <= 0) return
    end if
    r%rejection = 'non-finite candidate arithmetic'
    if (.not. all(ieee_is_finite([r%neutrality_residual, r%field_squared_residual, r%outer_residual, &
                                  r%edge_coefficient, r%minimum_curvature, r%escaping_flux, r%current]))) return
    r%rejection = 'closure residual exceeds tolerance'
    if (abs(r%neutrality_residual) > 2.e-7_dp .or. abs(r%outer_residual) > 2.e-7_dp .or. &
        abs(r%field_squared_residual) > 2.e-7_dp*max(1._dp, ctx%field**2)) return
    r%minimum_field_squared = 0
    do side = 1, sides
      lower = sides == 2 .and. side == 1
      lo = r%phi_min
      hi = 0
      if (lower .or. r%branch == 'B') hi = r%phi_h
      call bounded_grid(lo, hi, 129, grid)
      grid = [lo, grid, hi]
      call find_roots(rho_function, grid, 1.e-13_dp, zeros)
      do i = 1, size(zeros)
        if (zeros(i) <= lo + 1.e-8_dp*(hi - lo) .or. zeros(i) >= hi - 1.e-8_dp*(hi - lo)) cycle
        e2 = field_squared(ctx, r, zeros(i), lower)
        r%rejection = 'interior field zero blocks the connecting profile'
        if (e2 <= 1.e-11_dp*max(ctx%field**2, abs(endpoint_field_squared(ctx, r)), tiny(e2))) return
      end do
      grid = [grid, zeros]
      do i = 1, size(grid)
        value = grid(i)
        e2 = field_squared(ctx, r, value, lower)
        r%minimum_field_squared = min(r%minimum_field_squared, e2)
        scale = max(min(abs(value - r%phi_min), abs(value)), 1.e-14_dp)**2
        if (lower) scale = max(value - r%phi_min, 1.e-14_dp)
        r%rejection = 'non-finite candidate arithmetic'
        if (.not. ieee_is_finite(e2)) return
        r%rejection = 'negative field squared inside the profile'
        if (value > lo .and. value < hi .and. e2 < -2.e-9_dp*max(scale, 1.e-10_dp)) return
      end do
    end do
    if (r%boundary_turning) then
      curvature = charge(ctx, r, r%phi_h, .false.)
      r%rejection = 'wrong curvature at the zero-field boundary'
      if (r%branch == 'B' .and. curvature <= 0) return
      if (r%branch == 'C' .and. curvature >= 0) return
    end if
    r%accepted = .true.
    r%rejection = ''
  contains
    function rho_function(phi) result(y)
      real(dp), intent(in) :: phi
      real(dp) :: y
      y = charge(ctx, r, phi, lower)
    end function
  end subroutine

  subroutine depth_grid(lo, hi, points, grid)
    real(dp), intent(in) :: lo, hi
    integer, intent(in) :: points
    real(dp), allocatable, intent(out) :: grid(:)
    integer :: i, n, part
    real(dp) :: a, b
    allocate (grid(0))
    do part = 1, 3
      select case (part)
      case (1)
        a = lo
        b = min(1.e-12_dp, hi)
        n = max(65, points/3)
      case (2)
        a = 1.e-12_dp
        b = min(1.e4_dp, hi)
        n = points
      case (3)
        a = 1.e4_dp
        b = hi
        n = max(3, ceiling(log10(max(1._dp, hi/1.e4_dp))*points/16) + 1)
      end select
      if (b <= a) cycle
      do i = 0, n - 1
        grid = [grid, exp(log(a) + real(i, dp)/(n - 1)*log(b/a))]
      end do
    end do
    call sort_unique(grid)
  end subroutine

  subroutine bounded_grid(lo, hi, points, grid)
    real(dp), intent(in) :: lo, hi
    integer, intent(in) :: points
    real(dp), allocatable, intent(out) :: grid(:)
    integer :: i, n
    real(dp) :: f
    allocate (grid(0))
    n = max(24, points/2)
    do i = 1, n - 1
      f = real(i, dp)/n
      grid = [grid, lo + (hi - lo)*f]
    end do
    do i = 0, n - 1
      f = exp(log(1.e-12_dp) + real(i, dp)/(n - 1)*log(0.5e12_dp))
      grid = [grid, lo + (hi - lo)*f, hi - (hi - lo)*f]
    end do
    call sort_unique(grid)
  end subroutine

  subroutine add_boundary_neighbors(grid, boundaries)
    real(dp), allocatable, intent(inout) :: grid(:)
    real(dp), intent(in) :: boundaries(:)
    integer :: i, j
    do i = 1, size(boundaries)
      do j = 4, 8, 2
        grid = [grid, boundaries(i)*(1 - 10._dp**(-j)), boundaries(i)*(1 + 10._dp**(-j))]
      end do
    end do
    call sort_unique(grid)
  end subroutine

  subroutine find_roots(fun, samples, tolerance, roots)
    procedure(scalar_function) :: fun
    real(dp), intent(in) :: samples(:), tolerance
    real(dp), allocatable, intent(out) :: roots(:)
    real(dp), allocatable :: grid(:), values(:), extrema(:)
    real(dp) :: left, right, fl, fr, mid, fm, x1, x2, f1, f2, direction, ratio
    integer :: i, j
    grid = samples
    call sort_unique(grid)
    allocate (values(size(grid)), roots(0), extrema(0))
    do i = 1, size(grid)
      values(i) = fun(grid(i))
    end do
    ! Resolve sampled extrema to detect tangent roots and two crossings hidden
    ! in one interval. Finite sampling still gives no completeness guarantee.
    ratio = (sqrt(5._dp) - 1)/2
    do i = 2, size(grid) - 1
      if (.not. all(ieee_is_finite(values(i - 1:i + 1)))) cycle
      if ((values(i) - values(i - 1))*(values(i + 1) - values(i)) >= 0) cycle
      direction = 1
      if (values(i) > values(i - 1)) direction = -1
      left = grid(i - 1)
      right = grid(i + 1)
      do j = 1, 65
        x1 = right - ratio*(right - left)
        x2 = left + ratio*(right - left)
        f1 = direction*fun(x1)
        f2 = direction*fun(x2)
        if (.not. all(ieee_is_finite([f1, f2]))) exit
        if (f1 < f2) then
          right = x2
        else
          left = x1
        end if
      end do
      mid = (left + right)/2
      extrema = [extrema, mid]
      fm = fun(mid)
      if (ieee_is_finite(fm) .and. abs(fm) <= tolerance) roots = [roots, mid]
    end do
    grid = [grid, extrema]
    call sort_unique(grid)
    left = grid(1)
    fl = fun(left)
    do i = 2, size(grid)
      right = grid(i)
      fr = fun(right)
      if (fr == 0) roots = [roots, right]
      if (all(ieee_is_finite([fl, fr]))) then
        if (fl /= 0 .and. fr /= 0 .and. (fl > 0 .neqv. fr > 0)) then
          do j = 1, 100
            mid = (left + right)/2
            fm = fun(mid)
            if (.not. ieee_is_finite(fm)) exit
            if (fm == 0 .or. right - left <= 3.e-13_dp*max(abs(mid), 1.e-300_dp)) exit
            if (fm > 0 .eqv. fl > 0) then
              left = mid
              fl = fm
            else
              right = mid
            end if
          end do
          if (ieee_is_finite(fm)) roots = [roots, mid]
        end if
      end if
      left = grid(i)
      fl = fr
    end do
    call sort_unique(roots)
  end subroutine

  subroutine sort_unique(values)
    real(dp), allocatable, intent(inout) :: values(:)
    real(dp) :: value
    integer :: i, j, count
    do i = 2, size(values)
      value = values(i)
      j = i - 1
      do while (j >= 1)
        if (values(j) <= value) exit
        values(j + 1) = values(j)
        j = j - 1
      end do
      values(j + 1) = value
    end do
    count = 0
    do i = 1, size(values)
      if (count > 0) then
        if (values(i) == values(count)) cycle
      end if
      count = count + 1
      values(count) = values(i)
    end do
    values = values(:count)
  end subroutine

  function nan() result(value)
    real(dp) :: value
    value = ieee_value(0._dp, ieee_quiet_nan)
  end function
end module bem_source_kinetic_sheath
