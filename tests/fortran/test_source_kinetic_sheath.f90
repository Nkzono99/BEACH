program test_source_kinetic_sheath
  use bem_kinds, only: dp
  use bem_constants, only: pi
  use bem_source_kinetic_sheath
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none
  type(source_kinetic_result) :: answer
  type(source_kinetic_root) :: r
  real(dp) :: ni, ne, npe, phi, integral, t, weight, v, h
  integer :: i, j, ib, ineg

  call test_init(7)
  call test_begin('B density equals transported incoming velocity integral')
  r%branch = 'B'
  r%phi_h = 100._dp
  r%amplitude = 0.7_dp
  do i = 1, 4
    phi = 10._dp**(i - 2)
    call source_kinetic_density(100._dp, 0.2_dp, 0._dp, r, phi, .true., ni, ne, npe)
    integral = 0
    h = 12._dp/8000
    do j = 0, 8000
      t = j*h
      v = sqrt(2*phi + t*t)
      weight = 2
      if (mod(j, 2) == 1) weight = 4
      if (j == 0 .or. j == 8000) weight = 1
      integral = integral + weight*exp(-t*t/2)*t/v
    end do
    integral = r%amplitude/sqrt(2*pi)*integral*h/3
    call assert_close_dp(ne, integral, 2.e-11_dp, 'accelerated ambient density disagrees with velocity integral')
  end do
  call test_end()

  call test_begin('issue 38 coexisting B and N roots and nonzero current')
  call solve_source_kinetic(10._dp, 0.2_dp, 0.3_dp, 0.1_dp, 1836.15267343_dp, answer)
  call show_answer(answer)
  call assert_true(answer%status == source_kinetic_ok, 'coexistence search failed')
  call assert_true(size(answer%roots) == 2, 'both B and N must be retained')
  ib = 0
  ineg = 0
  do i = 1, size(answer%roots)
    if (answer%roots(i)%branch == 'B') ib = i
    if (answer%roots(i)%branch == 'N') ineg = i
  end do
  call assert_true(min(ib, ineg) > 0, 'branch classification changed')
  call assert_close_dp(answer%roots(ib)%phi_h, 0.0228867313_dp, 8.e-10_dp, 'B wall fixture')
  call assert_close_dp(answer%roots(ib)%amplitude, 0.5003210848_dp, 8.e-10_dp, 'B amplitude fixture')
  call assert_close_dp(answer%roots(ineg)%phi_h, -0.1329182235_dp, 8.e-10_dp, 'N wall fixture')
  call assert_close_dp(answer%roots(ineg)%phi_min, -0.1480124772_dp, 8.e-10_dp, 'N minimum fixture')
  call assert_close_dp(answer%roots(ineg)%amplitude, 0.8974844341_dp, 8.e-10_dp, 'N amplitude fixture')
  call assert_close_dp(answer%roots(ineg)%escaping_flux, 0.2781919119_dp, 8.e-10_dp, 'N escape fixture')
  call assert_close_dp(answer%roots(ineg)%current, 0.2027773700_dp, 8.e-10_dp, 'current must remain an output')
  call assert_true(.not. answer%search_complete, 'finite root search cannot certify completeness')
  call check_orbits(answer%roots(ineg), 0.3_dp, 0.1_dp)
  call test_end()

  call test_begin('Type A fixture satisfies independent orbit and Poisson integrals')
  call solve_source_kinetic(10._dp, 0.2_dp, 3.3_dp, 2.440435545625_dp, 1836.15267343_dp, answer)
  ib = 0
  do i = 1, size(answer%roots)
    if (answer%roots(i)%branch == 'A') ib = i
  end do
  call assert_true(ib > 0, 'zero-drift Type A fixture missing')
  if (ib > 0) then
    r = answer%roots(ib)
    call assert_close_dp(r%phi_h, 0.538231013617_dp, 2.e-10_dp, 'Type A wall potential')
    call assert_close_dp(r%phi_min, -0.0380442797254_dp, 2.e-11_dp, 'Type A minimum')
    call check_orbits(r, 3.3_dp, 2.440435545625_dp)
  end if
  call test_end()

  call test_begin('analytical exclusion differs from finite depth truncation')
  call solve_source_kinetic(10._dp, 0.2_dp, 0._dp, 0.1_dp, 1836._dp, answer)
  call assert_true(answer%status == source_kinetic_excluded .and. answer%search_complete, 'no-PE exclusion')
  call solve_source_kinetic(10._dp, 0.2_dp, 10._dp, 0.3_dp, 1836._dp, answer, maximum_depth=1.e4_dp)
  call assert_true(answer%status == source_kinetic_unresolved .and. .not. answer%search_complete, 'depth truncation')
  call test_end()

  call test_begin('zero field includes homogeneous and nonuniform endpoint states')
  call solve_source_kinetic(10._dp, 0.2_dp, 0._dp, 0._dp, 1836._dp, answer)
  call assert_true(size(answer%roots) == 1, 'dark homogeneous root count')
  call assert_true(answer%roots(1)%branch == 'B0', 'homogeneous state must be diagnosed as B0')
  call assert_close_dp(answer%roots(1)%amplitude, 2._dp, 0._dp, 'B0 neutrality')
  call solve_source_kinetic(10._dp, 0.2_dp, 1._dp, 0._dp, 1836._dp, answer)
  call show_answer(answer)
  call assert_true(size(answer%roots) > 0, 'nonuniform zero-field endpoint missing')
  call assert_true(answer%roots(1)%branch == 'C' .and. answer%roots(1)%boundary_turning, 'C endpoint classification')
  call assert_true(answer%roots(1)%phi_h < 0, 'zero field does not imply zero potential')
  call test_end()

  call test_begin('deep N roots retain neutral-edge and barrier diagnostics')
  call solve_source_kinetic(10._dp, 0.2_dp, 10._dp, 0.3_dp, 1836._dp, answer)
  call show_answer(answer)
  call assert_true(size(answer%roots) > 0, 'deep N root missing')
  r = answer%roots(1)
  call assert_true(r%branch == 'N' .and. r%phi_min < -1.e4_dp, 'deep branch classification')
  call assert_true(abs(r%outer_residual) < 2.e-7_dp, 'deep outer neutrality/first integral')
  call test_end()

  call test_begin('negative field dark C response is nonzero and source-free')
  call solve_source_kinetic(10._dp, 0.2_dp, 0._dp, -0.1_dp, 1836._dp, answer)
  call assert_true(size(answer%roots) > 0, 'dark negative response missing')
  call assert_true(answer%roots(1)%branch == 'C', 'dark negative branch')
  call assert_close_dp(answer%roots(1)%escaping_flux, 0._dp, 0._dp, 'dark escape')
  call test_end()
  call test_summary()
contains
  subroutine check_orbits(root, emission, field)
    type(source_kinetic_root), intent(in) :: root
    real(dp), intent(in) :: emission, field
    real(dp), parameter :: tau = 0.2_dp
    real(dp) :: potential, gap, upper, u, w, ne_integral, pe_integral, ni, ne, npe, du, integral
    integer :: side, k, j
    logical :: lower
    ! Direct local-velocity integrals of the transported source VDFs. The
    ! reflected population occupies velocities below the minimum barrier.
    do side = 1, 2
      lower = side == 1
      upper = 0
      if (lower) upper = root%phi_h
      do k = 1, 3
        potential = root%phi_min + (upper - root%phi_min)*k/4
        gap = potential - root%phi_min
        ne_integral = 0
        pe_integral = 0
        du = 12._dp/8000
        do j = 0, 8000
          u = j*du
          w = simpson_weight(j, 8000)*du/3
          if (lower) then
            ne_integral = ne_integral + w*root%amplitude*exp(root%phi_min - u*u/2)/sqrt(2*pi) &
                          *u/sqrt(u*u + 2*gap)
            pe_integral = pe_integral + w*emission/tau*exp((potential - root%phi_h)/tau - u*u/(2*tau))
          else
            ne_integral = ne_integral + w*root%amplitude*exp(potential - u*u/2)/sqrt(2*pi)
            pe_integral = pe_integral + w*emission/tau*exp((root%phi_min - root%phi_h)/tau - u*u/(2*tau)) &
                          *u/sqrt(u*u + 2*gap)
          end if
        end do
        du = sqrt(2*gap)/2000
        do j = 0, 2000
          u = j*du
          w = simpson_weight(j, 2000)*du/3
          if (lower) then
            pe_integral = pe_integral + w*emission/tau*exp((potential - root%phi_h)/tau - u*u/(2*tau))
          else
            ne_integral = ne_integral + w*root%amplitude*exp(potential - u*u/2)/sqrt(2*pi)
          end if
        end do
        call source_kinetic_density(10._dp, tau, emission, root, potential, lower, ni, ne, npe)
        call assert_close_dp(ne, ne_integral, 3.e-8_dp, 'A/N ambient source orbit integral')
        call assert_close_dp(npe, pe_integral, 3.e-8_dp, 'A/N emitted source orbit integral')
      end do
      integral = 0
      do j = 0, 4000
        u = real(j, dp)/4000
        potential = root%phi_min + (upper - root%phi_min)*u*u
        call source_kinetic_density(10._dp, tau, emission, root, potential, lower, ni, ne, npe)
        integral = integral - simpson_weight(j, 4000)*4*(upper - root%phi_min)*u*(ni - ne - npe)/12000
      end do
      if (lower) then
        call assert_close_dp(integral, field**2, 2.e-9_dp, 'inner Poisson first integral')
      else
        call assert_close_dp(integral, 0._dp, 2.e-9_dp, 'outer Poisson first integral')
      end if
    end do
  end subroutine

  real(dp) function simpson_weight(i, n) result(w)
    integer, intent(in) :: i, n
    w = 2
    if (mod(i, 2) == 1) w = 4
    if (i == 0 .or. i == n) w = 1
  end function

  subroutine show_answer(result)
    type(source_kinetic_result), intent(in) :: result
    integer :: k
    print *, trim(result%message)
    do k = 1, size(result%roots)
      print *, result%roots(k)%branch, result%roots(k)%phi_h, result%roots(k)%phi_min, result%roots(k)%amplitude
    end do
    do k = 1, size(result%rejected)
      print *, 'rejected ', result%rejected(k)%branch, result%rejected(k)%phi_h, trim(result%rejected(k)%rejection)
    end do
  end subroutine
end program
