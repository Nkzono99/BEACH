!> Zhao の初期推定、減衰 Newton 法、差分 Jacobian、小規模線形解法。
!! 残差の定義と物理的な許容条件は physics、最終的な根の選択は roots が担当する。
submodule(bem_matching_plane_zhao) bem_matching_plane_zhao_numerics
  implicit none

  integer, parameter :: root_max_iterations = 60
  integer, parameter :: root_max_backtracks = 24
  real(dp), parameter :: root_tolerance = 1.0e-9_dp

contains

  module procedure make_matching_branch_guesses

  real(dp) :: density, phi0, phi_m, ion, electron_free, electron_reflected, photo, captured
  real(dp) :: photo_scale, field_scale, ion_limit, positive_scale, negative_scale
  real(dp) :: potential(matching_root_seed_count), minimum(matching_root_seed_count)
  real(dp) :: source_ratio, total_flux, cumulative, median_energy
  integer :: index, bin
  logical :: valid
  character(len=9) :: side

  ! Every voltage scale follows the current query. For spectra, use the flux
  ! median as well as the mean so a small energetic tail cannot set every seed.
  photo_scale = params%t_phe_ev
  if (allocated(params%pe_spectrum%flux)) then
    total_flux = params%pe_spectrum%total_flux()
    cumulative = 0.0_dp
    do bin = 1, size(params%pe_spectrum%flux)
      cumulative = cumulative + params%pe_spectrum%flux(bin)
      if (total_flux <= 0.0_dp .or. cumulative < 0.5_dp*total_flux) cycle
      median_energy = 0.5_dp*(params%pe_spectrum%edge(bin - 1) + params%pe_spectrum%edge(bin))
      photo_scale = max(median_energy, 0.25_dp*params%t_phe_ev)
      exit
    end do
  end if
  field_scale = params%t_phe_ev*target_field_hat**2
  ion_limit = 0.5_dp*params%t_swe_ev*params%mach**2
  source_ratio = max(1.0_dp, params%n_phe0_m3/(2.0_dp*params%n_swi_inf_m3))
  positive_scale = max(photo_scale*(1.0_dp + log(source_ratio)), field_scale)
  negative_scale = max(params%t_swe_ev, photo_scale, field_scale)
  guesses = 0.0_dp
  count = 0
  select case (branch)
  case ('A')
    potential = positive_scale*[0.25_dp, 0.75_dp, 1.5_dp, 3.0_dp, 0.5_dp, 2.0_dp, -0.1_dp, -0.5_dp]
    minimum = -params%t_swe_ev*[0.02_dp, 0.1_dp, 0.3_dp, 1.0_dp, 0.5_dp, 0.02_dp, 0.3_dp, 1.0_dp]
    potential(7:8) = -negative_scale*[0.25_dp, 0.5_dp]
    minimum(7) = potential(7) - max(field_scale, 0.01_dp*photo_scale)
    minimum(8) = potential(8) - photo_scale
    side = 'upper'
  case ('B')
    potential = positive_scale*[0.002_dp, 0.02_dp, 0.1_dp, 0.3_dp, 0.7_dp, 1.5_dp, 3.0_dp, 6.0_dp]
    minimum = potential
    side = 'monotonic'
  case ('C')
    potential = -negative_scale*[0.002_dp, 0.02_dp, 0.1_dp, 0.3_dp, 0.7_dp, 1.5_dp, 3.0_dp, 6.0_dp]
    minimum = potential
    side = 'monotonic'
  case default
    return
  end select
  do index = 1, matching_root_seed_count
    phi0 = min(potential(index), 0.9_dp*ion_limit)
    phi_m = minimum(index)
    if (branch == 'A') phi_m = min(phi_m, phi0 - 0.1_dp*photo_scale)
    if (branch /= 'A') phi_m = phi0
    call evaluate_zhao_density_hat(params, branch, trim(side), 0.0_dp, phi0/params%t_phe_ev, &
                                   phi_m/params%t_phe_ev, 1.0_dp, ion, electron_free, &
                                   electron_reflected, photo, captured)
    if (.not. all(ieee_is_finite([ion, electron_free, electron_reflected, photo, captured]))) cycle
    if (electron_free + electron_reflected <= 0.0_dp .or. ion <= photo + captured) cycle
    density = params%n_phe_ref_m3*(ion - photo - captured)/(electron_free + electron_reflected)
    call encode_matching_unknowns(params, branch, phi0, phi_m, density, guesses(:, count + 1), valid)
    if (valid) count = count + 1
  end do
  end procedure make_matching_branch_guesses

  module procedure newton_matching_branch

  real(dp) :: y(3), f(3), jac(3, 3), delta(3), trial(3), trial_f(3)
  real(dp) :: norm, trial_norm, step
  integer :: n, iteration, backtrack
  logical :: valid, jacobian_ok, linear_ok, trial_valid

  n = merge(3, 2, branch == 'A')
  y = y0
  call evaluate_charge_residual(params, branch, target_field_hat, y, f, valid)
  if (.not. valid) then
    y_out = y
    final_norm = huge(1.0_dp)
    iterations = 0
    success = .false.
    return
  end if
  norm = maxval(abs(f(1:n)))
  success = .false.
  do iteration = 0, root_max_iterations
    if (norm <= root_tolerance) then
      success = .true.
      exit
    end if
    if (iteration == root_max_iterations) exit
    call matching_numerical_jacobian( &
      params, branch, target_field_hat, y, f, n, jac, jacobian_ok &
      )
    if (.not. jacobian_ok) exit
    call solve_matching_small_system(jac, -f, n, delta, linear_ok)
    if (.not. linear_ok) exit
    step = 1.0_dp
    do backtrack = 1, root_max_backtracks
      trial = y + step*delta
      call evaluate_charge_residual( &
        params, branch, target_field_hat, trial, trial_f, trial_valid &
        )
      if (trial_valid) then
        trial_norm = maxval(abs(trial_f(1:n)))
        if (trial_norm < norm) then
          y = trial
          f = trial_f
          norm = trial_norm
          exit
        end if
      end if
      step = 0.5_dp*step
    end do
    if (backtrack > root_max_backtracks) exit
  end do
  y_out = y
  final_norm = norm
  iterations = iteration
  end procedure newton_matching_branch

  subroutine matching_numerical_jacobian( &
    params, branch, target_field_hat, y, f0, n, jac, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y(3), f0(3)
    integer, intent(in) :: n
    real(dp), intent(out) :: jac(3, 3)
    logical, intent(out) :: success

    real(dp) :: yp(3), ym(3), fp(3), fm(3), h
    integer :: column
    logical :: plus_valid, minus_valid

    jac = 0.0_dp
    success = .true.
    do column = 1, n
      h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(y(column)))
      yp = y
      ym = y
      yp(column) = yp(column) + h
      ym(column) = ym(column) - h
      call evaluate_charge_residual(params, branch, target_field_hat, yp, fp, plus_valid)
      call evaluate_charge_residual(params, branch, target_field_hat, ym, fm, minus_valid)
      if (plus_valid .and. minus_valid) then
        jac(1:n, column) = (fp(1:n) - fm(1:n))/(2.0_dp*h)
      else if (plus_valid) then
        jac(1:n, column) = (fp(1:n) - f0(1:n))/h
      else if (minus_valid) then
        jac(1:n, column) = (f0(1:n) - fm(1:n))/h
      else
        success = .false.
        return
      end if
    end do
    success = all(ieee_is_finite(jac(1:n, 1:n)))
  end subroutine matching_numerical_jacobian

  subroutine solve_matching_small_system(a_in, b_in, n, x, success)
    real(dp), intent(in) :: a_in(3, 3), b_in(3)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(3)
    logical, intent(out) :: success

    real(dp) :: a(3, 3), b(3), factor, pivot_value, tmp
    integer :: i, j, k, pivot

    a = a_in
    b = b_in
    x = 0.0_dp
    success = .false.
    do k = 1, n
      pivot = k
      do i = k + 1, n
        if (abs(a(i, k)) > abs(a(pivot, k))) pivot = i
      end do
      if (.not. ieee_is_finite(a(pivot, k)) .or. abs(a(pivot, k)) <= 1.0e-14_dp) return
      if (pivot /= k) then
        do j = k, n
          tmp = a(k, j)
          a(k, j) = a(pivot, j)
          a(pivot, j) = tmp
        end do
        tmp = b(k)
        b(k) = b(pivot)
        b(pivot) = tmp
      end if
      pivot_value = a(k, k)
      do i = k + 1, n
        factor = a(i, k)/pivot_value
        a(i, k:n) = a(i, k:n) - factor*a(k, k:n)
        b(i) = b(i) - factor*b(k)
      end do
    end do
    do i = n, 1, -1
      x(i) = b(i)
      do j = i + 1, n
        x(i) = x(i) - a(i, j)*x(j)
      end do
      x(i) = x(i)/a(i, i)
    end do
    success = all(ieee_is_finite(x(1:n)))
  end subroutine solve_matching_small_system

end submodule bem_matching_plane_zhao_numerics
