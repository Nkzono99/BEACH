!> Charge-driven quasi-steady Zhao outer-plasma closure.
module bem_outer_plasma_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi, qe
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid, &
                                    outer_plasma_no_physical_solution, outer_plasma_numerical_failure
  use bem_sheath_model_core, only: zhao_params_type, evaluate_zhao_rho_hat, zhao_residuals_type_a, &
                                   zhao_residuals_type_b, zhao_residuals_type_c, swe_free_current_term
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer, parameter :: root_max_iterations = 60
  integer, parameter :: root_max_backtracks = 24
  integer, parameter :: rho_quadrature_panels = 256
  real(dp), parameter :: root_tolerance = 1.0e-9_dp
  real(dp), parameter :: profile_phi_end_hat = 1.0e-4_dp
  real(dp), parameter :: zero_field_tolerance_hat = 1.0e-12_dp

  type, public :: zhao_charge_root_type
    character(len=1) :: branch = ' '
    character(len=16) :: interface_side = 'lower'
    real(dp) :: phi0_v = 0.0_dp
    real(dp) :: phi_m_v = 0.0_dp
    real(dp) :: n_swe_inf_m3 = 0.0_dp
    real(dp) :: interface_field_v_m = 0.0_dp
    real(dp) :: net_current_density_a_m2 = 0.0_dp
    real(dp) :: residual_norm = 0.0_dp
    integer(i32) :: nonlinear_iterations = 0_i32
  end type zhao_charge_root_type

  public :: solve_zhao_charge_root
  public :: evaluate_zhao_interface_field
  public :: build_zhao_outer_profile
  public :: solve_outer_plasma_zhao
  public :: zhao_net_current_density

contains

  subroutine evaluate_zhao_interface_field(params, root, interface_field_v_m, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    real(dp), intent(out) :: interface_field_v_m
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: phi0_hat, phi_m_hat, density_hat, integral, field_squared, sign_field
    logical :: integral_ok

    interface_field_v_m = 0.0_dp
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. root%n_swe_inf_m3 <= 0.0_dp) then
      message = 'invalid Zhao root for interface-field evaluation'
      return
    end if
    if (root%branch == '0') then
      status = outer_plasma_ok
      return
    end if
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    density_hat = root%n_swe_inf_m3/params%n_phe_ref_m3
    select case (root%branch)
    case ('A')
      call integrate_rho_hat( &
        params, 'A', 'lower', phi_m_hat, phi0_hat, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      field_squared = -2.0_dp*integral
      sign_field = 1.0_dp
    case ('B', 'C')
      call integrate_rho_hat( &
        params, root%branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      field_squared = 2.0_dp*integral
      if (root%branch == 'B') then
        sign_field = 1.0_dp
      else
        sign_field = -1.0_dp
      end if
    case default
      message = 'unknown Zhao branch for interface-field evaluation'
      return
    end select
    if (.not. integral_ok .or. field_squared < -1.0e-10_dp) then
      status = outer_plasma_no_physical_solution
      message = 'Zhao root has no real interface field'
      return
    end if
    interface_field_v_m = sign_field*sqrt(max(0.0_dp, field_squared))* &
                          params%t_phe_ev/params%lambda_d_phe_ref_m
    if (.not. ieee_is_finite(interface_field_v_m)) then
      status = outer_plasma_numerical_failure
      message = 'non-finite Zhao interface field'
      return
    end if
    status = outer_plasma_ok
  end subroutine evaluate_zhao_interface_field

  subroutine solve_outer_plasma_zhao( &
    model, params, interface_field_v_m, grid_points, state, root, status, message, initial_root &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    integer(i32), intent(in) :: grid_points
    type(outer_plasma_state_type), intent(out) :: state
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root

    if (present(initial_root)) then
      call solve_zhao_charge_root( &
        model, params, interface_field_v_m, root, status, message, initial_root=initial_root &
        )
    else
      call solve_zhao_charge_root(model, params, interface_field_v_m, root, status, message)
    end if
    if (status /= outer_plasma_ok) then
      state = outer_plasma_state_type()
      state%model = 'kinetic_1d'
      state%applicability_status = status
      return
    end if
    call build_zhao_outer_profile(params, root, grid_points, state, status, message)
  end subroutine solve_outer_plasma_zhao

  subroutine solve_zhao_charge_root( &
    model, params, interface_field_v_m, root, status, message, initial_root &
    )
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root

    character(len=16) :: normalized_model
    character(len=1) :: order(3), candidate
    type(zhao_charge_root_type) :: trial_root
    real(dp) :: field_scale, target_field_hat, degenerate_density_hat
    integer :: candidate_count, candidate_index, warm_index
    character(len=1) :: first_candidate

    root = zhao_charge_root_type()
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. .not. ieee_is_finite(interface_field_v_m)) then
      message = 'invalid Zhao charge-driven parameters'
      return
    end if
    field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
    target_field_hat = interface_field_v_m/field_scale
    normalized_model = trim(lower_ascii(model))

    if (abs(target_field_hat) <= zero_field_tolerance_hat) then
      degenerate_density_hat = ( &
                               2.0_dp*params%n_swi_inf_m3/params%n_phe_ref_m3 - &
                               params%n_phe0_m3/params%n_phe_ref_m3 &
                               )/(1.0_dp + erf(params%u))
      root%phi0_v = 0.0_dp
      root%phi_m_v = 0.0_dp
      root%interface_field_v_m = 0.0_dp
      root%residual_norm = 0.0_dp
      root%nonlinear_iterations = 0_i32
      if (degenerate_density_hat > 0.0_dp .and. ieee_is_finite(degenerate_density_hat)) then
        root%branch = 'B'
        root%interface_side = 'monotonic'
        root%n_swe_inf_m3 = degenerate_density_hat*params%n_phe_ref_m3
        message = 'zero-field degenerate Zhao-B state'
      else
        root%branch = '0'
        root%interface_side = 'bootstrap'
        root%n_swe_inf_m3 = params%n_swi_inf_m3
        message = 'zero-field transient bootstrap; photoelectron outer population is not yet formed'
      end if
      root%net_current_density_a_m2 = zhao_net_current_density(params, root)
      status = outer_plasma_ok
      return
    end if

    call branch_order(normalized_model, params, target_field_hat, order, candidate_count, status, message)
    if (status /= outer_plasma_ok) return
    if (present(initial_root)) then
      if ((normalized_model == 'auto' .or. normalized_model == 'zhao_auto') .and. &
          index('ABC', initial_root%branch) > 0) then
        do warm_index = 1, candidate_count
          if (order(warm_index) /= initial_root%branch) cycle
          first_candidate = order(1)
          order(1) = order(warm_index)
          order(warm_index) = first_candidate
          exit
        end do
      end if
    end if
    do candidate_index = 1, candidate_count
      candidate = order(candidate_index)
      if (present(initial_root)) then
        call solve_one_branch( &
          params, candidate, interface_field_v_m, trial_root, status, message, initial_root=initial_root &
          )
      else
        call solve_one_branch(params, candidate, interface_field_v_m, trial_root, status, message)
      end if
      if (status == outer_plasma_ok) then
        root = trial_root
        return
      end if
    end do
    status = outer_plasma_no_physical_solution
    message = 'no charge-driven Zhao branch satisfies the prescribed interface field'
  end subroutine solve_zhao_charge_root

  subroutine branch_order(model, params, target_field_hat, order, count, status, message)
    character(len=*), intent(in) :: model
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: target_field_hat
    character(len=1), intent(out) :: order(3)
    integer, intent(out) :: count
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    order = ' '
    count = 0
    status = outer_plasma_ok
    message = ''
    select case (trim(model))
    case ('zhao_a', 'a')
      order(1) = 'A'
      count = 1
    case ('zhao_b', 'b')
      order(1) = 'B'
      count = 1
    case ('zhao_c', 'c')
      order(1) = 'C'
      count = 1
    case ('zhao_auto', 'auto')
      if (target_field_hat > 0.0_dp) then
        order = ['A', 'B', 'C']
      else
        order = ['C', 'A', 'B']
      end if
      if (params%alpha_rad*180.0_dp/pi < 20.0_dp .and. target_field_hat > 0.0_dp) then
        order = ['B', 'A', 'C']
      end if
      count = 3
    case default
      status = outer_plasma_invalid
      message = 'unknown Zhao charge-driven model'
    end select
  end subroutine branch_order

  subroutine solve_one_branch( &
    params, branch, interface_field_v_m, root, status, message, initial_root &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: interface_field_v_m
    type(zhao_charge_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(zhao_charge_root_type), intent(in), optional :: initial_root

    real(dp) :: target_field_hat, field_scale
    real(dp) :: guesses(3, 8), y(3), best_y(3), norm, best_norm
    integer :: guess_count, guess_index, iterations, best_iterations
    logical :: guess_valid, success

    root = zhao_charge_root_type()
    root%branch = branch
    root%interface_field_v_m = interface_field_v_m
    status = outer_plasma_no_physical_solution
    message = ''
    field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
    target_field_hat = interface_field_v_m/field_scale
    if ((branch == 'A' .or. branch == 'B') .and. target_field_hat <= 0.0_dp) then
      message = 'Zhao A/B lower interface requires a positive normal field'
      return
    end if
    if (branch == 'C' .and. target_field_hat >= 0.0_dp) then
      message = 'Zhao C interface requires a negative normal field'
      return
    end if

    call make_branch_guesses(params, branch, guesses, guess_count)
    if (present(initial_root)) then
      if (initial_root%branch == branch) then
        call encode_unknowns( &
          params, branch, initial_root%phi0_v, initial_root%phi_m_v, &
          initial_root%n_swe_inf_m3, y, guess_valid &
          )
        if (guess_valid) then
          guesses(:, 2:guess_count + 1) = guesses(:, 1:guess_count)
          guesses(:, 1) = y
          guess_count = guess_count + 1
        end if
      end if
    end if

    best_norm = huge(1.0_dp)
    best_y = guesses(:, 1)
    best_iterations = 0
    do guess_index = 1, guess_count
      call newton_branch( &
        params, branch, target_field_hat, guesses(:, guess_index), y, norm, iterations, success &
        )
      if (norm < best_norm) then
        best_norm = norm
        best_y = y
        best_iterations = iterations
      end if
      if (success) exit
    end do
    if (.not. success) then
      status = outer_plasma_numerical_failure
      message = 'charge-driven Zhao root Newton solve did not converge'
      return
    end if

    call decode_unknowns( &
      params, branch, y, root%phi0_v, root%phi_m_v, root%n_swe_inf_m3, guess_valid &
      )
    if (.not. guess_valid) then
      status = outer_plasma_numerical_failure
      message = 'charge-driven Zhao root decoded to an invalid state'
      return
    end if
    if (branch == 'A') then
      root%interface_side = 'lower'
    else
      root%interface_side = 'monotonic'
    end if
    root%residual_norm = norm
    root%nonlinear_iterations = int(iterations, i32)
    root%net_current_density_a_m2 = zhao_net_current_density(params, root)
    status = outer_plasma_ok
    message = ''
  end subroutine solve_one_branch

  subroutine make_branch_guesses(params, branch, guesses, count)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(out) :: guesses(3, 8)
    integer, intent(out) :: count

    real(dp) :: phi0, phi_m, density
    logical :: valid

    guesses = 0.0_dp
    select case (branch)
    case ('A')
      count = 5
      call encode_unknowns(params, branch, 3.6_dp, -0.5_dp, 8.2e6_dp, guesses(:, 1), valid)
      call encode_unknowns(params, branch, 2.8_dp, -0.3_dp, 8.0e6_dp, guesses(:, 2), valid)
      call encode_unknowns(params, branch, 4.5_dp, -0.8_dp, 8.4e6_dp, guesses(:, 3), valid)
      phi0 = 1.3_dp*params%t_phe_ev
      phi_m = -0.3_dp*params%t_phe_ev
      density = max(0.9_dp*params%n_swi_inf_m3, tiny(1.0_dp))
      call encode_unknowns(params, branch, phi0, phi_m, density, guesses(:, 4), valid)
      call encode_unknowns( &
        params, branch, 0.8_dp*params%t_phe_ev, -0.1_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 5), valid &
        )
    case ('B')
      count = 5
      call encode_unknowns(params, branch, 1.3_dp, 1.3_dp, 7.0e6_dp, guesses(:, 1), valid)
      call encode_unknowns(params, branch, 0.8_dp, 0.8_dp, 6.5e6_dp, guesses(:, 2), valid)
      call encode_unknowns(params, branch, 2.0_dp, 2.0_dp, 7.8e6_dp, guesses(:, 3), valid)
      call encode_unknowns( &
        params, branch, 0.6_dp*params%t_phe_ev, 0.6_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 4), valid &
        )
      call encode_unknowns( &
        params, branch, 0.2_dp*params%t_phe_ev, 0.2_dp*params%t_phe_ev, &
        0.8_dp*params%n_swi_inf_m3, guesses(:, 5), valid &
        )
    case ('C')
      count = 7
      call encode_unknowns(params, branch, -0.5_dp, -0.5_dp, 6.0e6_dp, guesses(:, 1), valid)
      call encode_unknowns(params, branch, -2.0_dp, -2.0_dp, 7.0e6_dp, guesses(:, 2), valid)
      call encode_unknowns(params, branch, -5.0_dp, -5.0_dp, 8.0e6_dp, guesses(:, 3), valid)
      call encode_unknowns(params, branch, -10.0_dp, -10.0_dp, 8.2e6_dp, guesses(:, 4), valid)
      call encode_unknowns(params, branch, -15.0_dp, -15.0_dp, 8.5e6_dp, guesses(:, 5), valid)
      call encode_unknowns( &
        params, branch, -1.0_dp*params%t_phe_ev, -1.0_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 6), valid &
        )
      call encode_unknowns( &
        params, branch, -3.0_dp*params%t_phe_ev, -3.0_dp*params%t_phe_ev, &
        0.9_dp*params%n_swi_inf_m3, guesses(:, 7), valid &
        )
    case default
      count = 0
    end select
  end subroutine make_branch_guesses

  subroutine encode_unknowns(params, branch, phi0_v, phi_m_v, density_m3, y, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, density_m3
    real(dp), intent(out) :: y(3)
    logical, intent(out) :: valid

    y = 0.0_dp
    valid = density_m3 > 0.0_dp .and. params%n_phe_ref_m3 > 0.0_dp .and. params%t_phe_ev > 0.0_dp
    if (.not. valid) return
    select case (branch)
    case ('A')
      valid = phi0_v > 0.0_dp .and. phi_m_v < 0.0_dp
      if (.not. valid) return
      y(1) = log(phi0_v/params%t_phe_ev)
      y(2) = log(-phi_m_v/params%t_phe_ev)
      y(3) = log(density_m3/params%n_phe_ref_m3)
    case ('B')
      valid = phi0_v > 0.0_dp
      if (.not. valid) return
      y(1) = log(phi0_v/params%t_phe_ev)
      y(2) = log(density_m3/params%n_phe_ref_m3)
    case ('C')
      valid = phi0_v < 0.0_dp
      if (.not. valid) return
      y(1) = log(-phi0_v/params%t_phe_ev)
      y(2) = log(density_m3/params%n_phe_ref_m3)
    case default
      valid = .false.
    end select
    valid = valid .and. all(ieee_is_finite(y))
  end subroutine encode_unknowns

  subroutine decode_unknowns(params, branch, y, phi0_v, phi_m_v, density_m3, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: y(3)
    real(dp), intent(out) :: phi0_v, phi_m_v, density_m3
    logical, intent(out) :: valid

    phi0_v = 0.0_dp
    phi_m_v = 0.0_dp
    density_m3 = 0.0_dp
    valid = all(ieee_is_finite(y))
    if (.not. valid) return
    if (y(1) < -50.0_dp .or. y(1) > log(200.0_dp)) then
      valid = .false.
      return
    end if
    select case (branch)
    case ('A')
      if (y(2) < -50.0_dp .or. y(2) > log(200.0_dp) .or. &
          y(3) < -30.0_dp .or. y(3) > log(1.0e6_dp)) then
        valid = .false.
        return
      end if
      phi0_v = params%t_phe_ev*exp(y(1))
      phi_m_v = -params%t_phe_ev*exp(y(2))
      density_m3 = params%n_phe_ref_m3*exp(y(3))
    case ('B')
      if (y(2) < -30.0_dp .or. y(2) > log(1.0e6_dp)) then
        valid = .false.
        return
      end if
      phi0_v = params%t_phe_ev*exp(y(1))
      phi_m_v = phi0_v
      density_m3 = params%n_phe_ref_m3*exp(y(2))
    case ('C')
      if (y(2) < -30.0_dp .or. y(2) > log(1.0e6_dp)) then
        valid = .false.
        return
      end if
      phi0_v = -params%t_phe_ev*exp(y(1))
      phi_m_v = phi0_v
      density_m3 = params%n_phe_ref_m3*exp(y(2))
    case default
      valid = .false.
    end select
    valid = valid .and. all(ieee_is_finite([phi0_v, phi_m_v, density_m3]))
  end subroutine decode_unknowns

  subroutine newton_branch(params, branch, target_field_hat, y0, y_out, final_norm, iterations, success)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y0(3)
    real(dp), intent(out) :: y_out(3), final_norm
    integer, intent(out) :: iterations
    logical, intent(out) :: success

    real(dp) :: y(3), f(3), jac(3, 3), delta(3), trial(3), trial_f(3)
    real(dp) :: norm, trial_norm, step
    integer :: n, iteration, backtrack
    logical :: valid, jacobian_ok, linear_ok, trial_valid

    if (branch == 'A') then
      n = 3
    else
      n = 2
    end if
    y = y0
    call eval_charge_residual(params, branch, target_field_hat, y, f, valid)
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
      call numerical_jacobian(params, branch, target_field_hat, y, f, n, jac, jacobian_ok)
      if (.not. jacobian_ok) exit
      call solve_small_system(jac, -f, n, delta, linear_ok)
      if (.not. linear_ok) exit
      step = 1.0_dp
      do backtrack = 1, root_max_backtracks
        trial = y + step*delta
        call eval_charge_residual(params, branch, target_field_hat, trial, trial_f, trial_valid)
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
  end subroutine newton_branch

  subroutine numerical_jacobian(params, branch, target_field_hat, y, f0, n, jac, success)
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
      call eval_charge_residual(params, branch, target_field_hat, yp, fp, plus_valid)
      call eval_charge_residual(params, branch, target_field_hat, ym, fm, minus_valid)
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
  end subroutine numerical_jacobian

  subroutine solve_small_system(a_in, b_in, n, x, success)
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
  end subroutine solve_small_system

  subroutine eval_charge_residual(params, branch, target_field_hat, y, residual, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y(3)
    real(dp), intent(out) :: residual(3)
    logical, intent(out) :: valid

    real(dp) :: phi0_v, phi_m_v, density_m3, phi0_hat, phi_m_hat, density_hat
    real(dp) :: raw(3), integral, field_squared, field_scale
    real(dp) :: x3(3), x2(2)
    logical :: integral_ok

    residual = 0.0_dp
    call decode_unknowns(params, branch, y, phi0_v, phi_m_v, density_m3, valid)
    if (.not. valid) return
    phi0_hat = phi0_v/params%t_phe_ev
    phi_m_hat = phi_m_v/params%t_phe_ev
    density_hat = density_m3/params%n_phe_ref_m3
    if (.not. ion_accessible(params, max(phi0_hat, 0.0_dp))) then
      valid = .false.
      return
    end if

    select case (branch)
    case ('A')
      x3 = [phi0_v, phi_m_v, density_m3]
      call zhao_residuals_type_a(params, x3, raw)
      call integrate_rho_hat( &
        params, branch, 'lower', phi_m_hat, phi0_hat, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) then
        valid = .false.
        return
      end if
      field_squared = -2.0_dp*integral
      if (field_squared < -1.0e-10_dp) then
        valid = .false.
        return
      end if
      field_scale = max(1.0_dp, target_field_hat*target_field_hat)
      residual(1) = raw(1)/params%n_phe_ref_m3
      residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_scale
      residual(3) = raw(3)
    case ('B')
      x2 = [phi0_v, density_m3]
      call zhao_residuals_type_b(params, x2, raw(1:2))
      call integrate_rho_hat( &
        params, branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) then
        valid = .false.
        return
      end if
      field_squared = 2.0_dp*integral
      if (field_squared < -1.0e-10_dp) then
        valid = .false.
        return
      end if
      field_scale = max(1.0_dp, target_field_hat*target_field_hat)
      residual(1) = raw(1)/params%n_phe_ref_m3
      residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_scale
    case ('C')
      x2 = [phi0_v, density_m3]
      call zhao_residuals_type_c(params, x2, raw(1:2))
      call integrate_rho_hat( &
        params, branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) then
        valid = .false.
        return
      end if
      field_squared = 2.0_dp*integral
      if (field_squared < -1.0e-10_dp) then
        valid = .false.
        return
      end if
      field_scale = max(1.0_dp, target_field_hat*target_field_hat)
      residual(1) = raw(1)/params%n_phe_ref_m3
      residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_scale
    case default
      valid = .false.
      return
    end select
    valid = all(ieee_is_finite(residual))
  end subroutine eval_charge_residual

  subroutine integrate_rho_hat( &
    params, branch, side, lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat, integral, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat
    real(dp), intent(out) :: integral
    logical, intent(out) :: success

    real(dp) :: t, phi_hat, jacobian, rho_hat, summand, weight, h
    integer :: point

    integral = 0.0_dp
    success = .false.
    if (rho_quadrature_panels < 2 .or. mod(rho_quadrature_panels, 2) /= 0) return
    if (.not. all(ieee_is_finite([ &
                                 lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat &
                                 ])) .or. density_hat <= 0.0_dp) return
    h = 1.0_dp/real(rho_quadrature_panels, dp)
    do point = 0, rho_quadrature_panels
      t = real(point, dp)*h
      phi_hat = lower_phi_hat + (upper_phi_hat - lower_phi_hat)*sin(0.5_dp*pi*t)**2
      jacobian = (upper_phi_hat - lower_phi_hat)*0.5_dp*pi*sin(pi*t)
      if (.not. ion_accessible(params, phi_hat)) return
      call evaluate_zhao_rho_hat( &
        params, branch, side, phi_hat, phi0_hat, phi_m_hat, density_hat, rho_hat &
        )
      if (.not. ieee_is_finite(rho_hat)) return
      summand = rho_hat*jacobian
      if (point == 0 .or. point == rho_quadrature_panels) then
        weight = 1.0_dp
      else if (mod(point, 2) == 0) then
        weight = 2.0_dp
      else
        weight = 4.0_dp
      end if
      integral = integral + weight*summand
    end do
    integral = integral*h/3.0_dp
    success = ieee_is_finite(integral)
  end subroutine integrate_rho_hat

  subroutine build_zhao_outer_profile(params, root, grid_points, state, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    integer(i32), intent(in) :: grid_points
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp), allocatable :: phi_hat(:), field_hat(:), rho_hat(:)
    real(dp) :: phi0_hat, phi_m_hat, density_hat, electron_current, ion_current, photo_current
    integer :: point
    logical :: success

    state = outer_plasma_state_type()
    state%model = 'kinetic_1d'
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_zhao_params(params) .or. grid_points < 5_i32 .or. &
        .not. ieee_is_finite(root%interface_field_v_m) .or. root%n_swe_inf_m3 <= 0.0_dp) then
      state%applicability_status = status
      message = 'invalid Zhao profile request'
      return
    end if

    allocate (phi_hat(grid_points), field_hat(grid_points), rho_hat(grid_points))
    allocate (state%z(grid_points), state%potential(grid_points), &
              state%field(grid_points), state%charge_density(grid_points))
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    density_hat = root%n_swe_inf_m3/params%n_phe_ref_m3
    if (root%branch == '0' .or. &
        (abs(phi0_hat) <= zero_field_tolerance_hat .and. &
         abs(root%interface_field_v_m) <= &
         zero_field_tolerance_hat*params%t_phe_ev/params%lambda_d_phe_ref_m)) then
      phi_hat = 0.0_dp
      field_hat = 0.0_dp
      rho_hat = 0.0_dp
      do point = 1, grid_points
        state%z(point) = params%lambda_d_phe_ref_m*real(point - 1, dp)/real(grid_points - 1, dp)
      end do
    else
      select case (root%branch)
      case ('A')
        call build_type_a_profile( &
          params, phi0_hat, phi_m_hat, density_hat, phi_hat, field_hat, rho_hat, state%z, success &
          )
      case ('B', 'C')
        call build_monotonic_profile( &
          params, root%branch, phi0_hat, density_hat, phi_hat, field_hat, rho_hat, state%z, success &
          )
      case default
        success = .false.
      end select
      if (.not. success) then
        status = outer_plasma_no_physical_solution
        state%applicability_status = status
        message = 'failed to reconstruct the charge-driven Zhao profile'
        return
      end if
    end if

    state%profile_n = grid_points
    state%kinetic_closure = 'zhao_charge_driven'
    state%zhao_branch = root%branch
    state%interface_z = 0.0_dp
    state%interface_potential = root%phi0_v
    state%zhao_phi0 = root%phi0_v
    state%zhao_phi_minimum = root%phi_m_v
    state%zhao_electron_density_infinity = root%n_swe_inf_m3
    state%infinity_potential = 0.0_dp
    state%debye_length = params%lambda_d_phe_ref_m
    state%interface_field = root%interface_field_v_m
    state%nonlinear_iterations = root%nonlinear_iterations
    state%nonlinear_residual = root%residual_norm
    state%potential = params%t_phe_ev*phi_hat
    state%field = (params%t_phe_ev/params%lambda_d_phe_ref_m)*field_hat
    state%field(1) = root%interface_field_v_m
    state%charge_density = qe*params%n_phe_ref_m3*rho_hat
    state%integrated_charge_per_area = trapz(state%z, state%charge_density)
    call zhao_current_components(params, root, electron_current, ion_current, photo_current)
    state%electron_current_density = electron_current
    state%ion_current_density = ion_current
    state%photoelectron_current_density = photo_current
    state%total_current_density = electron_current + ion_current + photo_current
    state%ready = .true.
    state%applicability_status = outer_plasma_ok
    status = outer_plasma_ok
  end subroutine build_zhao_outer_profile

  subroutine build_monotonic_profile( &
    params, branch, phi0_hat, density_hat, phi_hat, field_hat, rho_hat, z, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_hat, density_hat
    real(dp), intent(out) :: phi_hat(:), field_hat(:), rho_hat(:), z(:)
    logical, intent(out) :: success

    real(dp) :: x, phi_end_hat, integral, field_squared, phi_m_hat
    integer :: point, n
    logical :: integral_ok

    success = .false.
    n = size(phi_hat)
    if (size(field_hat) /= n .or. size(rho_hat) /= n .or. size(z) /= n) return
    phi_m_hat = phi0_hat
    if (branch == 'B') then
      phi_end_hat = profile_phi_end_hat
      if (phi_end_hat >= phi0_hat) phi_end_hat = 0.5_dp*phi0_hat
    else if (branch == 'C') then
      phi_end_hat = -profile_phi_end_hat
      if (phi_end_hat <= phi0_hat) phi_end_hat = 0.5_dp*phi0_hat
    else
      return
    end if
    do point = 1, n
      x = real(point - 1, dp)/real(n - 1, dp)
      phi_hat(point) = phi0_hat + (phi_end_hat - phi0_hat)*sin(0.5_dp*pi*x)**2
      call integrate_rho_hat( &
        params, branch, 'monotonic', phi_hat(point), 0.0_dp, phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) return
      field_squared = 2.0_dp*integral
      if (field_squared < -1.0e-9_dp) return
      if (branch == 'B') then
        field_hat(point) = sqrt(max(0.0_dp, field_squared))
      else
        field_hat(point) = -sqrt(max(0.0_dp, field_squared))
      end if
      if (.not. ion_accessible(params, phi_hat(point))) return
      call evaluate_zhao_rho_hat( &
        params, branch, 'monotonic', phi_hat(point), phi0_hat, phi_m_hat, density_hat, rho_hat(point) &
        )
    end do
    call integrate_profile_z(params%lambda_d_phe_ref_m, phi_hat, field_hat, z, success)
  end subroutine build_monotonic_profile

  subroutine build_type_a_profile( &
    params, phi0_hat, phi_m_hat, density_hat, phi_hat, field_hat, rho_hat, z, success &
    )
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: phi0_hat, phi_m_hat, density_hat
    real(dp), intent(out) :: phi_hat(:), field_hat(:), rho_hat(:), z(:)
    logical, intent(out) :: success

    real(dp) :: x, phi_end_hat, integral, field_squared
    integer :: point, upper_point, n, n_lower, n_upper
    logical :: integral_ok

    success = .false.
    n = size(phi_hat)
    if (size(field_hat) /= n .or. size(rho_hat) /= n .or. size(z) /= n) return
    n_lower = max(3, (n + 1)/3)
    n_upper = n - n_lower + 1
    if (n_upper < 3) return
    phi_end_hat = -profile_phi_end_hat
    if (phi_end_hat <= phi_m_hat) phi_end_hat = 0.5_dp*phi_m_hat

    do point = 1, n_lower
      x = real(point - 1, dp)/real(n_lower - 1, dp)
      phi_hat(point) = phi_m_hat + (phi0_hat - phi_m_hat)*cos(0.5_dp*pi*x)**2
      call integrate_rho_hat( &
        params, 'A', 'lower', phi_m_hat, phi_hat(point), phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) return
      field_squared = -2.0_dp*integral
      if (field_squared < -1.0e-9_dp) return
      field_hat(point) = sqrt(max(0.0_dp, field_squared))
      call evaluate_zhao_rho_hat( &
        params, 'A', 'lower', phi_hat(point), phi0_hat, phi_m_hat, density_hat, rho_hat(point) &
        )
    end do
    do upper_point = 2, n_upper
      point = n_lower + upper_point - 1
      x = real(upper_point - 1, dp)/real(n_upper - 1, dp)
      phi_hat(point) = phi_m_hat + (phi_end_hat - phi_m_hat)*sin(0.5_dp*pi*x)**2
      call integrate_rho_hat( &
        params, 'A', 'upper', phi_m_hat, phi_hat(point), phi0_hat, phi_m_hat, density_hat, &
        integral, integral_ok &
        )
      if (.not. integral_ok) return
      field_squared = -2.0_dp*integral
      if (field_squared < -1.0e-9_dp) return
      field_hat(point) = -sqrt(max(0.0_dp, field_squared))
      call evaluate_zhao_rho_hat( &
        params, 'A', 'upper', phi_hat(point), phi0_hat, phi_m_hat, density_hat, rho_hat(point) &
        )
    end do
    call integrate_profile_z(params%lambda_d_phe_ref_m, phi_hat, field_hat, z, success)
  end subroutine build_type_a_profile

  subroutine integrate_profile_z(lambda_d, phi_hat, field_hat, z, success)
    real(dp), intent(in) :: lambda_d, phi_hat(:), field_hat(:)
    real(dp), intent(out) :: z(:)
    logical, intent(out) :: success

    real(dp) :: field_average, delta_phi
    integer :: point

    success = .false.
    if (lambda_d <= 0.0_dp .or. size(phi_hat) /= size(field_hat) .or. size(z) /= size(phi_hat)) return
    z(1) = 0.0_dp
    do point = 2, size(phi_hat)
      delta_phi = abs(phi_hat(point) - phi_hat(point - 1))
      field_average = 0.5_dp*(abs(field_hat(point)) + abs(field_hat(point - 1)))
      if (field_average <= tiny(1.0_dp) .and. delta_phi > tiny(1.0_dp)) return
      z(point) = z(point - 1) + lambda_d*delta_phi/max(field_average, tiny(1.0_dp))
    end do
    success = all(ieee_is_finite(z)) .and. all(z(2:) > z(:size(z) - 1))
  end subroutine integrate_profile_z

  real(dp) function zhao_net_current_density(params, root) result(current_density)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root

    real(dp) :: electron_current, ion_current, photo_current

    call zhao_current_components(params, root, electron_current, ion_current, photo_current)
    current_density = electron_current + ion_current + photo_current
  end function zhao_net_current_density

  subroutine zhao_current_components(params, root, electron_current, ion_current, photo_current)
    type(zhao_params_type), intent(in) :: params
    type(zhao_charge_root_type), intent(in) :: root
    real(dp), intent(out) :: electron_current, ion_current, photo_current

    real(dp) :: a_swe, electron_term, ion_term, photo_term, scale
    real(dp) :: phi0_hat, phi_m_hat
    character(len=1) :: branch

    electron_current = 0.0_dp
    ion_current = 0.0_dp
    photo_current = 0.0_dp
    if (.not. valid_zhao_params(params) .or. root%n_swe_inf_m3 <= 0.0_dp) return
    phi0_hat = root%phi0_v/params%t_phe_ev
    phi_m_hat = root%phi_m_v/params%t_phe_ev
    branch = root%branch
    if (branch == '0') branch = 'B'
    ion_term = params%n_swi_inf_m3*sqrt( &
               2.0_dp*pi*params%t_swe_ev/params%t_phe_ev*params%m_e_kg/params%m_i_kg &
               )*params%mach
    select case (branch)
    case ('A')
      a_swe = sqrt(max(0.0_dp, -phi_m_hat/params%tau)) - params%u
      electron_term = swe_free_current_term(params, root%n_swe_inf_m3, a_swe)
      photo_term = params%n_phe0_m3*exp(phi_m_hat - phi0_hat)
    case ('B')
      electron_term = swe_free_current_term(params, root%n_swe_inf_m3, -params%u)
      photo_term = params%n_phe0_m3*exp(-phi0_hat)
    case ('C')
      a_swe = sqrt(max(0.0_dp, -phi0_hat/params%tau)) - params%u
      electron_term = swe_free_current_term(params, root%n_swe_inf_m3, a_swe)
      photo_term = params%n_phe0_m3
    case default
      return
    end select
    scale = qe*params%v_phe_th_mps/(2.0_dp*sqrt(pi))
    electron_current = -scale*electron_term
    ion_current = scale*ion_term
    photo_current = scale*photo_term
  end subroutine zhao_current_components

  pure logical function ion_accessible(params, phi_hat) result(accessible)
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: phi_hat

    accessible = params%tau > 0.0_dp .and. params%mach > 0.0_dp .and. &
                 1.0_dp - 2.0_dp*phi_hat/(params%tau*params%mach*params%mach) > 0.0_dp
  end function ion_accessible

  pure logical function valid_zhao_params(params) result(valid)
    type(zhao_params_type), intent(in) :: params

    valid = all(ieee_is_finite([ &
                               params%n_swi_inf_m3, params%n_phe_ref_m3, params%n_phe0_m3, &
                               params%t_swe_ev, params%t_phe_ev, params%m_i_kg, params%m_e_kg, &
                               params%v_phe_th_mps, params%mach, params%u, params%tau, &
                               params%lambda_d_phe_ref_m &
                               ])) .and. &
            params%n_swi_inf_m3 > 0.0_dp .and. params%n_phe_ref_m3 > 0.0_dp .and. &
            params%t_swe_ev > 0.0_dp .and. params%t_phe_ev > 0.0_dp .and. &
            params%m_i_kg > 0.0_dp .and. params%m_e_kg > 0.0_dp .and. &
            params%v_phe_th_mps > 0.0_dp .and. params%mach > 0.0_dp .and. &
            params%tau > 0.0_dp .and. params%lambda_d_phe_ref_m > 0.0_dp
  end function valid_zhao_params

  pure real(dp) function trapz(x, y) result(integral)
    real(dp), intent(in) :: x(:), y(:)
    integer :: point

    integral = 0.0_dp
    do point = 1, size(x) - 1
      integral = integral + 0.5_dp*(y(point) + y(point + 1))*(x(point + 1) - x(point))
    end do
  end function trapz

end module bem_outer_plasma_zhao
