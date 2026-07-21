module bem_outer_plasma_kinetic
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: pi, eps0, qe
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid, &
                                    outer_plasma_no_physical_solution, outer_plasma_numerical_failure
  use bem_outer_plasma_grid, only: outer_plasma_grid_type, init_outer_plasma_grid
  use bem_outer_plasma_zhao, only: zhao_charge_root_type, solve_outer_plasma_zhao, &
                                   solve_outer_plasma_zhao_column
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  type, public :: kinetic_outer_plasma_options_type
    character(len=32) :: kinetic_closure = 'absorbing_maxwellian'
    character(len=16) :: zhao_branch = 'auto'
    integer(i32) :: grid_points = 128_i32
    integer(i32) :: max_iterations = 40_i32
    real(dp) :: domain_length = 0.0_dp
    real(dp) :: grid_stretch = 2.0_dp
    real(dp) :: tail_length = 0.0_dp
    real(dp) :: interface_field = 0.0_dp
    real(dp) :: electron_charge = 0.0_dp
    real(dp) :: electron_mass = 0.0_dp
    real(dp) :: electron_density_infinity = 0.0_dp
    real(dp) :: electron_temperature_j = 0.0_dp
    real(dp) :: electron_drift_infinity = 0.0_dp
    real(dp) :: ion_charge = 0.0_dp
    real(dp) :: ion_mass = 0.0_dp
    real(dp) :: ion_density_infinity = 0.0_dp
    real(dp) :: ion_temperature_j = 0.0_dp
    real(dp) :: ion_gamma = 1.0_dp
    real(dp) :: ion_drift_infinity = 0.0_dp
    real(dp) :: photoelectron_charge = 0.0_dp
    real(dp) :: photoelectron_mass = 0.0_dp
    real(dp) :: photoelectron_temperature_j = 0.0_dp
    real(dp) :: photoelectron_emission_flux = 0.0_dp
    real(dp) :: photoelectron_reference_density = 0.0_dp
    real(dp) :: photoelectron_population_fraction = 1.0_dp
    real(dp) :: photoelectron_source_scale = 1.0_dp
    logical :: photoelectron_column_closure_enabled = .false.
    real(dp) :: photoelectron_column_target_m2 = 0.0_dp
    real(dp) :: zhao_alpha_deg = 60.0_dp
    real(dp) :: residual_tolerance = 1.0e-8_dp
    real(dp) :: external_current_density = 0.0_dp
  end type kinetic_outer_plasma_options_type

  public :: eval_absorbing_maxwellian_density
  public :: eval_cold_drift_density
  public :: kinetic_bohm_speed
  public :: eval_photoelectron_escape_return
  public :: eval_emitted_maxwellian_density
  public :: eval_kinetic_current_balance
  public :: eval_kinetic_residual_jacobian_action
  public :: solve_outer_plasma_kinetic

contains

  subroutine solve_outer_plasma_kinetic( &
    options, state, status, message, initial_potential, continuation_steps, initial_state &
    )
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: initial_potential(:)
    integer(i32), intent(out), optional :: continuation_steps
    type(outer_plasma_state_type), intent(in), optional :: initial_state
    type(kinetic_outer_plasma_options_type) :: step_options
    type(outer_plasma_state_type) :: step_state
    type(outer_plasma_grid_type) :: grid
    real(dp), allocatable :: current_profile(:)
    real(dp) :: current_field, target_field, remaining, increment, next_field
    real(dp) :: field_tolerance, minimum_increment
    integer(i32) :: successful_steps, attempts

    if (present(continuation_steps)) continuation_steps = 0_i32
    select case (trim(lower_ascii(options%kinetic_closure)))
    case ('zhao_charge_driven')
      if (present(initial_state)) then
        call solve_outer_plasma_zhao_closure(options, state, status, message, initial_state)
      else
        call solve_outer_plasma_zhao_closure(options, state, status, message)
      end if
      if (status == outer_plasma_ok .and. present(continuation_steps)) continuation_steps = 1_i32
      return
    case ('absorbing_maxwellian')
      continue
    case default
      state = outer_plasma_state_type()
      state%model = 'kinetic_1d'
      state%kinetic_closure = options%kinetic_closure
      state%applicability_status = outer_plasma_invalid
      status = outer_plasma_invalid
      message = 'unknown kinetic outer-plasma closure'
      return
    end select
    if (.not. present(initial_potential)) then
      call solve_outer_plasma_kinetic_fixed(options, state, status, message)
      return
    end if
    if (options%grid_points < 3_i32 .or. size(initial_potential) /= options%grid_points .or. &
        options%domain_length <= 0.0_dp .or. options%grid_stretch < 0.0_dp) then
      call solve_outer_plasma_kinetic_fixed( &
        options, state, status, message, initial_potential=initial_potential &
        )
      return
    end if

    call init_outer_plasma_grid(options%grid_points, options%domain_length, options%grid_stretch, grid)
    allocate (current_profile(grid%n))
    current_profile = initial_potential
    current_field = -(current_profile(2) - current_profile(1))/grid%dz(1)
    target_field = options%interface_field
    field_tolerance = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(current_field), abs(target_field))
    remaining = target_field - current_field
    if (abs(remaining) <= field_tolerance) then
      call solve_outer_plasma_kinetic_fixed( &
        options, state, status, message, initial_potential=current_profile &
        )
      if (status == outer_plasma_ok .and. present(continuation_steps)) continuation_steps = 1_i32
      return
    end if

    minimum_increment = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(current_field), abs(target_field))
    increment = remaining
    successful_steps = 0_i32
    attempts = 0_i32
    do while (attempts < 128_i32)
      attempts = attempts + 1_i32
      remaining = target_field - current_field
      if (abs(increment) >= abs(remaining)) then
        next_field = target_field
      else
        next_field = current_field + increment
      end if
      step_options = options
      step_options%interface_field = next_field
      call solve_outer_plasma_kinetic_fixed( &
        step_options, step_state, status, message, initial_potential=current_profile &
        )
      if (status == outer_plasma_ok) then
        current_profile = step_state%potential
        current_field = next_field
        successful_steps = successful_steps + 1_i32
        if (abs(target_field - current_field) <= field_tolerance) then
          state = step_state
          if (present(continuation_steps)) continuation_steps = successful_steps
          return
        end if
        remaining = target_field - current_field
        increment = sign(min(2.0_dp*abs(increment), abs(remaining)), remaining)
      else if (status == outer_plasma_numerical_failure) then
        increment = 0.5_dp*increment
        if (abs(increment) < minimum_increment) then
          state = step_state
          message = 'kinetic interface-field continuation exhausted its minimum increment: '//trim(message)
          return
        end if
      else
        state = step_state
        return
      end if
    end do

    state = step_state
    status = outer_plasma_numerical_failure
    state%applicability_status = status
    message = 'kinetic interface-field continuation reached its attempt limit'
  end subroutine solve_outer_plasma_kinetic

  subroutine solve_outer_plasma_zhao_closure(options, state, status, message, initial_state)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(outer_plasma_state_type), intent(in), optional :: initial_state

    type(zhao_params_type) :: params
    type(zhao_charge_root_type) :: root, initial_root
    real(dp) :: electron_temperature_ev, photoelectron_temperature_ev, charge_tolerance
    logical :: has_initial_root

    state = outer_plasma_state_type()
    state%model = 'kinetic_1d'
    state%kinetic_closure = 'zhao_charge_driven'
    status = outer_plasma_invalid
    message = ''
    charge_tolerance = 0.1_dp*qe
    if (options%grid_points < 5_i32 .or. options%ion_density_infinity <= 0.0_dp .or. &
        options%photoelectron_reference_density <= 0.0_dp .or. &
        options%electron_temperature_j <= 0.0_dp .or. options%photoelectron_temperature_j <= 0.0_dp .or. &
        options%electron_drift_infinity <= 0.0_dp .or. options%ion_drift_infinity <= 0.0_dp .or. &
        options%electron_mass <= 0.0_dp .or. options%ion_mass <= 0.0_dp .or. &
        options%photoelectron_mass <= 0.0_dp .or. &
        .not. ieee_is_finite(options%photoelectron_population_fraction) .or. &
        options%photoelectron_population_fraction < 0.0_dp .or. &
        (options%photoelectron_column_closure_enabled .and. &
         (.not. ieee_is_finite(options%photoelectron_column_target_m2) .or. &
          options%photoelectron_column_target_m2 < 0.0_dp .or. &
          .not. ieee_is_finite(options%domain_length) .or. options%domain_length <= 0.0_dp))) then
      state%applicability_status = status
      message = 'Zhao charge-driven closure parameters are invalid'
      return
    end if
    if (abs(abs(options%electron_charge) - qe) > charge_tolerance .or. &
        abs(options%ion_charge - qe) > charge_tolerance .or. &
        abs(abs(options%photoelectron_charge) - qe) > charge_tolerance) then
      state%applicability_status = status
      message = 'Zhao charge-driven closure requires singly charged physical species'
      return
    end if
    if (abs(options%external_current_density) > tiny(1.0_dp)) then
      state%applicability_status = status
      message = 'Zhao charge-driven closure does not support external_current_density'
      return
    end if

    electron_temperature_ev = options%electron_temperature_j/qe
    photoelectron_temperature_ev = options%photoelectron_temperature_j/qe
    call build_zhao_params( &
      options%zhao_alpha_deg, options%ion_density_infinity, options%photoelectron_reference_density, &
      electron_temperature_ev, photoelectron_temperature_ev, options%electron_drift_infinity, &
      options%ion_drift_infinity, options%ion_mass, options%electron_mass, params, &
      photoelectron_population_fraction=options%photoelectron_population_fraction, &
      photoelectron_source_scale=options%photoelectron_source_scale &
      )
    has_initial_root = .false.
    if (present(initial_state)) then
      has_initial_root = initial_state%ready .and. &
                         trim(lower_ascii(initial_state%model)) == 'kinetic_1d' .and. &
                         trim(lower_ascii(initial_state%kinetic_closure)) == 'zhao_charge_driven' .and. &
                         index('ABC0', initial_state%zhao_branch) > 0 .and. &
                         initial_state%zhao_electron_density_infinity > 0.0_dp
      if (has_initial_root) then
        initial_root%branch = initial_state%zhao_branch
        initial_root%phi0_v = initial_state%zhao_phi0
        initial_root%phi_m_v = initial_state%zhao_phi_minimum
        initial_root%n_swe_inf_m3 = initial_state%zhao_electron_density_infinity
        initial_root%interface_field_v_m = initial_state%interface_field
        initial_root%photoelectron_population_fraction = initial_state%photoelectron_population_fraction
      end if
    end if
    if (options%photoelectron_column_closure_enabled) then
      if (has_initial_root) then
        call solve_outer_plasma_zhao_column( &
          options%zhao_branch, params, options%interface_field, options%grid_points, &
          options%domain_length, options%photoelectron_column_target_m2, state, root, status, message, &
          initial_root=initial_root &
          )
      else
        call solve_outer_plasma_zhao_column( &
          options%zhao_branch, params, options%interface_field, options%grid_points, &
          options%domain_length, options%photoelectron_column_target_m2, state, root, status, message &
          )
      end if
    else if (has_initial_root) then
      call solve_outer_plasma_zhao( &
        options%zhao_branch, params, options%interface_field, options%grid_points, state, root, status, message, &
        initial_root=initial_root &
        )
    else
      call solve_outer_plasma_zhao( &
        options%zhao_branch, params, options%interface_field, options%grid_points, state, root, status, message, &
        allow_transient_bootstrap=.true. &
        )
    end if
    state%model = 'kinetic_1d'
    state%kinetic_closure = 'zhao_charge_driven'
    state%zhao_branch = root%branch
    state%zhao_phi0 = root%phi0_v
    state%zhao_phi_minimum = root%phi_m_v
    state%zhao_electron_density_infinity = root%n_swe_inf_m3
    state%photoelectron_population_fraction = root%photoelectron_population_fraction
    state%applicability_status = status
  end subroutine solve_outer_plasma_zhao_closure

  subroutine solve_outer_plasma_kinetic_fixed(options, state, status, message, initial_potential)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: initial_potential(:)
    type(outer_plasma_grid_type) :: grid
    real(dp), allocatable :: phi(:), residual(:), trial_residual(:), delta(:), trial(:)
    real(dp), allocatable :: lower(:), diagonal(:), regularized_diagonal(:), upper(:), border(:)
    real(dp) :: residual_norm, trial_norm, bohm_speed, seed_scale, pseudo_time
    real(dp) :: previous_interface_field, field_correction
    integer(i32) :: iteration, solve_status, recovery_attempt
    logical :: accepted

    state = outer_plasma_state_type()
    state%model = 'kinetic_1d'
    status = outer_plasma_invalid
    message = ''
    if (.not. valid_kinetic_options(options, message)) then
      state%applicability_status = status
      return
    end if
    bohm_speed = kinetic_bohm_speed( &
                 options%electron_temperature_j, options%ion_temperature_j, options%ion_mass, options%ion_gamma &
                 )
    if (options%ion_density_infinity > 0.0_dp .and. options%ion_drift_infinity < bohm_speed) then
      status = outer_plasma_no_physical_solution
      state%applicability_status = status
      message = 'ion infinity drift violates the kinetic Bohm entry condition'
      return
    end if

    call init_outer_plasma_grid(options%grid_points, options%domain_length, options%grid_stretch, grid)
    allocate (phi(grid%n), residual(grid%n), trial_residual(grid%n), delta(grid%n), trial(grid%n), &
              lower(grid%n), diagonal(grid%n), regularized_diagonal(grid%n), upper(grid%n), border(grid%n))
    if (present(initial_potential)) then
      if (size(initial_potential) /= grid%n) then
        message = 'kinetic initial profile size does not match the grid'
        state%applicability_status = status
        return
      end if
      phi = initial_potential
      previous_interface_field = -(phi(2) - phi(1))/grid%dz(1)
      field_correction = options%interface_field - previous_interface_field
      phi = phi + field_correction*(options%tail_length + options%domain_length - grid%z)
      if (.not. monotonic_electron_repelling_branch(options, phi)) then
        phi = options%interface_field*(options%tail_length + options%domain_length - grid%z)
      end if
    else
      phi = options%interface_field*(options%tail_length + options%domain_length - grid%z)
      if (abs(options%interface_field) <= 256.0_dp*epsilon(1.0_dp) .and. &
          options%photoelectron_emission_flux > 0.0_dp) then
        seed_scale = 1.0e-3_dp*options%photoelectron_temperature_j/abs(options%photoelectron_charge)
        phi = -seed_scale*(options%tail_length + options%domain_length - grid%z)/ &
              (options%tail_length + options%domain_length)
      end if
    end if
    call kinetic_residual(options, grid, phi, residual, solve_status)
    if (solve_status /= outer_plasma_ok) then
      status = solve_status
      state%applicability_status = status
      message = 'initial kinetic profile is outside the monotonic physical branch'
      return
    end if
    residual_norm = maxval(abs(residual))
    pseudo_time = max(minval(grid%dz)**2, tiny(1.0_dp))

    do iteration = 0_i32, options%max_iterations
      if (residual_norm <= options%residual_tolerance) exit
      if (iteration == options%max_iterations) then
        status = outer_plasma_numerical_failure
        state%applicability_status = status
        message = 'kinetic Newton solve reached max_iterations'
        return
      end if
      call kinetic_residual_jacobian(options, grid, phi, residual, lower, diagonal, upper, border, solve_status)
      if (solve_status /= outer_plasma_ok) then
        status = solve_status
        state%applicability_status = status
        message = 'kinetic Jacobian evaluation failed'
        return
      end if
      call solve_bordered_tridiagonal(lower, diagonal, upper, border, -residual, delta, accepted)
      if (accepted) then
        call accept_kinetic_step( &
          options, grid, phi, residual_norm, delta, trial, trial_residual, trial_norm, accepted &
          )
      end if
      if (.not. accepted) then
        do recovery_attempt = 0_i32, 16_i32
          regularized_diagonal = diagonal
          regularized_diagonal(2:grid%n - 1_i32) = &
            regularized_diagonal(2:grid%n - 1_i32) - 1.0_dp/pseudo_time
          call solve_bordered_tridiagonal( &
            lower, regularized_diagonal, upper, border, -residual, delta, accepted &
            )
          if (accepted) then
            call accept_kinetic_step( &
              options, grid, phi, residual_norm, delta, trial, trial_residual, trial_norm, accepted &
              )
          end if
          if (accepted) then
            pseudo_time = 4.0_dp*pseudo_time
            exit
          end if
          pseudo_time = 0.25_dp*pseudo_time
        end do
        if (.not. accepted) then
          status = outer_plasma_numerical_failure
          state%applicability_status = status
          message = 'kinetic Newton and pseudo-transient recovery failed on the monotonic branch'
          return
        end if
      end if
      phi = trial
      residual = trial_residual
      residual_norm = trial_norm
    end do

    call populate_kinetic_state(options, grid, phi, iteration, residual_norm, state, solve_status)
    if (solve_status /= outer_plasma_ok) then
      status = solve_status
      state%applicability_status = status
      message = 'converged kinetic profile failed density reconstruction'
      return
    end if
    status = outer_plasma_ok
    state%applicability_status = status
    state%ready = .true.
  end subroutine solve_outer_plasma_kinetic_fixed

  subroutine accept_kinetic_step( &
    options, grid, phi, residual_norm, delta, trial, trial_residual, trial_norm, accepted &
    )
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_grid_type), intent(in) :: grid
    real(dp), intent(in) :: phi(:), residual_norm, delta(:)
    real(dp), intent(out) :: trial(:), trial_residual(:), trial_norm
    logical, intent(out) :: accepted
    real(dp) :: step
    integer(i32) :: solve_status

    step = 1.0_dp
    accepted = .false.
    trial_norm = huge(1.0_dp)
    do while (step >= 2.0_dp**(-20))
      trial = phi + step*delta
      if (monotonic_electron_repelling_branch(options, trial)) then
        call kinetic_residual(options, grid, trial, trial_residual, solve_status)
        if (solve_status == outer_plasma_ok) then
          trial_norm = maxval(abs(trial_residual))
          if (trial_norm < residual_norm) then
            accepted = .true.
            return
          end if
        end if
      end if
      step = 0.5_dp*step
    end do
  end subroutine accept_kinetic_step

  logical function valid_kinetic_options(options, message) result(valid)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    character(len=*), intent(out) :: message
    real(dp) :: charge_scale

    valid = .false.
    message = ''
    if (options%grid_points < 3_i32 .or. options%max_iterations < 1_i32 .or. &
        options%domain_length <= 0.0_dp .or. options%tail_length <= 0.0_dp .or. &
        options%grid_stretch < 0.0_dp .or. options%residual_tolerance <= 0.0_dp) then
      message = 'kinetic grid and solver controls must be positive'
      return
    end if
    if (options%electron_charge >= 0.0_dp .or. options%electron_mass <= 0.0_dp .or. &
        options%electron_density_infinity < 0.0_dp .or. options%electron_temperature_j <= 0.0_dp .or. &
        options%ion_charge <= 0.0_dp .or. options%ion_mass <= 0.0_dp .or. &
        options%ion_density_infinity < 0.0_dp .or. options%ion_temperature_j < 0.0_dp .or. &
        options%ion_gamma < 0.0_dp .or. options%ion_drift_infinity <= 0.0_dp) then
      message = 'kinetic ambient species parameters are invalid'
      return
    end if
    if (options%photoelectron_emission_flux < 0.0_dp .or. &
        (options%photoelectron_emission_flux > 0.0_dp .and. &
         (options%photoelectron_charge >= 0.0_dp .or. options%photoelectron_mass <= 0.0_dp .or. &
          options%photoelectron_temperature_j <= 0.0_dp))) then
      message = 'kinetic photoelectron source parameters are invalid'
      return
    end if
    charge_scale = max(abs(options%electron_charge*options%electron_density_infinity), &
                       abs(options%ion_charge*options%ion_density_infinity), tiny(1.0_dp))
    if (abs(options%electron_charge*options%electron_density_infinity + &
            options%ion_charge*options%ion_density_infinity) > 1.0e-10_dp*charge_scale) then
      message = 'kinetic infinity state must be quasineutral for the Robin tail'
      return
    end if
    valid = .true.
  end function valid_kinetic_options

  subroutine kinetic_residual(options, grid, phi, residual, status)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_grid_type), intent(in) :: grid
    real(dp), intent(in) :: phi(:)
    real(dp), intent(out) :: residual(:)
    integer(i32), intent(out) :: status
    real(dp) :: electron_density, ion_density, photoelectron_density, susceptibility, speed, left_h, right_h
    integer(i32) :: j, closure_status

    status = outer_plasma_no_physical_solution
    residual = 0.0_dp
    if (.not. monotonic_electron_repelling_branch(options, phi)) return
    residual(1) = (phi(2) - phi(1))/grid%dz(1) + options%interface_field
    do j = 2_i32, grid%n - 1_i32
      call eval_absorbing_maxwellian_density( &
        phi(j), phi(1), options%electron_charge, options%electron_temperature_j, &
        options%electron_density_infinity, electron_density, susceptibility, closure_status &
        )
      if (closure_status /= outer_plasma_ok) return
      call eval_cold_drift_density( &
        phi(j), options%ion_charge, options%ion_mass, options%ion_density_infinity, &
        options%ion_drift_infinity, ion_density, speed, susceptibility, closure_status &
        )
      if (closure_status /= outer_plasma_ok) return
      photoelectron_density = 0.0_dp
      if (options%photoelectron_emission_flux > 0.0_dp) then
        call eval_emitted_maxwellian_density( &
          phi(j), phi(1), 0.0_dp, options%photoelectron_charge, options%photoelectron_mass, &
          options%photoelectron_temperature_j, options%photoelectron_emission_flux, &
          photoelectron_density, closure_status &
          )
        if (closure_status /= outer_plasma_ok) return
      end if
      left_h = grid%dz(j - 1_i32)
      right_h = grid%dz(j)
      residual(j) = 2.0_dp*((phi(j + 1_i32) - phi(j))/right_h - &
                            (phi(j) - phi(j - 1_i32))/left_h)/(left_h + right_h) + &
                    (options%electron_charge*electron_density + options%ion_charge*ion_density + &
                     options%photoelectron_charge*photoelectron_density)/eps0
    end do
    residual(grid%n) = (phi(grid%n) - phi(grid%n - 1_i32))/grid%dz(grid%n - 1_i32) + &
                       phi(grid%n)/options%tail_length
    if (.not. all(ieee_is_finite(residual))) then
      status = outer_plasma_numerical_failure
      return
    end if
    status = outer_plasma_ok
  end subroutine kinetic_residual

  subroutine kinetic_residual_jacobian(options, grid, phi, residual, lower, diagonal, upper, border, status)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_grid_type), intent(in) :: grid
    real(dp), intent(in) :: phi(:)
    real(dp), intent(out) :: residual(:), lower(:), diagonal(:), upper(:), border(:)
    integer(i32), intent(out) :: status
    real(dp) :: electron_density, ion_density, photoelectron_density, electron_local, electron_interface
    real(dp) :: ion_local, photoelectron_local, photoelectron_interface, speed, left_h, right_h, stencil_scale
    integer(i32) :: j, closure_status

    status = outer_plasma_no_physical_solution
    residual = 0.0_dp
    lower = 0.0_dp
    diagonal = 0.0_dp
    upper = 0.0_dp
    border = 0.0_dp
    if (.not. monotonic_electron_repelling_branch(options, phi)) return

    residual(1) = (phi(2) - phi(1))/grid%dz(1) + options%interface_field
    diagonal(1) = -1.0_dp/grid%dz(1)
    upper(1) = 1.0_dp/grid%dz(1)
    do j = 2_i32, grid%n - 1_i32
      call eval_absorbing_maxwellian_density( &
        phi(j), phi(1), options%electron_charge, options%electron_temperature_j, &
        options%electron_density_infinity, electron_density, electron_local, closure_status, electron_interface &
        )
      if (closure_status /= outer_plasma_ok) return
      call eval_cold_drift_density( &
        phi(j), options%ion_charge, options%ion_mass, options%ion_density_infinity, &
        options%ion_drift_infinity, ion_density, speed, ion_local, closure_status &
        )
      if (closure_status /= outer_plasma_ok) return
      photoelectron_density = 0.0_dp
      photoelectron_local = 0.0_dp
      photoelectron_interface = 0.0_dp
      if (options%photoelectron_emission_flux > 0.0_dp) then
        call eval_emitted_maxwellian_density( &
          phi(j), phi(1), 0.0_dp, options%photoelectron_charge, options%photoelectron_mass, &
          options%photoelectron_temperature_j, options%photoelectron_emission_flux, &
          photoelectron_density, closure_status, photoelectron_local, photoelectron_interface &
          )
        if (closure_status /= outer_plasma_ok) return
      end if
      left_h = grid%dz(j - 1_i32)
      right_h = grid%dz(j)
      stencil_scale = 2.0_dp/(left_h + right_h)
      residual(j) = stencil_scale*((phi(j + 1_i32) - phi(j))/right_h - &
                                   (phi(j) - phi(j - 1_i32))/left_h) + &
                    (options%electron_charge*electron_density + options%ion_charge*ion_density + &
                     options%photoelectron_charge*photoelectron_density)/eps0
      lower(j) = stencil_scale/left_h
      diagonal(j) = -stencil_scale*(1.0_dp/left_h + 1.0_dp/right_h) + &
                    (options%electron_charge*electron_local + options%ion_charge*ion_local + &
                     options%photoelectron_charge*photoelectron_local)/eps0
      upper(j) = stencil_scale/right_h
      border(j) = (options%electron_charge*electron_interface + &
                   options%photoelectron_charge*photoelectron_interface)/eps0
    end do
    residual(grid%n) = (phi(grid%n) - phi(grid%n - 1_i32))/grid%dz(grid%n - 1_i32) + &
                       phi(grid%n)/options%tail_length
    lower(grid%n) = -1.0_dp/grid%dz(grid%n - 1_i32)
    diagonal(grid%n) = 1.0_dp/grid%dz(grid%n - 1_i32) + 1.0_dp/options%tail_length
    if (.not. all(ieee_is_finite(residual)) .or. .not. all(ieee_is_finite(lower)) .or. &
        .not. all(ieee_is_finite(diagonal)) .or. .not. all(ieee_is_finite(upper)) .or. &
        .not. all(ieee_is_finite(border))) then
      status = outer_plasma_numerical_failure
      return
    end if
    status = outer_plasma_ok
  end subroutine kinetic_residual_jacobian

  subroutine eval_kinetic_residual_jacobian_action(options, phi, direction, residual, action, status)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    real(dp), intent(in) :: phi(:), direction(:)
    real(dp), intent(out) :: residual(:), action(:)
    integer(i32), intent(out) :: status
    type(outer_plasma_grid_type) :: grid
    real(dp), allocatable :: lower(:), diagonal(:), upper(:), border(:)
    integer(i32) :: j

    status = outer_plasma_invalid
    residual = 0.0_dp
    action = 0.0_dp
    if (size(phi) /= options%grid_points .or. size(direction) /= size(phi) .or. &
        size(residual) /= size(phi) .or. size(action) /= size(phi)) return
    call init_outer_plasma_grid(options%grid_points, options%domain_length, options%grid_stretch, grid)
    allocate (lower(grid%n), diagonal(grid%n), upper(grid%n), border(grid%n))
    call kinetic_residual_jacobian(options, grid, phi, residual, lower, diagonal, upper, border, status)
    if (status /= outer_plasma_ok) return
    action(1) = diagonal(1)*direction(1) + upper(1)*direction(2) + border(1)*direction(1)
    do j = 2_i32, grid%n - 1_i32
      action(j) = lower(j)*direction(j - 1_i32) + diagonal(j)*direction(j) + &
                  upper(j)*direction(j + 1_i32) + border(j)*direction(1)
    end do
    action(grid%n) = lower(grid%n)*direction(grid%n - 1_i32) + diagonal(grid%n)*direction(grid%n) + &
                     border(grid%n)*direction(1)
    if (.not. all(ieee_is_finite(action))) status = outer_plasma_numerical_failure
  end subroutine eval_kinetic_residual_jacobian_action

  logical function monotonic_electron_repelling_branch(options, phi) result(valid)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    real(dp), intent(in) :: phi(:)
    real(dp) :: tolerance

    tolerance = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(phi)))
    valid = all(ieee_is_finite(phi))
    if (.not. valid) return
    if (options%interface_field <= tolerance) then
      valid = options%electron_charge*phi(1) >= -tolerance .and. &
              all(phi(2:) >= phi(:size(phi) - 1) - tolerance)
    else
      valid = options%electron_charge*phi(1) <= tolerance .and. &
              all(phi(2:) <= phi(:size(phi) - 1) + tolerance)
    end if
  end function monotonic_electron_repelling_branch

  subroutine solve_bordered_tridiagonal(lower, diagonal, upper, border, rhs, solution, success)
    real(dp), intent(in) :: lower(:), diagonal(:), upper(:), border(:), rhs(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: success
    real(dp) :: base_solution(size(rhs)), border_solution(size(rhs)), denominator
    logical :: rhs_success, border_success

    call solve_tridiagonal(lower, diagonal, upper, rhs, base_solution, rhs_success)
    call solve_tridiagonal(lower, diagonal, upper, border, border_solution, border_success)
    success = rhs_success .and. border_success
    if (.not. success) return
    denominator = 1.0_dp + border_solution(1)
    if (.not. ieee_is_finite(denominator) .or. &
        abs(denominator) <= sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(border_solution(1)))) then
      success = .false.
      return
    end if
    solution = base_solution - border_solution*base_solution(1)/denominator
    success = all(ieee_is_finite(solution))
  end subroutine solve_bordered_tridiagonal

  subroutine solve_tridiagonal(lower, diagonal, upper, rhs, solution, success)
    real(dp), intent(in) :: lower(:), diagonal(:), upper(:), rhs(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: success
    real(dp) :: modified_upper(size(rhs)), modified_rhs(size(rhs)), pivot
    integer(i32) :: j, n

    n = size(rhs)
    solution = 0.0_dp
    modified_upper = 0.0_dp
    modified_rhs = 0.0_dp
    success = .false.
    pivot = diagonal(1)
    if (.not. ieee_is_finite(pivot) .or. abs(pivot) <= tiny(1.0_dp)) return
    modified_upper(1) = upper(1)/pivot
    modified_rhs(1) = rhs(1)/pivot
    do j = 2_i32, n
      pivot = diagonal(j) - lower(j)*modified_upper(j - 1_i32)
      if (.not. ieee_is_finite(pivot) .or. abs(pivot) <= tiny(1.0_dp)) return
      if (j < n) modified_upper(j) = upper(j)/pivot
      modified_rhs(j) = (rhs(j) - lower(j)*modified_rhs(j - 1_i32))/pivot
    end do
    solution(n) = modified_rhs(n)
    do j = n - 1_i32, 1_i32, -1_i32
      solution(j) = modified_rhs(j) - modified_upper(j)*solution(j + 1_i32)
    end do
    success = all(ieee_is_finite(solution))
  end subroutine solve_tridiagonal

  subroutine populate_kinetic_state(options, grid, phi, iterations, residual_norm, state, status)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_grid_type), intent(in) :: grid
    real(dp), intent(in) :: phi(:), residual_norm
    integer(i32), intent(in) :: iterations
    type(outer_plasma_state_type), intent(inout) :: state
    integer(i32), intent(out) :: status
    real(dp) :: electron_density, ion_density, photoelectron_density, susceptibility, speed
    real(dp) :: electron_flux, ion_flux, escape_fraction, return_fraction, barrier, normalization
    integer(i32) :: j, closure_status

    state%profile_n = grid%n
    state%kinetic_closure = options%kinetic_closure
    state%interface_z = 0.0_dp
    state%interface_potential = phi(1)
    state%infinity_potential = 0.0_dp
    state%debye_length = options%tail_length
    state%interface_field = options%interface_field
    state%nonlinear_iterations = iterations
    state%nonlinear_residual = residual_norm
    allocate (state%z(grid%n), state%potential(grid%n), state%field(grid%n), state%charge_density(grid%n))
    state%z = grid%z
    state%potential = phi
    state%field(1) = options%interface_field
    do j = 2_i32, grid%n - 1_i32
      state%field(j) = -(phi(j + 1_i32) - phi(j - 1_i32))/(grid%z(j + 1_i32) - grid%z(j - 1_i32))
    end do
    state%field(grid%n) = -(phi(grid%n) - phi(grid%n - 1_i32))/grid%dz(grid%n - 1_i32)
    do j = 1_i32, grid%n
      call eval_absorbing_maxwellian_density( &
        phi(j), phi(1), options%electron_charge, options%electron_temperature_j, &
        options%electron_density_infinity, electron_density, susceptibility, closure_status &
        )
      if (closure_status /= outer_plasma_ok) then
        status = closure_status
        return
      end if
      call eval_cold_drift_density( &
        phi(j), options%ion_charge, options%ion_mass, options%ion_density_infinity, &
        options%ion_drift_infinity, ion_density, speed, susceptibility, closure_status &
        )
      if (closure_status /= outer_plasma_ok) then
        status = closure_status
        return
      end if
      photoelectron_density = 0.0_dp
      if (options%photoelectron_emission_flux > 0.0_dp) then
        call eval_emitted_maxwellian_density( &
          phi(j), phi(1), 0.0_dp, options%photoelectron_charge, options%photoelectron_mass, &
          options%photoelectron_temperature_j, options%photoelectron_emission_flux, &
          photoelectron_density, closure_status &
          )
        if (closure_status /= outer_plasma_ok) then
          status = closure_status
          return
        end if
      end if
      state%charge_density(j) = options%electron_charge*electron_density + options%ion_charge*ion_density + &
                                options%photoelectron_charge*photoelectron_density
    end do
    state%integrated_charge_per_area = trapz(grid%z, state%charge_density)
    barrier = options%electron_charge*phi(1)/options%electron_temperature_j
    if (barrier >= 0.0_dp) then
      normalization = 1.0_dp + erf(sqrt(barrier))
      electron_flux = 2.0_dp*options%electron_density_infinity* &
                      sqrt(options%electron_temperature_j/(2.0_dp*pi*options%electron_mass))* &
                      exp(-barrier)/normalization
    else
      electron_flux = 2.0_dp*options%electron_density_infinity* &
                      sqrt(options%electron_temperature_j/(2.0_dp*pi*options%electron_mass))
    end if
    ion_flux = options%ion_density_infinity*options%ion_drift_infinity
    escape_fraction = 1.0_dp
    if (options%photoelectron_emission_flux > 0.0_dp) then
      call eval_photoelectron_escape_return( &
        phi(1), 0.0_dp, options%photoelectron_charge, options%photoelectron_temperature_j, &
        escape_fraction, return_fraction, closure_status &
        )
      if (closure_status /= outer_plasma_ok) then
        status = closure_status
        return
      end if
    end if
    call eval_kinetic_current_balance( &
      electron_flux, ion_flux, options%photoelectron_emission_flux, escape_fraction, &
      options%electron_charge, options%ion_charge, options%external_current_density, &
      state%electron_current_density, state%ion_current_density, state%photoelectron_current_density, &
      state%total_current_density &
      )
    status = outer_plasma_ok
  end subroutine populate_kinetic_state

  pure real(dp) function trapz(x, y) result(integral)
    real(dp), intent(in) :: x(:), y(:)
    integer(i32) :: j

    integral = 0.0_dp
    do j = 1_i32, size(x) - 1_i32
      integral = integral + 0.5_dp*(y(j) + y(j + 1_i32))*(x(j + 1_i32) - x(j))
    end do
  end function trapz

  subroutine eval_absorbing_maxwellian_density( &
    phi, phi_interface, charge, temperature_j, density_infinity, density, susceptibility, status, &
    derivative_interface &
    )
    real(dp), intent(in) :: phi, phi_interface, charge, temperature_j, density_infinity
    real(dp), intent(out) :: density, susceptibility
    integer(i32), intent(out) :: status
    real(dp), intent(out), optional :: derivative_interface
    real(dp) :: barrier_interface, barrier_local, denominator, reflected_factor
    real(dp) :: boltzmann_factor, derivative_reflected, derivative_denominator, acceleration

    density = 0.0_dp
    susceptibility = 0.0_dp
    if (present(derivative_interface)) derivative_interface = 0.0_dp
    status = outer_plasma_invalid
    if (.not. all(ieee_is_finite([phi, phi_interface, charge, temperature_j, density_infinity])) .or. &
        charge == 0.0_dp .or. temperature_j <= 0.0_dp .or. density_infinity < 0.0_dp) return

    barrier_interface = charge*phi_interface/temperature_j
    barrier_local = charge*(phi_interface - phi)/temperature_j
    if (barrier_interface < 0.0_dp) then
      acceleration = -charge*phi/temperature_j
      if (acceleration < -64.0_dp*epsilon(1.0_dp)) then
        status = outer_plasma_no_physical_solution
        return
      end if
      acceleration = max(0.0_dp, acceleration)
      density = density_infinity*exp(acceleration)*erfc(sqrt(acceleration))
      if (acceleration > sqrt(epsilon(1.0_dp))) then
        susceptibility = density_infinity*(-charge/temperature_j)*( &
                         exp(acceleration)*erfc(sqrt(acceleration)) - &
                         1.0_dp/sqrt(pi*acceleration) &
                         )
      end if
      status = outer_plasma_ok
      return
    else if (barrier_local < -64.0_dp*epsilon(1.0_dp)) then
      status = outer_plasma_no_physical_solution
      return
    end if
    barrier_local = max(0.0_dp, barrier_local)
    denominator = 1.0_dp + erf(sqrt(barrier_interface))
    reflected_factor = 1.0_dp + erf(sqrt(barrier_local))
    boltzmann_factor = exp(-charge*phi/temperature_j)
    density = density_infinity*boltzmann_factor*reflected_factor/denominator

    derivative_reflected = 0.0_dp
    if (barrier_local > sqrt(epsilon(1.0_dp))) then
      derivative_reflected = (-charge/temperature_j)*exp(-barrier_local)/sqrt(pi*barrier_local)
    end if
    susceptibility = density_infinity*boltzmann_factor/denominator*( &
                     (-charge/temperature_j)*reflected_factor + derivative_reflected &
                     )
    if (present(derivative_interface)) then
      derivative_reflected = 0.0_dp
      if (barrier_local > sqrt(epsilon(1.0_dp))) then
        derivative_reflected = (charge/temperature_j)*exp(-barrier_local)/sqrt(pi*barrier_local)
      end if
      derivative_denominator = 0.0_dp
      if (barrier_interface > sqrt(epsilon(1.0_dp))) then
        derivative_denominator = (charge/temperature_j)*exp(-barrier_interface)/sqrt(pi*barrier_interface)
      end if
      derivative_interface = density_infinity*boltzmann_factor*( &
                             derivative_reflected/denominator - &
                             reflected_factor*derivative_denominator/(denominator*denominator) &
                             )
    end if
    if (.not. all(ieee_is_finite([density, susceptibility]))) then
      density = 0.0_dp
      susceptibility = 0.0_dp
      status = outer_plasma_invalid
      return
    end if
    status = outer_plasma_ok
  end subroutine eval_absorbing_maxwellian_density

  subroutine eval_cold_drift_density( &
    phi, charge, mass, density_infinity, drift_infinity, density, speed, susceptibility, status &
    )
    real(dp), intent(in) :: phi, charge, mass, density_infinity, drift_infinity
    real(dp), intent(out) :: density, speed, susceptibility
    integer(i32), intent(out) :: status
    real(dp) :: speed_squared

    density = 0.0_dp
    speed = 0.0_dp
    susceptibility = 0.0_dp
    status = outer_plasma_invalid
    if (.not. all(ieee_is_finite([phi, charge, mass, density_infinity, drift_infinity])) .or. &
        charge == 0.0_dp .or. mass <= 0.0_dp .or. density_infinity < 0.0_dp .or. drift_infinity <= 0.0_dp) return
    speed_squared = drift_infinity*drift_infinity - 2.0_dp*charge*phi/mass
    if (.not. ieee_is_finite(speed_squared) .or. speed_squared <= 0.0_dp) then
      status = outer_plasma_no_physical_solution
      return
    end if
    speed = sqrt(speed_squared)
    density = density_infinity*drift_infinity/speed
    susceptibility = density*charge/(mass*speed_squared)
    status = outer_plasma_ok
  end subroutine eval_cold_drift_density

  pure real(dp) function kinetic_bohm_speed(electron_temperature_j, ion_temperature_j, ion_mass, ion_gamma) result(speed)
    real(dp), intent(in) :: electron_temperature_j, ion_temperature_j, ion_mass, ion_gamma

    if (electron_temperature_j <= 0.0_dp .or. ion_temperature_j < 0.0_dp .or. &
        ion_mass <= 0.0_dp .or. ion_gamma < 0.0_dp) then
      speed = 0.0_dp
    else
      speed = sqrt((electron_temperature_j + ion_gamma*ion_temperature_j)/ion_mass)
    end if
  end function kinetic_bohm_speed

  subroutine eval_photoelectron_escape_return( &
    phi_interface, phi_infinity, charge, temperature_j, escape_fraction, return_fraction, status &
    )
    real(dp), intent(in) :: phi_interface, phi_infinity, charge, temperature_j
    real(dp), intent(out) :: escape_fraction, return_fraction
    integer(i32), intent(out) :: status
    real(dp) :: barrier

    escape_fraction = 0.0_dp
    return_fraction = 0.0_dp
    status = outer_plasma_invalid
    if (.not. all(ieee_is_finite([phi_interface, phi_infinity, charge, temperature_j])) .or. &
        charge == 0.0_dp .or. temperature_j <= 0.0_dp) return
    barrier = max(0.0_dp, charge*(phi_infinity - phi_interface))
    escape_fraction = exp(-barrier/temperature_j)
    return_fraction = 1.0_dp - escape_fraction
    status = outer_plasma_ok
  end subroutine eval_photoelectron_escape_return

  subroutine eval_emitted_maxwellian_density( &
    phi, phi_interface, phi_infinity, charge, mass, temperature_j, emission_flux, density, status, &
    derivative_local, derivative_interface &
    )
    real(dp), intent(in) :: phi, phi_interface, phi_infinity, charge, mass, temperature_j, emission_flux
    real(dp), intent(out) :: density
    integer(i32), intent(out) :: status
    real(dp), intent(out), optional :: derivative_local, derivative_interface
    real(dp) :: barrier_total, barrier_local, acceleration, prefactor, return_barrier
    real(dp) :: derivative_kernel, returned_factor
    logical :: derivatives_finite

    density = 0.0_dp
    if (present(derivative_local)) derivative_local = 0.0_dp
    if (present(derivative_interface)) derivative_interface = 0.0_dp
    status = outer_plasma_invalid
    if (.not. all(ieee_is_finite([ &
                                 phi, phi_interface, phi_infinity, charge, mass, temperature_j, emission_flux &
                                 ])) .or. charge == 0.0_dp .or. mass <= 0.0_dp .or. &
        temperature_j <= 0.0_dp .or. emission_flux < 0.0_dp) return
    if (emission_flux == 0.0_dp) then
      status = outer_plasma_ok
      return
    end if
    barrier_total = charge*(phi_infinity - phi_interface)/temperature_j
    barrier_local = charge*(phi - phi_interface)/temperature_j
    prefactor = emission_flux*sqrt(pi*mass/(2.0_dp*temperature_j))
    if (barrier_total >= 0.0_dp) then
      if (barrier_local < -64.0_dp*epsilon(1.0_dp) .or. &
          barrier_local > barrier_total + 64.0_dp*epsilon(1.0_dp)) then
        status = outer_plasma_no_physical_solution
        return
      end if
      barrier_local = max(0.0_dp, min(barrier_total, barrier_local))
      return_barrier = max(0.0_dp, barrier_total - barrier_local)
      returned_factor = 1.0_dp + erf(sqrt(return_barrier))
      density = prefactor*exp(-barrier_local)*returned_factor
      if (present(derivative_local)) then
        derivative_kernel = 0.0_dp
        if (return_barrier > sqrt(epsilon(1.0_dp))) then
          derivative_kernel = exp(-return_barrier)/sqrt(pi*return_barrier)
        end if
        derivative_local = prefactor*exp(-barrier_local)*(-charge/temperature_j)*( &
                           returned_factor + derivative_kernel &
                           )
      end if
      if (present(derivative_interface)) derivative_interface = density*charge/temperature_j
    else
      acceleration = -barrier_local
      if (acceleration < -64.0_dp*epsilon(1.0_dp)) then
        status = outer_plasma_no_physical_solution
        return
      end if
      acceleration = max(0.0_dp, acceleration)
      density = prefactor*exp(acceleration)*erfc(sqrt(acceleration))
      derivative_kernel = 0.0_dp
      if (acceleration > sqrt(epsilon(1.0_dp))) then
        derivative_kernel = prefactor*(-charge/temperature_j)*( &
                            exp(acceleration)*erfc(sqrt(acceleration)) - &
                            1.0_dp/sqrt(pi*acceleration) &
                            )
      end if
      if (present(derivative_local)) derivative_local = derivative_kernel
      if (present(derivative_interface)) derivative_interface = -derivative_kernel
    end if
    derivatives_finite = .true.
    if (present(derivative_local)) derivatives_finite = derivatives_finite .and. ieee_is_finite(derivative_local)
    if (present(derivative_interface)) then
      derivatives_finite = derivatives_finite .and. ieee_is_finite(derivative_interface)
    end if
    if (.not. ieee_is_finite(density) .or. .not. derivatives_finite) then
      density = 0.0_dp
      if (present(derivative_local)) derivative_local = 0.0_dp
      if (present(derivative_interface)) derivative_interface = 0.0_dp
      return
    end if
    status = outer_plasma_ok
  end subroutine eval_emitted_maxwellian_density

  pure subroutine eval_kinetic_current_balance( &
    electron_absorption_flux, ion_absorption_flux, photoelectron_emission_flux, &
    photoelectron_escape_fraction, electron_charge, ion_charge, external_current_density, &
    electron_current, ion_current, photoelectron_current, total_current &
    )
    real(dp), intent(in) :: electron_absorption_flux, ion_absorption_flux, photoelectron_emission_flux
    real(dp), intent(in) :: photoelectron_escape_fraction, electron_charge, ion_charge, external_current_density
    real(dp), intent(out) :: electron_current, ion_current, photoelectron_current, total_current

    electron_current = electron_charge*electron_absorption_flux
    ion_current = ion_charge*ion_absorption_flux
    photoelectron_current = -electron_charge*photoelectron_emission_flux*photoelectron_escape_fraction
    total_current = electron_current + ion_current + photoelectron_current + external_current_density
  end subroutine eval_kinetic_current_balance

end module bem_outer_plasma_kinetic
