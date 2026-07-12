module bem_outer_plasma_kinetic
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: pi, eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid, &
                                    outer_plasma_no_physical_solution, outer_plasma_numerical_failure
  use bem_outer_plasma_grid, only: outer_plasma_grid_type, init_outer_plasma_grid
  implicit none
  private

  type, public :: kinetic_outer_plasma_options_type
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
    real(dp) :: residual_tolerance = 1.0e-8_dp
    real(dp) :: external_current_density = 0.0_dp
  end type kinetic_outer_plasma_options_type

  public :: eval_absorbing_maxwellian_density
  public :: eval_cold_drift_density
  public :: kinetic_bohm_speed
  public :: eval_photoelectron_escape_return
  public :: eval_emitted_maxwellian_density
  public :: eval_kinetic_current_balance
  public :: solve_outer_plasma_kinetic

contains

  subroutine solve_outer_plasma_kinetic(options, state, status, message, initial_potential)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: initial_potential(:)
    type(outer_plasma_grid_type) :: grid
    real(dp), allocatable :: phi(:), residual(:), trial_residual(:), jacobian(:, :), delta(:), trial(:)
    real(dp) :: residual_norm, trial_norm, step, bohm_speed, seed_scale
    real(dp) :: previous_interface_field, field_correction
    integer(i32) :: iteration, solve_status
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
    allocate (phi(grid%n), residual(grid%n), trial_residual(grid%n), jacobian(grid%n, grid%n), &
              delta(grid%n), trial(grid%n))
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

    do iteration = 0_i32, options%max_iterations
      if (residual_norm <= options%residual_tolerance) exit
      if (iteration == options%max_iterations) then
        status = outer_plasma_numerical_failure
        state%applicability_status = status
        message = 'kinetic Newton solve reached max_iterations'
        return
      end if
      call numerical_kinetic_jacobian(options, grid, phi, residual, jacobian, solve_status)
      if (solve_status /= outer_plasma_ok) then
        status = outer_plasma_numerical_failure
        state%applicability_status = status
        message = 'kinetic Jacobian evaluation failed'
        return
      end if
      call solve_dense_system(jacobian, -residual, delta, accepted)
      if (.not. accepted) then
        status = outer_plasma_numerical_failure
        state%applicability_status = status
        message = 'kinetic Newton Jacobian is singular'
        return
      end if
      step = 1.0_dp
      accepted = .false.
      do while (step >= 2.0_dp**(-20))
        trial = phi + step*delta
        if (monotonic_electron_repelling_branch(options, trial)) then
          call kinetic_residual(options, grid, trial, trial_residual, solve_status)
          if (solve_status == outer_plasma_ok) then
            trial_norm = maxval(abs(trial_residual))
            if (trial_norm < residual_norm) then
              accepted = .true.
              exit
            end if
          end if
        end if
        step = 0.5_dp*step
      end do
      if (.not. accepted) then
        status = outer_plasma_numerical_failure
        state%applicability_status = status
        message = 'kinetic Newton line search failed on the monotonic branch'
        return
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
  end subroutine solve_outer_plasma_kinetic

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

  subroutine numerical_kinetic_jacobian(options, grid, phi, residual, jacobian, status)
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(outer_plasma_grid_type), intent(in) :: grid
    real(dp), intent(in) :: phi(:), residual(:)
    real(dp), intent(out) :: jacobian(:, :)
    integer(i32), intent(out) :: status
    real(dp) :: perturbed(size(phi)), shifted_residual(size(phi)), increment
    integer(i32) :: column, closure_status

    status = outer_plasma_ok
    do column = 1_i32, size(phi)
      increment = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(phi(column)))
      perturbed = phi
      if (column == 1_i32 .and. abs(phi(column)) < increment) increment = -increment
      perturbed(column) = perturbed(column) + increment
      call kinetic_residual(options, grid, perturbed, shifted_residual, closure_status)
      if (closure_status /= outer_plasma_ok) then
        increment = -increment
        perturbed = phi
        perturbed(column) = perturbed(column) + increment
        call kinetic_residual(options, grid, perturbed, shifted_residual, closure_status)
      end if
      if (closure_status /= outer_plasma_ok) then
        status = closure_status
        return
      end if
      jacobian(:, column) = (shifted_residual - residual)/increment
    end do
  end subroutine numerical_kinetic_jacobian

  subroutine solve_dense_system(matrix, rhs, solution, success)
    real(dp), intent(in) :: matrix(:, :), rhs(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: success
    real(dp) :: work(size(rhs), size(rhs)), b(size(rhs)), factor, row_buffer(size(rhs)), pivot_value
    integer(i32) :: n, column, row, pivot

    n = size(rhs)
    work = matrix
    b = rhs
    success = .false.
    do column = 1_i32, n
      pivot = column - 1_i32 + maxloc(abs(work(column:n, column)), dim=1)
      pivot_value = work(pivot, column)
      if (.not. ieee_is_finite(pivot_value) .or. abs(pivot_value) <= tiny(1.0_dp)) return
      if (pivot /= column) then
        row_buffer = work(column, :)
        work(column, :) = work(pivot, :)
        work(pivot, :) = row_buffer
        factor = b(column)
        b(column) = b(pivot)
        b(pivot) = factor
      end if
      do row = column + 1_i32, n
        factor = work(row, column)/work(column, column)
        work(row, column:n) = work(row, column:n) - factor*work(column, column:n)
        b(row) = b(row) - factor*b(column)
      end do
    end do
    do row = n, 1_i32, -1_i32
      solution(row) = (b(row) - dot_product(work(row, row + 1_i32:n), solution(row + 1_i32:n)))/work(row, row)
    end do
    success = all(ieee_is_finite(solution))
  end subroutine solve_dense_system

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
    phi, phi_interface, charge, temperature_j, density_infinity, density, susceptibility, status &
    )
    real(dp), intent(in) :: phi, phi_interface, charge, temperature_j, density_infinity
    real(dp), intent(out) :: density, susceptibility
    integer(i32), intent(out) :: status
    real(dp) :: barrier_interface, barrier_local, denominator, reflected_factor
    real(dp) :: boltzmann_factor, derivative_reflected, acceleration

    density = 0.0_dp
    susceptibility = 0.0_dp
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
    phi, phi_interface, phi_infinity, charge, mass, temperature_j, emission_flux, density, status &
    )
    real(dp), intent(in) :: phi, phi_interface, phi_infinity, charge, mass, temperature_j, emission_flux
    real(dp), intent(out) :: density
    integer(i32), intent(out) :: status
    real(dp) :: barrier_total, barrier_local, acceleration, prefactor, return_barrier

    density = 0.0_dp
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
      density = prefactor*exp(-barrier_local)*(1.0_dp + erf(sqrt(return_barrier)))
    else
      acceleration = -barrier_local
      if (acceleration < -64.0_dp*epsilon(1.0_dp)) then
        status = outer_plasma_no_physical_solution
        return
      end if
      acceleration = max(0.0_dp, acceleration)
      density = prefactor*exp(acceleration)*erfc(sqrt(acceleration))
    end if
    if (.not. ieee_is_finite(density)) then
      density = 0.0_dp
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
