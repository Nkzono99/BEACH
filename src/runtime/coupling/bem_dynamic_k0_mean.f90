!> tracked粒子の局所depositと独立に、平面平均k=0電荷を陰的に進めるruntime adapter。
module bem_dynamic_k0_mean
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use bem_kinds, only: dp, i32
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_invalid
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type, &
                                      eval_photoelectron_escape_return
  use bem_mean_charge_integrator, only: &
    mean_charge_step_ok, mean_charge_step_result_type, &
    planar_debye_capacitance_per_area, advance_mean_charge_backward_euler
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: dynamic_k0_ok = 0_i32
  integer(i32), parameter, public :: dynamic_k0_invalid = 1_i32
  integer(i32), parameter, public :: dynamic_k0_numerical_failure = 2_i32

  integer(i32), parameter :: current_parameter_count = 6_i32
  integer(i32), parameter :: parameter_electron_current = 1_i32
  integer(i32), parameter :: parameter_ion_current = 2_i32
  integer(i32), parameter :: parameter_photoelectron_outward_current = 3_i32
  integer(i32), parameter :: parameter_photoelectron_charge = 4_i32
  integer(i32), parameter :: parameter_photoelectron_temperature = 5_i32
  integer(i32), parameter :: parameter_other_current = 6_i32

  !> 1 batchの陰的平均更新結果。正のcurrentはsurfaceを正に帯電させる。
  type, public :: dynamic_k0_step_type
    integer(i32) :: status = dynamic_k0_invalid
    integer(i32) :: iterations = 0_i32
    integer(i32) :: bracket_expansions = 0_i32
    real(dp) :: surface_charge_before_c = 0.0_dp
    real(dp) :: surface_charge_after_c = 0.0_dp
    real(dp) :: interface_potential_before_v = 0.0_dp
    real(dp) :: interface_potential_after_v = 0.0_dp
    real(dp) :: interface_field_after_v_m = 0.0_dp
    real(dp) :: electron_current_density_a_m2 = 0.0_dp
    real(dp) :: ion_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_current_density_a_m2 = 0.0_dp
    real(dp) :: additional_tracked_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_escape_fraction = 1.0_dp
    real(dp) :: photoelectron_return_fraction = 0.0_dp
    real(dp) :: total_current_density_a_m2 = 0.0_dp
    real(dp) :: backward_euler_residual_v = 0.0_dp
    real(dp) :: explicit_predictor_potential_v = 0.0_dp
    real(dp) :: returned_outer_flight_time_mean_s = 0.0_dp
    real(dp) :: estimated_returning_photoelectron_column_charge_per_area_c_m2 = 0.0_dp
    real(dp) :: backward_euler_residual_charge_c = 0.0_dp
    real(dp) :: photoelectron_source_charge_c = 0.0_dp
    real(dp) :: photoelectron_barrier_energy_j = 0.0_dp
    real(dp) :: photoelectron_energy_max_j = 0.0_dp
    real(dp) :: marginal_photoelectron_energy_j = -1.0_dp
    real(dp) :: marginal_photoelectron_escape_fraction = -1.0_dp
    real(dp) :: zhao_effective_source_scale = 0.0_dp
    real(dp) :: photoelectron_recross_charge_fraction = 0.0_dp
    real(dp) :: photoelectron_terminal_mismatch_charge_fraction = 0.0_dp
    character(len=1) :: zhao_branch = ' '
  end type dynamic_k0_step_type

  public :: advance_dynamic_k0_mean

contains

  !> ambient_linear_debyeの容量とspecies電流を用いてmean surface chargeを1 batch進める。
  !!
  !! e_bottom_zeroでは sigma=eps0*E、symmetric_vacuumではsigma=2*eps0*E。
  !! profileは phi=lambda_D*E なので、後者の単位面積容量は前者の2倍になる。
  subroutine advance_dynamic_k0_mean( &
    kinetic_options, lower_boundary_model, area_xy_m2, surface_charge_before_c, time_step_s, step, message, &
    tracked_electron_current_density_a_m2, tracked_ion_current_density_a_m2, &
    tracked_photoelectron_outward_current_density_a_m2, tracked_additional_current_density_a_m2 &
    )
    type(kinetic_outer_plasma_options_type), intent(in) :: kinetic_options
    character(len=*), intent(in) :: lower_boundary_model
    real(dp), intent(in) :: area_xy_m2, surface_charge_before_c, time_step_s
    type(dynamic_k0_step_type), intent(out) :: step
    character(len=*), intent(out) :: message
    real(dp), intent(in) :: tracked_electron_current_density_a_m2
    real(dp), intent(in) :: tracked_ion_current_density_a_m2
    real(dp), intent(in) :: tracked_photoelectron_outward_current_density_a_m2
    real(dp), intent(in) :: tracked_additional_current_density_a_m2

    type(mean_charge_step_result_type) :: integrated
    real(dp) :: parameters(current_parameter_count)
    real(dp) :: capacitance_per_area, boundary_factor, surface_charge_density
    integer(i32) :: current_status

    step = dynamic_k0_step_type()
    message = ''
    if (trim(lower_ascii(kinetic_options%kinetic_closure)) /= 'ambient_linear_debye') then
      message = 'dynamic k0 mean currently requires kinetic_closure=ambient_linear_debye'
      return
    end if
    select case (trim(lower_ascii(lower_boundary_model)))
    case ('e_bottom_zero')
      boundary_factor = 1.0_dp
    case ('symmetric_vacuum')
      boundary_factor = 2.0_dp
    case default
      message = 'dynamic k0 mean received an unsupported lower boundary model'
      return
    end select
    if (.not. all(ieee_is_finite([ &
                                 area_xy_m2, surface_charge_before_c, time_step_s, kinetic_options%tail_length, &
                                 tracked_electron_current_density_a_m2, tracked_ion_current_density_a_m2, &
                                 tracked_photoelectron_outward_current_density_a_m2, &
                                 tracked_additional_current_density_a_m2, kinetic_options%external_current_density &
                                 ])) .or. area_xy_m2 <= 0.0_dp .or. time_step_s <= 0.0_dp .or. &
        kinetic_options%tail_length <= 0.0_dp .or. &
        tracked_photoelectron_outward_current_density_a_m2 < 0.0_dp) then
      message = 'dynamic k0 mean received invalid geometry, time step, or tracked current'
      return
    end if
    if (tracked_photoelectron_outward_current_density_a_m2 > 0.0_dp) then
      if (.not. ieee_is_finite(kinetic_options%photoelectron_charge) .or. &
          kinetic_options%photoelectron_charge >= 0.0_dp .or. &
          .not. ieee_is_finite(kinetic_options%photoelectron_temperature_j) .or. &
          kinetic_options%photoelectron_temperature_j <= 0.0_dp) then
        message = 'dynamic k0 photoelectron barrier requires negative charge and positive temperature'
        return
      end if
    end if

    capacitance_per_area = boundary_factor* &
                           planar_debye_capacitance_per_area(kinetic_options%tail_length)
    surface_charge_density = surface_charge_before_c/area_xy_m2
    parameters(parameter_electron_current) = tracked_electron_current_density_a_m2
    parameters(parameter_ion_current) = tracked_ion_current_density_a_m2
    parameters(parameter_photoelectron_outward_current) = &
      tracked_photoelectron_outward_current_density_a_m2
    parameters(parameter_photoelectron_charge) = kinetic_options%photoelectron_charge
    parameters(parameter_photoelectron_temperature) = kinetic_options%photoelectron_temperature_j
    parameters(parameter_other_current) = &
      kinetic_options%external_current_density + tracked_additional_current_density_a_m2
    if (.not. ieee_is_finite(parameters(parameter_other_current))) then
      message = 'dynamic k0 combined external current is not finite'
      return
    end if
    call advance_mean_charge_backward_euler( &
      surface_charge_density, time_step_s, capacitance_per_area, tracked_current_model, parameters, integrated &
      )
    if (integrated%status /= mean_charge_step_ok) then
      step%status = dynamic_k0_numerical_failure
      write (message, '(a,i0)') 'dynamic k0 backward-Euler solve failed with status ', integrated%status
      return
    end if

    step%surface_charge_before_c = surface_charge_before_c
    step%surface_charge_after_c = integrated%surface_charge_density_c_m2*area_xy_m2
    step%interface_potential_before_v = surface_charge_density/capacitance_per_area
    step%interface_potential_after_v = integrated%potential_v
    step%interface_field_after_v_m = integrated%potential_v/kinetic_options%tail_length
    step%total_current_density_a_m2 = integrated%current_density_a_m2
    step%backward_euler_residual_v = integrated%backward_euler_residual_v
    step%explicit_predictor_potential_v = integrated%explicit_predictor_potential_v
    step%additional_tracked_current_density_a_m2 = tracked_additional_current_density_a_m2
    step%iterations = integrated%iterations
    step%bracket_expansions = integrated%bracket_expansions
    call eval_dynamic_currents( &
      step%interface_potential_after_v, &
      parameters(parameter_electron_current), parameters(parameter_ion_current), &
      parameters(parameter_photoelectron_outward_current), &
      parameters(parameter_photoelectron_charge), parameters(parameter_photoelectron_temperature), &
      parameters(parameter_other_current), step%electron_current_density_a_m2, &
      step%ion_current_density_a_m2, step%photoelectron_current_density_a_m2, &
      step%total_current_density_a_m2, step%photoelectron_escape_fraction, &
      step%photoelectron_return_fraction, current_status &
      )
    if (current_status /= outer_plasma_ok .or. &
        .not. all(ieee_is_finite([ &
                                 step%surface_charge_after_c, step%interface_potential_after_v, &
                                 step%interface_field_after_v_m, step%electron_current_density_a_m2, &
                                 step%ion_current_density_a_m2, step%photoelectron_current_density_a_m2, &
                                 step%photoelectron_escape_fraction, step%photoelectron_return_fraction, &
                                 step%total_current_density_a_m2 &
                                 ]))) then
      step%status = dynamic_k0_numerical_failure
      message = 'dynamic k0 current diagnostics are invalid'
      return
    end if
    step%status = dynamic_k0_ok
    message = ''
  end subroutine advance_dynamic_k0_mean

  pure function tracked_current_model(interface_potential_v, parameters) result(current_density_a_m2)
    real(dp), intent(in) :: interface_potential_v
    real(dp), intent(in) :: parameters(:)
    real(dp) :: current_density_a_m2
    real(dp) :: electron_current, ion_current, photoelectron_current
    real(dp) :: escape_fraction, return_fraction
    integer(i32) :: status

    if (size(parameters) /= current_parameter_count) then
      current_density_a_m2 = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    call eval_dynamic_currents( &
      interface_potential_v, &
      parameters(parameter_electron_current), parameters(parameter_ion_current), &
      parameters(parameter_photoelectron_outward_current), &
      parameters(parameter_photoelectron_charge), parameters(parameter_photoelectron_temperature), &
      parameters(parameter_other_current), electron_current, ion_current, photoelectron_current, &
      current_density_a_m2, escape_fraction, return_fraction, status &
      )
    if (status /= outer_plasma_ok) current_density_a_m2 = ieee_value(0.0_dp, ieee_quiet_nan)
  end function tracked_current_model

  pure subroutine eval_dynamic_currents( &
    interface_potential_v, tracked_electron_current, tracked_ion_current, &
    tracked_photoelectron_outward_current, photoelectron_charge, photoelectron_temperature, &
    other_current, electron_current, ion_current, photoelectron_current, total_current, &
    escape_fraction, return_fraction, status &
    )
    real(dp), intent(in) :: interface_potential_v
    real(dp), intent(in) :: tracked_electron_current, tracked_ion_current
    real(dp), intent(in) :: tracked_photoelectron_outward_current
    real(dp), intent(in) :: photoelectron_charge, photoelectron_temperature, other_current
    real(dp), intent(out) :: electron_current, ion_current, photoelectron_current, total_current
    real(dp), intent(out) :: escape_fraction, return_fraction
    integer(i32), intent(out) :: status

    electron_current = tracked_electron_current
    ion_current = tracked_ion_current
    photoelectron_current = 0.0_dp
    total_current = 0.0_dp
    escape_fraction = 1.0_dp
    return_fraction = 0.0_dp
    status = outer_plasma_invalid
    if (.not. all(ieee_is_finite([ &
                                 interface_potential_v, electron_current, ion_current, &
                                 tracked_photoelectron_outward_current, photoelectron_charge, &
                                 photoelectron_temperature, other_current &
                                 ])) .or. tracked_photoelectron_outward_current < 0.0_dp) return
    if (tracked_photoelectron_outward_current > 0.0_dp) then
      call eval_photoelectron_escape_return( &
        interface_potential_v, 0.0_dp, photoelectron_charge, photoelectron_temperature, &
        escape_fraction, return_fraction, status &
        )
      if (status /= outer_plasma_ok) return
    end if
    photoelectron_current = tracked_photoelectron_outward_current*escape_fraction
    total_current = electron_current + ion_current + photoelectron_current + other_current
    if (.not. all(ieee_is_finite([electron_current, ion_current, photoelectron_current, total_current]))) then
      electron_current = 0.0_dp
      ion_current = 0.0_dp
      photoelectron_current = 0.0_dp
      total_current = 0.0_dp
      status = outer_plasma_invalid
      return
    end if
    status = outer_plasma_ok
  end subroutine eval_dynamic_currents

end module bem_dynamic_k0_mean
