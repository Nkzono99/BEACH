!> 外部モデルから species 別の固定表面電流を解決する。
module bem_surface_current_model
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann, pi, qe
  use bem_app_config_types, only: app_config
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_config_helpers, only: species_number_density_m3, species_temperature_k
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params, try_solve_zhao_unknowns, &
                                   swe_free_current_term
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  type, extends(surface_closure_contract_type), public :: surface_current_model_result_type
    character(len=32) :: model = 'none'
    character(len=1) :: zhao_branch = ' '
    integer(i32) :: electron_species_idx = 0_i32
    integer(i32) :: ion_species_idx = 0_i32
    integer(i32) :: photoelectron_species_idx = 0_i32
    logical :: photoelectron_active = .false.
    real(dp) :: reference_area_m2 = 0.0_dp
    real(dp) :: phi0_v = 0.0_dp
    real(dp) :: phi_m_v = 0.0_dp
    real(dp) :: ambient_electron_density_m3 = 0.0_dp
    real(dp) :: electron_current_density_a_m2 = 0.0_dp
    real(dp) :: ion_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_emission_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_escape_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_return_current_density_a_m2 = 0.0_dp
    real(dp) :: net_current_density_a_m2 = 0.0_dp
    real(dp) :: photoelectron_budget_residual_current_density_a_m2 = 0.0_dp
    real(dp) :: surface_budget_residual_current_density_a_m2 = 0.0_dp
    character(len=32) :: kinetic_contract = 'none'
  end type surface_current_model_result_type

  public :: evaluate_surface_current_model
  public :: evaluate_surface_closure

contains

  !> モデル固有の診断値をシミュレータへ漏らさず、境界契約だけを返す。
  subroutine evaluate_surface_closure(app, contract)
    type(app_config), intent(in) :: app
    type(surface_closure_contract_type), intent(out) :: contract
    type(surface_current_model_result_type) :: detailed_result

    call evaluate_surface_current_model(app, detailed_result)
    contract = detailed_result%surface_closure_contract_type
  end subroutine evaluate_surface_closure

  !> 設定されたmodelをdispatchし、固定電流closure用のtarget配列を返す。
  subroutine evaluate_surface_current_model(app, result)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(out) :: result

    allocate ( &
      result%has_absorbed_target(app%n_particle_species), &
      result%has_emission_target(app%n_particle_species), &
      result%has_escape_target(app%n_particle_species), &
      result%has_inflow_kinetic_map(app%n_particle_species), &
      result%has_outflow_kinetic_barrier(app%n_particle_species), &
      result%absorbed_current_a(app%n_particle_species), &
      result%emission_current_a(app%n_particle_species), &
      result%escaped_particle_current_a(app%n_particle_species), &
      result%inflow_reservoir_potential_v(app%n_particle_species), &
      result%inflow_access_potential_v(app%n_particle_species), &
      result%inflow_kinetic_face(app%n_particle_species), &
      result%outflow_barrier_potential_v(app%n_particle_species), &
      result%outflow_barrier_face(app%n_particle_species) &
      )
    result%has_absorbed_target = .false.
    result%has_emission_target = .false.
    result%has_escape_target = .false.
    result%has_inflow_kinetic_map = .false.
    result%has_outflow_kinetic_barrier = .false.
    result%absorbed_current_a = 0.0_dp
    result%emission_current_a = 0.0_dp
    result%escaped_particle_current_a = 0.0_dp
    result%inflow_reservoir_potential_v = 0.0_dp
    result%inflow_access_potential_v = 0.0_dp
    result%inflow_kinetic_face = 0_i32
    result%outflow_barrier_potential_v = 0.0_dp
    result%outflow_barrier_face = 0_i32
    result%model = trim(lower_ascii(app%surface_current%model))

    select case (trim(result%model))
    case ('none')
      return
    case ('zhao_stationary')
      call evaluate_zhao_stationary_current(app, result)
    case default
      error stop 'Unknown surface current model dispatch.'
    end select
  end subroutine evaluate_surface_current_model

  subroutine evaluate_zhao_stationary_current(app, result)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(inout) :: result
    type(zhao_params_type) :: params
    integer(i32) :: electron_idx, ion_idx, photo_idx
    real(dp) :: electron_temperature_ev, photo_temperature_ev, solar_elevation_deg, photoelectron_ref_density_m3
    real(dp) :: electron_drift_mps, ion_drift_mps, a_swe, electron_term, ion_term, photo_escape_term
    real(dp) :: scale, area, budget_scale, budget_tolerance, electron_bottleneck_potential_v
    character(len=16) :: solver_name
    logical :: success, photoelectron_active

    electron_idx = species_index(app, app%surface_current%electron_species)
    ion_idx = species_index(app, app%surface_current%ion_species)
    photoelectron_active = app%surface_current%photoelectron_source_scale > 0.0_dp
    photo_idx = 0_i32
    if (photoelectron_active) photo_idx = species_index(app, app%surface_current%photoelectron_species)
    electron_temperature_ev = species_temperature_k(app%particle_species(electron_idx))*k_boltzmann/qe
    if (photoelectron_active) then
      photo_temperature_ev = species_temperature_k(app%particle_species(photo_idx))*k_boltzmann/qe
      solar_elevation_deg = app%surface_current%solar_elevation_deg
      photoelectron_ref_density_m3 = app%surface_current%photoelectron_ref_density_m3
    else
      ! source_scale=0ではPE正規化量は解に寄与しない。正のambient量でcore契約だけを満たす。
      photo_temperature_ev = electron_temperature_ev
      solar_elevation_deg = 90.0_dp
      photoelectron_ref_density_m3 = species_number_density_m3(app%particle_species(ion_idx))
    end if
    electron_drift_mps = -app%particle_species(electron_idx)%drift_velocity(3)
    ion_drift_mps = -app%particle_species(ion_idx)%drift_velocity(3)
    area = (app%sim%box_max(1) - app%sim%box_min(1))*(app%sim%box_max(2) - app%sim%box_min(2))
    if (app%surface_current%has_reference_area_m2) area = app%surface_current%reference_area_m2
    if (.not. ieee_is_finite(area) .or. area <= 0.0_dp) then
      error stop 'Zhao stationary surface-current reference area must be finite and positive.'
    end if

    call build_zhao_params( &
      solar_elevation_deg, &
      species_number_density_m3(app%particle_species(ion_idx)), &
      photoelectron_ref_density_m3, &
      electron_temperature_ev, photo_temperature_ev, electron_drift_mps, ion_drift_mps, &
      app%particle_species(ion_idx)%m_particle, app%particle_species(electron_idx)%m_particle, params, &
      photoelectron_source_scale=app%surface_current%photoelectron_source_scale &
      )
    if (.not. photoelectron_active) then
      select case (trim(lower_ascii(app%surface_current%zhao_branch)))
      case ('auto', 'c')
        solver_name = 'zhao_c'
      case default
        error stop 'photoelectron_source_scale=0 requires surface_current_model.zhao_branch="auto" or "c".'
      end select
    else
      solver_name = 'zhao_'//trim(lower_ascii(app%surface_current%zhao_branch))
    end if
    call try_solve_zhao_unknowns( &
      trim(solver_name), params, result%phi0_v, result%phi_m_v, result%ambient_electron_density_m3, &
      result%zhao_branch, success &
      )
    if (.not. success) error stop 'Zhao stationary surface-current root solve failed.'

    ion_term = params%n_swi_inf_m3*sqrt( &
               2.0_dp*pi*params%t_swe_ev/params%t_phe_ev*params%m_e_kg/params%m_i_kg &
               )*params%mach
    select case (result%zhao_branch)
    case ('A')
      a_swe = sqrt(max(0.0_dp, -result%phi_m_v/params%t_swe_ev)) - params%u
      electron_term = swe_free_current_term(params, result%ambient_electron_density_m3, a_swe)
      photo_escape_term = params%n_phe0_m3*exp((result%phi_m_v - result%phi0_v)/params%t_phe_ev)
    case ('B')
      electron_term = swe_free_current_term(params, result%ambient_electron_density_m3, -params%u)
      photo_escape_term = params%n_phe0_m3*exp(-result%phi0_v/params%t_phe_ev)
    case ('C')
      a_swe = sqrt(max(0.0_dp, -result%phi0_v/params%t_swe_ev)) - params%u
      electron_term = swe_free_current_term(params, result%ambient_electron_density_m3, a_swe)
      photo_escape_term = params%n_phe0_m3
    case default
      error stop 'Zhao stationary current returned an unknown branch.'
    end select
    scale = qe*params%v_phe_th_mps/(2.0_dp*sqrt(pi))
    result%electron_current_density_a_m2 = -scale*electron_term
    result%ion_current_density_a_m2 = scale*ion_term
    result%photoelectron_emission_current_density_a_m2 = scale*params%n_phe0_m3
    result%photoelectron_escape_current_density_a_m2 = scale*photo_escape_term
    result%photoelectron_return_current_density_a_m2 = &
      result%photoelectron_escape_current_density_a_m2 - result%photoelectron_emission_current_density_a_m2
    result%net_current_density_a_m2 = result%electron_current_density_a_m2 + &
                                      result%ion_current_density_a_m2 + &
                                      result%photoelectron_escape_current_density_a_m2
    result%photoelectron_budget_residual_current_density_a_m2 = &
      result%photoelectron_emission_current_density_a_m2 + &
      result%photoelectron_return_current_density_a_m2 - &
      result%photoelectron_escape_current_density_a_m2
    result%surface_budget_residual_current_density_a_m2 = &
      result%electron_current_density_a_m2 + result%ion_current_density_a_m2 + &
      result%photoelectron_emission_current_density_a_m2 + &
      result%photoelectron_return_current_density_a_m2
    if (.not. all(ieee_is_finite([ &
                                 result%electron_current_density_a_m2, result%ion_current_density_a_m2, &
                                 result%photoelectron_emission_current_density_a_m2, &
                                 result%photoelectron_escape_current_density_a_m2, &
                                 result%photoelectron_return_current_density_a_m2, result%net_current_density_a_m2, &
                                 result%photoelectron_budget_residual_current_density_a_m2, &
                                 result%surface_budget_residual_current_density_a_m2 &
                                 ]))) then
      error stop 'Zhao stationary surface-current evaluation produced non-finite currents.'
    end if
    if (result%electron_current_density_a_m2 >= 0.0_dp .or. &
        result%ion_current_density_a_m2 <= 0.0_dp .or. &
        result%photoelectron_return_current_density_a_m2 > 0.0_dp .or. &
        result%photoelectron_escape_current_density_a_m2 < 0.0_dp) then
      error stop 'Zhao stationary surface-current evaluation produced invalid channel signs.'
    end if
    if (photoelectron_active) then
      if (result%photoelectron_emission_current_density_a_m2 <= 0.0_dp) then
        error stop 'Zhao stationary photoelectron closure requires a positive emission current.'
      end if
    else if (any([ &
                 result%photoelectron_emission_current_density_a_m2, &
                 result%photoelectron_escape_current_density_a_m2, &
                 result%photoelectron_return_current_density_a_m2 &
                 ] /= 0.0_dp)) then
      error stop 'Zhao stationary no-photoelectron closure produced a nonzero photoelectron current.'
    end if
    budget_scale = max( &
                   abs(result%electron_current_density_a_m2), abs(result%ion_current_density_a_m2), &
                   abs(result%photoelectron_emission_current_density_a_m2), &
                   abs(result%photoelectron_escape_current_density_a_m2), &
                   abs(result%photoelectron_return_current_density_a_m2), tiny(1.0_dp) &
                   )
    budget_tolerance = sqrt(epsilon(1.0_dp))*budget_scale
    if (abs(result%photoelectron_budget_residual_current_density_a_m2) > budget_tolerance) then
      error stop 'Zhao stationary PE current budget does not close.'
    end if
    if (abs(result%surface_budget_residual_current_density_a_m2) > budget_tolerance) then
      error stop 'Zhao stationary surface current budget does not close.'
    end if

    result%active = .true.
    result%reference_area_m2 = area
    result%electron_species_idx = electron_idx
    result%ion_species_idx = ion_idx
    result%photoelectron_species_idx = photo_idx
    result%photoelectron_active = photoelectron_active
    result%has_absorbed_target([electron_idx, ion_idx]) = .true.
    result%absorbed_current_a(electron_idx) = &
      checked_area_current(area, result%electron_current_density_a_m2)
    result%absorbed_current_a(ion_idx) = checked_area_current(area, result%ion_current_density_a_m2)
    if (photoelectron_active) then
      result%has_absorbed_target(photo_idx) = .true.
      result%has_emission_target(photo_idx) = .true.
      result%has_escape_target(photo_idx) = .true.
      result%absorbed_current_a(photo_idx) = &
        checked_area_current(area, result%photoelectron_return_current_density_a_m2)
      result%emission_current_a(photo_idx) = &
        checked_area_current(area, result%photoelectron_emission_current_density_a_m2)
      ! escaped_to_infinity は粒子電荷の外向きfluxなので、正の表面帯電電流とは符号が逆。
      result%escaped_particle_current_a(photo_idx) = &
        -checked_area_current(area, result%photoelectron_escape_current_density_a_m2)
    end if

    ! Zhao の1-D外部シースを、z-high interfaceに対するkinetic boundary mapへ縮約する。
    ! Type Aの電子はphi_mがaccess bottleneckであり、Type B/Cはphi_infinity=0を使う。
    electron_bottleneck_potential_v = 0.0_dp
    if (result%zhao_branch == 'A') electron_bottleneck_potential_v = result%phi_m_v
    result%kinetic_contract = 'zhao_barrier_v1'
    result%has_inflow_kinetic_map([electron_idx, ion_idx]) = .true.
    result%inflow_reservoir_potential_v([electron_idx, ion_idx]) = 0.0_dp
    result%inflow_access_potential_v(electron_idx) = electron_bottleneck_potential_v
    result%inflow_access_potential_v(ion_idx) = 0.0_dp
    result%inflow_kinetic_face([electron_idx, ion_idx]) = 6_i32
    result%has_outflow_kinetic_barrier([electron_idx, ion_idx]) = .true.
    result%outflow_barrier_potential_v(electron_idx) = electron_bottleneck_potential_v
    result%outflow_barrier_potential_v(ion_idx) = 0.0_dp
    result%outflow_barrier_face([electron_idx, ion_idx]) = 6_i32
    if (photoelectron_active) then
      result%has_outflow_kinetic_barrier(photo_idx) = .true.
      result%outflow_barrier_potential_v(photo_idx) = electron_bottleneck_potential_v
      result%outflow_barrier_face(photo_idx) = 6_i32
    end if
  end subroutine evaluate_zhao_stationary_current

  integer(i32) function species_index(app, species_key) result(index_value)
    type(app_config), intent(in) :: app
    character(len=*), intent(in) :: species_key
    integer(i32) :: idx

    index_value = 0_i32
    do idx = 1_i32, app%n_particle_species
      if (.not. app%particle_species(idx)%enabled) cycle
      if (trim(app%particle_species(idx)%species_key) /= trim(species_key)) cycle
      index_value = idx
      return
    end do
    error stop 'Surface current model species resolution failed: '//trim(species_key)
  end function species_index

  real(dp) function checked_area_current(area_m2, current_density_a_m2) result(current_a)
    real(dp), intent(in) :: area_m2, current_density_a_m2

    if (.not. all(ieee_is_finite([area_m2, current_density_a_m2])) .or. area_m2 <= 0.0_dp) then
      error stop 'Zhao stationary surface-current target conversion received invalid input.'
    end if
    if (area_m2 > 1.0_dp .and. abs(current_density_a_m2) > huge(current_a)/area_m2) then
      error stop 'Zhao stationary surface-current target conversion overflowed.'
    end if
    current_a = area_m2*current_density_a_m2
    if (.not. ieee_is_finite(current_a)) then
      error stop 'Zhao stationary surface-current target conversion produced a non-finite current.'
    end if
    if (current_density_a_m2 /= 0.0_dp .and. current_a == 0.0_dp) then
      error stop 'Zhao stationary surface-current target conversion underflowed.'
    end if
  end function checked_area_current

end module bem_surface_current_model
