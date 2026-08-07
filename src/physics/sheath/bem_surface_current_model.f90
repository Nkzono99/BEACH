!> 外部モデルから species 別の固定表面電流を解決する。
module bem_surface_current_model
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann, pi, qe
  use bem_app_config_types, only: app_config
  use bem_config_helpers, only: species_number_density_m3, species_temperature_k
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params, try_solve_zhao_unknowns, &
                                   swe_free_current_term
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  type, public :: surface_current_model_result_type
    logical :: active = .false.
    character(len=32) :: model = 'none'
    character(len=1) :: zhao_branch = ' '
    integer(i32) :: electron_species_idx = 0_i32
    integer(i32) :: ion_species_idx = 0_i32
    integer(i32) :: photoelectron_species_idx = 0_i32
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
    logical, allocatable :: has_absorbed_target(:)
    logical, allocatable :: has_emission_target(:)
    logical, allocatable :: has_escape_target(:)
    real(dp), allocatable :: absorbed_current_a(:)
    real(dp), allocatable :: emission_current_a(:)
    real(dp), allocatable :: escaped_particle_current_a(:)
  end type surface_current_model_result_type

  public :: evaluate_surface_current_model

contains

  !> 設定されたmodelをdispatchし、固定電流closure用のtarget配列を返す。
  subroutine evaluate_surface_current_model(app, result)
    type(app_config), intent(in) :: app
    type(surface_current_model_result_type), intent(out) :: result

    allocate ( &
      result%has_absorbed_target(app%n_particle_species), &
      result%has_emission_target(app%n_particle_species), &
      result%has_escape_target(app%n_particle_species), &
      result%absorbed_current_a(app%n_particle_species), &
      result%emission_current_a(app%n_particle_species), &
      result%escaped_particle_current_a(app%n_particle_species) &
      )
    result%has_absorbed_target = .false.
    result%has_emission_target = .false.
    result%has_escape_target = .false.
    result%absorbed_current_a = 0.0_dp
    result%emission_current_a = 0.0_dp
    result%escaped_particle_current_a = 0.0_dp
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
    real(dp) :: electron_temperature_ev, photo_temperature_ev
    real(dp) :: electron_drift_mps, ion_drift_mps, a_swe, electron_term, ion_term, photo_escape_term
    real(dp) :: scale, area, budget_scale, budget_tolerance
    character(len=16) :: solver_name
    logical :: success

    electron_idx = species_index(app, app%surface_current%electron_species)
    ion_idx = species_index(app, app%surface_current%ion_species)
    photo_idx = species_index(app, app%surface_current%photoelectron_species)
    electron_temperature_ev = species_temperature_k(app%particle_species(electron_idx))*k_boltzmann/qe
    photo_temperature_ev = species_temperature_k(app%particle_species(photo_idx))*k_boltzmann/qe
    electron_drift_mps = -app%particle_species(electron_idx)%drift_velocity(3)
    ion_drift_mps = -app%particle_species(ion_idx)%drift_velocity(3)
    area = (app%sim%box_max(1) - app%sim%box_min(1))*(app%sim%box_max(2) - app%sim%box_min(2))
    if (app%surface_current%has_reference_area_m2) area = app%surface_current%reference_area_m2

    call build_zhao_params( &
      app%surface_current%solar_elevation_deg, &
      species_number_density_m3(app%particle_species(ion_idx)), &
      app%surface_current%photoelectron_ref_density_m3, &
      electron_temperature_ev, photo_temperature_ev, electron_drift_mps, ion_drift_mps, &
      app%particle_species(ion_idx)%m_particle, app%particle_species(electron_idx)%m_particle, params, &
      photoelectron_source_scale=app%surface_current%photoelectron_source_scale &
      )
    solver_name = 'zhao_'//trim(lower_ascii(app%surface_current%zhao_branch))
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
    if (result%photoelectron_return_current_density_a_m2 > 0.0_dp .or. &
        result%photoelectron_escape_current_density_a_m2 < 0.0_dp) then
      error stop 'Zhao stationary surface-current evaluation produced invalid PE current signs.'
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
    result%has_absorbed_target([electron_idx, ion_idx, photo_idx]) = .true.
    result%has_emission_target(photo_idx) = .true.
    result%has_escape_target(photo_idx) = .true.
    result%absorbed_current_a(electron_idx) = area*result%electron_current_density_a_m2
    result%absorbed_current_a(ion_idx) = area*result%ion_current_density_a_m2
    result%absorbed_current_a(photo_idx) = area*result%photoelectron_return_current_density_a_m2
    result%emission_current_a(photo_idx) = area*result%photoelectron_emission_current_density_a_m2
    ! escaped_to_infinity は粒子電荷の外向きfluxなので、正の表面帯電電流とは符号が逆。
    result%escaped_particle_current_a(photo_idx) = -area*result%photoelectron_escape_current_density_a_m2
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

end module bem_surface_current_model
