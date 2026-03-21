!> シース数値モデルと app_config / 注入ランタイムの橋渡しを行うモジュール。
module bem_sheath_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: sim_config
  use bem_app_config_types, only: app_config, particle_species_spec
  use bem_app_config_parser, only: &
    lower, resolve_inject_face, resolve_inward_normal, species_number_density_m3, species_temperature_k
  use bem_sheath_model_core, only: &
    sheath_model_species, zhao_params_type, zhao_local_state_type, &
    solve_no_photo_floating_potential, build_zhao_params, solve_zhao_unknowns, sample_zhao_state_at_z, &
    resolve_species_drift_speed, zhao_electron_vmin_normal, zhao_photo_vmin_normal, zhao_photo_emit_current_density
  implicit none

  real(dp), parameter :: qe = 1.602176634d-19

  type :: sheath_injection_context
    logical :: enabled = .false.
    logical :: has_photo_species = .false.
    logical :: has_local_reservoir_profile = .false.
    character(len=32) :: model = 'none'
    character(len=1) :: branch = ' '
    integer(i32) :: electron_species = 0_i32
    integer(i32) :: ion_species = 0_i32
    integer(i32) :: photo_species = 0_i32
    character(len=16) :: reference_face = ''
    integer(i32) :: reference_axis = 0_i32
    real(dp) :: reference_coordinate = 0.0d0
    real(dp) :: reference_inward_normal(3) = 0.0d0
    real(dp) :: reservoir_plane_distance_m = 0.0d0
    real(dp) :: reservoir_phi_v = 0.0d0
    real(dp) :: phi0_v = 0.0d0
    real(dp) :: phi_m_v = 0.0d0
    real(dp) :: n_swe_inf_m3 = 0.0d0
    real(dp) :: electron_number_density_m3 = 0.0d0
    real(dp) :: electron_vmin_normal = 0.0d0
    real(dp) :: ion_number_density_m3 = 0.0d0
    real(dp) :: ion_normal_speed_mps = 0.0d0
    real(dp) :: photo_emit_current_density_a_m2 = 0.0d0
    real(dp) :: photo_vmin_normal = 0.0d0
  end type sheath_injection_context

  public :: sheath_injection_context
  public :: resolve_sheath_injection_context

contains

  subroutine resolve_sheath_injection_context(cfg, ctx)
    type(app_config), intent(in) :: cfg
    type(sheath_injection_context), intent(out) :: ctx

    integer(i32) :: electron_idx, ion_idx, photo_idx
    type(particle_species_spec) :: spec_e_cfg, spec_i_cfg, spec_p_cfg
    type(sheath_model_species) :: spec_e, spec_i, spec_p
    type(zhao_params_type) :: params
    type(zhao_local_state_type) :: local_state
    real(dp) :: phi0_v, phi_m_v, n_swe_inf_m3
    real(dp) :: n_e_inf_m3, n_i_inf_m3, n_phe_ref_m3
    real(dp) :: t_swe_ev, t_phe_ev, v_d_electron_mps, v_d_ion_mps
    real(dp) :: reference_coordinate, reference_inward_normal(3)
    integer :: reference_axis
    character(len=32) :: model
    character(len=1) :: branch

    ctx = sheath_injection_context()
    model = trim(lower(cfg%sim%sheath_injection_model))
    ctx%model = model
    if (model == 'none') return

    call detect_sheath_species(cfg, electron_idx, ion_idx, photo_idx)
    spec_e_cfg = cfg%particle_species(electron_idx)
    spec_i_cfg = cfg%particle_species(ion_idx)
    spec_e = to_sheath_model_species(spec_e_cfg)
    spec_i = to_sheath_model_species(spec_i_cfg)

    ctx%enabled = .true.
    ctx%electron_species = electron_idx
    ctx%ion_species = ion_idx
    ctx%reference_face = trim(lower(spec_e_cfg%inject_face))
    call resolve_sheath_reference_plane( &
      cfg%sim, ctx%reference_face, reference_axis, reference_coordinate, reference_inward_normal &
      )
    ctx%reference_axis = int(reference_axis, i32)
    ctx%reference_coordinate = reference_coordinate
    ctx%reference_inward_normal = reference_inward_normal

    n_e_inf_m3 = spec_e%number_density_m3
    n_i_inf_m3 = spec_i%number_density_m3
    t_swe_ev = spec_e%temperature_k*k_boltzmann/qe
    v_d_electron_mps = resolve_species_drift_speed(spec_e, trim(lower(cfg%sim%sheath_electron_drift_mode)), reference_inward_normal)
    v_d_ion_mps = resolve_species_drift_speed(spec_i, trim(lower(cfg%sim%sheath_ion_drift_mode)), reference_inward_normal)

    if (v_d_electron_mps <= 0.0d0) error stop 'sheath electron drift must point inward.'
    if (v_d_ion_mps <= 0.0d0) error stop 'sheath ion drift must point inward.'

    select case (model)
    case ('floating_no_photo')
      call solve_no_photo_floating_potential(spec_e, spec_i, reference_inward_normal, phi0_v)
      ctx%branch = 'N'
      ctx%phi0_v = phi0_v
      ctx%phi_m_v = phi0_v
      ctx%n_swe_inf_m3 = n_e_inf_m3
      ctx%electron_number_density_m3 = n_e_inf_m3
      ctx%electron_vmin_normal = sqrt(max(0.0d0, 2.0d0*qe*max(0.0d0, -phi0_v))/abs(spec_e%m_particle))
      return
    case ('zhao_auto', 'zhao_a', 'zhao_b', 'zhao_c')
      if (photo_idx <= 0_i32) then
        error stop 'Zhao sheath injection requires one enabled negative-q photo_raycast species.'
      end if
      spec_p_cfg = cfg%particle_species(photo_idx)
      spec_p = to_sheath_model_species(spec_p_cfg)
      if (abs(spec_p_cfg%q_particle) <= 0.0d0) error stop 'Zhao photoelectron species must have non-zero q_particle.'
      ctx%photo_species = photo_idx
      ctx%has_photo_species = .true.
      t_phe_ev = spec_p%temperature_k*k_boltzmann/qe
      n_phe_ref_m3 = cfg%sim%sheath_photoelectron_ref_density_cm3*1.0d6

      call build_zhao_params( &
        cfg%sim%sheath_alpha_deg, n_i_inf_m3, n_phe_ref_m3, t_swe_ev, t_phe_ev, v_d_electron_mps, v_d_ion_mps, &
        spec_i%m_particle, abs(spec_e%m_particle), params &
        )
      call solve_zhao_unknowns(model, params, phi0_v, phi_m_v, n_swe_inf_m3, branch)

      ctx%branch = branch
      ctx%phi0_v = phi0_v
      ctx%phi_m_v = phi_m_v
      ctx%n_swe_inf_m3 = n_swe_inf_m3
      ctx%electron_number_density_m3 = n_swe_inf_m3
      ctx%electron_vmin_normal = zhao_electron_vmin_normal(branch, phi0_v, phi_m_v, abs(spec_e%m_particle))
      ctx%photo_vmin_normal = zhao_photo_vmin_normal(branch, phi0_v, phi_m_v, abs(spec_p%m_particle))
      ctx%photo_emit_current_density_a_m2 = zhao_photo_emit_current_density(branch, phi0_v, phi_m_v, spec_p_cfg%q_particle, params)

      if (cfg%sim%has_sheath_reference_coordinate) then
        call sample_zhao_reservoir_state(cfg%sim, spec_e_cfg, ctx, params, branch, phi0_v, phi_m_v, n_swe_inf_m3, local_state)
        ctx%has_local_reservoir_profile = .true.
        ctx%reservoir_phi_v = local_state%phi_v
        ctx%electron_number_density_m3 = local_state%electron_source_density_m3
        if (local_state%swe_reflected_active) then
          ctx%electron_vmin_normal = 0.0d0
        else
          ctx%electron_vmin_normal = local_state%vcut_swe_mps
        end if
        ctx%ion_number_density_m3 = local_state%n_swi_m3
        ctx%ion_normal_speed_mps = local_state%v_i_mps
      end if
    case default
      error stop 'Unknown sim.sheath_injection_model in runtime.'
    end select
  end subroutine resolve_sheath_injection_context

  subroutine detect_sheath_species(cfg, electron_idx, ion_idx, photo_idx)
    type(app_config), intent(in) :: cfg
    integer(i32), intent(out) :: electron_idx, ion_idx, photo_idx

    integer(i32) :: s
    character(len=16) :: mode

    electron_idx = 0_i32
    ion_idx = 0_i32
    photo_idx = 0_i32
    do s = 1_i32, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      mode = trim(lower(cfg%particle_species(s)%source_mode))
      select case (mode)
      case ('reservoir_face')
        if (cfg%particle_species(s)%q_particle < 0.0d0) then
          if (electron_idx == 0_i32) electron_idx = s
        else if (cfg%particle_species(s)%q_particle > 0.0d0) then
          if (ion_idx == 0_i32) ion_idx = s
        end if
      case ('photo_raycast')
        if (cfg%particle_species(s)%q_particle < 0.0d0 .and. photo_idx == 0_i32) photo_idx = s
      end select
    end do

    if (electron_idx <= 0_i32) error stop 'sheath injection requires one enabled negative-q reservoir_face species.'
    if (ion_idx <= 0_i32) error stop 'sheath injection requires one enabled positive-q reservoir_face species.'
    if (trim(lower(cfg%particle_species(electron_idx)%inject_face)) /= trim(lower(cfg%particle_species(ion_idx)%inject_face))) then
      error stop 'sheath electron/ion reservoir species must share the same inject_face.'
    end if
  end subroutine detect_sheath_species

  subroutine sample_zhao_reservoir_state(sim, spec_e, ctx, p, branch, phi0_v, phi_m_v, n_swe_inf_m3, state)
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec_e
    type(sheath_injection_context), intent(inout) :: ctx
    type(zhao_params_type), intent(in) :: p
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, n_swe_inf_m3
    type(zhao_local_state_type), intent(out) :: state

    integer :: axis
    real(dp) :: boundary_value, distance_m, axis_sign

    call resolve_inject_face(sim%box_min, sim%box_max, spec_e%inject_face, axis, boundary_value)
    axis_sign = -ctx%reference_inward_normal(axis)
    distance_m = (boundary_value - ctx%reference_coordinate)*axis_sign
    if (distance_m < -1.0d-12) then
      error stop 'sheath_reference_coordinate must lie on or inside the shared reservoir_face boundary.'
    end if
    ctx%reservoir_plane_distance_m = max(0.0d0, distance_m)
    call sample_zhao_state_at_z(p, branch, phi0_v, phi_m_v, n_swe_inf_m3, ctx%reservoir_plane_distance_m, state)
  end subroutine sample_zhao_reservoir_state

  subroutine resolve_sheath_reference_plane(sim, inject_face, axis, reference_coordinate, inward_normal)
    type(sim_config), intent(in) :: sim
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis
    real(dp), intent(out) :: reference_coordinate
    real(dp), intent(out) :: inward_normal(3)
    real(dp) :: boundary_value

    call resolve_inward_normal(inject_face, inward_normal)
    call resolve_inject_face(sim%box_min, sim%box_max, inject_face, axis, boundary_value)
    if (sim%has_sheath_reference_coordinate) then
      reference_coordinate = sim%sheath_reference_coordinate
    else
      reference_coordinate = boundary_value
    end if
  end subroutine resolve_sheath_reference_plane

  type(sheath_model_species) function to_sheath_model_species(spec_cfg) result(spec)
    type(particle_species_spec), intent(in) :: spec_cfg

    spec%q_particle = spec_cfg%q_particle
    spec%m_particle = spec_cfg%m_particle
    if (spec_cfg%has_number_density_m3 .or. spec_cfg%has_number_density_cm3) then
      spec%number_density_m3 = species_number_density_m3(spec_cfg)
    end if
    spec%temperature_k = species_temperature_k(spec_cfg)
    spec%drift_velocity = spec_cfg%drift_velocity
  end function to_sheath_model_species

end module bem_sheath_runtime
