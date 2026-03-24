!> app_config 非依存のシース core を検証するテスト。
program test_sheath_model_core
  use bem_kinds, only: dp
  use bem_sheath_model_core, only: &
    sheath_model_species, zhao_params_type, zhao_local_state_type, &
    solve_no_photo_floating_potential, build_zhao_params, solve_zhao_unknowns, sample_zhao_state_at_z, &
    resolve_species_drift_speed, zhao_electron_vmin_normal, zhao_photo_emit_current_density
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_close_dp
  implicit none

  real(dp), parameter :: ev_to_k = 1.160451812d4

  type(sheath_model_species) :: spec_e, spec_i, spec_p
  type(zhao_params_type) :: params
  type(zhao_local_state_type) :: state
  real(dp) :: inward_normal(3), phi0_v, phi_m_v, n_swe_inf_m3
  real(dp) :: v_d_electron_mps, v_d_ion_mps
  character(len=1) :: branch

  call test_init(3)

  inward_normal = [0.0d0, 0.0d0, -1.0d0]

  call test_begin('zhao_type_a')
  call configure_zhao_fixture(60.0d0, spec_e, spec_i, spec_p)
  v_d_electron_mps = resolve_species_drift_speed(spec_e, 'normal', inward_normal)
  v_d_ion_mps = resolve_species_drift_speed(spec_i, 'normal', inward_normal)
  call build_zhao_params( &
    60.0d0, spec_i%number_density_m3, 64.0d0*1.0d6, spec_e%temperature_k/ev_to_k, spec_p%temperature_k/ev_to_k, &
    v_d_electron_mps, v_d_ion_mps, spec_i%m_particle, abs(spec_e%m_particle), params &
    )
  call solve_zhao_unknowns('zhao_auto', params, phi0_v, phi_m_v, n_swe_inf_m3, branch)
  call assert_true(branch == 'A', 'core alpha=60 should select Type A')
  call assert_close_dp(phi0_v, 2.9712182827319435d0, 5.0d-6, 'core Type-A phi0 mismatch')
  call assert_close_dp(phi_m_v, -0.8169121871620854d0, 5.0d-6, 'core Type-A phi_m mismatch')
  call assert_close_dp(n_swe_inf_m3, 7.819215729579456d6, 5.0d1, 'core Type-A n_swe_inf mismatch')
  call assert_close_dp( &
    zhao_photo_emit_current_density(branch, phi0_v, phi_m_v, spec_p%q_particle, params), &
    3.9386846806723257d-7, 5.0d-13, 'core Type-A photo current mismatch' &
    )
  call assert_close_dp( &
    zhao_electron_vmin_normal(branch, phi0_v, phi_m_v, abs(spec_e%m_particle)), &
    5.360599783278604d5, 5.0d-1, 'core Type-A electron cutoff mismatch' &
    )
  call sample_zhao_state_at_z(params, branch, phi0_v, phi_m_v, n_swe_inf_m3, 4.0d0, state)
  call assert_close_dp(state%phi_v, -4.047310136961665d-2, 5.0d-8, 'core Type-A local phi mismatch')
  call assert_close_dp(state%electron_source_density_m3, 7.792887827459829d6, 5.0d1, 'core Type-A local electron density mismatch')
  call test_end()

  call test_begin('zhao_type_c')
  call configure_zhao_fixture(10.0d0, spec_e, spec_i, spec_p)
  v_d_electron_mps = resolve_species_drift_speed(spec_e, 'normal', inward_normal)
  v_d_ion_mps = resolve_species_drift_speed(spec_i, 'normal', inward_normal)
  call build_zhao_params( &
    10.0d0, spec_i%number_density_m3, 64.0d0*1.0d6, spec_e%temperature_k/ev_to_k, spec_p%temperature_k/ev_to_k, &
    v_d_electron_mps, v_d_ion_mps, spec_i%m_particle, abs(spec_e%m_particle), params &
    )
  call solve_zhao_unknowns('zhao_auto', params, phi0_v, phi_m_v, n_swe_inf_m3, branch)
  call assert_true(branch == 'C', 'core alpha=10 should select Type C')
  call sample_zhao_state_at_z(params, branch, phi0_v, phi_m_v, n_swe_inf_m3, 4.0d0, state)
  call assert_true(state%swe_reflected_active, 'core Type-C should activate reflected SWE branch')
  call assert_close_dp(state%phi_v, -3.43310064584225d0, 1.0d-4, 'core Type-C local phi mismatch')
  call assert_close_dp(state%electron_source_density_m3, 6.136203306856032d6, 1.0d2, 'core Type-C local electron density mismatch')
  call assert_close_dp(state%v_i_mps, 8.521786009033145d4, 1.0d-1, 'core Type-C local ion speed mismatch')
  call test_end()

  call test_begin('no_photo_floating')
  call configure_no_photo_fixture(spec_e, spec_i)
  call solve_no_photo_floating_potential(spec_e, spec_i, inward_normal, phi0_v)
  call assert_true(phi0_v < 0.0d0, 'core floating_no_photo phi0 should be negative')
  call test_end()

  call test_summary()

contains

  subroutine configure_zhao_fixture(alpha_deg, spec_e, spec_i, spec_p)
    real(dp), intent(in) :: alpha_deg
    type(sheath_model_species), intent(out) :: spec_e, spec_i, spec_p
    real(dp) :: speed

    speed = 4.68d5*sin(alpha_deg*acos(-1.0d0)/180.0d0)

    spec_e%q_particle = -1.602176634d-19
    spec_e%m_particle = 9.1093837015d-31
    spec_e%number_density_m3 = 8.7d0*1.0d6
    spec_e%temperature_k = 12.0d0*ev_to_k
    spec_e%drift_velocity = [0.0d0, 0.0d0, -speed]

    spec_i%q_particle = 1.602176634d-19
    spec_i%m_particle = 1.67262192369d-27
    spec_i%number_density_m3 = 8.7d0*1.0d6
    spec_i%temperature_k = 0.0d0
    spec_i%drift_velocity = [0.0d0, 0.0d0, -speed]

    spec_p%q_particle = -1.602176634d-19
    spec_p%m_particle = 9.1093837015d-31
    spec_p%temperature_k = 2.2d0*ev_to_k
  end subroutine configure_zhao_fixture

  subroutine configure_no_photo_fixture(spec_e, spec_i)
    type(sheath_model_species), intent(out) :: spec_e, spec_i

    spec_e%q_particle = -1.602176634d-19
    spec_e%m_particle = 9.1093837015d-31
    spec_e%number_density_m3 = 5.0d0*1.0d6
    spec_e%temperature_k = 10.0d0*ev_to_k
    spec_e%drift_velocity = [0.0d0, 0.0d0, -4.0d5]

    spec_i%q_particle = 1.602176634d-19
    spec_i%m_particle = 1.67262192369d-27
    spec_i%number_density_m3 = 5.0d0*1.0d6
    spec_i%temperature_k = 0.0d0
    spec_i%drift_velocity = [0.0d0, 0.0d0, -4.0d5]
  end subroutine configure_no_photo_fixture

end program test_sheath_model_core
