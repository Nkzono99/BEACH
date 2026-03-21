!> Zhao/無光電子シース注入コンテキストの解決を検証するテスト。
program test_sheath_injection_model
  use bem_kinds, only: dp
  use bem_app_config, only: app_config, default_app_config
  use bem_sheath_injection_model, only: sheath_injection_context, resolve_sheath_injection_context
  use bem_injection, only: compute_inflow_flux_from_drifting_maxwellian
  use bem_app_config_parser, only: resolve_inward_normal
  use test_support, only: assert_true, assert_close_dp, assert_allclose_1d
  implicit none

  type(app_config) :: cfg
  type(sheath_injection_context) :: ctx
  real(dp) :: gamma_e, gamma_i, inward_normal(3)

  call configure_zhao_fixture(cfg, 60.0d0)
  call resolve_sheath_injection_context(cfg, ctx)

  call assert_true(ctx%enabled, 'zhao_auto should enable sheath context')
  call assert_true(ctx%branch == 'A', 'alpha=60 Zhao auto should select Type A')
  call assert_close_dp(ctx%phi0_v, 2.9712182827319435d0, 5.0d-6, 'Type-A phi0 mismatch')
  call assert_close_dp(ctx%phi_m_v, -0.8169121871620854d0, 5.0d-6, 'Type-A phi_m mismatch')
  call assert_close_dp(ctx%n_swe_inf_m3, 7.819215729579456d6, 5.0d1, 'Type-A n_swe_inf mismatch')
  call assert_close_dp( &
    ctx%photo_emit_current_density_a_m2, 3.9386846806723257d-7, 5.0d-13, 'Type-A photo current mismatch' &
    )
  call assert_true(trim(ctx%reference_face) == 'z_high', 'default sheath reference face mismatch')
  call assert_true(ctx%reference_axis == 3, 'default sheath reference axis mismatch')
  call assert_close_dp(ctx%reference_coordinate, 4.0d0, 1.0d-12, 'default sheath reference coordinate mismatch')
  call assert_allclose_1d( &
    ctx%reference_inward_normal, [0.0d0, 0.0d0, -1.0d0], 1.0d-12, 'default sheath reference normal mismatch' &
    )

  call configure_zhao_fixture(cfg, 60.0d0)
  cfg%sim%sheath_reference_coordinate = 0.0d0
  cfg%sim%has_sheath_reference_coordinate = .true.
  call resolve_sheath_injection_context(cfg, ctx)
  call assert_close_dp(ctx%reference_coordinate, 0.0d0, 1.0d-12, 'explicit sheath reference coordinate mismatch')
  call assert_true(ctx%has_local_reservoir_profile, 'explicit reference should enable local reservoir profile')
  call assert_close_dp(ctx%reservoir_plane_distance_m, 4.0d0, 1.0d-12, 'Type-A reservoir plane distance mismatch')
  call assert_close_dp(ctx%reservoir_phi_v, -4.047310136961665d-2, 5.0d-8, 'Type-A local phi mismatch')
  call assert_close_dp(ctx%electron_number_density_m3, 7.792887827459829d6, 5.0d1, 'Type-A local electron source density mismatch')
  call assert_close_dp(ctx%electron_vmin_normal, 5.2261201693750557d5, 5.0d-1, 'Type-A local electron cutoff mismatch')
  call assert_close_dp(ctx%ion_number_density_m3, 8.699794680592455d6, 5.0d1, 'Type-A local ion density mismatch')
  call assert_close_dp(ctx%ion_normal_speed_mps, 4.053094542466367d5, 5.0d-1, 'Type-A local ion speed mismatch')

  call configure_zhao_fixture(cfg, 10.0d0)
  call resolve_sheath_injection_context(cfg, ctx)
  call assert_true(ctx%branch == 'C', 'alpha=10 Zhao auto should select Type C')
  call assert_close_dp(ctx%phi0_v, -4.798298347150015d0, 5.0d-6, 'Type-C phi0 mismatch')
  call assert_close_dp(ctx%n_swe_inf_m3, 8.168603119546532d6, 5.0d1, 'Type-C n_swe_inf mismatch')
  call assert_true(ctx%photo_emit_current_density_a_m2 > 0.0d0, 'Type-C photo current should be positive')

  call configure_zhao_fixture(cfg, 10.0d0)
  cfg%sim%sheath_reference_coordinate = 0.0d0
  cfg%sim%has_sheath_reference_coordinate = .true.
  call resolve_sheath_injection_context(cfg, ctx)
  call assert_true(ctx%has_local_reservoir_profile, 'Type-C explicit reference should enable local reservoir profile')
  call assert_close_dp(ctx%reservoir_phi_v, -3.43310064584225d0, 1.0d-4, 'Type-C local phi mismatch')
  call assert_close_dp(ctx%electron_number_density_m3, 6.136203306856032d6, 1.0d2, 'Type-C local electron source density mismatch')
  call assert_close_dp(ctx%electron_vmin_normal, 0.0d0, 1.0d-12, 'Type-C reflected electrons should fill the low-speed range')
  call assert_close_dp(ctx%ion_number_density_m3, 8.296687096334287d6, 2.0d1, 'Type-C local ion density mismatch')
  call assert_close_dp(ctx%ion_normal_speed_mps, 8.521786009033145d4, 1.0d-1, 'Type-C local ion speed mismatch')

  call configure_no_photo_fixture(cfg)
  call resolve_sheath_injection_context(cfg, ctx)
  call assert_true(ctx%branch == 'N', 'floating_no_photo should return synthetic branch N')
  call assert_true(ctx%phi0_v < 0.0d0, 'floating_no_photo phi0 should be negative')
  call assert_close_dp(ctx%photo_emit_current_density_a_m2, 0.0d0, 1.0d-18, 'floating_no_photo photo current mismatch')

  call resolve_inward_normal(cfg%particle_species(1)%inject_face, inward_normal)
  gamma_e = compute_inflow_flux_from_drifting_maxwellian( &
            cfg%particle_species(1)%number_density_cm3*1.0d6, cfg%particle_species(1)%temperature_ev*1.160451812d4, &
            cfg%particle_species(1)%m_particle, cfg%particle_species(1)%drift_velocity, inward_normal, &
            vmin_normal=ctx%electron_vmin_normal &
            )
  gamma_i = compute_inflow_flux_from_drifting_maxwellian( &
            cfg%particle_species(2)%number_density_cm3*1.0d6, cfg%particle_species(2)%temperature_ev*1.160451812d4, &
            cfg%particle_species(2)%m_particle, cfg%particle_species(2)%drift_velocity, inward_normal &
            )
  call assert_close_dp(gamma_e, gamma_i, 1.0d-6*max(1.0d0, gamma_i), 'floating_no_photo current balance mismatch')

contains

  subroutine configure_zhao_fixture(cfg, alpha_deg)
    type(app_config), intent(out) :: cfg
    real(dp), intent(in) :: alpha_deg
    real(dp) :: speed

    call default_app_config(cfg)
    cfg%sim%sheath_injection_model = 'zhao_auto'
    cfg%sim%sheath_alpha_deg = alpha_deg
    cfg%sim%sheath_photoelectron_ref_density_cm3 = 64.0d0
    cfg%sim%sheath_electron_drift_mode = 'normal'
    cfg%sim%sheath_ion_drift_mode = 'normal'
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0d0, 0.0d0, 0.0d0]
    cfg%sim%box_max = [1.0d0, 1.0d0, 4.0d0]
    cfg%n_particle_species = 3

    speed = 4.68d5*sin(alpha_deg*acos(-1.0d0)/180.0d0)

    cfg%particle_species(1)%enabled = .true.
    cfg%particle_species(1)%source_mode = 'reservoir_face'
    cfg%particle_species(1)%q_particle = -1.602176634d-19
    cfg%particle_species(1)%m_particle = 9.1093837015d-31
    cfg%particle_species(1)%number_density_cm3 = 8.7d0
    cfg%particle_species(1)%has_number_density_cm3 = .true.
    cfg%particle_species(1)%temperature_ev = 12.0d0
    cfg%particle_species(1)%has_temperature_ev = .true.
    cfg%particle_species(1)%inject_face = 'z_high'
    cfg%particle_species(1)%drift_velocity = [0.0d0, 0.0d0, -speed]

    cfg%particle_species(2)%enabled = .true.
    cfg%particle_species(2)%source_mode = 'reservoir_face'
    cfg%particle_species(2)%q_particle = 1.602176634d-19
    cfg%particle_species(2)%m_particle = 1.67262192369d-27
    cfg%particle_species(2)%number_density_cm3 = 8.7d0
    cfg%particle_species(2)%has_number_density_cm3 = .true.
    cfg%particle_species(2)%temperature_ev = 0.0d0
    cfg%particle_species(2)%has_temperature_ev = .true.
    cfg%particle_species(2)%inject_face = 'z_high'
    cfg%particle_species(2)%drift_velocity = [0.0d0, 0.0d0, -speed]

    cfg%particle_species(3)%enabled = .true.
    cfg%particle_species(3)%source_mode = 'photo_raycast'
    cfg%particle_species(3)%q_particle = -1.602176634d-19
    cfg%particle_species(3)%m_particle = 9.1093837015d-31
    cfg%particle_species(3)%temperature_ev = 2.2d0
    cfg%particle_species(3)%has_temperature_ev = .true.
    cfg%particle_species(3)%inject_face = 'z_high'
  end subroutine configure_zhao_fixture

  subroutine configure_no_photo_fixture(cfg)
    type(app_config), intent(out) :: cfg

    call default_app_config(cfg)
    cfg%sim%sheath_injection_model = 'floating_no_photo'
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0d0, 0.0d0, 0.0d0]
    cfg%sim%box_max = [1.0d0, 1.0d0, 4.0d0]
    cfg%n_particle_species = 2

    cfg%particle_species(1)%enabled = .true.
    cfg%particle_species(1)%source_mode = 'reservoir_face'
    cfg%particle_species(1)%q_particle = -1.602176634d-19
    cfg%particle_species(1)%m_particle = 9.1093837015d-31
    cfg%particle_species(1)%number_density_cm3 = 5.0d0
    cfg%particle_species(1)%has_number_density_cm3 = .true.
    cfg%particle_species(1)%temperature_ev = 10.0d0
    cfg%particle_species(1)%has_temperature_ev = .true.
    cfg%particle_species(1)%inject_face = 'z_high'
    cfg%particle_species(1)%drift_velocity = [0.0d0, 0.0d0, -4.0d5]

    cfg%particle_species(2)%enabled = .true.
    cfg%particle_species(2)%source_mode = 'reservoir_face'
    cfg%particle_species(2)%q_particle = 1.602176634d-19
    cfg%particle_species(2)%m_particle = 1.67262192369d-27
    cfg%particle_species(2)%number_density_cm3 = 5.0d0
    cfg%particle_species(2)%has_number_density_cm3 = .true.
    cfg%particle_species(2)%temperature_ev = 0.0d0
    cfg%particle_species(2)%has_temperature_ev = .true.
    cfg%particle_species(2)%inject_face = 'z_high'
    cfg%particle_species(2)%drift_velocity = [0.0d0, 0.0d0, -4.0d5]
  end subroutine configure_no_photo_fixture

end program test_sheath_injection_model
