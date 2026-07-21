program test_outer_plasma_kinetic_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe, k_boltzmann
  use bem_app_config_types, only: app_config, default_app_config
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_not_applicable
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type
  use bem_outer_plasma_kinetic_runtime, only: resolve_kinetic_outer_options
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32, assert_true
  implicit none

  type(app_config) :: app
  type(kinetic_outer_plasma_options_type) :: options
  integer(i32) :: status
  character(len=256) :: message

  call test_init(13)

  call test_begin('runtime adapter resolves the ambient reservoir VDFs')
  call init_fixture(app)
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'kinetic runtime resolution status mismatch')
  call assert_close_dp(options%interface_field, -0.5_dp, 0.0_dp, 'interface field mismatch')
  call assert_close_dp(options%electron_density_infinity, 2.0e6_dp, 0.0_dp, 'electron density mismatch')
  call assert_close_dp(options%electron_temperature_j, 3.0_dp*qe, 1.0e-32_dp, 'electron temperature mismatch')
  call assert_close_dp(options%ion_density_infinity, 2.0e6_dp, 0.0_dp, 'ion density mismatch')
  call assert_close_dp(options%ion_temperature_j, 500.0_dp*k_boltzmann, 1.0e-32_dp, 'ion temperature mismatch')
  call assert_close_dp(options%ion_drift_infinity, 4.0e4_dp, 0.0_dp, 'ion inward drift mismatch')
  call assert_close_dp(options%domain_length, 10.0_dp*app%outer_plasma%debye_length, 0.0_dp, 'domain length mismatch')
  call test_end()

  call test_begin('runtime adapter fails closed without an ion reservoir')
  app%particle_species(2)%enabled = .false.
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, 'missing ion reservoir must be not_applicable')
  call test_end()

  call test_begin('kinetic mean requires an explicit photoelectron source VDF')
  app%particle_species(2)%enabled = .true.
  app%outer_plasma%photoelectron_density_model = 'kinetic_mean'
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, 'missing photoelectron VDF must be not_applicable')
  call test_end()

  call test_begin('runtime adapter resolves Zhao charge-driven inputs without replacing emission current')
  call init_fixture(app)
  app%outer_plasma%kinetic_closure = 'zhao_charge_driven'
  app%outer_plasma%zhao_branch = 'c'
  app%sim%sheath_alpha_deg = 10.0_dp
  app%sim%sheath_photoelectron_ref_density_cm3 = 64.0_dp
  app%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.5e5_dp]
  app%n_particle_species = 3_i32
  app%particle_species(3)%enabled = .true.
  app%particle_species(3)%source_mode = 'photo_raycast'
  app%particle_species(3)%q_particle = -qe
  app%particle_species(3)%m_particle = 9.1093837139e-31_dp
  app%particle_species(3)%temperature_ev = 2.2_dp
  app%particle_species(3)%has_temperature_ev = .true.
  app%particle_species(3)%emit_current_density_a_m2 = &
    qe*64.0e6_dp*sin(10.0_dp*acos(-1.0_dp)/180.0_dp)* &
    sqrt(2.0_dp*2.2_dp*qe/app%particle_species(3)%m_particle)/(2.0_dp*sqrt(acos(-1.0_dp)))
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'Zhao runtime resolution status mismatch: '//trim(message))
  call assert_true(trim(options%kinetic_closure) == 'zhao_charge_driven', 'Zhao closure was not propagated')
  call assert_true(trim(options%zhao_branch) == 'c', 'Zhao branch was not propagated')
  call assert_close_dp(options%electron_drift_infinity, 1.5e5_dp, 0.0_dp, 'Zhao electron drift mismatch')
  call assert_close_dp(options%photoelectron_reference_density, 64.0e6_dp, 0.0_dp, 'Zhao reference density mismatch')
  call assert_close_dp(options%photoelectron_source_scale, 1.0_dp, 0.0_dp, 'Zhao source scale mismatch')
  call assert_close_dp(options%photoelectron_temperature_j, 2.2_dp*qe, 1.0e-32_dp, 'Zhao photoelectron temperature mismatch')
  call assert_close_dp(options%photoelectron_emission_flux, 0.0_dp, 0.0_dp, 'Zhao closure must not replace tracked emission')
  call test_end()

  call test_begin('runtime adapter scales the tracked Zhao source independently')
  app%outer_plasma%photoelectron_source_scale = 0.5_dp
  app%particle_species(3)%emit_current_density_a_m2 = 0.5_dp*app%particle_species(3)%emit_current_density_a_m2
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'scaled Zhao source resolution failed: '//trim(message))
  call assert_close_dp(options%photoelectron_source_scale, 0.5_dp, 0.0_dp, 'scaled Zhao source mismatch')
  call test_end()

  call test_begin('runtime adapter rejects inconsistent tracked Zhao emission')
  app%particle_species(3)%emit_current_density_a_m2 = 2.0_dp*app%particle_species(3)%emit_current_density_a_m2
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, &
                        'inconsistent Zhao tracked emission must be not_applicable')
  call test_end()

  call test_begin('runtime adapter derives the no-photo Zhao scale from upstream plasma')
  call init_fixture(app)
  app%outer_plasma%kinetic_closure = 'zhao_charge_driven'
  app%outer_plasma%zhao_branch = 'auto'
  app%outer_plasma%photoelectron_source_scale = 0.0_dp
  app%sim%sheath_photoelectron_ref_density_cm3 = 0.0_dp
  app%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.5e5_dp]
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'no-photo Zhao runtime resolution failed: '//trim(message))
  call assert_close_dp(options%photoelectron_source_scale, 0.0_dp, 0.0_dp, 'no-photo source scale mismatch')
  call assert_close_dp( &
    options%photoelectron_reference_density, options%ion_density_infinity, 0.0_dp, &
    'no-photo Zhao normalization density mismatch' &
    )
  call assert_close_dp( &
    options%photoelectron_temperature_j, options%electron_temperature_j, 0.0_dp, &
    'no-photo Zhao normalization temperature mismatch' &
    )
  call assert_close_dp(options%photoelectron_population_fraction, 0.0_dp, 0.0_dp, &
                       'no-photo Zhao population must vanish')
  call test_end()

  call test_begin('runtime adapter rejects an enabled source in no-photo Zhao mode')
  app%n_particle_species = 3_i32
  app%particle_species(3)%enabled = .true.
  app%particle_species(3)%source_mode = 'photo_raycast'
  app%particle_species(3)%q_particle = -qe
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, &
                        'no-photo Zhao must reject an enabled photoelectron source')
  call test_end()

  call test_begin('runtime adapter enforces upstream Zhao quasineutrality')
  app%n_particle_species = 2_i32
  app%particle_species(2)%number_density_m3 = 2.1e6_dp
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, &
                        'non-neutral Zhao upstream state must be not_applicable')
  call test_end()

  call test_begin('runtime adapter rejects a shifted Zhao photoelectron VDF')
  call init_fixture(app)
  app%outer_plasma%kinetic_closure = 'zhao_charge_driven'
  app%outer_plasma%zhao_branch = 'c'
  app%sim%sheath_alpha_deg = 10.0_dp
  app%sim%sheath_photoelectron_ref_density_cm3 = 64.0_dp
  app%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.5e5_dp]
  app%n_particle_species = 3_i32
  app%particle_species(3)%enabled = .true.
  app%particle_species(3)%source_mode = 'photo_raycast'
  app%particle_species(3)%q_particle = -qe
  app%particle_species(3)%m_particle = 9.1093837139e-31_dp
  app%particle_species(3)%temperature_ev = 2.2_dp
  app%particle_species(3)%has_temperature_ev = .true.
  app%particle_species(3)%emit_current_density_a_m2 = &
    qe*64.0e6_dp*sin(10.0_dp*acos(-1.0_dp)/180.0_dp)* &
    sqrt(2.0_dp*2.2_dp*qe/app%particle_species(3)%m_particle)/(2.0_dp*sqrt(acos(-1.0_dp)))
  app%particle_species(3)%emit_current_density_a_m2 = 0.5_dp*app%particle_species(3)%emit_current_density_a_m2
  app%particle_species(3)%normal_drift_speed = 1.0_dp
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, &
                        'shifted Zhao photoelectron VDF must be not_applicable')
  call test_end()

  call test_begin('runtime adapter enforces the Zhao cold-ion approximation')
  app%particle_species(3)%normal_drift_speed = 0.0_dp
  app%particle_species(2)%temperature_k = 0.2_dp*3.0_dp*qe/k_boltzmann
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, &
                        'warm-ion Zhao input must be not_applicable')
  call test_end()

  call test_begin('runtime adapter rejects duplicate Zhao populations')
  app%particle_species(2)%temperature_k = 500.0_dp
  app%n_particle_species = 4_i32
  app%particle_species(4) = app%particle_species(2)
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, &
                        'duplicate Zhao ion population must be not_applicable')
  call test_end()

  call test_begin('runtime adapter rejects full-vector Zhao drift mode')
  app%n_particle_species = 3_i32
  app%sim%sheath_electron_drift_mode = 'full'
  call resolve_kinetic_outer_options(app, -0.5_dp, options, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, &
                        'full-vector Zhao drift mode must be not_applicable')
  call test_end()

  call test_summary()

contains

  subroutine init_fixture(cfg)
    type(app_config), intent(out) :: cfg

    call default_app_config(cfg)
    cfg%outer_plasma%model = 'kinetic_1d'
    cfg%outer_plasma%debye_length = 2.0_dp
    cfg%outer_plasma%thermal_voltage = 3.0_dp
    cfg%n_particle_species = 2_i32
    cfg%particle_species(1)%enabled = .true.
    cfg%particle_species(1)%source_mode = 'reservoir_face'
    cfg%particle_species(1)%inject_face = 'z_high'
    cfg%particle_species(1)%q_particle = -qe
    cfg%particle_species(1)%m_particle = 9.1093837139e-31_dp
    cfg%particle_species(1)%number_density_m3 = 2.0e6_dp
    cfg%particle_species(1)%has_number_density_m3 = .true.
    cfg%particle_species(1)%temperature_ev = 3.0_dp
    cfg%particle_species(1)%has_temperature_ev = .true.
    cfg%particle_species(2)%enabled = .true.
    cfg%particle_species(2)%source_mode = 'reservoir_face'
    cfg%particle_species(2)%inject_face = 'z_high'
    cfg%particle_species(2)%q_particle = qe
    cfg%particle_species(2)%m_particle = 1836.0_dp*9.1093837139e-31_dp
    cfg%particle_species(2)%number_density_m3 = 2.0e6_dp
    cfg%particle_species(2)%has_number_density_m3 = .true.
    cfg%particle_species(2)%temperature_k = 500.0_dp
    cfg%particle_species(2)%has_temperature_k = .true.
    cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e4_dp]
  end subroutine init_fixture

end program test_outer_plasma_kinetic_runtime
