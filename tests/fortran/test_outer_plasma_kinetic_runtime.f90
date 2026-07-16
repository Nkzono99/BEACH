program test_outer_plasma_kinetic_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe, k_boltzmann
  use bem_app_config_types, only: app_config, default_app_config
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_not_applicable
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type
  use bem_outer_plasma_kinetic_runtime, only: resolve_kinetic_outer_options
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32
  implicit none

  type(app_config) :: app
  type(kinetic_outer_plasma_options_type) :: options
  integer(i32) :: status
  character(len=256) :: message

  call test_init(3)

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
