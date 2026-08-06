!> 自動表面電流modelのdispatchとZhao電流分解を検証する。
program test_surface_current_model
  use bem_kinds, only: dp, i32
  use bem_app_config, only: app_config, default_app_config
  use bem_app_config_types, only: particle_species_spec
  use bem_surface_current_model, only: surface_current_model_result_type, evaluate_surface_current_model
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  type(app_config) :: cfg
  type(surface_current_model_result_type) :: result
  real(dp) :: inward_speed

  call test_init(4)

  call test_begin('none_dispatch')
  call default_app_config(cfg)
  call evaluate_surface_current_model(cfg, result)
  call assert_true(.not. result%active, 'none current model must stay inactive')
  call test_end()

  call test_begin('zhao_stationary_type_b_channels')
  call configure_zhao_fixture(cfg)
  cfg%surface_current%zhao_branch = 'b'
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%zhao_branch == 'B', 'explicit Zhao Type B was not selected')
  call assert_current_decomposition(result, 'Type B')
  call test_end()

  call test_begin('zhao_stationary_type_c_channels')
  call configure_zhao_fixture(cfg, alpha_deg=10.0_dp)
  cfg%surface_current%zhao_branch = 'c'
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%zhao_branch == 'C', 'explicit Zhao Type C was not selected')
  call assert_current_decomposition(result, 'Type C')
  call test_end()

  call test_begin('zhao_stationary_type_a_channels')
  call configure_zhao_fixture(cfg)
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%active, 'Zhao current model must be active')
  call assert_true(result%zhao_branch == 'A', 'alpha=60 Zhao model must select Type A')
  call assert_close_dp(result%phi0_v, 2.9712182827319435_dp, 5.0e-6_dp, 'Zhao phi0 mismatch')
  call assert_close_dp(result%phi_m_v, -0.8169121871620854_dp, 5.0e-6_dp, 'Zhao phi_m mismatch')
  call assert_close_dp( &
    result%ambient_electron_density_m3, 7.819215729579456e6_dp, 5.0e1_dp, &
    'Zhao ambient electron density mismatch' &
    )
  call assert_close_dp( &
    result%photoelectron_escape_current_density_a_m2, 3.9386846806723257e-7_dp, 5.0e-13_dp, &
    'Zhao PE escape current mismatch' &
    )
  call assert_close_dp( &
    result%photoelectron_emission_current_density_a_m2 + &
    result%photoelectron_return_current_density_a_m2, &
    result%photoelectron_escape_current_density_a_m2, 1.0e-18_dp, &
    'PE emission/return/escape identity mismatch' &
    )
  call assert_close_dp( &
    result%net_current_density_a_m2, 0.0_dp, 1.0e-12_dp, &
    'stationary Zhao current must close to zero' &
    )
  call assert_close_dp( &
    result%absorbed_current_a(3) + result%emission_current_a(3), &
    2.0_dp*result%photoelectron_escape_current_density_a_m2, 1.0e-18_dp, &
    'PE total target current must use the configured reference area' &
    )
  call test_end()

  call test_summary()

contains

  subroutine configure_zhao_fixture(app, alpha_deg)
    type(app_config), intent(inout) :: app
    real(dp), intent(in), optional :: alpha_deg
    real(dp) :: resolved_alpha

    resolved_alpha = 60.0_dp
    if (present(alpha_deg)) resolved_alpha = alpha_deg
    inward_speed = 4.68e5_dp*sin(resolved_alpha*acos(-1.0_dp)/180.0_dp)
    if (.not. allocated(app%particle_species)) allocate (app%particle_species(3))
    app%n_particle_species = 3_i32
    app%particle_species(1:3) = particle_species_spec()
    app%particle_species(1:3)%enabled = .true.
    app%particle_species(1)%species_key = 'electron'
    app%particle_species(1)%q_particle = -1.602176634e-19_dp
    app%particle_species(1)%m_particle = 9.1093837015e-31_dp
    app%particle_species(1)%temperature_ev = 12.0_dp
    app%particle_species(1)%has_temperature_ev = .true.
    app%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -inward_speed]
    app%particle_species(2)%species_key = 'ion'
    app%particle_species(2)%q_particle = 1.602176634e-19_dp
    app%particle_species(2)%m_particle = 1.67262192369e-27_dp
    app%particle_species(2)%number_density_m3 = 8.7e6_dp
    app%particle_species(2)%has_number_density_m3 = .true.
    app%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -inward_speed]
    app%particle_species(3)%species_key = 'photoelectron'
    app%particle_species(3)%q_particle = -1.602176634e-19_dp
    app%particle_species(3)%m_particle = 9.1093837015e-31_dp
    app%particle_species(3)%temperature_ev = 2.2_dp
    app%particle_species(3)%has_temperature_ev = .true.
    app%surface_current%model = 'zhao_stationary'
    app%surface_current%zhao_branch = 'auto'
    app%surface_current%electron_species = 'electron'
    app%surface_current%ion_species = 'ion'
    app%surface_current%photoelectron_species = 'photoelectron'
    app%surface_current%solar_elevation_deg = resolved_alpha
    app%surface_current%photoelectron_ref_density_m3 = 64.0e6_dp
    app%surface_current%reference_area_m2 = 2.0_dp
    app%surface_current%has_reference_area_m2 = .true.
  end subroutine configure_zhao_fixture

  subroutine assert_current_decomposition(current, label)
    type(surface_current_model_result_type), intent(in) :: current
    character(len=*), intent(in) :: label

    call assert_true(current%electron_current_density_a_m2 < 0.0_dp, trim(label)//' electron current sign mismatch')
    call assert_true(current%ion_current_density_a_m2 > 0.0_dp, trim(label)//' ion current sign mismatch')
    call assert_true( &
      current%photoelectron_emission_current_density_a_m2 > 0.0_dp .and. &
      current%photoelectron_escape_current_density_a_m2 >= 0.0_dp .and. &
      current%photoelectron_return_current_density_a_m2 <= 0.0_dp, &
      trim(label)//' PE current signs mismatch' &
      )
    call assert_close_dp( &
      current%photoelectron_emission_current_density_a_m2 + &
      current%photoelectron_return_current_density_a_m2, &
      current%photoelectron_escape_current_density_a_m2, 1.0e-18_dp, &
      trim(label)//' PE emission/return/escape identity mismatch' &
      )
    call assert_close_dp( &
      current%net_current_density_a_m2, 0.0_dp, 1.0e-12_dp, &
      trim(label)//' stationary net current mismatch' &
      )
  end subroutine assert_current_decomposition

end program test_surface_current_model
