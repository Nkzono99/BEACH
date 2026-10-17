!> 自動表面電流modelのdispatchとZhao電流分解を検証する。
program test_surface_current_model
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_types, only: bc_open, bc_periodic
  use bem_app_config, only: app_config, default_app_config
  use bem_app_config_types, only: particle_species_spec
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_surface_current_model, only: &
    surface_current_model_result_type, evaluate_surface_current_model, evaluate_surface_closure, &
    solve_zhao_outflow_closure
  use bem_outflow_spectrum, only: outflow_spectrum_bins, outflow_spectrum_edges_ev, outflow_spectrum_bin
  use sheath_model, only: log_energy_edges
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, &
                          assert_equal_i32
  implicit none

  type(app_config) :: cfg
  type(surface_current_model_result_type) :: result
  type(surface_closure_contract_type) :: closure
  type(surface_current_model_result_type) :: refreshed
  character(len=256) :: message
  real(dp) :: barrier_v, expected_escape, cold_phi_m
  real(dp) :: edges(outflow_spectrum_bins + 1), spectrum(outflow_spectrum_bins), tail_spectrum(outflow_spectrum_bins)
  type(surface_current_model_result_type) :: maxwell_root
  integer(i32) :: j
  logical :: success

  call test_init(14)

  call test_begin('none_dispatch')
  call default_app_config(cfg)
  call evaluate_surface_current_model(cfg, result)
  call assert_true(.not. result%active, 'none current model must stay inactive')
  call evaluate_surface_closure(cfg, closure)
  call assert_true(.not. closure%active, 'none surface closure must stay inactive')
  call assert_true(allocated(closure%has_inflow_number_flux), 'surface closure flux flags must be allocated')
  call assert_true(allocated(closure%inflow_number_flux_m2_s), 'surface closure flux storage must be allocated')
  if (allocated(closure%has_inflow_number_flux)) then
    call assert_true( &
      size(closure%has_inflow_number_flux) == cfg%n_particle_species, &
      'surface closure flux flag capacity mismatch' &
      )
  end if
  if (allocated(closure%inflow_number_flux_m2_s)) then
    call assert_true( &
      size(closure%inflow_number_flux_m2_s) == cfg%n_particle_species, &
      'surface closure flux storage capacity mismatch' &
      )
  end if
  call test_end()

  call test_begin('zhao_stationary_type_b_channels')
  call configure_zhao_fixture(cfg)
  cfg%surface_current%zhao_branch = 'b'
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%zhao_branch == 'B', 'explicit Zhao Type B was not selected')
  call assert_current_decomposition(result, 'Type B')
  call assert_true(ieee_is_finite(result%top_field_v_m) .and. result%top_field_v_m > 0.0_dp, &
                   'Type B top field sign and finiteness')
  call assert_kinetic_contract(result, 0.0_dp, 'Type B')
  call test_end()

  call test_begin('zhao_stationary_type_c_channels')
  call configure_zhao_fixture(cfg, alpha_deg=10.0_dp)
  cfg%surface_current%zhao_branch = 'c'
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%zhao_branch == 'C', 'explicit Zhao Type C was not selected')
  call assert_current_decomposition(result, 'Type C')
  call assert_true(ieee_is_finite(result%top_field_v_m) .and. result%top_field_v_m < 0.0_dp, &
                   'Type C top field sign and finiteness')
  call assert_kinetic_contract(result, 0.0_dp, 'Type C')
  call test_end()

  call test_begin('zhao_stationary_type_a_known_root_and_channels')
  call configure_zhao_fixture(cfg)
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%zhao_branch == 'A', 'alpha=60 Zhao model must select Type A')
  ! sheath-model の J=0 根（太陽高度60度、無drift電子、冷たいイオン）。
  call assert_close_dp(result%phi0_v, 3.8401890937376293_dp, 5.0e-6_dp, 'Zhao phi0 mismatch')
  call assert_close_dp(result%phi_m_v, -0.33165814284053974_dp, 5.0e-6_dp, 'Zhao phi_m mismatch')
  call assert_close_dp( &
    result%ambient_electron_density_m3, 9.916982425098775e6_dp, 5.0e1_dp, &
    'Zhao ambient electron density mismatch' &
    )
  call assert_close_dp( &
    result%photoelectron_escape_current_density_a_m2, 3.308285627995705e-7_dp, 5.0e-13_dp, &
    'Zhao PE escape current mismatch' &
    )
  call assert_current_decomposition(result, 'Type A')
  call assert_true(ieee_is_finite(result%top_field_v_m) .and. result%top_field_v_m > 0.0_dp, &
                   'Type A top field sign and finiteness')
  call assert_close_dp(result%reference_area_m2, 2.0_dp, 0.0_dp, 'configured reference area mismatch')
  call assert_close_dp( &
    result%escaped_particle_current_a(3), &
    -2.0_dp*result%photoelectron_escape_current_density_a_m2, 1.0e-18_dp, &
    'Zhao PE escape target did not use the configured reference area' &
    )
  call assert_kinetic_contract(result, result%phi_m_v, 'Type A')
  call test_end()

  call test_begin('zhao_warm_ions_use_species_temperature_and_pressure_factor')
  ! 太陽風では入口エネルギー (約 1 keV) が T_i よりはるかに大きく、流体圧力項は根をほとんど動かさない。
  ! それでも温度と係数は sheath-model へ渡り、収支は閉じたままでなければならない。
  call configure_zhao_fixture(cfg)
  call evaluate_surface_current_model(cfg, result)
  cold_phi_m = result%phi_m_v
  cfg%particle_species(2)%temperature_ev = 10.0_dp
  cfg%surface_current%zhao_ion_pressure_factor = 3.0_dp
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%active .and. result%zhao_branch == 'A', 'warm-ion Zhao root must stay Type A')
  call assert_true(result%phi_m_v /= cold_phi_m, 'warm ions must enter the outer root')
  call assert_close_dp(result%phi_m_v, cold_phi_m, 1.0e-3_dp, 'warm ions at 1 keV entry must shift phi_m by < 1 mV')
  call assert_current_decomposition(result, 'warm-ion Type A')
  call assert_kinetic_contract(result, result%phi_m_v, 'warm-ion Type A')
  call test_end()

  call test_begin('zhao_stationary_without_photoelectrons')
  call configure_zhao_fixture(cfg)
  cfg%n_particle_species = 2_i32
  cfg%surface_current%photoelectron_species = ''
  cfg%surface_current%solar_elevation_deg = 0.0_dp
  cfg%surface_current%photoelectron_ref_density_m3 = 0.0_dp
  cfg%surface_current%photoelectron_source_scale = 0.0_dp
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%active, 'no-PE Zhao current model must be active')
  call assert_true(.not. result%photoelectron_active, 'no-PE Zhao result must disable photoelectron channels')
  call assert_true(result%photoelectron_species_idx == 0_i32, 'no-PE Zhao result must not resolve a PE species')
  call assert_true(result%zhao_branch == 'C', 'no-PE Zhao model must select Type C')
  call assert_true(all(result%has_absorbed_target), 'no-PE Zhao must target electron and ion absorption')
  call assert_true( &
    result%absorbed_current_a(1) < 0.0_dp .and. result%absorbed_current_a(2) > 0.0_dp, &
    'no-PE Zhao ambient target signs mismatch' &
    )
  call assert_true(.not. any(result%has_emission_target), 'no-PE Zhao must not create emission targets')
  call assert_true(.not. any(result%has_escape_target), 'no-PE Zhao must not create PE escape targets')
  call assert_close_dp(result%photoelectron_emission_current_density_a_m2, 0.0_dp, 0.0_dp, 'no-PE emission')
  call assert_close_dp(result%photoelectron_escape_current_density_a_m2, 0.0_dp, 0.0_dp, 'no-PE escape')
  call assert_close_dp(result%photoelectron_return_current_density_a_m2, 0.0_dp, 0.0_dp, 'no-PE return')
  call assert_close_dp(result%net_current_density_a_m2, 0.0_dp, 1.0e-12_dp, 'no-PE stationary net current')
  call assert_kinetic_contract(result, 0.0_dp, 'no-PE Zhao')
  call test_end()

  call test_begin('zhao_drifting_electrons_reject_type_a')
  ! 内向きdriftの電子で反射集団を持つA/Cは無限遠へ接続できない。autoは成立するBへ進む。
  call configure_zhao_fixture(cfg)
  cfg%particle_species(1)%drift_velocity = cfg%particle_species(2)%drift_velocity
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%zhao_branch == 'B', 'drifting electrons must fall back to the admissible Type B')
  call assert_close_dp(result%phi0_v, 4.058764424794152_dp, 5.0e-6_dp, 'drifting-electron Type B phi0 mismatch')
  call assert_current_decomposition(result, 'drifting-electron Type B')
  cfg%surface_current%zhao_branch = 'a'
  call solve_zhao_outflow_closure( &
    cfg, result%outer_photoelectron_flux_m2_s, 2.2_dp, ' ', refreshed, success, message &
    )
  call assert_true(.not. success, 'drifting electrons with reflected slow electrons have no Type A profile')
  call assert_true(len_trim(message) > 0, 'rejected Type A must report a reason')
  ! Accepting a narrow negative-E^2 band next to upstream recovers the algebraic A root
  ! (reference solved independently in Python: phi0, phi_m, and a 42 mV band).
  cfg%surface_current%zhao_upstream_band_tolerance = 0.05_dp
  call solve_zhao_outflow_closure( &
    cfg, result%outer_photoelectron_flux_m2_s, 2.2_dp, ' ', refreshed, success, message &
    )
  call assert_true(success, 'tolerated drifting Type A failed: '//trim(message))
  call assert_true(refreshed%zhao_branch == 'A', 'tolerated drifting root must stay Type A')
  call assert_close_dp(refreshed%phi0_v, 2.654426813939951_dp, 1.0e-6_dp, 'drifting Type A phi0 mismatch')
  call assert_close_dp(refreshed%phi_m_v, -1.189723811591865_dp, 1.0e-6_dp, 'drifting Type A phi_m mismatch')
  call assert_close_dp(refreshed%upstream_negative_band_v, 0.04203_dp, 2.0e-3_dp, 'drifting Type A band mismatch')
  call assert_current_decomposition(refreshed, 'tolerated drifting Type A')
  ! 場の下限でも同じ根を受理する。帯の内側の場の不足は mV/m の桁で、太陽風の運動電場 |u x B| と同程度。
  cfg%surface_current%zhao_upstream_band_tolerance = 0.0_dp
  cfg%surface_current%zhao_upstream_field_floor_v_m = 0.1_dp
  call solve_zhao_outflow_closure( &
    cfg, result%outer_photoelectron_flux_m2_s, 2.2_dp, ' ', refreshed, success, message &
    )
  call assert_true(success, 'field-floor drifting Type A failed: '//trim(message))
  call assert_close_dp(refreshed%phi0_v, 2.654426813939951_dp, 1.0e-6_dp, &
                       'field-floor Type A phi0 must match the width-tolerated root')
  call assert_true(refreshed%upstream_field_deficit_v_m > 1.0e-4_dp .and. refreshed%upstream_field_deficit_v_m < 0.1_dp, &
                   'field deficit inside the upstream band must be of order mV/m')
  cfg%surface_current%zhao_upstream_field_floor_v_m = 1.0e-5_dp
  call solve_zhao_outflow_closure( &
    cfg, result%outer_photoelectron_flux_m2_s, 2.2_dp, ' ', refreshed, success, message &
    )
  call assert_true(.not. success, 'a field deficit above the floor must be rejected')
  cfg%surface_current%zhao_upstream_field_floor_v_m = 0.0_dp
  cfg%surface_current%zhao_upstream_band_tolerance = 0.02_dp
  call solve_zhao_outflow_closure( &
    cfg, result%outer_photoelectron_flux_m2_s, 2.2_dp, ' ', refreshed, success, message &
    )
  call assert_true(.not. success, 'a band wider than the tolerance must be rejected')
  call test_end()

  call test_begin('zhao_outflow_refresh_with_emission_source_reproduces_root')
  call configure_zhao_fixture(cfg)
  cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  call evaluate_surface_current_model(cfg, result)
  call solve_zhao_outflow_closure( &
    cfg, result%outer_photoelectron_flux_m2_s, result%outer_photoelectron_mean_energy_ev, ' ', refreshed, &
    success, message &
    )
  call assert_true(success, 'emission-source outer solve failed: '//trim(message))
  call assert_true(refreshed%zhao_branch == result%zhao_branch, 'emission-source outer branch mismatch')
  call assert_close_dp(refreshed%phi0_v, result%phi0_v, 1.0e-8_dp, 'emission-source phi0 mismatch')
  call assert_close_dp(refreshed%phi_m_v, result%phi_m_v, 1.0e-8_dp, 'emission-source phi_m mismatch')
  call assert_close_dp( &
    refreshed%photoelectron_escape_current_density_a_m2, result%photoelectron_escape_current_density_a_m2, &
    1.0e-9_dp*result%photoelectron_escape_current_density_a_m2, 'emission-source PE escape mismatch' &
    )
  call assert_close_dp( &
    refreshed%photoelectron_emission_current_density_a_m2, result%photoelectron_emission_current_density_a_m2, &
    0.0_dp, 'outer solve changed the surface emission' &
    )
  call assert_current_decomposition(refreshed, 'emission-source outer solve')
  call assert_true(result%has_plane_gauge .and. refreshed%has_plane_gauge, 'Zhao closure must publish its H gauge')
  call assert_close_dp(result%plane_gauge_potential_v, result%phi0_v, 0.0_dp, 'H gauge must be the wall potential')
  call assert_true(result%outer_return_cell_uniform, 'x/y periodic Zhao closure must return particles cell-uniformly')
  call test_end()

  call test_begin('zhao_outflow_refresh_uses_observed_outflow_as_outer_source')
  call configure_zhao_fixture(cfg)
  call evaluate_surface_current_model(cfg, result)
  call assert_true(.not. result%outer_return_cell_uniform, 'open x/y Zhao closure must keep crossing-point returns')
  call solve_zhao_outflow_closure( &
    cfg, 0.3_dp*result%outer_photoelectron_flux_m2_s, 4.0_dp, result%zhao_branch, refreshed, success, message &
    )
  call assert_true(success, 'observed-outflow outer solve failed: '//trim(message))
  select case (refreshed%zhao_branch)
  case ('A')
    barrier_v = refreshed%phi0_v - refreshed%phi_m_v
  case ('B')
    barrier_v = refreshed%phi0_v
  case default
    barrier_v = 0.0_dp
  end select
  ! Half-Maxwellian flux above the outer barrier, independent of the Zhao root equations.
  expected_escape = qe*0.3_dp*result%outer_photoelectron_flux_m2_s*exp(-barrier_v/4.0_dp)
  call assert_close_dp( &
    refreshed%photoelectron_escape_current_density_a_m2, expected_escape, 1.0e-9_dp*expected_escape, &
    'outer PE escape must be the observed source above the barrier' &
    )
  call assert_close_dp( &
    refreshed%photoelectron_emission_current_density_a_m2, result%photoelectron_emission_current_density_a_m2, &
    0.0_dp, 'observed outflow must not replace the surface emission' &
    )
  call assert_close_dp( &
    refreshed%outer_photoelectron_flux_m2_s, 0.3_dp*result%outer_photoelectron_flux_m2_s, 0.0_dp, &
    'outer source flux receipt mismatch' &
    )
  call assert_close_dp( &
    refreshed%net_current_density_a_m2, 0.0_dp, 1.0e-6_dp*refreshed%ion_current_density_a_m2, &
    'refreshed outer root must float' &
    )
  call assert_current_decomposition(refreshed, 'observed-outflow outer solve')
  call assert_close_dp(refreshed%plane_gauge_potential_v, refreshed%phi0_v, 0.0_dp, 'refreshed H gauge mismatch')
  call test_end()

  call test_begin('outflow_spectrum_bins_cover_energy_axis')
  edges = outflow_spectrum_edges_ev()
  call assert_close_dp(edges(1), 0.0_dp, 0.0_dp, 'first spectrum edge must be 0 eV')
  call assert_close_dp(edges(2), 1.0e-3_dp, 1.0e-15_dp, 'second spectrum edge must be 1 meV')
  call assert_close_dp(edges(outflow_spectrum_bins + 1), 1.0e3_dp, 1.0e-9_dp, 'last spectrum edge must be 1 keV')
  call assert_true(all(edges(2:) > edges(:outflow_spectrum_bins)), 'spectrum edges must increase')
  ! The observed spectrum becomes a sheath-model source on the library's own outflow grid.
  call assert_true(all(edges == log_energy_edges()), 'spectrum edges must equal the sheath-model outflow grid')
  do j = 2_i32, outflow_spectrum_bins
    call assert_equal_i32( &
      outflow_spectrum_bin(sqrt(edges(j)*edges(j + 1))), j, 'bin centre must fall in its own bin' &
      )
  end do
  call assert_equal_i32(outflow_spectrum_bin(0.0_dp), 1_i32, 'zero energy must fall in the first bin')
  call assert_equal_i32(outflow_spectrum_bin(5.0e-4_dp), 1_i32, 'sub-meV energy must fall in the first bin')
  call assert_equal_i32(outflow_spectrum_bin(-1.0_dp), 1_i32, 'negative energy must fall in the first bin')
  call assert_equal_i32( &
    outflow_spectrum_bin(1.0e5_dp), outflow_spectrum_bins, 'energy above 1 keV must fall in the last bin' &
    )
  call test_end()

  call test_begin('zhao_outflow_refresh_binned_maxwellian_matches_analytic_source')
  ! 半Maxwell源をbinへ積分した外部源は、解析的な半Maxwell源と同じ根を与える（bin内一定の離散化誤差だけ）。
  call configure_zhao_fixture(cfg)
  call evaluate_surface_current_model(cfg, result)
  call solve_zhao_outflow_closure( &
    cfg, 0.3_dp*result%outer_photoelectron_flux_m2_s, 4.0_dp, result%zhao_branch, maxwell_root, success, message &
    )
  call assert_true(success, 'analytic outer solve failed: '//trim(message))
  call bin_shifted_half_maxwellian(0.3_dp*result%outer_photoelectron_flux_m2_s, 4.0_dp, 0.0_dp, spectrum)
  call solve_zhao_outflow_closure( &
    cfg, 0.3_dp*result%outer_photoelectron_flux_m2_s, 4.0_dp, result%zhao_branch, refreshed, success, message, &
    outer_spectrum_m2_s=spectrum &
    )
  call assert_true(success, 'binned outer solve failed: '//trim(message))
  call assert_true(refreshed%zhao_branch == maxwell_root%zhao_branch, 'binned outer branch mismatch')
  call assert_close_dp(refreshed%phi0_v, maxwell_root%phi0_v, 0.02_dp, 'binned phi0 mismatch')
  call assert_close_dp(refreshed%phi_m_v, maxwell_root%phi_m_v, 0.02_dp, 'binned phi_m mismatch')
  call assert_close_dp( &
    refreshed%photoelectron_escape_current_density_a_m2, maxwell_root%photoelectron_escape_current_density_a_m2, &
    0.01_dp*maxwell_root%photoelectron_escape_current_density_a_m2, 'binned PE escape mismatch' &
    )
  call assert_true( &
    all(refreshed%outer_photoelectron_spectrum_m2_s == spectrum), 'outer spectrum receipt must keep the input' &
    )
  call assert_true( &
    all(maxwell_root%outer_photoelectron_spectrum_m2_s == 0.0_dp), 'analytic source must not publish a spectrum' &
    )
  call assert_current_decomposition(refreshed, 'binned outer solve')
  call test_end()

  call test_begin('zhao_outflow_refresh_keeps_outflow_tail_shape')
  ! 三層レゴリスの流出に似た、2.2 eV の本体（85%）と加速された 40 eV の裾（15%）。平均 7.87 eV は裾で決まり、
  ! 同じ平均の半Maxwellへ縮約すると障壁を越える電子を過大に数えて A の根を失う。スペクトルのままなら A の根がある。
  call configure_zhao_fixture(cfg)
  cfg%particle_species(1)%temperature_ev = 10.0_dp
  cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
  cfg%particle_species(2)%number_density_m3 = 5.0e6_dp
  cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
  cfg%surface_current%zhao_branch = 'a'
  cfg%surface_current%zhao_upstream_band_tolerance = 0.3_dp
  call solve_zhao_outflow_closure(cfg, 2.84e12_dp, 7.87_dp, 'A', maxwell_root, success, message)
  call assert_true(.not. success, 'the mean-energy half-Maxwellian proxy unexpectedly kept a Type A root')
  call bin_shifted_half_maxwellian(0.85_dp*2.84e12_dp, 2.2_dp, 0.0_dp, spectrum)
  call bin_shifted_half_maxwellian(0.15_dp*2.84e12_dp, 40.0_dp, 0.0_dp, tail_spectrum)
  spectrum = spectrum + tail_spectrum
  call solve_zhao_outflow_closure( &
    cfg, 2.84e12_dp, 7.87_dp, 'A', refreshed, success, message, outer_spectrum_m2_s=spectrum &
    )
  call assert_true(success, 'body-and-tail spectrum lost its Type A root: '//trim(message))
  call assert_true(refreshed%zhao_branch == 'A', 'body-and-tail root must stay Type A')
  ! sheath-model を直接使った同じ源の根（φ0 = 1.376 V、上流の帯 41 mV）。
  call assert_close_dp(refreshed%phi0_v, 1.376_dp, 5.0e-3_dp, 'body-and-tail phi0 mismatch')
  barrier_v = refreshed%phi0_v - refreshed%phi_m_v
  ! bin内で一定の dGamma/dK のまま障壁より上を数える。Zhaoの根の式とは独立に、源の分割だけを検証する。
  expected_escape = 0.0_dp
  edges = outflow_spectrum_edges_ev()
  do j = 1_i32, outflow_spectrum_bins
    if (edges(j) >= barrier_v) then
      expected_escape = expected_escape + spectrum(j)
    else if (edges(j + 1) > barrier_v) then
      expected_escape = expected_escape + spectrum(j)*(edges(j + 1) - barrier_v)/(edges(j + 1) - edges(j))
    end if
  end do
  expected_escape = qe*expected_escape
  call assert_close_dp( &
    refreshed%photoelectron_escape_current_density_a_m2, expected_escape, 1.0e-6_dp*expected_escape, &
    'outer PE escape must be the binned source above the barrier' &
    )
  call assert_current_decomposition(refreshed, 'body-and-tail outer solve')
  call test_end()

  call test_begin('zhao_outflow_refresh_rejects_empty_source')
  call configure_zhao_fixture(cfg)
  call solve_zhao_outflow_closure(cfg, 0.0_dp, 2.0_dp, 'A', refreshed, success, message)
  call assert_true(.not. success, 'empty observed outflow must not produce an outer root')
  call assert_true(len_trim(message) > 0, 'empty observed outflow must report a reason')
  spectrum = 0.0_dp
  call solve_zhao_outflow_closure( &
    cfg, 1.0e12_dp, 2.0_dp, 'A', refreshed, success, message, outer_spectrum_m2_s=spectrum &
    )
  call assert_true(.not. success, 'an empty outer spectrum must not produce an outer root')
  spectrum(10) = -1.0_dp
  spectrum(20) = 2.0e12_dp
  call solve_zhao_outflow_closure( &
    cfg, 1.0e12_dp, 2.0_dp, 'A', refreshed, success, message, outer_spectrum_m2_s=spectrum &
    )
  call assert_true(.not. success, 'a negative outer spectrum bin must be rejected')
  call test_end()

  call test_begin('secondary_electron_barrier_role_and_phase_one_pe_target_split')
  call configure_zhao_fixture(cfg)
  cfg%n_particle_species = 4_i32
  cfg%particle_species(4) = particle_species_spec()
  cfg%particle_species(4)%enabled = .true.
  cfg%particle_species(4)%species_key = 'secondary'
  cfg%particle_species(4)%source_mode = 'secondary_emission'
  cfg%surface_current%secondary_electron_species = 'secondary'
  cfg%sim%bc_low(1:2) = bc_periodic
  cfg%sim%bc_high(1:2) = bc_periodic
  call evaluate_surface_current_model(cfg, result)
  call assert_true(result%has_outflow_kinetic_barrier(4), 'SEE must have the PE kinetic barrier role')
  call assert_close_dp(result%outflow_barrier_potential_v(4), result%outflow_barrier_potential_v(3), &
                       0.0_dp, 'SEE and PE barrier potentials must agree')
  call assert_equal_i32(result%outflow_barrier_face(4), 6_i32, 'SEE barrier face')
  call assert_true(result%outer_return_cell_uniform, 'SEE uses the periodic-cell outer return map')
  call assert_true(.not. result%has_absorbed_target(4) .and. .not. result%has_emission_target(4) .and. &
                   .not. result%has_escape_target(4), 'SEE must stay fully tracked')
  call solve_zhao_outflow_closure(cfg, result%outer_photoelectron_flux_m2_s, 2.2_dp, 'A', &
                                  refreshed, success, message, secondary_share=0.25_dp)
  call assert_true(success, 'mixed emission source solve: '//trim(message))
  call assert_close_dp(refreshed%escaped_particle_current_a(3), &
                       -0.75_dp*refreshed%photoelectron_escape_current_density_a_m2*refreshed%reference_area_m2, &
                       1.0e-20_dp, 'PE escape target must use its outflow share')
  call assert_close_dp(refreshed%absorbed_current_a(3), &
                       -refreshed%escaped_particle_current_a(3) - refreshed%emission_current_a(3), &
                       1.0e-20_dp, 'PE return target must remain escape minus emission')
  call test_end()

  call test_summary()

contains

  !> 法線エネルギーで K_c だけずらした半Maxwell流束 dGamma/dK ~ exp(-(K-K_c)/T) を出力binへ厳密に積分する。
  subroutine bin_shifted_half_maxwellian(total_flux, temperature_ev, shift_ev, bins)
    real(dp), intent(in) :: total_flux, temperature_ev, shift_ev
    real(dp), intent(out) :: bins(outflow_spectrum_bins)
    real(dp) :: bin_edges(outflow_spectrum_bins + 1), upper_tail
    integer(i32) :: k

    bin_edges = outflow_spectrum_edges_ev()
    do k = 1_i32, outflow_spectrum_bins
      ! 最後のbinは1 keV以上の裾も数える。
      upper_tail = 0.0_dp
      if (k < outflow_spectrum_bins) upper_tail = exp(-max(bin_edges(k + 1) - shift_ev, 0.0_dp)/temperature_ev)
      bins(k) = total_flux*(exp(-max(bin_edges(k) - shift_ev, 0.0_dp)/temperature_ev) - upper_tail)
    end do
  end subroutine bin_shifted_half_maxwellian

  subroutine configure_zhao_fixture(app, alpha_deg)
    type(app_config), intent(out) :: app
    real(dp), intent(in), optional :: alpha_deg
    real(dp) :: resolved_alpha, inward_speed

    call default_app_config(app)
    resolved_alpha = 60.0_dp
    if (present(alpha_deg)) resolved_alpha = alpha_deg
    inward_speed = 4.68e5_dp*sin(resolved_alpha*acos(-1.0_dp)/180.0_dp)
    app%n_particle_species = 3_i32
    app%particle_species(1:3) = particle_species_spec()
    app%particle_species(1:3)%enabled = .true.
    app%particle_species(1)%species_key = 'electron'
    app%particle_species(1)%q_particle = -1.602176634e-19_dp
    app%particle_species(1)%m_particle = 9.1093837015e-31_dp
    app%particle_species(1)%temperature_ev = 12.0_dp
    app%particle_species(1)%has_temperature_ev = .true.
    app%particle_species(1)%drift_velocity = 0.0_dp
    app%particle_species(2)%species_key = 'ion'
    app%particle_species(2)%q_particle = 1.602176634e-19_dp
    app%particle_species(2)%m_particle = 1.67262192369e-27_dp
    app%particle_species(2)%number_density_m3 = 8.7e6_dp
    app%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -inward_speed]
    ! 既知の根は冷たいイオンで解いたもの。温かいイオンは別のテストで調べる。
    app%particle_species(2)%temperature_ev = 0.0_dp
    app%particle_species(2)%has_temperature_ev = .true.
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

    call assert_true(current%active .and. current%photoelectron_active, trim(label)//' model activity mismatch')
    call assert_true(current%electron_current_density_a_m2 < 0.0_dp, trim(label)//' electron current sign mismatch')
    call assert_true(current%ion_current_density_a_m2 > 0.0_dp, trim(label)//' ion current sign mismatch')
    call assert_true( &
      current%photoelectron_emission_current_density_a_m2 > 0.0_dp .and. &
      current%photoelectron_escape_current_density_a_m2 >= 0.0_dp .and. &
      current%photoelectron_return_current_density_a_m2 <= 0.0_dp, &
      trim(label)//' PE current signs mismatch' &
      )
    call assert_true( &
      all(current%has_absorbed_target) .and. &
      count(current%has_emission_target) == 1 .and. current%has_emission_target(3) .and. &
      count(current%has_escape_target) == 1 .and. current%has_escape_target(3), &
      trim(label)//' current target routing mismatch' &
      )
    call assert_close_dp( &
      current%photoelectron_emission_current_density_a_m2 + &
      current%photoelectron_return_current_density_a_m2, &
      current%photoelectron_escape_current_density_a_m2, 1.0e-18_dp, &
      trim(label)//' PE emission/return/escape identity mismatch' &
      )
    call assert_close_dp( &
      current%photoelectron_budget_residual_current_density_a_m2, 0.0_dp, 1.0e-18_dp, &
      trim(label)//' PE budget residual mismatch' &
      )
    call assert_close_dp( &
      current%surface_budget_residual_current_density_a_m2, 0.0_dp, 1.0e-12_dp, &
      trim(label)//' surface budget residual mismatch' &
      )
    call assert_close_dp( &
      current%net_current_density_a_m2, 0.0_dp, 1.0e-12_dp, &
      trim(label)//' stationary net current mismatch' &
      )
    call assert_true( &
      current%absorbed_current_a(1) < 0.0_dp .and. current%absorbed_current_a(2) > 0.0_dp .and. &
      current%absorbed_current_a(3) <= 0.0_dp .and. current%emission_current_a(3) > 0.0_dp .and. &
      current%escaped_particle_current_a(3) <= 0.0_dp, &
      trim(label)//' total-current channel signs mismatch' &
      )
  end subroutine assert_current_decomposition

  subroutine assert_kinetic_contract(current, electron_bottleneck_v, label)
    type(surface_current_model_result_type), intent(in) :: current
    real(dp), intent(in) :: electron_bottleneck_v
    character(len=*), intent(in) :: label

    call assert_true(trim(current%kinetic_contract) == 'zhao_barrier_v1', trim(label)//' kinetic contract mismatch')
    call assert_true( &
      count(current%has_inflow_kinetic_map) == 2 .and. all(current%has_inflow_kinetic_map(1:2)), &
      trim(label)//' ambient inflow map routing mismatch' &
      )
    call assert_true( &
      all(current%inflow_reservoir_potential_v(1:2) == 0.0_dp), &
      trim(label)//' ambient reservoir gauge mismatch' &
      )
    call assert_true( &
      all(current%inflow_kinetic_face(1:2) == 6_i32), &
      trim(label)//' inflow map face must be z-high' &
      )
    ! 電子源はZhaoの準中性から解いた上流Maxwellian密度、イオン源は設定の太陽風密度。
    call assert_close_dp( &
      current%inflow_reservoir_density_m3(1), current%ambient_electron_density_m3, &
      1.0e-12_dp*current%ambient_electron_density_m3, trim(label)//' electron reservoir density mismatch' &
      )
    call assert_close_dp( &
      current%inflow_reservoir_density_m3(2), 8.7e6_dp, 1.0e-6_dp, &
      trim(label)//' ion reservoir density mismatch' &
      )
    call assert_close_dp( &
      current%inflow_access_potential_v(1), electron_bottleneck_v, 1.0e-12_dp, &
      trim(label)//' electron access bottleneck mismatch' &
      )
    call assert_close_dp( &
      current%inflow_access_potential_v(2), 0.0_dp, 1.0e-12_dp, &
      trim(label)//' ion access bottleneck mismatch' &
      )
    call assert_true(all(current%has_outflow_kinetic_barrier), trim(label)//' outflow barrier routing mismatch')
    call assert_close_dp( &
      current%outflow_barrier_potential_v(1), electron_bottleneck_v, 1.0e-12_dp, &
      trim(label)//' electron barrier potential mismatch' &
      )
    call assert_close_dp( &
      current%outflow_barrier_potential_v(2), 0.0_dp, 1.0e-12_dp, &
      trim(label)//' ion barrier potential mismatch' &
      )
    if (current%photoelectron_active) then
      call assert_close_dp( &
        current%outflow_barrier_potential_v(3), electron_bottleneck_v, 1.0e-12_dp, &
        trim(label)//' PE barrier potential mismatch' &
        )
    end if
    call assert_true(all(current%outflow_barrier_face == 6_i32), trim(label)//' barrier face must be z-high')
  end subroutine assert_kinetic_contract

end program test_surface_current_model
