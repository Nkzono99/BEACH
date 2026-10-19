!> Flat periodic plateでzhao_stationaryの外部障壁ゲージと観測PE流出による弱連成を検証する。
program test_zhao_outflow_refresh
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats, injection_state, bc_open, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            particle_inflow_reservoir
  use bem_charge_ledger, only: charge_ledger_type, finite_charge_sum
  use bem_surface_current_model, only: surface_current_model_result_type, evaluate_surface_current_model
  use bem_zhao_outflow_refresh, only: zhao_outflow_refresh_type
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_mpi, only: mpi_context
  use bem_outflow_spectrum, only: outflow_spectrum_bins, outflow_spectrum_edges_ev
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          assert_close_dp, delete_file_if_exists
  implicit none

  real(dp), parameter :: cell_width = 1.0e-4_dp
  real(dp), parameter :: box_height = 1.0e-3_dp
  character(len=*), parameter :: fixed_history_path = 'test_zhao_outflow_refresh_fixed_history_tmp.csv'
  character(len=*), parameter :: spectrum_history_path = 'test_zhao_outflow_refresh_spectrum_history_tmp.csv'
  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats, resumed_stats
  type(injection_state) :: inject_state
  type(charge_ledger_type) :: ledger
  type(surface_current_model_result_type) :: static_root
  real(dp) :: emission_flux, predicted_escape_fraction, observed_escape_fraction
  real(dp) :: edges(outflow_spectrum_bins + 1), split_edge, slow_fraction, last_batch_flux
  real(dp) :: row_time, row_low, row_high, row_flux
  integer :: fixed_history_unit, read_ios, fixed_history_rows, spectrum_history_unit, split_bin
  integer(i32) :: row_batch, row_bin, last_row_batch
  logical :: spectrum_rows_valid
  character(len=512) :: line

  call test_init(5)

  call test_begin('stationary_barrier_uses_outer_wall_potential_gauge')
  call configure_fixture(mesh, cfg, inject_state, 0_i32)
  cfg%periodic2%lower_boundary_model = 'e_top_from_sheath'
  call evaluate_surface_current_model(cfg, static_root)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, stats, inject_state=inject_state, charge_ledger=ledger)
  call assert_true(.not. stats%matching_plane_state_valid, 'a fixed Zhao root must not publish refresh state')
  call assert_close_dp(stats%surface_current_top_field_v_m, static_root%top_field_v_m, &
                       1.0e-12_dp*abs(static_root%top_field_v_m), 'static root top-field receipt')
  call assert_close_dp(stats%matching_plane_top_cell_mean_v_m, static_root%top_field_v_m, &
                       1.0e-12_dp*abs(static_root%top_field_v_m), 'simulator passes sheath field to snapshot')
  ! A flat plate returns every PE that reaches H below phi0-phi_m. With the H gauge at phi0, the tracked
  ! escape fraction is the half-Maxwellian tail exp(-(phi0-phi_m)/T) of the same root.
  predicted_escape_fraction = static_root%photoelectron_escape_current_density_a_m2/ &
                              static_root%photoelectron_emission_current_density_a_m2
  observed_escape_fraction = ledger%escaped_to_infinity(3)/ledger%emitted_from_surface(3)
  call assert_close_dp( &
    observed_escape_fraction, predicted_escape_fraction, 0.03_dp, &
    'tracked PE escape fraction must follow the outer barrier above phi0' &
    )
  call assert_close_dp( &
    finite_charge_sum(mesh%q_elem, 'fixed-root surface charge'), 0.0_dp, &
    1.0e-6_dp*abs(ledger%emitted_from_surface(3)), 'zero-current targets must keep the cell floating' &
    )
  ! 上流電子をZhao根の密度で注入すれば、固定電流の倍率はMC雑音の範囲で1になる。
  call assert_close_dp( &
    ledger%fixed_absorbed_weight_scale(1), 1.0_dp, 0.15_dp, &
    'electron inflow must carry the Zhao upstream density' &
    )
  call test_end()

  call test_begin('flat_plate_outflow_refresh_reproduces_emission_source')
  call configure_fixture(mesh, cfg, inject_state, 1_i32)
  call seed_particles_from_config(cfg)
  open (newunit=fixed_history_unit, file=fixed_history_path, status='replace', action='write')
  open (newunit=spectrum_history_unit, file=spectrum_history_path, status='replace', action='write')
  call run_absorption_insulator( &
    mesh, cfg, stats, history_stride=1_i32, inject_state=inject_state, fixed_current_history_unit=fixed_history_unit, &
    matching_plane_spectrum_history_unit=spectrum_history_unit &
    )
  close (fixed_history_unit)
  close (spectrum_history_unit)
  open (newunit=fixed_history_unit, file=fixed_history_path, status='old', action='read')
  fixed_history_rows = 0
  do
    read (fixed_history_unit, '(a)', iostat=read_ios) line
    if (read_ios /= 0) exit
    fixed_history_rows = fixed_history_rows + 1
  end do
  close (fixed_history_unit)
  call delete_file_if_exists(fixed_history_path)
  ! electron/ion吸収、PE吸収・放出・脱出の5 channelを毎batch書く。
  call assert_equal_i32( &
    int(fixed_history_rows, i32), 5_i32*cfg%sim%batch_count, 'fixed-current history must write every active channel' &
    )
  emission_flux = static_root%photoelectron_emission_current_density_a_m2/qe
  call assert_true(stats%matching_plane_state_valid, 'outflow refresh state was not published')
  call assert_equal_i32( &
    stats%matching_plane_iterations, cfg%sim%batch_count + 1_i32, &
    'every one-batch window must refresh the outer root' &
    )
  call assert_close_dp( &
    stats%matching_plane_feedback(1)/emission_flux, 1.0_dp, 0.02_dp, &
    'every PE leaving a flat plate must cross H once' &
    )
  call assert_close_dp( &
    stats%matching_plane_feedback(2), 2.2_dp, 0.15_dp, &
    'flat-plate outflow mean normal energy must equal the emission temperature' &
    )
  call assert_close_dp(stats%matching_plane_phi_v, static_root%phi0_v, 0.3_dp, 'refreshed phi0 drifted')
  ! 平板から出たPEの法線エネルギーは放出の半Maxwell分布のまま: K < K_s の割合は 1 - exp(-K_s/T)。
  call assert_close_dp( &
    sum(stats%matching_plane_photoelectron_spectrum_m2_s), stats%matching_plane_feedback(1), &
    1.0e-12_dp*stats%matching_plane_feedback(1), 'outer spectrum must carry the outer source flux' &
    )
  edges = outflow_spectrum_edges_ev()
  split_bin = minloc(abs(edges - 2.2_dp), dim=1)
  split_edge = edges(split_bin)
  slow_fraction = sum(stats%matching_plane_photoelectron_spectrum_m2_s(:split_bin - 1))/ &
                  sum(stats%matching_plane_photoelectron_spectrum_m2_s)
  call assert_close_dp( &
    slow_fraction, 1.0_dp - exp(-split_edge/2.2_dp), 0.04_dp, 'flat-plate outflow spectrum must stay half-Maxwellian' &
    )
  ! スペクトル history は外部根を更新した各 batch で、正の bin だけを書く。
  open (newunit=spectrum_history_unit, file=spectrum_history_path, status='old', action='read')
  last_row_batch = 0_i32
  last_batch_flux = 0.0_dp
  spectrum_rows_valid = .true.
  do
    read (spectrum_history_unit, '(a)', iostat=read_ios) line
    if (read_ios /= 0) exit
    read (line, *, iostat=read_ios) row_batch, row_time, row_bin, row_low, row_high, row_flux
    spectrum_rows_valid = spectrum_rows_valid .and. read_ios == 0
    if (read_ios /= 0) cycle
    spectrum_rows_valid = spectrum_rows_valid .and. row_flux > 0.0_dp .and. row_batch >= last_row_batch .and. &
                          row_low == edges(row_bin) .and. row_high == edges(row_bin + 1)
    if (row_batch /= last_row_batch) last_batch_flux = 0.0_dp
    last_row_batch = row_batch
    last_batch_flux = last_batch_flux + row_flux
  end do
  close (spectrum_history_unit)
  call delete_file_if_exists(spectrum_history_path)
  call assert_true(spectrum_rows_valid, 'spectrum history rows must be positive bins in batch order')
  call assert_equal_i32(last_row_batch, cfg%sim%batch_count, 'spectrum history must end at the last refresh')
  call assert_close_dp( &
    last_batch_flux, stats%matching_plane_feedback(1), 1.0e-12_dp*stats%matching_plane_feedback(1), &
    'last spectrum history rows must carry the outer source flux' &
    )
  call assert_close_dp( &
    stats%matching_plane_photoelectron_return_flux_m2_s + stats%matching_plane_photoelectron_escape_flux_m2_s, &
    stats%matching_plane_feedback(1), 1.0e-12_dp*stats%matching_plane_feedback(1), &
    'outer source return/escape partition mismatch' &
    )
  call assert_close_dp( &
    stats%matching_plane_displacement_c_m2, 0.0_dp, 1.0e-6_dp*emission_flux*qe*cfg%sim%batch_duration, &
    'refreshed zero-current targets must keep the cell floating' &
    )
  call test_end()

  call test_begin('outflow_refresh_state_reconstructs_on_resume')
  cfg%sim%batch_count = stats%batches + 1_i32
  call run_absorption_insulator(mesh, cfg, resumed_stats, initial_stats=stats, inject_state=inject_state)
  call assert_true(resumed_stats%matching_plane_state_valid, 'resumed outflow refresh state was not valid')
  call assert_equal_i32( &
    resumed_stats%matching_plane_iterations, stats%matching_plane_iterations + 1_i32, &
    'resume must continue the accepted outer-solve count' &
    )
  call assert_close_dp(resumed_stats%matching_plane_phi_v, stats%matching_plane_phi_v, 0.3_dp, 'resumed phi0 jump')
  call assert_true( &
    any(resumed_stats%matching_plane_photoelectron_spectrum_m2_s > 0.0_dp), 'resume must keep refreshing the spectrum' &
    )
  call test_end()

  call test_begin('combined_pe_see_spectrum_uses_distinct_flux_normalizations')
  call check_emission_group_refresh(.false.)
  call test_end()
  call test_begin('secondary_only_refresh_can_move_from_c_to_a')
  call check_emission_group_refresh(.true.)
  call test_end()

  call test_summary()

contains

  subroutine check_emission_group_refresh(night)
    logical, intent(in) :: night
    type(mesh_type) :: plate
    type(app_config) :: app
    type(injection_state) :: state
    type(sim_stats) :: receipts, restored
    type(surface_current_model_result_type) :: initial
    type(surface_closure_contract_type) :: closure, restored_closure
    type(zhao_outflow_refresh_type) :: refresh, resumed_refresh
    type(mpi_context) :: mpi
    real(dp) :: moments(4 + outflow_spectrum_bins, 4, 1), raw_charge(12)
    real(dp) :: spectrum_pe(outflow_spectrum_bins), spectrum_see(outflow_spectrum_bins)
    real(dp) :: flux_pe, flux_see, raw_out, batch_time, area, see_temperature, root_escape_charge
    real(dp) :: closed_charge(4), unclosed_charge(4)
    integer :: batch, k
    logical :: refreshed

    call configure_fixture(plate, app, state, 2_i32)
    call evaluate_surface_current_model(app, initial)
    flux_pe = 0.5_dp*initial%photoelectron_emission_current_density_a_m2/qe
    flux_see = 0.25_dp*flux_pe
    app%n_particle_species = 4_i32
    app%surface_current%secondary_electron_species = 'secondary'
    app%particle_species(4) = species_from_defaults()
    app%particle_species(4)%species_key = 'secondary'
    app%particle_species(4)%source_mode = 'secondary_emission'
    app%particle_species(4)%q_particle = -qe
    app%particle_species(4)%m_particle = app%particle_species(1)%m_particle
    if (night) then
      flux_see = 2.0_dp*flux_pe
      flux_pe = 0.0_dp
      app%surface_current%photoelectron_source_scale = 0.0_dp
      app%surface_current%photoelectron_species = ''
      app%surface_current%solar_elevation_deg = 0.0_dp
      app%surface_current%photoelectron_ref_density_m3 = 0.0_dp
      app%particle_species(3)%enabled = .false.
    end if
    see_temperature = 1.1_dp
    if (night) see_temperature = 2.2_dp
    area = cell_width**2
    edges = outflow_spectrum_edges_ev()
    do k = 1, outflow_spectrum_bins
      spectrum_pe(k) = exp(-edges(k)/2.2_dp) - exp(-edges(k + 1)/2.2_dp)
      spectrum_see(k) = exp(-edges(k)/see_temperature) - exp(-edges(k + 1)/see_temperature)
    end do
    receipts = sim_stats()
    call refresh%initialize(app, plate, receipts, closure)
    if (night) then
      call assert_true(receipts%matching_plane_phi_v < 0.0_dp, 'night must start from a zero-source Type C root')
      restored = receipts
      call resumed_refresh%initialize(app, plate, restored, restored_closure)
      call assert_close_dp(restored%matching_plane_phi_v, receipts%matching_plane_phi_v, 1.0e-12_dp, &
                           'zero-source checkpoint before the first SEE refresh must be resumable')
    end if
    do batch = 1, 2
      batch_time = 1.0e-3_dp*real(2*batch - 1, dp)
      app%sim%batch_duration = batch_time
      moments = 0.0_dp
      raw_charge = 0.0_dp
      if (.not. night) then
        raw_out = 100.0_dp*real(batch, dp)
        moments(1, 3, 1) = raw_out
        moments(2, 3, 1) = raw_out*2.2_dp*qe
        moments(5:, 3, 1) = raw_out*spectrum_pe
        raw_charge(7) = 2.0_dp*raw_out*qe
      end if
      raw_out = flux_see*area*batch_time
      moments(1, 4, 1) = raw_out
      moments(2, 4, 1) = raw_out*see_temperature*qe
      moments(5:, 4, 1) = raw_out*spectrum_see
      closed_charge = 0.0_dp
      unclosed_charge = 0.0_dp
      call refresh%commit_batch(app, plate, mpi, moments, raw_charge, closed_charge, unclosed_charge, int(batch, i32), &
                                receipts, closure, refreshed)
      call assert_true(refreshed .eqv. (batch == 2), 'refresh must wait for the full accepted-batch window')
    end do
    call assert_close_dp(receipts%matching_plane_feedback(2), (2.2_dp*flux_pe + see_temperature*flux_see)/(flux_pe + flux_see), &
                         1.0e-12_dp, 'mean normal energy must be flux weighted')
    call assert_close_dp(receipts%matching_plane_feedback(1), flux_pe + flux_see, &
                         1.0e-12_dp*(flux_pe + flux_see), 'PE transmission plus SEE weight/area/time normalization')
    ! 根の脱出束は PE+SEE 合算の値で、PE の目標はその全量を持つ。追跡した SEE 脱出は各バッチの closure が差し引く。
    root_escape_charge = qe*receipts%matching_plane_photoelectron_escape_flux_m2_s*area*batch_time
    if (night) then
      call assert_equal_i32(closure%escape_partner_species(3), 0_i32, 'night has no PE to absorb the SEE escape')
      call assert_true(closure%has_absorbed_target(4) .and. closure%has_emission_target(4) .and. &
                       closure%has_escape_target(4) .and. closure%emission_target_from_tracked_source(4), &
                       'SEE-only closure must pin SEE return, emission and escape to the root')
      call assert_close_dp(-closure%escaped_particle_current_a(4)*batch_time, root_escape_charge, &
                           1.0e-12_dp*root_escape_charge, 'SEE escape target must carry the root escape')
      call assert_close_dp(closure%emission_current_a(4) + closure%absorbed_current_a(4) + &
                           closure%escaped_particle_current_a(4), 0.0_dp, 1.0e-12_dp*closure%emission_current_a(4), &
                           'SEE emission minus return must equal the root escape')
    else
      call assert_equal_i32(closure%escape_partner_species(3), 4_i32, 'PE must subtract the tracked SEE escape')
      call assert_close_dp(-closure%escaped_particle_current_a(3)*batch_time, root_escape_charge, &
                           1.0e-12_dp*root_escape_charge, 'PE escape target must carry the combined root escape')
    end if
    ! 窓の途中のバッチ: PE の脱出目標が差し引いた SEE 脱出とその超過が、そのバッチの share と残差になる。
    moments = 0.0_dp
    raw_charge = 0.0_dp
    closed_charge = 0.0_dp
    unclosed_charge = 0.0_dp
    raw_charge(12) = -0.5_dp*root_escape_charge
    closed_charge(3) = 0.25_dp*root_escape_charge
    unclosed_charge(3) = 0.1_dp*root_escape_charge
    closed_charge(4) = 0.8_dp*root_escape_charge
    unclosed_charge(4) = 0.2_dp*root_escape_charge
    call refresh%commit_batch(app, plate, mpi, moments, raw_charge, closed_charge, unclosed_charge, 3_i32, &
                              receipts, closure, refreshed)
    call assert_true(.not. refreshed, 'a partial window must not refresh the outer root')
    if (night) then
      call assert_close_dp(receipts%matching_plane_secondary_share, 1.0_dp, 0.0_dp, 'SEE-only outflow owns the escape')
      call assert_close_dp(receipts%secondary_budget_residual_a_m2, 0.2_dp*root_escape_charge/(area*batch_time), &
                           1.0e-12_dp*root_escape_charge/(area*batch_time), &
                           'without PE the residual is the root escape the closure could not place')
    else
      call assert_close_dp(receipts%matching_plane_secondary_share, 0.25_dp, 1.0e-12_dp, &
                           'share must be the SEE escape the PE targets absorbed')
      call assert_close_dp(receipts%secondary_budget_residual_a_m2, 0.1_dp*root_escape_charge/(area*batch_time), &
                           1.0e-12_dp*root_escape_charge/(area*batch_time), 'residual must be the unabsorbed SEE escape')
    end if
    call assert_close_dp(sum(receipts%matching_plane_photoelectron_spectrum_m2_s), flux_pe + flux_see, &
                         1.0e-12_dp*(flux_pe + flux_see), 'combined spectrum flux conservation')
    call assert_true(all(abs(receipts%matching_plane_photoelectron_spectrum_m2_s - &
                             flux_pe*spectrum_pe - flux_see*spectrum_see) < 1.0e-12_dp*(flux_pe + flux_see)), &
                     'combined spectrum must retain both bin distributions')
    if (night) call assert_true(receipts%matching_plane_response(6) < 0.0_dp, 'SEE-only outflow must permit C to A')
    restored = receipts
    call resumed_refresh%initialize(app, plate, restored, restored_closure)
    call assert_close_dp(restored%matching_plane_secondary_share, receipts%matching_plane_secondary_share, 0.0_dp, &
                         'restart must retain the last share receipt')
    call assert_close_dp(restored_closure%top_field_v_m, closure%top_field_v_m, &
                         1.0e-10_dp*abs(closure%top_field_v_m), 'restart reconstructs external top field')
    call assert_close_dp(restored_closure%escaped_particle_current_a(3), closure%escaped_particle_current_a(3), &
                         1.0e-20_dp, 'restart PE target split')
    if (.not. night) then
      ! Empty outflow keeps the finite root; without tracked SEE escape the PE targets carry the whole root escape.
      moments = 0.0_dp
      raw_charge = 0.0_dp
      closed_charge = 0.0_dp
      unclosed_charge = 0.0_dp
      restored = receipts
      do batch = 4, 5
        call refresh%commit_batch(app, plate, mpi, moments, raw_charge, closed_charge, unclosed_charge, int(batch, i32), &
                                  receipts, closure, refreshed)
        call assert_true(.not. refreshed, 'an empty window must keep the previous root')
      end do
      call assert_close_dp(receipts%matching_plane_phi_v, restored%matching_plane_phi_v, 0.0_dp, 'empty-window root')
      call assert_close_dp(receipts%matching_plane_secondary_share, 0.0_dp, 0.0_dp, 'no SEE escape gives PE share one')
      call assert_close_dp(receipts%secondary_budget_residual_a_m2, 0.0_dp, 0.0_dp, 'no SEE escape leaves no residual')
      call assert_close_dp(-closure%escaped_particle_current_a(3), &
                           receipts%matching_plane_photoelectron_escape_flux_m2_s*qe*area, 1.0e-24_dp, &
                           'retained-root PE target must keep the root escape')
    end if
  end subroutine check_emission_group_refresh

  subroutine configure_fixture(fixture_mesh, fixture_cfg, state, refresh_batches)
    type(mesh_type), intent(out) :: fixture_mesh
    type(app_config), intent(out) :: fixture_cfg
    type(injection_state), intent(out) :: state
    integer(i32), intent(in) :: refresh_batches
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
    real(dp), parameter :: plate_z = 2.0e-6_dp
    real(dp), parameter :: inward_speed = 4.0529988897e5_dp

    v0(:, 1) = [0.0_dp, 0.0_dp, plate_z]
    v1(:, 1) = [cell_width, 0.0_dp, plate_z]
    v2(:, 1) = [0.0_dp, cell_width, plate_z]
    v0(:, 2) = [cell_width, 0.0_dp, plate_z]
    v1(:, 2) = [cell_width, cell_width, plate_z]
    v2(:, 2) = [0.0_dp, cell_width, plate_z]
    call init_mesh(fixture_mesh, v0, v1, v2, q0=[0.0_dp, 0.0_dp])
    fixture_mesh%elem_vacuum_sign = 1_i32
    fixture_mesh%vacuum_normals = fixture_mesh%normals

    call default_app_config(fixture_cfg)
    fixture_cfg%sim%rng_seed = 1357_i32
    fixture_cfg%sim%batch_count = 2_i32
    fixture_cfg%sim%dt = 1.0e-11_dp
    fixture_cfg%sim%batch_duration = 1.0e-3_dp
    fixture_cfg%sim%max_step = 20000_i32
    fixture_cfg%sim%q_floor = 1.0e-30_dp
    fixture_cfg%sim%field_solver = 'direct'
    fixture_cfg%sim%field_bc_mode = 'periodic2'
    fixture_cfg%sim%use_box = .true.
    fixture_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    fixture_cfg%sim%box_max = [cell_width, cell_width, box_height]
    fixture_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    fixture_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    fixture_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    fixture_cfg%periodic2%reference_mode_layers = 1_i32
    fixture_cfg%periodic2%panel_quadrature_order = 4_i32
    fixture_cfg%surface_current%model = 'zhao_stationary'
    fixture_cfg%surface_current%zhao_branch = 'auto'
    fixture_cfg%surface_current%electron_species = 'electron'
    fixture_cfg%surface_current%ion_species = 'ion'
    fixture_cfg%surface_current%photoelectron_species = 'photoelectron'
    fixture_cfg%surface_current%solar_elevation_deg = 60.0_dp
    fixture_cfg%surface_current%photoelectron_ref_density_m3 = 64.0e6_dp
    fixture_cfg%surface_current%outflow_refresh_batches = refresh_batches
    fixture_cfg%n_particle_species = 3_i32

    call configure_ambient_species(fixture_cfg, 1_i32, 'electron', -qe, 9.1093837015e-31_dp, 12.0_dp, 0.35_dp)
    call configure_ambient_species(fixture_cfg, 2_i32, 'ion', qe, 1.67262192369e-27_dp, 0.1_dp, 0.175_dp)
    fixture_cfg%particle_species(3) = species_from_defaults()
    fixture_cfg%particle_species(3)%species_key = 'photoelectron'
    fixture_cfg%particle_species(3)%source_mode = 'photo_raycast'
    fixture_cfg%particle_species(3)%rays_per_batch = 2000_i32
    fixture_cfg%particle_species(3)%emit_current_density_a_m2 = 1.0e-4_dp
    fixture_cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    fixture_cfg%particle_species(3)%surface_charge_closure = 'fixed_current'
    fixture_cfg%particle_species(3)%q_particle = -qe
    fixture_cfg%particle_species(3)%m_particle = 9.1093837015e-31_dp
    fixture_cfg%particle_species(3)%temperature_ev = 2.2_dp
    fixture_cfg%particle_species(3)%has_temperature_ev = .true.
    fixture_cfg%particle_species(3)%normal_drift_speed = 0.0_dp
    fixture_cfg%particle_species(3)%inject_face = 'z_high'
    fixture_cfg%particle_species(3)%pos_low = [0.0_dp, 0.0_dp, box_height]
    fixture_cfg%particle_species(3)%pos_high = [cell_width, cell_width, box_height]
    fixture_cfg%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    ! 太陽風のイオンビームだけがdriftする。driftのある電子ではType A/Cのシースが成立しない。
    fixture_cfg%particle_species(1)%drift_velocity = 0.0_dp
    fixture_cfg%particle_species(2)%drift_velocity(3) = -inward_speed

    allocate (state%macro_residual(3), state%boundary_macro_residual(6, 3))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp
    call prepare_periodic2_collision_mesh(fixture_mesh, fixture_cfg%sim)
  end subroutine configure_fixture

  subroutine configure_ambient_species(fixture_cfg, species_idx, species_key, charge, mass, temperature_ev, weight)
    type(app_config), intent(inout) :: fixture_cfg
    integer(i32), intent(in) :: species_idx
    character(len=*), intent(in) :: species_key
    real(dp), intent(in) :: charge, mass, temperature_ev, weight

    fixture_cfg%particle_species(species_idx) = species_from_defaults()
    fixture_cfg%particle_species(species_idx)%species_key = species_key
    fixture_cfg%particle_species(species_idx)%source_mode = 'volume_seed'
    fixture_cfg%particle_species(species_idx)%npcls_per_step = 0_i32
    fixture_cfg%particle_species(species_idx)%boundary_inflow_high(3) = particle_inflow_reservoir
    fixture_cfg%particle_species(species_idx)%surface_charge_closure = 'fixed_current'
    fixture_cfg%particle_species(species_idx)%q_particle = charge
    fixture_cfg%particle_species(species_idx)%m_particle = mass
    fixture_cfg%particle_species(species_idx)%w_particle = weight
    fixture_cfg%particle_species(species_idx)%number_density_m3 = 8.7e6_dp
    fixture_cfg%particle_species(species_idx)%temperature_ev = temperature_ev
    fixture_cfg%particle_species(species_idx)%has_temperature_ev = .true.
  end subroutine configure_ambient_species

end program test_zhao_outflow_refresh
