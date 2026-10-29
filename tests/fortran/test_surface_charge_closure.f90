!> 粒子の表面電荷 closure（neutral_return、fixed_current とその外部根・相手種・反射電子の扱い）を検証する。
program test_surface_charge_closure
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator, apply_fixed_surface_current_closure
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_mpi, only: mpi_context
  use bem_simulator_workspace, only: simulator_batch_workspace_type
  use bem_particles, only: init_particles
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_types, only: mesh_type, particles_soa, sim_stats, bc_open, bc_reflect
  use bem_charge_ledger, only: charge_ledger_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i64, assert_close_dp, assert_allclose_1d, delete_file_if_exists, &
                          run_self_probe, file_contains
  implicit none

  character(len=*), parameter :: fixed_current_nonfinite_path = 'test_surface_charge_closure_nonfinite_tmp.log'
  character(len=*), parameter :: fixed_current_empty_absorbed_path = 'test_surface_charge_closure_empty_absorbed_tmp.log'
  character(len=64) :: run_mode

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--fixed-current-nonfinite-probe') then
    call run_fixed_current_target_failure_probe(.false.)
    error stop 'fixed-current nonfinite target probe unexpectedly completed'
  else if (trim(run_mode) == '--fixed-current-empty-absorbed-probe') then
    call run_fixed_current_target_failure_probe(.true.)
    error stop 'fixed-current empty absorbed channel probe unexpectedly completed'
  end if

  call test_init(11)

  call test_begin('neutral_return_closes_mean_charge_and_preserves_redistribution')
  call test_neutral_return_closure()
  call test_end()

  call test_begin('fixed_absorbed_current_preserves_spatial_mapping')
  call test_fixed_absorbed_current_closure()
  call test_end()

  call test_begin('fixed_current_corrections_land_on_hit_and_emitter_elements')
  call test_fixed_current_element_distribution()
  call test_end()

  call test_begin('fixed_current_nonfinite_target_fails_closed')
  call test_fixed_current_nonfinite_target_fails_closed()
  call test_end()

  call test_begin('fixed_current_empty_absorbed_channel_fails_closed')
  call test_fixed_current_empty_absorbed_channel_fails_closed()
  call test_end()

  call test_begin('fixed_photo_currents_scale_emission_and_return_separately')
  call test_fixed_photo_current_closure()
  call test_end()

  call test_begin('zhao_budget_applies_escape_closure_separately_from_raw_trajectory')
  call test_zhao_fixed_current_budget()
  call test_end()

  call test_begin('fixed_current_escape_partner_subtracts_tracked_secondary_escape')
  call test_fixed_current_escape_partner()
  call test_end()

  call test_begin('fixed_current_secondary_only_return_target_is_emission_minus_escape')
  call test_fixed_current_secondary_only_closure()
  call test_end()

  call test_begin('fixed_current_reflected_electrons_credit_the_electron_and_shrink_the_emitter')
  call test_fixed_current_reflected_electrons()
  call test_end()

  call test_begin('fixed_current_reflected_only_pins_the_electron_to_the_root')
  call test_fixed_current_reflected_only()
  call test_end()

  call test_summary()

contains

  subroutine test_neutral_return_closure()
    type(mesh_type) :: neutral_mesh
    type(app_config) :: neutral_cfg
    type(sim_stats) :: neutral_stats
    type(charge_ledger_type) :: neutral_ledger
    real(dp) :: tri_v0(3, 4), tri_v1(3, 4), tri_v2(3, 4)

    ! 左98%の面は z=0.75 にあり、dt=0.6 のうちに z-high 反射後同じ面へ帰還する。
    ! 右2%の面は z=0.25 にあり、同じ dt では帰還せず max-step survivor になる。
    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 1) = [0.98_dp, 0.0_dp, 0.75_dp]
    tri_v2(:, 1) = [0.98_dp, 1.0_dp, 0.75_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 2) = [0.98_dp, 1.0_dp, 0.75_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.75_dp]
    tri_v0(:, 3) = [0.98_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 3) = [1.0_dp, 0.0_dp, 0.25_dp]
    tri_v2(:, 3) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v0(:, 4) = [0.98_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 4) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v2(:, 4) = [0.98_dp, 1.0_dp, 0.25_dp]
    call init_mesh(neutral_mesh, tri_v0, tri_v1, tri_v2)
    neutral_mesh%elem_vacuum_sign = 1_i32
    neutral_mesh%vacuum_normals = neutral_mesh%normals

    call default_app_config(neutral_cfg)
    neutral_cfg%sim%rng_seed = 998_i32
    neutral_cfg%sim%batch_count = 1_i32
    neutral_cfg%sim%batch_duration = 1.0_dp
    neutral_cfg%sim%dt = 0.6_dp
    neutral_cfg%sim%max_step = 1_i32
    neutral_cfg%sim%q_floor = 1.0e-30_dp
    neutral_cfg%sim%use_box = .true.
    neutral_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    neutral_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    neutral_cfg%sim%bc_high(3) = bc_open
    neutral_cfg%n_particle_species = 1_i32
    neutral_cfg%particle_species(1) = species_from_defaults()
    neutral_cfg%particle_species(1)%source_mode = 'photo_raycast'
    neutral_cfg%particle_species(1)%rays_per_batch = 1024_i32
    neutral_cfg%particle_species(1)%emit_current_density_a_m2 = 1.0_dp
    neutral_cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    neutral_cfg%particle_species(1)%boundary_high(3) = bc_reflect
    neutral_cfg%particle_species(1)%surface_charge_closure = 'neutral_return'
    neutral_cfg%particle_species(1)%q_particle = -1.0_dp
    neutral_cfg%particle_species(1)%m_particle = 1.0_dp
    neutral_cfg%particle_species(1)%temperature_k = 0.0_dp
    neutral_cfg%particle_species(1)%normal_drift_speed = 1.0_dp
    neutral_cfg%particle_species(1)%inject_face = 'z_high'
    neutral_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    neutral_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    neutral_cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    neutral_cfg%particle_species(1)%has_ray_direction = .true.

    call seed_particles_from_config(neutral_cfg)
    call run_absorption_insulator( &
      neutral_mesh, neutral_cfg, neutral_stats, charge_ledger=neutral_ledger &
      )

    call assert_true(neutral_ledger%absorbed_count(1) > 0_i64, 'neutral-return fixture needs resolved returns')
    call assert_true( &
      neutral_ledger%discarded_unresolved_count(1) > 0_i64, &
      'neutral-return fixture needs unresolved survivors' &
      )
    call assert_true( &
      neutral_ledger%neutral_return_unresolved_fraction(1) <= 0.05_dp, &
      'neutral-return fixture must stay inside the fixed applicability limit' &
      )
    call assert_equal_i64(neutral_stats%escaped_boundary, 0_i64, 'neutral-return photoelectrons must not escape')
    call assert_close_dp( &
      neutral_ledger%emitted_from_surface(1), &
      neutral_ledger%absorbed_on_surface(1) + neutral_ledger%discarded_unresolved(1), &
      1.0e-12_dp, 'raw neutral-return particle charge must close' &
      )
    call assert_close_dp( &
      neutral_ledger%neutral_return_correction(1), neutral_ledger%discarded_unresolved(1), &
      1.0e-12_dp, 'neutral-return correction must represent the unresolved raw charge' &
      )
    call assert_close_dp( &
      neutral_ledger%neutral_return_weight_scale(1), &
      neutral_ledger%emitted_from_surface(1)/neutral_ledger%absorbed_on_surface(1), &
      1.0e-12_dp, 'neutral-return weight scale mismatch' &
      )
    call assert_close_dp(sum(neutral_mesh%q_elem), 0.0_dp, 1.0e-12_dp, 'neutral-return mean charge must be zero')
    call assert_true(sum(abs(neutral_mesh%q_elem)) > 0.0_dp, 'neutral-return must preserve local charge redistribution')
    call assert_close_dp(neutral_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'neutral-return ledger residual mismatch')
  end subroutine test_neutral_return_closure

  subroutine test_fixed_absorbed_current_closure()
    type(mesh_type) :: fixed_mesh
    type(app_config) :: fixed_cfg
    type(sim_stats) :: fixed_stats
    type(charge_ledger_type) :: fixed_ledger
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    tri_v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(fixed_mesh, tri_v0, tri_v1, tri_v2)
    fixed_mesh%elem_vacuum_sign = 1_i32
    fixed_mesh%vacuum_normals = fixed_mesh%normals

    call default_app_config(fixed_cfg)
    fixed_cfg%sim%batch_count = 1_i32
    fixed_cfg%sim%batch_duration = 1.0_dp
    fixed_cfg%sim%dt = 1.0_dp
    fixed_cfg%sim%max_step = 1_i32
    fixed_cfg%sim%q_floor = 1.0e-30_dp
    fixed_cfg%n_particle_species = 1_i32
    fixed_cfg%particle_species(1) = species_from_defaults()
    fixed_cfg%particle_species(1)%source_mode = 'volume_seed'
    fixed_cfg%particle_species(1)%npcls_per_step = 1_i32
    fixed_cfg%particle_species(1)%q_particle = 1.0_dp
    fixed_cfg%particle_species(1)%m_particle = 1.0_dp
    fixed_cfg%particle_species(1)%w_particle = 2.0_dp
    fixed_cfg%particle_species(1)%pos_low = [0.2_dp, 0.2_dp, 0.5_dp]
    fixed_cfg%particle_species(1)%pos_high = fixed_cfg%particle_species(1)%pos_low
    fixed_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    fixed_cfg%particle_species(1)%temperature_k = 0.0_dp
    fixed_cfg%particle_species(1)%surface_charge_closure = 'fixed_current'
    fixed_cfg%particle_species(1)%target_absorbed_current_a = 5.0_dp
    fixed_cfg%particle_species(1)%has_target_absorbed_current_a = .true.

    call seed_particles_from_config(fixed_cfg)
    call run_absorption_insulator(fixed_mesh, fixed_cfg, fixed_stats, charge_ledger=fixed_ledger)
    call assert_equal_i64( &
      fixed_ledger%absorbed_count(1), 1_i64, &
      'single-sample fixed-current mapping must remain visible in raw diagnostics' &
      )
    call assert_close_dp( &
      fixed_ledger%absorbed_on_surface(1), 2.0_dp, 1.0e-12_dp, &
      'single-sample fixed-current mapping must preserve raw absorbed charge' &
      )
    call assert_close_dp(fixed_mesh%q_elem(1), 5.0_dp, 1.0e-12_dp, 'fixed absorbed target charge mismatch')
    call assert_close_dp( &
      fixed_ledger%fixed_absorbed_weight_scale(1), 2.5_dp, 1.0e-12_dp, &
      'fixed absorbed current scale mismatch' &
      )
    call assert_close_dp( &
      fixed_ledger%fixed_current_correction(1), 3.0_dp, 1.0e-12_dp, &
      'fixed absorbed current correction mismatch' &
      )
    call assert_close_dp(fixed_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'fixed absorbed ledger residual mismatch')
  end subroutine test_fixed_absorbed_current_closure

  !> 吸収と放出の補正が、それぞれの raw 電荷が載った要素へ raw 電荷に比例して配られることを要素ごとに確かめる。
  !! 正しい総量を 1 要素へまとめる、または hit/放出要素を取り違える実装を検出する。
  subroutine test_fixed_current_element_distribution()
    type(app_config) :: app
    type(surface_closure_contract_type) :: model
    type(particles_soa) :: batch
    type(simulator_batch_workspace_type) :: workspace
    type(mpi_context) :: mpi
    real(dp) :: x(3, 4), v(3, 4)

    call default_app_config(app)
    app%sim%batch_duration = 0.5_dp
    app%n_particle_species = 1_i32
    app%particle_species(1) = species_from_defaults()
    app%particle_species(1)%q_particle = -1.0_dp
    app%particle_species(1)%surface_charge_closure = 'fixed_current'
    ! Targets -2.625 C absorbed and 2 C emission reaction against raw -3.5 C and 4 C: scales 0.75 and 0.5.
    app%particle_species(1)%target_absorbed_current_a = -5.25_dp
    app%particle_species(1)%has_target_absorbed_current_a = .true.
    app%particle_species(1)%target_emission_current_a = 4.0_dp
    app%particle_species(1)%has_target_emission_current_a = .true.
    allocate (model%has_absorbed_target(1), model%has_emission_target(1), model%has_escape_target(1))
    model%has_absorbed_target = .false.
    model%has_emission_target = .false.
    model%has_escape_target = .false.

    ! Particles 1-3 are absorbed on elements 1, 3 and 3 with weights 1, 2 and 0.5; particle 4 escapes.
    x = 0.0_dp
    v = 0.0_dp
    call init_particles(batch, x, v, [-1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
                        [1.0_dp, 2.0_dp, 0.5_dp, 1.0_dp], species_id=[1_i32, 1_i32, 1_i32, 1_i32], &
                        source_element=[-1_i32, -1_i32, -1_i32, -1_i32])
    call workspace%init(3_i32, 1_i32, 1_i32)
    call workspace%prepare_particle_flags(4_i32)
    workspace%absorbed_flag(1:3) = .true.
    workspace%absorbed_element(1:3) = [1_i32, 3_i32, 3_i32]
    workspace%escaped_boundary_flag(4) = .true.
    ! Raw emission reaction charge sits on elements 2 and 3 only.
    workspace%photo_emission_dq(:, 1) = [0.0_dp, 3.0_dp, 1.0_dp]
    call apply_fixed_surface_current_closure(app, model, batch, 4_i32, workspace, mpi)

    call assert_close_dp(workspace%fixed_absorbed_weight_scale(1), 0.75_dp, 1.0e-15_dp, 'absorbed scale')
    call assert_close_dp(workspace%fixed_emission_weight_scale(1), 0.5_dp, 1.0e-15_dp, 'emission scale')
    ! Absorbed correction (scale - 1) q w on each hit element: element 1 gets +0.25, element 3 gets +0.625.
    call assert_allclose_1d(workspace%dq_thread(:, 1), [0.25_dp, 0.0_dp, 0.625_dp], 1.0e-15_dp, &
                            'absorbed correction must land on the hit elements in proportion to their raw charge')
    call assert_allclose_1d(workspace%photo_emission_dq(:, 1), [0.0_dp, 1.5_dp, 0.5_dp], 1.0e-15_dp, &
                            'emission reaction must stay on the emitting elements, scaled')
  end subroutine test_fixed_current_element_distribution

  subroutine test_fixed_current_nonfinite_target_fails_closed()
    call assert_fixed_current_failure_probe( &
      '--fixed-current-nonfinite-probe', fixed_current_nonfinite_path, &
      'fixed_current absorbed target charge overflowed', 'fixed-current nonfinite target' &
      )
  end subroutine test_fixed_current_nonfinite_target_fails_closed

  subroutine test_fixed_current_empty_absorbed_channel_fails_closed()
    call assert_fixed_current_failure_probe( &
      '--fixed-current-empty-absorbed-probe', fixed_current_empty_absorbed_path, &
      'fixed_current cannot map a nonzero absorbed target onto an empty raw channel', &
      'fixed-current empty absorbed channel' &
      )
  end subroutine test_fixed_current_empty_absorbed_channel_fails_closed

  subroutine assert_fixed_current_failure_probe(run_mode, output_path, expected_diagnostic, context)
    character(len=*), intent(in) :: run_mode, output_path, expected_diagnostic, context
    integer :: exit_status
    logical :: saw_expected_error

    call run_self_probe(run_mode, output_path, exit_status)
    call assert_true(exit_status /= 0, trim(context)//' probe should terminate with nonzero status')
    saw_expected_error = file_contains(output_path, expected_diagnostic)
    call delete_file_if_exists(output_path)
    call assert_true(saw_expected_error, trim(context)//' probe should report the rejected contract')
  end subroutine assert_fixed_current_failure_probe

  subroutine run_fixed_current_target_failure_probe(empty_absorbed_channel)
    logical, intent(in) :: empty_absorbed_channel
    type(mesh_type) :: fixed_mesh
    type(app_config) :: fixed_cfg
    type(sim_stats) :: fixed_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    tri_v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(fixed_mesh, tri_v0, tri_v1, tri_v2)
    fixed_mesh%elem_vacuum_sign = 1_i32
    fixed_mesh%vacuum_normals = fixed_mesh%normals

    call default_app_config(fixed_cfg)
    fixed_cfg%sim%batch_count = 1_i32
    fixed_cfg%sim%batch_duration = 1.5_dp
    fixed_cfg%sim%dt = 1.0_dp
    fixed_cfg%sim%max_step = 1_i32
    fixed_cfg%n_particle_species = 1_i32
    fixed_cfg%particle_species(1) = species_from_defaults()
    fixed_cfg%particle_species(1)%source_mode = 'volume_seed'
    fixed_cfg%particle_species(1)%npcls_per_step = 1_i32
    fixed_cfg%particle_species(1)%q_particle = 1.0_dp
    fixed_cfg%particle_species(1)%m_particle = 1.0_dp
    fixed_cfg%particle_species(1)%w_particle = 1.0_dp
    fixed_cfg%particle_species(1)%pos_low = [0.2_dp, 0.2_dp, 0.5_dp]
    fixed_cfg%particle_species(1)%pos_high = fixed_cfg%particle_species(1)%pos_low
    if (empty_absorbed_channel) then
      fixed_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, 1.0_dp]
    else
      fixed_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    end if
    fixed_cfg%particle_species(1)%temperature_k = 0.0_dp
    fixed_cfg%particle_species(1)%surface_charge_closure = 'fixed_current'
    if (empty_absorbed_channel) then
      fixed_cfg%particle_species(1)%target_absorbed_current_a = 1.0_dp
    else
      ! Division rounds this value up to the precheck boundary, while the product still overflows.
      fixed_cfg%particle_species(1)%target_absorbed_current_a = huge(1.0_dp)/fixed_cfg%sim%batch_duration
    end if
    fixed_cfg%particle_species(1)%has_target_absorbed_current_a = .true.

    call seed_particles_from_config(fixed_cfg)
    call run_absorption_insulator(fixed_mesh, fixed_cfg, fixed_stats)
  end subroutine run_fixed_current_target_failure_probe

  subroutine test_fixed_photo_current_closure()
    type(mesh_type) :: fixed_mesh
    type(app_config) :: fixed_cfg
    type(sim_stats) :: fixed_stats
    type(charge_ledger_type) :: fixed_ledger
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.75_dp]
    tri_v2(:, 1) = [1.0_dp, 1.0_dp, 0.75_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 2) = [1.0_dp, 1.0_dp, 0.75_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.75_dp]
    call init_mesh(fixed_mesh, tri_v0, tri_v1, tri_v2)
    fixed_mesh%elem_vacuum_sign = 1_i32
    fixed_mesh%vacuum_normals = fixed_mesh%normals

    call default_app_config(fixed_cfg)
    fixed_cfg%sim%rng_seed = 999_i32
    fixed_cfg%sim%batch_count = 1_i32
    fixed_cfg%sim%batch_duration = 1.0_dp
    fixed_cfg%sim%dt = 0.6_dp
    fixed_cfg%sim%max_step = 1_i32
    fixed_cfg%sim%q_floor = 1.0e-30_dp
    fixed_cfg%sim%use_box = .true.
    fixed_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    fixed_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    fixed_cfg%sim%bc_high(3) = bc_open
    fixed_cfg%n_particle_species = 1_i32
    fixed_cfg%particle_species(1) = species_from_defaults()
    fixed_cfg%particle_species(1)%source_mode = 'photo_raycast'
    fixed_cfg%particle_species(1)%rays_per_batch = 256_i32
    fixed_cfg%particle_species(1)%emit_current_density_a_m2 = 1.0_dp
    fixed_cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    fixed_cfg%particle_species(1)%boundary_high(3) = bc_reflect
    fixed_cfg%particle_species(1)%surface_charge_closure = 'fixed_current'
    fixed_cfg%particle_species(1)%target_absorbed_current_a = -0.25_dp
    fixed_cfg%particle_species(1)%has_target_absorbed_current_a = .true.
    fixed_cfg%particle_species(1)%target_emission_current_a = 2.0_dp
    fixed_cfg%particle_species(1)%has_target_emission_current_a = .true.
    fixed_cfg%particle_species(1)%q_particle = -1.0_dp
    fixed_cfg%particle_species(1)%m_particle = 1.0_dp
    fixed_cfg%particle_species(1)%temperature_k = 0.0_dp
    fixed_cfg%particle_species(1)%normal_drift_speed = 1.0_dp
    fixed_cfg%particle_species(1)%inject_face = 'z_high'
    fixed_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    fixed_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    fixed_cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    fixed_cfg%particle_species(1)%has_ray_direction = .true.

    call seed_particles_from_config(fixed_cfg)
    call run_absorption_insulator(fixed_mesh, fixed_cfg, fixed_stats, charge_ledger=fixed_ledger)
    call assert_equal_i64(fixed_stats%absorbed, 256_i64, 'fixed PE fixture must resolve every return')
    call assert_close_dp(sum(fixed_mesh%q_elem), 1.75_dp, 1.0e-12_dp, 'fixed PE net surface charge mismatch')
    call assert_close_dp( &
      fixed_ledger%fixed_absorbed_target_charge(1), -0.25_dp, 1.0e-12_dp, &
      'fixed PE return target mismatch' &
      )
    call assert_close_dp( &
      fixed_ledger%absorbed_on_surface(1)*fixed_ledger%fixed_absorbed_weight_scale(1), &
      fixed_ledger%fixed_absorbed_target_charge(1), 1.0e-12_dp, &
      'fixed PE return raw charge and scale must reproduce the applied target' &
      )
    call assert_close_dp( &
      fixed_ledger%fixed_emission_target_charge(1), 2.0_dp, 1.0e-12_dp, &
      'fixed PE emission target mismatch' &
      )
    call assert_close_dp( &
      -fixed_ledger%emitted_from_surface(1)*fixed_ledger%fixed_emission_weight_scale(1), &
      fixed_ledger%fixed_emission_target_charge(1), 1.0e-12_dp, &
      'fixed PE emission raw charge and scale must reproduce the applied target' &
      )
    call assert_true( &
      abs(fixed_ledger%fixed_absorbed_weight_scale(1) - fixed_ledger%fixed_emission_weight_scale(1)) > 1.0e-6_dp, &
      'fixed PE emission and return channels must use independent scales' &
      )
    call assert_close_dp(fixed_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'fixed PE ledger residual mismatch')
  end subroutine test_fixed_photo_current_closure

  subroutine test_zhao_fixed_current_budget()
    real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
    real(dp), parameter :: proton_mass = 1.67262192369e-27_dp
    type(mesh_type) :: zhao_mesh
    type(app_config) :: zhao_cfg
    type(sim_stats) :: zhao_stats
    type(charge_ledger_type) :: zhao_ledger
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2), inward_speed
    integer(i32) :: photo_idx

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 7.5e5_dp]
    tri_v1(:, 1) = [1.0e7_dp, 0.0_dp, 7.5e5_dp]
    tri_v2(:, 1) = [1.0e7_dp, 1.0e7_dp, 7.5e5_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 7.5e5_dp]
    tri_v1(:, 2) = [1.0e7_dp, 1.0e7_dp, 7.5e5_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0e7_dp, 7.5e5_dp]
    call init_mesh(zhao_mesh, tri_v0, tri_v1, tri_v2)
    zhao_mesh%elem_vacuum_sign = 1_i32
    zhao_mesh%vacuum_normals = zhao_mesh%normals

    call default_app_config(zhao_cfg)
    zhao_cfg%sim%rng_seed = 1234_i32
    zhao_cfg%sim%batch_count = 1_i32
    zhao_cfg%sim%batch_duration = 1.0_dp
    zhao_cfg%sim%dt = 1.0_dp
    zhao_cfg%sim%max_step = 2_i32
    zhao_cfg%sim%q_floor = 1.0e-30_dp
    zhao_cfg%sim%use_box = .true.
    zhao_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    zhao_cfg%sim%box_max = [1.0e7_dp, 1.0e7_dp, 1.0e6_dp]
    zhao_cfg%sim%bc_high(3) = bc_open
    zhao_cfg%n_particle_species = 3_i32
    inward_speed = 4.68e5_dp*sin(60.0_dp*acos(-1.0_dp)/180.0_dp)

    zhao_cfg%particle_species(1) = species_from_defaults()
    zhao_cfg%particle_species(1)%species_key = 'electron'
    zhao_cfg%particle_species(1)%source_mode = 'volume_seed'
    zhao_cfg%particle_species(1)%npcls_per_step = 512_i32
    zhao_cfg%particle_species(1)%q_particle = -1.0_dp
    zhao_cfg%particle_species(1)%m_particle = electron_mass
    zhao_cfg%particle_species(1)%w_particle = 1.0_dp
    zhao_cfg%particle_species(1)%temperature_ev = 12.0_dp
    zhao_cfg%particle_species(1)%has_temperature_ev = .true.
    zhao_cfg%particle_species(1)%drift_velocity = 0.0_dp
    zhao_cfg%particle_species(1)%pos_low = [5.0e6_dp, 5.0e6_dp, 9.0e5_dp]
    zhao_cfg%particle_species(1)%pos_high = zhao_cfg%particle_species(1)%pos_low
    zhao_cfg%particle_species(1)%surface_charge_closure = 'fixed_current'

    zhao_cfg%particle_species(2) = species_from_defaults()
    zhao_cfg%particle_species(2)%species_key = 'ion'
    zhao_cfg%particle_species(2)%source_mode = 'volume_seed'
    zhao_cfg%particle_species(2)%npcls_per_step = 128_i32
    zhao_cfg%particle_species(2)%q_particle = 1.0_dp
    zhao_cfg%particle_species(2)%m_particle = proton_mass
    zhao_cfg%particle_species(2)%w_particle = 1.0_dp
    zhao_cfg%particle_species(2)%number_density_m3 = 8.7e6_dp
    zhao_cfg%particle_species(2)%temperature_ev = 0.1_dp
    zhao_cfg%particle_species(2)%has_temperature_ev = .true.
    zhao_cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -inward_speed]
    zhao_cfg%particle_species(2)%pos_low = [5.0e6_dp, 5.0e6_dp, 9.0e5_dp]
    zhao_cfg%particle_species(2)%pos_high = zhao_cfg%particle_species(2)%pos_low
    zhao_cfg%particle_species(2)%surface_charge_closure = 'fixed_current'

    photo_idx = 3_i32
    zhao_cfg%particle_species(photo_idx) = species_from_defaults()
    zhao_cfg%particle_species(photo_idx)%species_key = 'photoelectron'
    zhao_cfg%particle_species(photo_idx)%source_mode = 'photo_raycast'
    zhao_cfg%particle_species(photo_idx)%rays_per_batch = 256_i32
    zhao_cfg%particle_species(photo_idx)%emit_current_density_a_m2 = 1.0_dp
    zhao_cfg%particle_species(photo_idx)%deposit_opposite_charge_on_emit = .true.
    zhao_cfg%particle_species(photo_idx)%boundary_high(3) = bc_reflect
    zhao_cfg%particle_species(photo_idx)%surface_charge_closure = 'fixed_current'
    zhao_cfg%particle_species(photo_idx)%q_particle = -1.0_dp
    zhao_cfg%particle_species(photo_idx)%m_particle = electron_mass
    zhao_cfg%particle_species(photo_idx)%temperature_ev = 2.2_dp
    zhao_cfg%particle_species(photo_idx)%has_temperature_ev = .true.
    zhao_cfg%particle_species(photo_idx)%normal_drift_speed = 1.0e6_dp
    zhao_cfg%particle_species(photo_idx)%inject_face = 'z_high'
    zhao_cfg%particle_species(photo_idx)%pos_low = [0.0_dp, 0.0_dp, 1.0e6_dp]
    zhao_cfg%particle_species(photo_idx)%pos_high = [1.0e7_dp, 1.0e7_dp, 1.0e6_dp]
    zhao_cfg%particle_species(photo_idx)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    zhao_cfg%particle_species(photo_idx)%has_ray_direction = .true.

    zhao_cfg%surface_current%model = 'zhao_stationary'
    zhao_cfg%surface_current%zhao_branch = 'auto'
    zhao_cfg%surface_current%electron_species = 'electron'
    zhao_cfg%surface_current%ion_species = 'ion'
    zhao_cfg%surface_current%photoelectron_species = 'photoelectron'
    zhao_cfg%surface_current%solar_elevation_deg = 60.0_dp
    zhao_cfg%surface_current%photoelectron_ref_density_m3 = 64.0e6_dp
    zhao_cfg%surface_current%reference_area_m2 = 1.0_dp
    zhao_cfg%surface_current%has_reference_area_m2 = .true.

    call seed_particles_from_config(zhao_cfg)
    call run_absorption_insulator(zhao_mesh, zhao_cfg, zhao_stats, charge_ledger=zhao_ledger)
    call assert_true(zhao_ledger%absorbed_count(1) > 0_i64, 'Zhao fixture needs raw electron absorption')
    call assert_true(zhao_ledger%absorbed_count(2) > 0_i64, 'Zhao fixture needs raw ion absorption')
    call assert_true(zhao_ledger%absorbed_count(photo_idx) > 0_i64, 'Zhao fixture needs raw PE return')
    call assert_close_dp( &
      zhao_ledger%escaped_to_infinity(photo_idx), 0.0_dp, 0.0_dp, &
      'reflected raw PE trajectory must not escape' &
      )
    call assert_true( &
      zhao_ledger%fixed_escape_target_charge(photo_idx) < 0.0_dp, &
      'Zhao PE escape target must carry negative particle charge outward' &
      )
    call assert_close_dp( &
      zhao_ledger%fixed_escape_correction(photo_idx), &
      zhao_ledger%fixed_escape_target_charge(photo_idx), 1.0e-18_dp, &
      'Zhao PE escape correction must separate target from zero raw escape' &
      )
    call assert_close_dp( &
      zhao_ledger%fixed_absorbed_target_charge(photo_idx) + &
      zhao_ledger%fixed_emission_target_charge(photo_idx) + &
      zhao_ledger%fixed_escape_target_charge(photo_idx), &
      0.0_dp, 1.0e-12_dp, 'Zhao PE applied channel continuity mismatch' &
      )
    call assert_close_dp( &
      sum(zhao_ledger%fixed_absorbed_target_charge) + sum(zhao_ledger%fixed_emission_target_charge), &
      0.0_dp, 1.0e-12_dp, 'Zhao applied surface-current budget must close' &
      )
  end subroutine test_zhao_fixed_current_budget

  !> PE の脱出・帰還目標は、根の脱出束からそのバッチで追跡した SEE 脱出を差し引いた値になる。
  subroutine test_fixed_current_escape_partner()
    type(app_config) :: app
    type(surface_closure_contract_type) :: model
    type(particles_soa) :: batch
    type(simulator_batch_workspace_type) :: workspace
    type(mpi_context) :: mpi
    real(dp) :: x(3, 3), v(3, 3)
    real(dp), parameter :: escape_current = 2.0_dp, emission_current = 5.0_dp

    call default_app_config(app)
    app%sim%batch_duration = 1.0_dp
    app%n_particle_species = 2_i32
    app%particle_species(1) = species_from_defaults()
    app%particle_species(1)%species_key = 'photoelectron'
    app%particle_species(1)%q_particle = -1.0_dp
    app%particle_species(1)%surface_charge_closure = 'fixed_current'
    app%particle_species(2) = species_from_defaults()
    app%particle_species(2)%species_key = 'secondary'
    app%particle_species(2)%q_particle = -1.0_dp

    model%active = .true.
    allocate (model%has_absorbed_target(2), model%has_emission_target(2), model%has_escape_target(2), &
              model%absorbed_current_a(2), model%emission_current_a(2), model%escaped_particle_current_a(2), &
              model%escape_partner_species(2))
    model%has_absorbed_target = [.true., .false.]
    model%has_emission_target = [.true., .false.]
    model%has_escape_target = [.true., .false.]
    model%absorbed_current_a = [escape_current - emission_current, 0.0_dp]
    model%emission_current_a = [emission_current, 0.0_dp]
    model%escaped_particle_current_a = [-escape_current, 0.0_dp]
    model%escape_partner_species = [2_i32, 0_i32]

    ! 粒子1: 帰還した PE、粒子2: 脱出した PE、粒子3: 脱出した SEE（重み 0.5）。
    x = 0.0_dp
    v = 0.0_dp
    call init_particles(batch, x, v, [-1.0_dp, -1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], &
                        [1.0_dp, 1.0_dp, 0.5_dp], species_id=[1_i32, 1_i32, 2_i32], &
                        source_element=[1_i32, 1_i32, 1_i32])
    call workspace%init(1_i32, 2_i32, 1_i32)
    call prepare_partner_fixture(workspace)
    call apply_fixed_surface_current_closure(app, model, batch, 3_i32, workspace, mpi)
    ! SEE が 0.5 脱出したので、PE の脱出目標は 2 − 0.5 = 1.5、帰還目標は 1.5 − 5 = −3.5 になる。
    call assert_close_dp(workspace%fixed_escape_closed_charge(1), 0.5_dp, 0.0_dp, 'partner escape charge')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(1), 0.0_dp, 0.0_dp, 'no excess below the root escape')
    call assert_close_dp(workspace%fixed_escape_target_charge(1), -1.5_dp, 1.0e-12_dp, 'PE escape target minus SEE escape')
    call assert_close_dp(workspace%fixed_escape_correction(1), -0.5_dp, 1.0e-12_dp, 'PE escape correction vs raw escape')
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), -3.5_dp, 1.0e-12_dp, &
                         'PE return target must absorb the SEE escape')
    call assert_close_dp(workspace%fixed_absorbed_weight_scale(1), 3.5_dp, 1.0e-12_dp, 'PE return scale')
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1) + workspace%fixed_emission_target_charge(1) + &
                         workspace%fixed_escape_target_charge(1), 0.0_dp, 1.0e-12_dp, 'PE channel continuity')
    call assert_close_dp(-workspace%fixed_escape_target_charge(1) - workspace%fixed_current_charge_values(2, 3), &
                         escape_current, 1.0e-12_dp, 'PE and SEE escape must add up to the root escape')
    call assert_true(all(workspace%fixed_escape_closed_charge(2:) == 0.0_dp), 'SEE itself has no partner')

    ! SEE の脱出が根の脱出束を超えると、PE の脱出目標は 0、帰還目標は放出全量、超過は残差に残る。
    batch%w(3) = 3.0_dp
    call prepare_partner_fixture(workspace)
    call apply_fixed_surface_current_closure(app, model, batch, 3_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_escape_closed_charge(1), escape_current, 0.0_dp, 'partner charge clamps at root')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(1), 1.0_dp, 1.0e-12_dp, 'excess SEE escape')
    call assert_close_dp(workspace%fixed_escape_target_charge(1), 0.0_dp, 0.0_dp, 'PE escape target clamps at zero')
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), -emission_current, 1.0e-12_dp, &
                         'all PE return when SEE carries the whole root escape')
  end subroutine test_fixed_current_escape_partner

  subroutine prepare_partner_fixture(workspace)
    type(simulator_batch_workspace_type), intent(inout) :: workspace

    call workspace%reset_before_injection()
    call workspace%prepare_particle_flags(3_i32)
    workspace%absorbed_flag(1) = .true.
    workspace%absorbed_element(1) = 1_i32
    workspace%escaped_boundary_flag(2:3) = .true.
    ! raw の PE 反作用電荷（放出 2 個分）。
    workspace%photo_emission_dq(1, 1) = 2.0_dp
  end subroutine prepare_partner_fixture

  !> PE がない SEE だけの構成では、SEE の再吸収目標が「raw 放出 − 根の脱出」になり、
  !! 根の脱出の超過と載せられない目標は未閉の電荷として残る。
  subroutine test_fixed_current_secondary_only_closure()
    type(app_config) :: app
    type(surface_closure_contract_type) :: model
    type(particles_soa) :: batch
    type(simulator_batch_workspace_type) :: workspace
    type(mpi_context) :: mpi
    real(dp) :: x(3, 2), v(3, 2)

    call default_app_config(app)
    app%sim%batch_duration = 1.0_dp
    app%n_particle_species = 1_i32
    app%particle_species(1) = species_from_defaults()
    app%particle_species(1)%species_key = 'secondary'
    app%particle_species(1)%q_particle = -1.0_dp

    ! 根: 帰還 0.1、脱出 0.3、放出目標（根に与えた外側源）0.4。放出 − 帰還 = 脱出。
    model%active = .true.
    allocate (model%has_absorbed_target(1), model%has_emission_target(1), model%has_escape_target(1), &
              model%absorbed_current_a(1), model%emission_current_a(1), model%escaped_particle_current_a(1), &
              model%escape_partner_species(1), model%emission_target_from_tracked_source(1))
    model%has_absorbed_target = .true.
    model%has_emission_target = .true.
    model%has_escape_target = .true.
    model%absorbed_current_a = -0.1_dp
    model%emission_current_a = 0.4_dp
    model%escaped_particle_current_a = -0.3_dp
    model%escape_partner_species = 0_i32
    model%emission_target_from_tracked_source = .true.

    ! 粒子1: 再吸収した SEE（重み 0.5）、粒子2: 脱出した SEE（重み 0.5）。raw の放出反作用電荷は 1.0。
    x = 0.0_dp
    v = 0.0_dp
    call init_particles(batch, x, v, [-1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp], [0.5_dp, 0.5_dp], &
                        species_id=[1_i32, 1_i32], source_element=[1_i32, 1_i32])
    call workspace%init(1_i32, 1_i32, 1_i32)
    call prepare_secondary_only_fixture(workspace, .true., .true., 1.0_dp)
    call apply_fixed_surface_current_closure(app, model, batch, 2_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), -0.1_dp, 1.0e-12_dp, 'SEE return target = root return')
    call assert_close_dp(workspace%fixed_absorbed_weight_scale(1), 0.2_dp, 1.0e-12_dp, 'SEE return scale')
    call assert_close_dp(workspace%fixed_emission_target_charge(1), 0.4_dp, 1.0e-12_dp, &
                         'SEE emission reaction = root outward source')
    call assert_close_dp(workspace%fixed_emission_weight_scale(1), 0.4_dp, 1.0e-12_dp, 'SEE emission scale')
    call assert_close_dp(workspace%photo_emission_dq(1, 1), 0.4_dp, 1.0e-12_dp, 'applied emission reaction charge')
    call assert_close_dp(workspace%fixed_escape_target_charge(1), -0.3_dp, 1.0e-12_dp, 'SEE escape target = root escape')
    call assert_close_dp(workspace%fixed_escape_correction(1), 0.2_dp, 1.0e-12_dp, 'SEE escape correction vs raw escape')
    call assert_close_dp(workspace%fixed_current_correction(1), -0.2_dp, 1.0e-12_dp, 'return and emission corrections')
    call assert_close_dp(workspace%fixed_escape_closed_charge(1), 0.3_dp, 1.0e-12_dp, 'root escape placed in the budget')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(1), 0.0_dp, 0.0_dp, 'nothing left unclosed')
    call assert_close_dp(workspace%fixed_emission_target_charge(1) + workspace%fixed_absorbed_target_charge(1) + &
                         workspace%fixed_escape_target_charge(1), 0.0_dp, 1.0e-12_dp, &
                         'SEE emission, return and escape targets must close')

    ! 根の脱出がそのバッチの放出を超えても、放出の反作用電荷を根の外側源まで広げて閉じる（繰り越しなし）。
    model%absorbed_current_a = -0.1_dp
    model%emission_current_a = 1.6_dp
    model%escaped_particle_current_a = -1.5_dp
    call prepare_secondary_only_fixture(workspace, .true., .true., 1.0_dp)
    call apply_fixed_surface_current_closure(app, model, batch, 2_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_emission_weight_scale(1), 1.6_dp, 1.0e-12_dp, 'emission scales up to the root source')
    call assert_close_dp(workspace%photo_emission_dq(1, 1), 1.6_dp, 1.0e-12_dp, 'applied emission reaction charge')
    call assert_close_dp(workspace%fixed_escape_closed_charge(1), 1.5_dp, 1.0e-12_dp, 'whole root escape placed')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(1), 0.0_dp, 0.0_dp, 'nothing carried over')

    ! 生の帰還がないバッチでは帰還目標を放出へ畳む。放出 − 帰還 = 脱出は変わらない。
    model%absorbed_current_a = -0.1_dp
    model%emission_current_a = 0.4_dp
    model%escaped_particle_current_a = -0.3_dp
    call prepare_secondary_only_fixture(workspace, .false., .true., 1.0_dp)
    call apply_fixed_surface_current_closure(app, model, batch, 2_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), 0.0_dp, 0.0_dp, 'return target folds into the emission')
    call assert_close_dp(workspace%fixed_absorbed_weight_scale(1), 1.0_dp, 0.0_dp, 'empty raw return keeps the return scale one')
    call assert_close_dp(workspace%fixed_emission_target_charge(1), 0.3_dp, 1.0e-12_dp, &
                         'emission reaction = root source minus the folded return')
    call assert_close_dp(workspace%photo_emission_dq(1, 1), 0.3_dp, 1.0e-12_dp, 'applied emission reaction charge')
    call assert_close_dp(workspace%fixed_current_correction(1), -0.7_dp, 1.0e-12_dp, 'emission correction closes the budget')
    call assert_close_dp(workspace%fixed_escape_closed_charge(1), 0.3_dp, 1.0e-12_dp, 'root escape closed through the emission')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(1), 0.0_dp, 0.0_dp, 'nothing unclosed through the emission')

    ! 源 0 の初期根（帰還 0、脱出 0）では反作用電荷を置かず、セル総電荷を動かさない。
    model%absorbed_current_a = 0.0_dp
    model%emission_current_a = 0.0_dp
    model%escaped_particle_current_a = 0.0_dp
    call prepare_secondary_only_fixture(workspace, .false., .true., 1.0_dp)
    call apply_fixed_surface_current_closure(app, model, batch, 2_i32, workspace, mpi)
    call assert_close_dp(workspace%photo_emission_dq(1, 1), 0.0_dp, 0.0_dp, 'zero-source root leaves no SEE reaction charge')
    call assert_close_dp(workspace%fixed_emission_weight_scale(1), 0.0_dp, 0.0_dp, 'zero-source emission scale')
    call assert_close_dp(workspace%fixed_current_correction(1), -1.0_dp, 1.0e-12_dp, 'the whole emission is corrected away')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(1), 0.0_dp, 0.0_dp, 'zero-source batch stays closed')
    call assert_close_dp(workspace%fixed_escape_target_charge(1), 0.0_dp, 0.0_dp, 'zero-source root has no escape target')

    ! SEE が放出されなかったバッチには目標を載せられない。根の脱出分を未閉として記録し、止まらない。
    model%absorbed_current_a = -0.1_dp
    model%emission_current_a = 0.4_dp
    model%escaped_particle_current_a = -0.3_dp
    call prepare_secondary_only_fixture(workspace, .false., .false., 0.0_dp)
    call apply_fixed_surface_current_closure(app, model, batch, 2_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_emission_target_charge(1), 0.0_dp, 0.0_dp, 'no emission target without emission')
    call assert_close_dp(workspace%fixed_emission_weight_scale(1), 1.0_dp, 0.0_dp, 'empty raw emission keeps scale one')
    call assert_close_dp(workspace%fixed_current_correction(1), 0.0_dp, 0.0_dp, 'no correction without emission')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(1), 0.3_dp, 1.0e-12_dp, &
                         'root escape left unclosed when no SEE was emitted')
    call assert_close_dp(workspace%fixed_escape_closed_charge(1), 0.0_dp, 0.0_dp, 'nothing placed without emission')
    call assert_close_dp(workspace%fixed_escape_target_charge(1), -0.3_dp, 1.0e-12_dp, 'escape target is still recorded')
  end subroutine test_fixed_current_secondary_only_closure

  !> 反射して H を抜けた上流電子を外部源に含めたときの closure: 電子の吸収目標には追跡した脱出を足し戻し、
  !! 放出電子の目標からは追跡した反射電子の脱出（追跡した源の種では帰還と放出も）を差し引く。
  subroutine test_fixed_current_reflected_electrons()
    type(app_config) :: app
    type(surface_closure_contract_type) :: model
    type(particles_soa) :: batch
    type(simulator_batch_workspace_type) :: workspace
    type(mpi_context) :: mpi
    real(dp) :: x(3, 4), v(3, 4)

    call default_app_config(app)
    app%sim%batch_duration = 1.0_dp
    app%n_particle_species = 2_i32
    app%particle_species(1) = species_from_defaults()
    app%particle_species(1)%species_key = 'electron'
    app%particle_species(1)%q_particle = -1.0_dp
    app%particle_species(2) = species_from_defaults()
    app%particle_species(2)%species_key = 'photoelectron'
    app%particle_species(2)%q_particle = -1.0_dp

    model%active = .true.
    allocate (model%has_absorbed_target(2), model%has_emission_target(2), model%has_escape_target(2), &
              model%absorbed_current_a(2), model%emission_current_a(2), model%escaped_particle_current_a(2), &
              model%escape_partner_species(2), model%emission_target_from_tracked_source(2), &
              model%reflected_electron_species(2), model%absorbed_target_credits_tracked_escape(2))
    model%has_absorbed_target = .true.
    model%has_emission_target = [.false., .true.]
    model%has_escape_target = [.false., .true.]
    model%escape_partner_species = 0_i32
    model%emission_target_from_tracked_source = .false.
    model%reflected_electron_species = [0_i32, 1_i32]
    model%absorbed_target_credits_tracked_escape = [.true., .false.]
    ! 根: 電子流入 1.0、PE 放出 1.0、脱出 0.5（反射電子を含む合算）、PE 帰還 = 放出 − 脱出 = 0.5。
    model%absorbed_current_a = [-1.0_dp, -0.5_dp]
    model%emission_current_a = [0.0_dp, 1.0_dp]
    model%escaped_particle_current_a = [0.0_dp, -0.5_dp]

    ! 粒子1: 吸収された電子（重み 1.0）、2: 障壁を越えて脱出した電子（0.2）、3: 再吸収 PE（0.5）、4: 脱出 PE（0.5）。
    x = 0.0_dp
    v = 0.0_dp
    call init_particles(batch, x, v, [-1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
                        [1.0_dp, 0.2_dp, 0.5_dp, 0.5_dp], species_id=[1_i32, 1_i32, 2_i32, 2_i32], &
                        source_element=[0_i32, 0_i32, 1_i32, 1_i32])
    call workspace%init(1_i32, 2_i32, 1_i32)
    call prepare_reflected_fixture(workspace)
    call apply_fixed_surface_current_closure(app, model, batch, 4_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), -0.8_dp, 1.0e-12_dp, &
                         'electron target = root inflow minus the tracked reflected escape')
    call assert_close_dp(workspace%fixed_absorbed_weight_scale(1), 0.8_dp, 1.0e-12_dp, 'electron scale')
    call assert_close_dp(workspace%fixed_escape_target_charge(2), -0.3_dp, 1.0e-12_dp, &
                         'PE escape target = root escape minus the reflected escape')
    call assert_close_dp(workspace%fixed_absorbed_target_charge(2), -0.7_dp, 1.0e-12_dp, &
                         'PE return target grows by the reflected escape (emission minus escape)')
    call assert_close_dp(workspace%fixed_emission_target_charge(2), 1.0_dp, 1.0e-12_dp, 'PE emission stays the surface source')
    call assert_close_dp(workspace%fixed_escape_closed_charge(2), 0.2_dp, 1.0e-12_dp, 'reflected escape placed')
    call assert_close_dp(workspace%fixed_escape_unclosed_charge(2), 0.0_dp, 0.0_dp, 'nothing unclosed')
    ! 収支: 電子 −0.8 + イオン +0.5（根の J_i = J_e − J_esc）+ PE 正味 (1.0 − 0.7) = 0。
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1) + 0.5_dp + workspace%fixed_emission_target_charge(2) + &
                         workspace%fixed_absorbed_target_charge(2), 0.0_dp, 1.0e-12_dp, 'plate budget closes with the root')
    call assert_close_dp(workspace%fixed_surface_net_charge, -0.8_dp + 0.3_dp, 1.0e-12_dp, &
                         'surface net sums the electron and PE deposits')

    ! 追跡した反射電子の脱出 0.2 が根の脱出 0.1 を超えると、PE の脱出目標は 0 で止まり帰還が放出全量になる。
    ! 電子へは PE 側が実際に減らした 0.1 だけを足し戻すので、板の収支は閉じたまま（#93）。
    model%absorbed_current_a = [-1.0_dp, -0.9_dp]
    model%escaped_particle_current_a = [0.0_dp, -0.1_dp]
    call prepare_reflected_fixture(workspace)
    call apply_fixed_surface_current_closure(app, model, batch, 4_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_escape_target_charge(2), 0.0_dp, 1.0e-12_dp, 'PE escape target stops at zero')
    call assert_close_dp(workspace%fixed_absorbed_target_charge(2), -1.0_dp, 1.0e-12_dp, 'PE return is the whole emission')
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), -0.9_dp, 1.0e-12_dp, &
                         'electron is credited only with the reflected escape the PE targets took off')
    ! 収支: 電子 −0.9 + イオン +0.9（根の J_i = J_e − J_esc = 1.0 − 0.1）+ PE 正味 (1.0 − 1.0) = 0。
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1) + 0.9_dp + workspace%fixed_emission_target_charge(2) + &
                         workspace%fixed_absorbed_target_charge(2), 0.0_dp, 1.0e-12_dp, &
                         'plate budget closes when the reflected escape exceeds the root escape')

    ! 追跡した源の種（PE なしの SEE）では、反射電子の脱出・帰還・外向きを脱出・帰還・放出の目標から外す。
    model%emission_target_from_tracked_source(2) = .true.
    model%absorbed_current_a = [-1.0_dp, -0.1_dp]
    model%emission_current_a = [0.0_dp, 0.4_dp]
    model%escaped_particle_current_a = [0.0_dp, -0.3_dp]
    call prepare_reflected_fixture(workspace)
    call apply_fixed_surface_current_closure(app, model, batch, 4_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), -0.8_dp, 1.0e-12_dp, 'electron credit is unchanged')
    call assert_close_dp(workspace%fixed_escape_target_charge(2), -0.1_dp, 1.0e-12_dp, 'SEE escape = root escape minus reflected')
    call assert_close_dp(workspace%fixed_absorbed_target_charge(2), -0.05_dp, 1.0e-12_dp, &
                         'SEE return = root return minus the reflected barrier return')
    call assert_close_dp(workspace%fixed_emission_target_charge(2), 0.15_dp, 1.0e-12_dp, &
                         'SEE emission = root outward minus the reflected outward')
    call assert_close_dp(workspace%fixed_emission_target_charge(2) + workspace%fixed_absorbed_target_charge(2) + &
                         workspace%fixed_escape_target_charge(2), 0.0_dp, 1.0e-12_dp, 'SEE emission minus return equals escape')
    ! 収支: 電子 −0.8 + イオン +0.7 + SEE 正味 0.10 = 0。
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1) + 0.7_dp + workspace%fixed_emission_target_charge(2) + &
                         workspace%fixed_absorbed_target_charge(2), 0.0_dp, 1.0e-12_dp, 'SEE-only plate budget closes')

    ! 追跡した反射電子の脱出が根の流入を超えることはないが、超えれば吸収なしまで縮めて止まらない。
    model%absorbed_current_a = [-0.1_dp, -0.1_dp]
    call prepare_reflected_fixture(workspace)
    call apply_fixed_surface_current_closure(app, model, batch, 4_i32, workspace, mpi)
    call assert_close_dp(workspace%fixed_absorbed_target_charge(1), 0.0_dp, 0.0_dp, 'electron target clamps at zero')
  end subroutine test_fixed_current_reflected_electrons

  !> PE も SEE もなく反射電子だけを外部源に数えるとき、電子の吸収と脱出は根の値に固定され、追跡した脱出の
  !! 揺らぎで動かない。イオンと合わせて板の収支は毎バッチ閉じる。
  subroutine test_fixed_current_reflected_only()
    type(app_config) :: app
    type(surface_closure_contract_type) :: model
    type(particles_soa) :: batch
    type(simulator_batch_workspace_type) :: workspace
    type(mpi_context) :: mpi
    real(dp) :: x(3, 2), v(3, 2), tracked_escape
    integer :: k

    call default_app_config(app)
    app%sim%batch_duration = 1.0_dp
    app%n_particle_species = 2_i32
    app%particle_species(1) = species_from_defaults()
    app%particle_species(1)%species_key = 'electron'
    app%particle_species(1)%q_particle = -1.0_dp
    app%particle_species(2) = species_from_defaults()
    app%particle_species(2)%species_key = 'ion'
    app%particle_species(2)%q_particle = 1.0_dp

    model%active = .true.
    allocate (model%has_absorbed_target(2), model%has_emission_target(2), model%has_escape_target(2), &
              model%absorbed_current_a(2), model%emission_current_a(2), model%escaped_particle_current_a(2), &
              model%escape_partner_species(2), model%emission_target_from_tracked_source(2), &
              model%reflected_electron_species(2), model%absorbed_target_credits_tracked_escape(2))
    ! イオンの目標は根の値（+0.8）で、このバッチにイオン粒子は置かないので電子だけを検査する。
    model%has_absorbed_target = [.true., .false.]
    model%has_emission_target = .false.
    model%has_escape_target = [.true., .false.]
    model%escape_partner_species = 0_i32
    model%emission_target_from_tracked_source = .false.
    model%reflected_electron_species = 0_i32
    model%absorbed_target_credits_tracked_escape = .false.
    ! 根: 電子流入 1.0、反射電子の脱出 0.2、イオン 0.8（J=0）。電子の吸収目標は流入 − 脱出 = 0.8。
    model%absorbed_current_a = [-0.8_dp, 0.8_dp]
    model%emission_current_a = 0.0_dp
    model%escaped_particle_current_a = [-0.2_dp, 0.0_dp]

    x = 0.0_dp
    v = 0.0_dp
    call workspace%init(1_i32, 2_i32, 1_i32)
    ! 追跡した脱出が根と違っても（0.35 と 0.05）、吸収目標は根の値のまま。
    do k = 1, 2
      tracked_escape = merge(0.35_dp, 0.05_dp, k == 1)
      call init_particles(batch, x, v, [-1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp], [1.0_dp, tracked_escape], &
                          species_id=[1_i32, 1_i32], source_element=[0_i32, 0_i32])
      call workspace%reset_before_injection()
      call workspace%prepare_particle_flags(2_i32)
      workspace%absorbed_flag(1) = .true.
      workspace%absorbed_element(1) = 1_i32
      workspace%escaped_boundary_flag(2) = .true.
      call apply_fixed_surface_current_closure(app, model, batch, 2_i32, workspace, mpi)
      call assert_close_dp(workspace%fixed_absorbed_target_charge(1), -0.8_dp, 1.0e-12_dp, &
                           'electron target is the root inflow minus the root escape')
      call assert_close_dp(workspace%fixed_absorbed_weight_scale(1), 0.8_dp, 1.0e-12_dp, 'electron scale')
      call assert_close_dp(workspace%fixed_escape_target_charge(1), -0.2_dp, 1.0e-12_dp, 'electron escape = root escape')
      call assert_close_dp(workspace%fixed_escape_correction(1), -0.2_dp + tracked_escape, 1.0e-12_dp, &
                           'escape correction is root minus tracked')
      call assert_close_dp(workspace%fixed_absorbed_target_charge(1) + 0.8_dp, 0.0_dp, 1.0e-12_dp, &
                           'electron and ion targets close the plate budget')
    end do
  end subroutine test_fixed_current_reflected_only

  subroutine prepare_reflected_fixture(workspace)
    type(simulator_batch_workspace_type), intent(inout) :: workspace

    call workspace%reset_before_injection()
    call workspace%prepare_particle_flags(4_i32)
    workspace%absorbed_flag(1) = .true.
    workspace%absorbed_element(1) = 1_i32
    workspace%escaped_boundary_flag(2) = .true.
    workspace%absorbed_flag(3) = .true.
    workspace%absorbed_element(3) = 1_i32
    workspace%escaped_boundary_flag(4) = .true.
    workspace%photo_emission_dq(1, 2) = 1.0_dp
    ! 外部障壁で戻った電子 0.05 個分（raw weight）。
    workspace%matching_plane_moments_thread(3, 1, 1) = 0.05_dp
  end subroutine prepare_reflected_fixture

  subroutine prepare_secondary_only_fixture(workspace, absorbed, escaped, emission)
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    logical, intent(in) :: absorbed, escaped
    real(dp), intent(in) :: emission

    call workspace%reset_before_injection()
    call workspace%prepare_particle_flags(2_i32)
    workspace%absorbed_flag(1) = absorbed
    workspace%absorbed_element(1) = 1_i32
    workspace%escaped_boundary_flag(2) = escaped
    workspace%photo_emission_dq(1, 1) = emission
  end subroutine prepare_secondary_only_fixture
end program test_surface_charge_closure
