!> 再開用チェックポイントの保存/復元を検証するテスト。
program test_restart
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file, &
                         restart_rng_state_path, restart_macro_residual_path, validate_restart_contract, &
                         restart_contract_ok, restart_contract_mismatch
  use bem_output_writer, only: write_result_files
  use bem_periodic_checkpoint, only: maybe_write_periodic_checkpoint, resolve_latest_checkpoint_dir
  use bem_checkpoint_contract, only: begin_checkpoint_publish, checkpoint_schema_version_current, &
                                     inspect_checkpoint_directory, publish_checkpoint_manifest
  use bem_app_config, only: app_config, default_app_config, species_from_defaults
  use bem_mpi, only: mpi_context
  use bem_charge_ledger, only: charge_ledger_type
  use bem_types, only: mesh_type, sim_stats, injection_state
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, assert_allclose_1d, &
                          delete_file_if_exists, ensure_directory, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(sim_stats) :: stats
  type(injection_state) :: state
  type(app_config) :: cfg, cfg_changed
  type(charge_ledger_type) :: ledger, restored_ledger
  type(mpi_context) :: mpi
  logical :: has_restart, exists
  integer(i32) :: contract_status
  real(dp) :: rng_probe
  character(len=256) :: contract_message
  character(len=1024) :: rng_rank_path, residual_global_path
  character(len=1024) :: resolved_checkpoint_dir
  character(len=*), parameter :: out_dir = 'test_restart_tmp'
  character(len=*), parameter :: matching_probe_log = 'test_restart_matching_probe.log'
  character(len=64) :: run_mode

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--matching-summary-probe') then
    call build_test_mesh(mesh)
    call load_restart_checkpoint(out_dir, mesh, stats, has_restart)
    stop 0
  end if
  if (trim(run_mode) == '--world-size-probe') then
    call build_test_mesh(mesh)
    call load_restart_checkpoint(out_dir, mesh, stats, has_restart, mpi_rank=0_i32, mpi_size=1_i32)
    stop 0
  end if

  call cleanup_restart_fixture(out_dir)
  call test_init(9)

  call test_begin('missing_checkpoint')
  call build_test_mesh(mesh)
  call load_restart_checkpoint('test_restart_missing', mesh, stats, has_restart)
  call assert_true(.not. has_restart, 'missing checkpoint should not be reported as restart')
  call assert_equal_i32(stats%batches, 0_i32, 'missing checkpoint should keep stats at defaults')
  call test_end()

  call test_begin('manifest_required_scalar_keys')
  call test_manifest_required_scalar_keys()
  call test_end()

  call test_begin('mpi_checkpoint_completeness_layers')
  call test_mpi_checkpoint_completeness_layers()
  call test_end()

  call test_begin('schema_v9_matching_keys_are_required_and_unique')
  call test_matching_summary_key_contract()
  call test_end()

  call test_begin('checkpoint_contract_and_ledger_restore')
  call default_app_config(cfg)
  cfg%n_particle_species = 2_i32
  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%species_key = 'electron'
  cfg%particle_species(2) = species_from_defaults()
  cfg%particle_species(2)%species_key = 'ion'
  stats = sim_stats()
  stats%processed_particles = 10_i64
  stats%absorbed = 7_i64
  stats%escaped = 3_i64
  stats%batches = 2_i32
  stats%last_rel_change = 1.0e-3_dp
  call ledger%init(2_i32)
  call ledger%reset(2_i32)
  ledger%surface_charge_before = 1.0_dp
  ledger%surface_charge_after = 2.0_dp
  ledger%injected_from_remote = [-3.0_dp, 4.0_dp]
  ledger%injected_count = [3_i64, 4_i64]
  ledger%neutral_return_correction = [-0.25_dp, 0.0_dp]
  ledger%neutral_return_weight_scale = [1.25_dp, 1.0_dp]
  ledger%neutral_return_unresolved_fraction = [0.2_dp, 0.0_dp]
  ledger%fixed_absorbed_target_charge = [-2.0_dp, 3.0_dp]
  ledger%fixed_absorbed_weight_scale = [2.0_dp, 1.5_dp]
  ledger%fixed_emission_target_charge = [4.0_dp, 0.0_dp]
  ledger%fixed_emission_weight_scale = [1.25_dp, 1.0_dp]
  ledger%fixed_escape_target_charge = [-1.5_dp, 0.0_dp]
  ledger%fixed_escape_correction = [-0.5_dp, 0.0_dp]
  ledger%fixed_current_correction = [-0.5_dp, 0.25_dp]
  mesh%q_elem = [1.0e-12_dp, -2.0e-12_dp]
  allocate (state%macro_residual(2), state%boundary_macro_residual(6, 2))
  state%macro_residual = [0.25_dp, 0.75_dp]
  state%boundary_macro_residual = 0.0_dp
  state%boundary_macro_residual(2, 1) = 0.125_dp
  state%boundary_macro_residual(5, 2) = 0.625_dp
  call write_result_files(out_dir, mesh, stats, cfg, charge_ledger=ledger)
  call write_rng_state_file(out_dir)
  call write_macro_residuals_file(out_dir, state)
  call publish_checkpoint_manifest(out_dir, stats%batches, 1_i32, .true., .true.)

  call build_test_mesh(mesh)
  state%macro_residual = 0.0_dp
  state%boundary_macro_residual = 0.0_dp
  call load_restart_checkpoint( &
    out_dir, mesh, stats, has_restart, state, app=cfg, charge_ledger=restored_ledger &
    )
  call assert_true(has_restart, 'checkpoint should load')
  call assert_equal_i32(restored_ledger%batch_count, 2_i32, 'ledger batch count mismatch')
  call assert_close_dp(restored_ledger%surface_charge_before, 1.0_dp, 1.0e-12_dp, 'ledger stock mismatch')
  call assert_allclose_1d( &
    restored_ledger%injected_from_remote, [-3.0_dp, 4.0_dp], 1.0e-12_dp, 'ledger flux restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%neutral_return_correction, [-0.25_dp, 0.0_dp], 1.0e-12_dp, &
    'neutral-return correction restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%neutral_return_weight_scale, [1.25_dp, 1.0_dp], 1.0e-12_dp, &
    'neutral-return scale restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%neutral_return_unresolved_fraction, [0.2_dp, 0.0_dp], 1.0e-12_dp, &
    'neutral-return unresolved fraction restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%fixed_absorbed_target_charge, [-2.0_dp, 3.0_dp], 1.0e-12_dp, &
    'fixed absorbed target restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%fixed_emission_target_charge, [4.0_dp, 0.0_dp], 1.0e-12_dp, &
    'fixed emission target restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%fixed_escape_target_charge, [-1.5_dp, 0.0_dp], 1.0e-12_dp, &
    'fixed escape target restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%fixed_escape_correction, [-0.5_dp, 0.0_dp], 1.0e-12_dp, &
    'fixed escape correction restore mismatch' &
    )
  call assert_allclose_1d( &
    restored_ledger%fixed_current_correction, [-0.5_dp, 0.25_dp], 1.0e-12_dp, &
    'fixed current correction restore mismatch' &
    )
  call assert_allclose_1d(mesh%q_elem, [1.0e-12_dp, -2.0e-12_dp], 1.0e-24_dp, 'charge restore mismatch')
  call assert_allclose_1d(state%macro_residual, [0.25_dp, 0.75_dp], 1.0e-12_dp, 'macro residual restore mismatch')
  call assert_close_dp( &
    state%boundary_macro_residual(2, 1), 0.125_dp, 1.0e-12_dp, 'boundary residual x-high restore mismatch' &
    )
  call assert_close_dp( &
    state%boundary_macro_residual(5, 2), 0.625_dp, 1.0e-12_dp, 'boundary residual z-low restore mismatch' &
    )

  cfg_changed = cfg
  cfg_changed%sim%dt = 2.0_dp*cfg%sim%dt
  call validate_restart_contract( &
    out_dir//'/summary.txt', mesh, cfg_changed, contract_status, contract_message &
    )
  call assert_equal_i32(contract_status, restart_contract_mismatch, 'model mismatch should be rejected')
  call validate_restart_contract(out_dir//'/summary.txt', mesh, cfg, contract_status, contract_message)
  call assert_equal_i32(contract_status, restart_contract_ok, 'matching contract should be accepted')
  call test_end()

  call test_begin('periodic_checkpoint_slots')
  cfg%output_dir = out_dir
  cfg%checkpoint_stride = 4_i32
  stats%batches = 4_i32
  stats%processed_particles = 20_i64
  ledger%batch_count = 4_i32
  mesh%q_elem = [4.0e-12_dp, -1.0e-12_dp]
  call maybe_write_periodic_checkpoint(cfg, mesh, stats, state, mpi, ledger)
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot0', &
    'first periodic checkpoint should publish slot0' &
    )

  stats%batches = 8_i32
  stats%processed_particles = 40_i64
  ledger%batch_count = 8_i32
  mesh%q_elem = [8.0e-12_dp, -3.0e-12_dp]
  call maybe_write_periodic_checkpoint(cfg, mesh, stats, state, mpi, ledger)
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot1', &
    'second periodic checkpoint should publish slot1' &
    )

  ! complete manifestを持つslotは、indexが欠落、破損、古い場合も回収する。
  call delete_file_if_exists(out_dir//'/checkpoint_latest.txt')
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot1', &
    'missing index must not hide a complete periodic checkpoint' &
    )
  call write_checkpoint_index_fixture(out_dir, -1_i32, -1_i32, malformed=.true.)
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot1', &
    'malformed index must not hide a complete periodic checkpoint' &
    )
  call write_checkpoint_index_fixture(out_dir, 0_i32, 4_i32)
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot1', &
    'stale index must not hide a newer complete periodic checkpoint' &
    )

  call build_test_mesh(mesh)
  call load_restart_checkpoint( &
    trim(resolved_checkpoint_dir), mesh, stats, has_restart, state, app=cfg, charge_ledger=restored_ledger &
    )
  call assert_true(has_restart, 'resolved periodic checkpoint should load')
  call assert_equal_i32(stats%batches, 8_i32, 'latest periodic checkpoint batch mismatch')
  call assert_allclose_1d(mesh%q_elem, [8.0e-12_dp, -3.0e-12_dp], 1.0e-24_dp, 'periodic charge mismatch')

  ! stale indexがslot0を指していても、最新slot1をactiveとみなしてslot0へ次世代を書く。
  stats%batches = 12_i32
  stats%processed_particles = 60_i64
  ledger%batch_count = 12_i32
  mesh%q_elem = [12.0e-12_dp, -4.0e-12_dp]
  call maybe_write_periodic_checkpoint(cfg, mesh, stats, state, mpi, ledger)
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot0', &
    'recovered latest slot must be preserved while the other slot is updated' &
    )

  ! loadできないfuture schemaのroot/slotは、より古いload可能なslotを隠さない。
  call begin_checkpoint_publish(out_dir//'/checkpoints/slot1')
  call write_structural_summary_fixture( &
    out_dir//'/checkpoints/slot1', checkpoint_schema_version_current + 1_i32, 100_i32, 1_i32 &
    )
  call publish_checkpoint_manifest( &
    out_dir//'/checkpoints/slot1', 100_i32, 1_i32, .false., .false. &
    )
  call inspect_checkpoint_directory(out_dir//'/checkpoints/slot1', exists)
  call assert_true(.not. exists, 'future schema must not be reported as a loadable complete checkpoint')
  call begin_checkpoint_publish(out_dir)
  call write_structural_summary_fixture( &
    out_dir, checkpoint_schema_version_current + 1_i32, 99_i32, 1_i32 &
    )
  call publish_checkpoint_manifest(out_dir, 99_i32, 1_i32, .false., .false.)
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot0', &
    'unsupported newer root or slot schema must not hide a loadable periodic checkpoint' &
    )

  stats%batches = 16_i32
  ledger%batch_count = 16_i32
  call write_result_files(out_dir, mesh, stats, cfg, charge_ledger=ledger)
  call write_rng_state_file(out_dir)
  call write_macro_residuals_file(out_dir, state)
  call publish_checkpoint_manifest(out_dir, stats%batches, 1_i32, .true., .true.)
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir, &
    'newer loadable complete final output must win over both periodic slots' &
    )

  ! 同じbatch番号の新しい世代を書き始め、旧RNG/residualが残った時点を再現する。
  mesh%q_elem = [16.0e-12_dp, -5.0e-12_dp]
  state%macro_residual = [0.5_dp, 0.25_dp]
  call random_number(rng_probe)
  call write_result_files(out_dir, mesh, stats, cfg, charge_ledger=ledger)
  inquire (file=out_dir//'/rng_state.txt', exist=exists)
  call assert_true(exists, 'mixed-generation fixture must retain a stale RNG state')
  inquire (file=out_dir//'/macro_residuals.csv', exist=exists)
  call assert_true(exists, 'mixed-generation fixture must retain a stale macro residual')
  inquire (file=out_dir//'/charge_ledger.csv', exist=exists)
  call assert_true(exists, 'mixed-generation fixture must retain its conditional charge ledger')
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot0', &
    'stale-but-present final state must not hide the last complete periodic checkpoint' &
    )

  call write_rng_state_file(out_dir)
  call write_macro_residuals_file(out_dir, state)
  call publish_checkpoint_manifest(out_dir, stats%batches, 1_i32, .true., .true.)
  call delete_file_if_exists(out_dir//'/macro_residuals.csv')
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot0', &
    'missing residual declared by manifest must fall back to the complete periodic checkpoint' &
    )
  call begin_checkpoint_publish(out_dir)
  call write_macro_residuals_file(out_dir, state)
  call publish_checkpoint_manifest(out_dir, stats%batches, 1_i32, .true., .true.)
  call delete_file_if_exists(out_dir//'/charge_ledger.csv')
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot0', &
    'missing ledger declared by summary must fall back to the complete periodic checkpoint' &
    )
  call test_end()

  call test_begin('write_checkpoint')
  call ensure_directory(out_dir)
  call write_summary_fixture(out_dir)
  call write_charges_fixture(out_dir)
  call write_rng_state_file(out_dir)
  state%macro_residual = [0.25_dp, 0.75_dp]
  call write_macro_residuals_file(out_dir, state)
  inquire (file=trim(out_dir)//'/rng_state.txt', exist=exists)
  call assert_true(exists, 'rng_state.txt should be created')
  inquire (file=trim(out_dir)//'/macro_residuals.csv', exist=exists)
  call assert_true(exists, 'macro_residuals.csv should be created')
  call test_end()

  call test_begin('ranked_paths')
  rng_rank_path = restart_rng_state_path(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  residual_global_path = restart_macro_residual_path(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  call assert_true(trim(rng_rank_path) == trim(out_dir)//'/rng_state_rank00001.txt', 'ranked rng path mismatch')
  call assert_true( &
    trim(residual_global_path) == trim(out_dir)//'/macro_residuals.csv', 'MPI residual path must remain global' &
    )
  call write_rng_state_file(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  call write_macro_residuals_file(out_dir, state, mpi_rank=1_i32, mpi_size=4_i32)
  inquire (file=trim(rng_rank_path), exist=exists)
  call assert_true(exists, 'rng_state_rank00001.txt should be created')
  inquire (file=trim(out_dir)//'/macro_residuals_rank00001.csv', exist=exists)
  call assert_true(.not. exists, 'non-root must not create a ranked macro residual file')
  call test_end()

  call test_begin('load_checkpoint')
  call build_test_mesh(mesh)
  state%macro_residual = 0.0_dp
  state%boundary_macro_residual = 0.0_dp
  call load_restart_checkpoint(out_dir, mesh, stats, has_restart, state)
  call assert_true(has_restart, 'complete checkpoint should be detected')
  call assert_equal_i64(stats%processed_particles, 10_i64, 'processed_particles mismatch')
  call assert_equal_i64(stats%absorbed, 7_i64, 'absorbed mismatch')
  call assert_equal_i64(stats%escaped, 3_i64, 'escaped mismatch')
  call assert_equal_i32(stats%batches, 2_i32, 'batches mismatch')
  call assert_equal_i64(stats%escaped_boundary, 1_i64, 'escaped_boundary mismatch')
  call assert_equal_i64(stats%survived_max_step, 2_i64, 'survived_max_step mismatch')
  call assert_equal_i64(stats%multiple_box_events_retry_attempted, 6_i64, 'retry attempted count mismatch')
  call assert_equal_i64(stats%multiple_box_events_retry_resolved, 2_i64, 'retry resolved count mismatch')
  call assert_equal_i64( &
    stats%multiple_box_events_soft_discarded, 4_i64, 'soft discarded count mismatch' &
    )
  call assert_equal_i32( &
    stats%adaptive_nonzero_mode_omp_threads, 6_i32, 'adaptive OpenMP thread count mismatch' &
    )
  call assert_close_dp( &
    stats%multiple_box_events_soft_discarded_abs_charge, 2.5e-15_dp, 1.0e-27_dp, &
    'soft discarded absolute charge mismatch' &
    )
  call assert_close_dp(stats%last_rel_change, 1.0e-3_dp, 1.0e-12_dp, 'last_rel_change mismatch')
  call assert_allclose_1d(mesh%q_elem, [1.0e-12_dp, -2.0e-12_dp], 1.0e-24_dp, 'charge restore mismatch')
  call assert_allclose_1d(state%macro_residual, [0.25_dp, 0.75_dp], 1.0e-12_dp, 'macro residual restore mismatch')
  call assert_close_dp( &
    state%boundary_macro_residual(2, 1), 0.125_dp, 1.0e-12_dp, 'boundary residual reload mismatch' &
    )
  call test_end()

  call cleanup_restart_fixture(out_dir)
  call test_summary()

contains

  subroutine test_manifest_required_scalar_keys()
    character(len=32), parameter :: required_keys(6) = [ &
                                    character(len=32) :: 'schema_version', 'state', 'batches', 'mpi_world_size', &
                                                         'macro_residuals_present', 'charge_ledger_present' &
                                                         ]
    integer :: key_idx

    call cleanup_restart_fixture(out_dir)
    call ensure_directory(out_dir)
    call write_structural_summary_fixture(out_dir, checkpoint_schema_version_current, 2_i32, 1_i32)
    call write_charges_fixture(out_dir)
    call write_rng_state_file(out_dir)
    call write_manifest_fixture_omitting(out_dir, '')
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(exists, 'handwritten complete manifest fixture was rejected')
    do key_idx = 1, size(required_keys)
      call write_manifest_fixture_omitting(out_dir, trim(required_keys(key_idx)))
      call inspect_checkpoint_directory(out_dir, exists)
      call assert_true(.not. exists, 'manifest accepted missing '//trim(required_keys(key_idx)))
    end do
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(exists, 'complete manifest fixture was rejected')
    call write_structural_summary_fixture( &
      out_dir, checkpoint_schema_version_current, 2_i32, 1_i32, omit_batches=.true. &
      )
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(.not. exists, 'summary metadata accepted a missing batches key')
    call cleanup_restart_fixture(out_dir)
  end subroutine test_manifest_required_scalar_keys

  subroutine write_manifest_fixture_omitting(dir_path, omitted_key)
    character(len=*), intent(in) :: dir_path, omitted_key
    integer :: u, ios

    open ( &
      newunit=u, file=trim(dir_path)//'/checkpoint_complete.txt', &
      status='replace', action='write', iostat=ios &
      )
    if (ios /= 0) error stop 'failed to open checkpoint manifest fixture'
    if (trim(omitted_key) /= 'schema_version') write (u, '(a)') 'schema_version=1'
    if (trim(omitted_key) /= 'state') write (u, '(a)') 'state=complete'
    if (trim(omitted_key) /= 'batches') write (u, '(a)') 'batches=2'
    if (trim(omitted_key) /= 'mpi_world_size') write (u, '(a)') 'mpi_world_size=1'
    if (trim(omitted_key) /= 'macro_residuals_present') write (u, '(a)') 'macro_residuals_present=F'
    if (trim(omitted_key) /= 'charge_ledger_present') write (u, '(a)') 'charge_ledger_present=F'
    close (u)
  end subroutine write_manifest_fixture_omitting

  subroutine test_mpi_checkpoint_completeness_layers()
    character(len=1024) :: executable_path, command

    call cleanup_restart_fixture(out_dir)
    call ensure_directory(out_dir)
    call write_matching_summary_fixture(out_dir, mpi_world_size=2_i32)
    call write_charges_fixture(out_dir)
    call write_rng_state_file(out_dir, mpi_rank=0_i32, mpi_size=2_i32)
    call write_rng_state_file(out_dir, mpi_rank=1_i32, mpi_size=2_i32)
    call publish_checkpoint_manifest(out_dir, 2_i32, 2_i32, .false., .false.)
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(exists, 'complete two-rank checkpoint fixture was rejected')

    call delete_file_if_exists(out_dir//'/charges.csv')
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(.not. exists, 'checkpoint inspection accepted missing charges.csv')
    call write_charges_fixture(out_dir)

    call delete_file_if_exists(out_dir//'/rng_state_rank00001.txt')
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(.not. exists, 'checkpoint inspection accepted a missing rank RNG state')
    call write_rng_state_file(out_dir, mpi_rank=1_i32, mpi_size=2_i32)

    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(.not. exists, 'checkpoint inspection accepted mismatched summary/manifest world sizes')
    call publish_checkpoint_manifest(out_dir, 2_i32, 2_i32, .false., .false.)
    call inspect_checkpoint_directory(out_dir, exists)
    call assert_true(exists, 'repaired two-rank checkpoint fixture was rejected')

    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" --world-size-probe > '//matching_probe_log//' 2>&1'
    call assert_restart_probe_fails( &
      command, 'mpi_world_size does not match current MPI world size', 'runtime world-size mismatch probe' &
      )
    call cleanup_restart_fixture(out_dir)
  end subroutine test_mpi_checkpoint_completeness_layers

  subroutine test_matching_summary_key_contract()
    character(len=1024) :: executable_path, command
    integer :: exit_status, command_status

    call cleanup_restart_fixture(out_dir)
    call ensure_directory(out_dir)
    call write_matching_summary_fixture(out_dir, omit_ion_access=.true.)
    call write_charges_fixture(out_dir)
    call write_rng_state_file(out_dir)
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" --matching-summary-probe > '//matching_probe_log//' 2>&1'
    call assert_restart_probe_fails( &
      command, 'schema-v9 summary is missing required matching-plane keys', 'missing matching key probe' &
      )

    call write_matching_summary_fixture(out_dir, duplicate_residual=.true.)
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call assert_restart_probe_fails( &
      command, 'duplicate matching-plane key: matching_plane_residual', 'duplicate matching key probe' &
      )

    call write_matching_summary_fixture(out_dir, ion_access_token='/')
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call assert_restart_probe_fails( &
      command, 'invalid real token: matching_plane_ion_access_potential_V', 'matching real-token probe' &
      )

    call write_matching_summary_fixture(out_dir, iteration_token='2*0')
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call assert_restart_probe_fails( &
      command, 'invalid integer token: matching_plane_iterations', 'matching integer-token probe' &
      )

    call write_matching_summary_fixture(out_dir, state_token='1')
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call assert_restart_probe_fails( &
      command, 'invalid logical token: matching_plane_state_valid', 'matching logical-token probe' &
      )

    call write_matching_summary_fixture(out_dir, inconsistent_budget=.true.)
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call assert_restart_probe_fails( &
      command, 'matching-plane photoelectron budget is inconsistent', 'matching PE-budget probe' &
      )

    call write_matching_summary_fixture(out_dir)
    call publish_checkpoint_manifest(out_dir, 2_i32, 1_i32, .false., .false.)
    call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=command_status)
    call assert_equal_i32(int(command_status, i32), 0_i32, 'valid matching summary probe command status mismatch')
    call assert_equal_i32(int(exit_status, i32), 0_i32, 'valid schema-v9 matching summary was rejected')
    call delete_file_if_exists(matching_probe_log)
    call cleanup_restart_fixture(out_dir)
  end subroutine test_matching_summary_key_contract

  subroutine assert_restart_probe_fails(command, expected_fragment, label)
    character(len=*), intent(in) :: command, expected_fragment, label
    integer :: exit_status, command_status
    logical :: saw_expected_failure

    exit_status = -1
    command_status = -1
    call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=command_status)
    call assert_equal_i32(int(command_status, i32), 0_i32, trim(label)//' command status mismatch')
    call assert_true(exit_status /= 0, trim(label)//' unexpectedly succeeded')
    call text_file_contains(matching_probe_log, expected_fragment, saw_expected_failure)
    call assert_true(saw_expected_failure, trim(label)//' did not reach its intended restart guard')
  end subroutine assert_restart_probe_fails

  subroutine text_file_contains(path, expected_fragment, found)
    character(len=*), intent(in) :: path, expected_fragment
    logical, intent(out) :: found
    character(len=1024) :: line
    integer :: u, ios

    found = .false.
    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read (u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, trim(expected_fragment)) > 0) then
        found = .true.
        exit
      end if
    end do
    close (u)
  end subroutine text_file_contains

  subroutine write_matching_summary_fixture( &
    dir_path, omit_ion_access, duplicate_residual, ion_access_token, iteration_token, state_token, inconsistent_budget, &
    mpi_world_size &
    )
    character(len=*), intent(in) :: dir_path
    logical, intent(in), optional :: omit_ion_access, duplicate_residual
    character(len=*), intent(in), optional :: ion_access_token, iteration_token, state_token
    logical, intent(in), optional :: inconsistent_budget
    integer(i32), intent(in), optional :: mpi_world_size
    integer :: u, ios
    integer(i32) :: world_size
    logical :: omit_access, duplicate_key, break_budget
    character(len=32) :: access_value, iterations_value, state_value

    omit_access = .false.
    duplicate_key = .false.
    break_budget = .false.
    access_value = '0.0'
    iterations_value = '2'
    state_value = 'T'
    world_size = 1_i32
    if (present(omit_ion_access)) omit_access = omit_ion_access
    if (present(duplicate_residual)) duplicate_key = duplicate_residual
    if (present(ion_access_token)) access_value = ion_access_token
    if (present(iteration_token)) iterations_value = iteration_token
    if (present(state_token)) state_value = state_token
    if (present(inconsistent_budget)) break_budget = inconsistent_budget
    if (present(mpi_world_size)) world_size = mpi_world_size
    open (newunit=u, file=trim(dir_path)//'/summary.txt', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open matching summary fixture'
    write (u, '(a,i0)') 'checkpoint_schema_version=', checkpoint_schema_version_current
    write (u, '(a)') 'mesh_nelem=2'
    write (u, '(a,i0)') 'mpi_world_size=', world_size
    write (u, '(a)') 'processed_particles=3'
    write (u, '(a)') 'absorbed=0'
    write (u, '(a)') 'escaped=1'
    write (u, '(a)') 'batches=2'
    write (u, '(a)') 'last_rel_change=0.0'
    write (u, '(a)') 'matching_plane_state_valid='//trim(state_value)
    write (u, '(a)') 'matching_plane_displacement_C_m2=0.0'
    write (u, '(a)') 'matching_plane_phi_V=1.0'
    write (u, '(a)') 'matching_plane_electron_inward_flux_m2_s=1.25'
    write (u, '(a)') 'matching_plane_ion_inward_flux_m2_s=1.25'
    write (u, '(a)') 'matching_plane_electron_access_potential_V=0.0'
    if (.not. omit_access) write (u, '(a)') 'matching_plane_ion_access_potential_V='//trim(access_value)
    write (u, '(a)') 'matching_plane_photoelectron_barrier_potential_V=0.0'
    write (u, '(a)') 'matching_plane_photoelectron_outward_flux_m2_s=1.0'
    write (u, '(a)') 'matching_plane_photoelectron_mean_normal_energy_eV=1.0'
    write (u, '(a)') 'matching_plane_electron_outward_flux_m2_s=0.0'
    write (u, '(a)') 'matching_plane_ion_outward_flux_m2_s=0.0'
    write (u, '(a)') 'matching_plane_photoelectron_return_flux_m2_s=0.0'
    if (break_budget) then
      write (u, '(a)') 'matching_plane_photoelectron_escape_flux_m2_s=0.5'
    else
      write (u, '(a)') 'matching_plane_photoelectron_escape_flux_m2_s=1.0'
    end if
    write (u, '(a)') 'matching_plane_iterations='//trim(iterations_value)
    write (u, '(a)') 'matching_plane_residual=0.0'
    if (duplicate_key) write (u, '(a)') 'matching_plane_residual=0.0'
    close (u)
  end subroutine write_matching_summary_fixture

  subroutine write_checkpoint_index_fixture(dir_path, slot, batches, malformed)
    character(len=*), intent(in) :: dir_path
    integer(i32), intent(in) :: slot, batches
    logical, intent(in), optional :: malformed
    integer :: u, ios
    logical :: write_malformed

    write_malformed = .false.
    if (present(malformed)) write_malformed = malformed
    open ( &
      newunit=u, file=trim(dir_path)//'/checkpoint_latest.txt', &
      status='replace', action='write', iostat=ios &
      )
    if (ios /= 0) error stop 'failed to open checkpoint index fixture'
    if (write_malformed) then
      write (u, '(a)') 'malformed checkpoint index fixture'
    else
      write (u, '(a)') 'schema_version=1'
      write (u, '(a,i0)') 'slot=', slot
      write (u, '(a,i0)') 'batches=', batches
    end if
    close (u)
  end subroutine write_checkpoint_index_fixture

  subroutine write_structural_summary_fixture(dir_path, schema_version, batches, mpi_world_size, omit_batches)
    character(len=*), intent(in) :: dir_path
    integer(i32), intent(in) :: schema_version, batches, mpi_world_size
    logical, intent(in), optional :: omit_batches
    integer :: u, ios
    logical :: skip_batches

    skip_batches = .false.
    if (present(omit_batches)) skip_batches = omit_batches

    open ( &
      newunit=u, file=trim(dir_path)//'/summary.txt', &
      status='replace', action='write', iostat=ios &
      )
    if (ios /= 0) error stop 'failed to open structural summary fixture'
    write (u, '(a,i0)') 'checkpoint_schema_version=', schema_version
    if (.not. skip_batches) write (u, '(a,i0)') 'batches=', batches
    write (u, '(a,i0)') 'mpi_world_size=', mpi_world_size
    close (u)
  end subroutine write_structural_summary_fixture

  subroutine cleanup_restart_fixture(dir_path)
    character(len=*), intent(in) :: dir_path

    call delete_file_if_exists(matching_probe_log)
    call delete_file_if_exists(trim(dir_path)//'/summary.txt')
    call delete_file_if_exists(trim(dir_path)//'/charges.csv')
    call delete_file_if_exists(trim(dir_path)//'/rng_state.txt')
    call delete_file_if_exists(trim(dir_path)//'/macro_residuals.csv')
    call delete_file_if_exists(trim(dir_path)//'/rng_state_rank00000.txt')
    call delete_file_if_exists(trim(dir_path)//'/rng_state_rank00001.txt')
    call delete_file_if_exists(trim(dir_path)//'/macro_residuals_rank00001.csv')
    call delete_file_if_exists(trim(dir_path)//'/mesh_triangles.csv')
    call delete_file_if_exists(trim(dir_path)//'/mesh_sources.csv')
    call delete_file_if_exists(trim(dir_path)//'/charge_ledger.csv')
    call delete_file_if_exists(trim(dir_path)//'/checkpoint_latest.txt')
    call delete_file_if_exists(trim(dir_path)//'/checkpoint_latest.txt.tmp')
    call delete_file_if_exists(trim(dir_path)//'/checkpoint_complete.txt')
    call delete_file_if_exists(trim(dir_path)//'/checkpoint_complete.txt.tmp')
    call cleanup_checkpoint_slot(trim(dir_path)//'/checkpoints/slot0')
    call cleanup_checkpoint_slot(trim(dir_path)//'/checkpoints/slot1')
    call remove_empty_directory(trim(dir_path)//'/checkpoints')
    call remove_empty_directory(dir_path)
  end subroutine cleanup_restart_fixture

  subroutine cleanup_checkpoint_slot(slot_dir)
    character(len=*), intent(in) :: slot_dir

    call delete_file_if_exists(trim(slot_dir)//'/summary.txt')
    call delete_file_if_exists(trim(slot_dir)//'/charges.csv')
    call delete_file_if_exists(trim(slot_dir)//'/rng_state.txt')
    call delete_file_if_exists(trim(slot_dir)//'/macro_residuals.csv')
    call delete_file_if_exists(trim(slot_dir)//'/charge_ledger.csv')
    call delete_file_if_exists(trim(slot_dir)//'/checkpoint_complete.txt')
    call delete_file_if_exists(trim(slot_dir)//'/checkpoint_complete.txt.tmp')
    call remove_empty_directory(slot_dir)
  end subroutine cleanup_checkpoint_slot

  !> 2要素メッシュを初期化する。
  subroutine build_test_mesh(mesh)
    type(mesh_type), intent(out) :: mesh
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    v0(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
    v1(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp]
    v2(:, 2) = [0.0_dp, 1.0_dp, 1.0_dp]
    call init_mesh(mesh, v0, v1, v2)
  end subroutine build_test_mesh

  !> `summary.txt` のフィクスチャを書き出す。
  subroutine write_summary_fixture(dir_path)
    character(len=*), intent(in) :: dir_path
    integer :: u, ios

    open (newunit=u, file=trim(dir_path)//'/summary.txt', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open summary fixture'
    write (u, '(a)') 'mesh_nelem=2'
    write (u, '(a)') 'processed_particles=10'
    write (u, '(a)') 'absorbed=7'
    write (u, '(a)') 'escaped=3'
    write (u, '(a)') 'batches=2'
    write (u, '(a)') 'escaped_boundary=1'
    write (u, '(a)') 'survived_max_step=2'
    write (u, '(a)') 'multiple_box_events_retry_attempted=6'
    write (u, '(a)') 'multiple_box_events_retry_resolved=2'
    write (u, '(a)') 'multiple_box_events_soft_discarded=4'
    write (u, '(a)') 'multiple_box_events_soft_discarded_abs_charge_C=2.5e-15'
    write (u, '(a)') 'adaptive_nonzero_mode_omp_threads=6'
    write (u, '(a)') 'last_rel_change=1.0e-3'
    close (u)
  end subroutine write_summary_fixture

  !> `charges.csv` のフィクスチャを書き出す。
  subroutine write_charges_fixture(dir_path)
    character(len=*), intent(in) :: dir_path
    integer :: u, ios

    open (newunit=u, file=trim(dir_path)//'/charges.csv', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open charges fixture'
    write (u, '(a)') 'elem_idx,charge_C'
    write (u, '(a)') '1,1.0e-12'
    write (u, '(a)') '2,-2.0e-12'
    close (u)
  end subroutine write_charges_fixture

end program test_restart
