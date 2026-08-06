!> 再開用チェックポイントの保存/復元を検証するテスト。
program test_restart
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file, &
                         restart_rng_state_path, restart_macro_residual_path, validate_restart_contract, &
                         restart_contract_ok, restart_contract_mismatch
  use bem_output_writer, only: write_result_files
  use bem_periodic_checkpoint, only: maybe_write_periodic_checkpoint, resolve_latest_checkpoint_dir
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
  character(len=256) :: contract_message
  character(len=1024) :: rng_rank_path, residual_global_path
  character(len=1024) :: resolved_checkpoint_dir
  character(len=*), parameter :: out_dir = 'test_restart_tmp'

  call cleanup_restart_fixture(out_dir)
  call test_init(6)

  call test_begin('missing_checkpoint')
  call build_test_mesh(mesh)
  call load_restart_checkpoint('test_restart_missing', mesh, stats, has_restart)
  call assert_true(.not. has_restart, 'missing checkpoint should not be reported as restart')
  call assert_equal_i32(stats%batches, 0_i32, 'missing checkpoint should keep stats at defaults')
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
  mesh%q_elem = [1.0e-12_dp, -2.0e-12_dp]
  allocate (state%macro_residual(2), state%boundary_macro_residual(6, 2))
  state%macro_residual = [0.25_dp, 0.75_dp]
  state%boundary_macro_residual = 0.0_dp
  state%boundary_macro_residual(2, 1) = 0.125_dp
  state%boundary_macro_residual(5, 2) = 0.625_dp
  call write_result_files(out_dir, mesh, stats, cfg, charge_ledger=ledger)
  call write_rng_state_file(out_dir)
  call write_macro_residuals_file(out_dir, state)

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

  call build_test_mesh(mesh)
  call load_restart_checkpoint( &
    trim(resolved_checkpoint_dir), mesh, stats, has_restart, state, app=cfg, charge_ledger=restored_ledger &
    )
  call assert_true(has_restart, 'resolved periodic checkpoint should load')
  call assert_equal_i32(stats%batches, 8_i32, 'latest periodic checkpoint batch mismatch')
  call assert_allclose_1d(mesh%q_elem, [8.0e-12_dp, -3.0e-12_dp], 1.0e-24_dp, 'periodic charge mismatch')

  stats%batches = 12_i32
  call write_result_files(out_dir, mesh, stats, cfg, charge_ledger=ledger)
  call delete_file_if_exists(out_dir//'/rng_state.txt')
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot1', &
    'incomplete newer final output must not hide the last complete periodic checkpoint' &
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

  subroutine cleanup_restart_fixture(dir_path)
    character(len=*), intent(in) :: dir_path

    call delete_file_if_exists(trim(dir_path)//'/summary.txt')
    call delete_file_if_exists(trim(dir_path)//'/charges.csv')
    call delete_file_if_exists(trim(dir_path)//'/rng_state.txt')
    call delete_file_if_exists(trim(dir_path)//'/macro_residuals.csv')
    call delete_file_if_exists(trim(dir_path)//'/rng_state_rank00001.txt')
    call delete_file_if_exists(trim(dir_path)//'/macro_residuals_rank00001.csv')
    call delete_file_if_exists(trim(dir_path)//'/mesh_triangles.csv')
    call delete_file_if_exists(trim(dir_path)//'/mesh_sources.csv')
    call delete_file_if_exists(trim(dir_path)//'/charge_ledger.csv')
    call delete_file_if_exists(trim(dir_path)//'/checkpoint_latest.txt')
    call delete_file_if_exists(trim(dir_path)//'/checkpoint_latest.txt.tmp')
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
