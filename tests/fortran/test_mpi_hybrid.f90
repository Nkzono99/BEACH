!> MPI+OpenMPハイブリッド実行時の集約・rank別resumeファイルを検証する。
program test_mpi_hybrid
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: qe
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_is_root, mpi_select_lowest_rank_i32_values, &
                     mpi_allreduce_sum_i32_scalar, mpi_allreduce_sum_real_dp_array, &
                     mpi_allreduce_sum_i64_array, mpi_allreduce_min_i32_scalar, mpi_allreduce_max_i32_scalar, &
                     mpi_world_barrier
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_run_outputs, only: run_outputs_type
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            init_particle_batch_from_config, particle_inflow_reservoir
  use bem_restart, only: load_restart_checkpoint, write_macro_residuals_file, &
                         restart_rng_state_path, restart_macro_residual_path
  use bem_periodic_checkpoint, only: resolve_latest_checkpoint_dir
  use bem_app_config_runtime, only: particle_source_plan_type, build_particle_source_plan
  use bem_types, only: mesh_type, particles_soa, sim_stats, injection_state, bc_open, bc_reflect
  use bem_charge_ledger, only: charge_ledger_type
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_panel_plan, update_state, &
                                  eval_point, destroy_plan, destroy_state
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, &
                          assert_close_dp, delete_file_if_exists, &
                          ensure_directory, remove_empty_directory
  implicit none

  type(run_outputs_type) :: run_outputs

  type(mpi_context) :: mpi
  type(mesh_type) :: mesh, mesh_restart
  type(app_config) :: cfg, reservoir_cfg
  type(sim_stats) :: stats, stats_restart
  type(injection_state) :: state, state_restart, reservoir_state
  type(particles_soa) :: reservoir_particles
  type(particle_source_plan_type) :: reservoir_plan
  type(charge_ledger_type) :: ledger
  logical :: has_restart
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  integer :: u, ios
  character(len=256) :: line
  integer(i32) :: n_lines, history_batch, history_elem_idx
  integer(i32) :: selected_rank, expected_rank, batch_idx, global_reservoir_count
  integer(i64) :: history_processed_particles
  integer(i32) :: local_failure_values(4), selected_failure_values(4), reduced_min, reduced_max
  real(dp) :: history_rel_change, history_charge
  character(len=*), parameter :: history_path = 'test_mpi_hybrid_history_tmp.csv'
  character(len=*), parameter :: out_dir = 'test_mpi_hybrid_restart_tmp'
  character(len=1024) :: rng_path, residual_path
  integer(i32), parameter :: grid_nx = 6_i32, grid_ny = 6_i32, grid_nsrc = 2_i32*grid_nx*grid_ny
  type(fmm_options_type) :: fmm_options
  type(fmm_plan_type) :: fmm_plan
  type(fmm_state_type) :: replicated_state, distributed_state
  real(dp) :: grid_v0(3, grid_nsrc), grid_v1(3, grid_nsrc), grid_v2(3, grid_nsrc), grid_q(grid_nsrc)
  real(dp) :: replicated_field(3), distributed_field(3), probe(3)
  integer(i32) :: node_idx, coeff_idx, mismatches
  character(len=1024) :: periodic_rng_path, resolved_checkpoint_dir

  call mpi_initialize(mpi)

  if (mpi_is_root(mpi)) then
    call delete_file_if_exists(history_path)
    call delete_file_if_exists(out_dir//'/summary.txt')
    call delete_file_if_exists(out_dir//'/charges.csv')
    call delete_file_if_exists(out_dir//'/rng_state.txt')
    call delete_file_if_exists(out_dir//'/macro_residuals.csv')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt.tmp')
  end if
  call mpi_world_barrier(mpi)

  rng_path = restart_rng_state_path(out_dir, mpi=mpi)
  residual_path = restart_macro_residual_path(out_dir, mpi=mpi)
  call delete_file_if_exists(rng_path)
  if (mpi_is_root(mpi)) call delete_file_if_exists(residual_path)
  periodic_rng_path = restart_rng_state_path(out_dir//'/checkpoints/slot0', mpi=mpi)
  call delete_file_if_exists(periodic_rng_path)
  periodic_rng_path = restart_rng_state_path(out_dir//'/checkpoints/slot1', mpi=mpi)
  call delete_file_if_exists(periodic_rng_path)
  call mpi_world_barrier(mpi)
  if (mpi_is_root(mpi)) then
    call cleanup_periodic_slot(out_dir//'/checkpoints/slot0')
    call cleanup_periodic_slot(out_dir//'/checkpoints/slot1')
    call delete_file_if_exists(out_dir//'/checkpoint_latest.txt')
    call delete_file_if_exists(out_dir//'/checkpoint_latest.txt.tmp')
    call remove_empty_directory(out_dir//'/checkpoints')
  end if
  call mpi_world_barrier(mpi)

  v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
  v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  call init_mesh(mesh, v0, v1, v2)
  mesh%elem_vacuum_sign = 1_i32
  mesh%vacuum_normals = mesh%normals

  call default_app_config(cfg)
  cfg%sim%rng_seed = 2468_i32
  cfg%sim%batch_count = 1_i32
  cfg%sim%dt = 1.0d0
  cfg%sim%max_step = 1_i32
  cfg%sim%q_floor = 1.0d-30
  cfg%output_dir = out_dir
  cfg%checkpoint_stride = 1_i32
  cfg%sim%use_box = .true.
  cfg%sim%box_min = [-1.0d0, -1.0d0, -2.0d0]
  cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]

  cfg%n_particle_species = 1_i32
  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%source_mode = 'volume_seed'
  cfg%particle_species(1)%npcls_per_step = 4_i32
  cfg%particle_species(1)%q_particle = 1.0d0
  cfg%particle_species(1)%m_particle = 1.0d0
  cfg%particle_species(1)%w_particle = 1.0d0
  cfg%particle_species(1)%pos_low = [0.2d0, 0.2d0, 0.8d0]
  cfg%particle_species(1)%pos_high = [0.2d0, 0.2d0, 0.8d0]
  cfg%particle_species(1)%drift_velocity = [0.0d0, 0.0d0, -2.0d0]
  cfg%particle_species(1)%temperature_k = 0.0d0

  call seed_particles_from_config(cfg, mpi=mpi)

  call test_init(10)

  call test_begin('mpi_lowest_rank_metadata_selection')
  expected_rank = mpi%size - 1_i32
  ! photo collision failure payload: species, ray, bounce, status
  local_failure_values = 0_i32
  if (mpi%rank == expected_rank) local_failure_values = [2_i32, 17_i32, 3_i32, 1_i32]
  call mpi_select_lowest_rank_i32_values( &
    mpi, mpi%rank == expected_rank, local_failure_values, selected_rank, selected_failure_values &
    )
  call assert_equal_i32(selected_rank, expected_rank, 'single failure rank selection mismatch')
  call assert_true(all(selected_failure_values == [2_i32, 17_i32, 3_i32, 1_i32]), 'single failure metadata mismatch')

  local_failure_values = [1_i32, 20_i32 + mpi%rank, 4_i32, 2_i32]
  call mpi_select_lowest_rank_i32_values( &
    mpi, .true., local_failure_values, selected_rank, selected_failure_values &
    )
  call assert_equal_i32(selected_rank, 0_i32, 'lowest failure rank selection mismatch')
  call assert_true(all(selected_failure_values == [1_i32, 20_i32, 4_i32, 2_i32]), 'lowest rank metadata mismatch')
  reduced_min = mpi%rank + 1_i32
  reduced_max = reduced_min
  call mpi_allreduce_min_i32_scalar(mpi, reduced_min)
  call mpi_allreduce_max_i32_scalar(mpi, reduced_max)
  call assert_equal_i32(reduced_min, 1_i32, 'MPI i32 minimum reduction mismatch')
  call assert_equal_i32(reduced_max, mpi%size, 'MPI i32 maximum reduction mismatch')
  call test_end()

  call test_begin('mpi_simulation')
  if (mpi_is_root(mpi)) then
    open (newunit=u, file=history_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open MPI hybrid history fixture'
    run_outputs = run_outputs_type(history_unit=u, history_stride=1_i32)
    call run_absorption_insulator(mesh, cfg, stats, mpi=mpi, charge_ledger=ledger, outputs=run_outputs)
    close (u)
  else
    call run_absorption_insulator(mesh, cfg, stats, mpi=mpi, charge_ledger=ledger)
  end if

  call assert_equal_i64(stats%processed_particles, 4_i64, 'mpi processed_particles mismatch')
  call assert_equal_i64(stats%absorbed, 4_i64, 'mpi absorbed mismatch')
  call assert_equal_i64(stats%escaped, 0_i64, 'mpi escaped mismatch')
  call assert_equal_i32(stats%batches, 1_i32, 'mpi batches mismatch')
  call assert_close_dp(mesh%q_elem(1), 4.0d0, 1.0d-12, 'mpi deposited charge mismatch')
  call assert_equal_i64(ledger%injected_count(1), 4_i64, 'mpi ledger injected count mismatch')
  call assert_equal_i64(ledger%absorbed_count(1), 4_i64, 'mpi ledger absorbed count mismatch')
  call assert_close_dp(ledger%injected_from_remote(1), 4.0_dp, 1.0e-12_dp, 'mpi ledger injected charge mismatch')
  call assert_close_dp(ledger%residual(), 0.0_dp, 1.0e-12_dp, 'mpi ledger residual mismatch')
  call resolve_latest_checkpoint_dir(out_dir, resolved_checkpoint_dir)
  call assert_true( &
    trim(resolved_checkpoint_dir) == out_dir//'/checkpoints/slot0', &
    'MPI periodic checkpoint should publish slot0' &
    )
  call init_mesh(mesh_restart, v0, v1, v2)
  mesh_restart%elem_vacuum_sign = 1_i32
  mesh_restart%vacuum_normals = mesh_restart%normals
  call load_restart_checkpoint( &
    trim(resolved_checkpoint_dir), mesh_restart, stats_restart, has_restart, mpi=mpi, app=cfg &
    )
  call assert_true(has_restart, 'MPI periodic checkpoint should load')
  call assert_equal_i32(stats_restart%batches, 1_i32, 'MPI periodic checkpoint batch mismatch')
  call assert_close_dp(mesh_restart%q_elem(1), 4.0_dp, 1.0e-12_dp, 'MPI periodic checkpoint charge mismatch')
  call test_end()

  call test_begin('mpi_soft_discard_global_aggregation')
  call run_mpi_soft_discard_aggregation_test(mpi)
  call test_end()

  call test_begin('mpi_neutral_return_layout_invariance')
  call run_mpi_neutral_return_layout_test(mpi)
  call test_end()

  call test_begin('mpi_shared_pool_rank_independent')
  call run_mpi_shared_pool_test(mpi, 'shared')
  call test_end()

  call test_begin('mpi_static_slots_rank_independent')
  call run_mpi_shared_pool_test(mpi, 'static')
  call test_end()

  call test_begin('mpi_history')
  if (mpi_is_root(mpi)) then
    n_lines = 0_i32
    open (newunit=u, file=history_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to read MPI hybrid history fixture'
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      n_lines = n_lines + 1_i32
      read (line, *, iostat=ios) &
        history_batch, history_processed_particles, history_rel_change, history_elem_idx, history_charge
      if (ios /= 0) error stop 'failed to parse MPI hybrid history fixture'
    end do
    close (u)
    call assert_equal_i32(n_lines, 1_i32, 'mpi history snapshot line count mismatch')
    call assert_equal_i32(history_batch, 1_i32, 'MPI history batch mismatch')
    call assert_equal_i64(history_processed_particles, 4_i64, 'MPI history must contain the global processed count')
    call assert_true(history_rel_change >= 0.0_dp, 'MPI history relative change must be nonnegative')
    call assert_equal_i32(history_elem_idx, 1_i32, 'MPI history element index mismatch')
    call assert_close_dp(history_charge, 4.0_dp, 1.0e-12_dp, 'MPI history must contain globally reduced charge')
    call delete_file_if_exists(history_path)
  end if
  call mpi_world_barrier(mpi)
  call test_end()

  call test_begin('mpi_global_reservoir_count')
  call default_app_config(reservoir_cfg)
  reservoir_cfg%sim%batch_count = 4_i32
  reservoir_cfg%sim%batch_duration = 1.0_dp
  reservoir_cfg%sim%has_batch_duration = .true.
  reservoir_cfg%sim%use_box = .true.
  reservoir_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  reservoir_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  reservoir_cfg%n_particle_species = 1_i32
  reservoir_cfg%particle_species(1) = species_from_defaults()
  reservoir_cfg%particle_species(1)%source_mode = 'volume_seed'
  reservoir_cfg%particle_species(1)%npcls_per_step = 0_i32
  reservoir_cfg%particle_species(1)%number_density_m3 = 1.0_dp
  reservoir_cfg%particle_species(1)%temperature_k = 0.0_dp
  reservoir_cfg%particle_species(1)%q_particle = 0.0_dp
  reservoir_cfg%particle_species(1)%m_particle = 1.0_dp
  reservoir_cfg%particle_species(1)%w_particle = 4.0_dp
  reservoir_cfg%particle_species(1)%has_w_particle = .true.
  reservoir_cfg%particle_species(1)%boundary_inflow_low(3) = particle_inflow_reservoir
  reservoir_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, 1.0_dp]
  allocate (reservoir_state%macro_residual(1), reservoir_state%boundary_macro_residual(6, 1))
  reservoir_state%macro_residual = 0.0_dp
  reservoir_state%boundary_macro_residual = 0.0_dp
  call build_particle_source_plan(reservoir_cfg, reservoir_plan, mpi=mpi)
  do batch_idx = 1_i32, 4_i32
    call init_particle_batch_from_config( &
      reservoir_cfg, batch_idx, reservoir_particles, state=reservoir_state, mpi=mpi, source_plan=reservoir_plan &
      )
    global_reservoir_count = reservoir_particles%n
    call mpi_allreduce_sum_i32_scalar(mpi, global_reservoir_count)
    call assert_equal_i32( &
      global_reservoir_count, merge(1_i32, 0_i32, batch_idx == 4_i32), &
      'MPI global reservoir count sequence mismatch' &
      )
    call assert_close_dp( &
      reservoir_state%macro_residual(1), 0.0_dp, 1.0e-15_dp, &
      'MPI source residual must remain zero' &
      )
    call assert_close_dp( &
      reservoir_state%boundary_macro_residual(5, 1), modulo(0.25_dp*real(batch_idx, dp), 1.0_dp), 1.0e-15_dp, &
      'MPI global boundary reservoir residual mismatch' &
      )
  end do
  call test_end()

  call test_begin('mpi_restart')
  call ensure_directory(out_dir)
  if (mpi_is_root(mpi)) then
    call write_summary_fixture(out_dir, mpi%size)
    call write_charges_fixture(out_dir)
  end if
  call mpi_world_barrier(mpi)

  allocate (state%macro_residual(1), state%boundary_macro_residual(6, 1))
  state%macro_residual(1) = 0.25d0 + 0.01d0*real(mpi%rank, dp)
  state%boundary_macro_residual = 0.0_dp
  state%boundary_macro_residual(2, 1) = 0.5_dp + 0.01_dp*real(mpi%rank, dp)
  ! The fixture summary has no schema key, so it is a legacy checkpoint that still lists one RNG file per rank.
  rng_path = restart_rng_state_path(out_dir, mpi=mpi)
  open (newunit=u, file=trim(rng_path), status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open legacy MPI RNG fixture'
  write (u, '(a)') '0'
  close (u)
  call write_macro_residuals_file(out_dir, state, mpi=mpi)
  call mpi_world_barrier(mpi)

  call init_mesh(mesh_restart, v0, v1, v2)
  mesh_restart%elem_vacuum_sign = 1_i32
  mesh_restart%vacuum_normals = mesh_restart%normals
  allocate (state_restart%macro_residual(1), state_restart%boundary_macro_residual(6, 1))
  state_restart%macro_residual = 0.0d0
  state_restart%boundary_macro_residual = 0.0_dp
  call load_restart_checkpoint(out_dir, mesh_restart, stats_restart, has_restart, state_restart, mpi=mpi)
  call assert_true(has_restart, 'mpi restart should be detected')
  call assert_equal_i64(stats_restart%processed_particles, 8_i64, 'mpi restart processed_particles mismatch')
  call assert_equal_i32(stats_restart%batches, 2_i32, 'mpi restart batches mismatch')
  call assert_close_dp(mesh_restart%q_elem(1), 2.0d0, 1.0d-12, 'mpi restart charge mismatch')
  call assert_close_dp(state_restart%macro_residual(1), 0.25d0, 1.0d-12, 'MPI residual must be restored from root')
  call assert_close_dp( &
    state_restart%boundary_macro_residual(2, 1), 0.5_dp, 1.0e-12_dp, &
    'MPI boundary residual must be restored from root' &
    )
  call test_end()

  call delete_file_if_exists(rng_path)
  if (mpi_is_root(mpi)) call delete_file_if_exists(residual_path)
  periodic_rng_path = restart_rng_state_path(out_dir//'/checkpoints/slot0', mpi=mpi)
  call delete_file_if_exists(periodic_rng_path)
  call mpi_world_barrier(mpi)
  if (mpi_is_root(mpi)) then
    call delete_file_if_exists(out_dir//'/summary.txt')
    call delete_file_if_exists(out_dir//'/charges.csv')
    call delete_file_if_exists(out_dir//'/rng_state.txt')
    call delete_file_if_exists(out_dir//'/macro_residuals.csv')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt.tmp')
    call delete_file_if_exists(out_dir//'/checkpoint_latest.txt')
    call delete_file_if_exists(out_dir//'/checkpoint_latest.txt.tmp')
    call delete_file_if_exists(out_dir//'/checkpoints/slot0/summary.txt')
    call delete_file_if_exists(out_dir//'/checkpoints/slot0/charges.csv')
    call delete_file_if_exists(out_dir//'/checkpoints/slot0/charge_ledger.csv')
    call delete_file_if_exists(out_dir//'/checkpoints/slot0/checkpoint_complete.txt')
    call delete_file_if_exists(out_dir//'/checkpoints/slot0/checkpoint_complete.txt.tmp')
    call remove_empty_directory(out_dir//'/checkpoints/slot0')
    call remove_empty_directory(out_dir//'/checkpoints')
    call remove_empty_directory(out_dir)
  end if

  call test_begin('mpi_distributed_m2l_matches_replicated_update')
  ! update_state with the MPI context splits the M2L targets over ranks and gathers the local
  ! expansions; every rank must end with the same expansions as the replicated update, bit for bit.
  call make_mpi_panel_grid(grid_v0, grid_v1, grid_v2, grid_q)
  fmm_options%theta = 0.45_dp
  fmm_options%leaf_max = 4_i32
  fmm_options%order = 4_i32
  fmm_options%target_box_min = [-1.2_dp, -1.2_dp, -0.5_dp]
  fmm_options%target_box_max = [1.2_dp, 1.2_dp, 2.5_dp]
  call build_panel_plan(fmm_plan, grid_v0, grid_v1, grid_v2, fmm_options)
  call update_state(fmm_plan, replicated_state, grid_q)
  call update_state(fmm_plan, distributed_state, grid_q, mpi)
  call assert_true(fmm_plan%m2l_pair_count > 0_i32, 'MPI M2L test needs a plan with M2L pairs')
  mismatches = 0_i32
  do node_idx = 1_i32, size(replicated_state%local, 2)
    if (distributed_state%local_active(node_idx) /= replicated_state%local_active(node_idx)) mismatches = mismatches + 1_i32
    do coeff_idx = 1_i32, size(replicated_state%local, 1)
      if (distributed_state%local(coeff_idx, node_idx) /= replicated_state%local(coeff_idx, node_idx)) then
        mismatches = mismatches + 1_i32
      end if
    end do
  end do
  call assert_equal_i32(mismatches, 0_i32, 'distributed M2L local expansions differ from the replicated ones')
  probe = [0.42_dp, 0.31_dp, 0.45_dp]
  call eval_point(fmm_plan, replicated_state, probe, replicated_field)
  call eval_point(fmm_plan, distributed_state, probe, distributed_field)
  call assert_true(all(distributed_field == replicated_field), 'distributed M2L field differs from the replicated one')
  call destroy_state(replicated_state)
  call destroy_state(distributed_state)
  call destroy_plan(fmm_plan)
  call test_end()

  call test_summary()

  call mpi_shutdown(mpi)

contains

  !> Two triangles per cell on a square grid in z=0, with smoothly varying charges.
  subroutine make_mpi_panel_grid(v0_out, v1_out, v2_out, q_out)
    real(dp), intent(out) :: v0_out(3, grid_nsrc), v1_out(3, grid_nsrc), v2_out(3, grid_nsrc), q_out(grid_nsrc)
    integer(i32) :: ix, iy, element
    real(dp) :: x0, x1, y0, y1, dx, dy

    dx = 1.6_dp/real(grid_nx, dp)
    dy = 1.6_dp/real(grid_ny, dp)
    element = 0_i32
    do iy = 1_i32, grid_ny
      y0 = -0.8_dp + real(iy - 1_i32, dp)*dy
      y1 = y0 + dy
      do ix = 1_i32, grid_nx
        x0 = -0.8_dp + real(ix - 1_i32, dp)*dx
        x1 = x0 + dx
        element = element + 1_i32
        v0_out(:, element) = [x0, y0, 0.0_dp]
        v1_out(:, element) = [x1, y0, 0.0_dp]
        v2_out(:, element) = [x1, y1, 0.0_dp]
        q_out(element) = 1.0_dp + 0.2_dp*sin(real(element, dp))
        element = element + 1_i32
        v0_out(:, element) = [x0, y0, 0.0_dp]
        v1_out(:, element) = [x1, y1, 0.0_dp]
        v2_out(:, element) = [x0, y1, 0.0_dp]
        q_out(element) = 1.0_dp + 0.2_dp*cos(real(element, dp))
      end do
    end do
  end subroutine make_mpi_panel_grid

  subroutine run_mpi_soft_discard_aggregation_test(mpi)
    type(mpi_context), intent(in) :: mpi
    type(mesh_type) :: soft_mesh
    type(app_config) :: soft_cfg
    type(sim_stats) :: soft_stats
    real(dp) :: soft_v0(3, 1), soft_v1(3, 1), soft_v2(3, 1)

    soft_v0(:, 1) = [10.0_dp, -1.0_dp, -1.0_dp]
    soft_v1(:, 1) = [10.0_dp, 1.0_dp, -1.0_dp]
    soft_v2(:, 1) = [10.0_dp, 0.0_dp, 1.0_dp]
    call init_mesh(soft_mesh, soft_v0, soft_v1, soft_v2)
    soft_mesh%elem_vacuum_sign = 1_i32
    soft_mesh%vacuum_normals = soft_mesh%normals
    call default_app_config(soft_cfg)
    soft_cfg%sim%rng_seed = 5119_i32
    soft_cfg%sim%batch_count = 1_i32
    soft_cfg%sim%dt = 1.0_dp
    soft_cfg%sim%max_step = 1_i32
    soft_cfg%sim%use_box = .true.
    soft_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    soft_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    soft_cfg%sim%bc_low(1) = bc_reflect
    soft_cfg%sim%bc_high(1) = bc_reflect
    soft_cfg%sim%multiple_box_events_policy = 'soft_discard'
    soft_cfg%sim%multiple_box_events_retry_backend = 'upper_panel_fourier'
    soft_cfg%sim%multiple_box_events_soft_discard_count_grace = 100000_i32
    soft_cfg%sim%multiple_box_events_soft_discard_fraction_limit = 1.0_dp
    soft_cfg%sim%multiple_box_events_soft_discard_abs_charge_limit = 1.0e9_dp
    soft_cfg%n_particle_species = 1_i32
    soft_cfg%particle_species(1) = species_from_defaults()
    soft_cfg%particle_species(1)%source_mode = 'volume_seed'
    soft_cfg%particle_species(1)%npcls_per_step = mpi%size
    soft_cfg%particle_species(1)%q_particle = 2.0_dp
    soft_cfg%particle_species(1)%m_particle = 1.0_dp
    soft_cfg%particle_species(1)%w_particle = 3.0_dp
    soft_cfg%particle_species(1)%pos_low = [0.9_dp, 0.2_dp, 0.2_dp]
    soft_cfg%particle_species(1)%pos_high = soft_cfg%particle_species(1)%pos_low
    soft_cfg%particle_species(1)%drift_velocity = [9.0_dp, 0.0_dp, 0.0_dp]
    soft_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(soft_cfg, mpi=mpi)

    call run_absorption_insulator(soft_mesh, soft_cfg, soft_stats, mpi=mpi)
    call assert_equal_i64( &
      soft_stats%multiple_box_events_retry_attempted, int(mpi%size, i64), &
      'MPI retry attempt count was not globally reduced' &
      )
    call assert_equal_i64( &
      soft_stats%multiple_box_events_retry_resolved, 0_i64, &
      'ineligible MPI retry unexpectedly resolved' &
      )
    call assert_equal_i64( &
      soft_stats%multiple_box_events_soft_discarded, int(mpi%size, i64), &
      'MPI soft-discard count was not globally reduced' &
      )
    call assert_equal_i64( &
      soft_stats%processed_particles, int(mpi%size, i64), &
      'MPI soft-discard fraction denominator was not globally reduced' &
      )
    call assert_close_dp( &
      soft_stats%multiple_box_events_soft_discarded_abs_charge, 6.0_dp*real(mpi%size, dp), 1.0e-12_dp, &
      'MPI soft-discard absolute charge was not globally reduced' &
      )
  end subroutine run_mpi_soft_discard_aggregation_test

  !> 同じglobal ray集合を追跡するMPI実行とserial参照でneutral-return閉包を比較する。
  !!
  !! バッチの ray は global slot の counter-based 乱数で決まるので、MPI 側で rank に分けた ray 集合は
  !! serial 参照の ray 集合と同じになる。全rankの物理的な放出・帰還・未解決電荷とmesh depositはserial参照と一致する。
  subroutine run_mpi_neutral_return_layout_test(mpi)
    type(mpi_context), intent(in) :: mpi
    integer(i32), parameter :: rays_per_rank = 1024_i32
    integer(i32), parameter :: signature_size = 13_i32
    type(mesh_type) :: reference_mesh, distributed_mesh
    type(app_config) :: reference_cfg, distributed_cfg
    type(sim_stats) :: reference_stats, distributed_stats
    type(charge_ledger_type) :: reference_ledger, distributed_ledger
    real(dp) :: reference_signature(signature_size), distributed_signature(signature_size)
    integer(i64) :: reference_counts(3)
    integer(i32) :: value_index
    real(dp) :: tolerance
    character(len=128) :: assertion_message

    reference_signature = 0.0_dp
    reference_counts = 0_i64
    if (mpi_is_root(mpi)) then
      ! Batch slots are counter-based, so the root alone traces exactly the global ray set of the MPI run.
      call setup_mpi_neutral_return_fixture(reference_mesh, reference_cfg, rays_per_rank*mpi%size)
      call seed_particles_from_config(reference_cfg)
      call run_absorption_insulator( &
        reference_mesh, reference_cfg, reference_stats, charge_ledger=reference_ledger &
        )
      call pack_neutral_return_signature(reference_mesh, reference_ledger, reference_signature)
      reference_counts = [ &
                         reference_ledger%emitted_count(1), reference_ledger%absorbed_count(1), &
                         reference_ledger%discarded_unresolved_count(1) &
                         ]
    end if
    call mpi_allreduce_sum_real_dp_array(mpi, reference_signature)
    call mpi_allreduce_sum_i64_array(mpi, reference_counts)

    call setup_mpi_neutral_return_fixture( &
      distributed_mesh, distributed_cfg, rays_per_rank*mpi%size &
      )
    call seed_particles_from_config(distributed_cfg, mpi=mpi)
    call run_absorption_insulator( &
      distributed_mesh, distributed_cfg, distributed_stats, mpi=mpi, &
      charge_ledger=distributed_ledger &
      )
    call pack_neutral_return_signature(distributed_mesh, distributed_ledger, distributed_signature)

    do value_index = 1_i32, signature_size
      tolerance = 1.0e-11_dp*max( &
                  1.0_dp, abs(reference_signature(value_index)), &
                  abs(distributed_signature(value_index)) &
                  )
      write (assertion_message, '(a,i0)') &
        'neutral-return MPI signature mismatch at component ', value_index
      call assert_close_dp( &
        distributed_signature(value_index), reference_signature(value_index), tolerance, &
        trim(assertion_message) &
        )
    end do

    call assert_equal_i64( &
      distributed_ledger%emitted_count(1), reference_counts(1), &
      'neutral-return MPI emitted count was not reduced across ranks' &
      )
    call assert_equal_i64( &
      distributed_ledger%absorbed_count(1), reference_counts(2), &
      'neutral-return MPI absorbed count was not reduced across ranks' &
      )
    call assert_equal_i64( &
      distributed_ledger%discarded_unresolved_count(1), reference_counts(3), &
      'neutral-return MPI unresolved count was not reduced across ranks' &
      )
    call assert_true( &
      distributed_ledger%absorbed_count(1) > 0_i64 .and. &
      distributed_ledger%discarded_unresolved_count(1) > 0_i64, &
      'neutral-return MPI fixture must exercise resolved and unresolved returns' &
      )
    call assert_true( &
      distributed_ledger%neutral_return_unresolved_fraction(1) <= 0.05_dp, &
      'neutral-return MPI fixture must stay inside the fixed applicability limit' &
      )
    call assert_equal_i64( &
      distributed_stats%escaped_boundary, 0_i64, &
      'neutral-return MPI reflected photoelectrons must not escape' &
      )
    call assert_close_dp( &
      sum(distributed_mesh%q_elem), 0.0_dp, 1.0e-11_dp, &
      'neutral-return MPI mean surface charge must close' &
      )
    call assert_true( &
      sum(abs(distributed_mesh%q_elem)) > 0.0_dp, &
      'neutral-return MPI closure must preserve local redistribution' &
      )
    call assert_close_dp( &
      distributed_ledger%residual(), 0.0_dp, 1.0e-11_dp, &
      'neutral-return MPI ledger residual must close' &
      )
  end subroutine run_mpi_neutral_return_layout_test

  !> shared pool: 全 rank で slot を取り合った結果が、root だけで同じ pool を処理した結果と一致する。
  !! さらにその checkpoint を別の rank 数で読み戻せる。
  subroutine run_mpi_shared_pool_test(mpi, distribution)
    type(mpi_context), intent(in) :: mpi
    character(len=*), intent(in) :: distribution
    integer(i32), parameter :: rays_per_batch = 2048_i32
    integer(i32), parameter :: signature_size = 13_i32
    character(len=*), parameter :: pool_dir = 'test_mpi_hybrid_pool_tmp'
    type(mesh_type) :: reference_mesh, pool_mesh, resume_mesh
    type(app_config) :: reference_cfg, pool_cfg
    type(sim_stats) :: reference_stats, pool_stats, resume_stats
    type(charge_ledger_type) :: reference_ledger, pool_ledger
    real(dp) :: reference_signature(signature_size), pool_signature(signature_size), tolerance
    integer(i64) :: reference_counts(4), pool_counts(4)
    integer(i32) :: value_index
    logical :: pool_has_restart
    character(len=128) :: assertion_message
    character(len=1024) :: pool_checkpoint_dir

    call cleanup_pool_dir(mpi, pool_dir)
    if (mpi_is_root(mpi)) call ensure_directory(pool_dir)
    call mpi_world_barrier(mpi)

    reference_signature = 0.0_dp
    reference_counts = 0_i64
    if (mpi_is_root(mpi)) then
      call setup_mpi_neutral_return_fixture(reference_mesh, reference_cfg, rays_per_batch)
      reference_cfg%sim%particle_distribution = distribution
      reference_cfg%sim%particle_pool_chunk = 5_i32
      call seed_particles_from_config(reference_cfg)
      call run_absorption_insulator( &
        reference_mesh, reference_cfg, reference_stats, charge_ledger=reference_ledger &
        )
      call pack_neutral_return_signature(reference_mesh, reference_ledger, reference_signature)
      reference_counts = [ &
                         reference_ledger%emitted_count(1), reference_ledger%absorbed_count(1), &
                         reference_ledger%discarded_unresolved_count(1), reference_stats%processed_particles &
                         ]
    end if
    call mpi_allreduce_sum_real_dp_array(mpi, reference_signature)
    call mpi_allreduce_sum_i64_array(mpi, reference_counts)

    call setup_mpi_neutral_return_fixture(pool_mesh, pool_cfg, rays_per_batch)
    pool_cfg%sim%particle_distribution = distribution
    pool_cfg%sim%particle_pool_chunk = 5_i32
    pool_cfg%output_dir = pool_dir
    pool_cfg%checkpoint_stride = 1_i32
    call seed_particles_from_config(pool_cfg, mpi=mpi)
    call run_absorption_insulator(pool_mesh, pool_cfg, pool_stats, mpi=mpi, charge_ledger=pool_ledger)
    call pack_neutral_return_signature(pool_mesh, pool_ledger, pool_signature)
    pool_counts = [pool_ledger%emitted_count(1), pool_ledger%absorbed_count(1), pool_ledger%discarded_unresolved_count(1), &
                   pool_stats%processed_particles]

    do value_index = 1_i32, signature_size
      tolerance = 1.0e-11_dp*max(1.0_dp, abs(reference_signature(value_index)))
      write (assertion_message, '(a,i0)') 'shared pool signature must not depend on the rank count: index ', value_index
      call assert_close_dp(pool_signature(value_index), reference_signature(value_index), tolerance, &
                           trim(assertion_message))
    end do
    do value_index = 1_i32, 4_i32
      write (assertion_message, '(a,i0)') 'shared pool count must not depend on the rank count: index ', value_index
      call assert_equal_i64(pool_counts(value_index), reference_counts(value_index), trim(assertion_message))
    end do
    if (distribution == 'shared') then
      call assert_true(pool_stats%particle_pool_claims >= int(mpi%size, i64), 'every rank should claim pool slots')
    end if

    ! The pool keys every draw on the global slot, so the checkpoint may be resumed with another rank count.
    call resolve_latest_checkpoint_dir(pool_dir, pool_checkpoint_dir)
    call init_mesh(resume_mesh, pool_mesh%v0, pool_mesh%v1, pool_mesh%v2)
    resume_mesh%elem_vacuum_sign = 1_i32
    resume_mesh%vacuum_normals = resume_mesh%normals
    call load_restart_checkpoint( &
      trim(pool_checkpoint_dir), resume_mesh, resume_stats, pool_has_restart, mpi_rank=0_i32, mpi_size=mpi%size + 1_i32, &
      app=pool_cfg &
      )
    call assert_true(pool_has_restart, 'shared pool checkpoint should load with a different rank count')
    call assert_equal_i32(resume_stats%batches, 1_i32, 'shared pool resumed batch count mismatch')
    call assert_equal_i64(resume_stats%processed_particles, pool_stats%processed_particles, &
                          'shared pool resumed processed_particles mismatch')
    call mpi_world_barrier(mpi)
    call cleanup_pool_dir(mpi, pool_dir)
  end subroutine run_mpi_shared_pool_test

  subroutine cleanup_pool_dir(mpi, pool_dir)
    type(mpi_context), intent(in) :: mpi
    character(len=*), intent(in) :: pool_dir
    character(len=1024) :: pool_rng_path

    pool_rng_path = restart_rng_state_path(pool_dir//'/checkpoints/slot0', mpi=mpi)
    call delete_file_if_exists(pool_rng_path)
    pool_rng_path = restart_rng_state_path(pool_dir, mpi=mpi)
    call delete_file_if_exists(pool_rng_path)
    call mpi_world_barrier(mpi)
    if (mpi_is_root(mpi)) then
      call cleanup_periodic_slot(pool_dir//'/checkpoints/slot0')
      call delete_file_if_exists(pool_dir//'/checkpoint_latest.txt')
      call delete_file_if_exists(pool_dir//'/checkpoint_latest.txt.tmp')
      call remove_empty_directory(pool_dir//'/checkpoints')
      call remove_empty_directory(pool_dir)
    end if
    call mpi_world_barrier(mpi)
  end subroutine cleanup_pool_dir

  subroutine cleanup_periodic_slot(slot_dir)
    character(len=*), intent(in) :: slot_dir

    call delete_file_if_exists(trim(slot_dir)//'/summary.txt')
    call delete_file_if_exists(trim(slot_dir)//'/charges.csv')
    call delete_file_if_exists(trim(slot_dir)//'/macro_residuals.csv')
    call delete_file_if_exists(trim(slot_dir)//'/charge_ledger.csv')
    call delete_file_if_exists(trim(slot_dir)//'/checkpoint_complete.txt')
    call delete_file_if_exists(trim(slot_dir)//'/checkpoint_complete.txt.tmp')
    call remove_empty_directory(slot_dir)
  end subroutine cleanup_periodic_slot

  subroutine setup_mpi_neutral_return_fixture(mesh, cfg, global_rays_per_batch)
    type(mesh_type), intent(out) :: mesh
    type(app_config), intent(out) :: cfg
    integer(i32), intent(in) :: global_rays_per_batch
    real(dp) :: tri_v0(3, 4), tri_v1(3, 4), tri_v2(3, 4)

    if (global_rays_per_batch < 1_i32) then
      error stop 'neutral-return MPI fixture requires at least one ray'
    end if

    ! 左98%のz=0.75側は同じbatch内で帰還し、右2%のz=0.25側は未解決になる。
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
    call init_mesh(mesh, tri_v0, tri_v1, tri_v2)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals

    call default_app_config(cfg)
    cfg%sim%rng_seed = 998_i32
    cfg%sim%batch_count = 1_i32
    cfg%sim%batch_duration = 1.0_dp
    cfg%sim%dt = 0.6_dp
    cfg%sim%max_step = 1_i32
    cfg%sim%q_floor = 1.0e-30_dp
    cfg%sim%field_solver = 'direct'
    cfg%field%backend = 'direct'
    cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%sim%bc_high(3) = bc_open
    cfg%n_particle_species = 1_i32
    cfg%particle_species(1) = species_from_defaults()
    cfg%particle_species(1)%source_mode = 'photo_raycast'
    cfg%particle_species(1)%rays_per_batch = global_rays_per_batch
    cfg%particle_species(1)%emit_current_density_a_m2 = 1.0_dp
    cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    cfg%particle_species(1)%boundary_high(3) = bc_reflect
    cfg%particle_species(1)%surface_charge_closure = 'neutral_return'
    cfg%particle_species(1)%q_particle = -1.0_dp
    cfg%particle_species(1)%m_particle = 1.0_dp
    cfg%particle_species(1)%temperature_k = 0.0_dp
    cfg%particle_species(1)%normal_drift_speed = 1.0_dp
    cfg%particle_species(1)%inject_face = 'z_high'
    cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    cfg%particle_species(1)%has_ray_direction = .true.
  end subroutine setup_mpi_neutral_return_fixture

  subroutine pack_neutral_return_signature(mesh, ledger, signature)
    type(mesh_type), intent(in) :: mesh
    type(charge_ledger_type), intent(in) :: ledger
    real(dp), intent(out) :: signature(:)

    if (size(mesh%q_elem) /= 4 .or. size(signature) /= 13) then
      error stop 'neutral-return MPI signature storage mismatch'
    end if
    signature = [ &
                mesh%q_elem, ledger%surface_charge_before, ledger%surface_charge_after, &
                ledger%emitted_from_surface(1), ledger%absorbed_on_surface(1), &
                ledger%discarded_unresolved(1), ledger%neutral_return_correction(1), &
                ledger%neutral_return_weight_scale(1), &
                ledger%neutral_return_unresolved_fraction(1), ledger%residual() &
                ]
  end subroutine pack_neutral_return_signature

  subroutine write_summary_fixture(dir_path, mpi_world_size)
    character(len=*), intent(in) :: dir_path
    integer(i32), intent(in) :: mpi_world_size
    integer :: u, ios

    open (newunit=u, file=trim(dir_path)//'/summary.txt', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open MPI summary fixture'
    write (u, '(a)') 'mesh_nelem=1'
    write (u, '(a,i0)') 'mpi_world_size=', mpi_world_size
    write (u, '(a)') 'processed_particles=8'
    write (u, '(a)') 'absorbed=8'
    write (u, '(a)') 'escaped=0'
    write (u, '(a)') 'batches=2'
    write (u, '(a)') 'escaped_boundary=0'
    write (u, '(a)') 'survived_max_step=0'
    write (u, '(a)') 'last_rel_change=1.0e-4'
    close (u)
  end subroutine write_summary_fixture

  subroutine write_charges_fixture(dir_path)
    character(len=*), intent(in) :: dir_path
    integer :: u, ios

    open (newunit=u, file=trim(dir_path)//'/charges.csv', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open MPI charges fixture'
    write (u, '(a)') 'elem_idx,charge_C'
    write (u, '(a)') '1,2.0'
    close (u)
  end subroutine write_charges_fixture

end program test_mpi_hybrid
