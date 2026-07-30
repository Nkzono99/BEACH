!> MPI+OpenMPハイブリッド実行時の集約・rank別resumeファイルを検証する。
program test_mpi_hybrid
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: eps0, pi, qe
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_is_root, mpi_select_lowest_rank_i32_values, &
                     mpi_allreduce_sum_i32_scalar, mpi_allreduce_sum_real_dp_array, &
                     mpi_allreduce_sum_i64_array, &
                     mpi_bcast_i32_array, mpi_bcast_real_dp_array, mpi_gatherv_real_dp_array, &
                     mpi_world_barrier
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            init_particle_batch_from_config
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file, &
                         restart_rng_state_path, restart_macro_residual_path
  use bem_types, only: mesh_type, particles_soa, sim_stats, injection_state, bc_open, bc_periodic
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_diagnostics_type
  use bem_dynamic_k0_mean, only: dynamic_k0_ok, dynamic_k0_step_type
  use bem_dynamic_k0_zhao, only: &
    measured_interface_energy_distribution_type, build_measured_interface_energy_distribution, &
    advance_dynamic_k0_zhao, zhao_profile_barrier_energy, measured_sample_escape_fraction
  use bem_outer_plasma_kinetic, only: &
    kinetic_outer_plasma_options_type, solve_outer_plasma_zhao_stationary
  use bem_outer_plasma_kinetic_runtime, only: resolve_kinetic_outer_options
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_state_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, &
                          assert_close_dp, delete_file_if_exists, &
                          ensure_directory, remove_empty_directory
  implicit none

  type(mpi_context) :: mpi
  type(mesh_type) :: mesh, mesh_restart
  type(app_config) :: cfg, reservoir_cfg
  type(sim_stats) :: stats, stats_restart
  type(injection_state) :: state, state_restart, reservoir_state
  type(particles_soa) :: reservoir_particles
  type(charge_ledger_type) :: ledger
  logical :: has_restart
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  integer :: u, ios
  character(len=256) :: line
  integer(i32) :: n_lines, selected_rank, expected_rank, batch_idx, global_reservoir_count
  integer(i32) :: local_failure_values(4), selected_failure_values(4)
  character(len=*), parameter :: history_path = 'test_mpi_hybrid_history_tmp.csv'
  character(len=*), parameter :: out_dir = 'test_mpi_hybrid_restart_tmp'
  character(len=1024) :: rng_path, residual_path

  call mpi_initialize(mpi)

  if (mpi_is_root(mpi)) then
    call delete_file_if_exists(history_path)
    call delete_file_if_exists(out_dir//'/summary.txt')
    call delete_file_if_exists(out_dir//'/charges.csv')
    call delete_file_if_exists(out_dir//'/rng_state.txt')
    call delete_file_if_exists(out_dir//'/macro_residuals.csv')
  end if
  call mpi_world_barrier(mpi)

  rng_path = restart_rng_state_path(out_dir, mpi=mpi)
  residual_path = restart_macro_residual_path(out_dir, mpi=mpi)
  call delete_file_if_exists(rng_path)
  if (mpi_is_root(mpi)) call delete_file_if_exists(residual_path)
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

  call test_init(9)

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
  call test_end()

  call test_begin('mpi_simulation')
  if (mpi_is_root(mpi)) then
    open (newunit=u, file=history_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open MPI hybrid history fixture'
    call run_absorption_insulator( &
      mesh, cfg, stats, history_unit=u, history_stride=1_i32, mpi=mpi, charge_ledger=ledger &
      )
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
  call test_end()

  call test_begin('mpi_neutral_return_layout_invariance')
  call run_mpi_neutral_return_layout_test(mpi)
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
    end do
    close (u)
    call assert_equal_i32(n_lines, 1_i32, 'mpi history snapshot line count mismatch')
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
  reservoir_cfg%particle_species(1)%source_mode = 'reservoir_face'
  reservoir_cfg%particle_species(1)%number_density_m3 = 1.0_dp
  reservoir_cfg%particle_species(1)%has_number_density_m3 = .true.
  reservoir_cfg%particle_species(1)%temperature_k = 0.0_dp
  reservoir_cfg%particle_species(1)%has_temperature_k = .true.
  reservoir_cfg%particle_species(1)%q_particle = 0.0_dp
  reservoir_cfg%particle_species(1)%m_particle = 1.0_dp
  reservoir_cfg%particle_species(1)%w_particle = 4.0_dp
  reservoir_cfg%particle_species(1)%has_w_particle = .true.
  reservoir_cfg%particle_species(1)%inject_face = 'z_low'
  reservoir_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 0.0_dp]
  reservoir_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 0.0_dp]
  reservoir_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, 1.0_dp]
  allocate (reservoir_state%macro_residual(1))
  reservoir_state%macro_residual = 0.0_dp
  do batch_idx = 1_i32, 4_i32
    call init_particle_batch_from_config( &
      reservoir_cfg, batch_idx, reservoir_particles, state=reservoir_state, mpi=mpi &
      )
    global_reservoir_count = reservoir_particles%n
    call mpi_allreduce_sum_i32_scalar(mpi, global_reservoir_count)
    call assert_equal_i32( &
      global_reservoir_count, merge(1_i32, 0_i32, batch_idx == 4_i32), &
      'MPI global reservoir count sequence mismatch' &
      )
    call assert_close_dp( &
      reservoir_state%macro_residual(1), modulo(0.25_dp*real(batch_idx, dp), 1.0_dp), 1.0e-15_dp, &
      'MPI global reservoir residual mismatch' &
      )
  end do
  call test_end()

  call test_begin('mpi_kinetic_outer_root_broadcast_and_hold')
  call run_mpi_kinetic_outer_test(mpi)
  call test_end()

  call test_begin('mpi_implicit_mean_layout_invariance')
  call run_mpi_implicit_mean_layout_test(mpi)
  call test_end()

  call test_begin('mpi_dynamic_zhao_empirical_layout_invariance')
  call run_mpi_dynamic_zhao_layout_test(mpi)
  call test_end()

  call test_begin('mpi_restart')
  call ensure_directory(out_dir)
  if (mpi_is_root(mpi)) then
    call write_summary_fixture(out_dir, mpi%size)
    call write_charges_fixture(out_dir)
  end if
  call mpi_world_barrier(mpi)

  allocate (state%macro_residual(1))
  state%macro_residual(1) = 0.25d0 + 0.01d0*real(mpi%rank, dp)
  call write_rng_state_file(out_dir, mpi=mpi)
  call write_macro_residuals_file(out_dir, state, mpi=mpi)
  call mpi_world_barrier(mpi)

  call init_mesh(mesh_restart, v0, v1, v2)
  mesh_restart%elem_vacuum_sign = 1_i32
  mesh_restart%vacuum_normals = mesh_restart%normals
  allocate (state_restart%macro_residual(1))
  state_restart%macro_residual = 0.0d0
  call load_restart_checkpoint(out_dir, mesh_restart, stats_restart, has_restart, state_restart, mpi=mpi)
  call assert_true(has_restart, 'mpi restart should be detected')
  call assert_equal_i64(stats_restart%processed_particles, 8_i64, 'mpi restart processed_particles mismatch')
  call assert_equal_i32(stats_restart%batches, 2_i32, 'mpi restart batches mismatch')
  call assert_close_dp(mesh_restart%q_elem(1), 2.0d0, 1.0d-12, 'mpi restart charge mismatch')
  call assert_close_dp(state_restart%macro_residual(1), 0.25d0, 1.0d-12, 'MPI residual must be restored from root')
  call test_end()

  call delete_file_if_exists(rng_path)
  if (mpi_is_root(mpi)) call delete_file_if_exists(residual_path)
  call mpi_world_barrier(mpi)
  if (mpi_is_root(mpi)) then
    call delete_file_if_exists(out_dir//'/summary.txt')
    call delete_file_if_exists(out_dir//'/charges.csv')
    call delete_file_if_exists(out_dir//'/rng_state.txt')
    call delete_file_if_exists(out_dir//'/macro_residuals.csv')
    call remove_empty_directory(out_dir)
  end if

  call test_summary()

  call mpi_shutdown(mpi)

contains

  !> 同一ray集合を複製したMPI実行とserial参照でneutral-return閉包を比較する。
  !!
  !! MPI側は各rankがserial参照と同じray集合を追跡し、global raysを
  !! world-size倍にしてmacro chargeを1/world-sizeへ落とす。したがって
  !! 全rankの物理的な放出・帰還・未解決電荷とmesh depositはserial参照と一致する。
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
      call setup_mpi_neutral_return_fixture(reference_mesh, reference_cfg, rays_per_rank)
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
    ! Test-only common stream: every rank traces the same rays while each
    ! macro charge is divided by the global ray count.
    call seed_particles_from_config(distributed_cfg, mpi_rank=0_i32, mpi_size=1_i32)
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
      distributed_ledger%emitted_count(1), int(mpi%size, i64)*reference_counts(1), &
      'neutral-return MPI emitted count was not reduced across ranks' &
      )
    call assert_equal_i64( &
      distributed_ledger%absorbed_count(1), int(mpi%size, i64)*reference_counts(2), &
      'neutral-return MPI absorbed count was not reduced across ranks' &
      )
    call assert_equal_i64( &
      distributed_ledger%discarded_unresolved_count(1), int(mpi%size, i64)*reference_counts(3), &
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
    cfg%particle_species(1)%z_high_boundary = 'reflect'
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

  !> rank配置の異なる同一(E_n, |dq|)集合をglobal CDFからZhao rootへ通す。
  !!
  !! layout 1ではrank 0をzero-rayにし、layout 2では同じ4 rayを再配置する。
  !! 3 rayは同一energyかつ異なるweightを持つため、可変長gather、group化、
  !! nonlinear Q(Phi_I)、collective profile採用、charge ledgerを一度に検証する。
  subroutine run_mpi_dynamic_zhao_layout_test(mpi)
    type(mpi_context), intent(in) :: mpi
    integer(i32), parameter :: signature_size = 15_i32
    type(mesh_type) :: zhao_mesh
    type(app_config) :: zhao_cfg
    type(electrostatic_snapshot_type) :: zhao_snapshot
    type(kinetic_outer_plasma_options_type) :: zhao_options
    type(outer_plasma_state_type) :: stationary_outer, accepted_outer
    type(measured_interface_energy_distribution_type) :: distribution
    type(dynamic_k0_step_type) :: step
    type(charge_ledger_type) :: layout_ledger
    real(dp), allocatable :: local_energy_j(:), local_charge_c(:)
    real(dp), allocatable :: global_energy_j(:), global_charge_c(:)
    real(dp) :: signatures(signature_size, 2), global_signature(signature_size)
    real(dp) :: candidate_charge(2), sample_weight(4)
    real(dp) :: area_xy, time_step_s, barrier_j, reference_current_density
    real(dp) :: source_charge, anchor_charge, base_charge, high_fraction
    real(dp) :: effective_source_scale, grouped_low_charge, tolerance
    real(dp) :: real_metadata(1)
    integer(i32) :: layout_index, value_index, status, distribution_status
    integer(i32) :: group_count, root_zero_marker, branch_code_sum
    integer(i32) :: metadata(2)
    character(len=1) :: layout_branch(2)
    character(len=160) :: assertion_message
    character(len=256) :: message

    call configure_mpi_dynamic_zhao_options(zhao_options)
    call solve_outer_plasma_zhao_stationary(zhao_options, stationary_outer, status, message)
    call assert_equal_i32(status, outer_plasma_ok, &
                          'MPI dynamic Zhao stationary solve failed: '//trim(message))
    call assert_true(stationary_outer%ready .and. stationary_outer%zhao_branch == 'A', &
                     'MPI dynamic Zhao fixture did not resolve a virtual-cathode A branch')
    if (status /= outer_plasma_ok) return

    area_xy = (9.89949e-5_dp)**2
    time_step_s = 1.0e-8_dp
    barrier_j = zhao_profile_barrier_energy(stationary_outer, zhao_options%photoelectron_charge)
    high_fraction = exp(-barrier_j/zhao_options%photoelectron_temperature_j)
    reference_current_density = abs(zhao_options%photoelectron_charge)* &
                                zhao_options%photoelectron_reference_density* &
                                sin(zhao_options%zhao_alpha_deg*pi/180.0_dp)* &
                                sqrt(2.0_dp*zhao_options%photoelectron_temperature_j/ &
                                     zhao_options%photoelectron_mass)/(2.0_dp*sqrt(pi))
    source_charge = reference_current_density*area_xy*time_step_s
    anchor_charge = eps0*area_xy*stationary_outer%interface_field
    sample_weight(1) = high_fraction*source_charge
    sample_weight(2) = 0.2_dp*(source_charge - sample_weight(1))
    sample_weight(3) = 0.3_dp*(source_charge - sample_weight(1))
    sample_weight(4) = source_charge - sum(sample_weight(:3))
    base_charge = anchor_charge - sample_weight(1)
    zhao_options%interface_field = stationary_outer%interface_field
    call setup_mpi_dynamic_zhao_snapshot( &
      stationary_outer, zhao_options, zhao_mesh, zhao_cfg, zhao_snapshot, mpi &
      )

    signatures = 0.0_dp
    layout_branch = ' '
    do layout_index = 1_i32, 2_i32
      call assign_mpi_dynamic_zhao_layout( &
        layout_index, mpi, barrier_j, sample_weight, local_energy_j, local_charge_c &
        )

      root_zero_marker = merge(1_i32, 0_i32, mpi%rank == 0_i32 .and. size(local_energy_j) == 0)
      call mpi_allreduce_sum_i32_scalar(mpi, root_zero_marker)
      if (layout_index == 1_i32) then
        call assert_equal_i32( &
          root_zero_marker, merge(1_i32, 0_i32, mpi%size > 1_i32), &
          'MPI dynamic Zhao zero-ray root fixture was not constructed as requested' &
          )
      end if

      call mpi_gatherv_real_dp_array(mpi, local_energy_j, global_energy_j, 0_i32)
      call mpi_gatherv_real_dp_array(mpi, local_charge_c, global_charge_c, 0_i32)
      distribution_status = -1_i32
      group_count = 0_i32
      grouped_low_charge = 0.0_dp
      step = dynamic_k0_step_type()
      accepted_outer = outer_plasma_state_type()
      effective_source_scale = 0.0_dp
      message = ''
      if (mpi_is_root(mpi)) then
        call build_measured_interface_energy_distribution( &
          global_energy_j, global_charge_c, distribution, distribution_status, message &
          )
        if (distribution_status == dynamic_k0_ok) then
          group_count = distribution%group_count
          if (group_count >= 2_i32) grouped_low_charge = distribution%group_charge_c(2)
          call advance_dynamic_k0_zhao( &
            zhao_options, stationary_outer, 'e_bottom_zero', area_xy, anchor_charge, &
            base_charge, time_step_s, distribution, step, accepted_outer, &
            effective_source_scale, message &
            )
        else
          step%status = distribution_status
        end if
      end if
      metadata = [distribution_status, group_count]
      call mpi_bcast_i32_array(mpi, metadata, 0_i32)
      distribution_status = metadata(1)
      group_count = metadata(2)
      real_metadata = [grouped_low_charge]
      call mpi_bcast_real_dp_array(mpi, real_metadata, 0_i32)
      grouped_low_charge = real_metadata(1)
      call broadcast_test_dynamic_zhao_step(mpi, step, effective_source_scale)

      call assert_equal_i32(distribution_status, dynamic_k0_ok, &
                            'MPI dynamic Zhao global energy distribution build failed')
      call assert_equal_i32(group_count, 2_i32, &
                            'MPI dynamic Zhao equal-energy rays were not grouped globally')
      call assert_close_dp( &
        grouped_low_charge, sum(sample_weight(2:)), 1.0e-12_dp*source_charge, &
        'MPI dynamic Zhao unequal weights were not summed in the common energy group' &
        )
      call assert_equal_i32(step%status, dynamic_k0_ok, &
                            'MPI dynamic Zhao nonlinear update failed: '//trim(message))
      if (distribution_status /= dynamic_k0_ok .or. step%status /= dynamic_k0_ok) return

      zhao_snapshot%kinetic_options%photoelectron_source_scale = effective_source_scale
      zhao_snapshot%kinetic_options%photoelectron_population_fraction = 1.0_dp
      zhao_snapshot%kinetic_options%photoelectron_column_closure_enabled = .false.
      zhao_snapshot%kinetic_options%zhao_branch = 'a'
      candidate_charge = 0.5_dp*step%surface_charge_after_c
      call zhao_snapshot%adopt_mean_outer(zhao_mesh, candidate_charge, accepted_outer)

      call build_mpi_dynamic_zhao_ledger( &
        mpi, anchor_charge, base_charge, step, local_energy_j, local_charge_c, layout_ledger &
        )
      call assert_close_dp( &
        step%photoelectron_source_charge_c, -layout_ledger%interface_outward_gross(2), &
        4096.0_dp*epsilon(1.0_dp)*source_charge, &
        'MPI dynamic Zhao gathered source disagrees with the independently reduced interface charge' &
        )
      call assert_close_dp(effective_source_scale, 1.0_dp, 1.0e-13_dp, &
                           'MPI dynamic Zhao measured source normalization changed')
      call assert_close_dp( &
        layout_ledger%interface_outward_gross(2), &
        layout_ledger%interface_returned_gross(2) + layout_ledger%escaped_to_infinity(2), &
        1.0e-12_dp*source_charge, &
        'MPI dynamic Zhao gross interface ledger does not close' &
        )
      call assert_close_dp(layout_ledger%residual(), 0.0_dp, 1.0e-12_dp*source_charge, &
                           'MPI dynamic Zhao charge ledger residual does not close')
      call assert_close_dp(step%backward_euler_residual_charge_c, 0.0_dp, &
                           1.0e-12_dp*source_charge, &
                           'MPI dynamic Zhao nonlinear charge residual does not close')

      call pack_mpi_dynamic_zhao_signature( &
        zhao_snapshot, step, layout_ledger, signatures(:, layout_index) &
        )
      layout_branch(layout_index) = zhao_snapshot%outer%zhao_branch
      branch_code_sum = int(iachar(layout_branch(layout_index)), i32)
      call mpi_allreduce_sum_i32_scalar(mpi, branch_code_sum)
      call assert_equal_i32( &
        branch_code_sum, mpi%size*int(iachar('A'), i32), &
        'MPI dynamic Zhao collective adoption changed branch across ranks' &
        )

      global_signature = signatures(:, layout_index)
      call mpi_allreduce_sum_real_dp_array(mpi, global_signature)
      global_signature = global_signature/real(mpi%size, dp)
      do value_index = 1_i32, signature_size
        tolerance = mpi_dynamic_zhao_signature_tolerance( &
                    value_index, global_signature(value_index), &
                    signatures(value_index, layout_index), source_charge, &
                    zhao_options%photoelectron_temperature_j &
                    )
        write (assertion_message, '(a,i0,a,i0)') &
          'dynamic Zhao rank mismatch at layout ', layout_index, ', component ', value_index
        call assert_close_dp( &
          signatures(value_index, layout_index), global_signature(value_index), tolerance, &
          trim(assertion_message) &
          )
      end do
    end do

    call assert_true(all(layout_branch == 'A'), &
                     'MPI dynamic Zhao layout comparison left the committed branch')
    do value_index = 1_i32, signature_size
      tolerance = mpi_dynamic_zhao_signature_tolerance( &
                  value_index, signatures(value_index, 1), signatures(value_index, 2), &
                  source_charge, zhao_options%photoelectron_temperature_j &
                  )
      write (assertion_message, '(a,i0)') &
        'dynamic Zhao rank-layout signature mismatch at component ', value_index
      call assert_close_dp( &
        signatures(value_index, 2), signatures(value_index, 1), tolerance, &
        trim(assertion_message) &
        )
    end do
  end subroutine run_mpi_dynamic_zhao_layout_test

  subroutine configure_mpi_dynamic_zhao_options(options)
    type(kinetic_outer_plasma_options_type), intent(out) :: options
    real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
    real(dp), parameter :: proton_mass = 1.67262192369e-27_dp

    options = kinetic_outer_plasma_options_type()
    options%kinetic_closure = 'zhao_charge_driven'
    options%zhao_branch = 'a'
    options%grid_points = 257_i32
    options%domain_length = 2.0_dp
    options%electron_charge = -qe
    options%electron_mass = electron_mass
    options%electron_temperature_j = 12.0_dp*qe
    options%electron_drift_infinity = 4.68e5_dp*sin(60.0_dp*pi/180.0_dp)
    options%ion_charge = qe
    options%ion_mass = proton_mass
    options%ion_density_infinity = 8.7e6_dp
    options%ion_temperature_j = 0.1_dp*qe
    options%ion_drift_infinity = options%electron_drift_infinity
    options%photoelectron_charge = -qe
    options%photoelectron_mass = electron_mass
    options%photoelectron_temperature_j = 2.2_dp*qe
    options%photoelectron_reference_density = 64.0e6_dp
    options%photoelectron_population_fraction = 1.0_dp
    options%photoelectron_source_scale = 1.0_dp
    options%photoelectron_column_closure_enabled = .false.
    options%zhao_alpha_deg = 60.0_dp
  end subroutine configure_mpi_dynamic_zhao_options

  subroutine setup_mpi_dynamic_zhao_snapshot( &
    stationary_outer, options, mesh, cfg, snapshot, mpi &
    )
    type(outer_plasma_state_type), intent(in) :: stationary_outer
    type(kinetic_outer_plasma_options_type), intent(in) :: options
    type(mesh_type), intent(out) :: mesh
    type(app_config), intent(out) :: cfg
    type(electrostatic_snapshot_type), intent(out) :: snapshot
    type(mpi_context), intent(in) :: mpi
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)
    real(dp), parameter :: box_width = 9.89949e-5_dp
    real(dp), parameter :: support_z = 2.0e-6_dp
    real(dp), parameter :: interface_z = 2.0e-4_dp
    real(dp) :: initial_charge

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, support_z]
    tri_v1(:, 1) = [box_width, 0.0_dp, support_z]
    tri_v2(:, 1) = [box_width, box_width, support_z]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, support_z]
    tri_v1(:, 2) = [box_width, box_width, support_z]
    tri_v2(:, 2) = [0.0_dp, box_width, support_z]
    initial_charge = eps0*box_width**2*stationary_outer%interface_field
    call init_mesh(mesh, tri_v0, tri_v1, tri_v2, &
                   q0=[0.5_dp*initial_charge, 0.5_dp*initial_charge])
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals

    call default_app_config(cfg)
    cfg%sim%field_solver = 'direct'
    cfg%sim%field_bc_mode = 'periodic2'
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    cfg%sim%box_max = [box_width, box_width, interface_z]
    cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    cfg%field%backend = 'direct'
    cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    cfg%panel%surface_side_policy = 'per_element'
    cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    cfg%periodic2%zero_mode_policy = 'exclude_k0'
    cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    cfg%periodic2%reference_mode_layers = 2_i32
    cfg%periodic2%panel_quadrature_order = 4_i32
    cfg%outer_plasma%model = 'kinetic_1d'
    cfg%outer_plasma%kinetic_closure = 'zhao_charge_driven'
    cfg%outer_plasma%interface_z = interface_z
    cfg%outer_plasma%debye_length = stationary_outer%debye_length
    cfg%outer_plasma%thermal_voltage = 12.0_dp
    call snapshot%init( &
      mesh, cfg%sim, cfg%field, cfg%periodic2, cfg%panel, cfg%outer_plasma, &
      kinetic_options=options, mpi=mpi &
      )
  end subroutine setup_mpi_dynamic_zhao_snapshot

  subroutine assign_mpi_dynamic_zhao_layout( &
    layout_index, mpi, barrier_j, sample_weight, local_energy_j, local_charge_c &
    )
    integer(i32), intent(in) :: layout_index
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(in) :: barrier_j, sample_weight(4)
    real(dp), allocatable, intent(out) :: local_energy_j(:), local_charge_c(:)
    real(dp) :: sample_energy(4)

    sample_energy = [ &
                    barrier_j + 0.5_dp*qe, barrier_j - 0.5_dp*qe, &
                    barrier_j - 0.5_dp*qe, barrier_j - 0.5_dp*qe &
                    ]
    select case (layout_index)
    case (1_i32)
      if (mpi%size == 1_i32 .or. mpi%rank == 1_i32) then
        allocate (local_energy_j(4), local_charge_c(4))
        local_energy_j = sample_energy
        local_charge_c = sample_weight
      else
        allocate (local_energy_j(0), local_charge_c(0))
      end if
    case (2_i32)
      if (mpi%size == 1_i32) then
        allocate (local_energy_j(4), local_charge_c(4))
        local_energy_j = sample_energy([4, 1, 3, 2])
        local_charge_c = sample_weight([4, 1, 3, 2])
      else if (mpi%rank == 0_i32) then
        allocate (local_energy_j(2), local_charge_c(2))
        local_energy_j = sample_energy([4, 1])
        local_charge_c = sample_weight([4, 1])
      else if (mpi%rank == 1_i32) then
        allocate (local_energy_j(2), local_charge_c(2))
        local_energy_j = sample_energy([3, 2])
        local_charge_c = sample_weight([3, 2])
      else
        allocate (local_energy_j(0), local_charge_c(0))
      end if
    case default
      error stop 'unknown MPI dynamic Zhao sample layout'
    end select
  end subroutine assign_mpi_dynamic_zhao_layout

  subroutine broadcast_test_dynamic_zhao_step(mpi, step, effective_source_scale)
    type(mpi_context), intent(in) :: mpi
    type(dynamic_k0_step_type), intent(inout) :: step
    real(dp), intent(inout) :: effective_source_scale
    integer(i32) :: integers(3)
    real(dp) :: values(13)

    integers = 0_i32
    values = 0.0_dp
    if (mpi_is_root(mpi)) then
      integers = [step%status, step%iterations, int(iachar(step%zhao_branch), i32)]
      values = [ &
               step%surface_charge_before_c, step%surface_charge_after_c, &
               step%interface_potential_before_v, step%interface_potential_after_v, &
               step%interface_field_after_v_m, step%photoelectron_escape_fraction, &
               step%photoelectron_return_fraction, step%backward_euler_residual_charge_c, &
               step%photoelectron_barrier_energy_j, step%marginal_photoelectron_energy_j, &
               step%marginal_photoelectron_escape_fraction, step%photoelectron_source_charge_c, &
               effective_source_scale &
               ]
    end if
    call mpi_bcast_i32_array(mpi, integers, 0_i32)
    call mpi_bcast_real_dp_array(mpi, values, 0_i32)
    if (.not. mpi_is_root(mpi)) then
      step = dynamic_k0_step_type()
      step%status = integers(1)
      step%iterations = integers(2)
      step%zhao_branch = achar(integers(3))
      step%surface_charge_before_c = values(1)
      step%surface_charge_after_c = values(2)
      step%interface_potential_before_v = values(3)
      step%interface_potential_after_v = values(4)
      step%interface_field_after_v_m = values(5)
      step%photoelectron_escape_fraction = values(6)
      step%photoelectron_return_fraction = values(7)
      step%backward_euler_residual_charge_c = values(8)
      step%photoelectron_barrier_energy_j = values(9)
      step%marginal_photoelectron_energy_j = values(10)
      step%marginal_photoelectron_escape_fraction = values(11)
      step%photoelectron_source_charge_c = values(12)
      step%zhao_effective_source_scale = values(13)
    end if
    effective_source_scale = values(13)
    step%zhao_effective_source_scale = effective_source_scale
  end subroutine broadcast_test_dynamic_zhao_step

  subroutine build_mpi_dynamic_zhao_ledger( &
    mpi, surface_charge_before, surface_charge_base, step, local_energy_j, local_charge_c, ledger &
    )
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(in) :: surface_charge_before, surface_charge_base
    type(dynamic_k0_step_type), intent(in) :: step
    real(dp), intent(in) :: local_energy_j(:), local_charge_c(:)
    type(charge_ledger_type), intent(out) :: ledger
    real(dp) :: local_source_charge, local_escape_charge, local_return_charge
    real(dp) :: ledger_values(6)
    integer :: sample_index

    local_source_charge = sum(local_charge_c)
    local_escape_charge = 0.0_dp
    do sample_index = 1, size(local_energy_j)
      local_escape_charge = local_escape_charge + local_charge_c(sample_index)* &
                            measured_sample_escape_fraction(local_energy_j(sample_index), step)
    end do
    local_return_charge = local_source_charge - local_escape_charge
    ledger_values = 0.0_dp
    if (mpi_is_root(mpi)) ledger_values(1) = surface_charge_base - surface_charge_before
    ledger_values(2) = -local_source_charge
    ledger_values(3) = -local_return_charge
    ledger_values(4) = -local_escape_charge
    ledger_values(5) = -local_source_charge
    ledger_values(6) = -local_return_charge
    call mpi_allreduce_sum_real_dp_array(mpi, ledger_values)

    call ledger%init(2_i32)
    ledger%surface_charge_before = surface_charge_before
    ledger%surface_charge_after = step%surface_charge_after_c
    ledger%injected_from_remote(1) = ledger_values(1)
    ledger%emitted_from_surface(2) = ledger_values(2)
    ledger%absorbed_on_surface(2) = ledger_values(3)
    ledger%escaped_to_infinity(2) = ledger_values(4)
    ledger%interface_outward_gross(2) = ledger_values(5)
    ledger%interface_returned_gross(2) = ledger_values(6)
  end subroutine build_mpi_dynamic_zhao_ledger

  subroutine pack_mpi_dynamic_zhao_signature(snapshot, step, ledger, signature)
    type(electrostatic_snapshot_type), intent(in) :: snapshot
    type(dynamic_k0_step_type), intent(in) :: step
    type(charge_ledger_type), intent(in) :: ledger
    real(dp), intent(out) :: signature(:)

    if (size(signature) /= 15) error stop 'MPI dynamic Zhao signature storage mismatch'
    signature = [ &
                snapshot%outer%interface_potential, snapshot%outer%interface_field, &
                zhao_profile_barrier_energy( &
                snapshot%outer, snapshot%kinetic_options%photoelectron_charge &
                ), &
                snapshot%outer%photoelectron_source_scale, step%surface_charge_after_c, &
                step%photoelectron_escape_fraction, step%photoelectron_source_charge_c, &
                ledger%injected_from_remote(1), ledger%emitted_from_surface(2), &
                ledger%interface_outward_gross(2), ledger%interface_returned_gross(2), &
                ledger%escaped_to_infinity(2), ledger%residual(), snapshot%outer%zhao_phi_minimum, &
                snapshot%outer%zhao_electron_density_infinity &
                ]
  end subroutine pack_mpi_dynamic_zhao_signature

  pure real(dp) function mpi_dynamic_zhao_signature_tolerance( &
    component, first_value, second_value, source_charge, photoelectron_temperature_j &
    ) result(tolerance)
    integer(i32), intent(in) :: component
    real(dp), intent(in) :: first_value, second_value
    real(dp), intent(in) :: source_charge, photoelectron_temperature_j

    select case (component)
    case (3_i32)
      tolerance = 1.0e-10_dp*max( &
                  abs(first_value), abs(second_value), photoelectron_temperature_j, tiny(1.0_dp) &
                  )
    case (5_i32, 7_i32:13_i32)
      tolerance = 1.0e-10_dp*max( &
                  abs(first_value), abs(second_value), source_charge, tiny(1.0_dp) &
                  )
    case default
      tolerance = 1.0e-10_dp*max(abs(first_value), abs(second_value), tiny(1.0_dp))
    end select
  end function mpi_dynamic_zhao_signature_tolerance

  !> 1-ray直列参照と、各rankに同じ1-rayを複製したMPI実行を比較する。
  !!
  !! MPI側のglobal raysはworld sizeなので、各rayのmacro chargeは直列参照の
  !! 1/world-sizeになる。同じ乱数seedを全rankへ与えることで粒子軌道を揃え、
  !! stochastic sampling差を混ぜずにanalytic weight、return shadow moments、
  !! charge ledger reductionのrank配置不変性を検証する。
  subroutine run_mpi_implicit_mean_layout_test(mpi)
    type(mpi_context), intent(in) :: mpi
    integer(i32), parameter :: signature_size = 14_i32
    integer(i32), parameter :: photoelectron_index = 3_i32
    type(mesh_type) :: reference_mesh, distributed_mesh
    type(app_config) :: reference_cfg, distributed_cfg
    type(sim_stats) :: reference_stats, distributed_stats
    type(injection_state) :: reference_injection, distributed_injection
    type(charge_ledger_type) :: reference_ledger, distributed_ledger
    type(electrostatic_diagnostics_type) :: reference_diagnostics, distributed_diagnostics
    real(dp) :: reference_signature(signature_size), distributed_signature(signature_size)
    real(dp) :: tolerance, charge_scale
    integer(i32) :: value_index
    character(len=128) :: assertion_message

    reference_signature = 0.0_dp
    if (mpi_is_root(mpi)) then
      call setup_mpi_implicit_mean_fixture(reference_mesh, reference_cfg, 1_i32)
      allocate (reference_injection%macro_residual(reference_cfg%n_particle_species))
      reference_injection%macro_residual = 0.0_dp
      call seed_particles_from_config(reference_cfg)
      call run_absorption_insulator( &
        reference_mesh, reference_cfg, reference_stats, inject_state=reference_injection, &
        charge_ledger=reference_ledger, electrostatic_diagnostics=reference_diagnostics &
        )
      call pack_implicit_mean_signature( &
        reference_mesh, reference_ledger, reference_diagnostics, reference_signature &
        )
    end if
    call mpi_allreduce_sum_real_dp_array(mpi, reference_signature)

    call setup_mpi_implicit_mean_fixture(distributed_mesh, distributed_cfg, mpi%size)
    allocate (distributed_injection%macro_residual(distributed_cfg%n_particle_species))
    distributed_injection%macro_residual = 0.0_dp
    ! Test-only common stream: each rank traces the same ray while its macro
    ! charge is divided by the global rays_per_batch.
    call seed_particles_from_config(distributed_cfg, mpi_rank=0_i32, mpi_size=1_i32)
    call run_absorption_insulator( &
      distributed_mesh, distributed_cfg, distributed_stats, inject_state=distributed_injection, mpi=mpi, &
      charge_ledger=distributed_ledger, electrostatic_diagnostics=distributed_diagnostics &
      )
    call pack_implicit_mean_signature( &
      distributed_mesh, distributed_ledger, distributed_diagnostics, distributed_signature &
      )

    do value_index = 1_i32, signature_size
      tolerance = 1.0e-10_dp*max( &
                  abs(reference_signature(value_index)), abs(distributed_signature(value_index)), tiny(1.0_dp) &
                  )
      write (assertion_message, '(a,i0)') 'implicit-mean MPI signature mismatch at component ', value_index
      call assert_close_dp( &
        distributed_signature(value_index), reference_signature(value_index), tolerance, trim(assertion_message) &
        )
    end do

    call assert_equal_i64( &
      distributed_ledger%emitted_count(photoelectron_index), int(mpi%size, i64), &
      'implicit-mean MPI emitted count was not reduced across ranks' &
      )
    call assert_equal_i64( &
      distributed_ledger%absorbed_count(photoelectron_index), int(mpi%size, i64), &
      'implicit-mean MPI returned count was not reduced across ranks' &
      )
    call assert_equal_i64( &
      distributed_ledger%escaped_count(photoelectron_index), 0_i64, &
      'analytic escape charge must not invent a sampled escape event count' &
      )
    call assert_true(distributed_ledger%emitted_from_surface(photoelectron_index) < 0.0_dp, &
                     'implicit-mean MPI fixture did not emit photoelectron charge')
    call assert_true(distributed_ledger%escaped_to_infinity(photoelectron_index) < 0.0_dp, &
                     'implicit-mean analytic weighting did not retain its nonzero escape charge')
    call assert_true( &
      abs(distributed_ledger%escaped_to_infinity(photoelectron_index)) < &
      abs(distributed_ledger%emitted_from_surface(photoelectron_index)), &
      'implicit-mean analytic escape charge must be smaller than emitted charge' &
      )
    call assert_close_dp( &
      distributed_ledger%interface_outward_gross(photoelectron_index), &
      distributed_ledger%interface_returned_gross(photoelectron_index) + &
      distributed_ledger%escaped_to_infinity(photoelectron_index), &
      1.0e-12_dp*abs(distributed_ledger%emitted_from_surface(photoelectron_index)), &
      'implicit-mean MPI interface gross charge does not close' &
      )
    charge_scale = max( &
                   abs(distributed_ledger%surface_charge_before), &
                   abs(distributed_ledger%surface_charge_after), tiny(1.0_dp) &
                   )
    call assert_true( &
      abs(distributed_ledger%residual()) <= 1.0e-12_dp*charge_scale, &
      'implicit-mean MPI charge ledger residual is not closed' &
      )
    call assert_true(distributed_diagnostics%implicit_mean_shadow_diagnostics_available, &
                     'implicit-mean MPI shadow diagnostics are unavailable')
    call assert_true( &
      distributed_diagnostics%implicit_mean_last_returned_outer_flight_time_mean > 0.0_dp, &
      'implicit-mean MPI returned shadow flight time is not positive' &
      )
    call assert_true( &
      distributed_diagnostics%implicit_mean_last_returning_pe_column_charge_per_area > 0.0_dp, &
      'implicit-mean MPI returning shadow column is not positive' &
      )
  end subroutine run_mpi_implicit_mean_layout_test

  subroutine setup_mpi_implicit_mean_fixture(mesh, cfg, global_rays_per_batch)
    type(mesh_type), intent(out) :: mesh
    type(app_config), intent(out) :: cfg
    integer(i32), intent(in) :: global_rays_per_batch
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)
    real(dp), parameter :: electron_mass = 9.1093837139e-31_dp
    real(dp), parameter :: proton_mass = 1.67262192595e-27_dp
    real(dp), parameter :: box_width = 1.0e-4_dp
    real(dp), parameter :: interface_z = 2.0e-4_dp
    real(dp), parameter :: support_z = 2.0e-6_dp
    real(dp), parameter :: debye_length = 0.1_dp
    real(dp), parameter :: interface_potential = 35.0_dp
    real(dp) :: initial_charge

    if (global_rays_per_batch < 1_i32) error stop 'implicit-mean MPI fixture requires at least one ray'

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, support_z]
    tri_v1(:, 1) = [box_width, 0.0_dp, support_z]
    tri_v2(:, 1) = [box_width, box_width, support_z]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, support_z]
    tri_v1(:, 2) = [box_width, box_width, support_z]
    tri_v2(:, 2) = [0.0_dp, box_width, support_z]
    initial_charge = eps0*box_width**2*interface_potential/debye_length
    call init_mesh(mesh, tri_v0, tri_v1, tri_v2, q0=[0.5_dp*initial_charge, 0.5_dp*initial_charge])
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals

    call default_app_config(cfg)
    cfg%mesh_mode = 'obj'
    cfg%sim%rng_seed = 4321_i32
    cfg%sim%batch_count = 1_i32
    cfg%sim%dt = 3.0e-10_dp
    cfg%sim%batch_duration = 3.0e-10_dp
    cfg%sim%has_batch_duration = .true.
    cfg%sim%max_step = 10_i32
    cfg%sim%q_floor = 1.0e-40_dp
    cfg%sim%field_solver = 'direct'
    cfg%sim%field_bc_mode = 'periodic2'
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    cfg%sim%box_max = [box_width, box_width, interface_z]
    cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    cfg%field%backend = 'direct'
    cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    cfg%panel%surface_side_policy = 'per_element'
    cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    cfg%periodic2%zero_mode_policy = 'exclude_k0'
    cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    cfg%periodic2%reference_mode_layers = 2_i32
    cfg%periodic2%panel_quadrature_order = 4_i32
    cfg%periodic2%interface_sample_n = 2_i32
    cfg%periodic2%interface_phi_tolerance = 1.0e6_dp
    cfg%periodic2%interface_field_tolerance = 1.0e6_dp
    cfg%outer_plasma%model = 'kinetic_1d'
    cfg%outer_plasma%kinetic_closure = 'ambient_linear_debye'
    cfg%outer_plasma%photoelectron_density_model = 'none'
    cfg%outer_plasma%return_model = 'kinetic_1d_profile_return'
    cfg%outer_plasma%interface_z = interface_z
    cfg%outer_plasma%debye_length = debye_length
    cfg%outer_plasma%thermal_voltage = 10.0_dp
    cfg%coupling%update_mode = 'implicit_mean'
    cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    cfg%coupling%outer_update_stride = 1_i32
    cfg%coupling%field_evolution_timescale = 1.0_dp
    cfg%coupling%max_frozen_field_ratio = 1.0_dp
    cfg%coupling%outer_queue_enabled = .false.

    cfg%n_particle_species = 3_i32
    cfg%particle_species(1) = species_from_defaults()
    cfg%particle_species(1)%species_key = 'ambient_electron'
    cfg%particle_species(1)%source_mode = 'reservoir_face'
    cfg%particle_species(1)%inject_face = 'z_high'
    cfg%particle_species(1)%q_particle = -qe
    cfg%particle_species(1)%m_particle = electron_mass
    cfg%particle_species(1)%w_particle = 1.0e30_dp
    cfg%particle_species(1)%has_w_particle = .true.
    cfg%particle_species(1)%number_density_m3 = 1.0e6_dp
    cfg%particle_species(1)%has_number_density_m3 = .true.
    cfg%particle_species(1)%temperature_ev = 2.0_dp
    cfg%particle_species(1)%has_temperature_ev = .true.
    cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
    cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    cfg%particle_species(1)%pos_high = [box_width, box_width, interface_z]

    cfg%particle_species(2) = species_from_defaults()
    cfg%particle_species(2)%species_key = 'ambient_proton'
    cfg%particle_species(2)%source_mode = 'reservoir_face'
    cfg%particle_species(2)%inject_face = 'z_high'
    cfg%particle_species(2)%q_particle = qe
    cfg%particle_species(2)%m_particle = proton_mass
    cfg%particle_species(2)%w_particle = 1.0e30_dp
    cfg%particle_species(2)%has_w_particle = .true.
    cfg%particle_species(2)%number_density_m3 = 1.0e6_dp
    cfg%particle_species(2)%has_number_density_m3 = .true.
    cfg%particle_species(2)%temperature_ev = 0.0_dp
    cfg%particle_species(2)%has_temperature_ev = .true.
    cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
    cfg%particle_species(2)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    cfg%particle_species(2)%pos_high = [box_width, box_width, interface_z]

    cfg%particle_species(3) = species_from_defaults()
    cfg%particle_species(3)%species_key = 'photoelectron'
    cfg%particle_species(3)%source_mode = 'photo_raycast'
    cfg%particle_species(3)%inject_face = 'z_high'
    cfg%particle_species(3)%q_particle = -qe
    cfg%particle_species(3)%m_particle = electron_mass
    cfg%particle_species(3)%temperature_ev = 2.2_dp
    cfg%particle_species(3)%has_temperature_ev = .true.
    cfg%particle_species(3)%emit_current_density_a_m2 = 1.0e-3_dp
    cfg%particle_species(3)%rays_per_batch = global_rays_per_batch
    cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    cfg%particle_species(3)%has_deposit_opposite_charge_on_emit = .true.
    cfg%particle_species(3)%normal_drift_speed = 0.0_dp
    cfg%particle_species(3)%pos_low = [0.25_dp*box_width, 0.25_dp*box_width, interface_z]
    cfg%particle_species(3)%pos_high = [0.75_dp*box_width, 0.75_dp*box_width, interface_z]
    cfg%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    cfg%particle_species(3)%has_ray_direction = .true.

    call prepare_periodic2_collision_mesh(mesh, cfg%sim)
  end subroutine setup_mpi_implicit_mean_fixture

  subroutine pack_implicit_mean_signature(mesh, ledger, diagnostics, signature)
    type(mesh_type), intent(in) :: mesh
    type(charge_ledger_type), intent(in) :: ledger
    type(electrostatic_diagnostics_type), intent(in) :: diagnostics
    real(dp), intent(out) :: signature(:)
    integer(i32), parameter :: photoelectron_index = 3_i32

    if (size(mesh%q_elem) /= 2 .or. size(signature) /= 14) then
      error stop 'implicit-mean MPI signature storage mismatch'
    end if
    signature = [ &
                mesh%q_elem(1), mesh%q_elem(2), ledger%surface_charge_after, &
                ledger%emitted_from_surface(photoelectron_index), &
                ledger%absorbed_on_surface(photoelectron_index), &
                ledger%escaped_to_infinity(photoelectron_index), &
                ledger%interface_outward_gross(photoelectron_index), &
                ledger%interface_returned_gross(photoelectron_index), &
                diagnostics%implicit_mean_last_returned_outer_flight_time_mean, &
                diagnostics%implicit_mean_last_returning_pe_column_charge_per_area, &
                diagnostics%outer_photoelectron_current_density, diagnostics%outer_total_current_density, &
                diagnostics%interface_potential, diagnostics%interface_field &
                ]
  end subroutine pack_implicit_mean_signature

  subroutine run_mpi_kinetic_outer_test(mpi)
    type(mpi_context), intent(in) :: mpi
    type(mesh_type) :: kinetic_mesh
    type(app_config) :: kinetic_cfg
    type(electrostatic_snapshot_type) :: kinetic_snapshot
    type(kinetic_outer_plasma_options_type) :: kinetic_options
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)
    real(dp) :: local_values(11), global_values(11), held_potential(128)
    integer(i32) :: kinetic_status
    character(len=256) :: kinetic_message

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    tri_v2(:, 1) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 2) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh( &
      kinetic_mesh, tri_v0, tri_v1, tri_v2, q0=[0.01_dp*eps0, 0.01_dp*eps0] &
      )
    kinetic_mesh%elem_vacuum_sign = 1_i32
    kinetic_mesh%vacuum_normals = kinetic_mesh%normals

    call default_app_config(kinetic_cfg)
    kinetic_cfg%sim%field_solver = 'direct'
    kinetic_cfg%sim%field_bc_mode = 'periodic2'
    kinetic_cfg%sim%use_box = .true.
    kinetic_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    kinetic_cfg%sim%box_max = [1.0_dp, 1.0_dp, 0.8_dp]
    kinetic_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    kinetic_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    kinetic_cfg%field%backend = 'direct'
    kinetic_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    kinetic_cfg%panel%surface_side_policy = 'per_element'
    kinetic_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    kinetic_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    kinetic_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    kinetic_cfg%periodic2%reference_mode_layers = 2_i32
    kinetic_cfg%periodic2%panel_quadrature_order = 8_i32
    kinetic_cfg%outer_plasma%model = 'kinetic_1d'
    kinetic_cfg%outer_plasma%interface_z = kinetic_cfg%sim%box_max(3)
    kinetic_cfg%outer_plasma%debye_length = 0.8_dp
    kinetic_cfg%outer_plasma%thermal_voltage = 2.0_dp
    kinetic_cfg%n_particle_species = 2_i32
    kinetic_cfg%particle_species(1) = species_from_defaults()
    kinetic_cfg%particle_species(1)%source_mode = 'reservoir_face'
    kinetic_cfg%particle_species(1)%inject_face = 'z_high'
    kinetic_cfg%particle_species(1)%q_particle = -qe
    kinetic_cfg%particle_species(1)%m_particle = 9.1093837139e-31_dp
    kinetic_cfg%particle_species(1)%number_density_m3 = 1.0e6_dp
    kinetic_cfg%particle_species(1)%has_number_density_m3 = .true.
    kinetic_cfg%particle_species(1)%temperature_ev = 2.0_dp
    kinetic_cfg%particle_species(1)%has_temperature_ev = .true.
    kinetic_cfg%particle_species(2) = species_from_defaults()
    kinetic_cfg%particle_species(2)%source_mode = 'reservoir_face'
    kinetic_cfg%particle_species(2)%inject_face = 'z_high'
    kinetic_cfg%particle_species(2)%q_particle = qe
    kinetic_cfg%particle_species(2)%m_particle = 1.67262192595e-27_dp
    kinetic_cfg%particle_species(2)%number_density_m3 = 1.0e6_dp
    kinetic_cfg%particle_species(2)%has_number_density_m3 = .true.
    kinetic_cfg%particle_species(2)%temperature_ev = 0.0_dp
    kinetic_cfg%particle_species(2)%has_temperature_ev = .true.
    kinetic_cfg%particle_species(2)%drift_velocity = [ &
                                                     0.0_dp, 0.0_dp, &
                                                     -4.0_dp*sqrt( &
                                                     2.0_dp*qe/kinetic_cfg%particle_species(2)%m_particle &
                                                     ) &
                                                     ]

    call resolve_kinetic_outer_options( &
      kinetic_cfg, 0.0_dp, kinetic_options, kinetic_status, kinetic_message &
      )
    call assert_equal_i32(kinetic_status, outer_plasma_ok, &
                          'MPI kinetic options must resolve: '//trim(kinetic_message))
    call kinetic_snapshot%init( &
      kinetic_mesh, kinetic_cfg%sim, kinetic_cfg%field, kinetic_cfg%periodic2, kinetic_cfg%panel, &
      kinetic_cfg%outer_plasma, kinetic_options=kinetic_options, mpi=mpi &
      )
    call kinetic_snapshot%refresh(kinetic_mesh)
    local_values = [ &
                   kinetic_snapshot%outer%interface_potential, kinetic_snapshot%outer%interface_field, &
                   kinetic_snapshot%outer%nonlinear_residual, &
                   kinetic_snapshot%outer%integrated_charge_per_area, &
                   kinetic_snapshot%outer%electron_current_density, &
                   kinetic_snapshot%outer%ion_current_density, &
                   kinetic_snapshot%outer%photoelectron_current_density, &
                   kinetic_snapshot%outer%total_current_density, &
                   real(kinetic_snapshot%outer%nonlinear_iterations, dp), &
                   sum(kinetic_snapshot%outer%potential), sum(kinetic_snapshot%outer%charge_density) &
                   ]
    global_values = local_values
    call mpi_allreduce_sum_real_dp_array(mpi, global_values)
    call assert_true( &
      all(abs(global_values - real(mpi%size, dp)*local_values) <= &
          1.0e-12_dp*max(1.0_dp, abs(global_values))), &
      'kinetic root solve/broadcast must be rank invariant' &
      )
    call assert_true(maxval(abs(kinetic_snapshot%outer%potential)) > 0.0_dp, &
                     'kinetic MPI fixture must exercise a nonzero profile')

    held_potential = kinetic_snapshot%outer%potential
    call kinetic_snapshot%refresh(kinetic_mesh, update_outer=.false.)
    call assert_true(all(kinetic_snapshot%outer%potential == held_potential), &
                     'held kinetic profile must remain bitwise unchanged')
  end subroutine run_mpi_kinetic_outer_test

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
