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
  use bem_types, only: mesh_type, particles_soa, sim_stats, injection_state, bc_open, bc_reflect, bc_periodic
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_diagnostics_type
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

  call test_init(6)

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
