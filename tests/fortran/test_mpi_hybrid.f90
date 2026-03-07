!> MPI+OpenMPハイブリッド実行時の集約・rank別resumeファイルを検証する。
program test_mpi_hybrid
  use bem_kinds, only: dp, i32
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_is_root, mpi_world_barrier
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file
  use bem_types, only: mesh_type, sim_stats, injection_state
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists, ensure_directory, &
    remove_empty_directory
  implicit none

  type(mpi_context) :: mpi
  type(mesh_type) :: mesh, mesh_restart
  type(app_config) :: cfg
  type(sim_stats) :: stats, stats_restart
  type(injection_state) :: state, state_restart
  logical :: has_restart
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  integer :: u, ios
  character(len=256) :: line
  integer(i32) :: n_lines
  character(len=*), parameter :: history_path = 'test_mpi_hybrid_history_tmp.csv'
  character(len=*), parameter :: out_dir = 'test_mpi_hybrid_restart_tmp'
  character(len=1024) :: rng_path, residual_path

  call mpi_initialize(mpi)

  if (mpi_is_root(mpi)) then
    call delete_file_if_exists(history_path)
    call delete_file_if_exists(out_dir // '/summary.txt')
    call delete_file_if_exists(out_dir // '/charges.csv')
    call delete_file_if_exists(out_dir // '/rng_state.txt')
    call delete_file_if_exists(out_dir // '/macro_residuals.csv')
  end if
  call mpi_world_barrier(mpi)

  rng_path = restart_rng_path(out_dir, mpi%rank, mpi%size)
  residual_path = restart_residual_path(out_dir, mpi%rank, mpi%size)
  call delete_file_if_exists(rng_path)
  call delete_file_if_exists(residual_path)
  call mpi_world_barrier(mpi)

  v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
  v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  call init_mesh(mesh, v0, v1, v2)

  call default_app_config(cfg)
  cfg%sim%rng_seed = 2468_i32
  cfg%sim%batch_count = 1_i32
  cfg%sim%dt = 1.0d0
  cfg%sim%max_step = 1_i32
  cfg%sim%softening = 1.0d-6
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

  call seed_particles_from_config(cfg, mpi_rank=mpi%rank, mpi_size=mpi%size)

  if (mpi_is_root(mpi)) then
    open(newunit=u, file=history_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open MPI hybrid history fixture'
    call run_absorption_insulator(mesh, cfg, stats, history_unit=u, history_stride=1_i32, mpi=mpi)
    close(u)
  else
    call run_absorption_insulator(mesh, cfg, stats, mpi=mpi)
  end if

  call assert_equal_i32(stats%processed_particles, 4_i32, 'mpi processed_particles mismatch')
  call assert_equal_i32(stats%absorbed, 4_i32, 'mpi absorbed mismatch')
  call assert_equal_i32(stats%escaped, 0_i32, 'mpi escaped mismatch')
  call assert_equal_i32(stats%batches, 1_i32, 'mpi batches mismatch')
  call assert_close_dp(mesh%q_elem(1), 4.0d0, 1.0d-12, 'mpi deposited charge mismatch')

  if (mpi_is_root(mpi)) then
    n_lines = 0_i32
    open(newunit=u, file=history_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to read MPI hybrid history fixture'
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      n_lines = n_lines + 1_i32
    end do
    close(u)
    call assert_equal_i32(n_lines, 1_i32, 'mpi history snapshot line count mismatch')
    call delete_file_if_exists(history_path)
  end if
  call mpi_world_barrier(mpi)

  call ensure_directory(out_dir)
  if (mpi_is_root(mpi)) then
    call write_summary_fixture(out_dir, mpi%size)
    call write_charges_fixture(out_dir)
  end if
  call mpi_world_barrier(mpi)

  allocate(state%macro_residual(1))
  state%macro_residual(1) = 0.25d0 + real(mpi%rank, dp)
  call write_rng_state_file(out_dir, mpi_rank=mpi%rank, mpi_size=mpi%size)
  call write_macro_residuals_file(out_dir, state, mpi_rank=mpi%rank, mpi_size=mpi%size)
  call mpi_world_barrier(mpi)

  call init_mesh(mesh_restart, v0, v1, v2)
  allocate(state_restart%macro_residual(1))
  state_restart%macro_residual = 0.0d0
  call load_restart_checkpoint( &
    out_dir, mesh_restart, stats_restart, has_restart, state_restart, mpi_rank=mpi%rank, mpi_size=mpi%size &
  )
  call assert_true(has_restart, 'mpi restart should be detected')
  call assert_equal_i32(stats_restart%processed_particles, 8_i32, 'mpi restart processed_particles mismatch')
  call assert_equal_i32(stats_restart%batches, 2_i32, 'mpi restart batches mismatch')
  call assert_close_dp(mesh_restart%q_elem(1), 2.0d0, 1.0d-12, 'mpi restart charge mismatch')
  call assert_close_dp(state_restart%macro_residual(1), 0.25d0 + real(mpi%rank, dp), 1.0d-12, 'mpi residual mismatch')

  call delete_file_if_exists(rng_path)
  call delete_file_if_exists(residual_path)
  call mpi_world_barrier(mpi)
  if (mpi_is_root(mpi)) then
    call delete_file_if_exists(out_dir // '/summary.txt')
    call delete_file_if_exists(out_dir // '/charges.csv')
    call delete_file_if_exists(out_dir // '/rng_state.txt')
    call delete_file_if_exists(out_dir // '/macro_residuals.csv')
    call remove_empty_directory(out_dir)
  end if

  call mpi_shutdown(mpi)

contains

  subroutine write_summary_fixture(dir_path, mpi_world_size)
    character(len=*), intent(in) :: dir_path
    integer(i32), intent(in) :: mpi_world_size
    integer :: u, ios

    open(newunit=u, file=trim(dir_path) // '/summary.txt', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open MPI summary fixture'
    write(u, '(a)') 'mesh_nelem=1'
    write(u, '(a,i0)') 'mpi_world_size=', mpi_world_size
    write(u, '(a)') 'processed_particles=8'
    write(u, '(a)') 'absorbed=8'
    write(u, '(a)') 'escaped=0'
    write(u, '(a)') 'batches=2'
    write(u, '(a)') 'escaped_boundary=0'
    write(u, '(a)') 'survived_max_step=0'
    write(u, '(a)') 'last_rel_change=1.0e-4'
    close(u)
  end subroutine write_summary_fixture

  subroutine write_charges_fixture(dir_path)
    character(len=*), intent(in) :: dir_path
    integer :: u, ios

    open(newunit=u, file=trim(dir_path) // '/charges.csv', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open MPI charges fixture'
    write(u, '(a)') 'elem_idx,charge_C'
    write(u, '(a)') '1,2.0'
    close(u)
  end subroutine write_charges_fixture

  function restart_rng_path(dir_path, rank, size) result(path)
    character(len=*), intent(in) :: dir_path
    integer(i32), intent(in) :: rank, size
    character(len=1024) :: path

    if (size <= 1_i32) then
      path = trim(dir_path) // '/rng_state.txt'
    else
      write(path, '(a,a,i5.5,a)') trim(dir_path), '/rng_state_rank', rank, '.txt'
    end if
  end function restart_rng_path

  function restart_residual_path(dir_path, rank, size) result(path)
    character(len=*), intent(in) :: dir_path
    integer(i32), intent(in) :: rank, size
    character(len=1024) :: path

    if (size <= 1_i32) then
      path = trim(dir_path) // '/macro_residuals.csv'
    else
      write(path, '(a,a,i5.5,a)') trim(dir_path), '/macro_residuals_rank', rank, '.csv'
    end if
  end function restart_residual_path

end program test_mpi_hybrid
