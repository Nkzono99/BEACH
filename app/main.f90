!> 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。
program main
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_stats, mesh_type, injection_state
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_is_root, mpi_world_size
  use bem_performance_profile, only: perf_configure_from_env, perf_set_output_context, perf_region_begin, &
                                     perf_region_end, perf_write_outputs, perf_region_program_total, perf_region_load_or_init, &
                                     perf_region_history_open, perf_region_write_results, perf_region_write_checkpoint
  use bem_simulator, only: run_absorption_insulator
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file
  use bem_output_writer, only: open_history_writer, print_run_summary, write_result_files, ensure_output_dir
  use bem_app_config, only: app_config, default_app_config, load_app_config, build_mesh_from_config, &
                            seed_particles_from_config
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: app
  type(sim_stats) :: stats
  type(sim_stats) :: initial_stats
  type(injection_state) :: inject_state
  type(mpi_context) :: mpi
  integer :: history_unit
  logical :: history_opened, resumed
  real(dp) :: perf_t0, perf_program_t0

  call perf_configure_from_env()
  call perf_region_begin(perf_region_program_total, perf_program_t0)
  call mpi_initialize(mpi)
  call perf_region_begin(perf_region_load_or_init, perf_t0)
  call load_or_init_run_state(app, mesh, initial_stats, inject_state, resumed, mpi)
  call perf_region_end(perf_region_load_or_init, perf_t0)
  call perf_set_output_context(trim(app%output_dir), app%write_output)
  if (mpi_is_root(mpi)) then
    call perf_region_begin(perf_region_history_open, perf_t0)
    call open_history_writer(app, resumed, history_opened, history_unit)
    call perf_region_end(perf_region_history_open, perf_t0)
  else
    history_opened = .false.
    history_unit = -1
  end if

  if (history_opened) then
    call run_absorption_insulator( &
      mesh, app, stats, history_unit=history_unit, history_stride=app%history_stride, initial_stats=initial_stats, &
      inject_state=inject_state, mpi=mpi &
      )
    close (history_unit)
  else
    call run_absorption_insulator(mesh, app, stats, initial_stats=initial_stats, inject_state=inject_state, mpi=mpi)
  end if

  if (mpi_is_root(mpi)) call print_run_summary(mesh, stats)

  if (app%write_output) then
    call ensure_output_dir(app%output_dir)
    if (mpi_is_root(mpi)) then
      call perf_region_begin(perf_region_write_results, perf_t0)
      call write_result_files(trim(app%output_dir), mesh, stats, app, mpi_world_size=mpi_world_size(mpi))
      call perf_region_end(perf_region_write_results, perf_t0)
    end if
    call perf_region_begin(perf_region_write_checkpoint, perf_t0)
    call write_rng_state_file(trim(app%output_dir), mpi=mpi)
    call write_macro_residuals_file(trim(app%output_dir), inject_state, mpi=mpi)
    call perf_region_end(perf_region_write_checkpoint, perf_t0)
    if (mpi_is_root(mpi)) print '(a,a)', 'results written to ', trim(app%output_dir)
  end if
  call perf_region_end(perf_region_program_total, perf_program_t0)
  call perf_write_outputs(mpi)
  call mpi_shutdown(mpi)

contains

  !> 設定読込・メッシュ構築・再開判定・乱数初期化をまとめて行う。
  !! @param[out] app 読み込み・既定値適用後のアプリ設定。
  !! @param[out] mesh 構築した三角形メッシュ。
  !! @param[out] initial_stats 再開時に引き継ぐ初期統計（新規実行時はゼロ）。
  !! @param[out] inject_state 種別ごとの注入残差状態。
  !! @param[out] resumed チェックポイントから再開した場合に `.true.`。
  subroutine load_or_init_run_state(app, mesh, initial_stats, inject_state, resumed, mpi)
    type(app_config), intent(out) :: app
    type(mesh_type), intent(out) :: mesh
    type(sim_stats), intent(out) :: initial_stats
    type(injection_state), intent(out) :: inject_state
    logical, intent(out) :: resumed
    type(mpi_context), intent(in) :: mpi
    character(len=256) :: cfg_path
    logical :: has_config

    call default_app_config(app)
    call resolve_config_path(cfg_path, has_config)
    if (has_config) then
      call load_app_config(trim(cfg_path), app)
    end if

    call build_mesh_from_config(app, mesh)
    call initialize_injection_state(inject_state, app%n_particle_species)
    initial_stats = sim_stats()
    resumed = .false.
    if (app%resume_output) then
      if (.not. app%write_output) error stop 'output.resume requires output.write_files = true.'
      call load_restart_checkpoint( &
        trim(app%output_dir), mesh, initial_stats, resumed, inject_state, mpi=mpi &
        )
    end if

    if (resumed) then
      if (mpi_is_root(mpi)) then
        print '(a,i0)', 'resuming_from_batches=', initial_stats%batches
        print '(a,i0)', 'resuming_from_processed_particles=', initial_stats%processed_particles
      end if
    else
      call seed_particles_from_config(app, mpi=mpi)
    end if
  end subroutine load_or_init_run_state

  !> 実行時設定ファイルの読み込みパスを決定する。
  !! @param[out] path 読み込む設定ファイルパス（未検出時は空文字）。
  !! @param[out] found 設定ファイルが解決できた場合に `.true.`。
  subroutine resolve_config_path(path, found)
    character(len=*), intent(out) :: path
    logical, intent(out) :: found
    logical :: has_primary
    character(len=*), parameter :: primary_config = 'beach.toml'

    path = ''
    found = .false.
    if (command_argument_count() >= 1) then
      call get_command_argument(1, path)
      if (len_trim(path) == 0) error stop 'Config path argument is empty.'
      found = .true.
      return
    end if

    inquire (file=primary_config, exist=has_primary)
    if (has_primary) then
      path = primary_config
      found = .true.
    end if
  end subroutine resolve_config_path

  !> 種数に合わせて注入状態をゼロ初期化する。
  !! @param[out] state 初期化する注入状態。
  !! @param[in] n_species 粒子種数（`macro_residual` 配列長）。
  subroutine initialize_injection_state(state, n_species)
    type(injection_state), intent(out) :: state
    integer, intent(in) :: n_species

    allocate (state%macro_residual(n_species))
    state%macro_residual = 0.0d0
  end subroutine initialize_injection_state

end program main

