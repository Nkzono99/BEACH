!> 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。
program main
!$ use omp_lib, only: omp_sched_dynamic, omp_set_schedule
  use, intrinsic :: iso_fortran_env, only: error_unit
  use bem_kinds, only: dp, i32
  use bem_version, only: beach_build_info, beach_version
  use bem_types, only: sim_stats, mesh_type, injection_state
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_is_root, mpi_world_barrier, mpi_world_size
  use bem_performance_profile, only: perf_configure_from_env, perf_set_output_context, perf_region_begin, &
                                     perf_region_end, perf_write_outputs, perf_region_program_total, perf_region_load_or_init, &
                                     perf_region_history_open, perf_region_write_results, perf_region_write_checkpoint
  use bem_simulator, only: run_absorption_insulator
  use bem_restart, only: load_restart_checkpoint, write_macro_residuals_file
  use bem_periodic_checkpoint, only: resolve_latest_checkpoint_dir
  use bem_checkpoint_contract, only: publish_checkpoint_manifest
  use bem_output_writer, only: print_run_summary, write_result_files, ensure_output_dir
  use bem_run_outputs, only: run_outputs_type, open_run_outputs, close_run_outputs
  use bem_app_config, only: app_config, default_app_config, load_app_config, build_mesh_from_config, &
                            seed_particles_from_config
  use bem_mesh, only: prepare_built_periodic2_collision_mesh
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type, electrostatic_snapshot_type
  use bem_physics_config_types, only: field_physics_config, panel_kernel_config, derive_field_panel_config
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: app
  type(sim_stats) :: stats
  type(sim_stats) :: initial_stats
  type(injection_state) :: inject_state
  type(charge_ledger_type) :: charge_ledger
  type(electrostatic_diagnostics_type) :: electrostatic_diagnostics
  type(mpi_context) :: mpi
  type(run_outputs_type) :: run_outputs
  logical :: resumed
  real(dp) :: perf_t0, perf_program_t0
  logical :: prepare_cache_only = .false.
!$ integer :: schedule_env_status

  call handle_early_cli()

  ! Balance complete particle tracks by default; an explicit OpenMP schedule takes precedence.
!$ call get_environment_variable('OMP_SCHEDULE', status=schedule_env_status)
!$ if (schedule_env_status == 1) call omp_set_schedule(omp_sched_dynamic, 1)

  call perf_configure_from_env()
  call perf_region_begin(perf_region_program_total, perf_program_t0)
  call mpi_initialize(mpi)
  call perf_region_begin(perf_region_load_or_init, perf_t0)
  call load_or_init_run_state(app, mesh, initial_stats, inject_state, charge_ledger, resumed, mpi)
  call perf_region_end(perf_region_load_or_init, perf_t0)
  if (prepare_cache_only) then
    call prepare_periodic_cache_and_stop()
  end if
  call perf_set_output_context(trim(app%output_dir), app%write_output)
  if (mpi_is_root(mpi)) call perf_region_begin(perf_region_history_open, perf_t0)
  call open_run_outputs(app, resumed, initial_stats%batches, mpi_is_root(mpi), run_outputs)
  if (mpi_is_root(mpi)) call perf_region_end(perf_region_history_open, perf_t0)

  call run_absorption_insulator( &
    mesh, app, stats, initial_stats=initial_stats, inject_state=inject_state, mpi=mpi, charge_ledger=charge_ledger, &
    electrostatic_diagnostics=electrostatic_diagnostics, outputs=run_outputs &
    )
  call close_run_outputs(run_outputs)

  if (mpi_is_root(mpi)) call print_run_summary(mesh, stats)

  if (app%write_output) then
    call ensure_output_dir(app%output_dir)
    if (mpi_is_root(mpi)) then
      call perf_region_begin(perf_region_write_results, perf_t0)
      if (allocated(run_outputs%mesh_potential_v)) then
        call write_result_files( &
          trim(app%output_dir), mesh, stats, app, mpi_world_size=mpi_world_size(mpi), &
          mesh_potential_v=run_outputs%mesh_potential_v, &
          charge_ledger=charge_ledger, electrostatic_diagnostics=electrostatic_diagnostics, &
          identity=run_outputs%summary_identity &
          )
      else
        call write_result_files( &
          trim(app%output_dir), mesh, stats, app, mpi_world_size=mpi_world_size(mpi), charge_ledger=charge_ledger, &
          electrostatic_diagnostics=electrostatic_diagnostics, identity=run_outputs%summary_identity &
          )
      end if
      call perf_region_end(perf_region_write_results, perf_t0)
    end if
    ! rootが旧完了manifestを無効化してから、全rankでrestart stateの公開へ進む。
    call mpi_world_barrier(mpi)
    call perf_region_begin(perf_region_write_checkpoint, perf_t0)
    call write_macro_residuals_file(trim(app%output_dir), inject_state, mpi=mpi)
    call mpi_world_barrier(mpi)
    if (mpi_is_root(mpi)) then
      call publish_checkpoint_manifest( &
        trim(app%output_dir), stats%batches, mpi_world_size(mpi), allocated(inject_state%macro_residual), .true. &
        )
    end if
    call mpi_world_barrier(mpi)
    call perf_region_end(perf_region_write_checkpoint, perf_t0)
    if (mpi_is_root(mpi)) print '(a,a)', 'results written to ', trim(app%output_dir)
  end if
  call perf_region_end(perf_region_program_total, perf_program_t0)
  call perf_write_outputs(mpi)
  call mpi_shutdown(mpi)

contains

  !> MPI 初期化や計算状態の構築前に完結する CLI option を処理する。
  subroutine handle_early_cli()
    character(len=256) :: arg
    character(len=:), allocatable :: config_path
    integer :: path_length

    if (command_argument_count() < 1) return
    call get_command_argument(1, arg)
    select case (trim(arg))
    case ('--version', '-V')
      print '(a)', trim(beach_version)
      stop
    case ('--build-info')
      print '(a)', beach_build_info
      stop
    case ('--check-config')
      if (command_argument_count() /= 2) error stop 'usage: beach --check-config beach.toml'
      call get_command_argument(2, length=path_length)
      allocate (character(len=path_length) :: config_path)
      call get_command_argument(2, config_path)
      call default_app_config(app)
      call load_app_config(config_path, app)
      print '(a,a)', 'config=', config_path
      print '(a)', 'checks=toml,semantic'
      print '(a)', 'status=ok'
      stop
    case ('--prepare-periodic-cache')
      if (command_argument_count() /= 2) error stop 'usage: beach --prepare-periodic-cache beach.toml'
      prepare_cache_only = .true.
    case ('--help', '-h')
      print '(a)', 'usage: beach [beach.toml]'
      print '(a)', '       beach --check-config beach.toml'
      print '(a)', '       beach --prepare-periodic-cache beach.toml   (build/publish the periodic2 operator cache, then exit)'
      print '(a)', '       beach --version'
      print '(a)', '       beach --build-info'
      stop
    end select
  end subroutine handle_early_cli

  !> 設定読込・メッシュ構築・再開判定・乱数初期化をまとめて行う。
  !! @param[out] app 読み込み・既定値適用後のアプリ設定。
  !! @param[out] mesh 構築した三角形メッシュ。
  !! @param[out] initial_stats 再開時に引き継ぐ初期統計（新規実行時はゼロ）。
  !! @param[out] inject_state 種別ごとの注入残差状態。
  !! @param[out] resumed チェックポイントから再開した場合に `.true.`。
  subroutine load_or_init_run_state(app, mesh, initial_stats, inject_state, charge_ledger, resumed, mpi)
    type(app_config), intent(out) :: app
    type(mesh_type), intent(out) :: mesh
    type(sim_stats), intent(out) :: initial_stats
    type(injection_state), intent(out) :: inject_state
    type(charge_ledger_type), intent(out) :: charge_ledger
    logical, intent(out) :: resumed
    type(mpi_context), intent(in) :: mpi
    character(len=256) :: cfg_path
    character(len=256) :: restart_dir
    logical :: has_config

    call default_app_config(app)
    call resolve_config_path(cfg_path, has_config)
    if (has_config) then
      call load_app_config(trim(cfg_path), app)
    end if

    ! The cache-only mode never reads checkpoints; it only needs the mesh and field settings.
    if (prepare_cache_only) app%resume_output = .false.
    call build_mesh_from_config(app, mesh)
    call prepare_built_periodic2_collision_mesh(mesh, app%sim)
    call initialize_injection_state(inject_state, app%n_particle_species)
    initial_stats = sim_stats()
    resumed = .false.
    if (app%resume_output) then
      restart_dir = app%output_dir
      if (len_trim(app%output_restart_from) > 0) restart_dir = app%output_restart_from
      call resolve_latest_checkpoint_dir(trim(restart_dir), restart_dir)
      call load_restart_checkpoint( &
        trim(restart_dir), mesh, initial_stats, resumed, inject_state, &
        mpi=mpi, require_checkpoint=.true., app=app, charge_ledger=charge_ledger &
        )
    end if

    if (resumed) then
      if (mpi_is_root(mpi)) then
        print '(a,a)', 'resuming_from_dir=', trim(restart_dir)
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
    if (prepare_cache_only) then
      call get_command_argument(2, path)
      if (len_trim(path) == 0) error stop 'Config path argument is empty.'
      found = .true.
      return
    end if
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

  !> periodic2 の遠方補正 operator を作成（または cache から読込）して publish し、終了する。
  !! 本番 run より多い MPI rank で cold build だけを済ませる用途。target anchor が rank に分配されるので
  !! rank 数に応じて短くなる。結果ファイルや checkpoint は書かない。
  subroutine prepare_periodic_cache_and_stop()
    type(electrostatic_snapshot_type) :: snapshot
    type(field_physics_config) :: field_config
    type(panel_kernel_config) :: panel_config
    type(electrostatic_diagnostics_type) :: diagnostics
    integer(8) :: clock_start, clock_end, clock_rate

    call system_clock(clock_start, clock_rate)
    call derive_field_panel_config(app%sim, field_config, panel_config)
    call snapshot%init(mesh, app%sim, field_config, app%periodic2, panel_config, mpi=mpi)
    call snapshot%get_diagnostics(diagnostics)
    call mpi_world_barrier(mpi)
    call system_clock(clock_end)
    if (mpi_is_root(mpi)) then
      print '(a,es12.4)', 'prepare_elapsed_s=', real(clock_end - clock_start, dp)/real(max(clock_rate, 1_8), dp)
      print '(a,a)', 'periodic2_status=', trim(diagnostics%status)
      print '(a,l1)', 'periodic2_split_periodic_active=', diagnostics%split_periodic_active
      print '(a,l1)', 'periodic2_cache_hit=', diagnostics%periodic_cache_hit
      print '(a,i0)', 'periodic2_operator_build_count=', diagnostics%periodic_operator_build_count
      print '(a,a)', 'periodic2_cache_fingerprint=', trim(diagnostics%periodic_cache_fingerprint)
      print '(a,a)', 'periodic2_cache_path=', trim(diagnostics%periodic_cache_path)
      print '(a,i0)', 'mpi_world_size=', mpi_world_size(mpi)
    end if
    call mpi_shutdown(mpi)
    stop
  end subroutine prepare_periodic_cache_and_stop

  !> 種数に合わせて注入状態をゼロ初期化する。
  !! @param[out] state 初期化する注入状態。
  !! @param[in] n_species 粒子種数（source残差と6面boundary流入残差の種方向配列長）。
  subroutine initialize_injection_state(state, n_species)
    type(injection_state), intent(out) :: state
    integer, intent(in) :: n_species

    allocate (state%macro_residual(n_species))
    allocate (state%boundary_macro_residual(6, n_species))
    state%macro_residual = 0.0d0
    state%boundary_macro_residual = 0.0d0
  end subroutine initialize_injection_state

end program main
