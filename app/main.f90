!> 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。
program main
  use bem_types, only: sim_stats, mesh_type, injection_state
  use bem_simulator, only: run_absorption_insulator
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file
  use bem_app_config, only: app_config, default_app_config, load_app_config, build_mesh_from_config, &
    seed_particles_from_config
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: app
  type(sim_stats) :: stats
  type(sim_stats) :: initial_stats
  type(injection_state) :: inject_state
  integer :: history_unit
  logical :: history_opened, resumed

  call load_or_init_run_state(app, mesh, initial_stats, inject_state, resumed)
  call open_history_writer(app, resumed, history_opened, history_unit)

  if (history_opened) then
    call run_absorption_insulator( &
      mesh, app, stats, history_unit=history_unit, history_stride=app%history_stride, initial_stats=initial_stats, &
      inject_state=inject_state &
    )
    close(history_unit)
  else
    call run_absorption_insulator(mesh, app, stats, initial_stats=initial_stats, inject_state=inject_state)
  end if

  call print_run_summary(mesh, stats)

  if (app%write_output) then
    call write_result_files(trim(app%output_dir), mesh, stats)
    call write_rng_state_file(trim(app%output_dir))
    call write_macro_residuals_file(trim(app%output_dir), inject_state)
    print '(a,a)', 'results written to ', trim(app%output_dir)
  end if

contains

  !> 設定読込・メッシュ構築・再開判定・乱数初期化をまとめて行う。
  !! @param[out] app 読み込み・既定値適用後のアプリ設定。
  !! @param[out] mesh 構築した三角形メッシュ。
  !! @param[out] initial_stats 再開時に引き継ぐ初期統計（新規実行時はゼロ）。
  !! @param[out] inject_state 種別ごとの注入残差状態。
  !! @param[out] resumed チェックポイントから再開した場合に `.true.`。
  subroutine load_or_init_run_state(app, mesh, initial_stats, inject_state, resumed)
    type(app_config), intent(out) :: app
    type(mesh_type), intent(out) :: mesh
    type(sim_stats), intent(out) :: initial_stats
    type(injection_state), intent(out) :: inject_state
    logical, intent(out) :: resumed
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
      call load_restart_checkpoint(trim(app%output_dir), mesh, initial_stats, resumed, inject_state)
    end if

    if (resumed) then
      print '(a,i0)', 'resuming_from_batches=', initial_stats%batches
      print '(a,i0)', 'resuming_from_processed_particles=', initial_stats%processed_particles
    else
      call seed_particles_from_config(app)
    end if
  end subroutine load_or_init_run_state

  !> 実行時設定ファイルの読み込みパスを決定する。
  !! @param[out] path 読み込む設定ファイルパス（未検出時は空文字）。
  !! @param[out] found 設定ファイルが解決できた場合に `.true.`。
  subroutine resolve_config_path(path, found)
    character(len=*), intent(out) :: path
    logical, intent(out) :: found
    logical :: has_primary, has_legacy
    character(len=*), parameter :: primary_config = 'beach.toml'
    character(len=*), parameter :: legacy_config = 'fortran_config.toml'

    path = ''
    found = .false.
    if (command_argument_count() >= 1) then
      call get_command_argument(1, path)
      if (len_trim(path) == 0) error stop 'Config path argument is empty.'
      found = .true.
      return
    end if

    inquire(file=primary_config, exist=has_primary)
    if (has_primary) then
      path = primary_config
      found = .true.
      return
    end if

    inquire(file=legacy_config, exist=has_legacy)
    if (has_legacy) then
      path = legacy_config
      found = .true.
      print '(a)', 'warning: fortran_config.toml is deprecated; rename to beach.toml.'
    end if
  end subroutine resolve_config_path

  !> 種数に合わせて注入状態をゼロ初期化する。
  !! @param[out] state 初期化する注入状態。
  !! @param[in] n_species 粒子種数（`macro_residual` 配列長）。
  subroutine initialize_injection_state(state, n_species)
    type(injection_state), intent(out) :: state
    integer, intent(in) :: n_species

    allocate(state%macro_residual(n_species))
    state%macro_residual = 0.0d0
  end subroutine initialize_injection_state

  !> 履歴 CSV のオープンとヘッダ初期化を行う。
  !! @param[in] app 出力設定を含むアプリ設定。
  !! @param[in] resumed 再開実行かどうか。
  !! @param[out] history_opened 履歴ファイルを開けた場合に `.true.`。
  !! @param[out] history_unit 履歴CSVの出力ユニット番号（未使用時は `-1`）。
  subroutine open_history_writer(app, resumed, history_opened, history_unit)
    type(app_config), intent(in) :: app
    logical, intent(in) :: resumed
    logical, intent(out) :: history_opened
    integer, intent(out) :: history_unit
    character(len=1024) :: history_path
    integer :: ios
    logical :: history_exists

    history_opened = .false.
    history_unit = -1
    if (.not. app%write_output) return
    if (app%history_stride <= 0) return

    call ensure_output_dir(app%output_dir)

    history_path = trim(app%output_dir) // '/charge_history.csv'
    inquire(file=trim(history_path), exist=history_exists)
    if (resumed) then
      open(newunit=history_unit, file=trim(history_path), status='unknown', position='append', action='write', iostat=ios)
    else
      open(newunit=history_unit, file=trim(history_path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open charge history file.'

    if (.not. resumed .or. .not. history_exists) then
      ! 再開時は既存ファイルがある場合だけヘッダ追記を避ける。
      write(history_unit, '(a)') 'batch,processed_particles,rel_change,elem_idx,charge_C'
    end if
    history_opened = .true.
  end subroutine open_history_writer

  !> 実行結果の主要統計を標準出力へ表示する。
  !! @param[in] mesh 実行後のメッシュ情報。
  !! @param[in] stats 実行後の統計値。
  subroutine print_run_summary(mesh, stats)
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats

    print '(a,i0)', 'mesh nelem=', mesh%nelem
    print '(a,i0)', 'processed_particles=', stats%processed_particles
    print '(a,i0)', 'absorbed=', stats%absorbed
    print '(a,i0)', 'escaped=', stats%escaped
    print '(a,i0)', 'batches=', stats%batches
    print '(a,i0)', 'escaped_boundary=', stats%escaped_boundary
    print '(a,i0)', 'survived_max_step=', stats%survived_max_step
    print '(a,es12.4)', 'last_rel_change=', stats%last_rel_change
    print '(a,es12.4)', 'field_time_s=', stats%field_time_s
    print '(a,es12.4)', 'push_time_s=', stats%push_time_s
    print '(a,es12.4)', 'collision_time_s=', stats%collision_time_s
    print '(a,*(es12.4,1x))', 'mesh charges=', mesh%q_elem
  end subroutine print_run_summary

  !> 解析結果を `summary.txt` / `charges.csv` / `mesh_triangles.csv` として保存する。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 書き出し対象のメッシュ。
  !! @param[in] stats 書き出し対象の統計値。
  subroutine write_result_files(out_dir, mesh, stats)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats

    call ensure_output_dir(out_dir)
    call write_summary_file(out_dir, mesh, stats)
    call write_charges_file(out_dir, mesh)
    call write_mesh_file(out_dir, mesh)
  end subroutine write_result_files

  !> 出力ディレクトリを作成する。
  !! @param[in] out_dir 作成対象ディレクトリのパス。
  subroutine ensure_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir
    character(len=1024) :: cmd
    integer :: ios

    cmd = 'mkdir -p "' // trim(out_dir) // '"'
    call execute_command_line(trim(cmd), wait=.true., exitstat=ios)
    if (ios /= 0) error stop 'Failed to create output directory.'
  end subroutine ensure_output_dir

  !> 実行統計を `summary.txt` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh メッシュ情報（要素数を書き出す）。
  !! @param[in] stats 実行統計。
  subroutine write_summary_file(out_dir, mesh, stats)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    character(len=1024) :: summary_path
    integer :: u, ios

    summary_path = trim(out_dir) // '/summary.txt'
    open(newunit=u, file=trim(summary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary file.'
    write(u, '(a,i0)') 'mesh_nelem=', mesh%nelem
    write(u, '(a,i0)') 'processed_particles=', stats%processed_particles
    write(u, '(a,i0)') 'absorbed=', stats%absorbed
    write(u, '(a,i0)') 'escaped=', stats%escaped
    write(u, '(a,i0)') 'batches=', stats%batches
    write(u, '(a,i0)') 'escaped_boundary=', stats%escaped_boundary
    write(u, '(a,i0)') 'survived_max_step=', stats%survived_max_step
    write(u, '(a,es24.16)') 'last_rel_change=', stats%last_rel_change
    write(u, '(a,es24.16)') 'field_time_s=', stats%field_time_s
    write(u, '(a,es24.16)') 'push_time_s=', stats%push_time_s
    write(u, '(a,es24.16)') 'collision_time_s=', stats%collision_time_s
    close(u)
  end subroutine write_summary_file

  !> 要素電荷を `charges.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素電荷を含むメッシュ情報。
  subroutine write_charges_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    character(len=1024) :: charges_path
    integer :: u, ios, i

    charges_path = trim(out_dir) // '/charges.csv'
    open(newunit=u, file=trim(charges_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charges file.'
    write(u, '(a)') 'elem_idx,charge_C'
    do i = 1, mesh%nelem
      write(u, '(i0,a,es24.16)') i, ',', mesh%q_elem(i)
    end do
    close(u)
  end subroutine write_charges_file

  !> 三角形メッシュを `mesh_triangles.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 頂点座標と要素電荷を含むメッシュ情報。
  subroutine write_mesh_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    character(len=1024) :: mesh_path
    integer :: u, ios, i

    mesh_path = trim(out_dir) // '/mesh_triangles.csv'
    open(newunit=u, file=trim(mesh_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh file.'
    write(u, '(a)') 'elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C'
    do i = 1, mesh%nelem
      write(u, '(i0,10(a,es24.16))') i, ',', mesh%v0(1, i), ',', mesh%v0(2, i), ',', mesh%v0(3, i), &
                                      ',', mesh%v1(1, i), ',', mesh%v1(2, i), ',', mesh%v1(3, i), &
                                      ',', mesh%v2(1, i), ',', mesh%v2(2, i), ',', mesh%v2(3, i), ',', mesh%q_elem(i)
    end do
    close(u)
  end subroutine write_mesh_file

end program main
