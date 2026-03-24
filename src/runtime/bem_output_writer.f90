!> 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。
module bem_output_writer
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_stats
  use bem_app_config_types, only: app_config
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: open_history_writer
  public :: print_run_summary
  public :: write_result_files
  public :: ensure_output_dir

contains

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

    history_path = trim(app%output_dir)//'/charge_history.csv'
    inquire (file=trim(history_path), exist=history_exists)
    if (resumed) then
      open (newunit=history_unit, file=trim(history_path), status='unknown', position='append', action='write', iostat=ios)
    else
      open (newunit=history_unit, file=trim(history_path), status='replace', action='write', iostat=ios)
    end if
    if (ios /= 0) error stop 'Failed to open charge history file.'

    if (.not. resumed .or. .not. history_exists) then
      ! 再開時は既存ファイルがある場合だけヘッダ追記を避ける。
      write (history_unit, '(a)') 'batch,processed_particles,rel_change,elem_idx,charge_C'
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
    print '(a,*(es12.4,1x))', 'mesh charges=', mesh%q_elem
  end subroutine print_run_summary

  !> 解析結果を `summary.txt` / `charges.csv` / `mesh_triangles.csv` などとして保存する。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 書き出し対象のメッシュ。
  !! @param[in] stats 書き出し対象の統計値。
  !! @param[in] cfg 出力設定を含むアプリ設定。
  subroutine write_result_files(out_dir, mesh, stats, cfg, mpi_world_size, mesh_potential_v)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_world_size
    real(dp), intent(in), optional :: mesh_potential_v(:)

    call ensure_output_dir(out_dir)
    call write_summary_file(out_dir, mesh, stats, mpi_world_size=mpi_world_size)
    call write_charges_file(out_dir, mesh)
    if (cfg%write_mesh_potential) then
      if (.not. present(mesh_potential_v)) then
        error stop 'write_result_files: mesh_potential_v is required when write_mesh_potential is enabled.'
      end if
      call write_mesh_potential_file(out_dir, mesh, mesh_potential_v)
    end if
    call write_mesh_file(out_dir, mesh)
    call write_mesh_sources_file(out_dir, mesh, cfg)
  end subroutine write_result_files

  !> 出力ディレクトリを作成する。
  !! @param[in] out_dir 作成対象ディレクトリのパス。
  subroutine ensure_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir
    character(len=1024) :: cmd
    integer :: ios

    cmd = 'mkdir -p "'//trim(out_dir)//'"'
    call execute_command_line(trim(cmd), wait=.true., exitstat=ios)
    if (ios /= 0) error stop 'Failed to create output directory.'
  end subroutine ensure_output_dir

  !> 実行統計を `summary.txt` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh メッシュ情報（要素数を書き出す）。
  !! @param[in] stats 実行統計。
  subroutine write_summary_file(out_dir, mesh, stats, mpi_world_size)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    integer(i32), intent(in), optional :: mpi_world_size
    character(len=1024) :: summary_path
    integer :: u, ios
    integer(i32) :: world_size

    summary_path = trim(out_dir)//'/summary.txt'
    open (newunit=u, file=trim(summary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary file.'
    world_size = 1_i32
    if (present(mpi_world_size)) world_size = max(1_i32, mpi_world_size)
    write (u, '(a,i0)') 'mesh_nelem=', mesh%nelem
    write (u, '(a,i0)') 'mesh_count=', max(1_i32, maxval(mesh%elem_mesh_id))
    write (u, '(a,i0)') 'mpi_world_size=', world_size
    write (u, '(a,i0)') 'processed_particles=', stats%processed_particles
    write (u, '(a,i0)') 'absorbed=', stats%absorbed
    write (u, '(a,i0)') 'escaped=', stats%escaped
    write (u, '(a,i0)') 'batches=', stats%batches
    write (u, '(a,i0)') 'escaped_boundary=', stats%escaped_boundary
    write (u, '(a,i0)') 'survived_max_step=', stats%survived_max_step
    write (u, '(a,es24.16)') 'last_rel_change=', stats%last_rel_change
    close (u)
  end subroutine write_summary_file

  !> 要素電荷を `charges.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素電荷を含むメッシュ情報。
  subroutine write_charges_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    character(len=1024) :: charges_path
    integer :: u, ios, i

    charges_path = trim(out_dir)//'/charges.csv'
    open (newunit=u, file=trim(charges_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charges file.'
    write (u, '(a)') 'elem_idx,charge_C'
    do i = 1, mesh%nelem
      write (u, '(i0,a,es24.16)') i, ',', mesh%q_elem(i)
    end do
    close (u)
  end subroutine write_charges_file

  !> 事前計算済み電位を `mesh_potential.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素情報（要素数の検証用）。
  !! @param[in] potential_v 各要素重心での電位 [V]。
  subroutine write_mesh_potential_file(out_dir, mesh, potential_v)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: potential_v(:)
    character(len=1024) :: potential_path
    integer :: u, ios, i

    if (size(potential_v) /= mesh%nelem) error stop 'precomputed mesh potential size mismatch.'

    potential_path = trim(out_dir)//'/mesh_potential.csv'
    open (newunit=u, file=trim(potential_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh_potential.csv.'
    write (u, '(a)') 'elem_idx,potential_V'
    do i = 1, mesh%nelem
      write (u, '(i0,a,es24.16)') i, ',', potential_v(i)
    end do
    close (u)
  end subroutine write_mesh_potential_file

  !> 三角形メッシュを `mesh_triangles.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 頂点座標と要素電荷を含むメッシュ情報。
  subroutine write_mesh_file(out_dir, mesh)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    character(len=1024) :: mesh_path
    integer :: u, ios, i

    mesh_path = trim(out_dir)//'/mesh_triangles.csv'
    open (newunit=u, file=trim(mesh_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh file.'
    write (u, '(a)') 'elem_idx,v0x,v0y,v0z,v1x,v1y,v1z,v2x,v2y,v2z,charge_C,mesh_id'
    do i = 1, mesh%nelem
      write (u, '(i0,10(a,es24.16),a,i0)') i, ',', mesh%v0(1, i), ',', mesh%v0(2, i), ',', mesh%v0(3, i), &
        ',', mesh%v1(1, i), ',', mesh%v1(2, i), ',', mesh%v1(3, i), &
        ',', mesh%v2(1, i), ',', mesh%v2(2, i), ',', mesh%v2(3, i), ',', &
        mesh%q_elem(i), ',', mesh%elem_mesh_id(i)
    end do
    close (u)
  end subroutine write_mesh_file

  !> メッシュ識別情報を `mesh_sources.csv` に書き出す。
  !! @param[in] out_dir 出力先ディレクトリ。
  !! @param[in] mesh 要素ごとの `mesh_id` を含むメッシュ情報。
  !! @param[in] cfg 元の入力設定。
  subroutine write_mesh_sources_file(out_dir, mesh, cfg)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(app_config), intent(in) :: cfg
    character(len=1024) :: path
    character(len=16) :: mode_key, source_kind, template_kind
    logical :: has_obj
    integer :: u, ios, i, mesh_id
    integer(i32) :: elem_count

    path = trim(out_dir)//'/mesh_sources.csv'
    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open mesh_sources.csv.'
    write (u, '(a)') 'mesh_id,source_kind,template_kind,elem_count'

    mode_key = trim(lower_ascii(cfg%mesh_mode))
    if (mode_key == 'obj') then
      source_kind = 'obj'
    else if (mode_key == 'template') then
      source_kind = 'template'
    else
      inquire (file=trim(cfg%obj_path), exist=has_obj)
      if (has_obj) then
        source_kind = 'obj'
      else
        source_kind = 'template'
      end if
    end if

    if (source_kind == 'obj') then
      write (u, '(i0,a,a,a,a,a,i0)') 1, ',', 'obj', ',', 'obj', ',', mesh%nelem
      close (u)
      return
    end if

    mesh_id = 0
    do i = 1, size(cfg%templates)
      if (.not. cfg%templates(i)%enabled) cycle
      mesh_id = mesh_id + 1
      template_kind = trim(lower_ascii(cfg%templates(i)%kind))
      elem_count = int(count(mesh%elem_mesh_id == mesh_id), kind=i32)
      write (u, '(i0,a,a,a,a,a,i0)') mesh_id, ',', 'template', ',', trim(template_kind), ',', elem_count
    end do
    close (u)
  end subroutine write_mesh_sources_file

end module bem_output_writer
