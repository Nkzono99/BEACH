!> 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。
program main
  use bem_kinds, only: i32, dp
  use bem_types, only: sim_stats, mesh_type, particles_soa
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, load_app_config, build_mesh_from_config, init_particles_from_config
  implicit none

  type(mesh_type) :: mesh
  type(particles_soa) :: pcls
  type(app_config) :: app
  type(sim_stats) :: stats
  character(len=256) :: cfg_path
  character(len=1024) :: cmd, history_path
  integer :: history_unit, ios
  logical :: history_opened

  call default_app_config(app)

  if (command_argument_count() >= 1) then
    call get_command_argument(1, cfg_path)
    call load_app_config(trim(cfg_path), app)
  end if

  call build_mesh_from_config(app, mesh)
  call init_particles_from_config(app, pcls)

  history_opened = .false.
  if (app%write_output .and. app%history_stride > 0) then
    cmd = 'mkdir -p "' // trim(app%output_dir) // '"'
    call execute_command_line(trim(cmd), wait=.true., exitstat=ios)
    if (ios /= 0) error stop 'Failed to create output directory.'

    history_path = trim(app%output_dir) // '/charge_history.csv'
    open(newunit=history_unit, file=trim(history_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charge history file.'
    write(history_unit, '(a)') 'batch,processed_particles,rel_change,elem_idx,charge_C'
    history_opened = .true.
  end if

  if (history_opened) then
    call run_absorption_insulator(mesh, app%sim, pcls, stats, history_unit=history_unit, history_stride=app%history_stride)
    close(history_unit)
  else
    call run_absorption_insulator(mesh, app%sim, pcls, stats)
  end if

  print '(a,i0)', 'mesh nelem=', mesh%nelem
  print '(a,i0)', 'processed_particles=', stats%processed_particles
  print '(a,i0)', 'absorbed=', stats%absorbed
  print '(a,i0)', 'escaped=', stats%escaped
  print '(a,i0)', 'batches=', stats%batches
  print '(a,es12.4)', 'last_rel_change=', stats%last_rel_change
  print '(a,*(es12.4,1x))', 'mesh charges=', mesh%q_elem

  if (app%write_output) then
    call write_result_files(trim(app%output_dir), mesh, stats)
    print '(a,a)', 'results written to ', trim(app%output_dir)
  end if

contains

  !> 解析結果を `summary.txt` / `charges.csv` / `mesh_triangles.csv`（履歴は実行中に `charge_history.csv` へ逐次書込） として出力ディレクトリへ保存する。
  !! @param[in] out_dir 入力引数。
  !! @param[in] mesh 入力引数。
  !! @param[in] stats 入力引数。
  subroutine write_result_files(out_dir, mesh, stats)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    character(len=1024) :: cmd, summary_path, charges_path, mesh_path
    integer :: u, ios, i

    cmd = 'mkdir -p "' // trim(out_dir) // '"'
    call execute_command_line(trim(cmd), wait=.true., exitstat=ios)
    if (ios /= 0) error stop 'Failed to create output directory.'

    summary_path = trim(out_dir) // '/summary.txt'
    open(newunit=u, file=trim(summary_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open summary file.'
    write(u, '(a,i0)') 'mesh_nelem=', mesh%nelem
    write(u, '(a,i0)') 'processed_particles=', stats%processed_particles
    write(u, '(a,i0)') 'absorbed=', stats%absorbed
    write(u, '(a,i0)') 'escaped=', stats%escaped
    write(u, '(a,i0)') 'batches=', stats%batches
    write(u, '(a,es24.16)') 'last_rel_change=', stats%last_rel_change
    close(u)

    charges_path = trim(out_dir) // '/charges.csv'
    open(newunit=u, file=trim(charges_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Failed to open charges file.'
    write(u, '(a)') 'elem_idx,charge_C'
    do i = 1, mesh%nelem
      write(u, '(i0,a,es24.16)') i, ',', mesh%q_elem(i)
    end do
    close(u)

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

  end subroutine write_result_files

end program main
