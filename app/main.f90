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

  call default_app_config(app)

  if (command_argument_count() >= 1) then
    call get_command_argument(1, cfg_path)
    call load_app_config(trim(cfg_path), app)
  end if

  call build_mesh_from_config(app, mesh)
  call init_particles_from_config(app, pcls)
  call run_absorption_insulator(mesh, app%sim, pcls, stats)

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

  subroutine write_result_files(out_dir, mesh, stats)
    character(len=*), intent(in) :: out_dir
    type(mesh_type), intent(in) :: mesh
    type(sim_stats), intent(in) :: stats
    character(len=1024) :: cmd, summary_path, charges_path
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
  end subroutine write_result_files

end program main
