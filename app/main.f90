program main
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
end program main
