program main
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, sim_stats, mesh_type, particles_soa
  use bem_templates, only: make_plane
  use bem_importers, only: load_obj_mesh
  use bem_injection, only: seed_rng, init_random_beam_particles
  use bem_simulator, only: run_absorption_insulator
  implicit none

  type(mesh_type) :: mesh
  type(particles_soa) :: pcls
  type(sim_config) :: cfg
  type(sim_stats) :: stats

  logical :: has_obj
  inquire(file='examples/simple_plate.obj', exist=has_obj)
  if (has_obj) then
    call load_obj_mesh('examples/simple_plate.obj', mesh)
  else
    call make_plane(mesh, size_x=1.0d0, size_y=1.0d0, nx=1, ny=1, center=[0.0d0, 0.0d0, 0.0d0])
  end if

  call seed_rng([12345_i32])
  call init_random_beam_particles(
    pcls=pcls,
    n=256_i32,
    q_particle=-1.602176634d-19,
    m_particle=9.10938356d-31,
    w_particle=1.0d0,
    pos_low=[-0.4d0, -0.4d0, 0.2d0],
    pos_high=[0.4d0, 0.4d0, 0.5d0],
    drift_velocity=[0.0d0, 0.0d0, -8.0d5],
    temperature_k=2.0d4
  )

  cfg%dt = 1.0d-9
  cfg%npcls_per_step = 64
  cfg%max_step = 400
  cfg%tol_rel = 1.0d-8
  cfg%softening = 1.0d-6
  cfg%b0 = [0.0d0, 0.0d0, 0.0d0]

  call run_absorption_insulator(mesh, cfg, pcls, stats)

  print '(a,i0)', 'mesh nelem=', mesh%nelem
  print '(a,i0)', 'processed_particles=', stats%processed_particles
  print '(a,i0)', 'absorbed=', stats%absorbed
  print '(a,i0)', 'escaped=', stats%escaped
  print '(a,i0)', 'batches=', stats%batches
  print '(a,es12.4)', 'last_rel_change=', stats%last_rel_change
  print '(a,*(es12.4,1x))', 'mesh charges=', mesh%q_elem
end program main
