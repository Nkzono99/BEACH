program main
  use bem_kinds, only: dp
  use bem_types, only: sim_config, sim_stats, mesh_type, particles_soa
  use bem_templates, only: make_plane
  use bem_importers, only: load_obj_mesh
  use bem_particles, only: init_particles
  use bem_simulator, only: run_absorption_insulator
  implicit none

  type(mesh_type) :: mesh
  type(particles_soa) :: pcls
  type(sim_config) :: cfg
  type(sim_stats) :: stats

  real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:)

  logical :: has_obj
  inquire(file='examples/simple_plate.obj', exist=has_obj)
  if (has_obj) then
    call load_obj_mesh('examples/simple_plate.obj', mesh)
  else
    call make_plane(mesh, size_x=1.0d0, size_y=1.0d0, nx=1, ny=1, center=[0.0d0, 0.0d0, 0.0d0])
  end if

  allocate(x(3,4), v(3,4), q(4), m(4))

  x(:,1) = [0.0d0, 0.0d0, 0.2d0]
  x(:,2) = [0.2d0, 0.1d0, 0.2d0]
  x(:,3) = [-0.2d0, -0.1d0, 0.2d0]
  x(:,4) = [0.7d0, 0.0d0, 0.2d0]
  v(:,1:4) = reshape([0.0d0,0.0d0,-1.0d6,  0.0d0,0.0d0,-1.0d6,  0.0d0,0.0d0,-1.0d6,  0.0d0,0.0d0,-1.0d6], [3,4])
  q = -1.602176634d-19
  m = 9.10938356d-31

  call init_particles(pcls, x, v, q, m)

  cfg%dt = 1.0d-9
  cfg%npcls_per_step = 2
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
