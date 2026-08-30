!> Public field-solver adapter contract for triangle-P0 panel FMM.
program test_dynamics_panel_fmm
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config
  use bem_templates, only: make_plane
  use bem_field_solver, only: field_solver_type
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, panel_surface_side_ok
  use test_support, only: test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  type(mesh_type) :: mesh
  type(field_solver_type) :: direct_solver, fmm_solver
  type(sim_config) :: direct_sim, fmm_sim
  type(field_physics_config) :: direct_field, fmm_field
  type(periodic2_physics_config) :: periodic
  type(panel_kernel_config) :: direct_panel, fmm_panel
  real(dp), parameter :: fmm_relative_tolerance = 2.0e-3_dp
  real(dp) :: target(3), direct_e(3), fmm_e(3), direct_phi, fmm_phi, rel_e, rel_phi
  real(dp), allocatable :: direct_mesh_phi(:), fmm_mesh_phi(:)
  integer(i32) :: element, status
  character(len=128) :: message

  call test_begin('panel_fmm_adapter_matches_direct_and_refreshes_charge')
  call make_plane(mesh, size_x=1.6_dp, size_y=1.6_dp, nx=8_i32, ny=8_i32, center=[0.0_dp, 0.0_dp, 0.0_dp])
  do element = 1_i32, mesh%nelem
    mesh%q_elem(element) = 1.0e-12_dp*(1.0_dp + 0.2_dp*sin(real(element, dp)))
  end do
  call resolve_panel_surface_sides(mesh, 'normal_plus', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'panel side setup status')

  direct_sim = sim_config()
  direct_sim%field_solver = 'direct'
  direct_sim%field_bc_mode = 'free'
  direct_field = field_physics_config(backend='direct', normalization='si')
  direct_panel = panel_kernel_config( &
                 kernel_id='triangle_p0_exact_direct', surface_side_policy='per_element' &
                 )
  call direct_solver%init(mesh, direct_sim, direct_field, periodic, direct_panel)

  fmm_sim = direct_sim
  fmm_sim%field_solver = 'fmm'
  fmm_sim%use_box = .true.
  fmm_sim%box_min = [-1.2_dp, -1.2_dp, -0.5_dp]
  fmm_sim%box_max = [1.2_dp, 1.2_dp, 2.5_dp]
  fmm_sim%tree_theta = 0.45_dp
  fmm_sim%tree_leaf_max = 4_i32
  fmm_sim%has_tree_theta = .true.
  fmm_sim%has_tree_leaf_max = .true.
  fmm_field = field_physics_config(backend='fmm', normalization='si')
  fmm_panel = panel_kernel_config( &
              kernel_id='triangle_p0_exact_p2m_near', surface_side_policy='per_element' &
              )
  call fmm_solver%init(mesh, fmm_sim, fmm_field, periodic, fmm_panel)

  target = [0.17_dp, -0.21_dp, 0.42_dp]
  call direct_solver%eval_e(mesh, target, direct_e)
  call fmm_solver%eval_e(mesh, target, fmm_e)
  call direct_solver%eval_potential(mesh, direct_sim, target, direct_phi)
  call fmm_solver%eval_potential(mesh, fmm_sim, target, fmm_phi)
  rel_e = sqrt(sum((fmm_e - direct_e)**2))/max(sqrt(sum(direct_e*direct_e)), tiny(1.0_dp))
  rel_phi = abs(fmm_phi - direct_phi)/max(abs(direct_phi), tiny(1.0_dp))
  call assert_true(rel_e < fmm_relative_tolerance, 'panel FMM adapter field error exceeds tolerance')
  call assert_true(rel_phi < fmm_relative_tolerance, 'panel FMM adapter potential error exceeds tolerance')
  allocate (direct_mesh_phi(mesh%nelem), fmm_mesh_phi(mesh%nelem))
  call direct_solver%compute_mesh_potential(mesh, direct_sim, direct_mesh_phi)
  call fmm_solver%compute_mesh_potential(mesh, fmm_sim, fmm_mesh_phi)
  call assert_true( &
    maxval(abs(fmm_mesh_phi - direct_mesh_phi))/max(maxval(abs(direct_mesh_phi)), tiny(1.0_dp)) < &
    fmm_relative_tolerance, &
    'panel FMM mesh potential must not add a second self term' &
    )
  deallocate (direct_mesh_phi, fmm_mesh_phi)

  do element = 1_i32, mesh%nelem
    mesh%q_elem(element) = 1.0e-12_dp*(0.6_dp + 0.3_dp*cos(0.5_dp*real(element, dp)))
  end do
  call fmm_solver%refresh(mesh)
  call direct_solver%eval_e(mesh, target, direct_e)
  call fmm_solver%eval_e(mesh, target, fmm_e)
  call direct_solver%eval_potential(mesh, direct_sim, target, direct_phi)
  call fmm_solver%eval_potential(mesh, fmm_sim, target, fmm_phi)
  rel_e = sqrt(sum((fmm_e - direct_e)**2))/max(sqrt(sum(direct_e*direct_e)), tiny(1.0_dp))
  rel_phi = abs(fmm_phi - direct_phi)/max(abs(direct_phi), tiny(1.0_dp))
  call assert_true(rel_e < fmm_relative_tolerance, 'refreshed panel FMM field error exceeds tolerance')
  call assert_true(rel_phi < fmm_relative_tolerance, 'refreshed panel FMM potential error exceeds tolerance')
  call test_end()

  call test_begin('panel_auto_selects_fmm_adapter')
  fmm_sim = direct_sim
  fmm_sim%field_solver = 'auto'
  fmm_sim%tree_min_nelem = 64_i32
  fmm_field = field_physics_config(backend='auto', normalization='si')
  fmm_panel = panel_kernel_config( &
              kernel_id='triangle_p0_exact_auto', surface_side_policy='per_element' &
              )
  call fmm_solver%init(mesh, fmm_sim, fmm_field, periodic, fmm_panel)
  call assert_true(trim(fmm_solver%mode) == 'fmm', 'large panel auto case must select FMM')
  call test_end()
  call test_summary()
end program test_dynamics_panel_fmm
