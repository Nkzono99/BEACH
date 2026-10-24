!> cached_kneq0 snapshot test と adaptive batch test が共有する periodic2 fixture。
module test_periodic2_cached_fixture
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  implicit none
  private

  public :: configure_cached_periodic2_fixture

contains

  !> 2 枚の panel を置いた 1 m 立方の periodic2 cell。cached_kneq0 FMM、symmetric_vacuum、一様外部場つき。
  !! cache_dir は実行ごとに一意な directory を渡す。
  subroutine configure_cached_periodic2_fixture(mesh_out, sim_out, field_out, periodic_out, panel_out, a, b, c, cache_dir)
    type(mesh_type), intent(out) :: mesh_out
    type(sim_config), intent(out) :: sim_out
    type(field_physics_config), intent(out) :: field_out
    type(periodic2_physics_config), intent(out) :: periodic_out
    type(panel_kernel_config), intent(out) :: panel_out
    real(dp), intent(out) :: a(3, 2), b(3, 2), c(3, 2)
    character(len=*), intent(in) :: cache_dir

    a(:, 1) = [0.10_dp, 0.10_dp, 0.25_dp]
    b(:, 1) = [0.55_dp, 0.10_dp, 0.25_dp]
    c(:, 1) = [0.10_dp, 0.55_dp, 0.25_dp]
    a(:, 2) = [0.45_dp, 0.45_dp, 0.65_dp]
    b(:, 2) = [0.90_dp, 0.45_dp, 0.65_dp]
    c(:, 2) = [0.90_dp, 0.90_dp, 0.65_dp]
    ! Resolve the two-panel fixture inside the box-cell proxy envelope.
    b = a + 0.45_dp*(b - a)
    c = a + 0.45_dp*(c - a)
    call init_mesh(mesh_out, a, b, c, q0=[1.0e-12_dp, 2.0e-12_dp])
    mesh_out%elem_vacuum_sign = 1_i32
    mesh_out%vacuum_normals = mesh_out%normals

    sim_out = sim_config()
    sim_out%field_solver = 'fmm'
    sim_out%field_bc_mode = 'periodic2'
    sim_out%field_periodic_far_correction = 'cached_kneq0'
    sim_out%field_periodic_image_layers = 1_i32
    sim_out%field_periodic_ewald_layers = 3_i32
    sim_out%field_periodic_cache_dir = cache_dir
    sim_out%field_periodic_generation_tolerance = 1.0e-8_dp
    sim_out%multiple_box_events_retry_backend = 'upper_panel_fourier'
    sim_out%field_normalization = 'si'
    sim_out%tree_theta = 0.5_dp
    sim_out%has_tree_theta = .true.
    sim_out%tree_leaf_max = 8_i32
    sim_out%has_tree_leaf_max = .true.
    sim_out%e0 = [0.2_dp, -0.1_dp, 0.3_dp]
    sim_out%use_box = .true.
    sim_out%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    sim_out%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    sim_out%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim_out%bc_high = [bc_periodic, bc_periodic, bc_open]
    field_out = field_physics_config(backend='fmm', normalization='si')
    periodic_out = periodic2_physics_config( &
                   nonzero_mode_backend='cached_kneq0', zero_mode_policy='exclude_k0', &
                   lower_boundary_model='symmetric_vacuum' &
                   )
    panel_out = panel_kernel_config( &
                kernel_id='triangle_p0_exact_p2m_near', surface_side_policy='per_element' &
                )
  end subroutine configure_cached_periodic2_fixture

end module test_periodic2_cached_fixture
