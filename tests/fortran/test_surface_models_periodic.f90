!> 周期セル内の導体再配分を cached_kneq0 と独立した Fourier reference で検証する。
program test_surface_models_periodic
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_surface_models, only: floating_conductor_solver_type
  use bem_types, only: mesh_type, sim_config, surface_model_conductor, surface_model_insulator, bc_open, bc_periodic
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, &
                          delete_file_if_exists, remove_empty_directory
  implicit none

  call test_init(2)

  call test_begin('periodic2_cached_conductor_equalizes_and_conserves_charge')
  call test_periodic2_cached_conductor_equalizes_and_conserves_charge()
  call test_end()

  call test_begin('periodic2_reference_conductor_is_translation_invariant')
  call test_periodic2_reference_conductor_is_translation_invariant()
  call test_end()

  call test_summary()

contains

  !> 2 軸周期セル内の導体板 (製品構成の cached_kneq0): 等電位と総電荷保存。
  subroutine test_periodic2_cached_conductor_equalizes_and_conserves_charge()
    type(mesh_type) :: mesh
    type(sim_config) :: sim
    type(field_physics_config) :: field_config
    type(periodic2_physics_config) :: periodic_config
    type(panel_kernel_config) :: panel_config
    type(electrostatic_snapshot_type) :: snapshot
    type(floating_conductor_solver_type) :: solver
    real(dp) :: v0(3, 3), v1(3, 3), v2(3, 3), q0(3), phi1, phi2, phi_scale
    integer(i32) :: mesh_ids(3), surface_models(3)
    character(len=512) :: cache_dir, cache_path
    integer(i64) :: clock_count

    call system_clock(count=clock_count)
    write (cache_dir, '(a,i0)') 'test_surface_models_tmp_', clock_count
    call periodic_plate_fixture(v0, v1, v2, q0, mesh_ids, surface_models)
    call init_periodic_mesh(mesh, v0, v1, v2, q0, mesh_ids, surface_models)
    call configure_periodic_cell(sim, field_config, periodic_config, panel_config)
    sim%field_periodic_cache_dir = cache_dir
    call snapshot%init(mesh, sim, field_config, periodic_config, panel_config)
    cache_path = snapshot%diagnostics%periodic_cache_path

    call solver%init(mesh, sim, snapshot)
    call assert_true(solver%active(), 'periodic conductor mesh must activate the solver')
    call solver%relax(mesh, sim, snapshot)

    call snapshot%refresh(mesh)
    call snapshot%eval_local_phi(mesh, sim, mesh%centers(:, 1), phi1)
    call snapshot%eval_local_phi(mesh, sim, mesh%centers(:, 2), phi2)
    phi_scale = max(abs(phi1), abs(phi2), 1.0e-30_dp)
    call assert_close_dp(mesh%q_elem(1) + mesh%q_elem(2), 1.0e-12_dp, 1.0e-24_dp, 'periodic conductor total charge mismatch')
    call assert_close_dp(mesh%q_elem(3), q0(3), 0.0_dp, 'periodic insulator charge must remain unchanged')
    call assert_close_dp(phi1/phi_scale, phi2/phi_scale, 1.0e-8_dp, 'periodic conductor elements should be equipotential')
    call assert_true(mesh%q_elem(1) /= mesh%q_elem(2), 'asymmetric insulator charge must redistribute the plate charge')

    call delete_file_if_exists(cache_path)
    call delete_file_if_exists(trim(cache_path)//'.lock')
    call remove_empty_directory(cache_dir)
  end subroutine test_periodic2_cached_conductor_equalizes_and_conserves_charge

  !> 厳密な panel_spectral_reference では、幾何をセル内で平行移動しても再配分結果は変わらない。
  !! FMM (cached_kneq0) は theta 近似のため平行移動で 1e-3 程度ずれるので、この比較は reference で行う。
  subroutine test_periodic2_reference_conductor_is_translation_invariant()
    type(mesh_type) :: mesh, shifted_mesh
    type(sim_config) :: sim
    type(field_physics_config) :: field_config
    type(periodic2_physics_config) :: periodic_config
    type(panel_kernel_config) :: panel_config
    type(electrostatic_snapshot_type) :: snapshot, shifted_snapshot
    type(floating_conductor_solver_type) :: solver, shifted_solver
    real(dp) :: v0(3, 3), v1(3, 3), v2(3, 3), q0(3), shift(3)
    integer(i32) :: mesh_ids(3), surface_models(3)

    call periodic_plate_fixture(v0, v1, v2, q0, mesh_ids, surface_models)
    call init_periodic_mesh(mesh, v0, v1, v2, q0, mesh_ids, surface_models)
    call configure_periodic_cell(sim, field_config, periodic_config, panel_config)
    sim%field_solver = 'direct'
    sim%field_periodic_far_correction = 'none'
    field_config = field_physics_config(backend='direct', normalization='si')
    periodic_config = periodic2_physics_config( &
                      nonzero_mode_backend='panel_spectral_reference', zero_mode_policy='exclude_k0', &
                      lower_boundary_model='symmetric_vacuum', reference_mode_layers=6_i32, panel_quadrature_order=12_i32 &
                      )
    panel_config = panel_kernel_config(kernel_id='triangle_p0_exact_direct', surface_side_policy='per_element')
    call snapshot%init(mesh, sim, field_config, periodic_config, panel_config)
    call solver%init(mesh, sim, snapshot)
    call solver%relax(mesh, sim, snapshot)
    call assert_true(mesh%q_elem(1) /= mesh%q_elem(2), 'reference case must redistribute the plate charge')

    shift = [0.30_dp, 0.0_dp, 0.0_dp]
    call init_periodic_mesh( &
      shifted_mesh, v0 + spread(shift, 2, 3), v1 + spread(shift, 2, 3), v2 + spread(shift, 2, 3), &
      q0, mesh_ids, surface_models &
      )
    call shifted_snapshot%init(shifted_mesh, sim, field_config, periodic_config, panel_config)
    call shifted_solver%init(shifted_mesh, sim, shifted_snapshot)
    call shifted_solver%relax(shifted_mesh, sim, shifted_snapshot)
    call assert_close_dp( &
      shifted_mesh%q_elem(1), mesh%q_elem(1), 1.0e-6_dp*abs(mesh%q_elem(1)), &
      'shifted periodic conductor charge 1 mismatch' &
      )
    call assert_close_dp( &
      shifted_mesh%q_elem(2), mesh%q_elem(2), 1.0e-6_dp*abs(mesh%q_elem(2)), &
      'shifted periodic conductor charge 2 mismatch' &
      )
  end subroutine test_periodic2_reference_conductor_is_translation_invariant

  !> conductor 板 (2 三角形, mesh_id 1) を z=0.25 に、帯電した insulator 三角形 (mesh_id 2) を z=0.65 に置く。
  subroutine periodic_plate_fixture(v0, v1, v2, q0, mesh_ids, surface_models)
    real(dp), intent(out) :: v0(3, 3), v1(3, 3), v2(3, 3), q0(3)
    integer(i32), intent(out) :: mesh_ids(3), surface_models(3)

    v0(:, 1) = [0.20_dp, 0.20_dp, 0.25_dp]
    v1(:, 1) = [0.50_dp, 0.20_dp, 0.25_dp]
    v2(:, 1) = [0.50_dp, 0.50_dp, 0.25_dp]
    v0(:, 2) = [0.20_dp, 0.20_dp, 0.25_dp]
    v1(:, 2) = [0.50_dp, 0.50_dp, 0.25_dp]
    v2(:, 2) = [0.20_dp, 0.50_dp, 0.25_dp]
    v0(:, 3) = [0.10_dp, 0.30_dp, 0.65_dp]
    ! Keep the charged panel inside its fixed box-cell proxy envelope.
    v1(:, 3) = [0.30_dp, 0.30_dp, 0.65_dp]
    v2(:, 3) = [0.10_dp, 0.45_dp, 0.65_dp]
    q0 = [1.0e-12_dp, 0.0_dp, 2.0e-12_dp]
    mesh_ids = [1_i32, 1_i32, 2_i32]
    surface_models = [surface_model_conductor, surface_model_conductor, surface_model_insulator]
  end subroutine periodic_plate_fixture

  subroutine init_periodic_mesh(mesh, v0, v1, v2, q0, mesh_ids, surface_models)
    type(mesh_type), intent(out) :: mesh
    real(dp), intent(in) :: v0(:, :), v1(:, :), v2(:, :), q0(:)
    integer(i32), intent(in) :: mesh_ids(:), surface_models(:)

    call init_mesh(mesh, v0, v1, v2, q0=q0, elem_mesh_id0=mesh_ids, elem_surface_model0=surface_models)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals
  end subroutine init_periodic_mesh

  subroutine configure_periodic_cell(sim, field_config, periodic_config, panel_config)
    type(sim_config), intent(out) :: sim
    type(field_physics_config), intent(out) :: field_config
    type(periodic2_physics_config), intent(out) :: periodic_config
    type(panel_kernel_config), intent(out) :: panel_config

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_far_correction = 'cached_kneq0'
    sim%field_periodic_image_layers = 1_i32
    sim%field_periodic_ewald_layers = 3_i32
    sim%field_periodic_generation_tolerance = 1.0e-8_dp
    sim%field_normalization = 'si'
    sim%tree_theta = 0.5_dp
    sim%has_tree_theta = .true.
    sim%tree_leaf_max = 8_i32
    sim%has_tree_leaf_max = .true.
    sim%e0 = [0.2_dp, -0.1_dp, 0.3_dp]
    sim%use_box = .true.
    sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    field_config = field_physics_config(backend='fmm', normalization='si')
    periodic_config = periodic2_physics_config( &
                      nonzero_mode_backend='cached_kneq0', zero_mode_policy='exclude_k0', &
                      lower_boundary_model='symmetric_vacuum' &
                      )
    panel_config = panel_kernel_config(kernel_id='triangle_p0_exact_p2m_near', surface_side_policy='per_element')
  end subroutine configure_periodic_cell

end program test_surface_models_periodic
