!> surface_model に応じた電荷緩和を検証するテスト。
program test_surface_models
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_mesh, only: init_mesh
  use bem_surface_models, only: apply_surface_model_charge_relaxation
  use bem_types, only: mesh_type, surface_model_conductor, surface_model_insulator
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  call test_init(4)

  call test_begin('isolated_conductor_equalizes_symmetric_charge')
  call test_isolated_conductor_equalizes_symmetric_charge()
  call test_end()

  call test_begin('fixed_charge_induces_conductor_dipole')
  call test_fixed_charge_induces_conductor_dipole()
  call test_end()

  call test_begin('uniform_field_induces_conductor_dipole')
  call test_uniform_field_induces_conductor_dipole()
  call test_end()

  call test_begin('multiple_mesh_id_conductors_conserve_charge_independently')
  call test_multiple_mesh_id_conductors_conserve_charge_independently()
  call test_end()

  call test_summary()

contains

  subroutine test_isolated_conductor_equalizes_symmetric_charge()
    type(mesh_type) :: mesh
    real(dp) :: centers(3, 2), q0(2)
    integer(i32) :: mesh_ids(2), surface_models(2)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [1.0d0, 0.0d0, 0.0d0]
    q0 = [2.0d-12, 0.0d0]
    mesh_ids = [1_i32, 1_i32]
    surface_models = [surface_model_conductor, surface_model_conductor]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)

    call apply_surface_model_charge_relaxation(mesh, 0.2d0, [0.0d0, 0.0d0, 0.0d0])

    call assert_close_dp(sum(mesh%q_elem), 2.0d-12, 1.0d-24, 'conductor total charge must be conserved')
    call assert_close_dp(mesh%q_elem(1), 1.0d-12, 1.0d-24, 'symmetric conductor charge q1 mismatch')
    call assert_close_dp(mesh%q_elem(2), 1.0d-12, 1.0d-24, 'symmetric conductor charge q2 mismatch')
  end subroutine test_isolated_conductor_equalizes_symmetric_charge

  subroutine test_fixed_charge_induces_conductor_dipole()
    type(mesh_type) :: mesh
    real(dp) :: centers(3, 3), q0(3), phi1, phi2
    integer(i32) :: mesh_ids(3), surface_models(3)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    centers(:, 3) = [-1.0d0, 0.0d0, 0.0d0]
    q0 = [0.0d0, 0.0d0, 1.0d-10]
    mesh_ids = [1_i32, 1_i32, 2_i32]
    surface_models = [surface_model_conductor, surface_model_conductor, surface_model_insulator]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)

    call apply_surface_model_charge_relaxation(mesh, 0.2d0, [0.0d0, 0.0d0, 0.0d0])

    phi1 = scaled_potential_at_elem(mesh, 1_i32, 0.2d0, [0.0d0, 0.0d0, 0.0d0])
    phi2 = scaled_potential_at_elem(mesh, 2_i32, 0.2d0, [0.0d0, 0.0d0, 0.0d0])
    call assert_close_dp(mesh%q_elem(1) + mesh%q_elem(2), 0.0d0, 1.0d-22, 'induced conductor net charge mismatch')
    call assert_close_dp(phi1, phi2, 1.0d-20, 'conductor elements should be equipotential')
    call assert_true(mesh%q_elem(1) < 0.0d0, 'near fixed positive charge should induce negative conductor charge')
    call assert_true(mesh%q_elem(2) > 0.0d0, 'far conductor side should induce positive charge')
  end subroutine test_fixed_charge_induces_conductor_dipole

  subroutine test_uniform_field_induces_conductor_dipole()
    type(mesh_type) :: mesh
    real(dp) :: centers(3, 2), q0(2), external_e(3), phi1, phi2
    integer(i32) :: mesh_ids(2), surface_models(2)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    q0 = [0.0d0, 0.0d0]
    mesh_ids = [1_i32, 1_i32]
    surface_models = [surface_model_conductor, surface_model_conductor]
    external_e = [1.0d0, 0.0d0, 0.0d0]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)

    call apply_surface_model_charge_relaxation(mesh, 0.2d0, external_e)

    phi1 = scaled_potential_at_elem(mesh, 1_i32, 0.2d0, external_e)
    phi2 = scaled_potential_at_elem(mesh, 2_i32, 0.2d0, external_e)
    call assert_close_dp(mesh%q_elem(1) + mesh%q_elem(2), 0.0d0, 1.0d-24, 'uniform-field net charge mismatch')
    call assert_close_dp(phi1, phi2, 1.0d-20, 'uniform-field conductor should be equipotential')
    call assert_true(mesh%q_elem(1) < 0.0d0, 'left side should become negative for +x external field')
    call assert_true(mesh%q_elem(2) > 0.0d0, 'right side should become positive for +x external field')
  end subroutine test_uniform_field_induces_conductor_dipole

  subroutine test_multiple_mesh_id_conductors_conserve_charge_independently()
    type(mesh_type) :: mesh
    real(dp) :: centers(3, 4), q0(4)
    integer(i32) :: mesh_ids(4), surface_models(4)

    centers(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    centers(:, 2) = [1.0d0, 0.0d0, 0.0d0]
    centers(:, 3) = [4.0d0, 0.0d0, 0.0d0]
    centers(:, 4) = [5.0d0, 0.0d0, 0.0d0]
    q0 = [2.0d-12, 0.0d0, 4.0d-12, 0.0d0]
    mesh_ids = [1_i32, 1_i32, 2_i32, 2_i32]
    surface_models = [surface_model_conductor, surface_model_conductor, &
                      surface_model_conductor, surface_model_conductor]
    call init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)

    call apply_surface_model_charge_relaxation(mesh, 0.2d0, [0.0d0, 0.0d0, 0.0d0])

    call assert_close_dp(sum(mesh%q_elem(1:2)), 2.0d-12, 1.0d-23, &
                         'first conductor mesh_id total charge mismatch')
    call assert_close_dp(sum(mesh%q_elem(3:4)), 4.0d-12, 1.0d-23, &
                         'second conductor mesh_id total charge mismatch')
  end subroutine test_multiple_mesh_id_conductors_conserve_charge_independently

  subroutine init_test_mesh(mesh, centers, q0, mesh_ids, surface_models)
    type(mesh_type), intent(out) :: mesh
    real(dp), intent(in) :: centers(:, :)
    real(dp), intent(in) :: q0(:)
    integer(i32), intent(in) :: mesh_ids(:), surface_models(:)
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), parameter :: half_size = 0.05d0
    integer(i32) :: i, n

    n = int(size(centers, 2), kind=i32)
    allocate (v0(3, n), v1(3, n), v2(3, n))
    do i = 1, n
      v0(:, i) = centers(:, i) + [-half_size, -half_size, 0.0d0]
      v1(:, i) = centers(:, i) + [half_size, 0.0d0, 0.0d0]
      v2(:, i) = centers(:, i) + [0.0d0, half_size, 0.0d0]
    end do
    call init_mesh(mesh, v0, v1, v2, q0=q0, elem_mesh_id0=mesh_ids, elem_surface_model0=surface_models)
  end subroutine init_test_mesh

  real(dp) function scaled_potential_at_elem(mesh, elem_i, softening, external_e) result(phi)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elem_i
    real(dp), intent(in) :: softening
    real(dp), intent(in) :: external_e(3)
    integer(i32) :: elem_j

    phi = -dot_product(external_e, mesh%centers(:, elem_i))/k_coulomb
    do elem_j = 1, mesh%nelem
      phi = phi + mesh%q_elem(elem_j)*potential_coeff(mesh, elem_i, elem_j, softening)
    end do
  end function scaled_potential_at_elem

  real(dp) function potential_coeff(mesh, elem_i, elem_j, softening) result(coeff)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elem_i, elem_j
    real(dp), intent(in) :: softening
    real(dp), parameter :: pi_dp = acos(-1.0d0)
    real(dp) :: dx, dy, dz, r2

    if (elem_i == elem_j) then
      coeff = 1.0d0/softening
      if (softening <= 0.0d0) coeff = 2.0d0*sqrt(pi_dp)/max(mesh%h_elem(elem_i), sqrt(tiny(1.0d0)))
      return
    end if

    dx = mesh%center_x(elem_i) - mesh%center_x(elem_j)
    dy = mesh%center_y(elem_i) - mesh%center_y(elem_j)
    dz = mesh%center_z(elem_i) - mesh%center_z(elem_j)
    r2 = dx*dx + dy*dy + dz*dz + softening*softening
    coeff = 1.0d0/sqrt(max(r2, tiny(1.0d0)))
  end function potential_coeff

end program test_surface_models
