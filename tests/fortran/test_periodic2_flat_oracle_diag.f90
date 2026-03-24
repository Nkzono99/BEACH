!> periodic2 flat-plane に対する oracle far correction の診断。
program test_periodic2_flat_oracle_diag
  use iso_fortran_env, only: output_unit
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_templates, only: make_plane
  use bem_field_solver, only: field_solver_type
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_correction_all_sources
  use test_support, only: test_init, test_begin, test_end, test_summary
  implicit none

  type(mesh_type) :: mesh_fmm
  type(field_solver_type) :: solver_none = field_solver_type()
  type(field_solver_type) :: solver_oracle = field_solver_type()
  type(sim_config) :: sim
  integer(i32) :: i, max_none_idx, max_oracle_idx, max_corr_idx
  real(dp) :: r(3), e_none(3), e_oracle(3), e_none_ref(3), e_oracle_ref(3)
  real(dp) :: corr_fmm(3), corr_ref(3)
  real(dp) :: points(3, 8)
  real(dp) :: ref_norm, none_ref_norm, diff_norm
  real(dp) :: max_none_rel_err, max_oracle_rel_err, max_corr_rel_err, max_corr_vs_shell

  call test_init(1)

  call test_begin('periodic2_flat_oracle_diag')
  call make_plane(mesh_fmm, size_x=1.0d0, size_y=1.0d0, nx=20_i32, ny=20_i32, center=[0.5d0, 0.5d0, 0.02d0])
  call assign_flat_periodic_test_charges(mesh_fmm)

  sim = sim_config()
  sim%softening = 1.0d-6
  sim%field_solver = 'fmm'
  sim%field_bc_mode = 'periodic2'
  sim%field_periodic_far_correction = 'none'
  sim%field_periodic_image_layers = 1_i32
  sim%field_periodic_ewald_layers = 4_i32
  sim%tree_min_nelem = 64_i32
  sim%use_box = .true.
  sim%box_min = [0.0d0, 0.0d0, 0.0d0]
  sim%box_max = [1.0d0, 1.0d0, 10.0d0]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  call solver_none%init(mesh_fmm, sim)
  call solver_none%refresh(mesh_fmm)

  sim%field_periodic_far_correction = 'm2l_root_oracle'
  call solver_oracle%init(mesh_fmm, sim)
  call solver_oracle%refresh(mesh_fmm)

  points(:, 1) = [0.50d0, 0.50d0, 0.03d0]
  points(:, 2) = [0.50d0, 0.50d0, 0.10d0]
  points(:, 3) = [0.85d0, 0.50d0, 0.03d0]
  points(:, 4) = [0.95d0, 0.50d0, 0.03d0]
  points(:, 5) = [0.95d0, 0.95d0, 0.03d0]
  points(:, 6) = [0.95d0, 0.95d0, 0.10d0]
  points(:, 7) = [0.20d0, 0.80d0, 0.03d0]
  points(:, 8) = [0.80d0, 0.20d0, 0.10d0]

  max_none_rel_err = 0.0d0
  max_oracle_rel_err = 0.0d0
  max_corr_rel_err = 0.0d0
  max_corr_vs_shell = 0.0d0
  max_none_idx = 0_i32
  max_oracle_idx = 0_i32
  max_corr_idx = 0_i32
  do i = 1_i32, int(size(points, 2), i32)
    r = points(:, i)
    call solver_none%eval_e(mesh_fmm, r, e_none)
    call solver_oracle%eval_e(mesh_fmm, r, e_oracle)
    call electric_field_at_periodic2_reference(mesh_fmm, solver_none, r, e_none_ref)
    call electric_field_at_periodic2_reference(mesh_fmm, solver_oracle, r, e_oracle_ref)

    ref_norm = max(sqrt(sum(e_none_ref*e_none_ref)), 1.0d-30)
    diff_norm = sqrt(sum((e_none - e_none_ref)*(e_none - e_none_ref)))
    if (diff_norm/ref_norm > max_none_rel_err) then
      max_none_rel_err = diff_norm/ref_norm
      max_none_idx = i
    end if

    ref_norm = max(sqrt(sum(e_oracle_ref*e_oracle_ref)), 1.0d-30)
    diff_norm = sqrt(sum((e_oracle - e_oracle_ref)*(e_oracle - e_oracle_ref)))
    if (diff_norm/ref_norm > max_oracle_rel_err) then
      max_oracle_rel_err = diff_norm/ref_norm
      max_oracle_idx = i
    end if

    corr_fmm = e_oracle - e_none
    corr_ref = e_oracle_ref - e_none_ref
    ref_norm = max(sqrt(sum(corr_ref*corr_ref)), 1.0d-30)
    diff_norm = sqrt(sum((corr_fmm - corr_ref)*(corr_fmm - corr_ref)))
    if (diff_norm/ref_norm > max_corr_rel_err) then
      max_corr_rel_err = diff_norm/ref_norm
      max_corr_idx = i
    end if

    none_ref_norm = max(sqrt(sum(e_none_ref*e_none_ref)), 1.0d-30)
    max_corr_vs_shell = max(max_corr_vs_shell, diff_norm/none_ref_norm)

    write (output_unit, '(a,i0,a,3(1x,es12.4))') 'point_', i, '_xyz=', r
    write (output_unit, '(a,3(1x,es12.4))') '  none_ref   =', e_none_ref
    write (output_unit, '(a,3(1x,es12.4))') '  none_fmm   =', e_none
    write (output_unit, '(a,3(1x,es12.4))') '  oracle_ref =', e_oracle_ref
    write (output_unit, '(a,3(1x,es12.4))') '  oracle_fmm =', e_oracle
    write (output_unit, '(a,3(1x,es12.4))') '  corr_ref   =', corr_ref
    write (output_unit, '(a,3(1x,es12.4))') '  corr_fmm   =', corr_fmm
  end do

  write (output_unit, '(a,es12.4,a,i0)') 'diag_flat_none_rel_err=', max_none_rel_err, ' idx=', max_none_idx
  write (output_unit, '(a,es12.4,a,i0)') 'diag_flat_oracle_rel_err=', max_oracle_rel_err, ' idx=', max_oracle_idx
  write (output_unit, '(a,es12.4,a,i0)') 'diag_flat_corr_rel_err=', max_corr_rel_err, ' idx=', max_corr_idx
  write (output_unit, '(a,es12.4)') 'diag_flat_corr_vs_shell=', max_corr_vs_shell
  call test_end()

  call test_summary()

contains

  subroutine assign_flat_periodic_test_charges(mesh)
    type(mesh_type), intent(inout) :: mesh
    integer(i32) :: elem_idx
    real(dp), parameter :: two_pi = 2.0d0*acos(-1.0d0)
    real(dp) :: x, y

    do elem_idx = 1_i32, mesh%nelem
      x = mesh%center_x(elem_idx)
      y = mesh%center_y(elem_idx)
      mesh%q_elem(elem_idx) = -1.0d-13*( &
                              1.0d0 &
                              + 0.35d0*cos(two_pi*x) &
                              - 0.25d0*sin(two_pi*y) &
                              + 0.20d0*cos(two_pi*(x + y)) &
                              )
    end do
  end subroutine assign_flat_periodic_test_charges

  subroutine electric_field_at_periodic2_images(mesh, r, softening, box_min, box_max, periodic_axes, nimg, e)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: softening
    real(dp), intent(in) :: box_min(3), box_max(3)
    integer(i32), intent(in) :: periodic_axes(2)
    integer(i32), intent(in) :: nimg
    real(dp), intent(out) :: e(3)

    integer(i32) :: elem_idx, ix, iy
    integer(i32) :: axis1, axis2
    real(dp) :: l1, l2, soft2, shift1, shift2
    real(dp) :: src(3), d(3), r2, inv_r3
    real(dp) :: e_core(3)

    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    l1 = box_max(axis1) - box_min(axis1)
    l2 = box_max(axis2) - box_min(axis2)
    soft2 = softening*softening
    e_core = 0.0d0

    do elem_idx = 1_i32, mesh%nelem
      do ix = -nimg, nimg
        do iy = -nimg, nimg
          src = [mesh%center_x(elem_idx), mesh%center_y(elem_idx), mesh%center_z(elem_idx)]
          shift1 = real(ix, dp)*l1
          shift2 = real(iy, dp)*l2
          src(axis1) = src(axis1) + shift1
          src(axis2) = src(axis2) + shift2
          d = r - src
          r2 = sum(d*d) + soft2
          if (r2 <= tiny(1.0d0)) cycle
          inv_r3 = 1.0d0/(sqrt(r2)*r2)
          e_core = e_core + mesh%q_elem(elem_idx)*inv_r3*d
        end do
      end do
    end do

    e = k_coulomb*e_core
  end subroutine electric_field_at_periodic2_images

  subroutine electric_field_at_periodic2_reference(mesh, solver, r, e)
    type(mesh_type), intent(in) :: mesh
    type(field_solver_type), intent(in) :: solver
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: e(3)
    real(dp) :: wrapped_r(3), e_core(3)
    integer(i32) :: axis

    wrapped_r = r
    do axis = 1_i32, 2_i32
      wrapped_r(solver%periodic_axes(axis)) = solver%target_box_min(solver%periodic_axes(axis)) + &
                                              modulo( &
                                              wrapped_r(solver%periodic_axes(axis)) - &
                                              solver%target_box_min(solver%periodic_axes(axis)), &
                                              solver%periodic_len(axis) &
                                              )
    end do

    call electric_field_at_periodic2_images( &
      mesh, wrapped_r, solver%softening, solver%target_box_min, solver%target_box_max, solver%periodic_axes, &
      solver%periodic_image_layers, e &
      )
    if (trim(solver%periodic_far_correction) /= 'm2l_root_oracle') return

    e_core = e/k_coulomb
    call add_periodic2_exact_ewald_correction_all_sources(solver%fmm_core_plan, solver%fmm_core_state, wrapped_r, e_core)
    e = k_coulomb*e_core
  end subroutine electric_field_at_periodic2_reference

end program test_periodic2_flat_oracle_diag
