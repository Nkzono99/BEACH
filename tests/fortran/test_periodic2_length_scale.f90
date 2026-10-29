!> periodic2 の cached_kneq0 が µm 単位の形状で長さの尺度に従い、k≠0 成分に面平均を漏らさないことを確かめる。
program test_periodic2_length_scale
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb, eps0
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, &
                                  build_plan, build_panel_plan, update_state, eval_points, eval_potential_points, &
                                  destroy_plan, destroy_state
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_allclose_1d, assert_true, &
                          unique_test_path, delete_periodic_cache, remove_empty_directory
  implicit none
  type(fmm_options_type) :: options
  type(fmm_plan_type) :: plan
  type(fmm_state_type) :: state
  real(dp) :: source(3, 4), charge(4), target(3, 3), field(3, 3), potential(3)
  real(dp) :: reference_field(3, 3), reference_potential(3), length, tolerance
  character(len=512) :: cache_path
  character(len=:), allocatable :: cache_dir
  integer :: case_idx, point

  cache_dir = unique_test_path('test_periodic2_length_scale_tmp_')
  source(:, 1) = [0.15_dp, 0.20_dp, 0.10_dp]
  source(:, 2) = [0.75_dp, 0.25_dp, 0.80_dp]
  source(:, 3) = [0.30_dp, 0.80_dp, 1.20_dp]
  source(:, 4) = [0.82_dp, 0.72_dp, 1.60_dp]
  charge = [1.0_dp, -0.2_dp, 0.5_dp, 0.7_dp]
  target(:, 1) = [0.35_dp, 0.40_dp, 3.20_dp]
  target(:, 2) = [0.60_dp, 0.20_dp, 1.20_dp]
  target(:, 3) = [0.10_dp, 0.75_dp, 0.40_dp]
  options%theta = 0.45_dp
  options%leaf_max = 2_i32
  options%order = 4_i32
  options%use_periodic2 = .true.
  options%periodic_axes = [1_i32, 2_i32]
  options%periodic_image_layers = 1_i32
  options%periodic_far_correction = 'cached_kneq0'
  options%periodic_ewald_layers = 4_i32
  options%periodic_cache_dir = cache_dir
  options%periodic_generation_tolerance = 1.0e-8_dp
  options%target_box_min = 0.0_dp
  ! Scaling is independent of the operator partition; the root multipole keeps this case quick.
  options%periodic_operator_source_cells = 1_i32
  options%periodic_operator_anchor_cells = 2_i32
  call test_init(3)
  call test_begin('periodic Coulomb field obeys its physical length scaling')
  do case_idx = 1, 2
    length = 1.0_dp
    if (case_idx == 2) length = 1.0e-4_dp
    options%periodic_len = length
    options%target_box_max = length*[1.0_dp, 1.0_dp, 3.2_dp]
    call build_plan(plan, length*source, options)
    cache_path = plan%periodic_cache_path
    call update_state(plan, state, charge)
    call eval_points(plan, state, length*target, field)
    call eval_potential_points(plan, state, length*target, potential)
    field = field*length**2
    potential = potential*length
    if (case_idx == 1) then
      reference_field = field
      reference_potential = potential
    else
      tolerance = 1.0e-7_dp*max(1.0_dp, maxval(abs(reference_field)))
      do point = 1, 3
        write (*, '(a,i0,a,es12.5)') 'point=', point, ' scaled field error=', &
          maxval(abs(field(:, point) - reference_field(:, point)))
        call assert_allclose_1d(field(:, point), reference_field(:, point), tolerance, &
                                'field must scale as inverse length squared')
      end do
      write (*, '(a,es12.5)') 'scaled potential error=', maxval(abs(potential - reference_potential))
      tolerance = 1.0e-7_dp*max(1.0_dp, maxval(abs(reference_potential)))
      call assert_allclose_1d(potential, reference_potential, tolerance, 'potential must scale as inverse length')
    end if
    call destroy_state(state)
    call destroy_plan(plan)
    call delete_periodic_cache(cache_path)
  end do
  call test_end()

  call test_begin('uniform sheet leaves no plane mean in the nonzero modes')
  ! A uniform sheet has no k/=0 content, so the cached_kneq0 potential must average to zero on every plane.
  ! Uncharged points stretch the source tree over a particle-bed height, as in a regolith cell. A single
  ! root multipole then left a leak of about 1e-2; the analytic shell mean and the box-cell partition of
  ! the far correction remove it with one image layer.
  options%leaf_max = 8_i32
  options%periodic_operator_source_cells = -1_i32
  options%periodic_operator_anchor_cells = -1_i32
  call check_uniform_sheet_leak(2.0e-4_dp)
  call test_end()
  call check_panel_vacuum_reference()
  call remove_empty_directory(cache_dir)
  call test_summary()

contains

  subroutine check_panel_vacuum_reference()
    type(mesh_type) :: mesh
    type(fmm_options_type) :: config
    integer(i32), parameter :: side = 8_i32, n = side*side
    real(dp) :: a(3, 2), b(3, 2), c(3, 2), q(2), points(3, n), fields(3, n), potentials(n)
    real(dp) :: reference_phi, reference_e(3), errors(n), horizontal(n), mean_error, rms_error, field_error
    real(dp) :: z, scale
    integer(i32) :: i, j, k, idx

    call test_begin('box_panel_operator_independent_fourier_vacuum')
    config = options
    ! Isolate the far operator: use direct finite images instead of ordinary FMM truncation in this two-panel oracle.
    config%theta = 0.01_dp
    config%leaf_max = 1_i32
    config%periodic_operator_source_cells = -1_i32
    config%periodic_operator_anchor_cells = -1_i32
    config%periodic_len = 1.0_dp
    config%target_box_min = 0.0_dp
    config%target_box_max = 1.0_dp
    a(:, 1) = [0.19_dp, 0.21_dp, 0.15_dp]
    b(:, 1) = [0.29_dp, 0.21_dp, 0.15_dp]
    c(:, 1) = [0.19_dp, 0.31_dp, 0.18_dp]
    a(:, 2) = [0.61_dp, 0.63_dp, 0.27_dp]
    b(:, 2) = [0.71_dp, 0.63_dp, 0.30_dp]
    c(:, 2) = [0.61_dp, 0.73_dp, 0.30_dp]
    q = [1.0e-12_dp, -0.4e-12_dp]
    call init_mesh(mesh, a, b, c, q0=q)
    call build_panel_plan(plan, a, b, c, config)
    cache_path = plan%periodic_cache_path
    call update_state(plan, state, q)
    scale = sum(abs(q))/eps0
    do k = 1_i32, 3_i32
      z = 0.3_dp + 0.2_dp*real(k, dp)
      idx = 0_i32
      do j = 1_i32, side
        do i = 1_i32, side
          idx = idx + 1_i32
          points(:, idx) = [(real(i, dp) - 0.5_dp)/side, (real(j, dp) - 0.5_dp)/side, z]
        end do
      end do
      call eval_points(plan, state, points, fields)
      call eval_potential_points(plan, state, points, potentials)
      do idx = 1_i32, n
        call eval_periodic_nonzero_panel_reference(mesh, points(:, idx), 1.0_dp, 1.0_dp, 12_i32, 16_i32, &
                                                   reference_phi, reference_e)
        errors(idx) = k_coulomb*potentials(idx) - reference_phi
        horizontal(idx) = sum((k_coulomb*fields(1:2, idx) - reference_e(1:2))**2)
      end do
      mean_error = sum(errors)/real(n, dp)
      rms_error = sqrt(sum((errors - mean_error)**2)/real(n, dp))
      field_error = sqrt(sum(horizontal)/real(n, dp))
      write (*, '(a,f4.2,a,3(es14.6,1x))') 'panel vacuum z=', z, ' mean_V,lateral_rms_V,horizontal_rms_V_m=', &
        mean_error, rms_error, field_error
      call assert_true(abs(mean_error)/scale < 3.0e-5_dp, 'vacuum mean potential leakage')
      call assert_true(rms_error/scale < 5.0e-5_dp, 'vacuum lateral potential error')
      call assert_true(field_error/scale < 5.0e-4_dp, 'vacuum horizontal field error')
    end do
    call destroy_state(state)
    call destroy_plan(plan)
    call delete_periodic_cache(cache_path)
    call test_end()
  end subroutine check_panel_vacuum_reference

  !> Plane means of the k/=0 potential of a 1 pC point-charge sheet, relative to the sheet scale sigma*L/(2 eps0).
  !! The core works in raw 1/r units, where that scale is 2*pi*sigma*L.
  subroutine check_uniform_sheet_leak(bound)
    real(dp), intent(in) :: bound
    integer(i32), parameter :: nsheet = 4_i32, ncolumn = 4_i32, nplane = 16_i32
    real(dp), parameter :: cell = 1.0e-4_dp, sheet_charge = 1.0e-12_dp
    real(dp), parameter :: heights(3) = [0.4_dp, 1.0_dp, 1.9_dp]
    real(dp) :: sheet_source(3, nsheet*nsheet + ncolumn), sheet_q(nsheet*nsheet + ncolumn)
    real(dp) :: plane_points(3, nplane*nplane), plane_phi(nplane*nplane), scale, leak, worst
    integer(i32) :: i, j, k, n

    n = 0_i32
    do i = 1_i32, nsheet
      do j = 1_i32, nsheet
        n = n + 1_i32
        sheet_source(:, n) = cell*[(real(i, dp) - 0.5_dp)/nsheet, (real(j, dp) - 0.5_dp)/nsheet, 0.02_dp]
        sheet_q(n) = sheet_charge/real(nsheet*nsheet, dp)
      end do
    end do
    do k = 1_i32, ncolumn
      n = n + 1_i32
      sheet_source(:, n) = cell*[0.5_dp, 0.5_dp, 0.2_dp + 1.5_dp*real(k - 1_i32, dp)/real(ncolumn - 1_i32, dp)]
      sheet_q(n) = 0.0_dp
    end do
    options%periodic_len = cell
    options%periodic_image_layers = 1_i32
    options%target_box_max = cell*[1.0_dp, 1.0_dp, 2.0_dp]
    call build_plan(plan, sheet_source, options)
    cache_path = plan%periodic_cache_path
    call update_state(plan, state, sheet_q)
    scale = 2.0_dp*acos(-1.0_dp)*sheet_charge/(cell*cell)*cell
    worst = 0.0_dp
    do k = 1_i32, size(heights)
      n = 0_i32
      do i = 1_i32, nplane
        do j = 1_i32, nplane
          n = n + 1_i32
          plane_points(:, n) = cell*[(real(i, dp) - 0.5_dp)/nplane, (real(j, dp) - 0.5_dp)/nplane, heights(k)]
        end do
      end do
      call eval_potential_points(plan, state, plane_points, plane_phi)
      ! A NaN potential would make the plane mean NaN, which MAX may drop.
      call assert_true(all(ieee_is_finite(plane_phi)), 'uniform-sheet potential is not finite')
      leak = abs(sum(plane_phi))/real(nplane*nplane, dp)/scale
      write (*, '(a,f5.2,a,es10.3)') 'z/L=', heights(k), ' relative plane-mean leak=', leak
      worst = max(worst, leak)
    end do
    call assert_true(worst < bound, 'k/=0 potential of a uniform sheet must have no plane mean')
    call destroy_state(state)
    call destroy_plan(plan)
    call delete_periodic_cache(cache_path)
  end subroutine check_uniform_sheet_leak

end program test_periodic2_length_scale
