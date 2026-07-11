program test_outer_plasma_local_mean
  use bem_kinds, only: dp, i32
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_not_applicable
  use bem_outer_plasma_local_mean, only: &
    build_accessible_fraction_from_heights, combine_accessible_charge_density, &
    residence_histogram_density, relative_local_mean_mismatch, sample_plasma_facing_height_field
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32
  implicit none

  real(dp) :: fraction(4), full_density(4), histogram_density(2), mismatch
  integer(i32) :: status
  integer(i32) :: multiple_intersection_count
  type(mesh_type) :: mesh
  real(dp), allocatable :: height(:)

  call test_init(6)

  call test_begin('height-field samples define accessible area monotonically')
  call build_accessible_fraction_from_heights( &
    [0.5_dp, 1.5_dp, 2.5_dp, 3.5_dp], [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], fraction, status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'accessible-area status mismatch')
  call assert_close_dp(fraction(1), 0.25_dp, 0.0_dp, 'accessible fraction 1 mismatch')
  call assert_close_dp(fraction(2), 0.50_dp, 0.0_dp, 'accessible fraction 2 mismatch')
  call assert_close_dp(fraction(3), 0.75_dp, 0.0_dp, 'accessible fraction 3 mismatch')
  call assert_close_dp(fraction(4), 1.00_dp, 0.0_dp, 'accessible fraction 4 mismatch')
  call test_end()

  call test_begin('rough single-valued mesh is sampled at cell centers')
  call build_sloped_mesh(mesh, duplicate_top=.false.)
  call sample_plasma_facing_height_field( &
    mesh, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], 2_i32, 2_i32, height, &
    multiple_intersection_count, status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'rough height-field status mismatch')
  call assert_equal_i32(multiple_intersection_count, 0_i32, 'single-valued mesh reported an overhang')
  call assert_close_dp(height(1), 0.25_dp, 1.0e-14_dp, 'rough height sample 1 mismatch')
  call assert_close_dp(height(2), 0.75_dp, 1.0e-14_dp, 'rough height sample 2 mismatch')
  call assert_close_dp(height(3), 0.25_dp, 1.0e-14_dp, 'rough height sample 3 mismatch')
  call assert_close_dp(height(4), 0.75_dp, 1.0e-14_dp, 'rough height sample 4 mismatch')
  call test_end()

  call test_begin('multiple vertical intersections fail the 1D geometry contract')
  call build_sloped_mesh(mesh, duplicate_top=.true.)
  call sample_plasma_facing_height_field( &
    mesh, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], 2_i32, 2_i32, height, &
    multiple_intersection_count, status &
    )
  call assert_equal_i32(status, outer_plasma_not_applicable, 'overhang geometry must fail closed')
  call assert_equal_i32(multiple_intersection_count, 4_i32, 'overhang column count mismatch')
  call test_end()

  call test_begin('conditional closure density is weighted exactly once')
  call combine_accessible_charge_density(fraction, [-4.0_dp, -4.0_dp, -4.0_dp, -4.0_dp], full_density, status)
  call assert_equal_i32(status, outer_plasma_ok, 'accessible charge status mismatch')
  call assert_close_dp(full_density(1), -1.0_dp, 0.0_dp, 'accessible charge 1 mismatch')
  call assert_close_dp(full_density(4), -4.0_dp, 0.0_dp, 'accessible charge 4 mismatch')
  call test_end()

  call test_begin('residence histogram uses accessible volume and observation time')
  call residence_histogram_density( &
    weight_by_bin=[20.0_dp, 30.0_dp], area_xy=10.0_dp, accessible_fraction=[0.5_dp, 1.0_dp], &
    bin_width=[2.0_dp, 3.0_dp], observation_time=2.0_dp, density=histogram_density, status=status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'residence histogram status mismatch')
  call assert_close_dp(histogram_density(1), 1.0_dp, 1.0e-15_dp, 'residence density 1 mismatch')
  call assert_close_dp(histogram_density(2), 0.5_dp, 1.0e-15_dp, 'residence density 2 mismatch')
  call test_end()

  call test_begin('zero accessible volume fails closed and mismatch is symmetric')
  call residence_histogram_density( &
    weight_by_bin=[1.0_dp, 1.0_dp], area_xy=1.0_dp, accessible_fraction=[0.0_dp, 1.0_dp], &
    bin_width=[1.0_dp, 1.0_dp], observation_time=1.0_dp, density=histogram_density, status=status &
    )
  call assert_equal_i32(status, outer_plasma_not_applicable, 'zero accessible volume must fail closed')
  mismatch = relative_local_mean_mismatch([2.0_dp, 4.0_dp], [1.0_dp, 5.0_dp])
  call assert_close_dp(mismatch, 0.5_dp, 1.0e-15_dp, 'local mean mismatch metric mismatch')
  call test_end()

  call test_summary()

contains

  subroutine build_sloped_mesh(value, duplicate_top)
    type(mesh_type), intent(out) :: value
    logical, intent(in) :: duplicate_top
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    integer(i32) :: n

    n = merge(4_i32, 2_i32, duplicate_top)
    allocate (v0(3, n), v1(3, n), v2(3, n))
    v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 1.0_dp]
    v2(:, 1) = [1.0_dp, 1.0_dp, 1.0_dp]
    v0(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp]
    v1(:, 2) = [1.0_dp, 1.0_dp, 1.0_dp]
    v2(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp]
    if (duplicate_top) then
      v0(:, 3:4) = v0(:, 1:2)
      v1(:, 3:4) = v1(:, 1:2)
      v2(:, 3:4) = v2(:, 1:2)
      v0(3, 3:4) = v0(3, 3:4) + 1.0_dp
      v1(3, 3:4) = v1(3, 3:4) + 1.0_dp
      v2(3, 3:4) = v2(3, 3:4) + 1.0_dp
    end if
    call init_mesh(value, v0, v1, v2)
  end subroutine build_sloped_mesh
end program test_outer_plasma_local_mean
