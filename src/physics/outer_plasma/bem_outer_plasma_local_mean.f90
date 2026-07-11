module bem_outer_plasma_local_mean
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_invalid, outer_plasma_not_applicable
  use bem_types, only: mesh_type
  implicit none
  private

  public :: build_accessible_fraction_from_heights
  public :: combine_accessible_charge_density
  public :: residence_histogram_density
  public :: relative_local_mean_mismatch
  public :: sample_plasma_facing_height_field

contains

  subroutine sample_plasma_facing_height_field( &
    mesh, xy_min, xy_max, sample_nx, sample_ny, surface_height, &
    multiple_intersection_count, status &
    )
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: xy_min(2), xy_max(2)
    integer(i32), intent(in) :: sample_nx, sample_ny
    real(dp), allocatable, intent(out) :: surface_height(:)
    integer(i32), intent(out) :: multiple_intersection_count
    integer(i32), intent(out) :: status
    real(dp), allocatable :: intersections(:)
    real(dp) :: x, y, z_hit, tolerance
    integer(i32) :: ix, iy, sample, element, intersection_count
    logical :: hit

    status = outer_plasma_invalid
    multiple_intersection_count = 0_i32
    if (mesh%nelem < 1_i32 .or. sample_nx < 1_i32 .or. sample_ny < 1_i32 .or. &
        .not. all(ieee_is_finite([xy_min, xy_max])) .or. any(xy_max <= xy_min)) return
    allocate (surface_height(sample_nx*sample_ny), intersections(mesh%nelem))
    surface_height = 0.0_dp
    tolerance = 256.0_dp*epsilon(1.0_dp)*max( &
                1.0_dp, max(maxval(abs(mesh%v0)), max(maxval(abs(mesh%v1)), maxval(abs(mesh%v2)))) &
                )
    sample = 0_i32
    do iy = 1_i32, sample_ny
      y = xy_min(2) + (real(iy, dp) - 0.5_dp)*(xy_max(2) - xy_min(2))/real(sample_ny, dp)
      do ix = 1_i32, sample_nx
        x = xy_min(1) + (real(ix, dp) - 0.5_dp)*(xy_max(1) - xy_min(1))/real(sample_nx, dp)
        sample = sample + 1_i32
        intersection_count = 0_i32
        do element = 1_i32, mesh%nelem
          call vertical_triangle_intersection(mesh, element, x, y, tolerance, hit, z_hit)
          if (.not. hit) cycle
          if (intersection_count == 0_i32 .or. &
              all(abs(intersections(1:intersection_count) - z_hit) > tolerance)) then
            intersection_count = intersection_count + 1_i32
            intersections(intersection_count) = z_hit
          end if
        end do
        if (intersection_count == 0_i32) return
        surface_height(sample) = maxval(intersections(1:intersection_count))
        if (intersection_count > 1_i32) multiple_intersection_count = multiple_intersection_count + 1_i32
      end do
    end do
    if (multiple_intersection_count > 0_i32) then
      status = outer_plasma_not_applicable
    else
      status = outer_plasma_ok
    end if
  end subroutine sample_plasma_facing_height_field

  subroutine vertical_triangle_intersection(mesh, element, x, y, tolerance, hit, z)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: element
    real(dp), intent(in) :: x, y, tolerance
    logical, intent(out) :: hit
    real(dp), intent(out) :: z
    real(dp) :: x0, y0, x1, y1, x2, y2, denominator, w0, w1, w2

    x0 = mesh%v0(1, element)
    y0 = mesh%v0(2, element)
    x1 = mesh%v1(1, element)
    y1 = mesh%v1(2, element)
    x2 = mesh%v2(1, element)
    y2 = mesh%v2(2, element)
    denominator = (y1 - y2)*(x0 - x2) + (x2 - x1)*(y0 - y2)
    hit = .false.
    z = 0.0_dp
    if (abs(denominator) <= tolerance) return
    w0 = ((y1 - y2)*(x - x2) + (x2 - x1)*(y - y2))/denominator
    w1 = ((y2 - y0)*(x - x2) + (x0 - x2)*(y - y2))/denominator
    w2 = 1.0_dp - w0 - w1
    if (min(w0, w1, w2) < -tolerance .or. max(w0, w1, w2) > 1.0_dp + tolerance) return
    z = w0*mesh%v0(3, element) + w1*mesh%v1(3, element) + w2*mesh%v2(3, element)
    hit = ieee_is_finite(z)
  end subroutine vertical_triangle_intersection

  subroutine build_accessible_fraction_from_heights(z, surface_height, accessible_fraction, status)
    real(dp), intent(in) :: z(:), surface_height(:)
    real(dp), intent(out) :: accessible_fraction(:)
    integer(i32), intent(out) :: status
    integer(i32) :: level

    accessible_fraction = 0.0_dp
    status = outer_plasma_invalid
    if (size(z) /= size(accessible_fraction) .or. size(z) < 1 .or. size(surface_height) < 1) return
    if (.not. all(ieee_is_finite(z)) .or. .not. all(ieee_is_finite(surface_height))) return
    if (size(z) > 1 .and. any(z(2:) <= z(:size(z) - 1))) return
    do level = 1_i32, size(z)
      accessible_fraction(level) = real(count(surface_height < z(level)), dp)/real(size(surface_height), dp)
    end do
    status = outer_plasma_ok
  end subroutine build_accessible_fraction_from_heights

  subroutine combine_accessible_charge_density(accessible_fraction, conditional_density, full_cell_density, status)
    real(dp), intent(in) :: accessible_fraction(:), conditional_density(:)
    real(dp), intent(out) :: full_cell_density(:)
    integer(i32), intent(out) :: status

    full_cell_density = 0.0_dp
    status = outer_plasma_invalid
    if (size(accessible_fraction) /= size(conditional_density) .or. &
        size(full_cell_density) /= size(accessible_fraction)) return
    if (.not. all(ieee_is_finite(accessible_fraction)) .or. &
        .not. all(ieee_is_finite(conditional_density)) .or. &
        any(accessible_fraction < 0.0_dp) .or. any(accessible_fraction > 1.0_dp)) return
    full_cell_density = accessible_fraction*conditional_density
    status = outer_plasma_ok
  end subroutine combine_accessible_charge_density

  subroutine residence_histogram_density( &
    weight_by_bin, area_xy, accessible_fraction, bin_width, observation_time, density, status &
    )
    real(dp), intent(in) :: weight_by_bin(:), area_xy, accessible_fraction(:), bin_width(:), observation_time
    real(dp), intent(out) :: density(:)
    integer(i32), intent(out) :: status

    density = 0.0_dp
    status = outer_plasma_invalid
    if (size(weight_by_bin) /= size(accessible_fraction) .or. &
        size(bin_width) /= size(accessible_fraction) .or. size(density) /= size(accessible_fraction)) return
    if (.not. all(ieee_is_finite([weight_by_bin, accessible_fraction, bin_width])) .or. &
        .not. all(ieee_is_finite([area_xy, observation_time])) .or. area_xy <= 0.0_dp .or. &
        observation_time <= 0.0_dp .or. any(weight_by_bin < 0.0_dp) .or. any(bin_width <= 0.0_dp) .or. &
        any(accessible_fraction < 0.0_dp) .or. any(accessible_fraction > 1.0_dp)) return
    if (any(accessible_fraction <= 0.0_dp .and. weight_by_bin > 0.0_dp)) then
      status = outer_plasma_not_applicable
      return
    end if
    where (accessible_fraction > 0.0_dp)
      density = weight_by_bin/(area_xy*accessible_fraction*bin_width*observation_time)
    elsewhere
      density = 0.0_dp
    end where
    status = outer_plasma_ok
  end subroutine residence_histogram_density

  pure real(dp) function relative_local_mean_mismatch(closure_density, histogram_density) result(mismatch)
    real(dp), intent(in) :: closure_density(:), histogram_density(:)
    real(dp) :: scale
    integer(i32) :: bin

    if (size(closure_density) /= size(histogram_density) .or. size(closure_density) < 1) then
      mismatch = huge(1.0_dp)
      return
    end if
    mismatch = 0.0_dp
    do bin = 1_i32, size(closure_density)
      scale = max(abs(closure_density(bin)), abs(histogram_density(bin)), tiny(1.0_dp))
      mismatch = max(mismatch, abs(closure_density(bin) - histogram_density(bin))/scale)
    end do
  end function relative_local_mean_mismatch

end module bem_outer_plasma_local_mean
