module bem_outer_plasma_grid
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  implicit none
  private

  type, public :: outer_plasma_grid_type
    integer(i32) :: n = 0_i32
    real(dp) :: length = 0.0_dp
    real(dp) :: stretch = 0.0_dp
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: dz(:)
  end type outer_plasma_grid_type

  public :: init_outer_plasma_grid
  public :: interpolate_outer_profile

contains

  subroutine init_outer_plasma_grid(n, length, stretch, grid)
    integer(i32), intent(in) :: n
    real(dp), intent(in) :: length, stretch
    type(outer_plasma_grid_type), intent(out) :: grid
    integer(i32) :: j
    real(dp) :: coordinate, denominator

    if (n < 3_i32 .or. .not. all(ieee_is_finite([length, stretch])) .or. &
        length <= 0.0_dp .or. stretch < 0.0_dp) then
      error stop 'Outer-plasma grid requires n >= 3, positive length, and nonnegative stretch.'
    end if
    grid%n = n
    grid%length = length
    grid%stretch = stretch
    allocate (grid%z(n), grid%dz(n - 1_i32))
    if (stretch <= sqrt(epsilon(1.0_dp))) then
      do j = 1_i32, n
        grid%z(j) = length*real(j - 1_i32, dp)/real(n - 1_i32, dp)
      end do
    else
      denominator = exp(stretch) - 1.0_dp
      do j = 1_i32, n
        coordinate = real(j - 1_i32, dp)/real(n - 1_i32, dp)
        grid%z(j) = length*(exp(stretch*coordinate) - 1.0_dp)/denominator
      end do
    end if
    grid%z(1) = 0.0_dp
    grid%z(n) = length
    grid%dz = grid%z(2:n) - grid%z(1:n - 1_i32)
  end subroutine init_outer_plasma_grid

  subroutine interpolate_outer_profile(grid, profile, z, value)
    type(outer_plasma_grid_type), intent(in) :: grid
    real(dp), intent(in) :: profile(:), z
    real(dp), intent(out) :: value
    integer(i32) :: low, high, middle
    real(dp) :: weight

    if (size(profile) /= grid%n .or. grid%n < 2_i32 .or. .not. ieee_is_finite(z)) then
      error stop 'Invalid outer-plasma profile interpolation input.'
    end if
    if (z <= grid%z(1)) then
      value = profile(1)
      return
    else if (z >= grid%z(grid%n)) then
      value = profile(grid%n)
      return
    end if
    low = 1_i32
    high = grid%n
    do while (high - low > 1_i32)
      middle = (low + high)/2_i32
      if (grid%z(middle) <= z) then
        low = middle
      else
        high = middle
      end if
    end do
    weight = (z - grid%z(low))/grid%dz(low)
    value = (1.0_dp - weight)*profile(low) + weight*profile(high)
  end subroutine interpolate_outer_profile

end module bem_outer_plasma_grid
