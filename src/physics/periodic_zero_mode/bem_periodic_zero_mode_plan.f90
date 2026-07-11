module bem_periodic_zero_mode_plan
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer(i32), parameter, public :: periodic_zero_mode_ok = 0_i32
  integer(i32), parameter, public :: periodic_zero_mode_invalid = 1_i32

  type :: panel_height_projection_type
    real(dp) :: z(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    integer(i32) :: break_index(3) = [0_i32, 0_i32, 0_i32]
    logical :: horizontal = .false.
  end type panel_height_projection_type

  type, public :: periodic_zero_mode_plan_type
    integer(i32) :: nelem = 0_i32
    integer(i32) :: nbreak = 0_i32
    real(dp) :: area_xy = 0.0_dp
    real(dp), allocatable :: break_z(:)
    type(panel_height_projection_type), allocatable :: panel(:)
  end type periodic_zero_mode_plan_type

  type, public :: periodic_zero_mode_state_type
    real(dp) :: e_bottom = 0.0_dp
    real(dp) :: z_gauge = 0.0_dp
    real(dp) :: phi_gauge = 0.0_dp
    real(dp) :: total_charge = 0.0_dp
    real(dp), allocatable :: cumulative_charge_coeff(:, :)
    real(dp), allocatable :: primitive_at_break(:)
    real(dp), allocatable :: sheet_charge(:)
  end type periodic_zero_mode_state_type

  public :: build_periodic_zero_mode_plan
  public :: build_periodic_zero_mode_height_plan
  public :: refresh_periodic_zero_mode_state

contains

  subroutine build_periodic_zero_mode_plan(mesh, area_xy, plan, status, message)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: area_xy
    type(periodic_zero_mode_plan_type), intent(out) :: plan
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), allocatable :: heights(:, :)
    integer(i32) :: elem

    status = periodic_zero_mode_invalid
    message = ''
    if (mesh%nelem <= 0_i32) then
      message = 'periodic zero-mode mesh must contain at least one element'
      return
    end if

    allocate (heights(3, mesh%nelem))
    do elem = 1, mesh%nelem
      heights(:, elem) = [mesh%v0(3, elem), mesh%v1(3, elem), mesh%v2(3, elem)]
    end do
    call build_periodic_zero_mode_height_plan(heights, area_xy, plan, status, message)
  end subroutine build_periodic_zero_mode_plan

  subroutine build_periodic_zero_mode_height_plan(source_heights, area_xy, plan, status, message)
    real(dp), intent(in) :: source_heights(:, :)
    real(dp), intent(in) :: area_xy
    type(periodic_zero_mode_plan_type), intent(out) :: plan
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), allocatable :: heights(:), unique_heights(:)
    real(dp) :: scale, tolerance
    integer(i32) :: elem, vertex, count, nelem

    status = periodic_zero_mode_invalid
    message = ''
    if (.not. ieee_is_finite(area_xy) .or. area_xy <= 0.0_dp) then
      message = 'periodic zero-mode area_xy must be finite and positive'
      return
    end if
    if (size(source_heights, 1) /= 3 .or. size(source_heights, 2) <= 0) then
      message = 'periodic zero-mode heights must have shape (3, nsource)'
      return
    end if
    if (any(.not. ieee_is_finite(source_heights))) then
      message = 'periodic zero-mode panel heights must be finite'
      return
    end if
    nelem = int(size(source_heights, 2), i32)
    allocate (heights(3*nelem))
    do elem = 1_i32, nelem
      heights(3*elem - 2:3*elem) = source_heights(:, elem)
    end do
    call sort_real(heights)
    scale = max(1.0_dp, maxval(abs(heights)))
    tolerance = 128.0_dp*epsilon(1.0_dp)*scale
    allocate (unique_heights(size(heights)))
    count = 1_i32
    unique_heights(1) = heights(1)
    do vertex = 2, size(heights)
      if (abs(heights(vertex) - unique_heights(count)) > tolerance) then
        count = count + 1_i32
        unique_heights(count) = heights(vertex)
      end if
    end do

    plan%nelem = nelem
    plan%nbreak = count
    plan%area_xy = area_xy
    allocate (plan%break_z(count), plan%panel(nelem))
    plan%break_z = unique_heights(1:count)
    do elem = 1, nelem
      plan%panel(elem)%z = source_heights(:, elem)
      call sort_three(plan%panel(elem)%z)
      plan%panel(elem)%horizontal = &
        abs(plan%panel(elem)%z(3) - plan%panel(elem)%z(1)) <= tolerance
      do vertex = 1, 3
        plan%panel(elem)%break_index(vertex) = nearest_break(plan%break_z, plan%panel(elem)%z(vertex))
      end do
    end do
    status = periodic_zero_mode_ok
  end subroutine build_periodic_zero_mode_height_plan

  subroutine refresh_periodic_zero_mode_state(plan, charge, e_bottom, z_gauge, phi_gauge, state)
    type(periodic_zero_mode_plan_type), intent(in) :: plan
    real(dp), intent(in) :: charge(:)
    real(dp), intent(in) :: e_bottom, z_gauge, phi_gauge
    type(periodic_zero_mode_state_type), intent(out) :: state
    real(dp), allocatable :: difference(:, :)
    real(dp) :: z0, z1, z2, denominator, coefficient(3)
    integer(i32) :: elem, interval, ninterval, i0, i1, i2

    if (size(charge) /= plan%nelem) error stop 'periodic zero-mode charge size mismatch.'
    if (any(.not. ieee_is_finite(charge))) error stop 'periodic zero-mode charges must be finite.'
    ninterval = plan%nbreak + 1_i32
    allocate (difference(3, ninterval + 1_i32))
    allocate (state%cumulative_charge_coeff(3, ninterval))
    allocate (state%primitive_at_break(plan%nbreak), state%sheet_charge(plan%nbreak))
    difference = 0.0_dp
    state%sheet_charge = 0.0_dp
    state%e_bottom = e_bottom
    state%z_gauge = z_gauge
    state%phi_gauge = phi_gauge
    state%total_charge = sum(charge)

    do elem = 1, plan%nelem
      z0 = plan%panel(elem)%z(1)
      z1 = plan%panel(elem)%z(2)
      z2 = plan%panel(elem)%z(3)
      i0 = plan%panel(elem)%break_index(1)
      i1 = plan%panel(elem)%break_index(2)
      i2 = plan%panel(elem)%break_index(3)
      if (plan%panel(elem)%horizontal) then
        state%sheet_charge(i0) = state%sheet_charge(i0) + charge(elem)
      else
        if (i1 > i0) then
          denominator = (z1 - z0)*(z2 - z0)
          coefficient = charge(elem)*[z0*z0/denominator, -2.0_dp*z0/denominator, 1.0_dp/denominator]
          call add_interval_range(difference, i0 + 1_i32, i1, coefficient)
        end if
        if (i2 > i1) then
          denominator = (z2 - z0)*(z2 - z1)
          coefficient = charge(elem)*[ &
                        1.0_dp - z2*z2/denominator, 2.0_dp*z2/denominator, -1.0_dp/denominator &
                        ]
          call add_interval_range(difference, i1 + 1_i32, i2, coefficient)
        end if
      end if
      call add_interval_range(difference, i2 + 1_i32, ninterval, [charge(elem), 0.0_dp, 0.0_dp])
    end do

    state%cumulative_charge_coeff(:, 1) = difference(:, 1)
    do interval = 2, ninterval
      state%cumulative_charge_coeff(:, interval) = &
        state%cumulative_charge_coeff(:, interval - 1) + difference(:, interval)
    end do
    state%primitive_at_break(1) = 0.0_dp
    do interval = 1, plan%nbreak - 1
      state%primitive_at_break(interval + 1) = state%primitive_at_break(interval) + &
                                               integrate_polynomial( &
                                               state%cumulative_charge_coeff(:, interval + 1), &
                                               plan%break_z(interval), plan%break_z(interval + 1) &
                                               )
    end do
  end subroutine refresh_periodic_zero_mode_state

  subroutine add_interval_range(difference, first, last, coefficient)
    real(dp), intent(inout) :: difference(:, :)
    integer(i32), intent(in) :: first, last
    real(dp), intent(in) :: coefficient(3)

    if (first > last) return
    difference(:, first) = difference(:, first) + coefficient
    difference(:, last + 1_i32) = difference(:, last + 1_i32) - coefficient
  end subroutine add_interval_range

  pure real(dp) function integrate_polynomial(coefficient, lower, upper) result(integral)
    real(dp), intent(in) :: coefficient(3), lower, upper

    integral = coefficient(1)*(upper - lower) + &
               0.5_dp*coefficient(2)*(upper*upper - lower*lower) + &
               coefficient(3)/3.0_dp*(upper**3 - lower**3)
  end function integrate_polynomial

  pure integer(i32) function nearest_break(break_z, value) result(index)
    real(dp), intent(in) :: break_z(:), value
    integer(i32) :: i

    index = 1_i32
    do i = 2, size(break_z)
      if (abs(break_z(i) - value) < abs(break_z(index) - value)) index = i
    end do
  end function nearest_break

  subroutine sort_real(values)
    real(dp), intent(inout) :: values(:)
    real(dp) :: key
    integer :: i, j

    do i = 2, size(values)
      key = values(i)
      j = i - 1
      do while (j >= 1)
        if (values(j) <= key) exit
        values(j + 1) = values(j)
        j = j - 1
      end do
      values(j + 1) = key
    end do
  end subroutine sort_real

  pure subroutine sort_three(values)
    real(dp), intent(inout) :: values(3)
    real(dp) :: temporary

    if (values(1) > values(2)) then
      temporary = values(1); values(1) = values(2); values(2) = temporary
    end if
    if (values(2) > values(3)) then
      temporary = values(2); values(2) = values(3); values(3) = temporary
    end if
    if (values(1) > values(2)) then
      temporary = values(1); values(1) = values(2); values(2) = temporary
    end if
  end subroutine sort_three

end module bem_periodic_zero_mode_plan
