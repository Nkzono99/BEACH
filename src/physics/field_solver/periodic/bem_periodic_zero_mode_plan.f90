module bem_periodic_zero_mode_plan
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
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
    !> Geometric part of the height-coincidence tolerance (zero_mode_heights_coincide).
    real(dp) :: height_tolerance = 0.0_dp
    real(dp), allocatable :: break_z(:)
    type(panel_height_projection_type), allocatable :: panel(:)
  end type periodic_zero_mode_plan_type

  type, public :: periodic_zero_mode_state_type
    real(dp) :: e_bottom = 0.0_dp
    real(dp) :: z_gauge = 0.0_dp
    real(dp) :: phi_gauge = 0.0_dp
    real(dp) :: total_charge = 0.0_dp
    !> 区間 k の累積電荷多項式の係数。k>=2 は左端 break_z(k-1)、k=1 は break_z(1) を原点とする局所座標
    !! t=z-原点 の多項式で、絶対 z の冪の相殺を避ける（zero_mode_interval_origin）。
    real(dp), allocatable :: cumulative_charge_coeff(:, :)
    real(dp), allocatable :: primitive_at_break(:)
    real(dp), allocatable :: sheet_charge(:)
  end type periodic_zero_mode_state_type

  public :: build_periodic_zero_mode_plan
  public :: build_periodic_zero_mode_height_plan
  public :: refresh_periodic_zero_mode_state
  public :: symmetric_vacuum_bottom_field
  public :: zero_mode_interval_origin
  public :: integrate_local_polynomial
  public :: zero_mode_heights_coincide

contains

  pure real(dp) function symmetric_vacuum_bottom_field(plan, charge) result(field)
    type(periodic_zero_mode_plan_type), intent(in) :: plan
    real(dp), intent(in) :: charge(:)

    field = -sum(charge)/(2.0_dp*eps0*plan%area_xy)
  end function symmetric_vacuum_bottom_field

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
    real(dp) :: tolerance
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
    ! Geometric scale of the problem: the slab height and the lateral cell size. Neither depends on
    ! the absolute z origin or on a 1 m unit; sqrt(area_xy) keeps a flat slab from treating
    ! rounding-level tilts (eps times the lateral coordinates) as inclined panels.
    tolerance = 128.0_dp*epsilon(1.0_dp)*max(heights(size(heights)) - heights(1), sqrt(area_xy))
    allocate (unique_heights(size(heights)))
    count = 1_i32
    unique_heights(1) = heights(1)
    do vertex = 2, size(heights)
      if (.not. zero_mode_heights_coincide(heights(vertex), unique_heights(count), tolerance)) then
        count = count + 1_i32
        unique_heights(count) = heights(vertex)
      end if
    end do

    plan%nelem = nelem
    plan%nbreak = count
    plan%area_xy = area_xy
    plan%height_tolerance = tolerance
    allocate (plan%break_z(count), plan%panel(nelem))
    plan%break_z = unique_heights(1:count)
    do elem = 1, nelem
      plan%panel(elem)%z = source_heights(:, elem)
      call sort_three(plan%panel(elem)%z)
      plan%panel(elem)%horizontal = &
        zero_mode_heights_coincide(plan%panel(elem)%z(3), plan%panel(elem)%z(1), tolerance)
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
    real(dp) :: z0, z1, z2, denominator, offset, shift
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

    ! 傾斜 panel の寄与は、加える区間と外す区間それぞれの局所原点で表す。係数は panel 自身の高さ幅だけで
    ! 決まり、絶対 z に依らない（平行移動で値が変わらない）。
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
          ! [z0, z1]: q*(z-z0)^2/denominator。原点 a では z-z0 = t+offset, offset=a-z0。
          denominator = (z1 - z0)*(z2 - z0)
          offset = plan%break_z(i0) - z0
          difference(:, i0 + 1_i32) = difference(:, i0 + 1_i32) + &
                                      charge(elem)*[offset*offset, 2.0_dp*offset, 1.0_dp]/denominator
          offset = plan%break_z(i1) - z0
          difference(:, i1 + 1_i32) = difference(:, i1 + 1_i32) - &
                                      charge(elem)*[offset*offset, 2.0_dp*offset, 1.0_dp]/denominator
        end if
        if (i2 > i1) then
          ! [z1, z2]: q*(1-(z2-z)^2/denominator)。原点 a では z2-z = offset-t, offset=z2-a。
          denominator = (z2 - z0)*(z2 - z1)
          offset = z2 - plan%break_z(i1)
          difference(:, i1 + 1_i32) = difference(:, i1 + 1_i32) + charge(elem)*[ &
                                      1.0_dp - offset*offset/denominator, 2.0_dp*offset/denominator, &
                                      -1.0_dp/denominator]
          offset = z2 - plan%break_z(i2)
          difference(:, i2 + 1_i32) = difference(:, i2 + 1_i32) - charge(elem)*[ &
                                      1.0_dp - offset*offset/denominator, 2.0_dp*offset/denominator, &
                                      -1.0_dp/denominator]
        end if
      end if
      difference(1, i2 + 1_i32) = difference(1, i2 + 1_i32) + charge(elem)
    end do

    ! 前の区間の多項式を次の区間の原点へ移してから差分を足す。移す量は隣接 break の間隔だけ。
    state%cumulative_charge_coeff(:, 1) = difference(:, 1)
    do interval = 2, ninterval
      shift = zero_mode_interval_origin(plan%break_z, interval) - &
              zero_mode_interval_origin(plan%break_z, interval - 1_i32)
      state%cumulative_charge_coeff(:, interval) = &
        shift_polynomial(state%cumulative_charge_coeff(:, interval - 1), shift) + difference(:, interval)
    end do
    ! Outside the source slab the integrated charge is known exactly. Do not
    ! extrapolate cancellation residue from the inclined-panel polynomials.
    state%cumulative_charge_coeff(:, 1) = 0.0_dp
    state%cumulative_charge_coeff(:, ninterval) = [state%total_charge, 0.0_dp, 0.0_dp]
    state%primitive_at_break(1) = 0.0_dp
    do interval = 1, plan%nbreak - 1
      state%primitive_at_break(interval + 1) = state%primitive_at_break(interval) + &
                                               integrate_local_polynomial( &
                                               state%cumulative_charge_coeff(:, interval + 1), &
                                               plan%break_z(interval + 1) - plan%break_z(interval) &
                                               )
    end do
  end subroutine refresh_periodic_zero_mode_state

  !> Two heights are one breakpoint when they differ by at most the plan's geometric tolerance or by
  !! the few ULPs that binary64 cannot resolve at their magnitude. Only the second term sees the
  !! absolute coordinate, so a representable tilt (tens of ULPs) far from the origin stays inclined.
  pure logical function zero_mode_heights_coincide(a, b, geometric_tolerance) result(coincide)
    real(dp), intent(in) :: a, b, geometric_tolerance

    coincide = abs(a - b) <= max(geometric_tolerance, 8.0_dp*spacing(max(abs(a), abs(b))))
  end function zero_mode_heights_coincide

  !> 区間 interval の累積電荷多項式の原点。区間 1 は最下 break、それ以外は区間の左端 break。
  pure real(dp) function zero_mode_interval_origin(break_z, interval) result(origin)
    real(dp), intent(in) :: break_z(:)
    integer(i32), intent(in) :: interval

    origin = break_z(max(interval - 1_i32, 1_i32))
  end function zero_mode_interval_origin

  !> 局所座標の多項式を原点から t=0..t_upper まで積分する。
  pure real(dp) function integrate_local_polynomial(coefficient, t_upper) result(integral)
    real(dp), intent(in) :: coefficient(3), t_upper

    integral = t_upper*(coefficient(1) + t_upper*(0.5_dp*coefficient(2) + t_upper*coefficient(3)/3.0_dp))
  end function integrate_local_polynomial

  !> c0+c1*t+c2*t^2 を原点を shift だけ進めた座標 t'=t-shift の係数へ書き換える。
  pure function shift_polynomial(coefficient, shift) result(shifted)
    real(dp), intent(in) :: coefficient(3), shift
    real(dp) :: shifted(3)

    shifted(1) = coefficient(1) + shift*(coefficient(2) + shift*coefficient(3))
    shifted(2) = coefficient(2) + 2.0_dp*shift*coefficient(3)
    shifted(3) = coefficient(3)
  end function shift_polynomial

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
