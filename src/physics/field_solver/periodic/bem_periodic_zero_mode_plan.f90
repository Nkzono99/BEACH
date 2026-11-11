module bem_periodic_zero_mode_plan
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type
  use bem_sorting, only: sort_real_ascending
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer(i32), parameter, public :: periodic_zero_mode_ok = 0_i32
  integer(i32), parameter, public :: periodic_zero_mode_invalid = 1_i32

  ! 1 枚の panel が累積電荷多項式の区間差分へ足す寄与。番号の順が 1 枚の中で足す順。
  integer(i32), parameter :: zero_mode_term_lower_enter = 1_i32  ! [z0, z1] の多項式を区間 i0+1 から加える
  integer(i32), parameter :: zero_mode_term_lower_leave = 2_i32  ! 同じ多項式を区間 i1+1 から外す
  integer(i32), parameter :: zero_mode_term_upper_enter = 3_i32  ! [z1, z2] の多項式を区間 i1+1 から加える
  integer(i32), parameter :: zero_mode_term_upper_leave = 4_i32  ! 同じ多項式を区間 i2+1 から外す
  integer(i32), parameter :: zero_mode_term_total = 5_i32        ! panel の全電荷を区間 i2+1 から加える

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
    !> 区間ごとの寄与の一覧。区間 k の寄与は `term_start(k)` から `term_start(k+1)-1` で、panel 番号の昇順、
    !! 同じ panel の中では refresh が足す順に並ぶ。`term_kind` は寄与の種類（zero_mode_term_*）。
    integer(i32), allocatable :: term_start(:), term_elem(:), term_kind(:)
  end type periodic_zero_mode_plan_type

  type, public :: periodic_zero_mode_state_type
    real(dp) :: e_bottom = 0.0_dp
    real(dp) :: z_gauge = 0.0_dp
    real(dp) :: phi_gauge = 0.0_dp
    real(dp) :: total_charge = 0.0_dp
    !> Charge primitive at z_gauge; every evaluation subtracts it, so refresh computes it once.
    real(dp) :: gauge_primitive = 0.0_dp
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
  public :: locate_zero_mode_interval
  public :: zero_mode_charge_primitive

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
    call sort_real_ascending(heights)
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
    call build_zero_mode_interval_terms(plan)
    status = periodic_zero_mode_ok
  end subroutine build_periodic_zero_mode_height_plan

  !> panel の寄与を区間ごとの一覧に並べる。区間ごとに panel の昇順・panel の中では種類の順なので、
  !! 一覧の順に足すと、panel を順に回って足した場合と同じ順になる。
  subroutine build_zero_mode_interval_terms(plan)
    type(periodic_zero_mode_plan_type), intent(inout) :: plan
    integer(i32) :: elem, pass, k, nterm, interval, ninterval
    integer(i32) :: term_interval(5), term_kind(5)
    integer(i32), allocatable :: fill(:)

    ninterval = plan%nbreak + 1_i32
    allocate (plan%term_start(ninterval + 1_i32), fill(ninterval))
    plan%term_start = 0_i32
    do pass = 1, 2
      do elem = 1, plan%nelem
        call panel_interval_terms(plan%panel(elem), nterm, term_interval, term_kind)
        do k = 1, nterm
          interval = term_interval(k)
          if (pass == 1) then
            plan%term_start(interval + 1_i32) = plan%term_start(interval + 1_i32) + 1_i32
          else
            plan%term_elem(fill(interval)) = elem
            plan%term_kind(fill(interval)) = term_kind(k)
            fill(interval) = fill(interval) + 1_i32
          end if
        end do
      end do
      if (pass == 1) then
        plan%term_start(1) = 1_i32
        do interval = 1, ninterval
          plan%term_start(interval + 1_i32) = plan%term_start(interval) + plan%term_start(interval + 1_i32)
        end do
        allocate (plan%term_elem(plan%term_start(ninterval + 1_i32) - 1_i32), &
                  plan%term_kind(plan%term_start(ninterval + 1_i32) - 1_i32))
        fill = plan%term_start(1:ninterval)
      end if
    end do
  end subroutine build_zero_mode_interval_terms

  !> 1 枚の panel の寄与を、足す順に（区間、種類）で返す。
  pure subroutine panel_interval_terms(panel, nterm, term_interval, term_kind)
    type(panel_height_projection_type), intent(in) :: panel
    integer(i32), intent(out) :: nterm, term_interval(5), term_kind(5)
    integer(i32) :: i0, i1, i2

    i0 = panel%break_index(1)
    i1 = panel%break_index(2)
    i2 = panel%break_index(3)
    nterm = 0_i32
    if (.not. panel%horizontal) then
      if (i1 > i0) then
        term_interval(nterm + 1:nterm + 2) = [i0 + 1_i32, i1 + 1_i32]
        term_kind(nterm + 1:nterm + 2) = [zero_mode_term_lower_enter, zero_mode_term_lower_leave]
        nterm = nterm + 2_i32
      end if
      if (i2 > i1) then
        term_interval(nterm + 1:nterm + 2) = [i1 + 1_i32, i2 + 1_i32]
        term_kind(nterm + 1:nterm + 2) = [zero_mode_term_upper_enter, zero_mode_term_upper_leave]
        nterm = nterm + 2_i32
      end if
    end if
    nterm = nterm + 1_i32
    term_interval(nterm) = i2 + 1_i32
    term_kind(nterm) = zero_mode_term_total
  end subroutine panel_interval_terms

  subroutine refresh_periodic_zero_mode_state(plan, charge, e_bottom, z_gauge, phi_gauge, state)
    type(periodic_zero_mode_plan_type), intent(in) :: plan
    real(dp), intent(in) :: charge(:)
    real(dp), intent(in) :: e_bottom, z_gauge, phi_gauge
    !> Arrays already sized for this plan are reused; every component is overwritten.
    type(periodic_zero_mode_state_type), intent(inout) :: state
    real(dp) :: z0, z1, z2, denominator, offset, shift, interval_sum(3)
    integer(i32) :: elem, interval, ninterval, i0, term, gauge_interval, gauge_breakpoint

    if (size(charge) /= plan%nelem) error stop 'periodic zero-mode charge size mismatch.'
    if (any(.not. ieee_is_finite(charge))) error stop 'periodic zero-mode charges must be finite.'
    ninterval = plan%nbreak + 1_i32
    if (allocated(state%cumulative_charge_coeff)) then
      if (size(state%cumulative_charge_coeff, 2) /= ninterval) deallocate (state%cumulative_charge_coeff)
    end if
    if (.not. allocated(state%cumulative_charge_coeff)) allocate (state%cumulative_charge_coeff(3, ninterval))
    call size_break_array(state%primitive_at_break, plan%nbreak)
    call size_break_array(state%sheet_charge, plan%nbreak)
    ! The per-interval differences are stored in cumulative_charge_coeff and become prefix sums in place.
    state%sheet_charge = 0.0_dp
    state%e_bottom = e_bottom
    state%z_gauge = z_gauge
    state%phi_gauge = phi_gauge
    state%total_charge = sum(charge)

    do elem = 1, plan%nelem
      if (plan%panel(elem)%horizontal) then
        i0 = plan%panel(elem)%break_index(1)
        state%sheet_charge(i0) = state%sheet_charge(i0) + charge(elem)
      end if
    end do

    ! 傾斜 panel の寄与は、加える区間と外す区間それぞれの局所原点で表す。係数は panel 自身の高さ幅だけで
    ! 決まり、絶対 z に依らない（平行移動で値が変わらない）。
    ! 区間ごとに、その区間への寄与を panel の順に足す（plan の一覧の順）。区間どうしは独立なので並列に求められ、
    ! 足す順は panel を順に回る場合と同じなので、値は thread 数によらない。
    !$omp parallel do default(none) schedule(static) shared(plan, charge, state, ninterval) &
    !$omp private(interval, term, elem, z0, z1, z2, denominator, offset, interval_sum)
    do interval = 1, ninterval
      interval_sum = 0.0_dp
      do term = plan%term_start(interval), plan%term_start(interval + 1_i32) - 1_i32
        elem = plan%term_elem(term)
        z0 = plan%panel(elem)%z(1)
        z1 = plan%panel(elem)%z(2)
        z2 = plan%panel(elem)%z(3)
        select case (plan%term_kind(term))
        case (zero_mode_term_lower_enter)
          ! [z0, z1]: q*(z-z0)^2/denominator。原点 a では z-z0 = t+offset, offset=a-z0。
          denominator = (z1 - z0)*(z2 - z0)
          offset = plan%break_z(plan%panel(elem)%break_index(1)) - z0
          interval_sum = interval_sum + charge(elem)*[offset*offset, 2.0_dp*offset, 1.0_dp]/denominator
        case (zero_mode_term_lower_leave)
          denominator = (z1 - z0)*(z2 - z0)
          offset = plan%break_z(plan%panel(elem)%break_index(2)) - z0
          interval_sum = interval_sum - charge(elem)*[offset*offset, 2.0_dp*offset, 1.0_dp]/denominator
        case (zero_mode_term_upper_enter)
          ! [z1, z2]: q*(1-(z2-z)^2/denominator)。原点 a では z2-z = offset-t, offset=z2-a。
          denominator = (z2 - z0)*(z2 - z1)
          offset = z2 - plan%break_z(plan%panel(elem)%break_index(2))
          interval_sum = interval_sum + charge(elem)*[ &
                         1.0_dp - offset*offset/denominator, 2.0_dp*offset/denominator, -1.0_dp/denominator]
        case (zero_mode_term_upper_leave)
          denominator = (z2 - z0)*(z2 - z1)
          offset = z2 - plan%break_z(plan%panel(elem)%break_index(3))
          interval_sum = interval_sum - charge(elem)*[ &
                         1.0_dp - offset*offset/denominator, 2.0_dp*offset/denominator, -1.0_dp/denominator]
        case (zero_mode_term_total)
          interval_sum(1) = interval_sum(1) + charge(elem)
        end select
      end do
      state%cumulative_charge_coeff(:, interval) = interval_sum
    end do
    !$omp end parallel do

    ! 前の区間の多項式を次の区間の原点へ移してから差分を足す。移す量は隣接 break の間隔だけ。
    do interval = 2, ninterval
      shift = zero_mode_interval_origin(plan%break_z, interval) - &
              zero_mode_interval_origin(plan%break_z, interval - 1_i32)
      state%cumulative_charge_coeff(:, interval) = &
        shift_polynomial(state%cumulative_charge_coeff(:, interval - 1), shift) + state%cumulative_charge_coeff(:, interval)
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
    call locate_zero_mode_interval(plan%break_z, plan%height_tolerance, z_gauge, gauge_interval, gauge_breakpoint)
    state%gauge_primitive = zero_mode_charge_primitive(plan, state, z_gauge, gauge_interval)
  end subroutine refresh_periodic_zero_mode_state

  subroutine size_break_array(values, nbreak)
    real(dp), allocatable, intent(inout) :: values(:)
    integer(i32), intent(in) :: nbreak

    if (allocated(values)) then
      if (size(values) /= nbreak) deallocate (values)
    end if
    if (.not. allocated(values)) allocate (values(nbreak))
  end subroutine size_break_array

  pure real(dp) function zero_mode_charge_primitive(plan, state, z, interval) result(primitive)
    type(periodic_zero_mode_plan_type), intent(in) :: plan
    type(periodic_zero_mode_state_type), intent(in) :: state
    real(dp), intent(in) :: z
    integer(i32), intent(in) :: interval

    ! 区間 1 は最下 break から、それ以外は区間の左端 break（局所原点）から積分する。
    primitive = integrate_local_polynomial( &
                state%cumulative_charge_coeff(:, interval), z - zero_mode_interval_origin(plan%break_z, interval) &
                )
    if (interval > 1_i32) primitive = primitive + state%primitive_at_break(interval - 1_i32)
  end function zero_mode_charge_primitive

  !> A target on a breakpoint (same coincidence rule as the plan's breakpoint merge) takes the
  !! requested sheet trace; any other target lies strictly inside one interval.
  pure subroutine locate_zero_mode_interval(break_z, tolerance, z, interval, breakpoint)
    real(dp), intent(in) :: break_z(:), tolerance, z
    integer(i32), intent(out) :: interval, breakpoint
    integer(i32) :: low, high, middle

    low = 1_i32
    high = size(break_z)
    breakpoint = 0_i32
    do while (low <= high)
      middle = (low + high)/2_i32
      if (zero_mode_heights_coincide(z, break_z(middle), tolerance)) then
        breakpoint = middle
        interval = middle + 1_i32
        return
      else if (z < break_z(middle)) then
        high = middle - 1_i32
      else
        low = middle + 1_i32
      end if
    end do
    interval = low
  end subroutine locate_zero_mode_interval

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
    integer(i32) :: low, high, middle

    ! break_z is sorted and nonempty; only the two neighbors can be nearest.
    low = 1_i32
    high = size(break_z)
    do while (low < high)
      middle = low + (high - low)/2_i32
      if (break_z(middle) < value) then
        low = middle + 1_i32
      else
        high = middle
      end if
    end do
    index = low
    if (index > 1_i32) then
      ! Preserve the previous linear scan's preference for the lower index on ties.
      if (abs(break_z(index - 1_i32) - value) <= abs(break_z(index) - value)) index = index - 1_i32
    end if
  end function nearest_break

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
