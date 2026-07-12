!> Resolve physical vacuum sides without changing ordered triangle winding.
module bem_panel_surface_sides
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: panel_surface_side_ok = 0_i32
  integer(i32), parameter, public :: panel_surface_side_open = 1_i32
  integer(i32), parameter, public :: panel_surface_side_inconsistent = 2_i32
  integer(i32), parameter, public :: panel_surface_side_nonmanifold = 3_i32
  integer(i32), parameter, public :: panel_surface_side_invalid_policy = 4_i32

  public :: resolve_panel_surface_sides

contains

  subroutine resolve_panel_surface_sides(mesh, policy, status, message, mesh_id)
    type(mesh_type), intent(inout) :: mesh
    character(len=*), intent(in) :: policy
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer(i32), intent(in), optional :: mesh_id
    character(len=32) :: normalized
    integer(i32), allocatable :: signs(:), neighbor(:, :), queue(:)
    logical, allocatable :: visited(:), active(:)
    integer(i32) :: i, edge, face, next_face, head, tail, component_sign
    real(dp) :: signed_volume, volume_scale

    status = panel_surface_side_ok
    message = ''
    normalized = lower_ascii(trim(policy))
    allocate (signs(mesh%nelem), active(mesh%nelem))
    active = .true.
    if (present(mesh_id)) active = mesh%elem_mesh_id == mesh_id
    if (.not. any(active)) then
      status = panel_surface_side_invalid_policy
      message = 'surface side mesh_id selects no elements'
      return
    end if
    signs = mesh%elem_vacuum_sign
    select case (trim(normalized))
    case ('normal_plus')
      where (active) signs = 1_i32
    case ('normal_minus')
      where (active) signs = -1_i32
    case ('outward_closed')
      allocate (neighbor(3, mesh%nelem), queue(mesh%nelem), visited(mesh%nelem))
      call build_closed_neighbors(mesh, active, neighbor, status, message)
      if (status /= panel_surface_side_ok) return
      where (active) signs = 0_i32
      visited = .false.
      do i = 1, mesh%nelem
        if (.not. active(i)) cycle
        if (visited(i)) cycle
        head = 1_i32
        tail = 1_i32
        queue(1) = i
        visited(i) = .true.
        signed_volume = 0.0_dp
        volume_scale = 0.0_dp
        do while (head <= tail)
          face = queue(head)
          head = head + 1_i32
          signed_volume = signed_volume + dot_product( &
                          mesh%v0(:, face), cross_product(mesh%v1(:, face), mesh%v2(:, face)) &
                          )/6.0_dp
          volume_scale = volume_scale + mesh%panel_area(face)*max(1.0_dp, norm2(mesh%centers(:, face)))
          do edge = 1, 3
            next_face = neighbor(edge, face)
            if (.not. visited(next_face)) then
              tail = tail + 1_i32
              queue(tail) = next_face
              visited(next_face) = .true.
            end if
          end do
        end do
        if (abs(signed_volume) <= 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, volume_scale)) then
          status = panel_surface_side_inconsistent
          message = 'closed component has zero or indeterminate signed volume'
          return
        end if
        component_sign = merge(1_i32, -1_i32, signed_volume > 0.0_dp)
        do face = 1, tail
          signs(queue(face)) = component_sign
        end do
      end do
    case default
      status = panel_surface_side_invalid_policy
      message = 'surface side policy must be normal_plus, normal_minus, or outward_closed'
      return
    end select

    mesh%elem_vacuum_sign = signs
    do i = 1, mesh%nelem
      if (signs(i) /= 0_i32) mesh%vacuum_normals(:, i) = real(signs(i), dp)*mesh%normals(:, i)
    end do
  end subroutine resolve_panel_surface_sides

  subroutine build_closed_neighbors(mesh, active, neighbor, status, message)
    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: active(mesh%nelem)
    integer(i32), intent(out) :: neighbor(3, mesh%nelem)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp) :: a(3), b(3), c(3), d(3), tolerance, coordinate_scale
    integer(i32) :: face, edge, other_face, other_edge, matches, matched_face
    logical :: same_direction, opposite_direction

    neighbor = 0_i32
    status = panel_surface_side_ok
    message = ''
    coordinate_scale = max(1.0_dp, max(maxval(abs(mesh%v0)), max(maxval(abs(mesh%v1)), maxval(abs(mesh%v2)))))
    tolerance = 128.0_dp*epsilon(1.0_dp)*coordinate_scale
    do face = 1, mesh%nelem
      if (.not. active(face)) cycle
      do edge = 1, 3
        call edge_vertices(mesh, face, edge, a, b)
        matches = 0_i32
        matched_face = 0_i32
        do other_face = 1, mesh%nelem
          if (.not. active(other_face)) cycle
          if (other_face == face) cycle
          if (mesh%elem_mesh_id(other_face) /= mesh%elem_mesh_id(face)) cycle
          do other_edge = 1, 3
            call edge_vertices(mesh, other_face, other_edge, c, d)
            same_direction = same_point(a, c, tolerance) .and. same_point(b, d, tolerance)
            opposite_direction = same_point(a, d, tolerance) .and. same_point(b, c, tolerance)
            if (same_direction) then
              status = panel_surface_side_inconsistent
              message = 'closed surface contains inconsistent face winding'
              return
            end if
            if (opposite_direction) then
              matches = matches + 1_i32
              matched_face = other_face
            end if
          end do
        end do
        if (matches == 0_i32) then
          status = panel_surface_side_open
          message = 'outward_closed requires every edge to have a partner'
          return
        else if (matches > 1_i32) then
          status = panel_surface_side_nonmanifold
          message = 'outward_closed requires a two-manifold mesh'
          return
        end if
        neighbor(edge, face) = matched_face
      end do
    end do
  end subroutine build_closed_neighbors

  pure subroutine edge_vertices(mesh, face, edge, a, b)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: face, edge
    real(dp), intent(out) :: a(3), b(3)

    select case (edge)
    case (1)
      a = mesh%v0(:, face)
      b = mesh%v1(:, face)
    case (2)
      a = mesh%v1(:, face)
      b = mesh%v2(:, face)
    case default
      a = mesh%v2(:, face)
      b = mesh%v0(:, face)
    end select
  end subroutine edge_vertices

  pure logical function same_point(a, b, tolerance) result(same)
    real(dp), intent(in) :: a(3), b(3), tolerance

    same = maxval(abs(a - b)) <= tolerance
  end function same_point

  pure real(dp) function norm2(vector) result(norm)
    real(dp), intent(in) :: vector(3)

    norm = sqrt(sum(vector*vector))
  end function norm2

  pure function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross_product

end module bem_panel_surface_sides
