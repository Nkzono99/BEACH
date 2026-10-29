!> Resolve physical vacuum sides without changing ordered triangle winding.
module bem_panel_surface_sides
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
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
  public :: resolve_panel_surface_sides_by_mesh_id

  !> Bucket index of the half-edge start points used by outward_closed edge matching.
  type :: edge_bucket_index_type
    integer(i64), allocatable :: cell_key(:, :)
    integer(i32), allocatable :: table(:)
    integer(i32), allocatable :: cell_of(:)
    integer(i32), allocatable :: cell_start(:), cell_half_edge(:)
    integer(i32), allocatable :: near_start(:), near_cell(:)
    integer(i32) :: ncell = 0_i32
  end type edge_bucket_index_type

contains

  !> Apply one policy to every element, or to the elements of one mesh_id.
  subroutine resolve_panel_surface_sides(mesh, policy, status, message, mesh_id)
    type(mesh_type), intent(inout) :: mesh
    character(len=*), intent(in) :: policy
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer(i32), intent(in), optional :: mesh_id
    integer(i32), allocatable :: elements(:)
    integer(i32) :: i
    real(dp) :: coordinate_scale

    status = panel_surface_side_ok
    message = ''
    if (present(mesh_id)) then
      elements = pack([(i, i=1_i32, mesh%nelem)], mesh%elem_mesh_id == mesh_id)
    else
      elements = [(i, i=1_i32, mesh%nelem)]
    end if
    if (size(elements) == 0) then
      status = panel_surface_side_invalid_policy
      message = 'surface side mesh_id selects no elements'
      return
    end if
    coordinate_scale = 0.0_dp
    if (is_outward_closed(policy)) coordinate_scale = mesh_coordinate_scale(mesh)
    call resolve_element_sides(mesh, policy, elements, coordinate_scale, status, message)
  end subroutine resolve_panel_surface_sides

  !> Apply `policies(k)` to the elements whose mesh_id is k, for k = 1..size(policies), in order.
  !! Elements are grouped by mesh_id once, so each part costs only its own size. On failure the
  !! parts before `failed_mesh_id` stay resolved and the failing part is left unchanged.
  subroutine resolve_panel_surface_sides_by_mesh_id(mesh, policies, status, message, failed_mesh_id)
    type(mesh_type), intent(inout) :: mesh
    character(len=*), intent(in) :: policies(:)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer(i32), intent(out) :: failed_mesh_id
    integer(i32), allocatable :: part_start(:), part_fill(:), elements(:)
    integer(i32) :: i, k, nparts
    real(dp) :: coordinate_scale

    status = panel_surface_side_ok
    message = ''
    failed_mesh_id = 0_i32
    nparts = int(size(policies), i32)
    allocate (part_start(nparts + 1_i32), part_fill(nparts), elements(mesh%nelem))
    part_start = 0_i32
    do i = 1_i32, mesh%nelem
      k = mesh%elem_mesh_id(i)
      if (k >= 1_i32 .and. k <= nparts) part_start(k + 1_i32) = part_start(k + 1_i32) + 1_i32
    end do
    part_start(1) = 1_i32
    do k = 1_i32, nparts
      part_start(k + 1_i32) = part_start(k + 1_i32) + part_start(k)
    end do
    part_fill = part_start(1:nparts)
    do i = 1_i32, mesh%nelem
      k = mesh%elem_mesh_id(i)
      if (k < 1_i32 .or. k > nparts) cycle
      elements(part_fill(k)) = i
      part_fill(k) = part_fill(k) + 1_i32
    end do

    coordinate_scale = 0.0_dp
    do k = 1_i32, nparts
      if (is_outward_closed(policies(k))) then
        coordinate_scale = mesh_coordinate_scale(mesh)
        exit
      end if
    end do
    do k = 1_i32, nparts
      if (part_start(k + 1_i32) == part_start(k)) then
        status = panel_surface_side_invalid_policy
        message = 'surface side mesh_id selects no elements'
      else
        call resolve_element_sides( &
          mesh, policies(k), elements(part_start(k):part_start(k + 1_i32) - 1_i32), coordinate_scale, status, message &
          )
      end if
      if (status /= panel_surface_side_ok) then
        failed_mesh_id = k
        return
      end if
    end do
  end subroutine resolve_panel_surface_sides_by_mesh_id

  !> Resolve one policy for `elements` (ascending element indices). The mesh is changed only on success.
  !! `coordinate_scale` sets the outward_closed vertex tolerance and must cover the whole mesh.
  subroutine resolve_element_sides(mesh, policy, elements, coordinate_scale, status, message)
    type(mesh_type), intent(inout) :: mesh
    character(len=*), intent(in) :: policy
    integer(i32), intent(in) :: elements(:)
    real(dp), intent(in) :: coordinate_scale
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: normalized
    integer(i32), allocatable :: signs(:), neighbor(:, :), queue(:)
    logical, allocatable :: visited(:)
    integer(i32) :: i, m, edge, face, next_face, head, tail, component_sign, element
    real(dp) :: signed_volume, volume_scale, origin(3), a(3), b(3), c(3)

    status = panel_surface_side_ok
    message = ''
    normalized = lower_ascii(trim(policy))
    m = int(size(elements), i32)
    allocate (signs(m))
    select case (trim(normalized))
    case ('normal_plus')
      signs = 1_i32
    case ('normal_minus')
      signs = -1_i32
    case ('outward_closed')
      allocate (neighbor(3, m), queue(m), visited(m))
      call build_closed_neighbors(mesh, elements, coordinate_scale, neighbor, status, message)
      if (status /= panel_surface_side_ok) return
      visited = .false.
      do i = 1_i32, m
        if (visited(i)) cycle
        head = 1_i32
        tail = 1_i32
        queue(1) = i
        visited(i) = .true.
        ! Signed volume about a vertex of this component, so neither the world origin nor a
        ! 1 m unit enters: each tetrahedron term is a.(b x c)/6 with a, b, c relative to it,
        ! and its rounding error is bounded by a few eps*|a||b||c|/6 (length^3 scale).
        origin = mesh%v0(:, elements(i))
        signed_volume = 0.0_dp
        volume_scale = 0.0_dp
        do while (head <= tail)
          face = queue(head)
          head = head + 1_i32
          element = elements(face)
          a = mesh%v0(:, element) - origin
          b = mesh%v1(:, element) - origin
          c = mesh%v2(:, element) - origin
          signed_volume = signed_volume + dot_product(a, cross_product(b, c))/6.0_dp
          volume_scale = volume_scale + norm2(a)*norm2(b)*norm2(c)/6.0_dp
          do edge = 1, 3
            next_face = neighbor(edge, face)
            if (.not. visited(next_face)) then
              tail = tail + 1_i32
              queue(tail) = next_face
              visited(next_face) = .true.
            end if
          end do
        end do
        if (abs(signed_volume) <= 128.0_dp*epsilon(1.0_dp)*volume_scale) then
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

    do i = 1_i32, m
      element = elements(i)
      mesh%elem_vacuum_sign(element) = signs(i)
      mesh%vacuum_normals(:, element) = real(signs(i), dp)*mesh%normals(:, element)
    end do
  end subroutine resolve_element_sides

  !> Find, for every edge of `elements`, the one oppositely oriented edge of another element with
  !! the same mesh_id. Two vertices match when every coordinate differs by at most
  !! 128*eps*coordinate_scale. Edges are matched through a hash of start-point buckets four
  !! tolerances wide, so every point within tolerance lies in the same or an adjacent bucket and
  !! the exact pairwise test sees all candidates: an edge with no partner, two partners, or a
  !! same-direction partner is reported for the first such edge in element order, as a full search.
  subroutine build_closed_neighbors(mesh, elements, coordinate_scale, neighbor, status, message)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elements(:)
    real(dp), intent(in) :: coordinate_scale
    integer(i32), intent(out) :: neighbor(:, :)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(edge_bucket_index_type) :: buckets
    real(dp) :: a(3), b(3), c(3), d(3), tolerance
    integer(i32), allocatable :: stamp(:)
    integer(i32) :: m, face, edge, half_edge, side, bucket, near_bucket, near_pos, pos
    integer(i32) :: other_half_edge, other_face, other_edge
    integer(i32) :: element, other_element, matches, matched_face
    logical :: same_direction, opposite_direction

    neighbor = 0_i32
    status = panel_surface_side_ok
    message = ''
    m = int(size(elements), i32)
    do face = 1_i32, m
      element = elements(face)
      if (.not. (all(ieee_is_finite(mesh%v0(:, element))) .and. all(ieee_is_finite(mesh%v1(:, element))) .and. &
                 all(ieee_is_finite(mesh%v2(:, element))))) then
        status = panel_surface_side_invalid_policy
        message = 'outward_closed requires finite vertex coordinates'
        return
      end if
    end do
    tolerance = 128.0_dp*epsilon(1.0_dp)*coordinate_scale
    call build_edge_bucket_index(mesh, elements, 4.0_dp*tolerance, buckets)

    allocate (stamp(3_i32*m))
    stamp = 0_i32
    do face = 1_i32, m
      element = elements(face)
      do edge = 1_i32, 3_i32
        half_edge = 3_i32*(face - 1_i32) + edge
        call edge_vertices(mesh, element, edge, a, b)
        matches = 0_i32
        matched_face = 0_i32
        ! A same-direction partner starts near a; an opposite one starts near b, the start of the next edge.
        do side = 1, 2
          if (side == 1) then
            bucket = buckets%cell_of(half_edge)
          else
            bucket = buckets%cell_of(3_i32*(face - 1_i32) + modulo(edge, 3_i32) + 1_i32)
          end if
          do near_pos = buckets%near_start(bucket), buckets%near_start(bucket + 1_i32) - 1_i32
            near_bucket = buckets%near_cell(near_pos)
            do pos = buckets%cell_start(near_bucket), buckets%cell_start(near_bucket + 1_i32) - 1_i32
              other_half_edge = buckets%cell_half_edge(pos)
              if (stamp(other_half_edge) == half_edge) cycle
              stamp(other_half_edge) = half_edge
              other_face = (other_half_edge - 1_i32)/3_i32 + 1_i32
              if (other_face == face) cycle
              other_element = elements(other_face)
              if (mesh%elem_mesh_id(other_element) /= mesh%elem_mesh_id(element)) cycle
              other_edge = other_half_edge - 3_i32*(other_face - 1_i32)
              call edge_vertices(mesh, other_element, other_edge, c, d)
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

  !> Hash every half-edge start point of `elements` into cubic buckets of edge `bucket_size` and
  !! list, per bucket, its half-edges (3*(face-1)+edge, ascending) and its occupied neighbor
  !! buckets (itself included). A non-finite bucket size, which only an unbounded coordinate
  !! scale produces, puts every point into one bucket.
  subroutine build_edge_bucket_index(mesh, elements, bucket_size, buckets)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elements(:)
    real(dp), intent(in) :: bucket_size
    type(edge_bucket_index_type), intent(out) :: buckets
    integer(i32), allocatable :: fill(:), found(:)
    integer(i32) :: m, nhalf, half_edge, face, edge, bucket, table_size, slot, pass, nfound, total, dx, dy, dz
    integer(i64) :: key(3)
    real(dp) :: point(3)

    m = int(size(elements), i32)
    nhalf = 3_i32*m
    table_size = 16_i32
    do while (table_size < 2_i32*nhalf)
      table_size = 2_i32*table_size
    end do
    allocate (buckets%cell_key(3, nhalf), buckets%table(0:table_size - 1_i32), buckets%cell_of(nhalf))
    buckets%table = 0_i32
    buckets%ncell = 0_i32
    do face = 1_i32, m
      do edge = 1_i32, 3_i32
        half_edge = 3_i32*(face - 1_i32) + edge
        select case (edge)
        case (1)
          point = mesh%v0(:, elements(face))
        case (2)
          point = mesh%v1(:, elements(face))
        case default
          point = mesh%v2(:, elements(face))
        end select
        if (ieee_is_finite(bucket_size)) then
          key = floor(point/bucket_size, kind=i64)
        else
          key = 0_i64
        end if
        slot = hash_slot(key, table_size)
        do
          bucket = buckets%table(slot)
          if (bucket == 0_i32) then
            buckets%ncell = buckets%ncell + 1_i32
            bucket = buckets%ncell
            buckets%cell_key(:, bucket) = key
            buckets%table(slot) = bucket
            exit
          end if
          if (all(buckets%cell_key(:, bucket) == key)) exit
          slot = iand(slot + 1_i32, table_size - 1_i32)
        end do
        buckets%cell_of(half_edge) = bucket
      end do
    end do

    allocate (buckets%cell_start(buckets%ncell + 1_i32), buckets%cell_half_edge(nhalf), fill(buckets%ncell))
    buckets%cell_start = 0_i32
    do half_edge = 1_i32, nhalf
      bucket = buckets%cell_of(half_edge)
      buckets%cell_start(bucket + 1_i32) = buckets%cell_start(bucket + 1_i32) + 1_i32
    end do
    buckets%cell_start(1) = 1_i32
    do bucket = 1_i32, buckets%ncell
      buckets%cell_start(bucket + 1_i32) = buckets%cell_start(bucket + 1_i32) + buckets%cell_start(bucket)
    end do
    fill = buckets%cell_start(1:buckets%ncell)
    do half_edge = 1_i32, nhalf
      bucket = buckets%cell_of(half_edge)
      buckets%cell_half_edge(fill(bucket)) = half_edge
      fill(bucket) = fill(bucket) + 1_i32
    end do

    ! Two passes over the 27 neighbor keys: count, then fill the per-bucket neighbor list.
    allocate (buckets%near_start(buckets%ncell + 1_i32), found(27))
    buckets%near_start(1) = 1_i32
    total = 0_i32
    do pass = 1, 2
      if (pass == 2) allocate (buckets%near_cell(total))
      do bucket = 1_i32, buckets%ncell
        nfound = 0_i32
        do dz = -1_i32, 1_i32
          do dy = -1_i32, 1_i32
            do dx = -1_i32, 1_i32
              key = buckets%cell_key(:, bucket) + [int(dx, i64), int(dy, i64), int(dz, i64)]
              slot = hash_slot(key, table_size)
              do while (buckets%table(slot) /= 0_i32)
                if (all(buckets%cell_key(:, buckets%table(slot)) == key)) then
                  nfound = nfound + 1_i32
                  found(nfound) = buckets%table(slot)
                  exit
                end if
                slot = iand(slot + 1_i32, table_size - 1_i32)
              end do
            end do
          end do
        end do
        if (pass == 1) then
          buckets%near_start(bucket + 1_i32) = buckets%near_start(bucket) + nfound
        else
          buckets%near_cell(buckets%near_start(bucket):buckets%near_start(bucket + 1_i32) - 1_i32) = found(1:nfound)
        end if
      end do
      total = buckets%near_start(buckets%ncell + 1_i32) - 1_i32
    end do
  end subroutine build_edge_bucket_index

  !> Open-addressing slot of a bucket key; each component is reduced first so no product overflows.
  pure integer(i32) function hash_slot(key, table_size) result(slot)
    integer(i64), intent(in) :: key(3)
    integer(i32), intent(in) :: table_size
    integer(i64), parameter :: prime = 1000003_i64
    integer(i64) :: mixed

    mixed = modulo(key(1), prime)*73856093_i64 + modulo(key(2), prime)*19349663_i64 + &
            modulo(key(3), prime)*83492791_i64
    slot = int(iand(mixed, int(table_size - 1_i32, i64)), i32)
  end function hash_slot

  logical function is_outward_closed(policy)
    character(len=*), intent(in) :: policy
    character(len=32) :: normalized

    normalized = lower_ascii(trim(policy))
    is_outward_closed = trim(normalized) == 'outward_closed'
  end function is_outward_closed

  !> Largest vertex coordinate magnitude of the whole mesh, at least 1; scales the vertex tolerance.
  pure real(dp) function mesh_coordinate_scale(mesh) result(scale)
    type(mesh_type), intent(in) :: mesh

    scale = max(1.0_dp, max(maxval(abs(mesh%v0)), max(maxval(abs(mesh%v1)), maxval(abs(mesh%v2)))))
  end function mesh_coordinate_scale

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
