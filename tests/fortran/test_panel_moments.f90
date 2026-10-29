!> Triangle geometry, cached moments, cubature, and surface-side invariants.
program test_panel_moments
  use bem_kinds, only: dp, i32
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_quadrature, only: panel_quadrature_plan_type, build_panel_quadrature, build_panel_duffy_quadrature
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, resolve_panel_surface_sides_by_mesh_id, &
                                     panel_surface_side_ok, panel_surface_side_open, panel_surface_side_inconsistent, &
                                     panel_surface_side_nonmanifold, panel_surface_side_invalid_policy
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          assert_close_dp, assert_allclose_1d
  implicit none

  type(panel_geometry_type) :: geom, cyclic, reversed, scaled
  type(panel_quadrature_plan_type) :: quad
  type(mesh_type) :: mesh
  real(dp) :: v0(3), v1(3), v2(3), cubature_value
  real(dp) :: tetra_v0(3, 4), tetra_v1(3, 4), tetra_v2(3, 4)
  real(dp), parameter :: unit_moment1(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.0_dp]
  real(dp), parameter :: unit_moment2(9) = [ &
                         1.0_dp/12.0_dp, 1.0_dp/24.0_dp, 0.0_dp, &
                         1.0_dp/24.0_dp, 1.0_dp/12.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp &
                         ]
  character(len=128) :: message
  integer(i32) :: status
  integer :: x_power, y_power

  v0 = [0.0_dp, 0.0_dp, 0.0_dp]
  v1 = [1.0_dp, 0.0_dp, 0.0_dp]
  v2 = [0.0_dp, 1.0_dp, 0.0_dp]
  call test_init(11)

  call test_begin('unit_triangle_geometry_and_raw_moments')
  call init_panel_geometry(v0, v1, v2, geom, status)
  call assert_equal_i32(status, panel_geometry_ok, 'valid triangle status mismatch')
  call assert_close_dp(geom%area, 0.5_dp, 1.0e-15_dp, 'triangle area mismatch')
  call assert_allclose_1d(geom%centroid, [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp], 1.0e-15_dp, 'centroid mismatch')
  call assert_allclose_1d(geom%normal, [0.0_dp, 0.0_dp, 1.0_dp], 1.0e-15_dp, 'normal mismatch')
  call assert_close_dp(geom%moment0, 0.5_dp, 1.0e-15_dp, 'zeroth moment mismatch')
  call assert_allclose_1d(geom%moment1, unit_moment1, 1.0e-15_dp, 'first moment mismatch')
  call assert_allclose_1d(reshape(geom%moment2, [9]), unit_moment2, 1.0e-15_dp, 'second moment mismatch')
  call test_end()

  call test_begin('mesh_owns_panel_geometry')
  call init_mesh(mesh, reshape(v0, [3, 1]), reshape(v1, [3, 1]), reshape(v2, [3, 1]))
  call assert_close_dp(mesh%panel_area(1), 0.5_dp, 1.0e-15_dp, 'mesh panel area mismatch')
  call assert_allclose_1d(mesh%panel_moment1(:, 1), unit_moment1, 1.0e-15_dp, &
                          'mesh panel moment mismatch')
  call assert_allclose_1d(reshape(mesh%panel_moment2(:, :, 1), [9]), unit_moment2, 1.0e-15_dp, &
                          'mesh second moment mismatch')
  call assert_allclose_1d(mesh%panel_edge_length(:, 1), [1.0_dp, sqrt(2.0_dp), 1.0_dp], 1.0e-15_dp, &
                          'mesh edge cache mismatch')
  call assert_close_dp(sum(mesh%panel_quad_weight(:, 1)), 0.5_dp, 1.0e-15_dp, 'mesh quadrature weight mismatch')
  call test_end()

  call test_begin('open_surface_explicit_side')
  call init_mesh(mesh, reshape(v0, [3, 1]), reshape(v1, [3, 1]), reshape(v2, [3, 1]))
  call resolve_panel_surface_sides(mesh, 'normal_plus', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'normal_plus side resolution failed')
  call assert_equal_i32(mesh%elem_vacuum_sign(1), 1_i32, 'vacuum sign mismatch')
  call assert_allclose_1d(mesh%vacuum_normals(:, 1), [0.0_dp, 0.0_dp, 1.0_dp], 1.0e-15_dp, &
                          'vacuum normal mismatch')
  call resolve_panel_surface_sides(mesh, 'outward_closed', status, message)
  call assert_equal_i32(status, panel_surface_side_open, 'open surface must reject outward_closed')
  call test_end()

  call test_begin('closed_surface_outward_side_is_winding_independent')
  tetra_v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  tetra_v2(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
  tetra_v0(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
  tetra_v2(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
  tetra_v0(:, 3) = [0.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
  tetra_v2(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
  tetra_v0(:, 4) = [1.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 4) = [0.0_dp, 1.0_dp, 0.0_dp]
  tetra_v2(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
  call init_mesh(mesh, tetra_v0, tetra_v1, tetra_v2)
  call resolve_panel_surface_sides(mesh, 'outward_closed', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'closed tetrahedron side resolution failed')
  call assert_true(all(mesh%elem_vacuum_sign == 1_i32), 'outward tetrahedron signs must be positive')
  call assert_true(all(sum(mesh%vacuum_normals*(mesh%centers - 0.25_dp), dim=1) > 0.0_dp), &
                   'vacuum normals must face away from tetrahedron interior')
  call init_mesh(mesh, tetra_v0, tetra_v2, tetra_v1)
  call resolve_panel_surface_sides(mesh, 'outward_closed', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'inward-wound tetrahedron side resolution failed')
  call assert_true(all(mesh%elem_vacuum_sign == -1_i32), 'inward tetrahedron signs must be negative')
  call assert_true(all(sum(mesh%vacuum_normals*(mesh%centers - 0.25_dp), dim=1) > 0.0_dp), &
                   'inward-wound tetrahedron vacuum normals must face outward')
  call test_end()

  call test_begin('closed_surface_side_is_scale_and_translation_invariant')
  call check_closed_surface_scaling(tetra_v0, tetra_v1, tetra_v2)
  call test_end()

  call test_begin('closed_surface_edge_matching_contract')
  call check_closed_surface_edge_matching(tetra_v0, tetra_v1, tetra_v2)
  call test_end()

  call test_begin('surface_sides_by_mesh_id')
  call check_surface_sides_by_mesh_id(tetra_v0, tetra_v1, tetra_v2)
  call test_end()

  call test_begin('degree5_cubature_integrates_all_monomials')
  call build_panel_quadrature(geom, quad)
  do x_power = 0, 5
    do y_power = 0, 5 - x_power
      cubature_value = sum(quad%weight*quad%position(1, :)**x_power*quad%position(2, :)**y_power)
      call assert_close_dp(cubature_value, unit_triangle_monomial_integral(x_power, y_power), 1.0e-14_dp, &
                           'degree-5 cubature monomial mismatch')
    end do
  end do
  call test_end()

  call test_begin('moment_permutation_and_scaling_invariants')
  call init_panel_geometry(v1, v2, v0, cyclic, status)
  call assert_equal_i32(status, panel_geometry_ok, 'cyclic triangle status mismatch')
  call assert_close_dp(cyclic%area, geom%area, 1.0e-15_dp, 'cyclic area mismatch')
  call assert_allclose_1d(cyclic%moment1, geom%moment1, 1.0e-15_dp, 'cyclic first moment mismatch')
  call assert_allclose_1d(reshape(cyclic%moment2, [9]), reshape(geom%moment2, [9]), 1.0e-15_dp, 'cyclic moments mismatch')
  call assert_allclose_1d(cyclic%normal, geom%normal, 1.0e-15_dp, 'cyclic normal mismatch')
  call init_panel_geometry(v0, v2, v1, reversed, status)
  call assert_equal_i32(status, panel_geometry_ok, 'reversed triangle status mismatch')
  call assert_allclose_1d(reversed%normal, -geom%normal, 1.0e-15_dp, 'reversed normal mismatch')
  call assert_allclose_1d(reversed%moment1, geom%moment1, 1.0e-15_dp, 'reversed first moment mismatch')
  call assert_allclose_1d(reshape(reversed%moment2, [9]), reshape(geom%moment2, [9]), 1.0e-15_dp, 'reversed moments mismatch')
  call init_panel_geometry(2.0_dp*v0, 2.0_dp*v1, 2.0_dp*v2, scaled, status)
  call assert_equal_i32(status, panel_geometry_ok, 'scaled triangle status mismatch')
  call assert_close_dp(scaled%moment0, 4.0_dp*geom%moment0, 1.0e-14_dp, 'scaled zeroth moment mismatch')
  call assert_allclose_1d(scaled%moment1, 8.0_dp*geom%moment1, 1.0e-14_dp, 'scaled first moment mismatch')
  call assert_allclose_1d(reshape(scaled%moment2, [9]), 16.0_dp*reshape(geom%moment2, [9]), 1.0e-14_dp, &
                          'scaled second moment mismatch')
  call test_end()

  call test_begin('reused_cubature_replaces_geometry_and_rule')
  call build_panel_duffy_quadrature(geom, 3_i32, quad)
  call build_panel_quadrature(scaled, quad)
  call assert_equal_i32(quad%npoint, 7_i32, 'cubature must replace the prior Duffy rule')
  call assert_true(all(shape(quad%position) == [3, 7]) .and. size(quad%weight) == 7, 'cubature resize mismatch')
  call assert_close_dp(sum(quad%weight), 2.0_dp, 1.0e-14_dp, 'reused cubature area mismatch')
  call assert_close_dp(sum(quad%weight*quad%position(1, :)), 4.0_dp/3.0_dp, 1.0e-14_dp, &
                       'reused cubature first moment mismatch')
  call build_panel_quadrature(geom, quad)
  call assert_close_dp(sum(quad%weight), 0.5_dp, 1.0e-14_dp, 'same-size rebuild retained previous area')
  call assert_close_dp(sum(quad%weight*quad%position(1, :)), 1.0_dp/6.0_dp, 1.0e-14_dp, &
                       'same-size rebuild retained previous positions')
  call test_end()

  call test_begin('thin_triangle_retains_collision_geometry')
  call init_mesh(mesh, reshape(v0, [3, 1]), reshape(v1, [3, 1]), reshape([0.0_dp, 1.0e-16_dp, 0.0_dp], [3, 1]))
  call assert_close_dp(mesh%panel_area(1), 0.0_dp, 0.0_dp, 'thin triangle must remain excluded from panel integration')
  call assert_close_dp(mesh%h_elem(1), sqrt(0.5e-16_dp), 1.0e-23_dp, 'thin collision triangle lost its size')
  call assert_allclose_1d(mesh%normals(:, 1), [0.0_dp, 0.0_dp, 1.0_dp], 0.0_dp, 'thin triangle normal mismatch')
  call assert_allclose_1d(mesh%centers(:, 1), [1.0_dp, 1.0e-16_dp, 0.0_dp]/3.0_dp, 0.0_dp, &
                          'thin triangle centroid mismatch')
  call assert_true(all(mesh%panel_quad_weight == 0.0_dp), 'thin triangle retained quadrature weights')
  call test_end()

  call test_summary()

contains

  !> The same closed tetrahedron at 1 m, 1 mm, 1 um and 1 nm, near the origin and translated
  !! as far as binary64 still resolves its edges, must give the same side resolution: outward
  !! winding +1, inward winding -1, and a doubly covered (zero-volume) closed sheet rejected.
  subroutine check_closed_surface_scaling(unit_v0, unit_v1, unit_v2)
    real(dp), intent(in) :: unit_v0(3, 4), unit_v1(3, 4), unit_v2(3, 4)
    real(dp), parameter :: scales(4) = [1.0_dp, 1.0e-3_dp, 1.0e-6_dp, 1.0e-9_dp]
    real(dp), parameter :: offsets(3) = [0.0_dp, 0.75_dp, -1.0e3_dp]
    type(mesh_type) :: scaled_mesh
    real(dp) :: shift(3), p(3, 4), sheet_v0(3, 2), sheet_v1(3, 2), sheet_v2(3, 2)
    integer :: i_scale, i_offset, vertex
    character(len=128) :: side_message
    integer(i32) :: side_status

    do i_scale = 1, size(scales)
      do i_offset = 1, size(offsets)
        shift = offsets(i_offset)*[1.0_dp, -0.5_dp, 0.25_dp]
        ! Shared vertices are built once so that coincident corners stay bitwise equal.
        p(:, 1) = shift
        do vertex = 2, 4
          p(:, vertex) = shift
          p(vertex - 1, vertex) = shift(vertex - 1) + scales(i_scale)
        end do
        call init_mesh(scaled_mesh, map_corners(unit_v0, p), map_corners(unit_v1, p), map_corners(unit_v2, p))
        call resolve_panel_surface_sides(scaled_mesh, 'outward_closed', side_status, side_message)
        call assert_equal_i32(side_status, panel_surface_side_ok, 'scaled outward tetrahedron rejected')
        call assert_true(all(scaled_mesh%elem_vacuum_sign == 1_i32), 'scaled outward tetrahedron sign')
        call init_mesh(scaled_mesh, map_corners(unit_v0, p), map_corners(unit_v2, p), map_corners(unit_v1, p))
        call resolve_panel_surface_sides(scaled_mesh, 'outward_closed', side_status, side_message)
        call assert_equal_i32(side_status, panel_surface_side_ok, 'scaled inward tetrahedron rejected')
        call assert_true(all(scaled_mesh%elem_vacuum_sign == -1_i32), 'scaled inward tetrahedron sign')
        sheet_v0(:, 1) = p(:, 1)
        sheet_v1(:, 1) = p(:, 2)
        sheet_v2(:, 1) = p(:, 3)
        sheet_v0(:, 2) = p(:, 1)
        sheet_v1(:, 2) = p(:, 3)
        sheet_v2(:, 2) = p(:, 2)
        call init_mesh(scaled_mesh, sheet_v0, sheet_v1, sheet_v2)
        call resolve_panel_surface_sides(scaled_mesh, 'outward_closed', side_status, side_message)
        call assert_equal_i32(side_status, panel_surface_side_inconsistent, 'zero-volume closed sheet accepted')
      end do
    end do
  end subroutine check_closed_surface_scaling

  !> outward_closed edge matching on the unit tetrahedron: corners duplicated within the vertex
  !! tolerance still close the surface and corners further apart open it; a face flipped against
  !! its neighbors is inconsistent; a back-facing copy of a face gives its edges two partners; two
  !! disjoint components in one part get their own outward sides.
  subroutine check_closed_surface_edge_matching(t0, t1, t2)
    real(dp), intent(in) :: t0(3, 4), t1(3, 4), t2(3, 4)
    real(dp) :: v0(3, 8), v1(3, 8), v2(3, 8), center(3, 8)
    type(mesh_type) :: closed
    character(len=128) :: side_message
    integer(i32) :: side_status

    ! The vertex tolerance at unit scale is 128 eps (about 2.8e-14) per coordinate.
    v0(:, 1:4) = t0 + 1.0e-15_dp*reshape([1, -2, 3, 0, 2, -1, -3, 1, 2, 1, 0, -2], [3, 4])
    v1(:, 1:4) = t1 - 2.0e-15_dp*reshape([2, 1, -1, 3, 0, 1, 1, -2, 0, -1, 2, 1], [3, 4])
    v2(:, 1:4) = t2
    call init_mesh(closed, v0(:, 1:4), v1(:, 1:4), v2(:, 1:4))
    call resolve_panel_surface_sides(closed, 'outward_closed', side_status, side_message)
    call assert_equal_i32(side_status, panel_surface_side_ok, 'corners duplicated within tolerance must close')
    call assert_true(all(closed%elem_vacuum_sign == 1_i32), 'duplicated-corner tetrahedron sign')
    v1(1, 2) = v1(1, 2) + 1.0e-12_dp
    call init_mesh(closed, v0(:, 1:4), v1(:, 1:4), v2(:, 1:4))
    call resolve_panel_surface_sides(closed, 'outward_closed', side_status, side_message)
    call assert_equal_i32(side_status, panel_surface_side_open, 'corner beyond tolerance must open the surface')
    call assert_true(all(closed%elem_vacuum_sign == 0_i32), 'rejected surface must keep its unresolved sides')

    v1(:, 1:4) = t1
    v2(:, 1:4) = t2
    v1(:, 3) = t2(:, 3)
    v2(:, 3) = t1(:, 3)
    call init_mesh(closed, t0, v1(:, 1:4), v2(:, 1:4))
    call resolve_panel_surface_sides(closed, 'outward_closed', side_status, side_message)
    call assert_equal_i32(side_status, panel_surface_side_inconsistent, 'flipped face must be inconsistent')

    v0(:, 1:4) = t0
    v1(:, 1:4) = t1
    v2(:, 1:4) = t2
    v0(:, 5) = t0(:, 1)
    v1(:, 5) = t2(:, 1)
    v2(:, 5) = t1(:, 1)
    call init_mesh(closed, v0(:, 1:5), v1(:, 1:5), v2(:, 1:5))
    call resolve_panel_surface_sides(closed, 'outward_closed', side_status, side_message)
    call assert_equal_i32(side_status, panel_surface_side_nonmanifold, 'three faces on one edge must be nonmanifold')

    ! Second tetrahedron: shifted by 3 along x and wound inward.
    v0(:, 1:4) = t0
    v1(:, 1:4) = t1
    v2(:, 1:4) = t2
    v0(:, 5:8) = t0 + spread([3.0_dp, 0.0_dp, 0.0_dp], 2, 4)
    v1(:, 5:8) = t2 + spread([3.0_dp, 0.0_dp, 0.0_dp], 2, 4)
    v2(:, 5:8) = t1 + spread([3.0_dp, 0.0_dp, 0.0_dp], 2, 4)
    call init_mesh(closed, v0, v1, v2)
    call resolve_panel_surface_sides(closed, 'outward_closed', side_status, side_message)
    call assert_equal_i32(side_status, panel_surface_side_ok, 'two closed components in one part')
    call assert_true(all(closed%elem_vacuum_sign == [1, 1, 1, 1, -1, -1, -1, -1]), 'per-component signs')
    center(:, 1:4) = 0.25_dp
    center(:, 5:8) = spread([3.25_dp, 0.25_dp, 0.25_dp], 2, 4)
    call assert_true(all(sum(closed%vacuum_normals*(closed%centers - center), dim=1) > 0.0_dp), &
                     'each component must face away from its own interior')
  end subroutine check_closed_surface_edge_matching

  !> One policy per mesh_id, applied in id order to elements that need not be contiguous: a
  !! coincident tetrahedron in another id never pairs with this one, and a failing id is reported
  !! with the earlier ids resolved and its own elements left unchanged.
  subroutine check_surface_sides_by_mesh_id(t0, t1, t2)
    real(dp), intent(in) :: t0(3, 4), t1(3, 4), t2(3, 4)
    type(mesh_type) :: parts
    character(len=128) :: side_message
    integer(i32) :: side_status, failed_mesh_id
    integer(i32), parameter :: ids(9) = [1, 2, 1, 2, 3, 1, 2, 1, 2]
    integer, parameter :: face_of(9) = [1, 1, 2, 2, 1, 3, 3, 4, 4]
    real(dp) :: v0(3, 9), v1(3, 9), v2(3, 9)
    integer :: element

    ! ids 1 and 2 hold the same tetrahedron (id 2 inward wound), interleaved; id 3 is one triangle.
    do element = 1, 9
      v0(:, element) = t0(:, face_of(element))
      v1(:, element) = t1(:, face_of(element))
      v2(:, element) = t2(:, face_of(element))
      if (ids(element) == 2_i32) then
        v1(:, element) = t2(:, face_of(element))
        v2(:, element) = t1(:, face_of(element))
      end if
    end do
    call init_mesh(parts, v0, v1, v2, elem_mesh_id0=ids)
    call resolve_panel_surface_sides_by_mesh_id( &
      parts, [character(len=16) :: 'outward_closed', 'OUTWARD_CLOSED', 'normal_minus'], side_status, side_message, &
      failed_mesh_id)
    call assert_equal_i32(side_status, panel_surface_side_ok, 'interleaved parts must resolve')
    call assert_equal_i32(failed_mesh_id, 0_i32, 'no failed mesh_id on success')
    call assert_true(all(parts%elem_vacuum_sign == [1, -1, 1, -1, -1, 1, -1, 1, -1]), 'per-part signs')
    call assert_true(all(parts%vacuum_normals(:, 5) == -parts%normals(:, 5)), 'normal_minus vacuum normal')

    call init_mesh(parts, v0, v1, v2, elem_mesh_id0=ids)
    call resolve_panel_surface_sides_by_mesh_id( &
      parts, [character(len=16) :: 'normal_plus', 'normal_minus', 'outward_closed'], side_status, side_message, &
      failed_mesh_id)
    call assert_equal_i32(side_status, panel_surface_side_open, 'open part must fail')
    call assert_equal_i32(failed_mesh_id, 3_i32, 'failing mesh_id')
    call assert_true(all(parts%elem_vacuum_sign(1:4) == [1, -1, 1, -1]), 'earlier parts stay resolved')
    call assert_equal_i32(parts%elem_vacuum_sign(5), 0_i32, 'failing part stays unresolved')

    call init_mesh(parts, v0, v1, v2, elem_mesh_id0=ids)
    call resolve_panel_surface_sides_by_mesh_id( &
      parts, [character(len=16) :: 'normal_plus', 'normal_plus', 'normal_plus', 'normal_plus'], side_status, &
      side_message, failed_mesh_id)
    call assert_equal_i32(side_status, panel_surface_side_invalid_policy, 'mesh_id without elements')
    call assert_equal_i32(failed_mesh_id, 4_i32, 'empty mesh_id reported')
  end subroutine check_surface_sides_by_mesh_id

  !> Replace each unit-tetrahedron corner (origin or a unit axis point) by its scaled, shifted copy.
  pure function map_corners(unit_vertices, corners) result(mapped)
    real(dp), intent(in) :: unit_vertices(3, 4), corners(3, 4)
    real(dp) :: mapped(3, 4)
    integer :: face

    do face = 1, 4
      if (all(unit_vertices(:, face) == 0.0_dp)) then
        mapped(:, face) = corners(:, 1)
      else
        mapped(:, face) = corners(:, 1 + maxloc(unit_vertices(:, face), dim=1))
      end if
    end do
  end function map_corners

  pure real(dp) function unit_triangle_monomial_integral(px, py) result(integral)
    integer, intent(in) :: px, py

    integral = real(factorial(px)*factorial(py), dp)/real(factorial(px + py + 2), dp)
  end function unit_triangle_monomial_integral

  pure integer function factorial(n) result(value)
    integer, intent(in) :: n
    integer :: factor

    value = 1
    do factor = 2, n
      value = value*factor
    end do
  end function factorial

end program test_panel_moments
