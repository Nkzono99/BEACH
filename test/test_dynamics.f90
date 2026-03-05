!> 電場評価・Boris更新・衝突判定の基本動作を検証するテスト。
program test_dynamics
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, hit_info
  use bem_mesh, only: init_mesh
  use bem_field, only: electric_field_at
  use bem_pusher, only: boris_push
  use bem_collision, only: find_first_hit
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  type(mesh_type) :: mesh_field, mesh_hit
  type(hit_info) :: hit
  real(dp) :: v0_field(3, 1), v1_field(3, 1), v2_field(3, 1), q0_field(1)
  real(dp) :: v0_hit(3, 2), v1_hit(3, 2), v2_hit(3, 2)
  real(dp) :: e(3), x_new(3), v_new(3), speed0, speed1
  real(dp) :: inv_r3, expected_ex

  v0_field(:, 1) = [1.0d0, 0.0d0, 0.0d0]
  v1_field(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  v2_field(:, 1) = [-1.0d0, -1.0d0, 0.0d0]
  q0_field(1) = 2.0d-9
  call init_mesh(mesh_field, v0_field, v1_field, v2_field, q0=q0_field)

  call electric_field_at(mesh_field, [1.0d0, 0.0d0, 0.0d0], 0.5d0, e)
  inv_r3 = 1.0d0 / (sqrt(1.25d0) * 1.25d0)
  expected_ex = k_coulomb * q0_field(1) * inv_r3
  call assert_close_dp(e(1), expected_ex, abs(expected_ex) * 1.0d-12, 'electric field Ex mismatch')
  call assert_close_dp(e(2), 0.0d0, 1.0d-20, 'electric field Ey should be zero')
  call assert_close_dp(e(3), 0.0d0, 1.0d-20, 'electric field Ez should be zero')

  call boris_push( &
    [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], &
    2.0d0, 1.0d0, 0.5d0, [1.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], &
    x_new, v_new &
  )
  call assert_allclose_1d(v_new, [1.0d0, 0.0d0, 0.0d0], 1.0d-12, 'boris (E only) velocity mismatch')
  call assert_allclose_1d(x_new, [0.5d0, 0.0d0, 0.0d0], 1.0d-12, 'boris (E only) position mismatch')

  call boris_push( &
    [0.0d0, 0.0d0, 0.0d0], [1.0d0, 2.0d0, -0.5d0], &
    1.0d0, 1.0d0, 0.1d0, [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 2.0d0], &
    x_new, v_new &
  )
  speed0 = sqrt(1.0d0 + 4.0d0 + 0.25d0)
  speed1 = sqrt(sum(v_new * v_new))
  call assert_close_dp(speed1, speed0, 1.0d-12, 'boris should preserve speed when E=0')

  v0_hit(:, 1) = [0.0d0, 0.0d0, 1.0d0]
  v1_hit(:, 1) = [1.0d0, 0.0d0, 1.0d0]
  v2_hit(:, 1) = [0.0d0, 1.0d0, 1.0d0]
  v0_hit(:, 2) = [0.0d0, 0.0d0, 0.0d0]
  v1_hit(:, 2) = [1.0d0, 0.0d0, 0.0d0]
  v2_hit(:, 2) = [0.0d0, 1.0d0, 0.0d0]
  call init_mesh(mesh_hit, v0_hit, v1_hit, v2_hit)

  call find_first_hit(mesh_hit, [0.2d0, 0.2d0, 2.0d0], [0.2d0, 0.2d0, -1.0d0], hit)
  call assert_true(hit%has_hit, 'segment should hit the mesh')
  call assert_equal_i32(hit%elem_idx, 1_i32, 'first hit should be the nearer triangle')
  call assert_close_dp(hit%t, 1.0d0 / 3.0d0, 1.0d-12, 'hit parameter mismatch')
  call assert_close_dp(hit%pos(3), 1.0d0, 1.0d-12, 'hit position mismatch')

  call find_first_hit(mesh_hit, [0.9d0, 0.9d0, 2.0d0], [0.9d0, 0.9d0, -1.0d0], hit)
  call assert_true(.not. hit%has_hit, 'segment outside triangle should not hit')
end program test_dynamics
