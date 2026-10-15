!> 照射セルの設定解決と壁の幾何学的問い合わせを検証する。
program test_irradiation_cell_contract
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_periodic, bc_open, bc_reflect
  use bem_app_config_types, only: app_config, default_app_config
  use bem_irradiation_cell_contract, only: irradiation_cell_contract_type, &
                                           resolve_irradiation_cell_contract, irradiation_cell_ok, irradiation_cell_invalid
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp
  implicit none
  type(app_config) :: app
  type(irradiation_cell_contract_type) :: cell
  integer(i32) :: status, face_bc(6)
  character(len=256) :: message

  call test_init(5)
  call test_begin('resolved_product_geometry')
  call default_app_config(app)
  app%sim%use_box = .true.
  app%sim%field_bc_mode = 'periodic2'
  app%sim%box_min = [-4.0_dp, -5.0_dp, -6.0_dp]
  app%sim%box_max = [7.0_dp, 8.0_dp, 9.0_dp]
  app%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  app%sim%bc_high = app%sim%bc_low
  app%periodic2%lower_boundary_model = 'e_bottom_zero'
  call resolve_irradiation_cell_contract(app, cell, status, message)
  call assert_equal_i32(status, irradiation_cell_ok, 'valid product configuration')
  call assert_true(len_trim(message) == 0, 'successful resolution has no diagnostic')
  call assert_equal_i32(cell%wall_face, 6_i32, 'bem_boundary face numbering')
  call assert_equal_i32(cell%wall_axis, 3_i32, 'wall axis')
  call assert_equal_i32(cell%outward_normal_sign, 1_i32, 'outward sign')
  call assert_equal_i32(cell%wall_face_bit(), 32_i32, 'wall bit')
  call assert_close_dp(cell%wall_coordinate(), 9.0_dp, 0.0_dp, 'wall coordinate from domain')
  call assert_true(cell%wall_face_name() == 'z_high', 'wall face name')
  call assert_true(all(cell%periodic_axes .eqv. [.true., .true., .false.]), 'resolved periodic axes')
  call assert_true(cell%lower_boundary_model == 'e_bottom_zero', 'resolved lower boundary')
  call assert_close_dp(cell%wall_area(app%sim%box_min, app%sim%box_max), 143.0_dp, 0.0_dp, 'wall area')
  call assert_true(cell%wall_tangents_periodic(), 'periodic wall tangents')
  app%periodic2%lower_boundary_model = 'symmetric_vacuum'
  call resolve_irradiation_cell_contract(app, cell, status, message)
  call assert_equal_i32(status, irradiation_cell_ok, 'second supported lower model')
  call assert_true(cell%lower_boundary_model == 'symmetric_vacuum', 'lower model is not hardcoded')
  call test_end()

  call test_begin('outward_open_crossing_and_normal_velocity')
  face_bc = bc_open
  call assert_close_dp(cell%outward_normal_velocity([4.0_dp, -7.0_dp, 2.0_dp]), 2.0_dp, 0.0_dp, 'normal component')
  call assert_true(cell%crosses_wall_outward(32_i32, face_bc, [4.0_dp, -7.0_dp, 2.0_dp]), 'open outward crossing')
  call assert_true(cell%crosses_wall_outward(33_i32, face_bc, [4.0_dp, -7.0_dp, 2.0_dp]), 'corner includes wall')
  call assert_true(.not. cell%crosses_wall_outward(1_i32, face_bc, [0.0_dp, 0.0_dp, 2.0_dp]), 'other face')
  call assert_true(.not. cell%crosses_wall_outward(32_i32, face_bc, [0.0_dp, 0.0_dp, -2.0_dp]), 'inward')
  call assert_true(.not. cell%crosses_wall_outward(32_i32, face_bc, [1.0_dp, 2.0_dp, 0.0_dp]), 'tangent')
  face_bc(6) = bc_reflect
  call assert_true(.not. cell%crosses_wall_outward(32_i32, face_bc, [0.0_dp, 0.0_dp, 2.0_dp]), 'reflecting wall')
  call test_end()

  call test_begin('queries_follow_contract_orientation')
  ! 問い合わせは現在の製品方向を再び直書きしない。負向き x 面を幾何 oracle とする。
  cell%wall_face = 1_i32
  cell%wall_axis = 1_i32
  cell%outward_normal_sign = -1_i32
  cell%wall_position = -4.0_dp
  cell%periodic_axes = [.false., .true., .true.]
  face_bc = bc_open
  call assert_equal_i32(cell%wall_face_bit(), 1_i32, 'low-face bit')
  call assert_close_dp(cell%outward_normal_velocity([-2.0_dp, 7.0_dp, 99.0_dp]), 2.0_dp, 0.0_dp, 'low-face normal')
  call assert_true(cell%crosses_wall_outward(1_i32, face_bc, [-2.0_dp, 7.0_dp, 99.0_dp]), 'low-face outward')
  call assert_true(cell%wall_face_name() == 'x_low', 'low-face name')
  call assert_close_dp(cell%wall_coordinate(), -4.0_dp, 0.0_dp, 'low-face coordinate')
  call assert_close_dp(cell%wall_area(app%sim%box_min, app%sim%box_max), 195.0_dp, 0.0_dp, 'low-face area')
  call assert_true(cell%wall_tangents_periodic(), 'low-face tangents')
  call test_end()

  call test_begin('unexpected_product_configuration_returns_diagnostic')
  app%sim%bc_low(1) = bc_open
  call resolve_irradiation_cell_contract(app, cell, status, message)
  call assert_equal_i32(status, irradiation_cell_invalid, 'incomplete periodic axis rejected')
  call assert_true(index(message, 'x/y periodic') > 0, 'direction diagnostic')
  app%sim%bc_low(1) = bc_periodic
  app%sim%bc_high(3) = bc_periodic
  call resolve_irradiation_cell_contract(app, cell, status, message)
  call assert_equal_i32(status, irradiation_cell_invalid, 'periodic wall axis rejected')
  app%sim%bc_high(3) = bc_open
  app%sim%use_box = .false.
  call resolve_irradiation_cell_contract(app, cell, status, message)
  call assert_equal_i32(status, irradiation_cell_invalid, 'missing product domain rejected')
  app%sim%use_box = .true.
  app%periodic2%lower_boundary_model = 'unknown'
  call resolve_irradiation_cell_contract(app, cell, status, message)
  call assert_equal_i32(status, irradiation_cell_invalid, 'unknown lower model rejected')
  call assert_true(index(message, 'fields.periodic.lower_boundary_model') > 0, 'public key diagnostic')
  call test_end()

  call test_begin('free_space_reference_remains_supported')
  call default_app_config(app)
  call resolve_irradiation_cell_contract(app, cell, status, message)
  call assert_equal_i32(status, irradiation_cell_ok, 'reference configuration')
  call assert_true(.not. any(cell%periodic_axes), 'reference periodic axes')
  call assert_true(.not. cell%wall_tangents_periodic(), 'reference tangents')
  call assert_true(cell%lower_boundary_model == 'not_applicable', 'reference lower model')
  call test_end()
  call test_summary()
end program test_irradiation_cell_contract
