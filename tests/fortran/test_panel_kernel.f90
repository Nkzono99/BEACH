!> Analytic P0 triangle kernel against independent quadrature and side limits.
program test_panel_kernel
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: eps0, k_coulomb
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_quadrature, only: panel_oracle_potential_field, panel_singular_potential_oracle
  use bem_panel_kernel, only: panel_potential_field, panel_potential_field_translated, panel_field_sum_translated, &
                              panel_potential_sum_translated, panel_side_principal_value, panel_side_normal_plus, &
                              panel_side_normal_minus
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_equal_i32, assert_close_dp, &
                          assert_allclose_1d
  implicit none

  type(panel_geometry_type) :: geom, moved, small_geom, neighbor_geom, random_geom, reference_geom
  real(dp) :: v0(3), v1(3), v2(3), target(3), phi, field(3), phi_ref, field_ref(3)
  real(dp) :: phi_plus, phi_minus, field_plus(3), field_minus(3), phi_self, self_ref
  real(dp) :: rotation(3, 3), shift(3), moved_target(3), moved_phi, moved_field(3)
  real(dp) :: scale_value, offset_value, neighbor_field_plus(3), neighbor_field_minus(3)
  real(dp), parameter :: charge = 2.5e-12_dp
  integer(i32) :: status, scale_idx, offset_idx
  real(dp) :: random_vertex(3, 3), translated_phi, translated_field(3)
  integer :: trial, accepted_trials
  integer(i64) :: lcg_state
  integer(i32), parameter :: batch_panels = 24_i32, batch_entries = 150_i32
  type(panel_geometry_type) :: batch_geom(batch_panels)
  real(dp) :: batch_charge(batch_panels), batch_shift1(batch_entries), batch_shift2(batch_entries)
  real(dp) :: sum_ex, sum_ey, sum_ez, sum_phi, ref_ex, ref_ey, ref_ez, ref_phi
  integer(i32) :: batch_index(batch_entries), panel_idx, entry
  real(dp), parameter :: scales(5) = [1.0e-9_dp, 1.0e-6_dp, 1.0_dp, 1.0e3_dp, 1.0e6_dp]
  real(dp), parameter :: offsets(3) = [0.0_dp, 1.0e3_dp, 1.0e6_dp]

  v0 = [0.0_dp, 0.0_dp, 0.0_dp]
  v1 = [1.0_dp, 0.0_dp, 0.0_dp]
  v2 = [0.0_dp, 1.0_dp, 0.0_dp]
  call init_panel_geometry(v0, v1, v2, geom, status)

  call test_init(8)

  call test_begin('off_surface_analytic_matches_oracle')
  call assert_equal_i32(status, panel_geometry_ok, 'valid geometry status mismatch')
  target = [0.23_dp, 0.31_dp, 0.47_dp]
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi, field)
  call panel_oracle_potential_field(geom, charge, target, 32_i32, phi_ref, field_ref)
  call assert_close_dp(phi, phi_ref, 2.0e-11_dp*max(1.0_dp, abs(phi_ref)), 'off-surface potential mismatch')
  call assert_allclose_1d(field, field_ref, 2.0e-10_dp*max(1.0_dp, maxval(abs(field_ref))), &
                          'off-surface field mismatch')
  call test_end()

  call test_begin('geometry_status_is_scale_and_translation_invariant')
  do scale_idx = 1, size(scales)
    scale_value = scales(scale_idx)
    do offset_idx = 1, size(offsets)
      offset_value = offsets(offset_idx)
      shift = [offset_value, -0.5_dp*offset_value, 0.25_dp*offset_value]
      call init_panel_geometry( &
        shift, shift + [scale_value, 0.0_dp, 0.0_dp], shift + [0.0_dp, scale_value, 0.0_dp], moved, status &
        )
      call assert_equal_i32(status, panel_geometry_ok, 'translated/scaled triangle was classified as degenerate')
    end do
  end do
  call test_end()

  call test_begin('coplanar_neighbor_does_not_receive_normal_jump')
  scale_value = 1.0e-7_dp
  call init_panel_geometry( &
    [0.0_dp, 0.0_dp, 0.0_dp], [scale_value, 0.0_dp, 0.0_dp], [0.0_dp, scale_value, 0.0_dp], &
    small_geom, status &
    )
  call assert_equal_i32(status, panel_geometry_ok, 'small panel geometry status mismatch')
  call init_panel_geometry( &
    [scale_value, 0.0_dp, 0.0_dp], [scale_value, scale_value, 0.0_dp], &
    [0.0_dp, scale_value, 0.0_dp], neighbor_geom, status &
    )
  call assert_equal_i32(status, panel_geometry_ok, 'small neighbor geometry status mismatch')
  call panel_potential_field( &
    small_geom, charge, neighbor_geom%centroid, panel_side_normal_plus, phi_plus, neighbor_field_plus &
    )
  call panel_potential_field( &
    small_geom, charge, neighbor_geom%centroid, panel_side_normal_minus, phi_minus, neighbor_field_minus &
    )
  call assert_allclose_1d( &
    neighbor_field_plus, neighbor_field_minus, &
    1.0e-12_dp*max(1.0_dp, maxval(abs(neighbor_field_plus))), &
    'coplanar point outside the panel received a normal jump' &
    )
  call test_end()

  call test_begin('off_surface_potential_is_two_sided_symmetric')
  target = geom%centroid + 1.0e-8_dp*geom%normal
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi_plus, field_plus)
  target = geom%centroid - 1.0e-8_dp*geom%normal
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi_minus, field_minus)
  call assert_close_dp(phi_plus, phi_minus, 1.0e-12_dp*max(1.0_dp, abs(phi_plus)), 'two-sided potential mismatch')
  call test_end()

  call test_begin('on_surface_principal_value_and_jump')
  call panel_potential_field(geom, charge, geom%centroid, panel_side_principal_value, phi_self, field)
  call panel_potential_field(geom, charge, geom%centroid, panel_side_normal_plus, phi_plus, field_plus)
  call panel_potential_field(geom, charge, geom%centroid, panel_side_normal_minus, phi_minus, field_minus)
  call panel_singular_potential_oracle(geom, charge, geom%centroid, 32_i32, self_ref)
  call assert_close_dp(phi_self, self_ref, 2.0e-11_dp*max(1.0_dp, abs(self_ref)), 'self potential mismatch')
  call assert_close_dp(phi_plus, phi_self, 1.0e-14_dp, 'plus-side potential mismatch')
  call assert_close_dp(phi_minus, phi_self, 1.0e-14_dp, 'minus-side potential mismatch')
  call assert_allclose_1d( &
    field_plus - field_minus, charge/(geom%area*eps0)*geom%normal, &
    1.0e-12_dp*max(1.0_dp, abs(charge/(geom%area*eps0))), 'normal field jump mismatch' &
    )
  call assert_allclose_1d( &
    field, 0.5_dp*(field_plus + field_minus), &
    1.0e-12_dp*max(1.0_dp, maxval(abs(field))), 'principal-value field is not the side-trace average' &
    )
  call test_end()

  call test_begin('rigid_transform_and_charge_scaling')
  rotation(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  rotation(:, 2) = [-1.0_dp, 0.0_dp, 0.0_dp]
  rotation(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
  shift = [2.0_dp, -3.0_dp, 0.7_dp]
  call init_panel_geometry(matmul(rotation, v0) + shift, matmul(rotation, v1) + shift, &
                           matmul(rotation, v2) + shift, moved, status)
  target = [0.23_dp, 0.31_dp, 0.47_dp]
  moved_target = matmul(rotation, target) + shift
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi, field)
  call panel_potential_field(moved, 2.0_dp*charge, moved_target, panel_side_principal_value, moved_phi, moved_field)
  call assert_close_dp(moved_phi, 2.0_dp*phi, 1.0e-12_dp*max(1.0_dp, abs(phi)), 'transformed potential mismatch')
  call assert_allclose_1d(moved_field, 2.0_dp*matmul(rotation, field), &
                          1.0e-11_dp*max(1.0_dp, maxval(abs(field))), 'transformed field mismatch')
  call test_end()

  call test_begin('translated_kernel_matches_shifted_geometry_bitwise')
  ! The FMM near field evaluates panels translated by image shifts. The translated kernel must
  ! reproduce the shifted-geometry evaluation bit for bit, including the on-plane fallback.
  lcg_state = 20261022_i64
  accepted_trials = 0
  do trial = 1, 400
    random_vertex(:, 1) = [lcg_uniform(), lcg_uniform(), lcg_uniform()]
    random_vertex(:, 2) = [lcg_uniform(), lcg_uniform(), lcg_uniform()]
    random_vertex(:, 3) = [lcg_uniform(), lcg_uniform(), lcg_uniform()]
    call init_panel_geometry(random_vertex(:, 1), random_vertex(:, 2), random_vertex(:, 3), random_geom, status)
    if (status /= panel_geometry_ok) cycle
    accepted_trials = accepted_trials + 1
    shift = [3.0_dp*lcg_uniform() - 1.5_dp, 3.0_dp*lcg_uniform() - 1.5_dp, 0.0_dp]
    if (modulo(trial, 3) == 0) shift = 0.0_dp
    reference_geom = random_geom
    reference_geom%vertex = reference_geom%vertex + spread(shift, 2, 3)
    reference_geom%centroid = reference_geom%centroid + shift
    if (modulo(trial, 5) == 0) then
      target = reference_geom%centroid
    else
      target = [2.0_dp*lcg_uniform() - 0.5_dp, 2.0_dp*lcg_uniform() - 0.5_dp, 2.0_dp*lcg_uniform() - 0.5_dp]
    end if
    call panel_potential_field(reference_geom, charge, target, panel_side_principal_value, phi, field)
    call panel_potential_field_translated(random_geom, shift, charge, target, translated_phi, translated_field)
    call assert_close_dp(translated_phi, phi, 0.0_dp, 'translated potential is not bitwise identical')
    call assert_allclose_1d(translated_field, field, 0.0_dp, 'translated field is not bitwise identical')
  end do
  call assert_equal_i32(merge(1_i32, 0_i32, accepted_trials >= 300), 1_i32, 'too few non-degenerate random panels')
  call test_end()

  call test_begin('batched_near_field_sum_matches_sequential_translated_kernel')
  ! The FMM near field sums translated panels in list order; the batched sum must reproduce the
  ! sequential accumulation bit for bit, including zero charges and on-plane targets.
  lcg_state = 7_i64
  panel_idx = 0_i32
  do while (panel_idx < batch_panels)
    random_vertex(:, 1) = [lcg_uniform(), lcg_uniform(), lcg_uniform()]
    random_vertex(:, 2) = [lcg_uniform(), lcg_uniform(), lcg_uniform()]
    random_vertex(:, 3) = [lcg_uniform(), lcg_uniform(), lcg_uniform()]
    call init_panel_geometry(random_vertex(:, 1), random_vertex(:, 2), random_vertex(:, 3), random_geom, status)
    if (status /= panel_geometry_ok) cycle
    panel_idx = panel_idx + 1_i32
    batch_geom(panel_idx) = random_geom
    batch_charge(panel_idx) = charge*(2.0_dp*lcg_uniform() - 1.0_dp)
  end do
  batch_charge(5) = 0.0_dp
  do entry = 1_i32, batch_entries
    batch_index(entry) = 1_i32 + int(lcg_uniform()*real(batch_panels, dp), i32)
    batch_shift1(entry) = real(int(3.0_dp*lcg_uniform(), i32) - 1_i32, dp)*2.5_dp
    batch_shift2(entry) = real(int(3.0_dp*lcg_uniform(), i32) - 1_i32, dp)*2.5_dp
  end do
  do trial = 1, 3
    select case (trial)
    case (1)
      target = [0.4_dp, 0.3_dp, 1.7_dp]
    case (2)
      target = batch_geom(batch_index(2))%centroid + [batch_shift1(2), batch_shift2(2), 0.0_dp]
    case default
      target = [-1.2_dp, 2.9_dp, -0.6_dp]
    end select
    ref_ex = 0.0_dp
    ref_ey = 0.0_dp
    ref_ez = 0.0_dp
    ref_phi = 0.0_dp
    do entry = 1_i32, batch_entries
      if (abs(batch_charge(batch_index(entry))) <= tiny(1.0_dp)) cycle
      shift = [batch_shift1(entry), batch_shift2(entry), 0.0_dp]
      call panel_potential_field_translated( &
        batch_geom(batch_index(entry)), shift, batch_charge(batch_index(entry)), target, phi, field &
        )
      ref_ex = ref_ex + field(1)/k_coulomb
      ref_ey = ref_ey + field(2)/k_coulomb
      ref_ez = ref_ez + field(3)/k_coulomb
      ref_phi = ref_phi + phi/k_coulomb
    end do
    sum_ex = 0.0_dp
    sum_ey = 0.0_dp
    sum_ez = 0.0_dp
    sum_phi = 0.0_dp
    call panel_field_sum_translated( &
      batch_geom, batch_index, [1_i32, 2_i32], batch_shift1, batch_shift2, batch_charge, batch_entries, target, &
      sum_ex, sum_ey, sum_ez &
      )
    call panel_potential_sum_translated( &
      batch_geom, batch_index, [1_i32, 2_i32], batch_shift1, batch_shift2, batch_charge, batch_entries, target, sum_phi &
      )
    call assert_close_dp(sum_ex, ref_ex, 0.0_dp, 'batched ex is not bitwise identical')
    call assert_close_dp(sum_ey, ref_ey, 0.0_dp, 'batched ey is not bitwise identical')
    call assert_close_dp(sum_ez, ref_ez, 0.0_dp, 'batched ez is not bitwise identical')
    call assert_close_dp(sum_phi, ref_phi, 0.0_dp, 'batched potential is not bitwise identical')
  end do
  call test_end()

  call test_summary()

contains

  !> Deterministic uniform deviate in [0, 1) from a 64-bit linear congruential generator.
  real(dp) function lcg_uniform() result(value)
    lcg_state = modulo(6364136223846793005_i64*lcg_state + 1442695040888963407_i64, 2_i64**62)
    value = real(lcg_state/2_i64**9, dp)/real(2_i64**53, dp)
  end function lcg_uniform
end program test_panel_kernel
