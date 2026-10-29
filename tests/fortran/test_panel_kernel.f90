!> Analytic P0 triangle kernel against independent quadrature and side limits.
program test_panel_kernel
  use, intrinsic :: ieee_exceptions, only: ieee_invalid, ieee_divide_by_zero, ieee_set_flag, ieee_get_flag
  use bem_kinds, only: dp, sp, i32, i64
  use bem_constants, only: eps0, k_coulomb
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_quadrature, only: panel_oracle_potential_field, panel_singular_potential_oracle
  use bem_panel_kernel, only: panel_potential_field, panel_potential_field_translated, panel_field_sum_translated, &
                              panel_potential_sum_translated, panel_field_sum_translated_sp, &
                              panel_potential_sum_translated_sp, panel_side_principal_value, panel_side_normal_plus, &
                              panel_side_normal_minus
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_equal_i32, assert_close_dp, &
                          assert_allclose_1d, assert_true
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
  logical :: invalid_raised, divide_by_zero_raised
  integer(i32), parameter :: batch_panels = 24_i32, batch_entries = 150_i32
  type(panel_geometry_type) :: batch_geom(batch_panels)
  real(dp) :: batch_charge(batch_panels), batch_shift1(batch_entries), batch_shift2(batch_entries)
  real(dp) :: sum_ex, sum_ey, sum_ez, sum_phi, ref_ex, ref_ey, ref_ez, ref_phi
  integer(i32) :: batch_index(batch_entries), panel_idx, entry
  real(dp) :: tab_vx(3, batch_entries), tab_vy(3, batch_entries), tab_vz(3, batch_entries)
  real(dp) :: tab_nx(batch_entries), tab_ny(batch_entries), tab_nz(batch_entries), tab_edge(3, batch_entries)
  real(dp) :: tab_ox(3, batch_entries), tab_oy(3, batch_entries), tab_oz(3, batch_entries), tab_area(batch_entries)
  real(sp) :: sp_vx(3, batch_entries), sp_vy(3, batch_entries), sp_vz(3, batch_entries)
  real(sp) :: sp_nx(batch_entries), sp_ny(batch_entries), sp_nz(batch_entries), sp_edge(3, batch_entries)
  real(sp) :: sp_ox(3, batch_entries), sp_oy(3, batch_entries), sp_oz(3, batch_entries), sp_area(batch_entries)
  real(dp) :: leaf_center(3), sp_ex, sp_ey, sp_ez, sp_phi, field_scale, potential_scale
  real(dp), parameter :: scales(5) = [1.0e-9_dp, 1.0e-6_dp, 1.0_dp, 1.0e3_dp, 1.0e6_dp]
  real(dp), parameter :: offsets(3) = [0.0_dp, 1.0e3_dp, 1.0e6_dp]

  v0 = [0.0_dp, 0.0_dp, 0.0_dp]
  v1 = [1.0_dp, 0.0_dp, 0.0_dp]
  v2 = [0.0_dp, 1.0_dp, 0.0_dp]
  call init_panel_geometry(v0, v1, v2, geom, status)

  call test_init(9)

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
    if (entry == 1_i32) batch_index(entry) = 5_i32
    batch_shift1(entry) = real(int(3.0_dp*lcg_uniform(), i32) - 1_i32, dp)*2.5_dp
    batch_shift2(entry) = real(int(3.0_dp*lcg_uniform(), i32) - 1_i32, dp)*2.5_dp
    ! The FMM plan stores each near panel translated as `vertex + shift`, component by component.
    tab_vx(:, entry) = batch_geom(batch_index(entry))%vertex(1, :) + batch_shift1(entry)
    tab_vy(:, entry) = batch_geom(batch_index(entry))%vertex(2, :) + batch_shift2(entry)
    tab_vz(:, entry) = batch_geom(batch_index(entry))%vertex(3, :) + 0.0_dp
    tab_nx(entry) = batch_geom(batch_index(entry))%normal(1)
    tab_ny(entry) = batch_geom(batch_index(entry))%normal(2)
    tab_nz(entry) = batch_geom(batch_index(entry))%normal(3)
    tab_edge(:, entry) = batch_geom(batch_index(entry))%edge_length
    tab_ox(:, entry) = batch_geom(batch_index(entry))%edge_outward(1, :)
    tab_oy(:, entry) = batch_geom(batch_index(entry))%edge_outward(2, :)
    tab_oz(:, entry) = batch_geom(batch_index(entry))%edge_outward(3, :)
    tab_area(entry) = batch_geom(batch_index(entry))%area
  end do
  do trial = 1, 4
    select case (trial)
    case (1)
      target = [0.4_dp, 0.3_dp, 1.7_dp]
    case (2)
      target = batch_geom(batch_index(2))%centroid + [batch_shift1(2), batch_shift2(2), 0.0_dp]
    case (3)
      target = [-1.2_dp, 2.9_dp, -0.6_dp]
    case default
      ! A zero-charge panel contributes nothing even when its vertex coincides with the target.
      target = batch_geom(5)%vertex(:, 1) + [batch_shift1(1), batch_shift2(1), 0.0_dp]
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
    call ieee_set_flag(ieee_invalid, .false.)
    call ieee_set_flag(ieee_divide_by_zero, .false.)
    call panel_field_sum_translated( &
      batch_entries, tab_vx, tab_vy, tab_vz, tab_nx, tab_ny, tab_nz, tab_edge, tab_ox, tab_oy, tab_oz, tab_area, &
      batch_index, batch_charge, batch_geom, target, sum_ex, sum_ey, sum_ez &
      )
    call panel_potential_sum_translated( &
      batch_entries, tab_vx, tab_vy, tab_vz, tab_nx, tab_ny, tab_nz, tab_edge, tab_ox, tab_oy, tab_oz, tab_area, &
      batch_index, batch_charge, batch_geom, target, sum_phi &
      )
    call ieee_get_flag(ieee_invalid, invalid_raised)
    call ieee_get_flag(ieee_divide_by_zero, divide_by_zero_raised)
    call assert_true(.not. invalid_raised, 'batched near field raised an invalid operation')
    call assert_true(.not. divide_by_zero_raised, 'batched near field divided by zero')
    call assert_close_dp(sum_ex, ref_ex, 0.0_dp, 'batched ex is not bitwise identical')
    call assert_close_dp(sum_ey, ref_ey, 0.0_dp, 'batched ey is not bitwise identical')
    call assert_close_dp(sum_ez, ref_ez, 0.0_dp, 'batched ez is not bitwise identical')
    call assert_close_dp(sum_phi, ref_phi, 0.0_dp, 'batched potential is not bitwise identical')
  end do
  call test_end()

  call test_begin('single_precision_near_field_sum_stays_within_1e-5')
  ! The optional single-precision tables hold the panels relative to a leaf center; the sums must agree
  ! with the double-precision sequential kernel to the single-precision level and use the double
  ! kernel for on-plane targets.
  leaf_center = [0.3_dp, 0.4_dp, 0.5_dp]
  do entry = 1_i32, batch_entries
    sp_vx(:, entry) = real(tab_vx(:, entry) - leaf_center(1), sp)
    sp_vy(:, entry) = real(tab_vy(:, entry) - leaf_center(2), sp)
    sp_vz(:, entry) = real(tab_vz(:, entry) - leaf_center(3), sp)
    sp_nx(entry) = real(tab_nx(entry), sp)
    sp_ny(entry) = real(tab_ny(entry), sp)
    sp_nz(entry) = real(tab_nz(entry), sp)
    sp_edge(:, entry) = real(tab_edge(:, entry), sp)
    sp_ox(:, entry) = real(tab_ox(:, entry), sp)
    sp_oy(:, entry) = real(tab_oy(:, entry), sp)
    sp_oz(:, entry) = real(tab_oz(:, entry), sp)
    sp_area(entry) = real(tab_area(entry), sp)
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
    field_scale = 0.0_dp
    potential_scale = 0.0_dp
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
      ! Random-sign charges cancel in the sum; the precision bound applies to the contributions.
      field_scale = field_scale + sqrt(sum(field*field))/k_coulomb
      potential_scale = potential_scale + abs(phi)/k_coulomb
    end do
    sp_ex = 0.0_dp
    sp_ey = 0.0_dp
    sp_ez = 0.0_dp
    sp_phi = 0.0_dp
    call panel_field_sum_translated_sp( &
      batch_entries, sp_vx, sp_vy, sp_vz, sp_nx, sp_ny, sp_nz, sp_edge, sp_ox, sp_oy, sp_oz, sp_area, batch_index, &
      batch_charge, batch_geom, [1_i32, 2_i32], batch_shift1, batch_shift2, target - leaf_center, target, &
      sp_ex, sp_ey, sp_ez &
      )
    call panel_potential_sum_translated_sp( &
      batch_entries, sp_vx, sp_vy, sp_vz, sp_nx, sp_ny, sp_nz, sp_edge, sp_ox, sp_oy, sp_oz, sp_area, batch_index, &
      batch_charge, batch_geom, [1_i32, 2_i32], batch_shift1, batch_shift2, target - leaf_center, target, sp_phi &
      )
    write (*, '(a,i0,a,es10.2,a,es10.2)') 'single-precision sum trial ', trial, ': field rel err=', &
      sqrt((sp_ex - ref_ex)**2 + (sp_ey - ref_ey)**2 + (sp_ez - ref_ez)**2)/max(field_scale, 1.0e-300_dp), &
      ' potential rel err=', abs(sp_phi - ref_phi)/max(potential_scale, 1.0e-300_dp)
    call assert_close_dp(sqrt((sp_ex - ref_ex)**2 + (sp_ey - ref_ey)**2 + (sp_ez - ref_ez)**2), 0.0_dp, &
                         1.0e-5_dp*field_scale, 'single-precision field sum differs from the double kernel by more than 1e-5')
    call assert_close_dp(sp_phi, ref_phi, 1.0e-5_dp*potential_scale, &
                         'single-precision potential sum differs from the double kernel by more than 1e-5')
  end do
  call test_end()

  call test_summary()

contains

  !> Deterministic uniform deviate in [0, 1) from a 64-bit linear congruential generator.
  !> Park-Miller minimal standard generator (multiplier 48271, modulus 2**31 - 1). The product of a
  !! 16-bit multiplier and a state below 2**31 stays far inside int64, so the sample sequence never
  !! depends on signed overflow. Returns a value in (0, 1).
  real(dp) function lcg_uniform() result(value)
    integer(i64), parameter :: modulus = 2147483647_i64

    lcg_state = modulo(48271_i64*lcg_state, modulus)
    value = real(lcg_state, dp)/real(modulus, dp)
  end function lcg_uniform
end program test_panel_kernel
