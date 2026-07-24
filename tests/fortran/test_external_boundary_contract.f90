!> 旧設定語彙から解決する外部境界契約とface ownershipを検証する。
program test_external_boundary_contract
  use bem_kinds, only: i32
  use bem_types, only: bc_open, bc_periodic
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_boundary_ok, external_boundary_invalid, &
    external_inflow_none, external_inflow_scalar_barrier, &
    external_inflow_kinetic_profile, external_open_potential_barrier, external_transport_none, &
    external_transport_kinetic_1d, external_transport_unified_3d, &
    resolve_external_boundary_contract, external_boundary_owns_event
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  type(external_boundary_contract_type) :: contract
  integer(i32) :: status
  integer(i32) :: face_bc(6)
  character(len=256) :: message

  call test_init(8)

  call test_begin('default_contract_has_no_external_owner')
  call resolve( &
    'none', 'escape', 'none', 'absorbing_maxwellian', 'none', 'none', .false., &
    contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_ok, 'default contract resolution failed')
  call assert_equal_i32(contract%inflow_map, external_inflow_none, 'default inflow mismatch')
  call assert_equal_i32(contract%interface_transport, external_transport_none, 'default transport mismatch')
  call test_end()

  call test_begin('unknown_outer_model_fails_closed')
  call resolve( &
    'none', 'escape', 'unknown', 'absorbing_maxwellian', 'none', 'none', .false., &
    contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_invalid, 'unknown outer model must be rejected')
  call test_end()

  call test_begin('scalar_inflow_and_open_barrier_remain_separate')
  call resolve( &
    'infinity_barrier', 'potential_barrier', 'none', 'absorbing_maxwellian', 'none', 'none', .false., &
    contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_ok, 'scalar contract resolution failed')
  call assert_equal_i32(contract%inflow_map, external_inflow_scalar_barrier, 'scalar inflow mismatch')
  call assert_equal_i32(contract%ordinary_open_model, external_open_potential_barrier, 'open barrier mismatch')
  call assert_equal_i32(contract%interface_transport, external_transport_none, 'scalar contract must not own z-high')
  call test_end()

  call test_begin('kinetic_profile_contract_is_coherent')
  call resolve( &
    'none', 'potential_barrier', 'kinetic_1d', 'absorbing_maxwellian', &
    'kinetic_1d_profile_return', 'electrostatic_1d_instant_return', .false., contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_ok, 'kinetic contract resolution failed')
  call assert_equal_i32(contract%inflow_map, external_inflow_kinetic_profile, 'kinetic inflow mismatch')
  call assert_equal_i32(contract%interface_transport, external_transport_kinetic_1d, 'kinetic transport mismatch')
  call test_end()

  call test_begin('profile_transport_rejects_competing_inflow_owner')
  call resolve( &
    'infinity_barrier', 'escape', 'kinetic_1d', 'absorbing_maxwellian', &
    'kinetic_1d_profile_return', 'electrostatic_1d_instant_return', .false., contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_invalid, 'kinetic transport must reject scalar inflow')
  call test_end()

  call test_begin('kinetic_queue_contract_is_explicitly_delayed')
  call resolve( &
    'none', 'escape', 'kinetic_1d', 'zhao_charge_driven', 'kinetic_1d_profile_return', &
    'electrostatic_1d_instant_return', .true., contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_ok, 'queued kinetic contract resolution failed')
  call assert_equal_i32(contract%inflow_map, external_inflow_kinetic_profile, 'kinetic inflow mismatch')
  call assert_equal_i32(contract%interface_transport, external_transport_kinetic_1d, 'kinetic transport mismatch')
  call assert_true(contract%queue_enabled, 'queue contract flag mismatch')
  call test_end()

  call test_begin('unified_contract_uses_explicit_3d_transport')
  call resolve( &
    'none', 'escape', 'unified_linear_response', 'absorbing_maxwellian', &
    'electrostatic_3d_explicit_orbit', 'electrostatic_3d_explicit_orbit', .false., contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_ok, 'unified contract resolution failed')
  call assert_equal_i32(contract%interface_transport, external_transport_unified_3d, 'unified transport mismatch')
  call test_end()

  call test_begin('z_high_membership_owns_periodic_corner_event')
  call resolve( &
    'none', 'escape', 'kinetic_1d', 'absorbing_maxwellian', &
    'kinetic_1d_profile_return', 'electrostatic_1d_instant_return', .false., contract, status, message &
    )
  face_bc = [bc_periodic, bc_periodic, bc_periodic, bc_periodic, bc_open, bc_open]
  call assert_true( &
    external_boundary_owns_event(contract, ibset(ibset(0_i32, 1_i32), 5_i32), face_bc), &
    'z-high plus lateral periodic event must retain outer ownership' &
    )
  call assert_true( &
    .not. external_boundary_owns_event(contract, ibset(0_i32, 1_i32), face_bc), &
    'lateral event alone must remain an ordinary box event' &
    )
  call test_end()

  call test_summary()

contains

  subroutine resolve( &
    reservoir, open_model, outer, closure, return_model, transfer, queue, resolved, code, text &
    )
    character(len=*), intent(in) :: reservoir, open_model, outer, closure, return_model, transfer
    logical, intent(in) :: queue
    type(external_boundary_contract_type), intent(out) :: resolved
    integer(i32), intent(out) :: code
    character(len=*), intent(out) :: text

    call resolve_external_boundary_contract( &
      reservoir, open_model, outer, closure, return_model, transfer, queue, resolved, code, text &
      )
  end subroutine resolve

end program test_external_boundary_contract
