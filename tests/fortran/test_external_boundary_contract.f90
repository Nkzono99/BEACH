!> 局所 reservoir 流入の実行時契約を検証する。
program test_external_boundary_contract
  use bem_kinds, only: i32
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_boundary_ok, external_boundary_invalid, &
    external_inflow_none, external_inflow_scalar_barrier, resolve_external_boundary_contract
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_equal_i32
  implicit none

  type(external_boundary_contract_type) :: contract
  integer(i32) :: status
  character(len=256) :: message

  call test_init(3)

  call test_begin('default_contract_is_source_vdf')
  call resolve_external_boundary_contract('none', contract, status, message)
  call assert_equal_i32(status, external_boundary_ok, 'default contract resolution failed')
  call assert_equal_i32(contract%inflow_map, external_inflow_none, 'default inflow mismatch')
  call test_end()

  call test_begin('scalar_inflow_contract')
  call resolve_external_boundary_contract('infinity_barrier', contract, status, message)
  call assert_equal_i32(status, external_boundary_ok, 'scalar contract resolution failed')
  call assert_equal_i32(contract%inflow_map, external_inflow_scalar_barrier, 'scalar inflow mismatch')
  call test_end()

  call test_begin('unknown_local_model_fails_closed')
  call resolve_external_boundary_contract('unknown', contract, status, message)
  call assert_equal_i32(status, external_boundary_invalid, 'unknown inflow model must be rejected')
  call test_end()

  call test_summary()
end program test_external_boundary_contract
