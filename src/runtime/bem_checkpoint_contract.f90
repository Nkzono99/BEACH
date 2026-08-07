!> BEACH checkpoint metadata shared by output and restart code.
module bem_checkpoint_contract
  use bem_kinds, only: i32
  implicit none
  private

  integer(i32), parameter, public :: checkpoint_schema_version_current = 7_i32

end module bem_checkpoint_contract
