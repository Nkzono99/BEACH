!> Ordered mesh identity for mapping saved element charges.
module bem_mesh_identity
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: mesh_type
  implicit none
  private

  integer(i64), parameter :: hash_modulus = 2147483647_i64
  integer(i64), parameter :: hash_multiplier_a = 65599_i64
  integer(i64), parameter :: hash_multiplier_b = 131071_i64

  type :: hash_state
    integer(i64) :: a = 146959810_i64
    integer(i64) :: b = 109951162_i64
  end type hash_state

  public :: mesh_fingerprint

contains

  function mesh_fingerprint(mesh) result(fingerprint)
    type(mesh_type), intent(in) :: mesh
    character(len=16) :: fingerprint
    type(hash_state) :: hash
    integer(i32) :: elem

    call feed_integer(hash, mesh%nelem)
    do elem = 1, mesh%nelem
      call feed_real_vector(hash, mesh%v0(:, elem))
      call feed_real_vector(hash, mesh%v1(:, elem))
      call feed_real_vector(hash, mesh%v2(:, elem))
      if (allocated(mesh%elem_mesh_id)) call feed_integer(hash, mesh%elem_mesh_id(elem))
      if (allocated(mesh%elem_surface_model)) call feed_integer(hash, mesh%elem_surface_model(elem))
      if (allocated(mesh%elem_epsilon_r)) call feed_real(hash, mesh%elem_epsilon_r(elem))
      if (allocated(mesh%elem_vacuum_sign)) call feed_integer(hash, mesh%elem_vacuum_sign(elem))
    end do
    fingerprint = finish_hash(hash)
  end function mesh_fingerprint

  subroutine feed_string(hash, value)
    type(hash_state), intent(inout) :: hash
    character(len=*), intent(in) :: value
    integer :: index

    call feed_byte(hash, len_trim(value))
    do index = 1, len_trim(value)
      call feed_byte(hash, iachar(value(index:index)))
    end do
  end subroutine feed_string

  subroutine feed_integer(hash, value)
    type(hash_state), intent(inout) :: hash
    integer(i32), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(i0)') value
    call feed_string(hash, trim(encoded))
  end subroutine feed_integer

  subroutine feed_real(hash, value)
    type(hash_state), intent(inout) :: hash
    real(dp), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(es24.16e3)') value
    call feed_string(hash, trim(adjustl(encoded)))
  end subroutine feed_real

  subroutine feed_real_vector(hash, values)
    type(hash_state), intent(inout) :: hash
    real(dp), intent(in) :: values(:)
    integer :: index

    call feed_integer(hash, int(size(values), i32))
    do index = 1, size(values)
      call feed_real(hash, values(index))
    end do
  end subroutine feed_real_vector

  subroutine feed_byte(hash, value)
    type(hash_state), intent(inout) :: hash
    integer, intent(in) :: value

    hash%a = modulo(hash%a*hash_multiplier_a + int(value, i64) + 1_i64, hash_modulus)
    hash%b = modulo(hash%b*hash_multiplier_b + int(value, i64) + 1_i64, hash_modulus)
  end subroutine feed_byte

  function finish_hash(hash) result(fingerprint)
    type(hash_state), intent(in) :: hash
    character(len=16) :: fingerprint

    write (fingerprint, '(z8.8,z8.8)') hash%a, hash%b
  end function finish_hash

end module bem_mesh_identity
