!> Versioned stream codec for cached periodic root operators.
module bem_coulomb_fmm_periodic_cache
  use, intrinsic :: iso_fortran_env, only: int8, int64
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_kinds, only: dp, i32
  implicit none
  private

  character(len=16), parameter :: cache_magic = 'BEACH-P2-CACHE1'
  integer(i32), parameter, public :: periodic_cache_format_version = 1_i32
  integer(i32), parameter, public :: periodic_cache_ok = 0_i32
  integer(i32), parameter, public :: periodic_cache_miss = 1_i32
  integer(i32), parameter, public :: periodic_cache_invalid = 2_i32
  integer(i32), parameter, public :: periodic_cache_io_error = 3_i32
  integer(i32), parameter :: cache_max_ncoef = 4096_i32
  integer(i32), parameter :: cache_max_ntarget = 1000000_i32

  public :: write_periodic_operator_cache
  public :: read_periodic_operator_cache
  public :: periodic_operator_checksum
  public :: periodic_cache_fingerprint

contains

  subroutine write_periodic_operator_cache(path, fingerprint, target_nodes, operator, status)
    character(len=*), intent(in) :: path, fingerprint
    integer(i32), intent(in) :: target_nodes(:)
    real(dp), intent(in) :: operator(:, :, :)
    integer(i32), intent(out) :: status
    character(len=128) :: stored_fingerprint
    character(len=:), allocatable :: temporary_path
    integer(i32) :: ncoef, ntarget
    integer(int64) :: checksum
    integer :: unit, ios, rename_status

    status = periodic_cache_invalid
    ncoef = int(size(operator, 1), i32)
    ntarget = int(size(operator, 3), i32)
    if (ncoef <= 0_i32 .or. size(operator, 2) /= ncoef .or. ntarget /= size(target_nodes)) return
    if (len_trim(path) == 0 .or. len_trim(fingerprint) == 0 .or. len_trim(fingerprint) > len(stored_fingerprint)) return

    stored_fingerprint = ''
    stored_fingerprint = trim(fingerprint)
    checksum = periodic_operator_checksum(target_nodes, operator)
    temporary_path = trim(path)//'.tmp'
    open (newunit=unit, file=temporary_path, access='stream', form='unformatted', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      status = periodic_cache_io_error
      return
    end if
    write (unit, iostat=ios) cache_magic, periodic_cache_format_version, stored_fingerprint, ncoef, ntarget, checksum
    if (ios == 0) write (unit, iostat=ios) target_nodes
    if (ios == 0) write (unit, iostat=ios) operator
    close (unit)
    if (ios /= 0) then
      status = periodic_cache_io_error
      return
    end if
    call atomic_rename(temporary_path, trim(path), rename_status)
    if (rename_status /= filesystem_success) then
      status = periodic_cache_io_error
      return
    end if
    status = periodic_cache_ok
  end subroutine write_periodic_operator_cache

  subroutine read_periodic_operator_cache( &
    path, fingerprint, target_nodes, operator, status, expected_ncoef, expected_ntarget &
    )
    character(len=*), intent(in) :: path, fingerprint
    integer(i32), allocatable, intent(out) :: target_nodes(:)
    real(dp), allocatable, intent(out) :: operator(:, :, :)
    integer(i32), intent(out) :: status
    integer(i32), intent(in), optional :: expected_ncoef, expected_ntarget
    character(len=16) :: stored_magic
    character(len=128) :: stored_fingerprint
    integer(i32) :: stored_version, ncoef, ntarget
    integer(int64) :: stored_checksum, file_size, expected_size
    integer :: unit, ios
    logical :: exists

    inquire (file=trim(path), exist=exists)
    if (.not. exists) then
      status = periodic_cache_miss
      return
    end if
    open (newunit=unit, file=trim(path), access='stream', form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      status = periodic_cache_io_error
      return
    end if
    read (unit, iostat=ios) stored_magic, stored_version, stored_fingerprint, ncoef, ntarget, stored_checksum
    if (ios /= 0 .or. stored_magic /= cache_magic .or. stored_version /= periodic_cache_format_version .or. &
        trim(stored_fingerprint) /= trim(fingerprint) .or. ncoef <= 0_i32 .or. ntarget <= 0_i32 .or. &
        ncoef > cache_max_ncoef .or. ntarget > cache_max_ntarget) then
      close (unit)
      status = periodic_cache_invalid
      return
    end if
    if (present(expected_ncoef)) then
      if (ncoef /= expected_ncoef) then
        close (unit)
        status = periodic_cache_invalid
        return
      end if
    end if
    if (present(expected_ntarget)) then
      if (ntarget /= expected_ntarget) then
        close (unit)
        status = periodic_cache_invalid
        return
      end if
    end if
    inquire (file=trim(path), size=file_size, iostat=ios)
    expected_size = 16_int64 + 3_int64*int(storage_size(0_i32)/8, int64) + 128_int64 + &
                    int(storage_size(0_int64)/8, int64) + &
                    int(ntarget, int64)*int(storage_size(0_i32)/8, int64) + &
                    int(ncoef, int64)*int(ncoef, int64)*int(ntarget, int64)* &
                    int(storage_size(0.0_dp)/8, int64)
    if (ios /= 0 .or. file_size /= expected_size) then
      close (unit)
      status = periodic_cache_invalid
      return
    end if
    allocate (target_nodes(ntarget), operator(ncoef, ncoef, ntarget), stat=ios)
    if (ios /= 0) then
      close (unit)
      status = periodic_cache_io_error
      return
    end if
    read (unit, iostat=ios) target_nodes
    if (ios == 0) read (unit, iostat=ios) operator
    close (unit)
    if (ios /= 0 .or. periodic_operator_checksum(target_nodes, operator) /= stored_checksum) then
      deallocate (target_nodes, operator)
      status = periodic_cache_invalid
      return
    end if
    status = periodic_cache_ok
  end subroutine read_periodic_operator_cache

  pure integer(int64) function periodic_operator_checksum(target_nodes, operator) result(checksum)
    integer(i32), intent(in) :: target_nodes(:)
    real(dp), intent(in) :: operator(:, :, :)
    integer(int8), allocatable :: bytes(:)
    integer :: idx

    checksum = -3750763034362895579_int64
    bytes = transfer(target_nodes, bytes, size(target_nodes)*storage_size(0_i32)/8)
    do idx = 1, size(bytes)
      checksum = ieor(checksum, int(bytes(idx), int64))
      checksum = checksum*int(z'00000100000001B3', int64)
    end do
    bytes = transfer(operator, bytes, size(operator)*storage_size(0.0_dp)/8)
    do idx = 1, size(bytes)
      checksum = ieor(checksum, int(bytes(idx), int64))
      checksum = checksum*int(z'00000100000001B3', int64)
    end do
  end function periodic_operator_checksum

  pure function periodic_cache_fingerprint(tag, integer_values, real_values) result(fingerprint)
    character(len=*), intent(in) :: tag
    integer(i32), intent(in) :: integer_values(:)
    real(dp), intent(in) :: real_values(:)
    character(len=16) :: fingerprint
    integer(int8), allocatable :: bytes(:)
    integer(int64) :: checksum
    integer :: idx

    checksum = -3750763034362895579_int64
    bytes = transfer(tag, bytes, len(tag))
    call hash_bytes(bytes, checksum)
    bytes = transfer(integer_values, bytes, size(integer_values)*storage_size(0_i32)/8)
    call hash_bytes(bytes, checksum)
    bytes = transfer(real_values, bytes, size(real_values)*storage_size(0.0_dp)/8)
    call hash_bytes(bytes, checksum)
    write (fingerprint, '(z16.16)') checksum
  end function periodic_cache_fingerprint

  pure subroutine hash_bytes(bytes, checksum)
    integer(int8), intent(in) :: bytes(:)
    integer(int64), intent(inout) :: checksum
    integer :: idx

    do idx = 1, size(bytes)
      checksum = ieor(checksum, int(bytes(idx), int64))
      checksum = checksum*int(z'00000100000001B3', int64)
    end do
  end subroutine hash_bytes
end module bem_coulomb_fmm_periodic_cache
