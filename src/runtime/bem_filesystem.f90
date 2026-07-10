!> Minimal POSIX filesystem operations used by the runtime.
module bem_filesystem
  use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_int, c_null_char, c_ptr
  implicit none
  private

  integer, parameter, public :: filesystem_success = 0
  integer, parameter, public :: filesystem_empty_path = 1
  integer, parameter, public :: filesystem_not_directory = 2
  integer, parameter, public :: filesystem_os_error = 3

  integer(c_int), parameter :: directory_mode = 511_c_int

  public :: create_directories

  interface
    integer(c_int) function c_mkdir(path, mode) bind(C, name='mkdir')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: mode
    end function c_mkdir

    type(c_ptr) function c_opendir(path) bind(C, name='opendir')
      import :: c_char, c_ptr
      character(kind=c_char), intent(in) :: path(*)
    end function c_opendir

    integer(c_int) function c_closedir(directory) bind(C, name='closedir')
      import :: c_int, c_ptr
      type(c_ptr), value :: directory
    end function c_closedir
  end interface

contains

  !> Recursively create a directory path without invoking a shell.
  subroutine create_directories(path, status)
    character(len=*), intent(in) :: path
    integer, intent(out) :: status

    integer :: i, path_length

    status = filesystem_success
    path_length = len_trim(path)
    if (path_length == 0) then
      status = filesystem_empty_path
      return
    end if

    do i = 1, path_length
      if (path(i:i) /= '/') cycle
      if (i <= 1) cycle
      call create_one_directory(path(:i - 1), status)
      if (status /= filesystem_success) return
    end do

    call create_one_directory(path(:path_length), status)
  end subroutine create_directories

  !> Create one path prefix, accepting it only when it can be opened as a directory.
  subroutine create_one_directory(path, status)
    character(len=*), intent(in) :: path
    integer, intent(out) :: status

    character(kind=c_char), allocatable :: c_path(:)
    integer(c_int) :: mkdir_status
    logical :: path_exists

    status = filesystem_success
    if (len(path) == 0) return
    if (is_directory(path)) return

    call to_c_string(path, c_path)
    mkdir_status = c_mkdir(c_path, directory_mode)
    if (is_directory(path)) return

    if (mkdir_status == 0_c_int) then
      status = filesystem_os_error
      return
    end if

    inquire (file=path, exist=path_exists)
    if (path_exists) then
      status = filesystem_not_directory
    else
      status = filesystem_os_error
    end if
  end subroutine create_one_directory

  !> Return whether a path can be opened as a POSIX directory.
  logical function is_directory(path)
    character(len=*), intent(in) :: path

    character(kind=c_char), allocatable :: c_path(:)
    type(c_ptr) :: directory
    integer(c_int) :: close_status

    call to_c_string(path, c_path)
    directory = c_opendir(c_path)
    is_directory = c_associated(directory)
    if (is_directory) close_status = c_closedir(directory)
  end function is_directory

  !> Convert a Fortran character value to a null-terminated C character array.
  subroutine to_c_string(path, c_path)
    character(len=*), intent(in) :: path
    character(kind=c_char), allocatable, intent(out) :: c_path(:)

    integer :: i

    allocate (c_path(len(path) + 1))
    do i = 1, len(path)
      c_path(i) = path(i:i)
    end do
    c_path(len(path) + 1) = c_null_char
  end subroutine to_c_string

end module bem_filesystem
