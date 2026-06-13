program test_version
  use bem_version, only: beach_version, beach_version_mode
  implicit none

  call assert_true(len_trim(beach_version) > 0, 'beach_version must be non-empty')
  call assert_true(len_trim(beach_version_mode) > 0, 'beach_version_mode must be non-empty')

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine assert_true
end program test_version
