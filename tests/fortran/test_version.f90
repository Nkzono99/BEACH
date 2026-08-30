program test_version
  use bem_version, only: beach_build_id, beach_build_info, beach_source_commit, beach_version, beach_version_mode
  implicit none

  character(len=*), parameter :: expected_build_info = &
                                 'build_info_schema_version=1'//new_line('a')// &
                                 'build_version='//beach_version//new_line('a')// &
                                 'build_version_mode='//beach_version_mode//new_line('a')// &
                                 'build_source_commit='//beach_source_commit//new_line('a')// &
                                 'build_id='//beach_build_id

  call assert_true(len_trim(beach_version) > 0, 'beach_version must be non-empty')
  call assert_true(len_trim(beach_version_mode) > 0, 'beach_version_mode must be non-empty')
  call assert_true(len_trim(beach_source_commit) > 0, 'beach_source_commit must be non-empty')
  call assert_true(len_trim(beach_build_id) > 0, 'beach_build_id must be non-empty')
  call assert_true( &
    len(beach_build_info) == len(expected_build_info) .and. beach_build_info == expected_build_info, &
    'build-info schema, field order, or payload mismatch' &
    )

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine assert_true
end program test_version
