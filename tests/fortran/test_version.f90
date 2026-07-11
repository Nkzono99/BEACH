program test_version
  use bem_version, only: beach_build_id, beach_build_info, beach_source_commit, beach_version, beach_version_mode
  implicit none

  call assert_true(len_trim(beach_version) > 0, 'beach_version must be non-empty')
  call assert_true(len_trim(beach_version_mode) > 0, 'beach_version_mode must be non-empty')
  call assert_true(len_trim(beach_source_commit) > 0, 'beach_source_commit must be non-empty')
  call assert_true(len_trim(beach_build_id) > 0, 'beach_build_id must be non-empty')
  call assert_true(index(beach_build_info, 'build_info_schema_version=1') > 0, 'build-info schema missing')
  call assert_true(index(beach_build_info, 'build_version='//beach_version) > 0, 'build-info version mismatch')
  call assert_true(index(beach_build_info, 'build_version_mode='//beach_version_mode) > 0, &
                   'build-info version mode mismatch')
  call assert_true(index(beach_build_info, 'build_source_commit='//beach_source_commit) > 0, &
                   'build-info source commit mismatch')
  call assert_true(index(beach_build_info, 'build_id='//beach_build_id) > 0, 'build-info id mismatch')

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine assert_true
end program test_version
