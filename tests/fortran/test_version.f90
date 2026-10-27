!> build metadata の定数と build-info の書式を検証する。
program test_version
  use bem_version, only: beach_build_id, beach_build_info, beach_source_commit, beach_version, beach_version_mode
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true
  implicit none

  character(len=*), parameter :: expected_build_info = &
                                 'build_info_schema_version=1'//new_line('a')// &
                                 'build_version='//beach_version//new_line('a')// &
                                 'build_version_mode='//beach_version_mode//new_line('a')// &
                                 'build_source_commit='//beach_source_commit//new_line('a')// &
                                 'build_id='//beach_build_id

  call test_init(1)

  call test_begin('build_metadata_and_build_info_layout')
  call assert_true(len_trim(beach_version) > 0, 'beach_version must be non-empty')
  call assert_true(len_trim(beach_version_mode) > 0, 'beach_version_mode must be non-empty')
  call assert_true(len_trim(beach_source_commit) > 0, 'beach_source_commit must be non-empty')
  call assert_true(len_trim(beach_build_id) > 0, 'beach_build_id must be non-empty')
  call assert_true( &
    len(beach_build_info) == len(expected_build_info) .and. beach_build_info == expected_build_info, &
    'build-info schema, field order, or payload mismatch' &
    )
  call test_end()

  call test_summary()
end program test_version
