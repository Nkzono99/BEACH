!> Matching-plane response CSV v1 の検証、補間、immutable cacheを検証する。
program test_matching_plane_response
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use bem_kinds, only: dp, i32
  use bem_matching_plane_response, only: &
    matching_plane_response_table_type, matching_plane_response_csv_header, &
    matching_plane_response_ok, matching_plane_response_invalid_argument, &
    matching_plane_response_invalid_header, matching_plane_response_invalid_metadata, &
    matching_plane_response_invalid_row, matching_plane_response_invalid_grid, &
    matching_plane_response_out_of_range, get_matching_plane_response_snapshot, &
    reset_matching_plane_response_snapshot_cache
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, assert_close_dp, &
    assert_allclose_1d, delete_file_if_exists
  implicit none

  character(len=*), parameter :: full_path = 'test_matching_plane_response_full.csv'
  character(len=*), parameter :: singleton_path = 'test_matching_plane_response_singleton.csv'
  character(len=*), parameter :: range_path = 'test_matching_plane_response_range.csv'
  character(len=*), parameter :: header_path = 'test_matching_plane_response_header.csv'
  character(len=*), parameter :: metadata_path = 'test_matching_plane_response_metadata.csv'
  character(len=*), parameter :: feedback_path = 'test_matching_plane_response_feedback.csv'
  character(len=*), parameter :: negative_path = 'test_matching_plane_response_negative.csv'
  character(len=*), parameter :: incomplete_path = 'test_matching_plane_response_incomplete.csv'
  character(len=*), parameter :: duplicate_path = 'test_matching_plane_response_duplicate.csv'
  character(len=*), parameter :: cache_path = 'test_matching_plane_response_cache.csv'
  character(len=*), parameter :: lexical_path = 'test_matching_plane_response_lexical.csv'

  type(matching_plane_response_table_type) :: table, cached_table
  real(dp) :: full_input(5, 32), full_output(6, 32)
  real(dp) :: singleton_input(5, 1), singleton_output(6, 1)
  real(dp) :: query(5), actual(6), expected(6), matching_plane_z_m
  real(dp), allocatable :: fingerprint_axes(:), fingerprint_values(:, :)
  integer(i32), allocatable :: fingerprint_axis_sizes(:)
  integer(i32) :: status
  character(len=512) :: message

  call cleanup_files()
  call reset_matching_plane_response_snapshot_cache()
  call test_init(8)

  call test_begin('full_cartesian_multilinear_and_fingerprint')
  call build_full_grid(full_input, full_output)
  call write_response_csv(full_path, full_input, full_output, 1.25_dp)
  call get_matching_plane_response_snapshot(full_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'full response table load failed: '//trim(message))
  call assert_true(table%is_loaded(), 'loaded response table did not report loaded state')
  query = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  call affine_response(query, expected)
  call table%evaluate(query, actual, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'full response interpolation failed: '//trim(message))
  call assert_allclose_1d(actual, expected, 2.0e-12_dp, 'five-dimensional affine interpolation mismatch')
  call table%get_matching_plane_z(matching_plane_z_m, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'matching-plane z getter failed')
  call assert_close_dp(matching_plane_z_m, 1.25_dp, 0.0_dp, 'matching-plane z metadata mismatch')
  call table%get_fingerprint_data( &
    fingerprint_axis_sizes, fingerprint_axes, fingerprint_values, matching_plane_z_m, status, message &
    )
  call assert_equal_i32(status, matching_plane_response_ok, 'fingerprint getter failed')
  call assert_true(all(fingerprint_axis_sizes == 2_i32), 'full-grid fingerprint axis sizes mismatch')
  call assert_allclose_1d( &
    fingerprint_axes, [-1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 8.0_dp], &
    0.0_dp, 'fingerprint axes are not canonical' &
    )
  call affine_response([-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], expected)
  call assert_allclose_1d( &
    fingerprint_values(:, 1), expected, 0.0_dp, 'axis-1-fastest first response value mismatch' &
    )
  call affine_response([1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp], expected)
  call assert_allclose_1d( &
    fingerprint_values(:, size(fingerprint_values, 2)), expected, 0.0_dp, &
    'axis-1-fastest last response value mismatch' &
    )
  call test_end()

  call test_begin('response_csv_numeric_grammar_is_strict')
  call write_lexical_token_table(lexical_path, '/', metadata_token=.true.)
  call get_matching_plane_response_snapshot(lexical_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_metadata, 'metadata slash token was accepted')
  call write_lexical_token_table(lexical_path, '/')
  call get_matching_plane_response_snapshot(lexical_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_row, 'data slash token was accepted')
  call write_lexical_token_table(lexical_path, '2*0')
  call get_matching_plane_response_snapshot(lexical_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_row, 'data repeat token was accepted')
  call write_lexical_token_table(lexical_path, '')
  call get_matching_plane_response_snapshot(lexical_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_row, 'data null token was accepted')
  call write_lexical_token_table(lexical_path, '0,0')
  call get_matching_plane_response_snapshot(lexical_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_row, 'extra response column was accepted')
  call test_end()

  call test_begin('singleton_axes_are_inactive')
  singleton_input(:, 1) = 0.0_dp
  singleton_output(:, 1) = [3.0_dp, 4.0_dp, 5.0_dp, -1.0_dp, 2.0_dp, 6.0_dp]
  call write_response_csv(singleton_path, singleton_input, singleton_output, -2.0_dp)
  call get_matching_plane_response_snapshot(singleton_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'singleton response table load failed')
  query = [-9.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp]
  call table%evaluate(query, actual, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'inactive singleton query was rejected')
  call assert_allclose_1d(actual, singleton_output(:, 1), 0.0_dp, 'singleton response changed with query')
  call test_end()

  call test_begin('active_axis_range_and_finite_query_fail_closed')
  call build_full_grid(full_input, full_output)
  call write_response_csv(range_path, full_input, full_output, 1.25_dp)
  call get_matching_plane_response_snapshot(range_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'range response table load failed')
  query = [1.0_dp + 16.0_dp*epsilon(1.0_dp), 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  call affine_response([1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], expected)
  call table%evaluate(query, actual, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'endpoint roundoff guard rejected a near-bound query')
  call assert_allclose_1d(actual, expected, 2.0e-12_dp, 'near-bound query was not snapped to the endpoint')
  query = [2.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  call table%evaluate(query, actual, status, message)
  call assert_equal_i32(status, matching_plane_response_out_of_range, 'active-axis extrapolation was not rejected')
  query = 0.0_dp
  query(2) = ieee_value(0.0_dp, ieee_quiet_nan)
  call table%evaluate(query, actual, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_argument, 'non-finite query was not rejected')
  call test_end()

  call test_begin('exact_header_and_required_metadata')
  singleton_input = 0.0_dp
  singleton_output = 0.0_dp
  call write_response_csv( &
    header_path, singleton_input, singleton_output, 0.0_dp, &
    header='invalid_'//matching_plane_response_csv_header &
    )
  call get_matching_plane_response_snapshot(header_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_header, 'non-exact v1 header was accepted')
  call write_response_csv(metadata_path, singleton_input, singleton_output, 0.0_dp, include_metadata=.false.)
  call get_matching_plane_response_snapshot(metadata_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_metadata, 'missing matching-plane z was accepted')
  call write_response_csv( &
    metadata_path, singleton_input, singleton_output, 0.0_dp, duplicate_metadata=.true. &
    )
  call get_matching_plane_response_snapshot(metadata_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_metadata, 'duplicate matching-plane z was accepted')
  call test_end()

  call test_begin('feedback_zero_and_nonnegative_contract')
  call write_missing_zero_table(feedback_path)
  call get_matching_plane_response_snapshot(feedback_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_grid, 'feedback axis without zero was accepted')
  singleton_input = 0.0_dp
  singleton_output = 1.0_dp
  singleton_input(4, 1) = -1.0_dp
  call write_response_csv(negative_path, singleton_input, singleton_output, 0.0_dp)
  call get_matching_plane_response_snapshot(negative_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_row, 'negative outward flux was accepted')
  singleton_input = 0.0_dp
  singleton_output(2, 1) = -1.0_dp
  call write_response_csv(negative_path, singleton_input, singleton_output, 0.0_dp)
  call get_matching_plane_response_snapshot(negative_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_row, 'negative inward response flux was accepted')
  call test_end()

  call test_begin('complete_cartesian_and_unique_coordinate_contract')
  call write_incomplete_table(incomplete_path)
  call get_matching_plane_response_snapshot(incomplete_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_grid, 'incomplete Cartesian product was accepted')
  call write_duplicate_table(duplicate_path)
  call get_matching_plane_response_snapshot(duplicate_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_invalid_grid, 'duplicate input coordinate was accepted')
  call test_end()

  call test_begin('path_snapshot_is_immutable_until_reset')
  singleton_input = 0.0_dp
  singleton_output(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
  call write_response_csv(cache_path, singleton_input, singleton_output, 0.5_dp)
  call get_matching_plane_response_snapshot(cache_path, table, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'initial cache table load failed')
  singleton_output(1, 1) = 9.0_dp
  call write_response_csv(cache_path, singleton_input, singleton_output, 0.75_dp)
  call get_matching_plane_response_snapshot(cache_path, cached_table, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'warm cache lookup failed')
  call cached_table%evaluate([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], actual, status, message)
  call assert_close_dp(actual(1), 1.0_dp, 0.0_dp, 'warm cache observed a rewritten file')
  call cached_table%get_matching_plane_z(matching_plane_z_m)
  call assert_close_dp(matching_plane_z_m, 0.5_dp, 0.0_dp, 'warm cache metadata changed after rewrite')
  call reset_matching_plane_response_snapshot_cache()
  call get_matching_plane_response_snapshot(cache_path, cached_table, status, message)
  call assert_equal_i32(status, matching_plane_response_ok, 'cache reload after reset failed')
  call cached_table%evaluate([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], actual, status, message)
  call assert_close_dp(actual(1), 9.0_dp, 0.0_dp, 'cache reset did not expose rewritten response')
  call cached_table%get_matching_plane_z(matching_plane_z_m)
  call assert_close_dp(matching_plane_z_m, 0.75_dp, 0.0_dp, 'cache reset did not expose rewritten metadata')
  call test_end()

  call cleanup_files()
  call reset_matching_plane_response_snapshot_cache()
  call test_summary()

contains

  subroutine build_full_grid(input, output)
    real(dp), intent(out) :: input(5, 32), output(6, 32)

    real(dp), parameter :: axes(2, 5) = reshape( &
                           [-1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 8.0_dp], [2, 5] &
                           )
    integer :: i1, i2, i3, i4, i5, row

    row = 0
    do i5 = 2, 1, -1
      do i4 = 2, 1, -1
        do i3 = 2, 1, -1
          do i2 = 2, 1, -1
            do i1 = 2, 1, -1
              row = row + 1
              input(:, row) = [axes(i1, 1), axes(i2, 2), axes(i3, 3), axes(i4, 4), axes(i5, 5)]
              call affine_response(input(:, row), output(:, row))
            end do
          end do
        end do
      end do
    end do
  end subroutine build_full_grid

  pure subroutine affine_response(input, output)
    real(dp), intent(in) :: input(5)
    real(dp), intent(out) :: output(6)

    integer :: component

    do component = 1, 6
      output(component) = 100.0_dp*real(component, dp) + &
                          real(component, dp)*input(1) + &
                          real(component + 1, dp)*input(2) + &
                          real(component + 2, dp)*input(3) + &
                          real(component + 3, dp)*input(4) + &
                          real(component + 4, dp)*input(5)
    end do
  end subroutine affine_response

  subroutine write_response_csv(path, input, output, matching_z, include_metadata, header, duplicate_metadata)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: input(:, :), output(:, :), matching_z
    logical, intent(in), optional :: include_metadata
    character(len=*), intent(in), optional :: header
    logical, intent(in), optional :: duplicate_metadata

    real(dp) :: row_values(11)
    integer :: unit_id, row
    logical :: write_metadata, write_duplicate_metadata

    write_metadata = .true.
    if (present(include_metadata)) write_metadata = include_metadata
    write_duplicate_metadata = .false.
    if (present(duplicate_metadata)) write_duplicate_metadata = duplicate_metadata
    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') '# provenance=test_matching_plane_response'
    if (write_metadata) write (unit_id, '(a,es24.16)') '# matching_plane_z_m=', matching_z
    if (write_metadata .and. write_duplicate_metadata) then
      write (unit_id, '(a,es24.16)') '# matching_plane_z_m=', matching_z
    end if
    write (unit_id, '(a)') ''
    if (present(header)) then
      write (unit_id, '(a)') header
    else
      write (unit_id, '(a)') matching_plane_response_csv_header
    end if
    write (unit_id, '(a)') '# numeric rows may be separated by comments'
    do row = 1, size(input, 2)
      row_values = [input(:, row), output(:, row)]
      write (unit_id, '(11(es24.16,:,","))') row_values
    end do
    close (unit_id)
  end subroutine write_response_csv

  subroutine write_missing_zero_table(path)
    character(len=*), intent(in) :: path

    real(dp) :: input(5, 2), output(6, 2)

    input = 0.0_dp
    input(2, :) = [1.0_dp, 2.0_dp]
    output = 1.0_dp
    call write_response_csv(path, input, output, 0.0_dp)
  end subroutine write_missing_zero_table

  subroutine write_incomplete_table(path)
    character(len=*), intent(in) :: path

    real(dp) :: input(5, 3), output(6, 3)

    input = 0.0_dp
    input(1, :) = [0.0_dp, 1.0_dp, 0.0_dp]
    input(2, :) = [0.0_dp, 0.0_dp, 1.0_dp]
    output = 1.0_dp
    call write_response_csv(path, input, output, 0.0_dp)
  end subroutine write_incomplete_table

  subroutine write_duplicate_table(path)
    character(len=*), intent(in) :: path

    real(dp) :: input(5, 4), output(6, 4)

    input = 0.0_dp
    input(1, :) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    input(4, :) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    output = 1.0_dp
    call write_response_csv(path, input, output, 0.0_dp)
  end subroutine write_duplicate_table

  subroutine write_lexical_token_table(path, token, metadata_token)
    character(len=*), intent(in) :: path, token
    logical, intent(in), optional :: metadata_token
    integer :: unit_id
    logical :: invalid_metadata

    invalid_metadata = .false.
    if (present(metadata_token)) invalid_metadata = metadata_token
    open (newunit=unit_id, file=path, status='replace', action='write')
    if (invalid_metadata) then
      write (unit_id, '(a)') '# matching_plane_z_m='//token
    else
      write (unit_id, '(a)') '# matching_plane_z_m=0.0'
    end if
    write (unit_id, '(a)') matching_plane_response_csv_header
    write (unit_id, '(a)') token//',0,0,0,0,0,0,0,0,0,0'
    close (unit_id)
  end subroutine write_lexical_token_table

  subroutine cleanup_files()
    call delete_file_if_exists(full_path)
    call delete_file_if_exists(singleton_path)
    call delete_file_if_exists(range_path)
    call delete_file_if_exists(header_path)
    call delete_file_if_exists(metadata_path)
    call delete_file_if_exists(feedback_path)
    call delete_file_if_exists(negative_path)
    call delete_file_if_exists(incomplete_path)
    call delete_file_if_exists(duplicate_path)
    call delete_file_if_exists(cache_path)
    call delete_file_if_exists(lexical_path)
  end subroutine cleanup_files

end program test_matching_plane_response
