!> CSV 応答テーブルの読み込み・格子検証・内容 fingerprint の生成を担う。
submodule(bem_matching_plane_response) bem_matching_plane_response_io
  use bem_string_utils, only: is_decimal_real_token
  implicit none

  integer, parameter :: response_column_count = 11
  integer, parameter :: response_line_length = 4096
  integer, parameter :: initial_row_capacity = 64
  integer(i64), parameter :: content_hash_modulus = 2147483647_i64
  integer(i64), parameter :: content_hash_multiplier_a = 65599_i64
  integer(i64), parameter :: content_hash_multiplier_b = 131071_i64
  character(len=*), parameter :: matching_plane_z_prefix = '# matching_plane_z_m='

  type :: matching_plane_axis_type
    real(dp), allocatable :: values(:)
  end type matching_plane_axis_type

  type :: matching_plane_content_hash_state_type
    integer(i64) :: a = 146959810_i64
    integer(i64) :: b = 109951162_i64
  end type matching_plane_content_hash_state_type

contains

  module procedure read_matching_plane_response_csv

  real(dp), allocatable :: rows(:, :)
  real(dp) :: parsed_row(response_column_count), matching_plane_z_m
  character(len=response_line_length) :: line
  character(len=:), allocatable :: left_adjusted
  integer :: unit_id, ios, line_number, row_count, capacity
  logical :: header_found, matching_plane_z_found

  table = matching_plane_response_table_type()
  call accept(status, message)
  capacity = initial_row_capacity
  allocate (rows(response_column_count, capacity))
  row_count = 0
  line_number = 0
  header_found = .false.
  matching_plane_z_found = .false.
  matching_plane_z_m = 0.0_dp

  open (newunit=unit_id, file=trim(path), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    call reject( &
      matching_plane_response_io_error, 'could not open matching-plane response CSV: '//trim(path), status, message &
      )
    return
  end if

  do
    read (unit_id, '(a)', iostat=ios) line
    if (ios < 0) exit
    line_number = line_number + 1
    if (ios > 0) then
      close (unit_id)
      call reject_line( &
        matching_plane_response_io_error, line_number, 'failed to read matching-plane response CSV.', status, message &
        )
      return
    end if
    if (len_trim(line) == 0) cycle
    left_adjusted = adjustl(trim(line))

    if (left_adjusted(1:1) == '#') then
      if (index(left_adjusted, '# matching_plane_z_m') == 1) then
        if (header_found) then
          close (unit_id)
          call reject_line( &
            matching_plane_response_invalid_metadata, line_number, &
            'matching_plane_z_m metadata must precede the CSV header.', status, message &
            )
          return
        end if
        if (matching_plane_z_found) then
          close (unit_id)
          call reject_line( &
            matching_plane_response_invalid_metadata, line_number, &
            'matching_plane_z_m metadata must appear exactly once.', status, message &
            )
          return
        end if
        call parse_matching_plane_z(left_adjusted, matching_plane_z_m, status, message)
        if (status /= matching_plane_response_ok) then
          close (unit_id)
          call prefix_line_number(line_number, message)
          return
        end if
        matching_plane_z_found = .true.
      end if
      cycle
    end if

    if (.not. header_found) then
      if (trim(line) /= matching_plane_response_csv_header) then
        close (unit_id)
        call reject_line( &
          matching_plane_response_invalid_header, line_number, &
          'matching-plane response CSV header does not match the exact v1 contract.', status, message &
          )
        return
      end if
      if (.not. matching_plane_z_found) then
        close (unit_id)
        call reject_line( &
          matching_plane_response_invalid_metadata, line_number, &
          'matching_plane_z_m metadata is required before the CSV header.', status, message &
          )
        return
      end if
      header_found = .true.
      cycle
    end if

    call parse_numeric_csv_row(line, parsed_row, status, message)
    if (status /= matching_plane_response_ok) then
      close (unit_id)
      call prefix_line_number(line_number, message)
      return
    end if
    if (row_count >= capacity) call grow_row_buffer(rows, capacity)
    row_count = row_count + 1
    rows(:, row_count) = parsed_row
  end do
  close (unit_id)

  if (.not. header_found) then
    call reject( &
      matching_plane_response_invalid_header, 'matching-plane response CSV contains no v1 header.', status, message &
      )
    return
  end if
  if (row_count <= 0) then
    call reject( &
      matching_plane_response_invalid_grid, 'matching-plane response CSV contains no numeric rows.', status, message &
      )
    return
  end if

  call build_matching_plane_response_table( &
    path, matching_plane_z_m, rows(:, :row_count), table, status, message &
    )
  end procedure read_matching_plane_response_csv

  subroutine build_matching_plane_response_table(path, matching_plane_z_m, rows, table, status, message)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: matching_plane_z_m
    real(dp), intent(in) :: rows(:, :)
    type(matching_plane_response_table_type), intent(out) :: table
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(matching_plane_axis_type) :: axis_data(matching_plane_response_input_count)
    logical, allocatable :: filled(:)
    integer(i32) :: coordinate_index(matching_plane_response_input_count)
    integer(i64) :: expected_points, linear_index, stride
    integer :: axis, row, max_axis_size, point_count

    table = matching_plane_response_table_type()
    call accept(status, message)
    if (size(rows, 1) /= response_column_count .or. size(rows, 2) <= 0) then
      call reject( &
        matching_plane_response_invalid_grid, 'matching-plane response rows have an invalid shape.', status, message &
        )
      return
    end if
    if (.not. ieee_is_finite(matching_plane_z_m)) then
      call reject( &
        matching_plane_response_invalid_metadata, 'matching_plane_z_m metadata must be finite.', status, message &
        )
      return
    end if
    if (any(.not. ieee_is_finite(rows))) then
      call reject( &
        matching_plane_response_invalid_row, 'matching-plane response values must all be finite.', status, message &
        )
      return
    end if
    if (any(rows(2:5, :) < 0.0_dp) .or. any(rows(7:8, :) < 0.0_dp)) then
      call reject( &
        matching_plane_response_invalid_row, &
        'matching-plane response fluxes and normal energies must be non-negative.', status, message &
        )
      return
    end if

    expected_points = 1_i64
    do axis = 1, matching_plane_response_input_count
      call unique_sorted_values(rows(axis, :), axis_data(axis)%values)
      table%axis_sizes(axis) = int(size(axis_data(axis)%values), i32)
      if (axis >= matching_plane_input_photoelectron_outward_flux .and. &
          table%axis_sizes(axis) > 1_i32 .and. &
          .not. any(axis_data(axis)%values == 0.0_dp)) then
        call reject( &
          matching_plane_response_invalid_grid, &
          'every matching-plane feedback axis must include zero.', status, message &
          )
        return
      end if
      if (expected_points > huge(expected_points)/int(table%axis_sizes(axis), i64)) then
        call reject( &
          matching_plane_response_invalid_grid, 'matching-plane Cartesian grid size overflow.', status, message &
          )
        return
      end if
      expected_points = expected_points*int(table%axis_sizes(axis), i64)
    end do
    if (expected_points > int(huge(0), i64)) then
      call reject( &
        matching_plane_response_invalid_grid, 'matching-plane Cartesian grid is too large.', status, message &
        )
      return
    end if
    if (expected_points /= int(size(rows, 2), i64)) then
      call reject( &
        matching_plane_response_invalid_grid, &
        'matching-plane response rows do not form one complete Cartesian product.', status, message &
        )
      return
    end if

    max_axis_size = maxval(table%axis_sizes)
    point_count = int(expected_points)
    allocate (table%axes(max_axis_size, matching_plane_response_input_count))
    table%axes = 0.0_dp
    do axis = 1, matching_plane_response_input_count
      table%axes(:table%axis_sizes(axis), axis) = axis_data(axis)%values
    end do
    allocate (table%response_values(matching_plane_response_output_count, point_count))
    allocate (filled(point_count), source=.false.)
    table%response_values = 0.0_dp

    do row = 1, size(rows, 2)
      do axis = 1, matching_plane_response_input_count
        coordinate_index(axis) = exact_axis_index(axis_data(axis)%values, rows(axis, row))
        if (coordinate_index(axis) <= 0_i32) then
          call reject( &
            matching_plane_response_invalid_grid, &
            'matching-plane row could not be mapped to its canonical axis.', status, message &
            )
          return
        end if
      end do
      linear_index = 1_i64
      stride = 1_i64
      do axis = 1, matching_plane_response_input_count
        linear_index = linear_index + int(coordinate_index(axis) - 1_i32, i64)*stride
        stride = stride*int(table%axis_sizes(axis), i64)
      end do
      if (filled(int(linear_index))) then
        call reject( &
          matching_plane_response_invalid_grid, &
          'matching-plane response grid contains a duplicate input coordinate.', status, message &
          )
        return
      end if
      filled(int(linear_index)) = .true.
      table%response_values(:, int(linear_index)) = rows(6:11, row)
    end do
    if (size(rows, 2) /= point_count .or. .not. all(filled)) then
      call reject( &
        matching_plane_response_invalid_grid, &
        'matching-plane response rows must form one complete Cartesian product.', status, message &
        )
      return
    end if

    table%source_path = trim(path)
    table%matching_plane_z_m = matching_plane_z_m
    table%content_fingerprint = compute_matching_plane_content_fingerprint(table)
    table%loaded = .true.
  end subroutine build_matching_plane_response_table

  subroutine parse_matching_plane_z(line, matching_plane_z_m, status, message)
    character(len=*), intent(in) :: line
    real(dp), intent(out) :: matching_plane_z_m
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=:), allocatable :: value_text
    integer :: ios

    matching_plane_z_m = 0.0_dp
    call accept(status, message)
    if (index(line, matching_plane_z_prefix) /= 1) then
      call reject( &
        matching_plane_response_invalid_metadata, &
        'matching_plane_z_m metadata must use "# matching_plane_z_m=<finite>".', status, message &
        )
      return
    end if
    value_text = trim(adjustl(line(len(matching_plane_z_prefix) + 1:)))
    if (len(value_text) == 0 .or. scan(value_text, ' ,'//achar(9)) > 0 .or. &
        .not. is_decimal_real_token(value_text)) then
      call reject( &
        matching_plane_response_invalid_metadata, 'matching_plane_z_m metadata has an invalid value.', status, message &
        )
      return
    end if
    read (value_text, *, iostat=ios) matching_plane_z_m
    if (ios /= 0 .or. .not. ieee_is_finite(matching_plane_z_m)) then
      matching_plane_z_m = 0.0_dp
      call reject( &
        matching_plane_response_invalid_metadata, 'matching_plane_z_m metadata must be finite.', status, message &
        )
    end if
  end subroutine parse_matching_plane_z

  subroutine parse_numeric_csv_row(line, values, status, message)
    character(len=*), intent(in) :: line
    real(dp), intent(out) :: values(response_column_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=:), allocatable :: record, token
    integer :: column, comma, first, last, ios

    values = 0.0_dp
    call accept(status, message)
    record = trim(line)
    first = 1
    do column = 1, response_column_count
      if (column < response_column_count) then
        comma = index(record(first:), ',')
        if (comma == 0) then
          call reject( &
            matching_plane_response_invalid_row, &
            'matching-plane response row must contain exactly 11 comma-separated values.', status, message &
            )
          return
        end if
        last = first + comma - 2
      else
        if (index(record(first:), ',') /= 0) then
          call reject( &
            matching_plane_response_invalid_row, &
            'matching-plane response row must contain exactly 11 comma-separated values.', status, message &
            )
          return
        end if
        last = len(record)
      end if
      if (last < first) then
        call reject( &
          matching_plane_response_invalid_row, 'matching-plane response CSV values must not be empty.', status, message &
          )
        return
      end if
      token = trim(adjustl(record(first:last)))
      if (len(token) == 0 .or. scan(token, ' '//achar(9)) > 0 .or. &
          .not. is_decimal_real_token(token)) then
        call reject( &
          matching_plane_response_invalid_row, 'matching-plane response CSV contains an invalid numeric token.', &
          status, message &
          )
        return
      end if
      read (token, *, iostat=ios) values(column)
      if (ios /= 0) then
        call reject( &
          matching_plane_response_invalid_row, 'matching-plane response CSV contains an invalid numeric token.', &
          status, message &
          )
        return
      end if
      first = last + 2
    end do
    if (any(.not. ieee_is_finite(values))) then
      call reject( &
        matching_plane_response_invalid_row, 'matching-plane response values must all be finite.', status, message &
        )
    end if
  end subroutine parse_numeric_csv_row

  subroutine grow_row_buffer(rows, capacity)
    real(dp), allocatable, intent(inout) :: rows(:, :)
    integer, intent(inout) :: capacity

    real(dp), allocatable :: grown(:, :)
    integer :: new_capacity

    if (capacity > huge(capacity)/2) error stop 'matching-plane response row capacity overflow.'
    new_capacity = 2*capacity
    allocate (grown(response_column_count, new_capacity))
    grown(:, :capacity) = rows
    call move_alloc(grown, rows)
    capacity = new_capacity
  end subroutine grow_row_buffer

  subroutine unique_sorted_values(values, unique)
    real(dp), intent(in) :: values(:)
    real(dp), allocatable, intent(out) :: unique(:)

    real(dp), allocatable :: sorted(:), compact(:)
    integer :: i, count

    allocate (sorted, source=values)
    call merge_sort_real(sorted)
    allocate (compact(size(sorted)))
    count = 1
    compact(1) = sorted(1)
    do i = 2, size(sorted)
      if (sorted(i) /= compact(count)) then
        count = count + 1
        compact(count) = sorted(i)
      end if
    end do
    allocate (unique(count), source=compact(:count))
  end subroutine unique_sorted_values

  subroutine merge_sort_real(values)
    real(dp), intent(inout) :: values(:)

    real(dp), allocatable :: work(:)
    integer :: width, left, middle, right, i, j, k, n

    n = size(values)
    if (n <= 1) return
    allocate (work(n))
    width = 1
    do while (width < n)
      left = 1
      do while (left <= n)
        middle = min(left + width, n + 1)
        right = min(left + 2*width - 1, n)
        i = left
        j = middle
        do k = left, right
          if (i >= middle) then
            work(k) = values(j)
            j = j + 1
          else if (j > right) then
            work(k) = values(i)
            i = i + 1
          else if (values(i) <= values(j)) then
            work(k) = values(i)
            i = i + 1
          else
            work(k) = values(j)
            j = j + 1
          end if
        end do
        left = left + 2*width
      end do
      values = work
      if (width > n/2) exit
      width = 2*width
    end do
  end subroutine merge_sort_real

  pure integer(i32) function exact_axis_index(axis_values, value) result(index_value)
    real(dp), intent(in) :: axis_values(:), value

    integer :: low, high, middle

    index_value = 0_i32
    low = 1
    high = size(axis_values)
    do while (low <= high)
      middle = (low + high)/2
      if (value == axis_values(middle)) then
        index_value = int(middle, i32)
        return
      else if (value < axis_values(middle)) then
        high = middle - 1
      else
        low = middle + 1
      end if
    end do
  end function exact_axis_index

  subroutine reject_line(code, line_number, text, status, message)
    integer(i32), intent(in) :: code
    integer, intent(in) :: line_number
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=512) :: contextual

    write (contextual, '(a,i0,a,a)') 'line ', line_number, ': ', trim(text)
    call reject(code, trim(contextual), status, message)
  end subroutine reject_line

  subroutine prefix_line_number(line_number, message)
    integer, intent(in) :: line_number
    character(len=*), intent(inout) :: message

    character(len=512) :: contextual

    write (contextual, '(a,i0,a,a)') 'line ', line_number, ': ', trim(message)
    message = trim(contextual)
  end subroutine prefix_line_number

  function compute_matching_plane_content_fingerprint(table) result(fingerprint)
    type(matching_plane_response_table_type), intent(in) :: table
    character(len=16) :: fingerprint
    type(matching_plane_content_hash_state_type) :: hash
    integer :: axis, point

    call content_hash_feed_string(hash, 'matching_plane_response_csv_v1')
    call content_hash_feed_real(hash, table%matching_plane_z_m)
    call content_hash_feed_integer(hash, matching_plane_response_input_count)
    do axis = 1, matching_plane_response_input_count
      call content_hash_feed_integer(hash, table%axis_sizes(axis))
      call content_hash_feed_real_vector(hash, table%axes(:table%axis_sizes(axis), axis))
    end do
    call content_hash_feed_integer(hash, matching_plane_response_output_count)
    call content_hash_feed_integer(hash, int(size(table%response_values, 2), i32))
    do point = 1, size(table%response_values, 2)
      call content_hash_feed_real_vector(hash, table%response_values(:, point))
    end do
    write (fingerprint, '(z8.8,z8.8)') hash%a, hash%b
  end function compute_matching_plane_content_fingerprint

  subroutine content_hash_feed_string(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    character(len=*), intent(in) :: value
    integer :: index

    call content_hash_feed_byte(hash, len_trim(value))
    do index = 1, len_trim(value)
      call content_hash_feed_byte(hash, iachar(value(index:index)))
    end do
  end subroutine content_hash_feed_string

  subroutine content_hash_feed_integer(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    integer(i32), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(i0)') value
    call content_hash_feed_string(hash, trim(encoded))
  end subroutine content_hash_feed_integer

  subroutine content_hash_feed_real(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    real(dp), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(es24.16e3)') value
    call content_hash_feed_string(hash, trim(adjustl(encoded)))
  end subroutine content_hash_feed_real

  subroutine content_hash_feed_real_vector(hash, values)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    real(dp), intent(in) :: values(:)
    integer :: index

    call content_hash_feed_integer(hash, int(size(values), i32))
    do index = 1, size(values)
      call content_hash_feed_real(hash, values(index))
    end do
  end subroutine content_hash_feed_real_vector

  subroutine content_hash_feed_byte(hash, value)
    type(matching_plane_content_hash_state_type), intent(inout) :: hash
    integer, intent(in) :: value

    hash%a = modulo(hash%a*content_hash_multiplier_a + int(value, i64) + 1_i64, content_hash_modulus)
    hash%b = modulo(hash%b*content_hash_multiplier_b + int(value, i64) + 1_i64, content_hash_modulus)
  end subroutine content_hash_feed_byte

  module procedure accept

  status = matching_plane_response_ok
  message = ''
  end procedure accept

  module procedure reject

  status = code
  message = text
  end procedure reject

end submodule bem_matching_plane_response_io
