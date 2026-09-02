!> Online Zhao evaluator から既存形式の matching-plane 応答表を生成する。
module bem_matching_plane_response_generator
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_app_config_types, only: app_config
  use bem_matching_plane_response, only: matching_plane_response_input_count, &
                                         matching_plane_response_output_count, &
                                         matching_plane_response_query_csv_header, &
                                         matching_plane_response_csv_header
  use bem_matching_plane_response_provider, only: matching_plane_response_provider_type, &
                                                  matching_plane_provider_ok
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_mpi, only: mpi_context
  use bem_string_utils, only: lower_ascii, is_decimal_real_token
  implicit none
  private

  integer, parameter :: query_line_length = 4096
  integer, parameter :: initial_query_capacity = 64

  integer(i32), parameter, public :: matching_plane_generator_ok = 0_i32
  integer(i32), parameter, public :: matching_plane_generator_invalid_argument = 1_i32
  integer(i32), parameter, public :: matching_plane_generator_io_error = 2_i32
  integer(i32), parameter, public :: matching_plane_generator_invalid_grid = 3_i32
  integer(i32), parameter, public :: matching_plane_generator_evaluation_failure = 4_i32

  public :: generate_matching_plane_zhao_response_table

contains

  !> 5列query gridを評価し、既存table backend互換の11列CSVを書く。
  subroutine generate_matching_plane_zhao_response_table(cfg, query_path, output_path, status, message)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: query_path, output_path
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(matching_plane_response_provider_type) :: provider
    type(mpi_context) :: serial_mpi
    real(dp), allocatable :: queries(:, :), responses(:, :)
    integer(i32) :: provider_status
    integer :: row, output_unit, ios, close_ios, rename_status
    character(len=512) :: provider_message
    character(len=:), allocatable :: temporary_output_path

    call accept_generator(status, message)
    if (trim(lower_ascii(cfg%surface_current%model)) /= 'matching_plane_quasistatic' .or. &
        trim(lower_ascii(cfg%surface_current%response_backend)) /= 'zhao_online') then
      call reject_generator( &
        matching_plane_generator_invalid_argument, &
        'response generation requires matching_plane_quasistatic with response_backend="zhao_online".', &
        status, message &
        )
      return
    end if
    if (trim(lower_ascii(cfg%surface_current%zhao_root_selection)) == 'continuation') then
      call reject_generator( &
        matching_plane_generator_invalid_argument, &
        'response-table generation does not define an accepted endpoint for Zhao continuation.', &
        status, message &
        )
      return
    end if
    if (len_trim(query_path) == 0 .or. len_trim(output_path) == 0) then
      call reject_generator( &
        matching_plane_generator_invalid_argument, &
        'query and output paths must not be empty.', status, message &
        )
      return
    end if
    if (trim(query_path) == trim(output_path)) then
      call reject_generator( &
        matching_plane_generator_invalid_argument, &
        'query and output paths must be different.', status, message &
        )
      return
    end if
    temporary_output_path = trim(output_path)//'.beach-zhao-response.tmp'
    if (trim(query_path) == temporary_output_path) then
      call reject_generator( &
        matching_plane_generator_invalid_argument, &
        'query path conflicts with the response generator temporary path.', status, message &
        )
      return
    end if

    serial_mpi = mpi_context()
    call provider%initialize(cfg, serial_mpi, provider_status, provider_message)
    if (provider_status /= matching_plane_provider_ok) then
      call reject_generator( &
        matching_plane_generator_evaluation_failure, &
        'Zhao response initialization failed: '//trim(provider_message), status, message &
        )
      return
    end if

    call read_query_grid(query_path, queries, status, message)
    if (status /= matching_plane_generator_ok) return
    call validate_cartesian_query_grid(queries, status, message)
    if (status /= matching_plane_generator_ok) return

    allocate (responses(matching_plane_response_output_count, size(queries, 2)))
    do row = 1, size(queries, 2)
      call provider%evaluate( &
        queries(:, row), serial_mpi, responses(:, row), provider_status, provider_message &
        )
      if (provider_status /= matching_plane_provider_ok) then
        call reject_generator( &
          matching_plane_generator_evaluation_failure, &
          'Zhao response evaluation failed at query row '//trim(integer_text(row))//': '// &
          trim(provider_message), status, message &
          )
        return
      end if
    end do

    open (newunit=output_unit, file=temporary_output_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      call reject_generator( &
        matching_plane_generator_io_error, &
        'could not open matching-plane response output: '//trim(output_path), status, message &
        )
      return
    end if
    write (output_unit, '(a,es24.16)', iostat=ios) '# matching_plane_z_m=', cfg%sim%box_max(3)
    if (ios == 0) write (output_unit, '(a)', iostat=ios) matching_plane_response_csv_header
    do row = 1, size(queries, 2)
      if (ios /= 0) exit
      write (output_unit, '(11(es24.16,:,","))', iostat=ios) queries(:, row), responses(:, row)
    end do
    close (output_unit, iostat=close_ios)
    if (ios == 0) ios = close_ios
    if (ios /= 0) then
      call reject_generator( &
        matching_plane_generator_io_error, &
        'failed while writing matching-plane response output.', status, message &
        )
      return
    end if
    call atomic_rename(temporary_output_path, trim(output_path), rename_status)
    if (rename_status /= filesystem_success) then
      call reject_generator( &
        matching_plane_generator_io_error, &
        'failed to atomically publish matching-plane response output.', status, message &
        )
    end if
  end subroutine generate_matching_plane_zhao_response_table

  subroutine read_query_grid(path, queries, status, message)
    character(len=*), intent(in) :: path
    real(dp), allocatable, intent(out) :: queries(:, :)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp), allocatable :: buffer(:, :), grown(:, :)
    real(dp) :: query(matching_plane_response_input_count)
    character(len=query_line_length) :: line
    character(len=:), allocatable :: record
    integer :: unit_id, ios, line_number, row_count, capacity
    logical :: header_found

    call accept_generator(status, message)
    capacity = initial_query_capacity
    allocate (buffer(matching_plane_response_input_count, capacity))
    row_count = 0
    line_number = 0
    header_found = .false.
    open (newunit=unit_id, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      allocate (queries(matching_plane_response_input_count, 0))
      call reject_generator( &
        matching_plane_generator_io_error, &
        'could not open matching-plane query CSV: '//trim(path), status, message &
        )
      return
    end if

    do
      read (unit_id, '(a)', iostat=ios) line
      if (ios < 0) exit
      line_number = line_number + 1
      if (ios > 0) then
        close (unit_id)
        allocate (queries(matching_plane_response_input_count, 0))
        call reject_generator( &
          matching_plane_generator_io_error, &
          'failed to read matching-plane query CSV line '//trim(integer_text(line_number))//'.', &
          status, message &
          )
        return
      end if
      record = trim(adjustl(line))
      if (len(record) == 0) cycle
      if (record(1:1) == '#') cycle
      if (.not. header_found) then
        if (record /= matching_plane_response_query_csv_header) then
          close (unit_id)
          allocate (queries(matching_plane_response_input_count, 0))
          call reject_generator( &
            matching_plane_generator_invalid_grid, &
            'matching-plane query CSV header does not match the required five-column contract.', &
            status, message &
            )
          return
        end if
        header_found = .true.
        cycle
      end if
      if (count_character(record, ',') /= matching_plane_response_input_count - 1_i32) then
        close (unit_id)
        allocate (queries(matching_plane_response_input_count, 0))
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'matching-plane query row '//trim(integer_text(line_number))//' must contain exactly five values.', &
          status, message &
          )
        return
      end if
      call parse_query_record(record, query, ios)
      if (ios /= 0 .or. any(.not. ieee_is_finite(query)) .or. any(query(2:5) < 0.0_dp)) then
        close (unit_id)
        allocate (queries(matching_plane_response_input_count, 0))
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'matching-plane query row '//trim(integer_text(line_number))//' contains invalid values.', &
          status, message &
          )
        return
      end if
      if (row_count == capacity) then
        if (capacity > huge(capacity)/2) then
          close (unit_id)
          allocate (queries(matching_plane_response_input_count, 0))
          call reject_generator( &
            matching_plane_generator_invalid_grid, &
            'matching-plane query row capacity overflowed.', status, message &
            )
          return
        end if
        allocate (grown(matching_plane_response_input_count, 2*capacity))
        grown(:, :capacity) = buffer
        call move_alloc(grown, buffer)
        capacity = 2*capacity
      end if
      row_count = row_count + 1
      buffer(:, row_count) = query
    end do
    close (unit_id)

    if (.not. header_found .or. row_count == 0) then
      allocate (queries(matching_plane_response_input_count, 0))
      call reject_generator( &
        matching_plane_generator_invalid_grid, &
        'matching-plane query CSV must contain the required header and at least one row.', &
        status, message &
        )
      return
    end if
    allocate (queries(matching_plane_response_input_count, row_count), source=buffer(:, :row_count))
  end subroutine read_query_grid

  !> List-directed repeat/null/slash syntaxを拒否し、5個の十進実数だけを読む。
  subroutine parse_query_record(record, query, ios)
    character(len=*), intent(in) :: record
    real(dp), intent(out) :: query(matching_plane_response_input_count)
    integer, intent(out) :: ios

    character(len=:), allocatable :: token
    integer :: column, comma, first, last

    query = 0.0_dp
    ios = 0
    first = 1
    do column = 1, matching_plane_response_input_count
      if (column < matching_plane_response_input_count) then
        comma = index(record(first:), ',')
        if (comma <= 0) then
          ios = 1
          return
        end if
        last = first + comma - 2
      else
        if (index(record(first:), ',') /= 0) then
          ios = 1
          return
        end if
        last = len(record)
      end if
      if (last < first) then
        ios = 1
        return
      end if
      token = trim(adjustl(record(first:last)))
      if (len(token) == 0 .or. scan(token, ' '//achar(9)) > 0 .or. &
          .not. is_decimal_real_token(token)) then
        ios = 1
        return
      end if
      read (token, *, iostat=ios) query(column)
      if (ios /= 0) return
      first = last + 2
    end do
  end subroutine parse_query_record

  subroutine validate_cartesian_query_grid(queries, status, message)
    real(dp), intent(in) :: queries(:, :)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp), allocatable :: axes(:, :)
    logical, allocatable :: filled(:)
    integer(i32) :: axis_sizes(matching_plane_response_input_count)
    integer :: axis, row, candidate, coordinate_index
    integer(i64) :: point_count, linear_index, stride
    logical :: found

    call accept_generator(status, message)
    if (size(queries, 1) /= matching_plane_response_input_count .or. size(queries, 2) <= 0) then
      call reject_generator( &
        matching_plane_generator_invalid_grid, &
        'matching-plane query grid has invalid dimensions.', status, message &
        )
      return
    end if
    allocate (axes(size(queries, 2), matching_plane_response_input_count))
    axes = 0.0_dp
    axis_sizes = 0_i32
    do axis = 1, matching_plane_response_input_count
      do row = 1, size(queries, 2)
        found = .false.
        do candidate = 1, axis_sizes(axis)
          if (axes(candidate, axis) == queries(axis, row)) then
            found = .true.
            exit
          end if
        end do
        if (found) cycle
        axis_sizes(axis) = axis_sizes(axis) + 1_i32
        axes(axis_sizes(axis), axis) = queries(axis, row)
      end do
      if (axis == 2 .and. .not. any(axes(:axis_sizes(axis), axis) == 0.0_dp)) then
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'the matching-plane photoelectron-flux axis must include zero.', status, message &
          )
        return
      end if
      if (axis == 3 .and. axis_sizes(axis) /= 1_i32) then
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'online Zhao table generation requires a singleton photoelectron-energy axis.', &
          status, message &
          )
        return
      end if
      if (axis == 3 .and. any(axes(:axis_sizes(2), 2) > 0.0_dp) .and. &
          axes(1, axis) <= 0.0_dp) then
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'positive photoelectron flux requires a positive singleton energy node.', &
          status, message &
          )
        return
      end if
      if (axis >= 4 .and. &
          (axis_sizes(axis) /= 1_i32 .or. axes(1, axis) /= 0.0_dp)) then
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'online Zhao ambient outward axes must be singleton zero nodes.', status, message &
          )
        return
      end if
    end do

    point_count = 1_i64
    do axis = 1, matching_plane_response_input_count
      if (point_count > huge(point_count)/int(axis_sizes(axis), i64)) then
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'matching-plane Cartesian query size overflowed.', status, message &
          )
        return
      end if
      point_count = point_count*int(axis_sizes(axis), i64)
    end do
    if (point_count /= int(size(queries, 2), i64)) then
      call reject_generator( &
        matching_plane_generator_invalid_grid, &
        'matching-plane queries must form one complete Cartesian product.', status, message &
        )
      return
    end if

    allocate (filled(size(queries, 2)), source=.false.)
    do row = 1, size(queries, 2)
      linear_index = 1_i64
      stride = 1_i64
      do axis = 1, matching_plane_response_input_count
        coordinate_index = 0
        do candidate = 1, axis_sizes(axis)
          if (axes(candidate, axis) == queries(axis, row)) then
            coordinate_index = candidate
            exit
          end if
        end do
        if (coordinate_index == 0) then
          call reject_generator( &
            matching_plane_generator_invalid_grid, &
            'matching-plane query could not be mapped to a canonical axis.', status, message &
            )
          return
        end if
        linear_index = linear_index + int(coordinate_index - 1, i64)*stride
        stride = stride*int(axis_sizes(axis), i64)
      end do
      if (filled(int(linear_index))) then
        call reject_generator( &
          matching_plane_generator_invalid_grid, &
          'matching-plane query grid contains a duplicate coordinate.', status, message &
          )
        return
      end if
      filled(int(linear_index)) = .true.
    end do
    if (.not. all(filled)) then
      call reject_generator( &
        matching_plane_generator_invalid_grid, &
        'matching-plane queries must form one complete Cartesian product.', status, message &
        )
    end if
  end subroutine validate_cartesian_query_grid

  integer(i32) function count_character(text, target) result(count_value)
    character(len=*), intent(in) :: text
    character(len=1), intent(in) :: target
    integer :: index_value

    count_value = 0_i32
    do index_value = 1, len(text)
      if (text(index_value:index_value) == target) count_value = count_value + 1_i32
    end do
  end function count_character

  function integer_text(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write (text, '(i0)') value
  end function integer_text

  subroutine accept_generator(status, message)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = matching_plane_generator_ok
    message = ''
  end subroutine accept_generator

  subroutine reject_generator(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = code
    message = ''
    message = trim(text)
  end subroutine reject_generator

end module bem_matching_plane_response_generator
