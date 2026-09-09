!> Offline Zhao ツール共通の数値 query CSV 入力。物理条件と格子条件は呼出側が扱う。
module bem_matching_plane_query_io
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_string_utils, only: is_decimal_real_token
  implicit none
  private

  integer, parameter :: query_line_length = 4096
  integer, parameter :: initial_query_capacity = 64
  integer(i32), parameter, public :: matching_plane_query_ok = 0_i32
  integer(i32), parameter, public :: matching_plane_query_io_error = 2_i32
  integer(i32), parameter, public :: matching_plane_query_invalid_grid = 3_i32

  public :: read_matching_plane_query_csv

contains

  !> 指定 header の列数だけ有限の十進実数を読み、入力順の列ベクトルとして返す。
  subroutine read_matching_plane_query_csv(path, expected_header, queries, status, message)
    character(len=*), intent(in) :: path, expected_header
    real(dp), allocatable, intent(out) :: queries(:, :)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp), allocatable :: buffer(:, :), grown(:, :), query(:)
    character(len=query_line_length) :: line
    character(len=:), allocatable :: record
    integer :: unit_id, ios, line_number, row_count, capacity, column_count
    logical :: header_found

    status = matching_plane_query_ok
    message = ''
    column_count = count_character(expected_header, ',') + 1
    capacity = initial_query_capacity
    allocate (queries(column_count, 0), buffer(column_count, capacity), query(column_count))
    row_count = 0
    line_number = 0
    header_found = .false.
    open (newunit=unit_id, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      status = matching_plane_query_io_error
      message = 'could not open matching-plane query CSV: '//trim(path)
      return
    end if

    do
      read (unit_id, '(a)', iostat=ios) line
      if (ios < 0) exit
      line_number = line_number + 1
      if (ios > 0) then
        status = matching_plane_query_io_error
        message = 'failed to read matching-plane query CSV line '//trim(integer_text(line_number))//'.'
        exit
      end if
      record = trim(adjustl(line))
      if (len(record) == 0) cycle
      if (record(1:1) == '#') cycle
      if (.not. header_found) then
        if (record /= expected_header) then
          status = matching_plane_query_invalid_grid
          message = 'matching-plane query CSV header does not match the required columns.'
          exit
        end if
        header_found = .true.
        cycle
      end if
      if (count_character(record, ',') /= column_count - 1) then
        status = matching_plane_query_invalid_grid
        message = 'matching-plane query row '//trim(integer_text(line_number))// &
                  ' must contain exactly '//trim(integer_text(column_count))//' values.'
        exit
      end if
      call parse_query_record(record, query, ios)
      if (ios /= 0 .or. any(.not. ieee_is_finite(query))) then
        status = matching_plane_query_invalid_grid
        message = 'matching-plane query row '//trim(integer_text(line_number))//' contains invalid values.'
        exit
      end if
      if (row_count == capacity) then
        if (capacity > huge(capacity)/2) then
          status = matching_plane_query_invalid_grid
          message = 'matching-plane query row capacity overflowed.'
          exit
        end if
        allocate (grown(column_count, 2*capacity))
        grown(:, :capacity) = buffer
        call move_alloc(grown, buffer)
        capacity = 2*capacity
      end if
      row_count = row_count + 1
      buffer(:, row_count) = query
    end do
    close (unit_id)
    if (status /= matching_plane_query_ok) return

    if (.not. header_found .or. row_count == 0) then
      status = matching_plane_query_invalid_grid
      message = 'matching-plane query CSV must contain the required header and at least one row.'
      return
    end if
    queries = buffer(:, :row_count)
  end subroutine read_matching_plane_query_csv

  !> List-directed repeat/null/slash syntax を拒否し、十進実数 token だけを読む。
  subroutine parse_query_record(record, query, ios)
    character(len=*), intent(in) :: record
    real(dp), intent(out) :: query(:)
    integer, intent(out) :: ios

    character(len=:), allocatable :: token
    integer :: column, comma, first, last

    query = 0.0_dp
    ios = 0
    first = 1
    do column = 1, size(query)
      if (column < size(query)) then
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

  pure integer function count_character(text, target) result(count_value)
    character(len=*), intent(in) :: text
    character(len=1), intent(in) :: target
    integer :: index_value

    count_value = 0
    do index_value = 1, len(text)
      if (text(index_value:index_value) == target) count_value = count_value + 1
    end do
  end function count_character

  function integer_text(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write (text, '(i0)') value
  end function integer_text

end module bem_matching_plane_query_io
