!> In-place O(n log n) heapsorts shared by plan builders.
module bem_sorting
  use bem_kinds, only: dp, i32
  implicit none
  private

  public :: sort_real_ascending
  public :: order_by_i32_key

contains

  !> Sort finite real values ascending in place (heapsort; equal values may change order).
  !! @param[inout] values Values to sort.
  subroutine sort_real_ascending(values)
    real(dp), intent(inout) :: values(:)
    real(dp) :: key
    integer :: start, last

    last = size(values)
    do start = last/2, 1, -1
      call sift_down(start, last)
    end do
    do last = size(values), 2, -1
      key = values(last)
      values(last) = values(1)
      values(1) = key
      call sift_down(1, last - 1)
    end do
  contains
    subroutine sift_down(first, limit)
      integer, intent(in) :: first, limit
      integer :: root, child
      real(dp) :: key

      root = first
      key = values(root)
      do while (root <= limit/2)
        child = 2*root
        if (child < limit) then
          if (values(child) < values(child + 1)) child = child + 1
        end if
        if (key >= values(child)) exit
        values(root) = values(child)
        root = child
      end do
      values(root) = key
    end subroutine sift_down
  end subroutine sort_real_ascending

  !> Indices that list `keys` in ascending order (heapsort; equal keys may appear in any order).
  !! @param[in] keys Integer keys.
  !! @param[out] order Permutation with `keys(order)` ascending.
  subroutine order_by_i32_key(keys, order)
    integer(i32), intent(in) :: keys(:)
    integer(i32), intent(out) :: order(:)
    integer(i32) :: item
    integer :: start, last, i

    if (size(order) /= size(keys)) error stop 'order_by_i32_key size mismatch.'
    order = [(int(i, i32), i=1, size(keys))]
    last = size(order)
    do start = last/2, 1, -1
      call sift_down(start, last)
    end do
    do last = size(order), 2, -1
      item = order(last)
      order(last) = order(1)
      order(1) = item
      call sift_down(1, last - 1)
    end do
  contains
    subroutine sift_down(first, limit)
      integer, intent(in) :: first, limit
      integer :: root, child
      integer(i32) :: item

      root = first
      item = order(root)
      do while (root <= limit/2)
        child = 2*root
        if (child < limit) then
          if (keys(order(child)) < keys(order(child + 1))) child = child + 1
        end if
        if (keys(item) >= keys(order(child))) exit
        order(root) = order(child)
        root = child
      end do
      order(root) = item
    end subroutine sift_down
  end subroutine order_by_i32_key

end module bem_sorting
