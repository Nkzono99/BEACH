!> ASCII 文字列操作ユーティリティ。
module bem_string_utils
  implicit none
  private

  public :: lower_ascii

contains

  !> ASCII 英字を小文字化した文字列を返す。
  pure function lower_ascii(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i, code

    out = s
    do i = 1, len(s)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end function lower_ascii

end module bem_string_utils
