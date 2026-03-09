!> `bem_app_config_parser` の文字列パース補助手続きを実装する submodule。
submodule (bem_app_config_parser) bem_app_config_parser_parse_utils
  implicit none
contains

  !> `key = value` 形式の行を分割してキーと値を返す。
  module procedure split_key_value
    integer :: p

    p = index(line, '=')
    if (p <= 0) then
      key = ''
      value = ''
      return
    end if
    key = lower(trim(adjustl(line(:p - 1))))
    value = trim(adjustl(line(p + 1:)))
  end procedure split_key_value

  !> 値文字列を倍精度実数へ変換する。
  module procedure parse_real
    read (text, *) out
  end procedure parse_real

  !> 値文字列を `integer(i32)` へ変換する。
  module procedure parse_int
    read (text, *) out
  end procedure parse_int

  !> 値文字列を論理値へ変換する。
  module procedure parse_logical
    character(len=32) :: t

    t = lower(trim(adjustl(text)))
    select case (t)
    case ('true', '.true.')
      out = .true.
    case ('false', '.false.')
      out = .false.
    case default
      error stop 'Logical value must be true/false (or .true./.false.).'
    end select
  end procedure parse_logical

  !> 値文字列を引用符除去済みの文字列へ変換する。
  module procedure parse_string
    character(len=:), allocatable :: tmp

    tmp = trim(adjustl(text))
    if (len(tmp) >= 2 .and. tmp(1:1) == '"' .and. tmp(len(tmp):len(tmp)) == '"') then
      tmp = tmp(2:len(tmp) - 1)
    end if
    out = trim(tmp)
  end procedure parse_string

  !> `[x, y, z]` 形式の値文字列を3成分実数配列へ変換する。
  module procedure parse_real3
    character(len=256) :: t

    t = trim(adjustl(text))
    if (t(1:1) == '[') t = t(2:)
    if (t(len_trim(t):len_trim(t)) == ']') t = t(:len_trim(t) - 1)
    read (t, *) out(1), out(2), out(3)
  end procedure parse_real3

  !> 境界条件モード文字列を内部定数へ変換する。
  module procedure parse_boundary_mode
    character(len=64) :: mode

    call parse_string(text, mode)
    select case (trim(lower(mode)))
    case ('open', 'outflow', 'escape')
      out = bc_open
    case ('reflect', 'reflection')
      out = bc_reflect
    case ('periodic')
      out = bc_periodic
    case default
      error stop 'Unknown boundary condition mode in [sim].'
    end select
  end procedure parse_boundary_mode

  !> 行内コメント (`#` 以降) を除去した文字列を返す。
  module procedure strip_comment
    integer :: p

    p = index(line, '#')
    if (p > 0) then
      out = trim(line(:p - 1))
    else
      out = trim(line)
    end if
  end procedure strip_comment

  !> ASCII 英字を小文字化した文字列を返す。
  module procedure lower
    integer :: i, c

    o = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) o(i:i) = achar(c + 32)
    end do
  end procedure lower

  !> 文字列が指定した接尾辞で終わるかを判定する。
  module procedure ends_with
    integer :: ls, lf

    ls = len_trim(s)
    lf = len_trim(suffix)
    if (lf > ls) then
      ends_it = .false.
    else
      ends_it = (s(ls - lf + 1:ls) == suffix(1:lf))
    end if
  end procedure ends_with

end submodule bem_app_config_parser_parse_utils
